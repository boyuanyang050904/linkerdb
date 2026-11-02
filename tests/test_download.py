import gzip
import json
import random
import threading
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path

import pytest

from linkerdb import download as dl
from linkerdb.config import Settings
from linkerdb.download import (DownloadError, Downloader, NotFoundError, fetch_structure, fetch_structures,
                               is_valid_gzip, is_valid_mmcif_gz, is_valid_pdb_gz)

# incompressible payload (~300 kB) so that interrupted transfers leave real partial data
CIF = gzip.compress(b"data_TEST\n" + random.Random(0).randbytes(150_000).hex().encode(), compresslevel=1)
PDB = gzip.compress(b"HEADER    TEST\nATOM      1  CA  ALA A   1       0.000   0.000   0.000\n")


class Server:
    """Tiny HTTP server with Range/ETag/If-Range support and failure injection."""

    def __init__(self):
        self.files: dict[str, bytes] = {}
        self.etags: dict[str, str] = {}
        self.fail: dict[str, list[str]] = {}   # path -> queue of actions: "503", "truncate", "ignore-range"
        self.log: list[tuple[str, dict]] = []
        srv = self

        class H(BaseHTTPRequestHandler):
            def log_message(self, *a):
                pass

            def _send(self, code, body=b"", headers=None):
                self.send_response(code)
                for k, v in (headers or {}).items():
                    self.send_header(k, v)
                self.send_header("Content-Length", str(len(body)))
                self.end_headers()
                self.wfile.write(body)

            def do_GET(self):
                srv.log.append((self.path, dict(self.headers)))
                actions = srv.fail.get(self.path, [])
                action = actions.pop(0) if actions else None
                if action == "503":
                    return self._send(503)
                if self.path == "/json":
                    return self._send(200, b'{"ok": 1}')
                if self.path not in srv.files:
                    return self._send(404)
                body = srv.files[self.path]
                etag = srv.etags.get(self.path, '"v1"')
                rng = self.headers.get("Range")
                if rng and action != "ignore-range" and self.headers.get("If-Range", etag) == etag:
                    start = int(rng.split("=")[1].split("-")[0])
                    if start >= len(body):
                        return self._send(416, headers={"Content-Range": f"bytes */{len(body)}"})
                    part = body[start:]
                    return self._send(206, part, {"Content-Range": f"bytes {start}-{len(body) - 1}/{len(body)}",
                                                  "ETag": etag})
                if action == "truncate":
                    self.send_response(200)
                    self.send_header("ETag", etag)
                    self.send_header("Content-Length", str(len(body)))
                    self.end_headers()
                    self.wfile.write(body[: len(body) // 2])
                    self.wfile.flush()
                    self.close_connection = True
                    return
                return self._send(200, body, {"ETag": etag})

            def do_POST(self):
                srv.log.append((self.path, dict(self.headers)))
                actions = srv.fail.get(self.path, [])
                if actions and actions.pop(0) == "503":
                    return self._send(503)
                n = int(self.headers.get("Content-Length", 0))
                payload = json.loads(self.rfile.read(n))
                return self._send(200, json.dumps({"echo": payload}).encode())

        self.httpd = ThreadingHTTPServer(("127.0.0.1", 0), H)
        self.url = f"http://127.0.0.1:{self.httpd.server_address[1]}"
        self.thread = threading.Thread(target=self.httpd.serve_forever, daemon=True)
        self.thread.start()

    def close(self):
        self.httpd.shutdown()

    def hits(self, path):
        return [h for p, h in self.log if p == path]


@pytest.fixture
def server():
    s = Server()
    yield s
    s.close()


@pytest.fixture
def downloader():
    return Downloader(timeout=5, retries=2, backoff=0, sleep=lambda s: None)


def test_fetch_and_cache(server, downloader, tmp_path):
    server.files["/a.cif.gz"] = CIF
    dest = tmp_path / "a.cif.gz"
    assert downloader.fetch(server.url + "/a.cif.gz", dest, is_valid_mmcif_gz) == dest
    assert dest.read_bytes() == CIF
    downloader.fetch(server.url + "/a.cif.gz", dest, is_valid_mmcif_gz)
    assert len(server.hits("/a.cif.gz")) == 1  # second call served from cache
    assert not list(tmp_path.glob("*.part*"))


def test_resume_from_partial_file(server, downloader, tmp_path):
    url = server.url + "/a.cif.gz"
    server.files["/a.cif.gz"] = CIF
    dest = tmp_path / "a.cif.gz"
    import hashlib
    part = tmp_path / f"a.cif.gz.{hashlib.sha1(url.encode()).hexdigest()[:8]}.part"
    part.write_bytes(CIF[:100])
    Path(str(part) + ".meta").write_text(json.dumps({"url": url, "validator": '"v1"'}))
    downloader.fetch(url, dest, is_valid_mmcif_gz)
    hdr = server.hits("/a.cif.gz")[0]
    assert hdr["Range"] == "bytes=100-" and hdr["If-Range"] == '"v1"'
    assert dest.read_bytes() == CIF


def test_resume_restarts_when_remote_changed(server, downloader, tmp_path):
    url = server.url + "/a.cif.gz"
    server.files["/a.cif.gz"] = CIF
    server.etags["/a.cif.gz"] = '"v2"'
    import hashlib
    part = tmp_path / f"a.cif.gz.{hashlib.sha1(url.encode()).hexdigest()[:8]}.part"
    part.write_bytes(b"garbage-from-old-version")
    Path(str(part) + ".meta").write_text(json.dumps({"url": url, "validator": '"v1"'}))
    downloader.fetch(url, tmp_path / "a.cif.gz", is_valid_mmcif_gz)
    assert (tmp_path / "a.cif.gz").read_bytes() == CIF


def test_partial_from_other_url_is_discarded(server, downloader, tmp_path):
    url = server.url + "/a.cif.gz"
    server.files["/a.cif.gz"] = CIF
    import hashlib
    part = tmp_path / f"a.cif.gz.{hashlib.sha1(url.encode()).hexdigest()[:8]}.part"
    part.write_bytes(b"xx")
    Path(str(part) + ".meta").write_text(json.dumps({"url": "http://elsewhere/a", "validator": None}))
    downloader.fetch(url, tmp_path / "a.cif.gz", is_valid_mmcif_gz)
    assert "Range" not in server.hits("/a.cif.gz")[0]
    assert (tmp_path / "a.cif.gz").read_bytes() == CIF


def test_truncated_transfer_is_resumed(server, downloader, tmp_path):
    server.files["/a.cif.gz"] = CIF
    server.fail["/a.cif.gz"] = ["truncate"]
    downloader.fetch(server.url + "/a.cif.gz", tmp_path / "a.cif.gz", is_valid_mmcif_gz)
    hits = server.hits("/a.cif.gz")
    assert len(hits) == 2 and "Range" not in hits[0]
    start = int(hits[1]["Range"].split("=")[1].rstrip("-"))
    assert 0 < start <= len(CIF) // 2
    assert (tmp_path / "a.cif.gz").read_bytes() == CIF


def test_server_ignoring_range_restarts_cleanly(server, downloader, tmp_path):
    server.files["/a.cif.gz"] = CIF
    server.fail["/a.cif.gz"] = ["truncate", "ignore-range"]
    downloader.fetch(server.url + "/a.cif.gz", tmp_path / "a.cif.gz", is_valid_mmcif_gz)
    assert (tmp_path / "a.cif.gz").read_bytes() == CIF


def test_416_means_complete(server, downloader, tmp_path):
    url = server.url + "/a.cif.gz"
    server.files["/a.cif.gz"] = CIF
    import hashlib
    part = tmp_path / f"a.cif.gz.{hashlib.sha1(url.encode()).hexdigest()[:8]}.part"
    part.write_bytes(CIF)
    Path(str(part) + ".meta").write_text(json.dumps({"url": url, "validator": '"v1"'}))
    downloader.fetch(url, tmp_path / "a.cif.gz", is_valid_mmcif_gz)
    assert (tmp_path / "a.cif.gz").read_bytes() == CIF


def test_mirror_fallback_and_not_found(server, downloader, tmp_path):
    server.files["/m2/a.cif.gz"] = CIF
    urls = [server.url + "/m1/a.cif.gz", server.url + "/m2/a.cif.gz"]
    downloader.fetch(urls, tmp_path / "a.cif.gz", is_valid_mmcif_gz)
    assert (tmp_path / "a.cif.gz").read_bytes() == CIF
    with pytest.raises(NotFoundError):
        downloader.fetch([server.url + "/nope1", server.url + "/nope2"], tmp_path / "b.cif.gz")


def test_invalid_content_is_rejected(server, downloader, tmp_path):
    server.files["/bad.cif.gz"] = b"not gzip at all"
    with pytest.raises(DownloadError) as exc:
        downloader.fetch(server.url + "/bad.cif.gz", tmp_path / "bad.cif.gz", is_valid_mmcif_gz)
    assert not isinstance(exc.value, NotFoundError)
    assert not (tmp_path / "bad.cif.gz").exists()
    assert not list(tmp_path.glob("*.part"))


def test_plain_mirror_is_gzipped(server, downloader, tmp_path):
    server.files["/a.cif"] = gzip.decompress(CIF)
    dest = downloader.fetch(server.url + "/a.cif", tmp_path / "a.cif.gz", is_valid_mmcif_gz)
    assert gzip.decompress(dest.read_bytes()) == gzip.decompress(CIF)


def test_retry_on_503(server, downloader, tmp_path):
    server.files["/a.cif.gz"] = CIF
    server.fail["/a.cif.gz"] = ["503", "503"]
    downloader.fetch(server.url + "/a.cif.gz", tmp_path / "a.cif.gz", is_valid_mmcif_gz)
    assert len(server.hits("/a.cif.gz")) == 3
    server.fail["/b"] = ["503"] * 5
    server.files["/b"] = CIF
    with pytest.raises(DownloadError):
        downloader.fetch(server.url + "/b", tmp_path / "b")


def test_get_json(server, downloader):
    assert downloader.get_json(server.url + "/json") == {"ok": 1}
    assert downloader.get_json(server.url + "/q", method="POST", payload={"x": 1}) == {"echo": {"x": 1}}
    with pytest.raises(NotFoundError):
        downloader.get_json(server.url + "/missing")
    server.fail["/q"] = ["503"] * 5
    with pytest.raises(DownloadError):
        downloader.get_json(server.url + "/q", method="POST", payload={})


def test_fetch_many(server, downloader, tmp_path):
    server.files["/a.cif.gz"] = CIF
    res = downloader.fetch_many([([server.url + "/a.cif.gz"], tmp_path / "a.cif.gz", is_valid_mmcif_gz),
                                 ([server.url + "/x"], tmp_path / "x", None)], workers=2)
    assert res[tmp_path / "a.cif.gz"] is None
    assert isinstance(res[tmp_path / "x"], NotFoundError)


def test_fetch_structure_mmcif_then_pdb_fallback(server, downloader, tmp_path, monkeypatch):
    import linkerdb.config as cfg
    monkeypatch.setattr(cfg, "MMCIF_MIRRORS", (server.url + "/cif/{id}.cif.gz",))
    monkeypatch.setattr(cfg, "PDB_FORMAT_MIRRORS", (server.url + "/pdb/{lid}.pdb.gz",))
    settings = Settings(data_dir=tmp_path, retries=0)
    server.files["/cif/1ABC.cif.gz"] = CIF
    server.files["/pdb/2abc.pdb.gz"] = PDB
    path, fmt = fetch_structure("1abc", settings, downloader)
    assert fmt == "cif" and path == settings.structure_path("1abc")
    path, fmt = fetch_structure("2abc", settings, downloader)
    assert fmt == "pdb" and path.name == "2abc.pdb.gz"
    with pytest.raises(NotFoundError):
        fetch_structure("3abc", settings, downloader)
    res = fetch_structures(["1abc", "3abc"], settings, downloader, workers=2)
    assert res["1ABC"][1] == "cif" and isinstance(res["3ABC"], NotFoundError)


def test_validators(tmp_path):
    good = tmp_path / "g.gz"
    good.write_bytes(CIF)
    trunc = tmp_path / "t.gz"
    trunc.write_bytes(CIF[:-20])
    assert is_valid_gzip(good) and not is_valid_gzip(trunc)
    assert is_valid_mmcif_gz(good) and not is_valid_mmcif_gz(trunc)
    p = tmp_path / "p.gz"
    p.write_bytes(PDB)
    assert is_valid_pdb_gz(p) and not is_valid_pdb_gz(trunc)
    notcif = tmp_path / "n.gz"
    notcif.write_bytes(gzip.compress(b"hello"))
    assert not is_valid_mmcif_gz(notcif)
    assert dl._content_range_start("bytes 5-9/10") == 5 and dl._content_range_start("junk") is None
    assert dl._content_range_total("bytes */10") == 10 and dl._content_range_total("bytes 0-1/*") is None
    assert dl._content_range_total("junk") is None
