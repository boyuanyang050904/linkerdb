"""HTTP downloader with a local cache, byte-level resume and mirror fallback.

Resume protocol
---------------
Data is streamed into ``<dest>.<urlhash>.part`` (one per mirror URL). If a previous attempt left a partial file,
the next request carries ``Range: bytes=<size>-`` and, when an ETag/Last-Modified
validator was recorded, ``If-Range`` so that a changed remote file is re-sent in
full (HTTP 200) instead of being spliced onto stale bytes. A 206 response is
appended, a 200 response restarts the file. Only after the content validates is
the ``.part`` file atomically renamed to ``dest``; therefore any existing
``dest`` is complete by construction and is served from cache.
"""

from __future__ import annotations

import gzip
import hashlib
import json
import logging
import os
import shutil
import time
import zlib
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from typing import Callable, Iterable, Sequence

import requests

log = logging.getLogger(__name__)

CHUNK = 1 << 16
RETRY_STATUS = {408, 425, 429, 500, 502, 503, 504}
NOT_FOUND_STATUS = {404, 410}

Validator = Callable[[Path], bool]


class DownloadError(RuntimeError):
    """Raised when no mirror produced a valid file."""


class NotFoundError(DownloadError):
    """The resource does not exist on the server (404/410)."""


class _Retryable(Exception):
    pass


# ---------------------------------------------------------------- validators
def is_valid_gzip(path: Path) -> bool:
    """True if ``path`` is a complete gzip stream (CRC and length are checked)."""
    try:
        with gzip.open(path, "rb") as fh:
            while fh.read(CHUNK):
                pass
        return True
    except (OSError, EOFError, zlib.error):
        return False


def is_valid_mmcif_gz(path: Path) -> bool:
    if not is_valid_gzip(path):
        return False
    with gzip.open(path, "rt", errors="replace") as fh:
        head = fh.read(4096).lstrip()
    return head.startswith("data_")


def is_valid_pdb_gz(path: Path) -> bool:
    if not is_valid_gzip(path):
        return False
    with gzip.open(path, "rt", errors="replace") as fh:
        head = fh.read(4096)
    return any(tag in head for tag in ("HEADER", "ATOM  ", "REMARK", "CRYST1"))


def is_nonempty(path: Path) -> bool:
    return path.exists() and path.stat().st_size > 0


# ---------------------------------------------------------------- downloader
class Downloader:
    def __init__(
        self,
        timeout: float = 60.0,
        retries: int = 4,
        backoff: float = 1.0,
        session: requests.Session | None = None,
        sleep: Callable[[float], None] = time.sleep,
    ) -> None:
        self.timeout = timeout
        self.retries = retries
        self.backoff = backoff
        self.session = session or requests.Session()
        self.session.headers.setdefault("User-Agent", "linkerdb/0.1 (+https://github.com/boyuanyang050904/linkerdb)")
        self._sleep = sleep

    # public -------------------------------------------------------------
    def fetch(
        self,
        urls: str | Sequence[str],
        dest: Path | str,
        validator: Validator | None = None,
        force: bool = False,
    ) -> Path:
        """Download the first URL that yields a valid file into ``dest``.

        URLs whose path does not end in ``.gz`` while ``dest`` does are gzip-compressed
        locally after download (PDBe serves plain ``.cif``).
        """
        dest = Path(dest)
        if isinstance(urls, str):
            urls = [urls]
        validator = validator or is_nonempty
        if dest.exists() and not force:
            return dest
        dest.parent.mkdir(parents=True, exist_ok=True)
        errors: list[str] = []
        all_not_found = True
        for url in urls:
            needs_gzip = dest.name.endswith(".gz") and not url.split("?")[0].endswith(".gz")
            # one partial file per URL so that mirrors never overwrite each other's progress
            tag = hashlib.sha1(url.encode()).hexdigest()[:8]
            part = dest.with_name(f"{dest.name}.{tag}.part")
            try:
                self._fetch_with_retries(url, part)
            except NotFoundError as exc:
                errors.append(f"{url}: not found")
                log.debug("%s", exc)
                continue
            except DownloadError as exc:
                all_not_found = False
                errors.append(f"{url}: {exc}")
                continue
            self._meta_path(part).unlink(missing_ok=True)
            if needs_gzip:
                gz_part = dest.with_name(f"{dest.name}.{tag}.gz.part")
                with open(part, "rb") as src, gzip.open(gz_part, "wb") as dst:
                    shutil.copyfileobj(src, dst)
                part.unlink()
                part = gz_part
            if validator(part):
                os.replace(part, dest)
                for stale in dest.parent.glob(dest.name + ".*.part*"):
                    stale.unlink(missing_ok=True)
                return dest
            all_not_found = False
            errors.append(f"{url}: downloaded content failed validation")
            part.unlink(missing_ok=True)
        exc_type = NotFoundError if all_not_found else DownloadError
        raise exc_type(f"could not fetch {dest.name}: " + "; ".join(errors))

    def fetch_many(
        self,
        jobs: Iterable[tuple[Sequence[str], Path, Validator | None]],
        workers: int = 8,
    ) -> dict[Path, Exception | None]:
        """Run several fetches concurrently; returns ``{dest: None | exception}``."""
        jobs = list(jobs)
        results: dict[Path, Exception | None] = {}

        def run(job):
            urls, dest, validator = job
            try:
                self.fetch(urls, dest, validator)
                return Path(dest), None
            except Exception as exc:  # noqa: BLE001 - reported to caller
                return Path(dest), exc

        with ThreadPoolExecutor(max_workers=max(1, workers)) as pool:
            for dest, err in pool.map(run, jobs):
                results[dest] = err
        return results

    def get_json(self, url: str, method: str = "GET", payload: dict | None = None) -> dict | list:
        """Small helper with the same retry policy for JSON APIs."""
        last: Exception | None = None
        for attempt in range(self.retries + 1):
            try:
                if method == "POST":
                    r = self.session.post(url, json=payload, timeout=self.timeout)
                else:
                    r = self.session.get(url, timeout=self.timeout)
                if r.status_code in NOT_FOUND_STATUS:
                    raise NotFoundError(f"{url}: HTTP {r.status_code}")
                if r.status_code in RETRY_STATUS:
                    raise _Retryable(f"HTTP {r.status_code}")
                r.raise_for_status()
                return r.json()
            except NotFoundError:
                raise
            except (requests.ConnectionError, requests.Timeout, _Retryable, ValueError) as exc:
                last = exc
                if attempt < self.retries:
                    self._sleep(self.backoff * 2**attempt)
            except requests.HTTPError as exc:
                raise DownloadError(f"{url}: {exc}") from exc
        raise DownloadError(f"{url}: giving up after {self.retries + 1} attempts: {last}")

    # internals -----------------------------------------------------------
    @staticmethod
    def _meta_path(part: Path) -> Path:
        return part.with_name(part.name + ".meta")

    def _fetch_with_retries(self, url: str, part: Path) -> None:
        last: Exception | None = None
        for attempt in range(self.retries + 1):
            try:
                self._fetch_once(url, part)
                return
            except NotFoundError:
                raise
            except (requests.ConnectionError, requests.Timeout, requests.exceptions.ChunkedEncodingError,
                    _Retryable) as exc:
                last = exc
                log.info("retrying %s after error: %s", url, exc)
                if attempt < self.retries:
                    self._sleep(self.backoff * 2**attempt)
        raise DownloadError(f"giving up after {self.retries + 1} attempts: {last}")

    def _fetch_once(self, url: str, part: Path) -> None:
        offset = part.stat().st_size if part.exists() else 0
        headers = {"Accept-Encoding": "identity"}
        meta_path = self._meta_path(part)
        if offset:
            # Only resume bytes that came from the very same URL; a partial file from
            # another mirror (different compression) must never be spliced.
            try:
                meta = json.loads(meta_path.read_text()) if meta_path.exists() else {}
            except ValueError:
                meta = {}
            if meta.get("url") != url:
                part.unlink(missing_ok=True)
                offset = 0
            else:
                headers["Range"] = f"bytes={offset}-"
                if meta.get("validator"):
                    headers["If-Range"] = meta["validator"]
        with self.session.get(url, headers=headers, stream=True, timeout=self.timeout) as r:
            if r.status_code in NOT_FOUND_STATUS:
                raise NotFoundError(f"{url}: HTTP {r.status_code}")
            if r.status_code == 416 and offset:
                # Range not satisfiable: we already hold the whole object.
                total = _content_range_total(r.headers.get("Content-Range", ""))
                if total is None or total == offset:
                    return
                part.unlink(missing_ok=True)
                raise _Retryable("stale partial file (416)")
            if r.status_code in RETRY_STATUS:
                raise _Retryable(f"HTTP {r.status_code}")
            if r.status_code not in (200, 206):
                raise DownloadError(f"{url}: unexpected HTTP {r.status_code}")
            mode = "wb"
            if r.status_code == 206:
                start = _content_range_start(r.headers.get("Content-Range", ""))
                if start == offset:
                    mode = "ab"
                else:  # server ignored our offset: refuse to splice
                    part.unlink(missing_ok=True)
                    raise _Retryable(f"unexpected Content-Range start {start} (wanted {offset})")
            validator = r.headers.get("ETag") or r.headers.get("Last-Modified")
            if mode == "wb":
                meta_path.write_text(json.dumps({"url": url, "validator": validator}))
            expected = r.headers.get("Content-Length")
            written = 0
            with open(part, mode) as fh:
                for chunk in r.iter_content(CHUNK):
                    if chunk:
                        fh.write(chunk)
                        written += len(chunk)
            if expected is not None and written != int(expected):
                raise _Retryable(f"short read: {written} of {expected} bytes")


def _content_range_start(value: str) -> int | None:
    # "bytes 100-199/49212"
    try:
        return int(value.split()[1].split("-")[0])
    except (IndexError, ValueError):
        return None


def _content_range_total(value: str) -> int | None:
    try:
        total = value.split("/")[1]
        return None if total == "*" else int(total)
    except (IndexError, ValueError):
        return None


# ---------------------------------------------------------------- structures
def structure_urls(pdb_id: str, mirrors: Sequence[str]) -> list[str]:
    lid = pdb_id.lower()
    return [m.format(id=pdb_id.upper(), lid=lid, mid=lid[1:3]) for m in mirrors]


def fetch_structure(pdb_id: str, settings, downloader: Downloader | None = None) -> tuple[Path, str]:
    """Fetch one entry, mmCIF first, PDB format as last resort.

    Returns ``(path, fmt)`` with ``fmt`` in ``{"cif", "pdb"}``.
    """
    from .config import MMCIF_MIRRORS, PDB_FORMAT_MIRRORS

    downloader = downloader or Downloader(timeout=settings.timeout, retries=settings.retries)
    cif_path = settings.structure_path(pdb_id, "cif")
    try:
        return downloader.fetch(structure_urls(pdb_id, MMCIF_MIRRORS), cif_path, is_valid_mmcif_gz), "cif"
    except DownloadError as cif_exc:
        pdb_path = settings.structure_path(pdb_id, "pdb")
        try:
            return downloader.fetch(structure_urls(pdb_id, PDB_FORMAT_MIRRORS), pdb_path, is_valid_pdb_gz), "pdb"
        except DownloadError as pdb_exc:
            exc_type = NotFoundError if isinstance(cif_exc, NotFoundError) and isinstance(
                pdb_exc, NotFoundError) else DownloadError
            raise exc_type(f"{pdb_id}: mmCIF failed ({cif_exc}); PDB format failed ({pdb_exc})") from pdb_exc


def fetch_structures(pdb_ids: Iterable[str], settings, downloader: Downloader | None = None,
                     workers: int | None = None) -> dict[str, tuple[Path, str] | Exception]:
    downloader = downloader or Downloader(timeout=settings.timeout, retries=settings.retries)
    ids = list(dict.fromkeys(i.upper() for i in pdb_ids))
    out: dict[str, tuple[Path, str] | Exception] = {}

    def run(pid):
        try:
            return pid, fetch_structure(pid, settings, downloader)
        except Exception as exc:  # noqa: BLE001 - recorded per entry
            return pid, exc

    with ThreadPoolExecutor(max_workers=max(1, workers or settings.workers)) as pool:
        for pid, res in pool.map(run, ids):
            out[pid] = res
    return out
