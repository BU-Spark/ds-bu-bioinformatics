"""Immutable, checksummed HTTP inputs for cohort pipelines."""
from __future__ import annotations
import hashlib
import copy
try:
    import fcntl
except ImportError:  # Windows
    fcntl = None
    import msvcrt
from contextlib import contextmanager
import gzip
import http.client
import json
import shutil
import time
import urllib.error
import urllib.request
import zlib
import zipfile
from datetime import datetime, timezone
from pathlib import Path

GEO_URL = 'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE21nnn/GSE21034/matrix/GSE21034-GPL10264_series_matrix.txt.gz'
GPL_URL = 'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL10264&targ=self&form=text&view=full'
CBIO_BASE = 'https://www.cbioportal.org/api/studies/prad_mskcc'
MYGENE_URL = 'https://mygene.info/v3/query'
LFS_POINTER = b'version https://git-lfs'


class DownloadError(ValueError):
    """An incomplete or invalid response that must never enter the cache."""


def validate_download(path: Path, name: str) -> None:
    """Check the transfer's payload before publishing a cached input."""
    try:
        with path.open('rb') as handle:
            if handle.read(len(LFS_POINTER)) == LFS_POINTER:
                # GitHub can answer with the Git LFS pointer instead of the payload.
                raise ValueError('Received a Git LFS pointer, not the file contents')
        if name.endswith('.gz'):
            # Read through the trailer to check both gzip CRC and completeness.
            with gzip.open(path, 'rb') as handle:
                while handle.read(1024 * 1024):
                    pass
        elif name.endswith('.zip'):
            with zipfile.ZipFile(path) as archive:
                if archive.testzip() is not None:
                    raise ValueError('Corrupt ZIP member')
        elif name.endswith('.json'):
            json.loads(path.read_text(encoding='utf-8'))
        elif name.startswith('GPL') and name.endswith('.soft.txt'):
            content = path.read_text(encoding='utf-8')
            platform = name.split('.')[0]
            if f'^PLATFORM = {platform}\n' not in content or '!platform_table_end' not in content:
                raise ValueError(f'Expected a complete {platform} SOFT record')
    except (ValueError, OSError, EOFError, UnicodeError, zlib.error, zipfile.BadZipFile) as error:
        raise DownloadError(f'Invalid downloaded payload for {name}: {error}') from error


def sha256(path: Path) -> str:
    with path.open('rb') as handle:
        return hashlib.file_digest(handle, 'sha256').hexdigest()


def write_json(path: Path, value) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + '.part')
    temporary.write_text(json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + '\n')
    temporary.replace(path)


class SourceCache:
    @contextmanager
    def transaction(self):
        """Serialize snapshot builds; roll back newly supplied local inputs on failure.

        Valid downloads remain reusable. Rejected local bundles/configs/supplements
        do not become immutable accepted inputs. The lock also prevents stale
        manifests in concurrent shared-builder processes.
        """
        self.root.mkdir(parents=True, exist_ok=True)
        with (self.root / '.build.lock').open('a+b') as lock:
            try:
                if fcntl is not None:
                    fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
                else:
                    if lock.tell() == 0:
                        lock.write(b'0')
                        lock.flush()
                    lock.seek(0)
                    msvcrt.locking(lock.fileno(), msvcrt.LK_NBLCK, 1)
            except OSError as error:
                raise ValueError('Another build is using this source snapshot') from error
            if self.manifest_path.exists():
                current = json.loads(self.manifest_path.read_text())
                if current.get('accessions') != self.manifest['accessions']:
                    raise ValueError('Concurrent source accession change')
                self.manifest = current
            before = copy.deepcopy(self.manifest)
            try:
                yield
            except BaseException:
                for name in list(self.manifest['sources']):
                    if name not in before['sources'] and self.manifest['sources'][name]['method'] == 'LOCAL':
                        (self.root / name).unlink(missing_ok=True)
                        del self.manifest['sources'][name]
                if self.manifest['sources']:
                    write_json(self.manifest_path, self.manifest)
                else:
                    self.manifest_path.unlink(missing_ok=True)
                raise
            finally:
                if fcntl is not None:
                    fcntl.flock(lock, fcntl.LOCK_UN)
                else:
                    lock.seek(0)
                    msvcrt.locking(lock.fileno(), msvcrt.LK_UNLCK, 1)

    def __init__(self, root: Path, offline: bool = False, accessions=None,
                 user_agent: str = 'Taylor-GSE21034/2.0'):
        self.root = root
        self.offline = offline
        self.user_agent = user_agent
        self.manifest_path = root / 'manifest.json'
        self.manifest = json.loads(self.manifest_path.read_text()) if self.manifest_path.exists() else {
            'schema_version': 1,
            'accessions': list(accessions) if accessions is not None else ['GSE21034', 'GPL10264', 'prad_mskcc'],
            'sources': {}}
        if accessions is not None and self.manifest['accessions'] != list(accessions):
            raise ValueError('Cached cohort accessions differ; use a new snapshot')

    def fetch(self, name: str, url: str, body: bytes | None = None) -> Path:
        path = self.root / name
        if path.resolve().is_relative_to(self.root.resolve()) is False:
            raise ValueError('Cache path escapes snapshot')
        request_info = {'url': url, 'method': 'POST' if body is not None else 'GET',
                        'body': body.decode() if body is not None else None}
        entry = self.manifest['sources'].get(name)
        if entry:
            if any(entry[key] != value for key, value in request_info.items()):
                raise ValueError(f'Cached request changed: {name}; use a new snapshot')
            if not path.exists() or sha256(path) != entry['sha256'] or path.stat().st_size != entry['bytes']:
                raise ValueError(f'Missing/corrupt cached input: {path}; use a new snapshot')
            return path
        if self.offline:
            raise ValueError(f'Offline input absent from manifest: {path}')
        if path.exists():
            raise ValueError(f'Unmanifested cached input: {path}; use a new snapshot')
        path.parent.mkdir(parents=True, exist_ok=True)
        temporary = path.with_name(path.name + '.part')
        headers = {'User-Agent': self.user_agent, 'Accept': '*/*'}
        if body is not None:
            headers['Content-Type'] = 'application/x-www-form-urlencoded'
        request = urllib.request.Request(url, data=body, headers=headers)
        print(f'Acquiring {name}', flush=True)
        for attempt in range(3):
            try:
                with urllib.request.urlopen(request, timeout=120) as response, temporary.open('wb') as out:
                    while chunk := response.read(1024 * 1024):
                        out.write(chunk)
                    response_metadata = {key: response.headers.get(key) for key in
                                         ['Content-Type', 'Content-Length', 'ETag', 'Last-Modified', 'total-count']}
                    resolved_url = response.url
                size = temporary.stat().st_size
                if not size:
                    raise DownloadError(f'Empty HTTP response: {url}')
                expected = response_metadata['Content-Length']
                if expected is not None and size != int(expected):
                    raise DownloadError(f'Incomplete HTTP response: {url}: {size} of {expected} bytes')
                validate_download(temporary, name)
                temporary.replace(path)
                break
            except (urllib.error.URLError, TimeoutError, ConnectionError,
                    http.client.IncompleteRead, DownloadError):
                temporary.unlink(missing_ok=True)
                if attempt == 2:
                    raise
                time.sleep(attempt + 1)
            finally:
                temporary.unlink(missing_ok=True)
        self.manifest['sources'][name] = {
            **request_info, 'path': name, 'resolved_url': resolved_url,
            'retrieved_at_utc': datetime.now(timezone.utc).isoformat(),
            'sha256': sha256(path), 'bytes': path.stat().st_size,
            'response_metadata': response_metadata}
        write_json(self.manifest_path, self.manifest)
        return path

    def json(self, name: str, url: str, body: bytes | None = None):
        return json.loads(self.fetch(name, url, body).read_text())

    def ingest(self, name: str, source_path: Path | None, *, source: str) -> Path:
        """Freeze explicitly supplied research files; cached inputs replay without originals."""
        path = self.root / name
        if not path.resolve().is_relative_to(self.root.resolve()) or path == self.manifest_path:
            raise ValueError('Local input path escapes cache or replaces manifest')
        request = {'url': source, 'method': 'LOCAL', 'body': None}
        entry = self.manifest['sources'].get(name)
        if entry:
            if any(entry.get(key) != value for key, value in request.items()):
                raise ValueError(f'Cached local input identity changed: {name}; use a new snapshot')
            if not path.is_file() or sha256(path) != entry['sha256'] or path.stat().st_size != entry['bytes']:
                raise ValueError(f'Missing/corrupt cached input: {name}')
            if source_path is not None and source_path.exists() and sha256(source_path) != entry['sha256']:
                raise ValueError(f'Local input changed: {name}; use a new snapshot')
            return path
        if self.offline:
            raise ValueError(f'Offline local input absent from manifest: {name}')
        if source_path is None or not source_path.is_file():
            raise ValueError(f'Local input required: {name}')
        if path.exists():
            raise ValueError(f'Unmanifested cached input: {name}; use a new snapshot')
        path.parent.mkdir(parents=True, exist_ok=True)
        temporary = path.with_name(path.name + '.part')
        try:
            shutil.copyfile(source_path, temporary)
            if not temporary.stat().st_size:
                raise ValueError(f'Empty local input: {name}')
            validate_download(temporary, name)
            digest = sha256(temporary)
            temporary.replace(path)
        finally:
            temporary.unlink(missing_ok=True)
        self.manifest['sources'][name] = {
            **request, 'path': name, 'resolved_url': source,
            'retrieved_at_utc': datetime.now(timezone.utc).isoformat(),
            'sha256': digest, 'bytes': path.stat().st_size, 'response_metadata': {}}
        write_json(self.manifest_path, self.manifest)
        return path

    def verify_all(self):
        for name, entry in self.manifest['sources'].items():
            path = self.root / name
            if not path.is_file() or sha256(path) != entry['sha256']:
                raise ValueError(f'Missing/corrupt cached input: {path}')
