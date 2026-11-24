"""Run the existing cohort builder using only exact, verified cached HTTP responses.

Reuse requires identical URL, HTTP method and request body, and verified response
bytes/SHA-256. Responses retain their original acquisition metadata. Missing exact
responses fail offline rather than using a different annotation query.
"""
from pathlib import Path
import importlib
import shutil
import sys

from ds_bu_bioinformatics.sources import SourceCache, sha256, write_json


class ExactAnnotationCache(SourceCache):
    def fetch(self, name, url, body=None):
        if name not in self.manifest['sources'] and name.startswith('mygene/'):
            import json
            request = {'url':url,'method':'POST' if body is not None else 'GET',
                       'body':body.decode() if body is not None else None}
            for manifest in sorted(self.root.parent.glob('*/manifest.json')):
                if manifest.parent.resolve()==self.root.resolve():continue
                peer=json.loads(manifest.read_text())
                for peer_name,entry in peer['sources'].items():
                    if not all(entry.get(k)==v for k,v in request.items()):continue
                    source=manifest.parent/peer_name
                    if not source.is_file() or source.stat().st_size!=entry['bytes'] or sha256(source)!=entry['sha256']:continue
                    target=self.root/name;target.parent.mkdir(parents=True,exist_ok=True)
                    if target.exists():raise ValueError('Unmanifested target annotation file')
                    shutil.copyfile(source,target)
                    if sha256(target)!=entry['sha256']:raise ValueError('Reused response digest mismatch')
                    with self._manifest_lock:
                        self.manifest['sources'][name]={**entry,'path':name}
                        write_json(self.manifest_path,self.manifest)
                    print('Reused exact immutable response',name,'from',manifest.parent.name,flush=True)
                    return super().fetch(name,url,body)
        return super().fetch(name,url,body)


def main():
    key,snapshot,root=sys.argv[1:]
    module=importlib.import_module('ds_bu_bioinformatics.'+key+'_pipeline')
    cohort=module.cohort
    cache=ExactAnnotationCache(Path('data/raw')/snapshot,offline=True,accessions=[cohort.SERIES,cohort.PLATFORM])
    root=Path(root)
    module.build(cache,root/'processed',root/'interim',root/'anndata')


if __name__=='__main__':main()
