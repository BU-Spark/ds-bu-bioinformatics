"""Run real-source cohort CLIs in isolated output roots; preserve per-cohort logs."""
import argparse
import concurrent.futures
import json
import os
from pathlib import Path
import subprocess
import sys
import time

COHORTS = [
    ('taylor', 'taylor_pipeline', 'taylor-v2'), ('cambridge', 'cambridge_pipeline', 'cambridge-v2'),
    ('stockholm', 'stockholm_pipeline', 'stockholm-v1'), ('long', 'long_pipeline', 'long-v1'),
    ('mortensen', 'mortensen_pipeline', 'mortensen-v1'), ('tcga', 'tcga_pipeline', 'tcga-prad-v2'),
    ('sboner', 'sboner_pipeline', 'sboner-v2'), ('belfast', 'belfast_pipeline', 'belfast-v1'),
    ('sun', 'sun_pipeline', 'sun-v1'), ('cancermap', 'cancermap_pipeline', 'cancermap-v1'),
    ('cpc-gene', 'cpc_gene_pipeline', 'cpc-gene-v1'), ('moscow', 'moscow_pipeline', 'moscow-v1'),
    ('nakagawa', 'nakagawa_pipeline', 'nakagawa-v1'), ('cit', 'cit_pipeline', 'cit-v1'),
    ('dkfz', 'dkfz_eopc_pipeline', 'dkfz-eopc-v1'), ('cpgea', 'cpgea_pipeline', 'cpgea-v1'),
    ('erho', 'erho_pipeline', 'erho-v1'), ('karnes', 'karnes_pipeline', 'karnes-v1'),
    ('den', 'den_pipeline', 'den-v1'), ('klein', 'klein_pipeline', 'klein-v1'),
]
LARGE = {'erho', 'karnes', 'den', 'klein'}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=Path('data/full-validation'))
    parser.add_argument('--cohort', action='append')
    parser.add_argument('--workers', type=int, default=2)
    parser.add_argument('--offline', action='store_true')
    parser.add_argument('--public-sources', action='store_true',
                        help='Use the explicit public DKFZ subset and prepared public CPGEA derivative; see docs/authorized_rnaseq.md')
    parser.add_argument('--timeout', type=int, default=7200,
                        help='Maximum seconds per build/check (default: 7200; cold downloads can be large)')
    args = parser.parse_args()
    known = {key for key, _, _ in COHORTS}
    if args.cohort and set(args.cohort) - known:
        parser.error('Unknown cohort: ' + ', '.join(sorted(set(args.cohort) - known)))
    if args.workers < 1 or args.timeout < 1:
        parser.error('--workers and --timeout must be positive')
    repo = Path(__file__).resolve().parents[1]
    output = args.output.resolve()
    logs = output / 'logs'
    logs.mkdir(parents=True, exist_ok=True)
    env = {**os.environ, 'PYTHONPATH': str(repo / 'src'), 'MPLCONFIGDIR': '/tmp/ds-bio-matplotlib',
           'OPENBLAS_NUM_THREADS': '1', 'OMP_NUM_THREADS': '1', 'VECLIB_MAXIMUM_THREADS': '1'}
    jobs = [c for c in COHORTS if not args.cohort or c[0] in args.cohort]
    if args.public_sources:
        jobs = [(key, 'dkfz_public_pipeline', 'dkfz-eopc-public-recovery-v1') if key=='dkfz'
                else (key,module,'cpgea-public-v1') if key=='cpgea'
                else (key,module,'cit-fullsource-v1') if key=='cit' else (key,module,snapshot)
                for key,module,snapshot in jobs]
    results = []

    def run(spec):
        key, module, snapshot = spec
        log = logs / (key + '.log')
        command = [sys.executable, '-m', 'ds_bu_bioinformatics.' + module, '--data-dir', str(repo / 'data'),
                   '--snapshot', snapshot, '--output-root', str(output / 'processed'),
                   '--interim-root', str(output / 'interim'), '--anndata-root', str(output / 'anndata')]
        if module == 'tcga_pipeline':
            position = command.index('--interim-root')
            del command[position:position+2]
        if key=='cpgea' and args.public_sources:
            prepared=repo/'data/raw/cpgea-public-discovery/prepared'
            command.extend(['--expression-file',str(prepared/'expression.tsv'),
                            '--sample-file',str(prepared/'samples.tsv'),
                            '--source-config',str(prepared/'source-config.json')])
        task_env = dict(env)
        rscript = repo / '.codex-build/cit-r/bin/Rscript'
        if key == 'cit' and rscript.exists():
            task_env['CIT_RSCRIPT'] = str(rscript)
        if args.offline:
            command.append('--offline')
        started = time.monotonic()
        row = {'cohort': key, 'snapshot': snapshot, 'command': command, 'log': str(log),
               'started_utc': time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())}
        with log.open('w') as handle:
            try:
                proc = subprocess.run(command, cwd=repo, env=task_env, stdout=handle,
                                      stderr=subprocess.STDOUT, timeout=args.timeout)
                row['returncode'] = proc.returncode
                if proc.returncode == 0:
                    checked = subprocess.run([sys.executable, '-m', 'ds_bu_bioinformatics.output_check',
                        '--data-dir', str(output), '--snapshot', snapshot,
                        '--json', str(logs / (key + '-check.json'))], cwd=repo, env=env,
                        stdout=handle, stderr=subprocess.STDOUT, timeout=args.timeout)
                    row['checker_returncode'] = checked.returncode
            except subprocess.TimeoutExpired:
                row['returncode'] = 'timeout'
        row['seconds'] = round(time.monotonic() - started, 2)
        row['tail'] = log.read_text(errors='replace').splitlines()[-8:]
        return row

    def record(row):
        results.append(row)
        (output / 'attempts.json').write_text(json.dumps(results, indent=2) + '\n')
        print(json.dumps({k: row[k] for k in ['cohort', 'returncode', 'seconds', 'tail']}), flush=True)

    # Large exon matrices run alone to avoid multiplying multi-GB resident sets.
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.workers) as pool:
        for future in concurrent.futures.as_completed([pool.submit(run, c) for c in jobs if c[0] not in LARGE]):
            record(future.result())
    for spec in jobs:
        if spec[0] in LARGE:
            record(run(spec))


if __name__ == '__main__':
    main()
