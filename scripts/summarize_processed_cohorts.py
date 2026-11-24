"""Summarize observed processed exports without loading their expression matrices."""
import argparse
import csv
import json
from pathlib import Path

import pandas as pd


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root',type=Path,default=Path('data/full-validation-all'))
    args=parser.parse_args();root=args.root
    rows=[];sets=[]
    for folder in sorted((root/'processed').iterdir()):
        if not (folder/'validation.json').is_file():continue
        v=json.loads((folder/'validation.json').read_text());m=json.loads((folder/'metadata.json').read_text())
        if not v.get('structurally_valid'):continue
        genes=pd.read_csv(folder/'genes.csv',dtype=str,keep_default_na=False)
        samples=pd.read_csv(folder/'samples.csv',dtype=str,keep_default_na=False)
        clinical=pd.read_csv(folder/'clinical.csv',keep_default_na=False)
        endpoints=[]
        for col in clinical:
            if not col.endswith('_event') or '__' in col:continue
            e=col[:-6];eligible=clinical[e+'_eligible'].astype(str).str.lower().eq('true')
            events=pd.to_numeric(clinical[col],errors='coerce')
            endpoints.append(f'{e}: eligible={int(eligible.sum())}, events={int((eligible & events.eq(1)).sum())}, known_events={int(events.eq(1).sum())}')
        gs=set(genes.gene_id)-{''};sets.append(gs)
        rows.append({'snapshot':folder.name,'study':m.get('study',''),'profiles':len(samples),
            'features':len(genes),'unique_gene_ids':len(gs),'expression_cells':len(samples)*len(genes),
            'expression_csv_gz_bytes':(folder/'expression.csv.gz').stat().st_size,
            'expression_h5ad_bytes':(root/'anndata'/folder.name/'expression.h5ad').stat().st_size,
            'scale':m['expression'].get('scale',m['expression'].get('transformations',m['expression'].get('source_scale','unspecified'))),'endpoints':'; '.join(endpoints),
            'processed_directory':str(folder.resolve()),'anndata_file':str((root/'anndata'/folder.name/'expression.h5ad').resolve())})
    with (root/'cohort_summary.csv').open('w') as h:
        writer=csv.DictWriter(h,fieldnames=list(rows[0]) if rows else ['snapshot'],lineterminator='\n');writer.writeheader();writer.writerows(rows)
    summary={'exported_cohorts':len(rows),'cohort_profiles':sum(r['profiles'] for r in rows),
             'expression_cells':sum(r['expression_cells'] for r in rows),
             'interpretation':['Profiles are not globally deduplicated patients; cross-study overlap remains unresolved.',
                 'Gene IDs must be harmonized and scales normalized before joint analysis.',
                 'Klein has restricted deposited coverage; CPGEA is a secondary public derivative.',
                 'DKFZ uses the public processed subset of 95 donors, not the full controlled raw deposit.',
                 'Eligible endpoint records do not override cohort sampling/design restrictions.'],
             'cohorts':rows}
    (root/'cohort_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps({k:v for k,v in summary.items() if k!='cohorts'},indent=2))


if __name__=='__main__':main()
