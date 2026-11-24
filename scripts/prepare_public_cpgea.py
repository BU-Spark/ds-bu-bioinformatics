"""Prepare an explicitly experimental public CPGEA derivative for the existing CPGEA CLI.

Preserves deposited TPM values; does not infer follow-up or resolve FOXA1 labels.
The original files and every excluded non-Ensembl row remain in the source audit.
"""
import argparse
from pathlib import Path
import re
import urllib.parse
import xml.etree.ElementTree as ET
import zipfile

import pandas as pd

from ds_bu_bioinformatics.sources import SourceCache, sha256, write_json

COMMIT = '13a8137dce2e2f447702354502d4b903dec40080'
BASE = f'https://raw.githubusercontent.com/nationstrong/CPGEA/{COMMIT}/'
FILES = {'mutant_TPM.csv': 'LuXin_6176/FOXA1_mutation_expr_TPM.csv',
         'wildtype_TPM.csv': 'LuXin_6176/FOXA1_normal_expr_TPM.csv',
         'history.txt': 'LuXin_6176/.Rhistory', 'conversion.R': 'LuXin_6176/fpkm2tpm.R',
         'clinical.xlsx': 'supplementary data/Supplementary Data 2 Clinical and pathological information of specimen.xlsx'}


def clinical_table(path):
    ns = {'s': 'http://schemas.openxmlformats.org/spreadsheetml/2006/main'}
    with zipfile.ZipFile(path) as archive:
        strings = [''.join(t.itertext()) for t in ET.fromstring(archive.read('xl/sharedStrings.xml')).findall('s:si', ns)]
        rows = []
        for row in ET.fromstring(archive.read('xl/worksheets/sheet1.xml')).findall('s:sheetData/s:row', ns):
            record = {}
            for cell in row.findall('s:c', ns):
                if cell.find('s:f', ns) is not None:
                    raise ValueError('Clinical formula requires explicit source review')
                value = cell.find('s:v', ns)
                token = '' if value is None else value.text
                if cell.attrib.get('t') == 's':
                    token = strings[int(token)]
                elif cell.attrib.get('t') not in (None, 'n'):
                    raise ValueError('Unexpected clinical cell type')
                record[re.sub(r'\d+$', '', cell.attrib['r'])] = token
            rows.append(record)
    header = rows[0]
    frame = pd.DataFrame([{header[k]: v for k, v in r.items()} for r in rows[1:]]).fillna('')
    if len(frame) != 208 or not frame['Sequencing ID'].is_unique or not frame['Sequencing ID'].str.fullmatch(r'\d+').all():
        raise ValueError('Unexpected original clinical subject inventory')
    return frame


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--cache', type=Path, default=Path('data/raw/cpgea-public-discovery'))
    parser.add_argument('--offline', action='store_true')
    args = parser.parse_args()
    cache = SourceCache(args.cache, offline=args.offline, accessions=['PRJCA001124', 'nationstrong/CPGEA'])
    with cache.transaction():
        paths = {name: cache.fetch(name, BASE + urllib.parse.quote(path, safe='/')) for name, path in FILES.items()}
        frames = [pd.read_csv(paths[name], dtype={'Gene_id': str, 'SYMBOL': str}, keep_default_na=False,
                              float_precision='round_trip') for name in ('mutant_TPM.csv', 'wildtype_TPM.csv')]
        if any(f.columns[:2].tolist() != ['SYMBOL', 'Gene_id'] for f in frames):
            raise ValueError('Unexpected source expression annotation columns')
        pd.testing.assert_frame_equal(frames[0].iloc[:, :2], frames[1].iloc[:, :2])
        if not frames[0].Gene_id.is_unique:
            raise ValueError('Repeated original gene identities; no implicit aggregation')
        all_expression = pd.concat([f.set_index('Gene_id').iloc[:, 1:] for f in frames], axis=1)
        if all_expression.shape != (20615, 134) or not all_expression.columns.is_unique:
            raise ValueError('Unexpected deposited expression cohort inventory')
        if not all(re.fullmatch(r'T\d+', x) for x in all_expression.columns):
            raise ValueError('Unrecognized expression subject label')
        clinical = clinical_table(paths['clinical.xlsx'])
        clinical.insert(0, 'expression_id', 'T' + clinical['Sequencing ID'])
        if not set(all_expression.columns).issubset(clinical.expression_id):
            raise ValueError('Expression subject absent from original clinical supplement')
        clinical['specimen'] = 'primary_tumor'
        group = {name: label for f, label in zip(frames, ['repository_mutant_group', 'repository_wildtype_group']) for name in f.columns[2:]}
        clinical['repository_expression_group'] = clinical.expression_id.map(group).fillna('not_in_expression')
        clinical['repository_group_disagrees_with_clinical_FOXA1'] = (
            clinical.repository_expression_group.eq('repository_mutant_group') & clinical.FOXA1_mutation.eq('wide_type'))
        if not clinical.BCR.isin(['YES', 'NO', 'NA', 'Dead', '']).all():
            raise ValueError('Unreviewed source BCR label')
        # Death in a recurrence column establishes neither recurrence nor its absence.
        # Preserve the original label, with an explicitly missing derived event.
        clinical['bcr_status_for_export'] = clinical.BCR.where(clinical.BCR.isin(['YES', 'NO']), 'NA')
        annotations = frames[0][['Gene_id', 'SYMBOL']].copy()
        annotations['retained_input_namespace'] = annotations.Gene_id.str.fullmatch(r'ENSG\d+(?:\.\d+)?')
        annotations['reason'] = annotations.retained_input_namespace.map({True: 'recognized_Ensembl_gene', False: 'unsupported_source_namespace_preserved_in_audit'})
        expression = all_expression.loc[annotations.loc[annotations.retained_input_namespace, 'Gene_id']]
        output = cache.root / 'prepared'
        output.mkdir(exist_ok=True)
        expression.to_csv(output / 'expression.tsv', sep='\t', index_label='feature_id')
        clinical.to_csv(output / 'samples.tsv', sep='\t', index=False)
        annotations.to_csv(output / 'source_feature_audit.tsv', sep='\t', index=False)
        config = {
            'study_accession': 'PRJCA001124', 'source': f'Experimental public derivative: https://github.com/nationstrong/CPGEA/tree/{COMMIT}; clinical Supplementary Data 2; original study doi:10.1038/s41586-020-2135-x',
            'expression_kind': 'tpm', 'expression_scale': 'Unmodified repository-labelled TPM; repository derives TPM from log2 FPKM before this pipeline; excluded XLOC rows audited, no renormalization after exclusion',
            'sample_linkage': 'Exact T + original Sequencing ID string, all 134 expression labels present among 208 clinical subjects. This derived crosswalk and derivative provenance require scientific review.',
            'specimen_definition': 'Repository primary-tumor expression derivative; mutant and normal filenames describe FOXA1 analysis groups, not a tumor/benign pairing. Ten mutant-group labels disagree with clinical wide_type; both labels retained without adjudication.',
            'sample_columns': {'sample_id': 'expression_id', 'expression_id': 'expression_id', 'patient_id': 'Sequencing ID', 'tissue': 'specimen'},
            'tumor_values': ['primary_tumor'], 'excluded_tissue_values': [], 'replicate_policy': 'error', 'primary_endpoint': 'bcr',
            'endpoints': {'bcr': {'definition': 'Deposited BCR YES/NO only; Dead and NA remain unknown events with original labels preserved; precise endpoint definition unreconciled', 'time_origin': 'Unverified; deposited duration string deliberately not interpreted', 'censoring': 'Unverified; no survival-eligible durations exported', 'event_column': 'bcr_status_for_export', 'event_codes': {'YES': 1, 'NO': 0}, 'time_column': None, 'time_unit': 'months'}},
            'preparation_script_sha256': sha256(Path(__file__)), 'original_source_manifest_sha256': sha256(cache.manifest_path),
            'source_feature_audit_sha256': sha256(output / 'source_feature_audit.tsv'),
            'limitations': ['Public secondary derivative, not a verified replacement for the original processed deposit', 'FOXA1 groups and inferred patient crosswalk require review', 'No survival-duration interpretation', '1998 non-Ensembl rows preserved in original files and the preparation audit, excluded from mapping input']}
        write_json(output / 'source-config.json', config)
        cache.verify_all()
        print(f'Prepared {expression.shape[1]} profiles x {expression.shape[0]} Ensembl rows; {len(annotations)-len(expression)} unsupported rows audited; {int(clinical.repository_group_disagrees_with_clinical_FOXA1.sum())} group discrepancies retained')


if __name__ == '__main__':
    main()
