"""CIT E-MTAB-6128: Affymetrix CEL processing and SDRF clinical alignment."""
from __future__ import annotations

import csv
import concurrent.futures
import os
import re
import shutil
import subprocess
import tempfile
from pathlib import Path

import pandas as pd

from .geo_cohort_pipeline import add_outcomes, decode_bcr
from .sources import sha256

SERIES = 'E-MTAB-6128'
PLATFORM = 'A-GEOD-16686'
STUDY = 'CIT / Kamoun'
SNAPSHOT = 'cit-v1'
SCALE = 'Locally computed log2 RMA core transcript-cluster intensities'
BASE = 'https://www.ebi.ac.uk/biostudies/'
SDRF = SERIES + '.sdrf.txt'
STUDY_URL = BASE + 'api/v1/studies/' + SERIES
FILES_URL = BASE + 'files/' + SERIES + '/'
R_SCRIPT = Path(__file__).with_name('cit_rma.R')
BCR = 'Characteristics[bcr event]'
TIME = 'Characteristics[months to bcr or last news]'


def outcome(record):
    status = record.get('sdrf__bcr_event')
    time = record.get('sdrf__months_to_bcr_or_last_news')
    status = None if str(status) == 'not available' else status
    time = None if str(time) in ('not applicable', 'not available') else time
    return decode_bcr(status, time, {'1': 1, '0': 0})


def read_sdrf(path):
    with path.open(encoding='utf-8-sig', newline='') as f:
        reader = csv.reader(f, delimiter='\t')
        header = next(reader)
        required = ['Source Name', 'Assay Name', 'Array Data File', 'Array Design REF',
                    'Description', 'Characteristics[patient id]', BCR, TIME,
                    'Factor Value[bcr event]', 'Factor Value[months to bcr or last news]']
        if any(header.count(k) != 1 for k in required):
            raise ValueError('Missing or ambiguous CIT SDRF fields')
        positions = {k: header.index(k) for k in required}
        time_column = positions[TIME] + 1
        if header[time_column] != 'Unit[time unit]':
            raise ValueError('Missing CIT BCR time unit column')
        characteristics = [(i, 'sdrf__' + re.sub(r'[^a-z0-9]+', '_', h[16:-1].lower()).strip('_'))
                           for i, h in enumerate(header) if h.startswith('Characteristics[')]
        if len({k for _, k in characteristics}) != len(characteristics):
            raise ValueError('Duplicate CIT characteristics')
        records = []
        for row in reader:
            # Rows with no follow-up omit only the final factor-unit ontology fields.
            if (len(row) == len(header) - 3 and
                    header[-3:] == ['Unit[time unit]', 'Term Source REF', 'Term Accession Number'] and
                    row[-1] == 'not applicable'):
                row += ['', '', '']
            if len(row) != len(header):
                raise ValueError('Malformed CIT SDRF row')
            r = {k: row[i] for k, i in positions.items()}
            if (r['Array Design REF'] != PLATFORM or
                    r[BCR] != r['Factor Value[bcr event]'] or r[TIME] != r['Factor Value[months to bcr or last news]']):
                raise ValueError('CIT platform or clinical factor fields disagree')
            if r[TIME] not in ('not applicable', 'not available', '') and row[time_column] != 'month':
                raise ValueError('CIT BCR time is not recorded in months')
            cel = r['Array Data File']
            if not re.fullmatch(r'OCT_\d+_HuGene2\.CEL', cel) or cel != r['Assay Name'] + '.CEL':
                raise ValueError('Unexpected CIT assay filename')
            label = r['Description']
            if label not in {'Tumoral tissue', 'Non tumoral tissue', 'normal tissue adjacent to tumour'}:
                raise ValueError('Unknown CIT tissue description')
            records.append({'sample_id': r['Source Name'], 'patient_id': r['Characteristics[patient id]'],
                'assay_id': r['Assay Name'], 'cel_file': cel, 'platform_id': PLATFORM,
                'tissue_description': label, 'selected': label == 'Tumoral tissue',
                'selection_reason': 'primary_tumor' if label == 'Tumoral tissue' else 'non_tumoral_tissue',
                **{k: row[i] for i, k in characteristics}})
    audit = pd.DataFrame(records)
    if audit.empty or any(not audit[c].is_unique or audit[c].eq('').any() for c in ('sample_id','assay_id','cel_file')):
        raise ValueError('Missing or duplicate CIT assay identities')
    clinical = audit.loc[audit.selected].reset_index(drop=True)
    if clinical.empty or not clinical.patient_id.is_unique or not clinical.patient_id.str.fullmatch(r'Patient\d+').all():
        raise ValueError('Missing or duplicate CIT tumor patient identities')
    clinical = add_outcomes(clinical, outcome)
    samples = clinical[['sample_id','patient_id','assay_id','cel_file','platform_id','tissue_description']].copy()
    return clinical, samples, audit


def cel_sizes(record):
    """Read the published file inventory without treating paths as commands."""
    result = {}
    def visit(node):
        if isinstance(node, dict):
            name = node.get('path')
            if isinstance(name, str) and name.endswith('.CEL'):
                if name in result:
                    raise ValueError('Duplicate CEL in BioStudies inventory')
                result[name] = int(node['size'])
            for v in node.values():
                visit(v)
        elif isinstance(node, list):
            for v in node:
                visit(v)
    visit(record)
    return result


def read_annotations(path, features):
    table = pd.read_csv(path, sep='\t', dtype=str, keep_default_na=False)
    if list(table.columns) != ['PROBEID','ENTREZID','SYMBOL'] or set(table.PROBEID) != set(features):
        raise ValueError('CIT annotation keys do not cover the RMA rows')
    rows = []
    for feature, group in table.groupby('PROBEID', sort=False):
        pairs = set(zip(group.ENTREZID, group.SYMBOL))
        valid = len(pairs) == 1 and all(g.isdigit() and bool(s.strip()) for g,s in pairs)
        rows.append({'feature_id': feature, 'gene_id': ';'.join(sorted({g for g,_ in pairs if g})),
                     'gene_symbol': ';'.join(sorted({s for _,s in pairs if s})),
                     'annotation_eligible': valid,
                     'mapping_status': 'mapped' if valid else 'missing_or_ambiguous_annotation'})
    return pd.DataFrame(rows).set_index('feature_id').loc[features].reset_index()


def run_rma(paths):
    rscript = shutil.which(os.environ.get('CIT_RSCRIPT', 'Rscript'))
    if not rscript:
        raise ValueError('CIT requires the pinned R/Bioconductor environment; see docs/cit.md')
    with tempfile.TemporaryDirectory(prefix='cit-rma-') as temp:
        directory = Path(temp)
        listing = directory/'cels.txt'
        listing.write_text(''.join(str(p.resolve())+'\n' for p in paths))
        subprocess.run([rscript, '--vanilla', str(R_SCRIPT), str(listing), temp], check=True)
        frame = pd.read_csv(directory/'expression.tsv', sep='\t', index_col=0,
                            dtype={'feature_id':str}, float_precision='round_trip')
        qc = pd.read_csv(directory/'cel_qc.tsv', sep='\t', keep_default_na=False)
        if (qc.columns.tolist() != ['cel_file','cel_qc_pass','cel_qc_reason']
                or qc.cel_file.tolist() != [p.name for p in paths]
                or qc.cel_qc_pass.dtype != bool or qc.cel_qc_reason.eq('').any()):
            raise ValueError('R raw-array QC inventory differs from immutable CEL inputs')
        if frame.empty or not frame.index.is_unique or frame.columns.tolist() != qc.loc[qc.cel_qc_pass,'cel_file'].tolist():
            raise ValueError('RMA expression assay/feature identities differ')
        annotations = read_annotations(directory/'annotations.tsv', frame.index)
        versions = (directory/'versions.txt').read_text().splitlines()
    return frame, annotations, versions, qc


def prepare(cache):
    record = cache.json(SERIES+'.json', STUDY_URL)
    if record.get('accno') != SERIES:
        raise ValueError('Unexpected BioStudies accession')
    clinical, samples, audit = read_sdrf(cache.fetch(SDRF, FILES_URL+SDRF))
    sizes = cel_sizes(record)
    if set(sizes) != set(audit.cel_file):
        raise ValueError('SDRF and BioStudies CEL inventory differ')
    workers = int(os.environ.get('CIT_DOWNLOAD_WORKERS', '4'))
    if not 1 <= workers <= 8:
        raise ValueError('CIT_DOWNLOAD_WORKERS must be between 1 and 8')
    with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as pool:
        paths = list(pool.map(lambda name: cache.fetch(name, FILES_URL + name), audit.cel_file))
    if any(p.stat().st_size != sizes[p.name] for p in paths):
        raise ValueError('CEL size differs from BioStudies inventory')
    full, annotations, versions, qc = run_rma(paths)
    audit = audit.merge(qc, on='cel_file', how='left', validate='one_to_one', sort=False)
    if audit.cel_qc_pass.isna().any():
        raise ValueError('Missing raw-array QC decision')
    failed = ~audit.cel_qc_pass
    audit.loc[failed, 'selected'] = False
    audit.loc[failed, 'selection_reason'] = 'excluded_unreadable_or_invalid_CEL'
    selected = set(audit.loc[audit.selected,'sample_id'])
    clinical = clinical.loc[clinical.sample_id.isin(selected)].reset_index(drop=True)
    samples = samples.loc[samples.sample_id.isin(selected)].reset_index(drop=True)
    expression = full.loc[:, samples.cel_file].copy()
    expression.columns = samples.sample_id.tolist()
    details = {'source_repository':'EMBL-EBI BioStudies / ArrayExpress', 'bcr_source':'SDRF Characteristics fields, cross-checked against Factor Value fields',
               'sample_selection':'Select explicit Tumoral tissue; retain incomplete outcomes; unique deposited patient IDs',
               'normalization':'oligo::rma(target="core") jointly across readable CEL arrays with finite nonnegative raw signals, including non-tumor arrays, before tumor selection; unreadable/invalid source arrays excluded with individual QC audit',
               'rma_script_sha256':sha256(R_SCRIPT), 'r_versions':versions,
               'annotation':'Pinned hugene20sttranscriptcluster.db; require one Entrez/symbol pair per transcript cluster; no averaging across clusters',
               'clinical_time':'SDRF months to BCR or last news; preserve zero and missing outcomes but flag ineligible'}
    return clinical, samples, audit, expression, annotations, {'raw_array_qc.csv':qc}, versions, details


def summarize(clinical, audit, mapping, extra):
    return {'deposited_arrays':len(audit), 'deposited_subjects':int(audit.patient_id.nunique()),
            'rma_input_arrays':int(audit.cel_qc_pass.sum()), 'raw_arrays_excluded_for_qc':int((~audit.cel_qc_pass).sum())}
