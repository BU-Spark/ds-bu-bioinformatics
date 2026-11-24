"""Validate in-memory and reloaded final Taylor outputs."""
import numpy as np
import pandas as pd

BASELINE = {'primary_profiles': 131, 'recurrence_events': 27, 'disease_free_observations': 104, 'gene_rows': 19307}


def validate(expression, clinical, samples, genes, mapping, audit):
    if expression.empty or not np.isfinite(expression.to_numpy()).all():
        raise ValueError('Empty or non-finite final expression')
    if not expression.index.is_unique or not expression.columns.is_unique:
        raise ValueError('Duplicate expression labels')
    if expression.columns.tolist() != samples.sample_id.tolist() or expression.columns.tolist() != clinical.sample_id.tolist():
        raise ValueError('Expression/sample/clinical order differs')
    if expression.index.astype(str).tolist() != genes.feature_id.astype(str).tolist():
        raise ValueError('Expression/gene annotation order differs')
    if not genes.feature_id.is_unique or set(expression.index.astype(str)) != set(mapping.loc[mapping.selected, 'feature_id'].astype(str)):
        raise ValueError('Selected features lost or duplicated')
    if samples.patient_id.isna().any() or not samples.patient_id.is_unique:
        raise ValueError('Missing/duplicate selected patients')
    if not samples[['sample_id', 'patient_id', 'geo_accession', 'cbio_sample_id']].equals(clinical[['sample_id', 'patient_id', 'geo_accession', 'cbio_sample_id']]):
        raise ValueError('Clinical/sample identity mismatch')
    if mapping.feature_id.nunique() == 0 or mapping.mapping_status.isna().any():
        raise ValueError('Incomplete mapping audit')
    if audit.sample_id.duplicated().any() or set(audit.loc[audit.selected, 'sample_id']) != set(samples.sample_id):
        raise ValueError('Sample audit selection mismatch')
    if not clinical.dfs_event.dropna().isin([0, 1]).all():
        raise ValueError('Invalid DFS event codes')
    eligible = clinical.dfs_eligible
    if clinical.loc[eligible, ['dfs_event', 'dfs_months']].isna().any().any() or not clinical.loc[eligible, 'dfs_months'].gt(0).all():
        raise ValueError('Invalid eligible DFS times/status')
    current = {'primary_profiles': len(samples), 'recurrence_events': int(clinical.dfs_event.eq(1).sum()),
               'disease_free_observations': int(clinical.dfs_event.eq(0).sum()), 'gene_rows': genes.gene_id.nunique()}
    return {'structurally_valid': True, 'counts': current, 'legacy_baseline': BASELINE,
            'baseline_differences': {key: current[key] - value for key, value in BASELINE.items()},
            'retained_feature_rows': len(genes), 'geo_profiles': len(audit), 'source_features': mapping.feature_id.nunique(),
            'feature_dispositions': mapping.drop_duplicates('feature_id').mapping_status.value_counts().to_dict(),
            'eligible_dfs_profiles': int(eligible.sum()),
            'clinical_missingness': clinical.isna().sum().to_dict(),
            'expression_range': [float(expression.to_numpy().min()), float(expression.to_numpy().max())]}


def validate_export(path, expected_expression=None):
    expression = pd.read_csv(path / 'expression.csv.gz', index_col=0, dtype={'feature_id': 'string', 'gene_id': 'string'}, float_precision='round_trip')
    clinical = pd.read_csv(path / 'clinical.csv', dtype={'sample_id': 'string', 'patient_id': 'string', 'geo_accession': 'string', 'cbio_sample_id': 'string'})
    samples = pd.read_csv(path / 'samples.csv', dtype='string')
    # Keep identities in the same pandas dtype for strict comparison.
    for col in ['sample_id', 'patient_id', 'geo_accession', 'cbio_sample_id']:
        clinical[col] = clinical[col].astype('string')
    genes = pd.read_csv(path / 'genes.csv', dtype={'feature_id': 'string', 'gene_id': 'string'})
    mapping = pd.read_csv(path / 'feature_mapping.csv.gz', dtype={'feature_id': 'string'})
    audit = pd.read_csv(path / 'sample_audit.csv')
    if expected_expression is not None:
        if expression.index.astype(str).tolist() != expected_expression.index.astype(str).tolist() or expression.columns.tolist() != expected_expression.columns.tolist():
            raise ValueError('Export changed expression labels')
        if not np.array_equal(expression.to_numpy(), expected_expression.to_numpy()):
            raise ValueError('Export changed deposited expression values')
    return validate(expression, clinical, samples, genes, mapping, audit)
