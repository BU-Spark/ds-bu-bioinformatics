"""Public processed Gerhauser 2018 RNA subset; no controlled EGA raw-read acquisition."""
import csv

import pandas as pd

from .gene_identifiers import annotate
from .geo_cohort_pipeline import add_outcomes, decode_bcr

SERIES = 'EGAD00001004791'
STUDY = 'DKFZ / EOPC (public Gerhauser 2018 processed subset)'
PLATFORM = 'RNA-seq, deposited RPKM-labelled values'
SNAPSHOT = 'dkfz-eopc-public-v1'
CBIO_STUDY = 'prostate_dkfz_2018'
COMMIT = '86690e1ed9752b1dcd50b5657f5f05eafa4b6b78'
ACCESSIONS = [SERIES, CBIO_STUDY, 'public-processed-subset-v1']
SCALE = 'Unmodified cBioPortal RPKM-labelled RNA-seq values; no added normalization or regional averaging'
BASE = f'https://media.githubusercontent.com/media/cBioPortal/datahub/{COMMIT}/public/{CBIO_STUDY}/'
TEXT_BASE = f'https://raw.githubusercontent.com/cBioPortal/datahub/{COMMIT}/public/{CBIO_STUDY}/'


def outcome(record):
    # The deposited duration column has no declared unit. Retain it verbatim,
    # rather than turning plausible-looking numbers into survival durations.
    return decode_bcr(record.get('source__BCR_STATUS'), None, {'0':0, '1':1})


ENDPOINTS = {'bcr':outcome}


def prepare(cache):
    study_path=cache.fetch('meta_study.txt',TEXT_BASE+'meta_study.txt')
    study=study_path.read_text()
    if f'cancer_study_identifier: {CBIO_STUDY}' not in study or 'pmid: 30537516' not in study:
        raise ValueError('Wrong DKFZ public study; requires Gerhauser Cancer Cell 2018')
    profile_path=cache.fetch('meta_mrna_seq_rpkm.txt',TEXT_BASE+'meta_mrna_seq_rpkm.txt')
    profile=profile_path.read_text()
    if f'cancer_study_identifier: {CBIO_STUDY}' not in profile or 'datatype: CONTINUOUS' not in profile:
        raise ValueError('Wrong public RNA expression profile')
    paths={n:cache.fetch(n,BASE+n) for n in ['data_clinical_patient.txt','data_clinical_sample.txt','data_mrna_seq_rpkm.txt']}
    patients=pd.read_csv(paths['data_clinical_patient.txt'],sep='\t',comment='#',dtype=str,keep_default_na=False)
    specimens=pd.read_csv(paths['data_clinical_sample.txt'],sep='\t',comment='#',dtype=str,keep_default_na=False)
    if not patients.PATIENT_ID.is_unique or not specimens.SAMPLE_ID.is_unique or not set(specimens.PATIENT_ID).issubset(patients.PATIENT_ID):
        raise ValueError('Public patient/sample identifiers are missing or repeated')
    joined=specimens.merge(patients,on='PATIENT_ID',validate='many_to_one',sort=False)
    with paths['data_mrna_seq_rpkm.txt'].open() as h:
        header=next(csv.reader(h,delimiter='\t'))
    if header[:2]!=['Hugo_Symbol','Entrez_Gene_Id'] or len(header)<3 or len(header)!=len(set(header)):
        raise ValueError('Unexpected public expression source columns')
    raw=pd.read_csv(paths['data_mrna_seq_rpkm.txt'],sep='\t',dtype={'Hugo_Symbol':str,'Entrez_Gene_Id':str},
                    keep_default_na=False,na_values=['NA','N/A','NULL','NaN','nan'],float_precision='round_trip')
    if raw.columns.tolist()!=header:
        raise ValueError('Expression parser changed source headers')
    if not set(header[2:]).issubset(joined.SAMPLE_ID):
        raise ValueError('Unmatched public expression assay')
    symbols=raw.Hugo_Symbol.fillna('').astype(str)
    duplicate=symbols.duplicated(False)|symbols.eq('')
    features=[f'{symbol or "unlabelled"}_row{i+1}' if duplicate.iloc[i] else symbol for i,symbol in enumerate(symbols)]
    if len(features)!=len(set(features)):
        raise ValueError('Source feature-key collision')
    expression=raw.iloc[:,2:].astype(float)
    expression.index=pd.Index(features,name='feature_id')
    scope_mask=(joined.SOURCE.eq('ICGC') & joined.AGE_GROUP.eq('EO') & joined.SAMPLE_CLASS.eq('Tumor') & joined.LOCALIZED_TUMOR.eq('TRUE'))
    present=joined.SAMPLE_ID.isin(expression.columns)
    if not scope_mask[present].all():
        raise ValueError('RNA columns include an unverified study group or nonlocalized/nonprimary specimen')
    candidates=joined.loc[scope_mask & present]
    if candidates.empty:
        raise ValueError('No early-onset ICGC primary RNA specimens')
    winners=set(candidates.sort_values('SAMPLE_ID').drop_duplicates('PATIENT_ID').SAMPLE_ID)
    audit=joined.rename(columns={c:'source__'+c for c in joined}).copy()
    audit['sample_id']=joined.SAMPLE_ID
    audit['patient_id']=joined.PATIENT_ID
    audit['geo_accession']=CBIO_STUDY+':'+joined.SAMPLE_ID
    audit['platform_id']=PLATFORM
    audit['title']=joined.SAMPLE_ID
    audit['selected']=joined.SAMPLE_ID.isin(winners)
    audit['selection_reason']='outside_ICGC_early_onset_primary_scope'
    audit.loc[scope_mask & ~present,'selection_reason']='source_sample_absent_from_expression'
    audit.loc[scope_mask & present,'selection_reason']='other_tumor_region_same_patient'
    audit.loc[audit.selected,'selection_reason']='selected_primary_region_lowest_sample_id'
    clinical=add_outcomes(audit.loc[audit.selected].reset_index(drop=True),outcome)
    samples=clinical[['sample_id','patient_id','geo_accession','platform_id','title']].copy()
    source=expression.loc[:,samples.sample_id]
    if source.lt(0).any().any():raise ValueError('Negative deposited RNA expression')
    mappings,candidates_map=annotate(cache,pd.unique(symbols).tolist(),scope='symbol')
    annotations=mappings.set_index('feature_id').loc[symbols].reset_index(drop=True)
    annotations['feature_id']=features
    annotations['source_gene_symbol']=symbols.to_numpy()
    annotations['source_entrez_from_deposit']=raw.Entrez_Gene_Id.fillna('').astype(str).to_numpy()
    supplied=annotations.source_entrez_from_deposit.ne('')
    conflict=supplied & annotations.source_entrez_from_deposit.ne(annotations.gene_id)
    annotations.loc[conflict,'annotation_eligible']=False
    annotations.loc[conflict,'mapping_status']='source_entrez_identity_conflict'
    projection=pd.DataFrame({'feature_id':features,'source_symbol':symbols})
    candidates_map=projection.merge(candidates_map.rename(columns={'feature_id':'source_symbol'}),on='source_symbol',how='left',sort=False)
    details={'sample_selection':'Require exact ICGC, EO, Tumor and localized TRUE source labels; one region per donor by lowest original sample ID, independently of outcomes; never average',
        'patient_identity':'Original cBioPortal SAMPLE_ID to PATIENT_ID table; exact sample-column matching',
        'annotation':'Original Hugo symbols preserved, with row-qualified duplicate feature keys; cached exact human symbol-to-Entrez mapping; ambiguous/retired symbol mappings excluded and audited',
        'expression_kind':'rpkm','study_source':f'https://github.com/cBioPortal/datahub/tree/{COMMIT}/public/{CBIO_STUDY}',
        'source_profile_metadata':profile,'endpoints':{'bcr':{'definition':'Deposited BCR_STATUS 0/1', 'time_origin':'Source label indicates surgery, but duration unit unverified', 'censoring':'Source last follow-up field retained verbatim; no duration inferred', 'unit':'months','duration_available':False}},
        'limitations':['Public processed subset, not controlled raw EGA files: 118 tumor regions from 95 donors; original EGA description lists 96 donors plus controls',
            'The cBioPortal profile description incorrectly refers to nbl cases and mismatched feature/sample totals; actual matrix identities are checked against the Gerhauser 2018 study and primary EO/ICGC clinical metadata',
            'RPKM label retained; no equivalence to raw-count reprocessing or alternate FPKM derivations assumed',
            'No survival durations: deposited duration unit is not explicitly supplied',
            'The 23 additional tumor regions are audited rather than averaged; selected region sensitivity remains a research decision']}
    return clinical,samples,audit,source,annotations,{'gene_mapping_candidates.csv.gz':candidates_map},[study,profile],details


def summarize(clinical,audit,mapping,extra):
    return {'public_expression_regions':int(audit.selection_reason.isin(['selected_primary_region_lowest_sample_id','other_tumor_region_same_patient']).sum()),
            'other_tumor_regions_excluded':int(audit.selection_reason.eq('other_tumor_region_same_patient').sum())}
