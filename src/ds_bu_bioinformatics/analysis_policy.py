"""Separate source-record completeness from permission for ordinary Cox analysis."""

EXPERIMENTAL = {
    'GSE94767', 'GSE107299', 'GSE46691', 'GSE62116', 'GSE72291', 'GSE62667',
    'GSE229904', 'GSE116918', 'GSE25136', 'GSE10645', 'EGAD00001004791', 'PRJCA001124',
}
DESIGNS = {'GSE46691': 'nested case-control', 'GSE62116': 'case-cohort',
           'GSE10645': 'selected clinical outcome groups', 'GSE16560': 'outcome-selected extreme groups'}


def policy(accession):
    if accession == 'GSE16560':
        return {'validation_status': 'existing workflow; restricted analysis design',
                'sampling_design': DESIGNS[accession], 'ordinary_cox_allowed': False,
                'reason': 'Outcome-selected lethal/indolent groups require an analysis appropriate to sampling'}
    if accession not in EXPERIMENTAL:
        return {'validation_status': 'existing workflow', 'ordinary_cox_allowed': True}
    design = DESIGNS.get(accession, 'not independently verified')
    return {'validation_status': 'experimental; full-source validation pending',
            'sampling_design': design, 'ordinary_cox_allowed': False,
            'reason': 'Full-source sample, endpoint and sampling-design validation pending; '
                      'clinical record completeness does not authorize ordinary Cox regression',
            'patient_overlap': 'Not established across studies; do not pool as independent cohorts'}
