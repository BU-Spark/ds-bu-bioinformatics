"""Dedicated dkfz-eopc processed RNA-seq source pipeline; authorized local acquisition."""
from .authorized_rnaseq import AuthorizedRNASeq, main as _main
from .geo_cohort_pipeline import build as _build


def build(cache, output_root, interim_root, anndata_root, *, expression_file=None, sample_file=None, source_config=None, **kwargs):
    cohort = AuthorizedRNASeq('dkfz-eopc', expression_file=expression_file, sample_file=sample_file, source_config=source_config)
    return _build(cache, output_root, interim_root, anndata_root, cohort=cohort, **kwargs)


def main(argv=None):
    return _main('dkfz-eopc', argv)


if __name__ == '__main__':
    main()
