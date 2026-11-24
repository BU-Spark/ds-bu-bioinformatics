"""Dedicated dkfz-eopc processed RNA-seq source pipeline; authorized local acquisition."""
from .authorized_rnaseq import AuthorizedRNASeq, main as _main
from .geo_cohort_pipeline import build as _build


def build(cache, output_root, interim_root, anndata_root, *, expression_file=None, sample_file=None, source_config=None, **kwargs):
    cohort = AuthorizedRNASeq('dkfz-eopc', expression_file=expression_file, sample_file=sample_file, source_config=source_config)
    return _build(cache, output_root, interim_root, anndata_root, cohort=cohort, **kwargs)


def main(argv=None):
    import sys
    args=list(sys.argv[1:] if argv is None else argv)
    if '--public-source' in args:
        from .dkfz_public_pipeline import main as public_main
        return public_main([a for a in args if a!='--public-source'])
    return _main('dkfz-eopc', args)


if __name__ == '__main__':
    main()
