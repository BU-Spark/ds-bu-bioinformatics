"""CLI for Belfast GSE116918 diagnostic biopsies."""
from . import belfast
from .geo_cohort_pipeline import build as _build, main as _main


def build(cache, output_root, interim_root, anndata_root):
    return _build(cache, output_root, interim_root, anndata_root, cohort=belfast)


def main(argv=None):
    return _main(belfast, argv)


if __name__ == '__main__':
    main()
