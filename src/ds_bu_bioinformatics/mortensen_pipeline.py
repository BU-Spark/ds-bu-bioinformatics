"""CLI for the Mortensen GSE46602 cohort."""
from . import mortensen
from .geo_cohort_pipeline import build as _build, main as _main


def build(cache, output_root, interim_root, anndata_root):
    return _build(cache, output_root, interim_root, anndata_root, cohort=mortensen)


def main(argv=None):
    return _main(mortensen, argv)


if __name__ == '__main__':
    main()
