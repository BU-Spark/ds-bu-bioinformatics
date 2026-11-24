"""CLI for the Long GSE54460 FPKM cohort."""
from . import long_cohort
from .geo_cohort_pipeline import build as _build, main as _main


def build(cache, output_root, interim_root, anndata_root):
    return _build(cache, output_root, interim_root, anndata_root, cohort=long_cohort)


def main(argv=None):
    return _main(long_cohort, argv)


if __name__ == '__main__':
    main()
