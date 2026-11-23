"""Native source acquisition and export for cpc_gene."""
from . import cpc_gene as cohort
from .geo_cohort_pipeline import build as _build, main as _main


def build(cache, output_root, interim_root, anndata_root, **kwargs):
    return _build(cache, output_root, interim_root, anndata_root, cohort=cohort, **kwargs)


def main(argv=None):
    return _main(cohort, argv)


if __name__ == '__main__':
    main()
