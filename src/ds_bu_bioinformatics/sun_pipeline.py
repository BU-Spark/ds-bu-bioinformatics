"""CLI for Sun/Goodison GSE25136 recurrence-status-only export."""
from . import sun
from .geo_cohort_pipeline import build as _build, main as _main


def build(cache, output_root, interim_root, anndata_root):
    return _build(cache, output_root, interim_root, anndata_root, cohort=sun)


def main(argv=None):
    return _main(sun, argv)


if __name__ == '__main__':
    main()
