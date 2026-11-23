"""CLI for Nakagawa GSE10645 paired DASL gene-target panels."""
from . import nakagawa
from .geo_cohort_pipeline import build as _build, main as _main


def build(cache, output_root, interim_root, anndata_root, **kwargs):
    return _build(cache, output_root, interim_root, anndata_root, cohort=nakagawa, **kwargs)


def main(argv=None):
    return _main(nakagawa, argv)


if __name__ == '__main__':
    main()
