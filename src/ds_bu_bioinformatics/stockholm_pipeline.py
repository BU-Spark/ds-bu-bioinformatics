"""Stockholm GSE70769 CLI using the shared Illumina cohort export workflow."""
from . import stockholm
from .cambridge_pipeline import build as _build, main as _main


def build(cache, output_root, interim_root, anndata_root):
    return _build(cache, output_root, interim_root, anndata_root, cohort=stockholm)


def main(argv=None):
    return _main(argv, cohort=stockholm)


if __name__ == '__main__':
    main()
