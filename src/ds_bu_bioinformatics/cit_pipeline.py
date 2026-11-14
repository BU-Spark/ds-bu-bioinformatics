"""CLI for CIT E-MTAB-6128."""
from . import cit
from .geo_cohort_pipeline import main as _main


def main(argv=None):
    return _main(cit, argv)


if __name__ == '__main__':
    main()
