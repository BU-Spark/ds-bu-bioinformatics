"""The user-facing contract: every requested cohort has a callable CLI, not just a registry entry."""
import contextlib
import importlib
import io
from pathlib import Path
import tomllib
import unittest

COMMANDS = {
    'tcga-prad', 'taylor-gse21034', 'cambridge-gse70768', 'stockholm-gse70769',
    'cancermap-gse94767', 'cit-e-mtab-6128', 'cpc-gene-gse107299', 'long-gse54460',
    'dkfz-eopc', 'belfast-gse116918', 'sboner-gse16560', 'erho-gse46691',
    'nakagawa-gse10645', 'karnes-gse62116', 'den-gse72291', 'klein-gse62667',
    'cpgea', 'moscow-gse229904', 'mortensen-gse46602', 'sun-gse25136',
}


class PipelineCoverageTest(unittest.TestCase):
    def test_all_twenty_registered_commands_resolve_and_expose_help_without_acquisition(self):
        root = Path(__file__).resolve().parents[1]
        scripts = tomllib.loads((root / 'pyproject.toml').read_text())['project']['scripts']
        self.assertEqual(set(scripts) - {'import-cohort', 'module-scores', 'check-outputs', 'pooled-analysis',
                                         'decipher-comparison'}, COMMANDS)
        for command in sorted(COMMANDS):
            with self.subTest(command=command), contextlib.redirect_stdout(io.StringIO()), \
                 contextlib.redirect_stderr(io.StringIO()), self.assertRaises(SystemExit) as stopped:
                module, function = scripts[command].split(':')
                getattr(importlib.import_module(module), function)(['--help'])
            self.assertEqual(stopped.exception.code, 0, command)


if __name__ == '__main__':
    unittest.main()
