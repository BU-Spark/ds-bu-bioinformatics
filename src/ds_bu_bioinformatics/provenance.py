"""Hash reachable local module dependencies, excluding unrelated cohort adapters."""
import ast
from pathlib import Path

from .sources import sha256


def code_hashes(*modules):
    root = Path(__file__).parent
    pending = [getattr(m, '__name__', m.__class__.__module__).rsplit('.', 1)[-1]
               if not isinstance(m, str) else m.rsplit('.', 1)[-1] for m in modules]
    seen = set()
    while pending:
        name = pending.pop()
        path = root / (name + '.py')
        if name in seen or not path.is_file():
            continue
        seen.add(name)
        tree = ast.parse(path.read_text())
        for node in ast.walk(tree):
            if isinstance(node, ast.ImportFrom):
                if node.level == 1:
                    pending.extend([node.module.split('.')[0]] if node.module else [a.name for a in node.names])
                elif node.module and node.module.startswith('ds_bu_bioinformatics.'):
                    pending.append(node.module.split('.')[1])
    return {name + '.py': sha256(root / (name + '.py')) for name in sorted(seen)}
