#!/usr/bin/env python3
"""Embed the tested helper module into the self-contained Colab notebook."""
import argparse
import hashlib
import json
from pathlib import Path


def expected_source(root):
    code = (root / "workflow_support.py").read_text(encoding="utf-8")
    sha = hashlib.sha256(code.encode()).hexdigest()
    return ("#@title 0) Load workflow helpers (run once)\n"
            "# Generated from workflow_support.py by scripts/sync_notebook.py.\n"
            f"_WORKFLOW_SOURCE_SHA = {sha!r}\n" + code).splitlines(keepends=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    path = root / "workflow.ipynb"
    notebook = json.loads(path.read_text(encoding="utf-8"))
    cells = [cell for cell in notebook["cells"] if cell.get("id") == "workflow-support"]
    if len(cells) != 1:
        raise SystemExit("Expected one workflow-support cell.")
    expected = expected_source(root)
    if args.check:
        if cells[0]["source"] != expected:
            raise SystemExit("Helper cell is stale. Run python scripts/sync_notebook.py.")
        print("Embedded workflow helpers match the tested Python source.")
    else:
        cells[0]["source"] = expected
        path.write_text(json.dumps(notebook, indent=1, ensure_ascii=False) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
