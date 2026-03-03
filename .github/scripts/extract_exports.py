"""
.github/scripts/extract_exports.py

Walks R/, finds every function tagged with @export, and extracts its full
roxygen2 documentation block. Outputs a JSON blob to GITHUB_OUTPUT so that
upgrade_readme.py can consume it without re-parsing the R files.

Also deterministically picks a section focus for this week based on the
ISO week number, so the README improves a different area each week in rotation.
"""

import os
import re
import sys
import json
from datetime import date

SEARCH_DIR = "R"

SECTION_ROTATION = [
    "Getting Started / Installation",
    "Usage Examples",
    "Function Reference Table",
    "Motivation & Use Cases",
    "FAQ / Common Pitfalls",
    "Contributing & Development Setup",
]


def extract_function_docs(filepath: str) -> list[dict]:
    """
    Returns a list of dicts, one per @export-tagged function:
      { "name": str, "file": str, "roxygen": str }
    """
    with open(filepath, "r", errors="ignore") as f:
        lines = f.readlines()

    results = []
    for i, line in enumerate(lines):
        m = re.match(r"^(\w+)\s*(?:<-|=)\s*function\s*\(", line)
        if not m:
            continue
        fn_name = m.group(1)

        # Collect the contiguous roxygen2 block above this line
        roxygen_lines = []
        j = i - 1
        while j >= 0 and re.match(r"^#'", lines[j]):
            roxygen_lines.insert(0, lines[j].rstrip())
            j -= 1

        if any("@export" in l for l in roxygen_lines):
            results.append({
                "name": fn_name,
                "file": filepath,
                "roxygen": "\n".join(roxygen_lines),
            })

    return results


def main():
    all_exports = []

    for root, _dirs, files in os.walk(SEARCH_DIR):
        for fname in sorted(files):
            if not fname.endswith(".R"):
                continue
            fpath = os.path.join(root, fname)
            try:
                all_exports.extend(extract_function_docs(fpath))
            except Exception as e:
                print(f"⚠️  Skipping {fpath}: {e}", file=sys.stderr)

    if not all_exports:
        print(
            "❌ No @export-tagged functions found in R/.\n"
            "   Add #' @export above your exported function definitions.",
            file=sys.stderr,
        )
        sys.exit(1)

    print(f"✅ Found {len(all_exports)} exported functions:")
    for exp in all_exports:
        print(f"   • {exp['name']} ({exp['file']})")

    # Rotate section focus by ISO week number so it's deterministic
    week_number = date.today().isocalendar().week
    section_focus = SECTION_ROTATION[week_number % len(SECTION_ROTATION)]
    print(f"\n📌 This week's README focus: {section_focus}")

    exports_json = json.dumps(all_exports)

    github_output = os.environ.get("GITHUB_OUTPUT")
    if github_output:
        with open(github_output, "a") as fh:
            fh.write(f"exports_json={exports_json}\n")
            fh.write(f"section_focus={section_focus}\n")
    else:
        print(f"\nEXPORTS_JSON={exports_json}")
        print(f"SECTION_FOCUS={section_focus}")


if __name__ == "__main__":
    main()
