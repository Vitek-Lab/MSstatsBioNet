"""
.github/scripts/pick_target.py

Walks the R/ directory, finds all functions tagged with @export in their
roxygen2 block, picks one at random, and writes the result to GITHUB_OUTPUT
so subsequent workflow steps can read TARGET_FILE and TARGET_FUNCTION.
"""

import os
import re
import random
import sys

SEARCH_DIR = "R"


def find_exported_functions(filepath: str) -> list[str]:
    """
    Returns the names of all functions in an R file whose immediately
    preceding roxygen2 block (#' lines) contains an @export tag.
    Supports both <- and = assignment styles.
    """
    with open(filepath, "r", errors="ignore") as f:
        lines = f.readlines()

    found = []
    for i, line in enumerate(lines):
        m = re.match(r"^(\w+)\s*(?:<-|=)\s*function\s*\(", line)
        if not m:
            continue
        fn_name = m.group(1)
        # Walk backwards through contiguous roxygen2 comment lines
        j = i - 1
        while j >= 0 and re.match(r"^#'", lines[j]):
            if "@export" in lines[j]:
                found.append(fn_name)
                break
            j -= 1

    return found


def main():
    candidates = []

    for root, _dirs, files in os.walk(SEARCH_DIR):
        for fname in files:
            if not fname.endswith(".R"):
                continue
            fpath = os.path.join(root, fname)
            try:
                for fn in find_exported_functions(fpath):
                    candidates.append((fpath, fn))
            except Exception as e:
                print(f"⚠️  Skipping {fpath}: {e}", file=sys.stderr)

    if not candidates:
        print(
            "❌ No @export-tagged functions found in R/.\n"
            "   Make sure your roxygen2 blocks include @export, e.g.:\n"
            "   #' @export\n"
            "   my_fn <- function(...) { ... }",
            file=sys.stderr,
        )
        sys.exit(1)

    target_file, target_function = random.choice(candidates)

    print(f"✅ Selected : {target_file} :: {target_function}")
    print(f"   Pool size: {len(candidates)} exported functions found in R/")

    github_output = os.environ.get("GITHUB_OUTPUT")
    if github_output:
        with open(github_output, "a") as fh:
            fh.write(f"target_file={target_file}\n")
            fh.write(f"target_function={target_function}\n")
    else:
        # Useful when running the script locally
        print(f"TARGET_FILE={target_file}")
        print(f"TARGET_FUNCTION={target_function}")


if __name__ == "__main__":
    main()
