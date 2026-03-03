"""
.github/scripts/upgrade_readme.py

Reads the current README.md and the full list of exported functions (from
EXPORTS_JSON), then asks the AI to make one focused, cohesive improvement
to the README based on the week's section focus.

Required environment variables:
  OPENAI_API_KEY  - your OpenAI secret key
  EXPORTS_JSON    - JSON string produced by extract_exports.py
"""

import os
import sys
import json
from openai import OpenAI

README_PATH = "README.md"

PROMPT_TEMPLATE = """\
You are a technical writer helping improve the README for an R package.

## Current README

```markdown
{readme_content}
```

## Exported Functions

Below are all the exported functions in this package, along with their \
current roxygen2 documentation:

{exports_block}

## Your Task

Make a single, cohesive improvement to the README focused on: **{section_focus}**

Guidelines:
- Improve or add the "{section_focus}" section using the exported functions above as your source of truth.
- Write clean, idiomatic R code examples that actually work.
- Keep the tone friendly and practical — help someone get up and running quickly.
- Do not remove or substantially alter any existing sections outside your focus area.
- Do not add placeholder text like "coming soon" or "TODO".
- Return ONLY the full updated README.md contents, with no explanation or markdown fences around it.
"""


def format_exports_block(exports: list[dict]) -> str:
    blocks = []
    for exp in exports:
        blocks.append(f"### `{exp['name']}` ({exp['file']})\n\n{exp['roxygen']}")
    return "\n\n---\n\n".join(blocks)


def main():
    api_key      = os.environ.get("OPENAI_API_KEY")
    exports_json = os.environ.get("EXPORTS_JSON")

    missing = [k for k, v in {
        "OPENAI_API_KEY": api_key,
        "EXPORTS_JSON": exports_json,
    }.items() if not v]

    if missing:
        print(f"❌ Missing required environment variables: {', '.join(missing)}", file=sys.stderr)
        sys.exit(1)

    exports = json.loads(exports_json)
    section_focus = os.environ.get("SECTION_FOCUS", "Usage Examples")

    if not os.path.exists(README_PATH):
        print(f"❌ {README_PATH} not found in repo root.", file=sys.stderr)
        sys.exit(1)

    with open(README_PATH, "r") as f:
        readme_content = f.read()

    exports_block = format_exports_block(exports)

    prompt = PROMPT_TEMPLATE.format(
        readme_content=readme_content,
        exports_block=exports_block,
        section_focus=section_focus,
    )

    client = OpenAI(api_key=api_key)
    print(f"⏳ Calling OpenAI to improve README (focus: {section_focus}) ...")

    response = client.chat.completions.create(
        model="gpt-4o",
        messages=[{"role": "user", "content": prompt}],
        temperature=0.4,
    )

    updated_readme = response.choices[0].message.content.strip()

    with open(README_PATH, "w") as f:
        f.write(updated_readme)

    print(f"✅ README updated successfully.")


if __name__ == "__main__":
    main()
