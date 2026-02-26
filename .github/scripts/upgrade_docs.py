"""
.github/scripts/upgrade_docs.py

Reads TARGET_FILE and TARGET_FUNCTION from the environment (set by the
pick_target.py step), calls the OpenAI API to improve the roxygen2
documentation block for that function, and writes the result back to
the same file.

Required environment variables:
  OPENAI_API_KEY   - your OpenAI secret key (set as a GitHub Secret)
  TARGET_FILE      - path to the R file (e.g. R/utils.R)
  TARGET_FUNCTION  - name of the function to document (e.g. my_fn)
"""

import os
import sys
from openai import OpenAI

PROMPT_TEMPLATE = """\
You are a documentation specialist for R packages.

I have an R file called `{target_file}`. \
Please improve the roxygen2 documentation block for the function `{target_function}`.

Here is the full file:

```r
{file_content}
```

Rules:
- Only modify the roxygen2 comment block (#' lines) immediately above `{target_function}`. \
Do NOT change any code.
- Use proper roxygen2 tags: @description, @param (one per parameter, with type in \
\\code{{}} and a clear description), @return, @examples.
- Keep the existing @export tag exactly as-is.
- Examples should be simple, runnable, and realistic.
- Return ONLY the full updated file contents, with no markdown fences or explanation.
"""


def main():
    api_key         = os.environ.get("OPENAI_API_KEY")
    target_file     = os.environ.get("TARGET_FILE")
    target_function = os.environ.get("TARGET_FUNCTION")

    missing = [k for k, v in {
        "OPENAI_API_KEY": api_key,
        "TARGET_FILE": target_file,
        "TARGET_FUNCTION": target_function,
    }.items() if not v]

    if missing:
        print(f"❌ Missing required environment variables: {', '.join(missing)}", file=sys.stderr)
        sys.exit(1)

    with open(target_file, "r") as f:
        file_content = f.read()

    prompt = PROMPT_TEMPLATE.format(
        target_file=target_file,
        target_function=target_function,
        file_content=file_content,
    )

    client = OpenAI(api_key=api_key)
    print(f"⏳ Calling OpenAI for {target_file}::{target_function} ...")

    response = client.chat.completions.create(
        model="gpt-4o",
        messages=[{"role": "user", "content": prompt}],
        temperature=0.3,
    )

    updated_content = response.choices[0].message.content.strip()

    with open(target_file, "w") as f:
        f.write(updated_content)

    print(f"✅ Doc upgrade written to {target_file}")


if __name__ == "__main__":
    main()
