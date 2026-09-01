#!/usr/bin/env bash
set -euo pipefail

src="${1:?Usage: r-to-odt FILE.R}"
base="${src%.R}"

{
    printf '# %s\n\n' "$(basename "$src")"
    printf '```r\n'
    cat "$src"
    printf '\n```\n'
} > "${base}.md"

pandoc "${base}.md" \
    --from=markdown \
    --to=docx \
    --highlight-style=tango \
    --output="${base}.docx"

soffice --headless \
    --convert-to odt \
    --outdir "$(dirname "$base")" \
    "${base}.docx"

echo "Created ${base}.odt"
