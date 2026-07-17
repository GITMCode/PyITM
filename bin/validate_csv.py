#!/usr/bin/env python3
"""Validate variables.csv and render a visual table of all LaTeX prettynames.

Outputs: variables_prettynames.png
"""

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pyitm.fileio.variables import _load_variable_csv

rows, index = _load_variable_csv()

# ---- Validation checks ----

errors = []

# 1. Unique shortnames
shorts = [r['shortname'] for r in rows]
dupes = set(s for s in shorts if shorts.count(s) > 1)
if dupes:
    errors.append(f"Duplicate shortnames: {dupes}")

# 2. Every row has shortname and longname
for r in rows:
    if not r['shortname']:
        errors.append(f"Empty shortname in row: {r}")
    if not r['longname']:
        errors.append(f"Empty longname for: {r['shortname']}")

# 3. Check for alias collisions (two different rows claiming the same alias)
seen = {}
for r in rows:
    all_names = [r['shortname'], r['longname']] + r['modelnames']
    for name in all_names:
        key = name.lower()
        if key in seen and seen[key] != r['shortname']:
            errors.append(
                f"Alias collision: '{name}' maps to both "
                f"'{seen[key]}' and '{r['shortname']}'"
            )
        seen[key] = r['shortname']

if errors:
    print("ERRORS:")
    for e in errors:
        print(f"  {e}")
else:
    print(f"All {len(rows)} variables passed validation.")
    print(f"  {len(index)} total aliases in index.")

# ---- Render prettynames ----

has_pretty = [r for r in rows if r['prettyname']]
no_pretty = [r for r in rows if not r['prettyname']]

if no_pretty:
    print(f"\n{len(no_pretty)} variables WITHOUT prettynames:")
    for r in no_pretty:
        print(f"  {r['shortname']}")

# Build table data
ncols = 3
nrows_per_col = (len(has_pretty) + ncols - 1) // ncols

fig, axes = plt.subplots(1, ncols, figsize=(14, nrows_per_col * 0.35 + 1))
fig.suptitle(f'variables.csv: {len(has_pretty)} LaTeX prettynames', fontsize=14)

for col_idx, ax in enumerate(axes):
    ax.set_xlim(0, 1)
    ax.set_ylim(0, nrows_per_col + 1)
    ax.axis('off')

    start = col_idx * nrows_per_col
    end = min(start + nrows_per_col, len(has_pretty))

    for i, r in enumerate(has_pretty[start:end]):
        y = nrows_per_col - i
        shortname = r['shortname']
        pretty = f"${r['prettyname']}$"
        unit = f"({r['unit']})" if r['unit'] else ''

        ax.text(0.0, y, shortname, fontsize=9, fontfamily='monospace',
                va='center')
        ax.text(0.45, y, pretty, fontsize=11, va='center')
        ax.text(0.75, y, unit, fontsize=9, va='center', color='gray')

plt.tight_layout()
outfile = 'variables_prettynames.png'
plt.savefig(outfile, dpi=150, bbox_inches='tight')
print(f"\nRendered to: {outfile}")
