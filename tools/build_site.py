#!/usr/bin/env python3
"""Build docs/index.html (the GitHub Pages site) from tools/site_template.html and the saved results.

* numbers in the text come from report/numbers.tex (itself generated from results/data/*.csv);
* the two result tables come from simulation/results/summary_*.md;
* figures, the report PDF and the CSVs are copied into docs/assets so the site is self-contained.

Usage:  python tools/build_site.py
"""
from __future__ import annotations

import html
import re
import shutil
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
REPO = "https://github.com/ParthSarda5/mean-vs-median"      # <- edit if the repository is named differently
DOCS, ASSETS = ROOT / "docs", ROOT / "docs" / "assets"
RES = ROOT / "simulation" / "results"


def clean(v: str) -> str:
    """LaTeX-ish value from numbers.tex -> plain HTML text."""
    return v.replace("\\,", ",").replace("--", "–")


def read_numbers() -> dict[str, str]:
    txt = (ROOT / "report" / "numbers.tex").read_text()
    return {k: clean(v) for k, v in re.findall(r"\\newcommand\{\\(\w+)\}\{(.*)\}", txt)}


def read_md_table(path: Path) -> tuple[list[str], list[list[str]]]:
    lines = [l for l in path.read_text().splitlines() if l.startswith("|")]
    cells = lambda l: [c.strip() for c in l.strip().strip("|").split("|")]
    return cells(lines[0]), [cells(l) for l in lines[2:]]


def pretty_label(s: str) -> str:
    if s == "Cauchy":
        return "Cauchy (t<sub>1</sub>)"
    if re.fullmatch(r"t\d+", s):
        return f"t<sub>{s[1:]}</sub>"
    if s.startswith("Tukey"):
        return f"Tukey ε = {float(s.split('=')[1]):.2f}"
    if s.startswith("Shift"):
        return f"Outliers ε = {float(s.split('=')[1]):.2f}"
    return html.escape(s)


def cell(v: str) -> str:
    return {"inf": "∞", "--": "–"}.get(v, v)


def table_are() -> str:
    _, rows = read_md_table(RES / "summary_are.md")
    head = ["Model", "n·Var mean", "n·Var median", "median, theory", "ARE n = 100", "ARE n = 1000", "ARE theory",
            "eff. mean", "eff. median"]
    out = ['<table><caption><b>Table 1.</b> Efficiency of the mean and the median (n·Var at n = 1000; B = 100,000). '
           'Shaded rows: the median is more efficient. “eff.” is the Cramér–Rao bound divided by n·Var at n = 1000 '
           '(entries marginally above 1 are Monte Carlo noise); “–”: the mean has infinite variance, so the simulated ratio is meaningless.</caption>',
           "<thead><tr>" + "".join(f"<th>{h}</th>" for h in head) + "</tr></thead><tbody>"]
    for r in rows:
        theory = r[6]
        median_wins = theory == "inf" or float(theory) > 1
        out.append(f'<tr class="{"med" if median_wins else ""}"><td>{pretty_label(r[0])}</td>'
                   + "".join(f"<td>{cell(v)}</td>" for v in r[1:]) + "</tr>")
    out.append("</tbody></table>")
    return "".join(out)


def table_comp() -> str:
    _, rows = read_md_table(RES / "summary_comparators.md")
    head = ["Model", "mean", "median", "10% trimmed", "Huber", "Hodges–Lehmann"]
    out = ['<table><caption><b>Table 2.</b> n·MSE about θ = 0 at n = 100 (B = 50,000; Monte Carlo error about 1%). '
           'Lowest in bold. The mean has no finite MSE for Cauchy data.</caption>',
           "<thead><tr>" + "".join(f"<th>{h}</th>" for h in head) + "</tr></thead><tbody>"]
    for r in rows:
        if r[0].endswith("theory"):
            out.append('<tr class="theory"><td>N(0,1), theory (n → ∞)</td>' + "".join(f"<td>{v}</td>" for v in r[1:]) + "</tr>")
            continue
        vals = [float(v) for v in r[1:] if v not in ("inf", "--")]
        best = min(vals)
        tds = "".join(f"<td><b>{v}</b></td>" if v not in ("inf", "--") and float(v) == best else f"<td>{cell(v)}</td>"
                      for v in r[1:])
        out.append(f"<tr><td>{pretty_label(r[0])}</td>{tds}</tr>")
    out.append("</tbody></table>")
    return "".join(out)


def copy_assets() -> None:
    (ASSETS / "data").mkdir(parents=True, exist_ok=True)
    for name in ("fig_core.png", "fig_robust.png", "fig_core.pdf", "fig_robust.pdf"):
        shutil.copy2(RES / "figures" / name, ASSETS / name)
    for csv in (RES / "data").glob("*.csv"):
        shutil.copy2(csv, ASSETS / "data" / csv.name)
    pdf = ROOT / "report" / "simulation_report.pdf"
    if pdf.exists():
        shutil.copy2(pdf, ASSETS / "simulation_report.pdf")
    (DOCS / ".nojekyll").write_text("")                                 # serve the folder as-is


def main() -> None:
    page = (ROOT / "tools" / "site_template.html").read_text()
    nums = read_numbers()
    page = page.replace("@@TABLE_ARE@@", table_are()).replace("@@TABLE_COMP@@", table_comp()).replace("@@REPO@@", REPO)
    missing = []
    for key in set(re.findall(r"@@(\w+)@@", page)):
        if key in nums:
            page = page.replace(f"@@{key}@@", nums[key])
        else:
            missing.append(key)
    if missing:
        raise SystemExit(f"placeholders without a value in numbers.tex: {sorted(missing)}")
    copy_assets()
    (DOCS / "index.html").write_text(page)
    print(f"wrote {DOCS / 'index.html'} ({len(page)/1024:.0f} KiB) and copied assets to {ASSETS}")


if __name__ == "__main__":
    main()
