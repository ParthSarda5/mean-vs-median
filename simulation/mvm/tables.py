"""Build the LaTeX tables, the prose macros (numbers.tex) and Markdown summaries from results/data.

Nothing in the report is typed by hand: every table entry and every number quoted in the text
is generated here from the CSVs written by the experiments.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from . import config as C
from . import theory as T
from .distributions import Dist
from .estimators import ALL

INF = r"$\infty$"


def _d(kind, nu=0.0, eps=0.0):
    return Dist(kind, nu=nu, eps=eps, scale=C.TUKEY_SCALE, shift=C.SHIFT_LOC)


def f(x, dec=2):
    if x is None or (isinstance(x, float) and np.isnan(x)):
        return "--"
    if np.isinf(x):
        return INF
    return f"{x:.{dec}f}"


def tex_name(label: str) -> str:
    if label == "N(0,1)":
        return "$N(0,1)$"
    if label == "Cauchy":
        return r"Cauchy ($t_1$)"
    if label == "Laplace":
        return "Laplace"
    if label.startswith("t"):
        return f"$t_{{{label[1:]}}}$"
    if label.startswith("Tukey"):
        return rf"Tukey $\varepsilon{{=}}{float(label.split('=')[1]):.2f}$"
    if label.startswith("Shift"):
        return rf"Outliers $\varepsilon{{=}}{float(label.split('=')[1]):.2f}$"
    return label


def _row(df, dist, n, est):
    r = df[(df.dist == dist) & (df.n == n) & (df.estimator == est)]
    return r.iloc[0]


# ---------------------------------------------------------------------------------------------
TABLE1_ROWS = ["N(0,1)", "t10", "t5", "t4", "t3", "t2", "Cauchy", "Laplace", "Tukey eps=0.1"]


def _dist_of(main, label):
    r = main[main.dist == label].iloc[0]
    return _d(r.kind, r.nu, r.eps)


def table_are(main):
    """Table 1: variances, simulated vs theoretical ARE, and efficiency w.r.t. the CRLB."""
    lines = [r"\begin{tabular}{@{}lrrrrrrrr@{}}", r"\toprule",
             r" & \multicolumn{3}{c}{$n\cdot\mathrm{Var}$ at $n{=}1000$} & \multicolumn{3}{c}{ARE $=\mathrm{Var}(\bar X)/\mathrm{Var}(M)$}"
             r" & \multicolumn{2}{c}{eff.\ vs CRLB$^{\ddagger}$}\\",
             r"\cmidrule(lr){2-4}\cmidrule(lr){5-7}\cmidrule(l){8-9}",
             r"Model & mean & median & med.\ thy & $n{=}100$ & $n{=}1000$ & theory & mean & median\\", r"\midrule"]
    max_rel_se = 0.0
    md = []
    for lab in TABLE1_ROWS:
        d = _dist_of(main, lab)
        m1k, md1k = _row(main, lab, 1000, "mean"), _row(main, lab, 1000, "median")
        a100, a1k = _row(main, lab, 100, "median"), _row(main, lab, 1000, "median")
        inf_var = not np.isfinite(d.var)
        thy_are = T.are_median_vs_mean(d)
        crlb = T.crlb_n(d)
        nv_mean = np.inf if inf_var else m1k.n_var
        eff_mean = 0.0 if inf_var else crlb / m1k.n_var
        eff_med = crlb / md1k.n_var
        sim100 = "--" if inf_var else f(a100.var_ratio_vs_mean, 3)
        sim1k = "--" if inf_var else f(a1k.var_ratio_vs_mean, 3)
        if not inf_var:
            max_rel_se = max(max_rel_se, m1k.se_var / m1k["var"], md1k.se_var / md1k["var"])
        bold = (lambda s: rf"\textbf{{{s}}}") if (np.isinf(thy_are) or thy_are > 1) else (lambda s: s)
        cells = [tex_name(lab), f(nv_mean), f(md1k.n_var), f(T.avar_median(d)), bold(sim100) if not inf_var else sim100,
                 bold(sim1k) if not inf_var else sim1k, bold(f(thy_are, 3)), f(eff_mean), f(eff_med)]
        lines.append(" & ".join(cells) + r"\\")
        md.append([lab, f(nv_mean).replace(INF, "inf"), f(md1k.n_var), f(T.avar_median(d)),
                   sim100, sim1k, f(thy_are, 3).replace(INF, "inf"), f(eff_mean), f(eff_med)])
    lines += [r"\bottomrule", r"\end{tabular}"]
    (C.TABLE_DIR / "table_are.tex").write_text("\n".join(lines) + "\n")
    _write_md("summary_are.md", ["Model", "n·Var mean (n=1000)", "n·Var median (n=1000)", "median, theory",
                                 "ARE n=100", "ARE n=1000", "ARE theory", "eff. mean", "eff. median"], md)
    return max_rel_se


COMP_ROWS = ["N(0,1)", "t5", "t3", "Laplace", "Cauchy", "Tukey eps=0.1", "Shift eps=0.1"]


def table_comp(comp):
    """Table 2: n * MSE of all five estimators at n = 100 (lowest in bold)."""
    order = list(ALL)
    head = [r"Model", "mean", "median", r"10\% trim", "Huber", "HL"]
    lines = [r"\begin{tabular}{@{}lrrrrr@{}}", r"\toprule", " & ".join(head) + r"\\", r"\midrule"]
    thy = [T.normal_n_var(e) for e in order]
    lines.append(r"\emph{$N(0,1)$, theory ($n\to\infty$)} & " + " & ".join(f(x) for x in thy) + r"\\")
    md = [["N(0,1) theory"] + [f(x) for x in thy]]
    lines.append(r"\midrule")
    for lab in COMP_ROWS:
        vals = []
        for e in order:
            r = _row(comp, lab, 100, e)
            inf = (e == "mean" and lab == "Cauchy")
            vals.append(np.inf if inf else r.n_mse)
        shown = [f(v) for v in vals]
        best = min(float(x) for x in shown if x not in ("--", INF))        # ties at the printed precision share the bold
        cells = [rf"\textbf{{{x}}}" if x not in ("--", INF) and float(x) == best else x for x in shown]
        lines.append(tex_name(lab) + " & " + " & ".join(cells) + r"\\")
        md.append([lab] + [f(v).replace(INF, "inf") for v in vals])
    lines += [r"\bottomrule", r"\end{tabular}"]
    (C.TABLE_DIR / "table_comp.tex").write_text("\n".join(lines) + "\n")
    _write_md("summary_comparators.md", ["Model", "mean", "median", "10% trim", "Huber", "HL"], md)


def _write_md(name, header, rows):
    out = ["| " + " | ".join(header) + " |", "|" + "|".join(["---"] * len(header)) + "|"]
    out += ["| " + " | ".join(map(str, r)) + " |" for r in rows]
    (C.DATA_DIR.parent / name).write_text("\n".join(out) + "\n")


# ---------------------------------------------------------------------------------------------
def numbers(main, comp, brk, sens, max_rel_se):
    """Write numbers.tex: \\newcommand's for every figure quoted in the prose."""
    M = {}

    def put(name, val):
        M[name] = val

    # normal model
    for n, w in ((10, "Ten"), (25, "TwentyFive"), (1000, "Thousand")):
        put(f"nvarMedNormal{w}", f(_row(main, "N(0,1)", n, "median").n_var, 2))
    put("nvarMeanNormalThousand", f(_row(main, "N(0,1)", 1000, "mean").n_var, 2))
    put("areNormalSim", f(_row(main, "N(0,1)", 1000, "median").var_ratio_vs_mean, 3))
    put("areNormalThy", f(2 / np.pi, 3))
    # t family
    for lab, w in (("t3", "Tthree"), ("t4", "Tfour"), ("t5", "Tfive")):
        put(f"are{w}Sim", f(_row(main, lab, 1000, "median").var_ratio_vs_mean, 3))
        put(f"are{w}Thy", f(T.are_median_vs_mean(_dist_of(main, lab)), 3))
    put("tCross", f(T.t_crossover(), 2))
    # Laplace
    put("areLapHundred", f(_row(main, "Laplace", 100, "median").var_ratio_vs_mean, 2))
    put("areLapThousand", f(_row(main, "Laplace", 1000, "median").var_ratio_vs_mean, 2))
    # Cauchy
    for q, w in ((0.5, "Fifty"), (0.9, "Ninety"), (0.99, "NinetyNine")):
        put(f"cauchyQ{w}Thy", f(np.tan(np.pi * q / 2), 1 if q < 0.99 else 0))
    put("cauchyMeanQfiftyRange", "{}--{}".format(
        f(main[(main.dist == "Cauchy") & (main.estimator == "mean")].q50_abs.min(), 3),
        f(main[(main.dist == "Cauchy") & (main.estimator == "mean")].q50_abs.max(), 3)))
    put("cauchyMedQfiftyTen", f(_row(main, "Cauchy", 10, "median").q50_abs, 2))
    put("cauchyMedQfiftyThousand", f(_row(main, "Cauchy", 1000, "median").q50_abs, 3))
    nv = _row(main, "Cauchy", 1000, "median").n_var
    put("cauchyMedNvarThousand", f(nv, 2))
    put("cauchyMedNvarThy", f(np.pi**2 / 4, 2))
    put("cauchyEffMed", f(2.0 / nv, 2))
    put("cauchyEffThy", f(8 / np.pi**2, 2))
    put("cauchyMeanNvarTen", f"{_row(main, 'Cauchy', 10, 'mean').n_var:.1e}".replace("e+0", r"\times10^{").replace("e+", r"\times10^{") + "}")
    put("cauchyMeanNvarThousand", f"{_row(main, 'Cauchy', 1000, 'mean').n_var:.1e}".replace("e+0", r"\times10^{").replace("e+", r"\times10^{") + "}")
    # Tukey contamination
    for e, w in ((0.05, "Five"), (0.10, "Ten"), (0.20, "Twenty")):
        put(f"areTukey{w}", f(_row(main, f"Tukey eps={e:g}", 1000, "median").var_ratio_vs_mean, 2))
    # asymmetric outliers (n = 100): MSE ratio mean / median and levels
    for e, w in ((0.05, "Five"), (0.10, "Ten")):
        lab = f"Shift eps={e:g}"
        mm, mdn = _row(comp, lab, 100, "mean"), _row(comp, lab, 100, "median")
        put(f"mseShift{w}Mean", f(mm.mse, 2))
        put(f"mseShift{w}Median", f(mdn.mse, 3))
        put(f"ratioShift{w}", f(mm.mse / mdn.mse, 0))
        put(f"biasShift{w}Mean", f(mm.bias, 2))
        put(f"biasShift{w}Median", f(mdn.bias, 2))
    # comparators at n = 100
    word = {"N(0,1)": "Normal", "t5": "Tfive", "t3": "Tthree", "Laplace": "Laplace", "Cauchy": "Cauchy",
            "Tukey eps=0.1": "TukeyTen", "Shift eps=0.1": "ShiftTen", "Shift eps=0.05": "ShiftFive",
            "Shift eps=0.01": "ShiftOne"}
    ew = {"mean": "Mean", "median": "Median", "trim10": "Trim", "huber": "Huber", "hl": "HL"}
    for lab, w in word.items():
        for e in ALL:
            if lab == "Cauchy" and e == "mean":
                continue                                    # infinite variance: no meaningful MSE
            put(f"cmp{w}{ew[e]}", f(_row(comp, lab, 100, e).n_mse, 2))
    # breakdown experiment
    n = C.BREAKDOWN_N
    for e, w in (("mean", "Mean"), ("trim10", "Trim"), ("hl", "HL"), ("median", "Median"), ("huber", "Huber")):
        s = brk[brk.estimator == e]
        k = int(s[s.median_abs > C.BREAKDOWN_THRESHOLD].k.min())
        put(f"bdK{w}", str(k))
        put(f"bdPct{w}", f(100 * k / n, 1))
    # sensitivity plateaus
    for e, w in (("median", "Median"), ("huber", "Huber"), ("trim10", "Trim"), ("hl", "HL")):
        put(f"sc{w}", f(abs(sens[e].iloc[-1]), 2))
    put("icMedian", f(1 / (2 * 0.3989422804), 2))
    # study size / precision
    put("Bmain", f"{C.B_MAIN:,}".replace(",", r"\,"))
    put("Bcomp", f"{C.B_COMP:,}".replace(",", r"\,"))
    put("maxRelSE", f(100 * max_rel_se, 1))
    put("seedValue", str(C.SEED))

    body = "\n".join(f"\\newcommand{{\\{k}}}{{{v}}}" for k, v in sorted(M.items()))
    (C.TABLE_DIR.parent / "numbers.tex").write_text("% AUTO-GENERATED by simulation/mvm/tables.py -- do not edit\n" + body + "\n")
    return M


def make_all():
    C.TABLE_DIR.mkdir(parents=True, exist_ok=True)
    main = pd.read_csv(C.DATA_DIR / "main_grid.csv")
    comp = pd.read_csv(C.DATA_DIR / "comparators.csv")
    brk = pd.read_csv(C.DATA_DIR / "breakdown.csv")
    sens = pd.read_csv(C.DATA_DIR / "sensitivity.csv")
    max_rel_se = table_are(main)
    table_comp(comp)
    numbers(main, comp, brk, sens, max_rel_se)
