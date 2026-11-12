"""Publication figures for the simulation report (light surface, validated categorical palette).

Colour follows the entity: the same estimator keeps the same hue in every panel.  Aqua, yellow
and magenta are below 3:1 contrast on white, so every series also carries its own marker shape
and a legend (the relief rule), and all plotted values are in results/data/*.csv.
"""
from __future__ import annotations

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D

from . import config as C
from . import theory as T
from .distributions import Dist
from .estimators import ALL, SHORT

# ---- palette (reference categorical order; slots 1-5, validated on a white surface) ----------
COL = {"mean": "#2a78d6", "median": "#eb6834", "trim10": "#1baf7a", "huber": "#eda100", "hl": "#e87ba4"}
MARK = {"mean": "o", "median": "s", "trim10": "^", "huber": "D", "hl": "v"}
INK, INK2, MUTED, GRID, BASE = "#0b0b0b", "#52514e", "#898781", "#e1e0d9", "#c3c2b7"
FULL_W = 7.09          # A4 text width with 1.5 cm margins, in inches


def style():
    plt.rcParams.update({
        "figure.dpi": 150, "savefig.dpi": 300, "pdf.fonttype": 42, "ps.fonttype": 42,
        "font.family": "DejaVu Sans", "mathtext.fontset": "dejavusans",
        "font.size": 7.2, "axes.labelsize": 7.2, "axes.titlesize": 7.8, "axes.titleweight": "normal",
        "axes.titlelocation": "left", "axes.titlepad": 5,
        "xtick.labelsize": 6.6, "ytick.labelsize": 6.6, "legend.fontsize": 6.5,
        "text.color": INK, "axes.labelcolor": INK2, "xtick.color": INK2, "ytick.color": INK2,
        "axes.edgecolor": BASE, "axes.linewidth": 0.8, "xtick.major.width": 0.8, "ytick.major.width": 0.8,
        "xtick.major.size": 2.5, "ytick.major.size": 2.5,
        "axes.spines.top": False, "axes.spines.right": False,
        "axes.grid": True, "axes.grid.axis": "y", "grid.color": GRID, "grid.linewidth": 0.6, "grid.linestyle": "-",
        "axes.axisbelow": True, "lines.linewidth": 1.5, "lines.markersize": 3.6,
        "legend.frameon": False, "legend.handlelength": 1.8, "legend.borderaxespad": 0.2,
        "figure.facecolor": "white", "axes.facecolor": "white",
    })


def _panel(ax, letter, title):
    ax.set_title(rf"$\mathbf{{({letter})}}$ {title}")


def _save(fig, name):
    fig.savefig(C.FIG_DIR / f"{name}.pdf", bbox_inches="tight", pad_inches=0.02)
    fig.savefig(C.FIG_DIR / f"{name}.png", bbox_inches="tight", pad_inches=0.02, dpi=220)
    plt.close(fig)


def _ser(df, **kw):
    for k, v in kw.items():
        df = df[df[k] == v]
    return df


# =============================================================================================
def figure_core():
    """(a) normal model, (b) ARE vs t degrees of freedom, (c) Cauchy running estimates."""
    main = pd.read_csv(C.DATA_DIR / "main_grid.csv")
    sweep = pd.read_csv(C.DATA_DIR / "t_sweep.csv")
    paths = pd.read_csv(C.DATA_DIR / "cauchy_paths.csv")
    fig, axs = plt.subplots(1, 3, figsize=(FULL_W, 2.15), gridspec_kw={"wspace": 0.36})

    # (a) n*Var vs n at N(0,1) ----------------------------------------------------------------
    ax = axs[0]
    for est, lim, lab in (("mean", 1.0, "1"), ("median", np.pi / 2, r"$\pi/2$")):
        s = _ser(main, dist="N(0,1)", estimator=est).sort_values("n")
        ax.axhline(lim, color=COL[est], lw=0.9, ls=(0, (4, 3)), alpha=0.9, zorder=1)
        ax.errorbar(s["n"], s["n_var"], yerr=s["n"] * s["se_var"] * 1.96, color=COL[est], marker=MARK[est],
                    capsize=0, zorder=3)
        ax.text(1.05e3, lim + 0.025, lab, color=INK2, va="bottom", ha="right", fontsize=6.6)
    ax.set_xscale("log")
    ax.set_xlim(8, 1.15e3)
    ax.set_xticks([10, 25, 50, 100, 1000])
    ax.set_xticklabels(["10", "25", "50", "100", "1000"])
    ax.xaxis.set_minor_locator(matplotlib.ticker.NullLocator())
    ax.set_ylim(0.85, 1.95)
    ax.set_xlabel("sample size $n$")
    ax.set_ylabel(r"$n\cdot\mathrm{Var}$")
    ax.legend(handles=[Line2D([], [], color=COL[e], marker=MARK[e], lw=1.5, label=l)
                       for e, l in (("mean", "sample mean"), ("median", "sample median"))],
              loc="center", bbox_to_anchor=(0.55, 0.30), handlelength=1.6)     # in the empty band between the two curves
    _panel(ax, "a", "Normal: variance vs $n$")

    # (b) ARE vs nu -----------------------------------------------------------------------------
    ax = axs[1]
    nu = np.geomspace(2.35, 40, 400)
    ax.plot(nu, [T.are_t(v) for v in nu], color=INK, lw=1.5, zorder=2, label="theory (asymptotic)")
    ax.axhline(1.0, color=BASE, lw=0.9, zorder=1)
    ax.axhline(2 / np.pi, color=MUTED, lw=0.8, ls=(0, (4, 3)), zorder=1)
    ax.text(2.42, 2 / np.pi * 0.985, r"$2/\pi$", color=INK2, va="top", ha="left", fontsize=6.6)
    star = T.t_crossover()
    ax.plot([star], [1.0], marker="o", ms=4.5, color="white", mec=INK, mew=1.2, zorder=4)
    ax.annotate(rf"$\nu^*\!\approx\!{star:.2f}$", (star, 1.0), xytext=(7.5, 2.0), fontsize=6.8, color=INK,
                arrowprops=dict(arrowstyle="-", color=MUTED, lw=0.8, shrinkA=0, shrinkB=3))
    for n, fc, ms, z in ((1000, INK2, 3.2, 3), (100, "none", 6.0, 4)):
        s = _ser(sweep, n=n, estimator="median").sort_values("nu")
        ax.errorbar(s["nu"], s["var_ratio_vs_mean"], yerr=1.96 * s["se_var_ratio_vs_mean"], ls="none", marker="o",
                    mfc=fc, mec=INK2, color=INK2, ms=ms, mew=0.9, capsize=0, zorder=z, label=f"simulated, $n={n}$")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(2.3, 40)
    ax.set_ylim(0.5, 7)
    ax.set_xticks([3, 5, 10, 30])
    ax.set_xticklabels(["3", "5", "10", "30"])
    ax.xaxis.set_minor_locator(matplotlib.ticker.NullLocator())
    ax.set_yticks([0.5, 1, 2, 4])
    ax.set_yticklabels(["0.5", "1", "2", "4"])
    ax.yaxis.set_minor_locator(matplotlib.ticker.NullLocator())
    ax.text(38, 1.07, "median wins", fontsize=6.4, color=INK2, va="bottom", ha="right")
    ax.text(38, 0.93, "mean wins", fontsize=6.4, color=INK2, va="top", ha="right")
    ax.set_xlabel(r"degrees of freedom $\nu$")
    ax.set_ylabel("ARE (median vs. mean)")
    ax.legend(loc="upper right", handlelength=1.4, labelspacing=0.25)
    _panel(ax, "b", r"Student-$t_\nu$: efficiency")

    # (c) Cauchy running estimates ------------------------------------------------------------
    ax = axs[2]
    ax.axhline(0, color=BASE, lw=0.9, zorder=1)
    for est in ("mean", "median"):
        col = {"mean": "running_mean", "median": "running_median"}[est]
        for p, s in paths.groupby("path"):
            ax.plot(s["n"], s[col], color=COL[est], lw=0.95, alpha=0.95 if est == "median" else 0.85,
                    zorder=3 if est == "median" else 2)
    ax.set_xscale("log")
    ax.set_ylim(-9, 9)
    ax.set_xlim(1, paths["n"].max() * 1.05)
    ax.set_xlabel("number of observations")
    ax.set_ylabel("estimate (true centre $=0$)")
    ax.legend(handles=[Line2D([], [], color=COL["mean"], lw=1.5, label="running mean"),
                       Line2D([], [], color=COL["median"], lw=1.5, label="running median")],
              loc="lower right", handlelength=1.5, labelspacing=0.25, frameon=True, facecolor="white",
              edgecolor="none", framealpha=0.9)
    _panel(ax, "c", "Cauchy: 4 sample paths")
    _save(fig, "fig_core")


# =============================================================================================
def figure_robust():
    """(d,e) MSE vs contamination, (f) sensitivity curves, (g) breakdown."""
    comp = pd.read_csv(C.DATA_DIR / "comparators.csv")
    sens = pd.read_csv(C.DATA_DIR / "sensitivity.csv")
    brk = pd.read_csv(C.DATA_DIR / "breakdown.csv")
    n = C.N_SWEEP_CONTAM
    fig, axs = plt.subplots(2, 2, figsize=(FULL_W, 3.8), gridspec_kw={"wspace": 0.24, "hspace": 0.56, "top": 0.915})

    # shared legend row on top
    handles = [Line2D([], [], color=COL[e], marker=MARK[e], lw=1.5, label=l) for e, l in
               (("mean", "mean"), ("median", "median"), ("trim10", "10% trimmed mean"),
                ("huber", "Huber ($k{=}1.345$)"), ("hl", "Hodges–Lehmann"))]
    handles.append(Line2D([], [], color=INK2, lw=1.0, ls=(0, (1, 1.6)), label="first-order theory (mean, median)"))
    fig.legend(handles=handles, loc="upper center", ncol=6, bbox_to_anchor=(0.5, 1.0), columnspacing=1.2,
               handletextpad=0.5)

    def mse_panel(ax, kind, ylog, ylabel):
        for est in ALL:
            s = comp[(comp.kind == kind) & (comp.n == n) & (comp.estimator == est)].sort_values("eps")
            y = s["mse"] if ylog else s["n_mse"]
            ax.plot(s["eps"], y, color=COL[est], marker=MARK[est], zorder=3 if est in ("mean", "median") else 2)
        eps = np.linspace(0, 0.2, 81)
        for est in ("mean", "median"):
            ds = [Dist(kind, eps=e, scale=C.TUKEY_SCALE, shift=C.SHIFT_LOC) for e in eps]
            th = np.array([T.asymptotic_mse(d, est, n) for d in ds])
            ax.plot(eps, th if ylog else th * n, color=COL[est], lw=1.0, ls=(0, (1, 1.6)), zorder=4)
        if ylog:
            ax.set_yscale("log")
        ax.set_xlim(-0.005, 0.205)
        ax.set_xlabel(r"contamination fraction $\varepsilon$")
        ax.set_ylabel(ylabel)

    mse_panel(axs[0, 0], "shift", True, r"MSE about $0$ ($n=100$)")
    _panel(axs[0, 0], "d", r"Outliers: $(1\!-\!\varepsilon)N(0,1)+\varepsilon N(10,1)$")
    mse_panel(axs[0, 1], "tukey", False, r"$n\cdot$MSE ($n=100$)")
    _panel(axs[0, 1], "e", r"Tukey: $(1\!-\!\varepsilon)N(0,1)+\varepsilon N(0,3^2)$")

    # (f) sensitivity curves ---------------------------------------------------------------------
    ax = axs[1, 0]
    ax.axhline(0, color=BASE, lw=0.9, zorder=1)
    lw = {"trim10": 3.0, "huber": 1.2}                 # trimmed mean and Huber nearly coincide: draw one wide, one thin
    for est in ALL:
        ax.plot(sens["y"], sens[est], color=COL[est], lw=lw.get(est, 1.4), zorder=3 if est == "mean" else 2)
    ax.set_xlim(-C.SC_YMAX, C.SC_YMAX)
    ax.set_ylim(-4.2, 4.2)
    ax.set_xlabel(r"location of the added observation $y$")
    ax.set_ylabel(r"$\mathrm{SC}_{100}(y)$")
    ax.text(-7.7, 3.6, "mean: unbounded", color=INK2, fontsize=6.4, ha="left", va="top")
    ax.text(5.2, -2.7, "all others: bounded", color=INK2, fontsize=6.4, ha="center", va="center")
    _panel(ax, "f", "Sensitivity curves ($n=100$)")

    # (g) breakdown -----------------------------------------------------------------------------
    ax = axs[1, 1]
    for est in ALL:
        s = brk[brk.estimator == est].sort_values("k")
        ax.plot(s["frac"], np.maximum(s["median_abs"], 1e-2), color=COL[est], lw=1.4, marker=MARK[est],
                markevery=(0, 6), zorder=3 if est in ("mean", "median") else 2)
    for x, lab in ((0.0, "mean"), (0.10, "trim"), (0.293, "HL"), (0.5, "median,\nHuber")):
        ax.axvline(x, color=MUTED, lw=0.7, ls=(0, (1, 2)), zorder=1)
        ax.text(x + 0.008, 6e7, lab, fontsize=6.2, color=INK2, va="top", ha="left", linespacing=1.0)
    ax.set_yscale("log")
    ax.set_ylim(1e-2, 1e8)
    ax.set_xlim(-0.01, 0.6)
    ax.set_yticks([1e-2, 1, 1e2, 1e4, 1e6])
    ax.set_yticklabels([r"$10^{-2}$", "1", r"$10^{2}$", r"$10^{4}$", r"$10^{6}$"])
    ax.yaxis.set_minor_locator(matplotlib.ticker.NullLocator())
    ax.set_xlabel(r"fraction of observations replaced by $10^6$")
    ax.set_ylabel(r"median $|$estimate$|$ ($n=101$)")
    _panel(ax, "g", "Breakdown experiment")
    _save(fig, "fig_robust")


def make_all():
    style()
    figure_core()
    figure_robust()
