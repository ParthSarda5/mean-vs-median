#!/usr/bin/env python3
"""Run the whole simulation study:  python run_all.py [--quick] [--only main tsweep ...] [--workers 2]

Writes CSVs to results/data, figures to results/figures and LaTeX tables to ../report/tables
(--quick writes to results_quick/ instead and leaves the real results untouched).
"""
import argparse
import sys
import time

from mvm import config as C
from mvm.experiments import ALL_EXPERIMENTS

NEEDS_WORKERS = {"main", "tsweep", "comparators"}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--quick", action="store_true", help="smoke test: ~20x fewer replications")
    ap.add_argument("--only", nargs="*", choices=list(ALL_EXPERIMENTS), help="run only these experiments")
    ap.add_argument("--workers", type=int, default=2, help="processes for the Monte Carlo grids")
    ap.add_argument("--no-figures", action="store_true", help="skip figures and tables")
    args = ap.parse_args()

    if args.quick:                      # smoke test: never overwrite the real results
        qd = C.ROOT / "results_quick"
        C.DATA_DIR, C.FIG_DIR, C.TABLE_DIR = qd / "data", qd / "figures", qd / "tables"
    C.DATA_DIR.mkdir(parents=True, exist_ok=True)
    C.FIG_DIR.mkdir(parents=True, exist_ok=True)
    C.TABLE_DIR.mkdir(parents=True, exist_ok=True)
    for name in args.only or ALL_EXPERIMENTS:
        t0 = time.time()
        fn = ALL_EXPERIMENTS[name]
        kw = {"quick": args.quick}
        if name in NEEDS_WORKERS:
            kw["workers"] = args.workers
        fn(**kw)
        print(f"== {name} done in {time.time() - t0:.0f}s", file=sys.stderr, flush=True)

    if not args.no_figures:
        from mvm import plotting, tables
        plotting.make_all()
        tables.make_all()


if __name__ == "__main__":
    main()
