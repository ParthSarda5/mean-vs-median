"""Central configuration: every number that defines the study lives here."""
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]          # .../simulation
DATA_DIR = ROOT / "results" / "data"
FIG_DIR = ROOT / "results" / "figures"
TABLE_DIR = ROOT.parent / "report" / "tables"

for _d in (DATA_DIR, FIG_DIR, TABLE_DIR):                 # make sure output folders exist
    _d.mkdir(parents=True, exist_ok=True)

SEED = 20261110                                     # master seed; every cell derives its own stream from it

# --- main grid (mean & median only) ---------------------------------------------------------
N_GRID = (10, 25, 50, 100, 500, 1000)
B_MAIN = 100_000                                    # replications per (distribution, n) cell
T_NUS = (1, 2, 3, 4, 5, 10, 30)                     # Student-t degrees of freedom (t_1 = Cauchy)
EPS_GRID = (0.01, 0.05, 0.10, 0.20)                 # contamination fractions (eps = 0 is N(0,1))
TUKEY_SCALE = 3.0                                   # contaminating component N(0, 3^2)
SHIFT_LOC = 10.0                                    # contaminating component N(10, 1)

# --- t sweep for ARE(nu) ---------------------------------------------------------------------
T_SWEEP_NUS = (3, 3.5, 4, 4.5, 5, 6, 8, 10, 15, 30)
T_SWEEP_N = (100, 1000)
B_SWEEP = 50_000

# --- comparators (mean, median, 10% trimmed mean, Huber, Hodges-Lehmann) ---------------------
N_COMP = (25, 100)
B_COMP = 50_000
EPS_SWEEP = (0.0, 0.01, 0.02, 0.05, 0.10, 0.15, 0.20)
N_SWEEP_CONTAM = 100

# --- Cauchy running estimates ----------------------------------------------------------------
CAUCHY_PATH_LEN = 20_000
CAUCHY_N_PATHS = 4

# --- breakdown experiment --------------------------------------------------------------------
BREAKDOWN_N = 101
BREAKDOWN_REPS = 1_000
BREAKDOWN_VALUE = 1e6                               # gross-error value substituted for k observations
BREAKDOWN_THRESHOLD = 1e3                           # "broken down" once the median |estimate| exceeds this

# --- sensitivity curves ----------------------------------------------------------------------
SC_N = 100
SC_YMAX = 8.0
SC_NPOINTS = 801

N_BATCHES = 20                                      # batches per cell, used for Monte Carlo standard errors
QUANTILES = (0.5, 0.9, 0.99)                        # quantiles of |error| reported (needed when moments fail)

QUICK_FACTOR = 20                                   # --quick divides B, reps and path lengths by this
