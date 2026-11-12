# Robustness

**Owner:** _to be filled in by the group_  **Status:** section text not yet added to the repository

## Scope (from the project proposal)

* Influence functions: $\mathrm{IF}_{\bar X}(y)=y-\mu$ (unbounded) versus $\mathrm{IF}_{M}(y)=\operatorname{sign}(y-m)/(2f(m))$ (bounded; $\pm\sqrt{\pi/2}\approx\pm1.25$ at $N(0,1)$).
* Finite-sample breakdown points: $0$ for the mean, $1/2$ for the median.
* Tukey (1960) contamination model $(1-\varepsilon)F+\varepsilon G$; the even/odd-$n$ effect for the median (Hodges–Lehmann, 1967).

## Evidence the simulation component already provides

| Topic | Data | Figure / test |
|---|---|---|
| Finite-sample influence (sensitivity curves) | [`simulation/results/data/sensitivity.csv`](../simulation/results/data/sensitivity.csv) | Fig. 2(f) in the report; `test_sensitivity_curves_approach_influence_functions` |
| Breakdown (k of n = 101 points replaced by $10^6$) | [`simulation/results/data/breakdown.csv`](../simulation/results/data/breakdown.csv) | Fig. 2(g); `test_breakdown_mean_at_one_outlier_median_at_half` |
| MSE under contamination | [`simulation/results/data/comparators.csv`](../simulation/results/data/comparators.csv) | Fig. 2(d,e) |

Observed breakdown (n = 101): mean at k = 1, 10 % trimmed mean at k = 11, Hodges–Lehmann at k = 30 (theory $1-1/\sqrt2\approx29.3\%$), median and Huber (MAD scale) at k = 51 (50.5 %).

## Where to put the section

Add the text here (e.g. `robustness/robustness_section.tex`) and link it from the root `README.md` and from `docs/index.html`.
