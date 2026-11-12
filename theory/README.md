# Theory

**Owner:** _to be filled in by the group_  **Status:** section text not yet added to the repository

## Scope (from the project proposal)

* Normal model: $\operatorname{Var}(\bar X_n)=\sigma^2/n$, Fisher information $I(\mu)=1/\sigma^2$, so $\bar X_n$ attains the Cramér–Rao lower bound.
* Sample median: exact law through order statistics (Beta transform of $U_{(k+1)}$) and the limit
  $\sqrt n\,(M_n-m)\xrightarrow{d}N\big(0,[4f(m)^2]^{-1}\big)$, hence $\operatorname{ARE}(M_n,\bar X_n)=2/\pi\approx0.637$ at the normal.
* General location families: $\operatorname{ARE}=4\sigma^2 f(m)^2$; Student-$t_\nu$ (Cramér–Rao bound $\tfrac{\nu+3}{\nu+1}\cdot\tfrac1n$, crossover at $\nu^*\approx4.68$), Laplace (ARE $=2$), Cauchy
  ($\bar X_n\overset{d}{=}X_1$; median variance $\approx\pi^2/(4n)$ against bound $2/n$, efficiency $8/\pi^2\approx0.81$).

## Numerical cross-checks that already exist

Every formula above is implemented in [`simulation/mvm/theory.py`](../simulation/mvm/theory.py) and tested against numerical integration
and Monte Carlo in [`simulation/tests/test_theory.py`](../simulation/tests/test_theory.py) and
[`simulation/tests/test_simulation.py`](../simulation/tests/test_simulation.py). Cite them as verification of the derivations.

## Where to put the section

Add the derivations here (e.g. `theory/theory_section.tex`) and link them from the root `README.md` and from `docs/index.html`
(the Theory block of `tools/site_template.html` currently holds a short summary of the same results).
