# Application: predicting the winner on class data

**Owner:** _to be filled in by the group_  **Status:** analysis and results not yet added to the repository

## Plan (from the project proposal)

1. Collect **anonymous** class data on variables of different shape (e.g. height, hours of sleep, Instagram followers). Record no names; participation is optional.
2. For each variable compute the plug-in estimate
   $$\widehat{\operatorname{ARE}}(M_n,\bar X_n)\approx 4\,\hat\sigma^2\,\hat f(\hat m)^2,$$
   with $\hat\sigma^2$ the sample variance, $\hat m$ the sample median and $\hat f$ a Gaussian kernel density estimate.
   $\widehat{\operatorname{ARE}}<1$ predicts the mean is more precise, $>1$ the median.
3. Check the prediction with $10^4$ bootstrap resamples (smoothed bootstrap for rounded variables such as height): compare
   $\operatorname{Var}^*(\bar X)/\operatorname{Var}^*(M)$ with $\widehat{\operatorname{ARE}}$ and report percentile intervals.

## Things to keep in mind (found while preparing the proposal)

* A synthetic check with $n=40$ showed the plug-in estimate is **noisy and biased low** under normality (mean 0.55 against the true 0.637, s.d. 0.10; for lognormal data s.d. about 0.35),
  so report intervals and treat only values clearly away from 1 as a verdict.
* For **skewed** variables the mean and the median estimate different parameters; the ARE then compares precision only, and the write-up should say so.
* The simulation component supplies reusable pieces: Gaussian-kernel and bootstrap code can follow the structure of `simulation/mvm/engine.py`,
  and `simulation/mvm/theory.py::are_median_vs_mean` gives the theoretical ARE for the benchmark distributions.

## Where to put the work

`application/analyze.py` (analysis), `application/data/` (anonymous CSV), `application/README.md` (results). Then link it from the root `README.md` and `docs/index.html`.
