
## Robustness of the Sample Mean and Sample Median

Robustness refers to the ability of an estimator to remain stable when a small proportion of observations are contaminated or contain extreme values. The sample mean is sensitive to extreme observations because it uses the magnitude of every data point. The sample median, on the other hand, depends on the middle observation(s), making it less affected by extreme values.

To study this behaviour, we conducted three types of robustness analyses: contamination experiments, sensitivity curve analysis, and a breakdown experiment. We also compared the mean and median with three alternative estimators: the 10% trimmed mean, Huber M-estimator, and Hodges–Lehmann estimator.

### Effect of Contamination

We considered two types of contamination to understand how different forms of outliers affect the estimators. All comparisons were made using the mean squared error (MSE), with the true centre fixed at $\theta=0$.

**Asymmetric outlier contamination**

We considered the model

$$
X \sim (1-\varepsilon)N(0,1)+\varepsilon N(10,1),
$$

where $\varepsilon$ is the proportion of contaminated observations. The main distribution is $N(0,1)$, while the contaminated observations come from a distribution centred at 10.

The target remains the centre of the original distribution, $\theta=0$. Under this model, both the mean and median become biased. However, the mean is affected much more strongly because even a small number of observations near 10 can pull the average away from zero.

For example, at $\varepsilon=0.10$, the scaled MSE values ($n\cdot\mathrm{MSE}$, with $n=100$) are approximately:

| Estimator | $n\cdot\mathrm{MSE}$ |
|---|---:|
| Sample mean | 109.84 |
| Sample median | 3.94 |
| 10% trimmed mean | 15.86 |
| Huber M-estimator | 5.68 |
| Hodges–Lehmann estimator | 6.41 |

The sample mean has an MSE approximately 28 times that of the sample median. The median performs much better because the contaminated observations are concentrated far from the original centre. The Huber estimator also performs well, while the trimmed mean and Hodges–Lehmann estimator provide intermediate levels of protection.

As the contamination fraction increases, the MSE of the mean rises sharply. This demonstrates its sensitivity to observations that are far from the main data distribution.

**Tukey scale contamination**

We also considered symmetric contamination of the form

$$
X \sim (1-\varepsilon)N(0,1)
+\varepsilon N(0,3^2).
$$

Here, the contaminated observations have a larger spread but remain centred at zero. Thus, unlike the previous experiment, the contamination does not shift the centre; it introduces observations with unusually large magnitudes.

At $\varepsilon=0.10$, the approximate scaled MSE values are:

| Estimator | $n\cdot\mathrm{MSE}$ |
|---|---:|
| Sample mean | 1.80 |
| Sample median | 1.81 |
| 10% trimmed mean | 1.31 |
| Huber M-estimator | 1.31 |
| Hodges–Lehmann estimator | 1.32 |

In this case, the mean and median have almost identical MSE values. However, the three alternative estimators perform better. The trimmed mean reduces the effect of large observations by removing a fixed proportion from both tails. The Huber and Hodges–Lehmann estimators also reduce the influence of extreme observations.

**Observation:** The two contamination experiments show that robustness depends on the type of contamination. When the outliers shift the centre, the median performs much better than the mean. When the contamination only increases the spread symmetrically, the mean and median perform similarly, while the three compromise estimators achieve lower MSE.

### Sensitivity Curves

A sensitivity curve measures how much an estimator changes when one additional observation is introduced into an otherwise fixed dataset.

For this experiment, we considered $n=100$ observations generated from standard normal quantiles. We then added one observation at a time, varying its value $y$ over the interval $[-8,8]$.

The sensitivity curve is defined as

$$
SC_{100}(y)=100\left[T_{101}(x_1,\ldots,x_{100},y)-T_{100}(x_1,\ldots,x_{100})\right],
$$

where $T_n$ represents the estimator based on $n$ observations.

The sensitivity curves show the following behaviour:

- **Sample mean:** Its sensitivity curve is a straight line that continues to increase as the added observation becomes more extreme. Therefore, its sensitivity is unbounded.
- **Sample median:** Its sensitivity remains bounded. Once the added observation is sufficiently far from the centre, increasing its magnitude further has little effect on the median.
- **10% trimmed mean:** Its sensitivity is bounded because extreme observations are removed from the tails before calculating the mean.
- **Huber M-estimator:** Its sensitivity is bounded because the estimator reduces the influence of observations with large deviations from the centre.
- **Hodges–Lehmann estimator:** Its sensitivity is also bounded in this experiment.

Thus, the sensitivity curves illustrate why the mean is more vulnerable to extreme observations, whereas the other four estimators show bounded sensitivity.

### Breakdown Experiment

The breakdown point measures the smallest proportion of observations that must be replaced by extreme values before an estimator can be made arbitrarily large.

In our experiment, we generated samples of size $n=101$ from $N(0,1)$. We then replaced the first $k$ observations by $10^6$, increasing $k$ from 0 to 60.

The breakdown condition was defined as the smallest $k$ for which the median of the absolute estimated values over 1,000 replications exceeded $10^3$.

The results are:

| Estimator | Number of replacements required ($k$) | Approximate breakdown fraction |
|---|---:|---:|
| Sample mean | 1 | $1/101 \approx 0.99\%$ |
| 10% trimmed mean | 11 | $11/101 \approx 10.89\%$ |
| Hodges–Lehmann | 30 | $30/101 \approx 29.70\%$ |
| Sample median | 51 | $51/101 \approx 50.50\%$ |
| Huber M-estimator | 51 | $51/101 \approx 50.50\%$ |

The sample mean breaks down after replacing just one observation. This is because a single sufficiently large value can make the mean arbitrarily large.

The sample median can withstand almost half the observations being replaced by extreme values. It breaks down only when at least 51 of the 101 observations are replaced.

The trimmed mean and Hodges–Lehmann estimator have intermediate breakdown behaviour. The Huber estimator, with the fixed MAD scale used in this experiment, has a breakdown point comparable to that of the median.

**Observation:** The breakdown experiment confirms that the sample median and Huber estimator are much more resistant to extreme contamination than the sample mean.

### Overall Findings on Robustness

The three experiments provide different ways of understanding robustness:

1. **Contamination experiments** show how the estimators behave when observations are drawn from distributions containing outliers. The mean performs poorly under asymmetric contamination, while the median is considerably more stable.
2. **Sensitivity curves** show how much an estimator changes when a single extreme observation is added. The mean has unbounded sensitivity, whereas the other four estimators have bounded sensitivity in the experiment.
3. **Breakdown experiment** shows how many observations must be replaced by extreme values to make an estimator take an extremely large value. The mean breaks down with only one replacement, while the median and Huber estimator withstand approximately half the sample being replaced.

Overall, the sample mean is highly sensitive to extreme observations, whereas the sample median is considerably more robust. The trimmed mean, Huber M-estimator, and Hodges–Lehmann estimator provide alternative compromises between efficiency and robustness. The choice of estimator therefore depends not only on the underlying distribution but also on the possibility of contamination in the data.
