## Theoretical Background

### 1. Basic Concepts and Properties of Estimators

**Estimator.** An estimator is a statistic calculated from a random sample and used to estimate an unknown population parameter. Let $X_1,\ldots,X_n$ be an i.i.d. sample from a population with parameter $\theta$. An estimator of $\theta$ is denoted by $T=T(X_1,\ldots,X_n)$.

**Sample Mean and Sample Median.** The sample mean is defined as $\bar X_n=\frac{1}{n}\sum_{i=1}^{n}X_i$. The sample median is the middle value of the observations when they are arranged in non-decreasing order. Let $X_{(1)}\leq X_{(2)}\leq\cdots\leq X_{(n)}$ denote the ordered observations (order statistics). The sample median $M_n$ is defined as $M_n=X_{(k+1)}$ if $n=2k+1$, and $M_n=\frac{X_{(k)}+X_{(k+1)}}{2}$ if $n=2k$. Thus, for an odd sample size, the median is the single middle observation, whereas for an even sample size, it is the arithmetic mean of the two middle observations.

**Bias and Unbiasedness.** The bias of an estimator $T$ of $\theta$ is $\mathrm{Bias}(T)=E(T)-\theta$. An estimator is called unbiased if $E(T)=\theta$. For a symmetric population with centre $\theta$, both the sample mean (when its expectation exists) and the sample median are unbiased under the usual conditions.

**Variance.** The variance measures the variability of an estimator over repeated samples: $\mathrm{Var}(T)=E[(T-E(T))^2]$. A smaller variance indicates greater precision. For an i.i.d. sample with finite population variance $\sigma^2$, $\mathrm{Var}(\bar X_n)=\frac{\sigma^2}{n}$.

**Mean Squared Error (MSE).** The mean squared error measures the expected squared difference between an estimator and the true parameter: $\mathrm{MSE}(T)=E[(T-\theta)^2]$. It can be decomposed as $\mathrm{MSE}(T)=\mathrm{Var}(T)+[\mathrm{Bias}(T)]^2$. For an unbiased estimator, $\mathrm{MSE}(T)=\mathrm{Var}(T)$. Thus, variance and MSE can be used to compare unbiased estimators, while MSE also accounts for bias.

**Consistency.** An estimator $T_n$ is consistent for $\theta$ if it converges in probability to $\theta$ as the sample size increases. This means that, for every $\varepsilon>0$, $P(|T_n-\theta|>\varepsilon)\to 0$. Under suitable conditions, both the sample mean and sample median are consistent estimators of their corresponding population centres.

### 2. Asymptotic Behaviour and Efficiency

Asymptotic theory studies the behaviour of estimators as the sample size $n$ tends to infinity.

**Asymptotic Behaviour of the Sample Mean.** Suppose $X_1,\ldots,X_n$ are i.i.d. with finite mean $\mu$ and finite variance $\sigma^2$.

**Law of Large Numbers (LLN).** The sample mean converges in probability to the population mean: $\bar X_n\to\mu$. Thus, the sample mean is a consistent estimator of $\mu$.

**Central Limit Theorem (CLT).** Under the same finite-variance condition, the standardized sample mean converges in distribution to a standard normal random variable: $\frac{\sqrt n(\bar X_n-\mu)}{\sigma}\to N(0,1)$. Consequently, for large $n$, $\bar X_n$ is approximately distributed as $N\left(\mu,\frac{\sigma^2}{n}\right)$. Its asymptotic variance, in the $\sqrt n$ scaling, is $\sigma^2$.

**Asymptotic Behaviour of the Sample Median.** Let $m$ denote the population median and let $f(x)$ be the population density. If the density is continuous and positive at $m$, then the sample median has the asymptotic distribution $\sqrt n(M_n-m)\to N\left(0,\frac{1}{4f(m)^2}\right)$. Therefore, for large $n$, $M_n$ is approximately distributed as $N\left(m,\frac{1}{4nf(m)^2}\right)$. Its asymptotic variance, in the $\sqrt n$ scaling, is $\mathrm{AVar}(M_n)=\frac{1}{4f(m)^2}$. This result provides the theoretical basis for comparing the asymptotic efficiency of the sample mean and median.

**Relative Efficiency.** For two unbiased estimators $T_1$ and $T_2$ of the same parameter, the relative efficiency of $T_1$ with respect to $T_2$ is defined as $\mathrm{RE}(T_1,T_2)=\frac{\mathrm{Var}(T_2)}{\mathrm{Var}(T_1)}$. A value greater than one indicates that $T_1$ has smaller variance and is more efficient.

**Asymptotic Relative Efficiency (ARE).** The asymptotic relative efficiency compares the asymptotic variances of two estimators: $\mathrm{ARE}(T_1,T_2)=\frac{\mathrm{AVar}(T_2)}{\mathrm{AVar}(T_1)}$. For the sample median relative to the sample mean, when the population has finite variance and a positive density at its median, $\mathrm{ARE}(M_n,\bar X_n)=4\sigma^2 f(m)^2$. Therefore, $\mathrm{ARE}>1$ indicates that the median is more efficient; $\mathrm{ARE}<1$ indicates that the mean is more efficient; and $\mathrm{ARE}=1$ indicates that both have equal asymptotic efficiency.

**Normal Distribution: Reference Model.** Let $X\sim N(\mu,\sigma^2)$. Its probability density function is $f(x)=\frac{1}{\sigma\sqrt{2\pi}}\exp\left[-\frac{(x-\mu)^2}{2\sigma^2}\right]$. For this distribution, mean = median = mode = $\mu$. The sample mean is the maximum likelihood estimator (MLE) of $\mu$ when $\sigma^2$ is known. The Fisher information per observation for $\mu$ is $I(\mu)=\frac{1}{\sigma^2}$. Hence, the Cramér-Rao lower bound for an unbiased estimator of $\mu$ based on $n$ observations is $\mathrm{CRLB}(\mu)=\frac{1}{nI(\mu)}=\frac{\sigma^2}{n}$. Since $\mathrm{Var}(\bar X_n)=\frac{\sigma^2}{n}$, the sample mean attains the CRLB and is efficient.

At the median $m=\mu$, $f(\mu)=\frac{1}{\sigma\sqrt{2\pi}}$. Thus, $\mathrm{ARE}(M_n,\bar X_n)=4\sigma^2\left(\frac{1}{\sigma\sqrt{2\pi}}\right)^2=\frac{2}{\pi}\approx0.637$. Therefore, under normality, the sample mean is asymptotically more efficient than the sample median.

### 3. Comparison Under Heavy-Tailed Distributions

**Laplace Distribution.** The Laplace distribution, also called the double exponential distribution, is symmetric around its location parameter and has heavier tails than the normal distribution. If $X\sim\mathrm{Laplace}(\mu,b)$, its density is $f(x)=\frac{1}{2b}\exp\left(-\frac{|x-\mu|}{b}\right)$, where $b>0$. Its important properties are $E(X)=\mu$, $\mathrm{Median}(X)=\mu$, and $\mathrm{Var}(X)=2b^2$. The density at the median is $f(\mu)=\frac{1}{2b}$. The asymptotic variance of the sample mean is $\mathrm{AVar}(\bar X_n)=2b^2$, while that of the sample median is $\mathrm{AVar}(M_n)=\frac{1}{4f(\mu)^2}=b^2$. Therefore, $\mathrm{ARE}(M_n,\bar X_n)=\frac{2b^2}{b^2}=2$. Thus, the sample median is asymptotically twice as efficient as the sample mean under the Laplace model.

*Note:* The Laplace density has a cusp at its location parameter. The usual regularity conditions for the classical Cramér-Rao inequality do not hold in their standard form for this location model.

**Student's t-Distribution.** The Student's t-distribution is symmetric and has heavier tails than the normal distribution. Its tail heaviness depends on the degrees of freedom $\nu$. For $X\sim t_\nu$, the density is $f(x)=\frac{\Gamma\left(\frac{\nu+1}{2}\right)}{\sqrt{\nu\pi}\,\Gamma\left(\frac{\nu}{2}\right)}\left(1+\frac{x^2}{\nu}\right)^{-\frac{\nu+1}{2}}$, where $\nu>0$ is the degrees of freedom and $\Gamma$ is the Gamma function. Its median is zero. Its mean and variance are $E(X)=0$ for $\nu>1$, and $\mathrm{Var}(X)=\frac{\nu}{\nu-2}$ for $\nu>2$. The mean does not exist for $0<\nu\leq1$, and the variance is infinite for $1<\nu\leq2$.

The density at the median is $f_\nu(0)=\frac{\Gamma\left(\frac{\nu+1}{2}\right)}{\sqrt{\nu\pi}\,\Gamma\left(\frac{\nu}{2}\right)}$. For $\nu>2$, the asymptotic variance of the sample mean is $\mathrm{AVar}(\bar X_n)=\frac{\nu}{\nu-2}$. The asymptotic variance of the sample median is $\mathrm{AVar}(M_n)=\frac{1}{4f_\nu(0)^2}$. Hence, $\mathrm{ARE}(M_n,\bar X_n)=\frac{4\nu f_\nu(0)^2}{\nu-2}$.

The comparison depends on the degrees of freedom: for approximately $2<\nu<4.679$, the median is asymptotically more efficient, whereas for $\nu>4.679$, the mean is asymptotically more efficient. Thus, the relative efficiency changes with tail heaviness. As $\nu$ increases, the Student's t-distribution approaches the standard normal distribution.

### 4. Cauchy Distribution and Overall Comparison

**Cauchy Distribution.** The Cauchy distribution is a symmetric continuous distribution with extremely heavy tails. Its mean and variance do not exist. If $X\sim\mathrm{Cauchy}(x_0,\gamma)$, its density is $f(x)=\frac{1}{\pi\gamma}\left[1+\left(\frac{x-x_0}{\gamma}\right)^2\right]^{-1}$, where $\gamma>0$. Its median and mode are both equal to $x_0$, while $E(X)$ does not exist and $\mathrm{Var}(X)$ does not exist. Therefore, the sample mean does not estimate a finite population mean in this model.

**Asymptotic Behaviour of the Sample Median.** The density at the population median $x_0$ is $f(x_0)=\frac{1}{\pi\gamma}$. Using the asymptotic distribution of the sample median, $\mathrm{AVar}(M_n)=\frac{1}{4f(x_0)^2}$. Consequently, $\mathrm{AVar}(M_n)=\frac{\pi^2\gamma^2}{4}$. The sample median is consistent for $x_0$ and has a finite asymptotic variance. However, the sample mean has no finite variance. Therefore, the usual variance-based ARE between the sample mean and sample median is undefined. The Cauchy distribution illustrates a situation in which the median provides a meaningful estimator of the population centre, whereas the mean does not.

**Overall Theoretical Comparison.** The following table summarises the asymptotic comparison of the sample median with respect to the sample mean.

| Distribution | ARE of Median | More Efficient |
|:---|:---:|:---|
| Normal | $\frac{2}{\pi}\approx0.637$ | Mean |
| Laplace | $2$ | Median |
| Student's $t_\nu$ | $\frac{4\nu f_\nu(0)^2}{\nu-2}$ | Depends on $\nu$ |
| Cauchy | Not defined | Median suitable for location |

**Conclusion.** The theoretical comparison shows that the relative efficiency of the sample mean and sample median depends on the underlying population distribution. Under normality, the sample mean is more efficient and attains the Cramér-Rao lower bound. Under the Laplace distribution, the sample median is twice as efficient asymptotically as the mean. For the Student's t-distribution, their relative efficiency depends on the degrees of freedom. Under the Cauchy distribution, the population mean and variance do not exist, making the usual variance-based comparison impossible. Therefore, neither estimator is universally superior; their suitability depends on the distribution and its properties.
