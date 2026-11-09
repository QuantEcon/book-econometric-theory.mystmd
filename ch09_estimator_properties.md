---
title: "Properties of Estimators"
---

When the concept of estimators was first introduced in {ref}`s-stats` we mentioned that there is nothing in the definition that guarantees—or even alludes to—good performance. Let's now start to think about evaluating estimators.

## Sampling Distributions

Our first step is to view estimators as random entities dependent on the data, and to consider their distributions under different specifications of the data-generating process. The notation below follows that set out in {ref}`ss-tep`.

### Estimators as Random Elements

By definition, a statistic $T$ is a $\bB$-measurable transformation of the sample space $\Zdata$ into some feature space $S$. As a function of the sample $\zdata  = (\boldz_1, \ldots, \boldz_N)$, it can be viewed as a random element mapping outcomes in the underlying space $\Omega$ to $S$ via

```{math}
:enumerated: false

\omega \mapsto 
    T \left[ \zdata(\omega) \right] \in S
```

Here "random element\" stands for either a random vector (or scalar), or something more general, such as the random distribution $\hat P_N$, or a density estimate, which is a random function of the sample.

```{prf:example}

Consider a sample mean $\bar x_N := \frac{1}{N} \sum_{n=1}^N x_n$. The observations $x_1, \ldots, x_N$ are understood to be random variables on some probability space $(\Omega, \fF, \PP)$. The sample mean inherits this property. Formally, it is the random variable

$$
\bar x_N(\omega) := \frac{1}{N} \sum_{n=1}^N x_n(\omega)
        \qquad (\omega \in \Omega)
$$ (eq-sirv)
```


```{prf:example}

A histogram is a random estimate of the density of the unknown distribution $P$, mapping sample realizations into step functions on the outcome space.
```


The idea that estimators are random can be counterintuitive. When we begin to study econometrics, we tend to think of data as fixed numbers determined by some prior historical outcome and statistics as deterministic functions of these numbers. Hence we only observe one value of any particular statistic—one sample mean, one sample variance, and so on. However, if we want to build a theory of estimators, then we need to consider their performance prior to observing the data. Each data set is viewed as a random draw from some unknown distribution. We seek estimators that perform well in the face of this uncertainty.

Once we agree that estimators are random, we know they must have distributions. These distributions determine how effective an estimator is at estimating the intended feature. Let's now turn to investigating these distributions.

(ss-wasd)=
### Sampling Distributions
Suppose that $T$ is a statistic mapping data $\zdata$ into $S$. As just discussed, $T(\zdata )$ is a random element taking values in $S$. Its distribution $\lL(T)$ on $S$ can be expressed as

$$
\lL(T) (B) 
    = \PP   \left\{ \, T \left(\zdata \right) \in B \, \right\}
$$ (eq-dsampdist)

$\lL(T)$ is called the **sampling distribution** of $T$. Continuing to write $\Pdata$ for the joint distribution of the sample $\zdata$, we can also define the sampling distribution as

```{math}
:enumerated: false

\lL(T) (B) = \Pdata\setntn{\bolds \in \Zdata}{T(\bolds) \in B}
```

or even $\Pdata \circ T^{-1}$ if you like that kind of notation. Note that the sampling distribution of $T$ is fully determined by

(i) the statistic $T$ that maps data to feature space,

(ii) the joint distribution $\Pdata$ of the data.

Since $\Pdata$ is unknown, the sampling distribution is likewise unknown. It can, however, be estimated, as we'll see below.

Mathematically, there is no difference between the "sampling distribution of $T$\" and the "distribution of $T$.\" Nonetheless, the terminology is standard.

```{prf:example}
:label: eg-sdsm

Let $x_1,\ldots,x_N$ be IID on $\RR$ with $P = \lL(x_n) = \nN(\mu, \sigma^2)$ for all $n$. Let $T$ be the sample mean $\bar x_N$. By independence, the joint distribution $\Pdata$ of $\zdata = \boldx := (x_1,\ldots,x_N)$ is $\nN(\mu \boldone, \sigma^2 \boldI)$. Applying [fact %s](#fa-nipre) gives

$$
\lL(\bar x_N) 
        = \nN \left(\mu, \, \frac{\sigma^2}{N} \right)
$$ (eq-sdsmnc)
This is the sampling distribution of $\bar x_N$.
```


Sampling distributions are often complicated. Sometimes the best we can do is obtain the sampling distribution of some transformation of an estimator.

```{prf:example}
:label: eg-sdq

Let $\boldx := (x_1,\ldots,x_N)$ have the same distribution as [example %s](#eg-sdsm). Let $s_N^2$ be the sample variance. If $\sigma > 0$, then

$$
\lL(q_N) = \chi^2(N - 1)
        \quad \text{where} \quad
        q_N := \frac{N \, s^2_N}{\sigma^2}
$$ (eq-sdsmnc2)
```


To see why {eq}`eq-sdsmnc2` holds, recall from [example %s](#eg-pvonesr) that if $S := \Span\{\boldone\}$ and $\boldP = \proj S$, then the residual projection satisfies $\boldM \boldx = \boldx - \bar x \boldone$. Hence

```{math}
:enumerated: false

N s_N^2 = \sum_{n=1}^N (x_n - \bar x_N)^2 = \| \boldM \boldx \|^2
```

Since $\boldM \boldone = \boldzero$ (see [fact %s](#fa-opt3)) and $\boldM$ is idempotent and symmetric ([fact %s](#fa-pmsi)), setting $\boldxi := \sigma^{-1} (\boldx - \mu \boldone)$ gives

$$
q_N
    = \sigma^{-2} \| \boldM \boldx \|^2
    = \left\| \boldM \sigma^{-1} (\boldx - \mu \boldone) \right\|^2
    =  \| \boldM \boldxi \|^2
    = \boldxi^\T \boldM \boldxi
$$ (eq-sn2m)

To see that the right-hand side has distribution $\chi^2(N - 1)$, apply [fact %s](#fa-trop) to obtain $\trace \boldM = N-1$, followed by [fact %s](#fa-mcfi).

```{prf:property}
:label: fa-ismsv

For an IID sample $x_1, \ldots, x_N$ with common distribution $\nN(\mu, \sigma^2)$, the sample mean and sample variance are independent random variables.
```


Many long and involved proofs of [fact %s](#fa-ismsv) can be found in the literature, but a simple and informative proof can be obtained from [fact %s](#fa-poci). See [exercise %s](#ex-ismsv) and its solution.

In the preceding examples we were able to obtain an analytical expression for the sampling distribution of our estimator (or a simple transformation of the estimator) from the underlying distribution of the data. In reality this case is the exception.

```{prf:example}
:label: eg-lnmsm

A random variable $x$ is said to be **lognormally distributed** with parameters $\mu$, $\sigma$ if $\lL(\ln x) = \nN(\mu, \sigma^2)$. We write $\lL(x) = \lN(\mu, \sigma^2)$. No neat expression exists for $\lL(\bar x_N)$ when the data are IID lognormal.
```


When no closed-form exists, we can still approximate the sampling distribution using simulation, conditional on a specification of the joint distribution of the data. For example, [figure %s](#f-lognorm_sample_mean) shows a histogram approximation to $\lL(\bar x_N)$ when $N=20$ and $x_1, \ldots, x_N$ is IID $\lN(0, 1)$. The histogram is produced by generating 10,000 observations of $\bar x_N$. The loop that produced the observations is given in [listing %s](#list-gsm).

```{figure} figures/lognorm_sample_mean.svg
:name: f-lognorm_sample_mean

Sampling distribution of $\bar x_N$ when $N=20$ and data are IID $\lN(0, 1)$
```

```{code-block} python
:caption: Generates 10,000 observations of the sample mean (Python)
:name: list-gsm
import numpy as np

N = 20
num_reps = 10000
xbar_outcomes = np.empty(num_reps)  # Allocate memory

for i in range(num_reps):
    x = np.exp(np.random.randn(N))  # Generate N iid lognormal RVs
    xbar_outcomes[i] = x.mean()
```

(sss-dp)=
#### Desirable Properties
The quality or performance of a given estimator is characterized by its sampling distribution. Let $\hgamma$ be an estimator of some feature $\gamma$ based on a sample $\zdata$. What would be an ideal sampling distribution for $\hgamma$? The answer is obvious: a point mass that puts probability one on the single point $\gamma$.

This scenario rarely arises because our sample is *random*. So here's a more realistic objective:

> It is desirable that the sampling distribution of an estimator of $\gamma$ concentrate most of its probability mass in a small neighborhood around $\gamma$, regardless of the joint distribution $\Pdata$ of the sample.

This statement considers only what is desirable, not what is possible. For example, if we look back at [figure %s](#f-lognorm_sample_mean), we see that the sampling distribution of $\bar x_N$ does not cluster tightly around the mean. (The feature $\gamma$ that we are estimating is the mean of the $\lN(0, 1)$ distribution, which is represented by the vertical black line.) But this is not unexpected, since the sample size is only 20. We simply don't have enough information to estimate the feature precisely.

Nonetheless, these ideas give us a way to start comparing different estimators of a given feature from a given amount of data. For example, suppose that we are estimating the mean of a distribution from a random sample of size $20$. Suppose that each is drawn from the $\lN(0, 1)$ distribution, the mean of which is about 1.65. The sampling distributions of three estimators of the mean under the associated joint distribution have been approximated by simulation and represented as box plots in [figure %s](#f-sampling_distributions). The top sampling distribution is that of the mid-range estimator, the second is that of a maximum likelihood estimator to be discussed later, and the third is that of the sample mean. The last two put more probability mass in the region around 1.65. Hence they outperform the mid-range estimator in this setting.

```{figure} figures/sampling_distributions.svg
:name: f-sampling_distributions

Sampling distributions of three estimators of a lognormal mean
```

In what follows we will formalize the ideas presented above, including both finite sample and asymptotic theory. But first let's discuss a general way of estimating sampling distributions from the data.

### The Bootstrap

While understanding sampling distributions is central to evaluating the precision of estimators, sampling distributions cannot be observed in the wild because they depend on $\Pdata$. For example, suppose $x_1, \ldots, x_N$ are assumed to be IID $\nN(\mu, \sigma^2)$ with $\mu$ and $\sigma$ unknown. By {eq}`eq-sdsmnc`, the distribution of $\bar x_N$ is $\nN(\mu, \sigma^2/N)$. While this provides valuable information, it doesn't fully pin down the distribution because $\mu$ and $\sigma$ are unknown.

Things are even more complicated for the sample mean of IID lognormals in [example %s](#eg-lnmsm). There we can't easily write down the sampling distribution as a function of the parameters. Problems multiply further if we are not prepared to assume a parametric form for the distribution of the observations.

There is, however, a very general technique that allows us to approximate sampling distributions given data. We'll discuss the case of scalar data, although the vector case is similar. To begin, let $T$ be any statistic and consider its sampling distribution

```{math}
:enumerated: false

\lL_P(T) 
    := \text{ the distribution of }
    T(x_1, \ldots, x_N) 
    \text{ when } x_1, \ldots, x_N \iidsim P
```

The sampling distribution is unknown to us because $P$ is unknown. Suppose, however, that we have access to a sample $x^o_1, \ldots, x^o_N$ of IID draws from $P$. We then have a good estimator of $P$ given our sample: the empirical distribution $\hat P_N$ generated by $x^o_1, \ldots, x^o_N$. Given this empirical distribution, the **bootstrap distribution** of $T$ is defined as

```{math}
:enumerated: false

\lL_{\hat P_N}(T) 
    := \text{ the distribution of }
    T(x_1, \ldots, x_N) 
    \text{ when } x_1, \ldots, x_N \iidsim \hat P_N
```

Thus the bootstrap distribution is just the plug-in estimator of the sampling distribution.

Is $\lL_{\hat P_N}(T)$ a good approximation to $\lL_P(T)$? By the Glivenko–Cantelli theorem, $\hat P_N$ will be close to $P$ when $N$ is large. If $Q \mapsto \lL_Q(T)$ is suitably continuous, then $\lL_{\hat P_N}(T)$ will also be close to $\lL_P(T)$. For details, see §29.2 of {cite}`dasguptaasymptotic`.

#### Simulation

One remaining problem with $\lL_{\hat P_N}(T)$ is that it's not particularly tractable. We can, however, generate draws from this distribution easily on a computer. Here's an algorithm that returns $M$ draws from $\lL_{\hat P_N}(T)$:

- set $\hat P_N = $ the empirical distribution of observed data $x^o_1, \ldots, x^o_N$
- **for** $m$ in $1, \ldots, M$ **do**
  - draw $x_1^b, \ldots, x_N^b$ independently from $\hat P_N$
  - set $T^b_m = T(x^b_1, \ldots, x_N^b)$
- **end**
- return the sample $T_1^b, \ldots, T_M^b$

Using the definition of the empirical distribution, drawing samples from $\hat P_N$ is easy: we just make repeated draws from the set $x^o_1, \ldots, x^o_N$ with equal probability on each element. {numref}`list-jboot` contains some Julia code implementing the whole algorithm as a function that returns $M$ bootstrap samples. The other arguments to the function are `xo`, the array of observed data, and `stat`, which is a function representing $T$.

```{code-block} julia
:caption: Function to generate bootstrap samples (Julia)
:name: list-jboot
function bootstrap(xo, stat, M)
    N = length(xo)
    T_b = Vector{Float64}(undef, M)
    x_b = Vector{Float64}(undef, N)
    for m in 1:M
        for i in 1:N
            x_b[i] = xo[rand(1:N)]
        end
        T_b[m] = stat(x_b)
    end
    return T_b
end
```

Here's an example function call:

``` julia
julia> bootstrap([1, 2, 3], mean, 3)
3-element Vector{Float64}:
 1.6666666666666667
 1.0
 1.6666666666666667
```

Here is the same routine in Python:

``` python
def bootstrap(xo, stat, M):
    N = len(xo)
    T_b = np.empty(M)
    for m in range(M):
        x_b = np.random.choice(xo, size=N)
        T_b[m] = stat(x_b)
    return T_b
```

{numref}`f-bootstrap_hist` shows a histogram of observations from a bootstrap distribution generated using the Julia code above. The underlying sample $x^o_1, \ldots, x^o_N$ was created using $N=200$ IID draws from the $\lN(0, 1)$ distribution. The statistic used was the sample median. The density shown in the same figure is a numerical approximation of the true sampling distribution of the median over $200$ IID draws from the $\lN(0, 1)$ distribution. The (population) median of the $\lN(0, 1)$ distribution is 1.

```{figure} figures/bootstrap_hist.svg
:name: f-bootstrap_hist

Bootstrap draws and sampling distribution for the median
```

(sss-se)=
#### Standard Errors
As discussed in {ref}`sss-dp`, a good estimator is one that concentrates its probability mass close to the target feature. The extent to which probability mass concentrates is usually measured by the variance or standard deviation of the sampling distribution. The standard deviation is preferable because it comes in units that are easier to interpret.

Because the sampling distribution is typically unknown, so is its standard deviation. Hence we replace the true standard deviation with an estimate. An estimate of the standard deviation of a sampling distribution is usually referred to as the **standard error**. More generally,

```{math}
:enumerated: false

\se(\hat \gamma) 
    := \text{ an estimate of the standard deviation of } \hat \gamma
```

This isn't a proper definition because we haven't specified *which* estimate of the standard deviation we are talking about, but nevertheless the terminology is common. One way to compute a standard error of an estimator is to take the sample standard deviation of the bootstrap draws.

```{prf:example}

{numref}`f-firms_median` shows a histogram of draws from the bootstrap distribution of the median applied to 2,118 observations of US firm sizes by sales, as measured by 2013 Compustat data. Values are in millions of US dollars. The median of the sample itself is 269.9. The sample standard deviation of this particular set of bootstrap draws is 25.4.
```


```{figure} figures/firms_median.svg
:name: f-firms_median

Bootstrapped draws for the median of firm sizes by sales
```

(s-fsp)=
## Evaluating Estimators
As mentioned above, the standard definition of an estimator does not include any restrictions that imply sensible behavior (as opposed to, say, the definition of a probability). We now turn again to the problem of assessing estimators, following our informal discussion in {ref}`sss-dp`. We'll begin with finite sample properties.

### Bias

One common criterion for assessing estimators is bias, or lack thereof. Let $\hgamma$ be an estimator of a given feature $\gamma = \gamma(P)$. The **bias** of $\hgamma$ is defined as

$$
\bias_P(\hgamma, \gamma) = \EEP  \hgamma - \gamma(P)
$$ (eq-bias)

The estimator $\hgamma$ is called **unbiased** for $\gamma$ over the class of distributions $\pP$ if

```{math}
:enumerated: false

\bias_P(\hgamma, \gamma) = 0
    \quad \text{for all } P \in \pP
```

The definition is the same when the feature $\gamma$ and estimator $\hgamma$ are vector-valued.

The notation $\EEP \hgamma$ indicates that expectation is taken under the assumption that the data points are IID with common distribution $P$. In general, the bias depends on $P$. For example, in the IID setting, the mid-range estimator is unbiased as an estimator of the mean of a uniform distribution on $\RR$ with unknown end points, but biased for many other distributions, including the lognormal.

{numref}`f-sampling_distributions` already suggests bias of the mid-range estimator in the lognormal case. Here's a simulation in `R` that provides further confirmation.

``` r
> mr <- function(x) return((min(x) + max(x)) / 2)
> observations <- replicate(5000, mr(rlnorm(20)))
> mean(observations)   # Approximate sample mean of mid-range
[1] 3.800108
> exp(1/2)      # Mean of lognormal
[1] 1.648721
```

```{prf:property}
:label: fa-ubpi

Let $P$ be a distribution on $\RR^K$, let $\{\boldx_n\}$ be a sample with $\lL(\boldx_n) = P$ for all $n$, and let $\hat P_N$ be the empirical distribution. If

(i) $\gamma(P) = \int h(\bolds) P(\diff \bolds)$ for some integrable function $h$, and

(ii) $\hgamma$ is the plug-in estimator $\hgamma = \int h(\bolds) \hat P_N(\diff  \bolds) = \frac{1}{N} \sum_{n=1}^N h(\boldx_n)$,

then $\hgamma$ is unbiased for $\gamma$ over the set of distributions such that $\int h(\bolds) P(\diff \bolds)$ exists.
```


[Fact %s](#fa-ubpi) is immediate from linearity of expectations, since

```{math}
:enumerated: false

\EE \hgamma
    = \EE \left[ \frac{1}{N} \sum_{n=1}^N h(\boldx_n) \right]
    = \frac{1}{N} \sum_{n=1}^N \EE h(\boldx_n) 
    = \int h(\bolds) P(\diff  \bolds) = \gamma(P)
```

Note that independence isn't used here. [Fact %s](#fa-ubpi) depends only on the existence of a common marginal $P$.

```{prf:example}

For any identically distributed sample, the $k$th sample moment is an unbiased estimator of the $k$th moment, whenever the latter exists.
```


```{prf:example}
:label: eg-svb

Let $x_1, \ldots, x_N$ be an IID sample with finite variance $\sigma^2$. The sample variance $s_N^2$ is, in general, a biased estimator of $\sigma^2$, but not by much. In particular,

$$
\EE s^2_N = \sigma^2 \frac{N-1}{N}
$$ (eq-esv)
```


To obtain this expression, we can combine the definition of $q_N$ in {eq}`eq-sdsmnc2` and the expression in {eq}`eq-sn2m` to obtain

```{math}
:enumerated: false

\EE s_N^2  
    = \frac{1}{N} \sigma^2 \EE[\boldxi^\T \boldM \boldxi]
    = \frac{1}{N} \sigma^2 \trace \boldM
```

Here $\boldM$ is the residual projection for the subspace $\Span\{\boldone\}$ and the second equality is due to [fact %s](#fa-eoqf). Now apply [fact %s](#fa-trop).

(ss-voe)=
### Variance
The variance of a given estimator gives a measure of its dispersion. Sensible estimators almost always have the property that the variance goes to zero as $N \to \infty$.

```{prf:example}

If $x_1,\ldots, x_N$ are uncorrelated with common finite variance $\sigma^2$, then

$$
\var \bar x_N
        = \var \left[ \frac{1}{N} \sum_{n=1}^N x_n \right]
        = \frac{\sigma^2}{N}
$$ (eq-varsm)
We saw this already in {ref}`ss-lln`.
```


````{prf:example}
:label: eg-sdq2

Assume the IID Gaussian setting of [example %s](#eg-sdq), where it was shown that the sample variance satisfies $\lL(\sigma^{-2} N \, s^2_N) = \chi^2(N - 1)$. The variance of a $\chi^2(k)$ distribution is $2k$. Using this fact and some algebra gives

```{math}
:enumerated: false

\var  s^2_N  = \frac{2 \sigma^4}{N^2} (N - 1)
```
````


Is small variance always desirable for an estimator? If the estimator is biased, then perhaps not, since probability mass might be concentrated in the wrong place. But for an unbiased estimator, low variance is good. It means that probability mass is concentrated around the feature that we wish to estimate.

Low variance is, of course, relative. For example, is the variance in {eq}`eq-varsm` low or high? One way to approach this kind of question is to take the class of unbiased estimators of a given feature $\gamma$, and find the estimator in the class with the lowest variance. For given $\gamma$ and given data $\boldx_1,\ldots,\boldx_N$, the estimator in the set of unbiased estimators

```{math}
:enumerated: false

U_{\gamma} := \{\text{all statistics } \hgamma
        \text{ with } \EE \hgamma = \gamma\}
```

that has the lowest variance within this class is called the **minimum variance unbiased estimator**.

Of course, it's quite possible that no minimizer exists. Even if it does, it might be hard to determine in practice, or require strong assumptions on the unknown distribution of the data. For these reasons it's common to focus on smaller classes than $U_{\gamma}$. For example, the estimator in the set of *linear* unbiased estimators

```{math}
:enumerated: false

U^{\ell}_{\gamma} := \{\text{all linear statistics } \hgamma 
        \text{ with } \EE \hgamma = \gamma\}
```

with the lowest variance—if it exists—is called the **best linear unbiased estimator**.

````{prf:example}

Let $x_1,\ldots,x_N$ be IID with common distribution $P$, where $P$ has finite mean $\mu \neq 0$ and variance $\sigma^2$. The set of linear estimators of $\mu$ is given by

```{math}
:enumerated: false

\left\{ 
        \text{all statistics of the form } 
        \hat \mu = \sum_{n=1}^N \alpha_n x_n, \text{ where }
        \alpha_n \in \RR \text{ for } n = 1,\ldots,N
        \right\}
```

Hence the set of linear unbiased estimators of $\mu$ is

```{math}
:enumerated: false

U^{\ell}_{\mu} := 
        \left\{ 
        \text{all } 
        \hat \mu = \sum_{n=1}^N \alpha_n x_n \text{ with }
        \alpha_n \in \RR,\;\; n = 1,\ldots,N
        \text{ and } 
        \EE \left[ \sum_{n=1}^N \alpha_n x_n \right] = \mu
        \right\}
```

Using linearity of expectations, this set can be rewritten as

```{math}
:enumerated: false

U^{\ell}_{\mu} := 
        \left\{ 
        \text{all } 
        \hat \mu = \sum_{n=1}^N \alpha_n x_n \text{ with }
        \sum_{n=1}^N \alpha_n = 1
        \right\}
```

By independence and our rules for variance of sums , the variance of an element of this class is given by

```{math}
:enumerated: false

\var \left[ \sum_{n=1}^N \alpha_n x_n \right]
        =   \sum_{n=1}^N \alpha_n^2 \var x_n 
        =   \sigma^2 \sum_{n=1}^N \alpha_n^2
```

To find the best linear unbiased estimator, we need to solve

```{math}
:enumerated: false

\text{minimize } \; \sigma^2 \sum_{n=1}^N \alpha_n^2 \;
        \text{ over all } \alpha_1,\ldots,\alpha_N \text{ with }
        \sum_{n=1}^N \alpha_n = 1
```

To solve this constrained optimization problem, we can use the Lagrangian, setting

```{math}
:enumerated: false

L(\alpha_1,\ldots,\alpha_N ; \lambda) := \sigma^2 \sum_{n=1}^N \alpha_n^2 - 
        \lambda \left[\sum_{n=1}^N \alpha_n - 1 \right]
```

where $\lambda$ is the Lagrange multiplier. Differentiating with respect to $\alpha_n$ and setting the result equal to zero, the minimizer $\alpha_n^*$ satisfies $\alpha_n^* = \lambda  (2 \sigma^2)^{-1}$ for each $n$. In particular, each $\alpha_n^*$ takes the same value, and hence, from the constraint $\sum_n \alpha_n^* = 1$, we have $\alpha_n^* = 1/N$. Using these values, our estimator becomes

```{math}
:enumerated: false

\sum_{n=1}^N \alpha_n^* x_n  
        = \bar x_N
```

Thus, for IID data with finite variance, the sample mean is the best linear unbiased estimator of $\mu$.
````


The theory of minimum variance unbiased estimators and best linear unbiased estimators is traditionally regarded as important. However, there's no convincing reason to restrict ourselves to unbiased estimators. In fact, as we'll see, there are good reasons not to.

(ss-vvb)=
### Variance versus Bias
The mean squared error (MSE) of a given estimator $\hgamma$ of some feature $\gamma \in \RR^K$ is

$$
\mse (\hgamma, \gamma) := \EE \left\{ \, \| \hgamma - \gamma \|^2 \, \right\}
$$ (eq-mse)

As for bias, $\mse(\hgamma, \gamma)$ depends on the joint distribution of the data as well as the specification of $\hgamma$, although we haven't made that dependence explicit in the notation.

Mean squared error is an overall measure of the performance of an estimator that takes into account both bias and variance. In fact in the scalar case we can neatly decompose mean squared error into the sum of variance and squared bias:

$$
\mse (\hgamma, \gamma) = \var \hgamma + \bias(\hgamma, \gamma)^2
$$ (eq-dmse)

This is just {eq}`eq-dmse0` specialized to the current setting.

````{prf:example}
:label: eg-sdq21

Consider again the sample variance $s_N^2$ under the IID Gaussian setting of [example %s](#eg-sdq2). We saw there that $\var  s^2_N = 2 \sigma^4 N^{-2} (N - 1)$. Combining this result with {eq}`eq-dmse` and {eq}`eq-esv`, we have

```{math}
:enumerated: false

\mse \left( s^2_N, \sigma^2 \right) 
        = \frac{2 \sigma^4}{N^2} (N - 1) 
        + \left[ \sigma^2 \frac{N-1}{N} - \sigma^2 \right]^2
        = \frac{\sigma^4}{N^2}  (2N - 1)
```
````


As can be seen from both {eq}`eq-mse` and {eq}`eq-dmse`, when $\mse (\hat \gamma, \gamma)$ is small for a given estimator $\hgamma$, the estimator concentrates probability mass around the target feature $\gamma$. Hence minimizing MSE is a natural and reasonable criterion. In doing so, we are making our subjective loss function explicit, while choosing a loss function that is simple and highly tractable.

When minimizing MSE, we typically find a trade-off between bias and variance: lower variance costs more bias and vice-versa. Moreover solutions usually turn out to be interior. In particular, we can often reduce variance significantly by accepting a small amount of bias.

```{prf:example}
:label: eg-udmmse

Let $\hgamma$ be any unbiased estimator of a feature $\gamma \in \RR$. Let $v := \var \hgamma$. Consider the class of estimators $\setntn{ \lambda \hgamma }{\lambda \in \RR}$. As [exercise %s](#ex-udmmse) asks you to verify, the value of $\lambda$ that minimizes $\mse(\lambda \hgamma, \gamma)$ is

$$
\lambda^* := \frac{\gamma^2}{\gamma^2 + v}
$$ (eq-oeib)
In other words, the unbiased estimator does not minimize MSE unless $v = 0$.
```


For most sensible estimators, $v$ converges to zero as the sample size $\to \infty$. The message is that unbiased estimators only approach optimality when the data size is large. Here "large\" should be understood as relative to model complexity. More complex models need more data to be estimated effectively.

#### Stein's Example

There is a famous example in the theory of statistics due to {cite}`stein1956inadmissibility` that pertains to so-called admissible estimators. Here we'll borrow Stein's problem and proposed solution to illustrate several ideas, including the relationship between bias and mean squared error.

To make things concrete, consider the following application: In the field of public finance, one topic of interest is ability-based taxation, which has optimality properties in terms of welfare in certain models. Ability is not observable. Suppose, in particular, that we have a population of $K$ individuals, where the $k$th individual has idiosyncratic ability $\mu_k \in [0, 1]$. For each agent, we get to observe $N$ noisy signals of ability, each of the form $x_{kn} = \mu_k + \epsilon_{kn}$, $n=1, \ldots, N$. The noise terms $\epsilon_{kn}$ are all independent $\nN(0, \sigma^2)$. Our task is to estimate the parameter vector $\boldmu = (\mu_1, \ldots, \mu_K)$ given the data from the observations $x_{nk}$. To make things simpler we'll assume that only $\boldmu$ is unknown.

A natural estimator of $\boldmu$ is the vector sample mean

```{math}
:enumerated: false

\bar \boldx_N 
    =
    \begin{pmatrix}
        \frac{1}{N} \sum_{n=1}^N x_{1n}
        \\
        \vdots
        \\
        \frac{1}{N} \sum_{n=1}^N x_{Kn}
    \end{pmatrix}
```

This is the estimator produced by applying the sample analogue principle to our problem. It's also the maximum likelihood estimator and best linear unbiased.

What Stein showed is that, whenever $K > 2$, there exists an alternative estimator $\hboldmu_N$ such that

$$
\EE \left\{ \, \| \hboldmu_N - \boldmu \|^2 \,  \right\}
    < \EE \left\{ \, \| \bar \boldx_N - \boldmu \|^2 \,  \right\}
    \quad
    \text{for all } \boldmu \in \RR^K
$$ (eq-mlii)

In other words, $\hboldmu_N$ uniformly dominates $\bar \boldx_N$ in terms of MSE, for all possible values of $\boldmu$! The estimator in question is what is now known as the **James–Stein estimator**. It is defined as

```{math}
:enumerated: false

\hboldmu_N 
    = \kappa \,
    \bar \boldx_N
    \quad \text{where} \quad
    \kappa :=
        1 - \frac{K - 2 }{N} \frac{\sigma^2}{ \| \bar \boldx_N \|^2}
```

This estimator is biased. It achieves lower mean squared error by reducing variance at the expense of a small amount of bias. A relatively simple proof of {eq}`eq-mlii` can be found on p. 34 of {cite}`young2005essentials`.

We can get some of the intuition from simulation. In our simulation we set $K=200$, $N=5$ and $\sigma=1$. Thus we have 5 noisy observations on the ability of each of 200 agents. For each agent we choose their ability $\mu_k$ independently from a uniform distribution on $[0, 1]$ and next generate the disturbances $\epsilon_{kn}$ as independent standard normals. We then compute the estimators $\bar \boldx_N$ and $\hboldmu_N$ and the corresponding error vectors.

```{figure} figures/stein.svg
:name: f-stein

Error vectors for the two estimators
```

{numref}`f-stein` shows a typical realization. The bars in the top panel are the elements of the vector $\bar \boldx_N - \boldmu$, which records the errors for the sample mean estimator $\bar \boldx_N$. Those in the lower panel correspond to the errors of the James–Stein estimator. These are noticeably downward biased and also noticeably smaller than the errors of the sample mean estimator. This is reflected in the SSE values in the top left of each panel, which give the sum of squared errors (i.e., the squared norms of the error vectors). The SSE for the sample mean estimator is nearly twice as large.

(ss-ape)=
### Asymptotic Properties
Let's now turn to asymptotic properties of estimators. Asymptotic theory concerns whether or not the sampling distribution of a given estimator increasingly concentrates around the feature we wish to estimate, the rate at which this occurs, and the shape of the sampling distribution in large samples.

A fundamental concept in asymptotic statistics is consistency. Let $\{\hgamma_N\}$ be a sequence of estimators of a given feature $\gamma$, based on a sample of size $N$. We say that $\hgamma_N$ is **consistent** for $\gamma$ if

$$
\hgamma_N \toprob \gamma
    \quad \text{as} \quad 
    N \to \infty
$$ (eq-consistent)

The definition is the same whether the estimators are vector or scalar-valued. The statement that $\hgamma_N$ is consistent means more correctly that the sequence $\{\hgamma_N\}$ is consistent.

As for bias, whether consistency holds depends not just on $\hat \gamma_N$ but also on the joint distribution of the data. If we say that $\hgamma_N$ is consistent over a class of distributions $\pP$, we mean that {eq}`eq-consistent` holds true for any distribution in $\pP$.

````{prf:example}

The sample mean $\bar x_N$ of any IID sample is consistent for the mean over the class of distributions on $\RR$ with finite first moment. Indeed, if $P$ is such a distribution and $x_1, \ldots, x_N$ are IID draws from $P$, then the law of large numbers gives

```{math}
:enumerated: false

\frac{1}{N} \sum_{n=1}^N x_n \toprob \int s P(\diff s)
         \quad \text{ as } \quad N \to \infty
```
````


More generally, $\hgamma_N = \frac{1}{N} \sum_{n=1}^N h(\boldx_n)$ is consistent for $\gamma = \int h(\bolds) P(\diff \bolds)$ whenever

(i) $\boldx_1, \ldots, \boldx_N$ is an IID sample with $\lL(\boldx_n) = P$ for all $n$ and

(ii) $\int h(\bolds) P(\diff \bolds)$ exists.

(See {eq}`eq-lln0g`.) For example, the $k$th sample moment is consistent for the $k$th moment over the class of distributions with finite $k$th moment.

Here's an example of consistency that requires a bit more work.

```{prf:example}
:label: eg-svic

For any IID sample $x_1, \ldots, x_N$ with finite variance, the sample variance $s^2_N$ is consistent for the variance. This can be established from the expression

$$
s_N^2 = \frac{1}{N} \sum_{n=1}^N (x_n - \mu)^2 - (\bar x_N - \mu)^2
$$ (eq-sn2e)
where $\mu$ is the common mean (see [ex. %s](#ex-sn2e) and its solution). Applying the law of large numbers to {eq}`eq-sn2e`, the first term converges to $\sigma^2$, while $(\bar x_N - \mu) \toprob 0$, and hence $(\bar x_N - \mu)^2 \toprob 0$ by [fact %s](#fa-reconpro). Consistency of $s^2_N$ follows.
```


```{prf:property}
:label: fa-cpc

If $\hgamma_N$ is consistent for $\gamma$ and $g$ is any continuous function, then $g(\hgamma_N)$ is consistent for $g(\gamma)$.
```


[Fact %s](#fa-cpc) is immediate from [fact %s](#fa-reconpro).

```{prf:example}
:label: eg-ssdc

The sample standard deviation $s_N = \sqrt{s^2_N}$ is consistent for the standard deviation whenever $s^2_N$ is consistent for the variance.
```


Reasonable estimators should be consistent, at least for the IID case. In fact, as we saw in [theorem %s](#t-fts), the entire unknown distribution can be recovered asymptotically with repeated independent samples.

(sss-ad)=
#### Asymptotic Distributions
In the IID Gaussian setting of [example %s](#eg-sdsm), the distribution of the sample mean is $\nN(\mu, \sigma^2 / N)$. Beyond the information contained in the mean and variance, the fact that the sampling distribution is normal gives us important additional knowledge. We know that the distribution is symmetric, that it has thin tails, and that it puts around 95% of its mass within two standard deviations of the mean. We can determine probabilities of intervals or the expectation of any function of $\bar x_N$.

The only concern is our assumption that the observations are Gaussian. Economic theory rarely provides any convincing basis for this assumption. Fortunately, we have the central limit theorem to fall back on, which tells us (remarkably!) that the sample mean is always asymptotically Gaussian, provided that the underlying observations have finite second moment. Even better, the same approximation is valid for many estimators beyond the sample mean.

To clarify, let $\{\hgamma_N\}$ be a sequence of estimators for some feature $\gamma$. We say that $\hgamma_N$ is **asymptotically normal** if there exists a positive definite matrix $\Sigma$ such that

$$
\sqrt{N} (\hgamma_N - \gamma) \tod \nN(\boldzero, \Sigma)  
    \quad \text{as} \quad
    N \to \infty
$$ (eq-asn)

When {eq}`eq-asn` holds, $\Sigma$ is called the **asymptotic variance–covariance matrix** of $\hgamma_N$.

````{prf:example}
:label: eg-svcan

Let $x_1, \ldots, x_N$ be IID with common mean $\mu$ and variance $\sigma^2$. If $\EE[x_n^4] < \infty$, then the sample variance $s_N^2$ is asymptotically normal with

$$
\sqrt{N}(s_N^2 - \sigma^2)
        \tod \nN(0, m_4 - \sigma^4)
        \quad \text{where} \quad
        m_4 := \EE[(x_n - \mu)^4]
$$ (eq-avarsn2)
To show this, we can modify {eq}`eq-sn2e` to get

$$
\sqrt{N}(s_N^2 - \sigma^2)
        = \sqrt{N} 
        \left[ 
            \frac{1}{N} \sum_{n=1}^N (x_n - \mu)^2 - \sigma^2
        \right]
        -  \sqrt{N} (\bar x_N - \mu)^2
$$ (eq-sn2es)
The last term on the right-hand side of {eq}`eq-sn2es` converges to zero in probability, since

```{math}
:enumerated: false

\sqrt{N} (\bar x_N - \mu)^2
        = a_N b_N
        \quad \text{where} \quad
        a_N := \sqrt{N} (\bar x_N - \mu), \;
        b_N := \bar x_N - \mu
```

By the CLT, the LLN and [fact %s](#fa-sluti), we have $a_N b_N \toprob 0$.

To complete the proof, let $Y_n := (x_n - \mu)^2$. The first term in {eq}`eq-sn2es` can then be written as $\sqrt{N}(\bar Y_N - \EE [Y_n])$. Applying the CLT, this expression converges to a zero-mean normal with variance

```{math}
:enumerated: false

\var Y_n
        = \EE[ (Y_n - \sigma^2)^2 ]
        = \EE[ (x_n - \mu)^4 - 2(x_n - \mu)^2 \sigma^2 + \sigma^4 ]
```

The claim follows.
````


```{prf:example}
:label: eg-ssdcan

Under the same assumptions as [example %s](#eg-svcan), the sample standard deviation is also asymptotically normal, with

$$
\sqrt{N}(s_N - \sigma)
        \tod \nN\left(0, \frac{m_4 - \sigma^4}{4 \sigma^2} \right)
$$ (eq-avarsn)
{prf:ref}`ex-ssdcan` asks you to confirm this.
```


```{prf:property}

Asymptotic normality implies consistency. In particular, if $\{\hgamma_N\}$ satisfies {eq}`eq-asn`, then $\hgamma_N \toprob \gamma$ as $N \to \infty$ (see [ex. %s](#ex-anormic)).
```


Asymptotic normality also implies a rate at which $\hgamma_N$ converges to $\gamma$. Indeed, if $\hgamma_N$ is asymptotically normal, then $\sqrt{N} ( \hgamma_N - \gamma)$ does not diverge. This means that $\hgamma_N - \gamma$ goes to zero fast enough to offset the diverging term $\sqrt{N}$. To emphasize this point, we say that an asymptotically normal estimator $\hgamma_N$ is **$\sqrt{N}$-consistent** for $\gamma$.

Still more significant is the fact that asymptotic normality gives us an approximation to the entire sampling distribution. This will provides us with a means of performing inference. [Chapter 10](#c-ciht) gives details.

(ss-risk)=
### Decision Theory
We've now reviewed the most common methods of assessing estimators. Next we'll look more deeply at the theory of estimation and decisions under uncertainty. Our focus will be on finite sample properties. This is the most interesting case, since it captures our need to generalize on the basis of limited information.

Despite the efforts of many brilliant statisticians and econometricians, we still lack a convincing general theory of optimality of estimators that rests on purely objective criteria. For example, while the theory of minimum variance unbiased estimators initially seems attractive because it can potentially yield estimators without any subjective judgment on the part of the statistician, there are simple settings where minimum variance unbiased estimators are *inadmissible*, which means that they fail a minimal requirement for what constitutes a sensible estimator (see below).

Without a fully objective theory of "good\" estimators, there are competing notions of what good means, and we are forced to take a stand on what our criteria for assessing estimators are in different situations. This means that the estimation problem is not fully specified until we make our preferences explicit.

For economists it's natural to think of specifying preference in terms of utility or loss. To fit with modern statistical terminology, we will think in terms of loss. Thus a full description of an estimation problem includes specification of the loss incurred when an estimator performs poorly.

Let's frame these ideas in line with the seminal work of {cite}`wald1939contributions` and subsequent research on decision theory. A general decision problem consists of

(i) a sample space $\xX$,

(ii) an action space $\aA$,

(iii) a set $\ddD$ of decision rules, which are maps $d \colon \xX \to \aA$

(iv) a universe $\Theta$ of possible states of the world, with typical element $\theta$,

(v) a set of probability distributions $\{P_\theta\}$ over the sample space $\xX$, and

(vi) a loss function $L \colon \aA \times \Theta \to \RR$, with interpretation

```{math}
:enumerated: false

L(a, \theta) = \text{loss incurred by choosing action } a \text{ when the
                state of the world is } \theta.
```

The state of the world is treated as unknown. The problem is to choose a decision rule $d \in \ddD$ that generates low loss in some sense to be made precise. Central to this discussion is the notion of **risk**, which is defined as

$$
\risk(d, \theta) := \int L[d(x), \theta ] \, P_\theta(\diff x)
$$ (eq-waldrisk)

The integral is over all $x \in \xX$. The notation $\risk$ is used to distinguish this decision-theoretic concept of risk from $R(f)$, the prediction risk of policy $f$ as defined in {eq}`eq-rf`.

````{prf:example}

Let $\pi$ be an inflation rate and let $x$ be a noisy signal of inflation received by firms. Each firm responds by choosing a price $p = p(x)$ for their product. Let $\Pi(p, \pi)$ represent profits from choosing price $p$ when the inflation rate is $\pi$, and let loss be negative profits (see {cite}`diego` for details and motivation). The risk is then

```{math}
:enumerated: false

\risk(p, \theta) 
        = - \int \Pi(p(x), \pi) P_\theta(\diff x, \diff \pi)
```

If the distribution $P_\theta$ is known to the firm (think of $\{P_\theta\}$ as a singleton), then we have a standard problem of choosing price to maximize expected profits under known probabilities for outcomes. If, however, $P_\theta$ is not known, we have a problem of **Knightian uncertainty**. The risk is unknown because we have uncertainty over probabilities.
````


The problem of choosing estimators is a special case of the decision theory described above. Suppose that we observe data $\zdata$ taking values in $\Zdata$ with joint distribution $P_\theta$ indexed by $\theta \in \Theta$. (The set $\Theta$ can be infinite dimensional and $\theta$ is just an index.) We wish to estimate some $S$-valued feature $\gamma_\theta = \gamma(P_\theta)$ using this data. An estimator in this context is a decision rule $\hat \gamma$ mapping $\Zdata$ to the action space $\aA = S$. If we specify the loss $L(\hat \gamma, \theta)$ of choosing $\hat \gamma$ when the state of the world is $\theta$, our risk becomes

$$
\risk(\hat \gamma, \theta) 
    = \int_{\Zdata} L(\hat \gamma(\boldz), \gamma_\theta) P_\theta(\diff \boldz)
$$ (eq-norage)

```{prf:example}
:label: eg-ccdp

Suppose we want to estimate expected returns $\mu_\theta := \int s P_\theta(\diff s)$ to holding an asset based on IID observations $x_1, \ldots, x_N$ from $P_\theta$. With quadratic loss for prediction error and the sample mean as our estimator, the risk is

$$
\risk(\bar x_N, \theta) 
        = \EE_\theta \left[  (\bar x_N - \mu_\theta)^2 \right]
$$ (eq-rsmean)
Since $\bar x_N$ is unbiased for $\mu_\theta$, this evaluates to the variance of $\bar x_N$ under $P_\theta$. By {eq}`eq-varsm`, this is equal to $\sigma_\theta^2/N = \int (s - \mu_\theta)^2 P_\theta (\diff s)/N$.
```


In the preceding example, the risk of the estimator is its MSE. More generally, if $\hat \gamma$ and $\gamma$ in {eq}`eq-norage` are vectors and $L$ measures squared vector deviation, then $\risk(\hat \gamma, \theta) = \mse(\hgamma, \gamma)$. Thus, risk gives us a generalization of the idea of evaluating estimators based on MSE, as described in {ref}`ss-vvb`. Since loss is subjective, there is no reason we should confine ourselves to quadratic loss.

```{prf:example}
:label: eg-hp

As an estimate of central tendency of housing prices in a given location, it's common to use the sample median rather the sample mean. This is because large upper tails distort the price faced by a "typical buyer.\" The choice of median over mean can be understood as minimization of a different loss function when estimating price given pricing data. See [example %s](#eg-mmed).
```


#### Choosing Decision Rules

In the current setup, a natural idea for making optimal decisions—and for choosing optimal estimators, in particular—is to choose decision rules with low risk. The problem with this approach is that the risk depends, in general, on the unknown state of the world $\theta$, and hence the risk-minimizing decision rule also depends on the unknown value $\theta$. In other words, we face Knightian uncertainty. We can't evaluate risk because we don't know the probabilities.

````{prf:example}

Let's return to the setting of [example %s](#eg-ccdp). The risk of the sample mean is equal to $\sigma_\theta^2/N$ for every value of $\theta$. Consider now a new estimator $\hgamma$ that is identically equal to $1$. In other words, we are predicting 100% return on our investment. This is a great estimator when true mean returns are $1$, since

```{math}
:enumerated: false

\mu_\theta = 1
        \quad \implies \quad
        \risk(\hgamma, \theta) 
        = \EE \left[  (1 - \mu_\theta)^2 \right] 
        =   (1 - 1)^2  
        = 0
```

Compared with {eq}`eq-rsmean`, we have $\risk(\hgamma, \theta) < \risk(\bar x_N, \theta)$ whenever $\sigma_\theta > 0$. But if $\mu_\theta$ diverges sufficiently from $1$, then

```{math}
:enumerated: false

\risk(\hgamma, \theta) 
        =   (1 - \mu_\theta)^2
        >  \frac{\sigma_\theta^2}{N} 
        = \risk(\bar x_N, \theta)
```
````


This example illustrates the fact that we cannot, in general, choose estimators based on lowest risk in a classical setting. One useful thing we can do, however, is exclude estimators that are always dominated in terms of risk. These are the so-called inadmissible estimators. In particular, a decision rule $d$ is called **inadmissible** for the general decision problem if there exists a second rule $e$ such that

(i) $\risk(e, \theta) \leq \risk(d, \theta)$ for all $\theta \in \Theta$ and

(ii) $\risk(e, \theta) < \risk(d, \theta)$ for at least one $\theta \in \Theta$.

A decision rule that is not inadmissible is called **admissible**

```{prf:example}

Equation {eq}`eq-mlii` tells us that the vector sample mean is inadmissible as an estimator of the mean when $L$ is quadratic and the universe of distributions is the multivariate Gaussians.
```


Another idea that comes out of the desire to minimize risk is minimax. A decision rule $d_m$ is called a **minimax** rule if

```{math}
:enumerated: false

r(d_m) \leq r(d) 
    \text{ for all } d \in \ddD
    \quad \text{where} \quad
    r(d) := \sup_{\theta \in \Theta} \risk(d, \theta)
```

A minimax rule is one that performs well in the worst possible state of the world. Minimax decision rules are an important component of robust control and the study of policy under uncertainty. See, for example, {cite}`hansen2008robustness`. On the other hand, minimax decision rules can be inadmissible (see, e.g., {cite:author}`young2005essentials` {cite:year}`young2005essentials`, p. 9).

#### Bayes Rules and Bayes Risk

The Bayesian approach to statistics leads naturally to a coherent and attractive decision theory. The horse race is not entirely fair, however, since we are granted additional information expressed as a prior. Let's review the main ideas.

Consider again the generic decision problem discussed above, but now with a prior distribution $\pi = \pi(\theta)$ over states of the world. In addition each $P_\theta$ is now restricted to be a density $p(x \given \theta)$. As in {ref}`ss-be`, we obtain the posterior for $\theta$ from the prior given data $x$ via Bayes' law:

```{math}
:enumerated: false

\pi(\theta \given x) 
    = \frac{p(x \given \theta) \pi(\theta)} {p(x)}
    \qquad 
    (x \in \xX, \; \theta \in \Theta)
```

Given decision rule $d$ with risk function $\risk(d, \theta)$, the **posterior risk** or **Bayes risk** of $d$ is

```{math}
:enumerated: false

r_\pi(d) := \int \risk(d, \theta) \pi(\theta) \diff \theta
```

A decision rule $d$ that minimizes the Bayes risk is called a **Bayes rule**. If $d$ is an estimator, then a Bayes rule is called a **Bayes estimator**. Existence of such rules is not guaranteed, but in many practical instances they are well defined and unique. We give some simple examples below.

There are two particularly attractive features of Bayes rules, and Bayes estimators in particular. One is that Bayes rules are always admissible under mild regularity conditions. See, for example, theorem 10.4.2 of {cite}`casella2002statistical`. A second is that Bayes rules can be computed as a function of observed data, meaning that we don't have to concern ourselves with how to act in situations that potentially never occur.

To understand the second point, note that we can write the joint density of $(x, \theta)$ as either $p(x \given \theta) \pi(\theta)$ or $\pi(\theta \given x) p(x)$. Using a change in order of integration followed by this identity, the Bayes risk can be expressed as

```{math}
:enumerated: false

\begin{aligned}
r_\pi(d) 
    & = \int_\Theta 
        \left\{
            \int_\xX
            L(d(x), \theta) p(x \given \theta)
            \diff x
        \right\}
        \pi(\theta) 
        \diff \theta   
    \\
    & = \int_\xX 
        \left\{
            \int_\Theta
            L(d(x), \theta) p(x \given \theta) \pi(\theta) 
            \diff \theta   
        \right\}
            \diff x
    \\
    & = \int_\xX 
        \left\{
            \int_\Theta
            L(d(x), \theta) \pi(\theta \given x) 
            \diff \theta   
        \right\}
            p(x)
            \diff x
\end{aligned}
```

It follows from this last expression that $d$ will be a Bayes rule whenever

$$
d(x) \in \argmin_{a \in \aA} 
            \int_\Theta
            L(a, \theta) \pi(\theta \given x) 
            \diff \theta   
            \quad \text{for all } x \in \xX
$$ (eq-bmv)

This is akin to the scenario encountered in dynamic programming after applying Bellman's principle of optimality: once we know the value function, we can choose the optimal action just by responding to any observed state.

````{prf:example}
:label: eg-mmed

Let $\theta$ be scalar, let $L(a, \theta) = (a - \theta)^2$ and let $\pi(\theta \given x)$ be the posterior distribution of $\theta$ given data $x$. The Bayes estimator of $\theta$ is, by {eq}`eq-bmv`, the solution to

```{math}
:enumerated: false

\min_{a \in \RR} 
            \int (a - \theta)^2 \pi(\theta \given x) \diff \theta
```

As shown in a different context in [exercise %s](#ex-minqrisk), the minimizer is the mean of $\pi(\theta \given x)$. In other words, with quadratic loss, the Bayes estimator is the mean of the posterior distribution.
````


````{prf:example}
:label: eg-mmed2

If we repeat the setting of [example %s](#eg-mmed) but with the absolute loss function $L(a, \theta) = |a - \theta|$ instead of quadratic loss, the Bayes estimator of $\theta$ is

```{math}
:enumerated: false

\argmin_{a \in \RR} 
            \int |a - \theta| \pi(\theta \given x) \diff \theta
```

As [exercise %s](#ex-mmed2) asks you to show, the solution is the median of $\pi(\theta \given x)$.
````


## Further Reading

There are many excellent general texts on statistical inference, including {cite}`casella2002statistical`, {cite}`wasserman2013all`, and {cite}`young2005essentials`. For a discussion of estimation theory from the perspective of machine learning, see {cite}`friedman2001elements` or {cite}`abu2012learning`.

## Exercises

```{exercise}
:label: ex-sn2e

Prove {eq}`eq-sn2e` in [example %s](#eg-svic).[^1]
```


```{exercise}
:label: ex-ismsv

Confirm [fact %s](#fa-ismsv).[^2]
```


```{exercise}
:label: ex-udmmse

Verify the claim in {eq}`eq-oeib` of [example %s](#eg-udmmse).
```


```{exercise}
:label: ex-cwwn

Let $x_1, \ldots, x_N$ be IID random variables uniformly distributed on $[0,\theta]$. Consider the estimator of $\theta$ defined by $\hat \theta_N := \frac{2}{N} \sum_{n=1}^N x_n$.

(i) Calculate the bias of $\hat \theta_N$.

(ii) Calculate the mean squared error of $\hat \theta_N$.

(iii) Is $\hat \theta_N$ consistent? Why or why not?
```


````{exercise}
:label: ex-cwwnex

As in [ex. %s](#ex-cwwn), let $x_1,\ldots,x_N$ be IID draws from $U[0,\theta]$, where $\theta$ is an unknown positive parameter. Let $\mu := \EE[x_n]$. Show that the estimator

```{math}
:enumerated: false

\hat \mu := \frac{N+1}{2N} \cdot \max\{x_1, \ldots, x_N\}
```

is an unbiased estimator of $\mu$.[^3]
````


```{exercise}
:label: ex-fi

Let $x_1,\ldots,x_N$ be IID and standard normal. Let $\bar x$ be the sample mean. Show that $\bar x$ and $x_i - \bar x$ are independent for $i = 1,\ldots,N$.
```


```{exercise}
:label: ex-ch09_estimator_properties-auto-1

Following on from [ex. %s](#ex-fi), let $x_1,\ldots,x_N$ be IID and standard normal, let $\bar x$ be the sample mean, and $s^2$ be the sample variance. Show that $\bar x$ and $s^2$ are independent.
```


```{exercise}
:label: ex-dmse

Confirm {eq}`eq-dmse`: Show that, for any estimator $\hat \gamma$ of $\gamma$, we have $\mse (\hat \gamma, \gamma) = \var \hat \gamma + \bias(\hat \gamma, \gamma)^2$.
```


```{exercise}
:label: ex-rauc

Let $\hat \gamma_N$ be an estimator of $\gamma \in \RR$, based on a sample of size $N$.

(i) Is it true that if $\mse(\hat \gamma_N, \gamma) \to 0$ as $N \to \infty$, then $\hat \gamma_N$ is consistent for $\gamma$? Why or why not?

(ii) Is it true that if $\hat \gamma_N$ is asymptotically unbiased and $\var \hat \gamma_N \to 0$, then $\hat \gamma_N$ is consistent for $\gamma$? Why or why not?
```


```{exercise}
:label: ex-ssdcan

Assuming the conditions of [examples %s](#eg-svcan)–[%s](#eg-ssdcan), show that {eq}`eq-avarsn` is valid.
```


````{exercise}
:label: ex-ctchi

Let $x_1,\ldots,x_N$ be IID with mean $\mu$ and variance $\sigma^2$. Let $\bar x_N$ be the sample mean, and let $\sigma_N$ be a consistent estimator of $\sigma$. What is the limiting distribution of

```{math}
:enumerated: false

y_N := N \left( \frac{\bar x_N - \mu}{\sigma_N} \right)^2
```
````


```{exercise}
:label: ex-anormic

Let $\hat \theta_N$ be an estimator of $\theta$. Show that if $\hat \theta_N$ is asymptotically normal, then $\hat \theta_N$ is consistent for $\theta$. (You can restrict attention to the scalar case.)
```


````{exercise}
:label: ex-pdead

Consider the following estimation problem. Let $\{p_{\boldtheta}\}_{\boldtheta \in \Theta}$ be a parametric class of densities. Now let a particular $\boldtheta_0 \in \Theta$ be fixed, and suppose that an asymptotically normal sequence of estimators $\{\hat \boldtheta_N\}$ exists for it, in the sense that

```{math}
:enumerated: false

\sqrt{N} (\hat \boldtheta_N - \boldtheta_0) 
        \tod \nN(\boldzero, \Sigma)
```

for some positive definite $\Sigma$. Fix a point $y$ and consider $p(y, \hat \boldtheta_N)$ as an estimator of $p(y, \boldtheta_0)$. State conditions under which this estimator is consistent and asymptotically normal. Derive the asymptotic distribution under these conditions.
````


```{exercise}
:label: ex-mmed2

Confirm the claim in [example %s](#eg-mmed2) that the minimizer is the median of $\pi(\theta \given x)$. In other words, it is the value $a \in \RR$ that satisfies

$$
\int_{-\infty}^a \pi(\theta \given x)
            = \int_a^{\infty} \pi(\theta \given x) = \frac{1}{2}
$$ (eq-mmed2)
```


### Solutions to Selected Exercises

````{solution} ex-sn2e

Let $s^2_N$, $\boldx$ and $\mu$ be as stated in the exercise. Here's one way to prove {eq}`eq-sn2e`. Let $S = \Span\{\boldone\}$, let $\boldP = \proj S$, and let $\boldM$ be the residual projection. Recalling [example %s](#eg-pvonesr), when $\boldy := \boldx - \mu \boldone$, we have $\boldP \boldy = \bar y \boldone$ and $\boldM \boldy = \boldy - \bar y \boldone$. Applying $\| \boldy \|^2 = \| \boldP \boldy \|^2 + \| \boldM \boldy \|^2$ now gives

```{math}
:enumerated: false

\sum_{n=1}^N (x_n - \mu)^2 
        = N (\bar x_N - \mu)^2 + 
            \sum_{n=1}^N (x_n - \bar x_N)^2
```

Some minor rearrangements give {eq}`eq-sn2e`.
````


````{solution} ex-ismsv

Let $\lL(\boldx) = \nN(\mu \boldone, \sigma^2 \boldI)$ and let $\boldxi := \sigma^{-1}(\boldx - \mu \boldone)$. Since the images of a standard normal random vector under a projection and its residual are independent ([fact %s](#fa-poci)) and since arbitrary deterministic functions preserve independence ([fact %s](#fa-rviifi)), it is enough to find functions $f_1$ and $f_2$ and a subspace $S$ and its corresponding projection $\boldP$ and residual $\boldM$ such that

$$
\bar x_N = f_1(\boldP \boldxi)
        \quad \text{and} \quad
        s^2_N = f_2(\boldM \boldxi)
$$ (eq-ismtc)
To this end, recall from [example %s](#eg-pvonesr) that if $S := \Span\{\boldone\}$ and $\boldP = \proj S$, then $\boldP \boldx = \bar x_N \boldone$. In addition, $\boldP \boldone = \boldone$. From this you will be able to show that $\bar x_N = \sigma \bolde_1^\T \boldP \boldxi + \mu$, where $\bolde_1$ is the first canonical basis vector. This confirms the first part of {eq}`eq-ismtc`. For the second part, apply {eq}`eq-sn2m` to get

```{math}
:enumerated: false

s^2_N = \frac{\sigma^2 \boldxi^\T \boldM \boldxi }{N}
        
```
````


````{solution} ex-cwwnex

Let $z := \max\{x_1, \ldots, x_N\}$. By the definition of the maximum, we have

```{math}
:enumerated: false

z \leq s
        \quad \text{ if and only if } \quad
            x_1 \leq s \text{ and } x_2 \leq s \text{ and } \cdots
            \text{ and } x_N \leq s
```

```{math}
:enumerated: false

\fore
        \{z \leq s \} = \cap_{n=1}^N \{x_n \leq s \}
```

```{math}
:enumerated: false

\fore
        \PP \{z \leq s \} 
        = \PP \cap_{n=1}^N \{x_n \leq s \}
        = \prod_{n=1}^N \PP \{x_n \leq s \}
```

For $s \in [0, \theta]$, the CDF of $x_n$ is $F(s) = s / \theta$, and we obtain

```{math}
:enumerated: false

\PP \{z \leq s \} = \left( \frac{s}{\theta} \right)^N
```

Differentiating shows that density $p$ of $z$ on the interval $[0,\theta]$ is $p(s) = N \frac{s^{N-1}}{\theta^N}$. The expectation of $z$ is therefore $\int_0^{\theta} s p(s) \diff s = \frac{N}{N+1} \theta$. We conclude that

```{math}
:enumerated: false

\EE[ \hat \mu ] 
        = \frac{N+1}{2N} \EE[z]  
        = \frac{N+1}{2N} \frac{N}{N+1} \theta
        = \frac{\theta}{2} 
        = \mu
```

In other words, the estimator is unbiased.
````


````{solution} ex-dmse

Adding and subtracting $\EE[\hat \gamma]$, we get

```{math}
:enumerated: false

\mse(\hat \gamma, \gamma) := \EE[ (\hat \gamma - \gamma)^2 ]
            = \EE[ (\hat \gamma - \EE[\hat \gamma] + \EE[\hat \gamma] - \gamma)^2 ]
```

Expanding the square and minor manipulations yield the desired result.
````


````{solution} ex-ssdcan

The aim is to verify the claim

```{math}
:enumerated: false

\sqrt{N}(s_N - \sigma)
        \tod \nN\left(0, \frac{m_4 - \sigma^4}{4 \sigma^2} \right)
```

under the conditions of [example %s](#eg-ssdcan). Letting $g(x) = \sqrt{x}$ and applying [theorem %s](#t-dm) along with {eq}`eq-avarsn2`, we have

```{math}
:enumerated: false

\sqrt{N}(s_N - \sigma)
        = \sqrt{N}(g(s_N^2) - g(\sigma^2))
        \tod \nN\left(0, g'(\sigma^2)^2 (m_4 - \sigma^4) \right)
```

The result follows.
````


````{solution} ex-ctchi

Let

```{math}
:enumerated: false

w_N 
        := \sqrt{N}  \frac{\bar x_N - \mu}{\sigma_N} 
        = \frac{\sigma}{\sigma_N} \sqrt{N}  \frac{\bar x_N - \mu}{\sigma}
```

Since $\sigma_N \toprob \sigma$ by assumption, [fact %s](#fa-reconpro) yields $\sigma / \sigma_N \toprob \sigma / \sigma = 1$. Applying the central limit theorem and Slutsky's theorem ([fact %s](#fa-slut)) together, we then have $w_N \tod z \sim \nN(0,1)$. By the continuous mapping theorem ([fact %s](#fa-cmtetc)), $y_N = w_N^2$ converges in distribution to $z^2$. As discussed in {ref}`ss-pnd`, the distribution of $z^2$ is $\chi^2(1)$.
````


````{solution} ex-anormic

Fix $\delta > 0$. It suffices to show that for any positive number $\epsilon$ we have

$$
\lim_{N \to \infty} \PP\{ | \hat \theta_N - \theta| > \delta \} \leq \epsilon
$$ (eq-sfco)
(Note: If $a$ is any nonnegative real number and $a \leq \epsilon$ for any $\epsilon > 0$, then $a = 0$.) To establish {eq}`eq-sfco`, fix $\epsilon > 0$. Let $z$ be zero-mean Gaussian with variance equal to the asymptotic variance of $\hat \theta_N$. Choose $M$ such that $\PP\{|z| \geq M\} \leq \epsilon$. For $N$ such that $\sqrt{N} \delta \geq M$ we have

```{math}
:enumerated: false

\PP\{ | \hat \theta_N - \theta| > \delta \} 
        = \PP\{ \sqrt{N} | \hat \theta_N - \theta| > \sqrt{N} \delta \} 
        \leq \PP\{ \sqrt{N} | \hat \theta_N - \theta| > M \}
```

Taking $N \to \infty$, applying asymptotic normality, the continuous mapping theorem ([fact %s](#fa-cmtetc)) and the definition of $M$ gives {eq}`eq-sfco`.
````


````{solution} ex-pdead

First observe that, under the conditions of the exercise, $p(y, \hat \boldtheta_N)$ is consistent for $p(y, \boldtheta_0)$ whenever $\boldtheta \mapsto p(y, \boldtheta)$ is continuous at $\boldtheta_0$. Indeed, since $\hat \boldtheta_N$ is asymptotically normal, it is also consistent. Under this continuity assumption, in probability convergence of $p(y, \hat \boldtheta_N)$ to $p(y, \boldtheta_0)$ now follows from [fact %s](#fa-reconpro).

Regarding asymptotic normality, suppose that $\boldtheta \mapsto p(y, \boldtheta)$ is differentiable at $\boldtheta_0$, and that the vector of partial derivatives $\nabla_{\boldtheta} p(y, \boldtheta_0)$ is not the zero vector. Applying {eq}`eq-dmmv`, we then have

```{math}
:enumerated: false

\sqrt{n} \{ p(y, \hat \boldtheta_N) - p(y, \boldtheta_0) \}
        \tod 
        \nN(0, \nabla_{\boldtheta} p(y, \boldtheta_0)^\T \Sigma
            \nabla_{\boldtheta} p(y, \boldtheta_0))
        
```
````


````{solution} ex-mmed2

The minimization problem in [example %s](#eg-mmed2) can alternatively be written as

```{math}
:enumerated: false

\argmin_{a \in \RR} 
           \left \{
                   \int_{-\infty}^a (a - \theta) \pi(\theta \given x) \diff \theta   
                   +
                   \int_a^{\infty} (\theta - a) \pi(\theta \given x) \diff \theta   
            \right\}
```

Differentiating with respect to $a$, setting the result to zero, and using the fact that the posterior is a density gives {eq}`eq-mmed2`.
````


[^1]: Hint: Review [example %s](#eg-pvonesr).

[^2]: Hints: Look at [fact %s](#fa-poci) and the proofs following [example %s](#eg-sdq).

[^3]: Hint: If you're stuck, review [ex. %s](#ex-ostatu) and its solution.
