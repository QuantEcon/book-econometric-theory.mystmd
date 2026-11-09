---
title: "Estimators"
label: c-statlearn
---

In probability theory we aim to deduce the likelihood of different outcomes based on known probability distributions. Statistics is the inverse problem. We want to infer unknown probability distributions from outcomes we observe. In this chapter we begin to explore these ideas.

## The Estimation Problem

Let's start off by clarifying the nature of the problem we want to study.

(ss-tep)=
### Definitions
From a pedagogical point of view, it's best to start our study of statistics in a setting where successive observations are both independent and identically distributed. Later we'll extend to dependent data. In an IID setting, the fundamental problem of econometrics and statistics is this:

```{admonition} Problem
:label: p-fpstat3

We observe independent $Z$-valued draws $\boldz_1, \ldots, \boldz_N$ from a common but unknown distribution $P \in \pP$, where $\pP$ is a class of distributions on $Z$. We wish to infer some features of $P$ from this sample.
```


The set $\pP$ is the universe of distributions we are willing to consider. It can be anything, including the set of all distributions on the outcome space. One of the main tasks of economic theory is to restrict $\pP$ and thereby narrow down the set of distributions we have to search over.

```{prf:example}

{cite}`benhabib2015wealth` study the wealth distribution in a model with idiosyncratic capital income risk. The model predicts that the wealth distribution will have a Pareto right tail. Intuitively, this is because stochastic returns on assets are multiplicative, so positive shocks are self-reinforcing, and a sequence of positive shocks can lead to very high wealth relative to the median.
```


Let's set out some notation:

- $Z$ is the **outcome space** where each **observation** $\boldz_n$ takes values.

- $\zdata$ denotes the **sample** or **data set** $(\boldz_1, \ldots, \boldz_N)$.

- $\Zdata := \times_{n=1}^N Z$ is the **sample space** in which $\zdata$ takes values.

- $\Pdata := \lL(\zdata)$ is the **joint distribution of the sample**.

Although calling $\Zdata$ the sample space follows standard nomenclature, we are overloading the terminology first used. (Unfortunately, statisticians and probabilists don't always coordinate on terminology.)

As above, $P = \lL(\boldz_n)$ is the common distribution of the observations. For as long as we continue to focus on IID data, the joint distribution $\Pdata$ will be the $N$th product of $P$, as defined.

```{prf:example}
:label: eg-siie

Let $x_1, \ldots, x_N$ be observations of labor income from a given population. We model $x_1, \ldots, x_N$ as IID draws from a common univariate distribution $P$. In our terminology, $x_n$ is an observation, the outcome space $Z$ is $\RR$, the sample space $\Zdata$ is $\RR^N$, and the sample $\zdata$ is the vector $(x_1, \ldots, x_N)$. Features of $P$ that we might wish to learn about include
```


- the mean and higher moments of $P$,

- measures of dispersion, such as the variance, or properties of tails,

- the median and other quantiles, and

- $P$ itself, or the density of $P$ if it exists.

In [example %s](#eg-siie), the individual observations are scalar. In other scenarios, observations are vector-valued:

````{prf:example}
:label: eg-fprof

Suppose that we wish to learn about the relationship between profitability and R&D spending within a group of firms. Let $\boldz_n = (x_n, y_n)$ be an observation of these two quantities at the $n$th firm. We treat the observations as IID across firms, with common marginal distribution $P = \lL(\boldz_n)$. The outcome space is $Z = \RR^2$, and the sample space is

```{math}
:enumerated: false

\Zdata := \RR^2 \times \cdots \times \RR^2 = \RR^{2 \times N}
```
````


Because the marginal distribution $P$ of the observations is now multivariate, new features of the distribution come in to play, such as

- correlations across coordinates of $P$,

- the variance–covariance matrix associated with $P$, and

- parameters controlling dependence when, say, $P$ is modeled via a copula over certain marginals.

#### Features

We referred above to "features\" of $P$ that we might be interested in estimating. Let's define a **feature** of $P$ to be an object of the form

$$
\gamma(P)
    \quad \text{for some} \quad
    \gamma \colon \pP \to S
$$ (eq-dfeature)

The set $S$ is left arbitrary to accommodate all possible features. When $P$ is understood, we'll write $\gamma(P)$ as $\gamma$. Here are some examples for univariate $P$ routinely estimated in econometric studies:

- $\gamma(P) = \int s^k P(\diff s)$, the $k$th moment of $P$.

- $\gamma(P) = \inf\setntn{s \in \RR}{P(- \infty, s] \geq 1/2}$, the median of $P$.

- $\gamma(P) = P$, when we want to estimate $P$ itself.

- $\gamma(P) =$ the density of $P$ when $P$ is absolutely continuous.

If $P$ is multivariate over $\boldz = (\boldx, y)$, then one feature of interest is the regression function $f^*(\boldx) := \EE[ y \given \boldx]$. This function is uniquely determined by $P$ (see {ref}`ss-ce`).

#### Parametric versus Nonparametric Classes

In the statistical problems listed above, it is assumed that the unknown distribution belongs to some class $\pP$. We call $\pP$ a **parametric class** if it can be expressed as

```{math}
:enumerated: false

\pP = \{P_{\boldtheta}\}_{\boldtheta \in \Theta} :=:
    \setntn{P_{\boldtheta}}{\boldtheta \in \Theta}
    \quad \text{for some } \Theta \subset \RR^K
```

In other words, a class of distributions is parametric if it can be indexed by finitely many parameters. A class of distributions is called **nonparametric** if it is not parametric.

````{prf:example}
:label: eg-normpc

Let $\pP$ be the set of all univariate normal distributions with positive variance. Equivalently,

```{math}
:enumerated: false

\pP := 
        \left\{ 
        \text{all } p \; \st \;
        p(s) = \frac{1}{\sqrt{2 \pi} \sigma}
       \exp \left\{ - \frac{(s - \mu)^2}{2\sigma^2} \right\}
       \; \text{ for some }
       \; \mu \in \RR, \; \sigma > 0
       \right\}
```

The set $\pP$ is an example of a parametric class. The parameters are $\boldtheta = (\mu, \sigma)$. A particular choice of the parameters determines (parameterizes) an element of the class.
````


```{prf:example}

If the outcome space $Z$ is finite, containing $J$ elements, then every distribution $P$ on $Z$ can be represented by $J-1$ parameters (the probability $p_j$ of each outcome but the last) and hence any family $\pP$ of distributions on $Z$ is necessarily a parametric class.
```


```{prf:example}

The set of all distributions on $\RR$ cannot be parameterized by a finite vector of parameters because the space of distributions is infinite dimensional. Hence $\pP :=$ all distributions on $\RR$ is a nonparametric class.
```


````{prf:example}

Let $\pP$ be the set of all absolutely continuous distributions on $\RR$ with finite second moment:

```{math}
:enumerated: false

\pP := \left\{ \text{all } p \colon \RR \to \RR \st
        p \geq 0, \; \int p(s)\diff s = 1, \; \int s^2 p(s) \diff s < \infty \right\}
```

This set is nonparametric.
````


Many traditional methods of inference are parametric in nature: The data are assumed to be generated by an unknown element $P_{\boldtheta}$ of a parametric class $\pP$. The aim is to estimate $\boldtheta$ using the data. Once we have an estimate $\hboldtheta$ of $\boldtheta$, we can plug this estimate back into the parametric class to obtain an estimate $P_{\hboldtheta}$ of $P_{\boldtheta}$.

The concept of a feature defined in {eq}`eq-dfeature` is a generalization of the notion of a parameter or vector of parameters. In parametric settings, the feature $\gamma$ that we are most interested in estimating is the parameter vector $\boldtheta$. If we have a good estimate $\hboldtheta$ of $\boldtheta$, we can estimate any feature $\gamma = \gamma(P_{\boldtheta})$ via $\hat \gamma = \gamma(P_{\hboldtheta})$.

Following common usage, we will sometimes refer to the $\boldtheta$ associated with the $P_\boldtheta$ that generates the data as the **true value** of the parameter vector. However, it's important to remember that this is an assumption, not a "truth.\" We *assume* that the data are generated by some member of a parametric class. This assumption can be completely false.

In our initial formulation of the estimation problem we considered estimation of features rather than just parameters because it is suboptimal to always restrict ourselves to parametric assumptions on $\pP$. We'll discuss the relative merits of parametric and nonparametric estimation more in [chapter 14](#c-rms).

(s-stats)=
### Statistics and Estimators
A statistic is any observable function of the data. More precisely, a **statistic** is any $\bB$-measurable function

```{math}
:enumerated: false

T \colon \Zdata \to S
```

that can be evaluated once the data $\zdata$ are observed. As in {eq}`eq-dfeature`, the set $S$ is left arbitrary to accommodate all the possible features that we might wish to estimate. Sometimes we write $T$ as $T_N$ to emphasize dependence on the sample size. The assumption that a statistic is $\bB$-measurable is just a basic regularity condition.

An **estimator** is a statistic used to infer some feature $\gamma(P)$ of an unknown distribution $P$. Thus a statistic becomes an estimator when paired with and compared to a feature of the distribution. Note that there is nothing in the definition of an estimator that implies it will be a sensible estimator of the target feature, let alone a good one.

````{prf:example}
:label: eg-smean

If the feature $\gamma$ we wish to infer is the mean of the marginal distribution $P$ of IID data $x_1, \ldots, x_N$, then the most common estimator is the **sample mean**

```{math}
:enumerated: false

\bar x_N := \frac{1}{N} \sum_{n=1}^N x_n
```

Formally, $\bar x_N$ is the mapping from $\Zdata = \RR^N$ to $S = \RR$ defined by

```{math}
:enumerated: false

\zdata = (x_1, \ldots, x_N)
        \mapsto
        T(x_1, \ldots, x_N)
        = \frac{1}{N} \sum_{n=1}^N x_n  \in \RR
```

This mapping is regarded as an estimator of the unknown mean $\gamma(P) = \int s P(\diff s)$.
````


````{prf:example}

The sample mean is not the only way to estimate the mean. For example, we could also use the so-called **mid-range estimator**

```{math}
:enumerated: false

m_N := \frac{\min_n x_n + \max_n x_n }{2}
```

Another option is a **truncated sample mean**, where values $x_n$ with $|x_n| \geq r$ are truncated for some specified value of $r$. The truncated sample mean is often used to estimate location parameters in heavy tailed distributions.
````


```{prf:example}

Given sample $x_1, \ldots, x_N$, let $y_n$ be the $n$th largest observation of the sample. If $N$ is the odd number $2m+1$, the **sample median** is defined as $y_{m+1}$. If $N=2m$, the sample median is $0.5(y_m + y_{m+1})$. For example,
```


``` julia
julia> median([1, 3, 5])
3.0
julia> median([2, 4, 6, 8])
5.0
```

```{prf:example}
:label: eg-smoment

Following on from [example %s](#eg-smean), a common estimator of the $k$th moment $\int s^k P(\diff s)$ of $P$ is the $k$th **sample moment** $\frac{1}{N} \sum_{n=1}^N x_n^k$.
```


```{prf:example}
:label: eg-svariance

A common estimator of the variance of $P$ is the **sample variance**

$$
s_N^2 := \frac{1}{N} \sum_{n=1}^N (x_n - \bar x_N)^2
$$ (eq-sampvar)
The standard deviation is usually estimated using **sample standard deviation**

$$
s_N 
        := \sqrt{s^2_N} 
        = \left[ 
            \frac{1}{N} \sum_{n=1}^N (x_n - \bar x_N)^2
            \right]^{1/2}
$$ (eq-sampsd)
```


In some texts the sample variance is defined as $\frac{1}{N-1} \sum_{n=1}^N (x_n - \bar x_N)^2$ rather than {eq}`eq-sampvar` in order to produce an unbiased estimator. However, {eq}`eq-sampvar` fits better with the theory in this text and is perhaps the more common definition.

```{prf:example}

Given bivariate data $\zdata = ((x_1, y_1), \ldots, (x_N, y_N))$, the **sample covariance** is the statistic

$$
\frac{1}{N} \sum_{n=1}^N (x_n - \bar x_N)(y_n - \bar y_N)
$$ (eq-sampcov)
The **sample correlation** is the sample covariance divided by the product of the two sample standard deviations. With some rearranging, this becomes

$$
\frac{\sum_{n=1}^N (x_n - \bar x_N)(y_n - \bar y_N)}
        { 
        \sqrt{ \sum_{n=1}^N (x_n - \bar x_N)^2 \sum_{n=1}^N (y_n - \bar y_N)^2 }
        }
$$ (eq-sampcorr)
```


````{prf:example}

In the case where our observations are vectors $\boldx_1,\ldots,\boldx_N$ in $\RR^K$, the sample mean is the random vector defined by

```{math}
:enumerated: false

\bar \boldx_N := \frac{1}{N} \sum_{n=1}^N \boldx_n
```

The variance–covariance matrix defined in {eq}`eq-vcd2` is most often estimated with the **sample variance–covariance matrix**

$$
\hboldSigma_N := \frac{1}{N} \sum_{n=1}^N 
            [ (\boldx_n - \bar \boldx_N) (\boldx_n - \bar \boldx_N)^\T]
$$ (eq-svarcov)
````


(ss-eds)=
### Empirical Distributions
The **empirical distribution** of a given $Z$-valued sample $\zdata = (\boldz_1, \ldots, \boldz_N)$ is the discrete distribution on $Z$ that puts equal probability $1/N$ on each sample point $\boldz_n$. Another way to state this is that $\hat P_N$ assigns to each Borel set $B \subset Z$ the number

$$
\hat P_N(B) = \frac{1}{N} \sum_{n=1}^N \1\{\boldz_n \in B\}
$$ (eq-dempd)

This is just the fraction of the sample that falls in $B$. (Refer back to [figure %s](#f-additivity) if you wish to see a visual representation.) The expectation of a function $h$ with respect to $\hat P_N$ is

$$
\int h(\bolds) \hat P_N(\diff \bolds) = \frac{1}{N} \sum_{n=1}^N h(\boldz_n)
$$ (eq-expue)

as follows from {eq}`eq-hasp2`.

```{prf:example}
:label: eg-hued

Let $\boldz_n$ be the scalar $x_n$. The sample mean can be expressed in terms of the empirical distribution as

$$
\bar x_N 
        = \frac{1}{N} \sum_{n=1}^N x_n
        = \int s \hat P_N(\diff s)
$$ (eq-hued)
In other words, the sample mean is the mean of the empirical distribution.
```


The empirical distribution is a statistic, mapping observations $\boldz_1, \ldots, \boldz_N$ into $\hat P_N = \frac{1}{N} \sum_{n=1}^N \1\{\boldz_n \in \cdot \}$, a random element of the set of all distributions on $Z$. If we think of $\boldz_1, \ldots, \boldz_N$ as independent draws from common but unknown distribution $P$, then $\hat P_N$ becomes an estimator of $P$. In particular, if the feature of $P$ we want to infer is $P$ itself, the simplest natural estimator is the empirical distribution.[^1]

With scalar data $x_1, \ldots, x_N$ we can visualize the empirical distribution by plotting its CDF. The CDF of $\hat P_N$ will be denoted in what follows by $\hat F_N$. Specializing $B$ in {eq}`eq-dempd` to $(-\infty, s]$, we get

```{math}
:enumerated: false

\hat F_N(s) = \frac{1}{N} \sum_{n=1}^N \1\{x_n \leq s\} 
    \qquad (s \in \RR)
```

The CDF $\hat F_N$ is called the **empirical cumulative distribution function**, or **ECDF**, corresponding to the sample. Graphically, $\hat F_N$ is a step function, with an upward jump of $1/N$ at each data point. {numref}`f-ecdf_beta` shows an example where draws are from the CDF $F$ of the Beta$(5,5)$ distribution. Each panel shows an observation of the ECDF, which is in turn constructed from $N = 15$ observations from $F$.

```{figure} figures/ecdf_beta.svg
:name: f-ecdf_beta

$F$ and four observations of $\hat F_N$ when $N=15$
```

#### Convergence

The empirical distribution is asymptotically an excellent estimator of $P$. If $\boldz_1, \ldots, \boldz_N$ are IID with common distribution $P$ and $\hat P_N$ is the empirical distribution, then, by the law of large numbers,

$$
\hat P_N(B) = \frac{1}{N} \sum_{n=1}^N \1\{\boldz_n \in B\}
    \toprob 
    \PP\{\boldz_n \in B\} = P(B)
$$

for any Borel set $B$. See, in particular, {eq}`eq-lln0lg`.

Specializing to the scalar case with $B := (-\infty, s]$, we have $\hat F_N(s) \toprob F(s)$ for any $s \in \RR$, where $F$ is the CDF of $P$. In fact a much stronger statement is also true. It is sometimes called the fundamental theorem of statistics, or the **Glivenko–Cantelli** theorem:

```{prf:theorem}
:label: t-fts

Let $x_1, \ldots, x_N$ be IID draws from $F$. If $\hat F_N$ is the corresponding ECDF, then

$$
\| F - \hat F_N \|_{\infty}
        \toprob 0
        \quad \text{as} \quad
        N \to \infty
$$ (eq-dev)
```


The norm here is the supremum norm, meaning that

```{math}
:enumerated: false

\| F - \hat F_N \|_{\infty}
    := \sup_{s \in \RR} |\hat F_N(s) - F(s)|
```

Roughly speaking, this is the maximal deviation over the domain. (For a definition of the supremum, see {ref}`sss-sai`.) In fact the convergence occurs "almost surely," which is a stronger notion than in probability. See theorem 11.4.2 of {cite}`dudley2002real`.

{numref}`f-ecdf_lim` illustrates the convergence in {eq}`eq-dev`. Each ECDF $\hat F_N$ is based on a single sample $x_1, \ldots, x_N$ of independent draws from $F$, where $F$ is Beta$(5,5)$. Each subplot shows 8 realizations of the ECDF at a different value of $N$. Across the subplots, the sample size is stepped up from 10 to 500.

{prf:ref}`t-fts` tells us that, at least in the IID setting, if we have an infinite amount of data, then in the limit we can learn the underlying distribution without having to impose any assumptions at all. Moreover, once we know the distribution $P$, we know any feature $\gamma = \gamma(P)$.

```{figure} figures/ecdf_lim.svg
:name: f-ecdf_lim

Realizations of $\hat F_N$ with four different sample sizes
```

Does this offer a solution to the basic problem of estimation stated? The answer is no. In practice, we only ever have a finite amount of data. Inference is about generalization, as discussed in the introduction. There's no need to generalize if we know the full population.

Moreover, when we still have only finitely many observations, the empirical distribution puts all of its mass on these points. In other words, it treats the sample like the unknown distribution, without acknowledging the fact that we have limited information. This can be viewed as an extreme form of overfitting. We'll discuss related ideas in a number of contexts throughout the book.

#### Identification

A class of distributions $\pP = \{P_{\boldtheta}\}$ indexed by $\boldtheta \in \Theta$ is called **identifiable** if the map $\boldtheta \mapsto P_{\boldtheta}$ is one-to-one (see {ref}`ss-functions`) on $\Theta$.

```{prf:example}

The class of normal densities in [example %s](#eg-normpc) is identifiable. In particular, it can be shown that if $(\mu_a, \sigma_a)$ and $(\mu_b, \sigma_b)$ are distinct vectors, then the distributions $\nN(\mu_a, \sigma^2_a)$ and $\nN(\mu_b, \sigma^2_b)$ differ at at least one point.
```


Identifiability means that the parameter vector associated with the unknown distribution can eventually be distinguished from the data. To see this, suppose that $\{P_{\boldtheta}\}$ is identified on $\Theta$ and nature generates an infinite sequence of observations $\{\boldz_n\}$ from $P = P_{\boldtheta}$. Let $\boldtheta'$ be any other vector in $\Theta$ and let $P'= P_{\boldtheta'}$. By identifiability, there exists at least one Borel set $B$ with $P(B) \neq P'(B)$. Since the empirical distribution $\hat P_N(B)$ converges to $P(B)$, we can in the limit conclude that the data are not generated by $P'$.

## Estimation Principles

We started off with the problem of estimating a feature $\gamma = \gamma(P)$ of the unknown distribution $P$ from a sample. We've introduced a number of different statistics aimed at estimating various features. Where do these estimators come from and what, if anything, do they have in common?

### The Sample Analogue Principle

Most of the estimators defined above can be derived from a simple principle. In the modern statistical literature this is often called the **plug-in method**. In econometrics it goes by a variety of names, including "analogue estimation\" (see, in particular, {cite:author}`manski1994analog` {cite:year}`manski1994analog`). We'll call it the **sample analogue principle**. The principle is:

```{math}
:enumerated: false

\text{to estimate } \gamma(P),\;
    \text{ use } \,
    \gamma(\hat P_N)
```

Here $\hat P_N$ is the empirical distribution constructed from the sample. We replace the unknown distribution $P$ with the observable distribution $\hat P_N$ and then evaluate $\gamma(\hat P_N)$.

````{prf:example}

Let $x_1, \ldots, x_N$ be draws from unknown distribution $P$. Suppose we want to estimate the mean $\gamma(P) := \int sP(\diff s)$. The sample analogue principle tells us to replace $P$ with $\hat P_N$, which gives

```{math}
:enumerated: false

\gamma(\hat P_N) = \int s \hat P_N(\diff s) = \bar x_N
```

(The last equality is just {eq}`eq-hued`.) Thus the sample mean is the estimator of the mean produced by the sample analogue principle.
````


The $k$th sample moment applies the same principle to estimate the $k$th moment.

````{prf:example}

When it exists, the variance of $P$ can be written as

```{math}
:enumerated: false

\sigma^2 = \gamma(P) = \int \left[ t - \int s P(\diff s) \right]^2
        P(\diff t)
```

Applying the sample analogue principle leads to the estimator

```{math}
:enumerated: false

\int \left[ t - \int s \hat P_N(\diff s) \right]^2 \hat P_N(\diff t)
            = \frac{1}{N} \sum_{n=1}^N (x_n - \bar x_N)^2
```

This is precisely the sample variance from [example %s](#eg-svariance).
````


The estimators in these two examples are relatively straightforward, but the set of all applications of the sample analogue principle is very broad. As we'll see below, estimation methods that can be obtained as special cases include least squares regression, maximum likelihood, the method of moments, and the generalized method of moments.

(sss-blpsa)=
#### Best Linear Prediction
Here's an example that illustrates the generality of the sample analogue principle. Recall the best linear prediction problem from {ref}`sss-blp`, where $\alpha$ and $\beta$ are chosen to minimize $\EE[ (y - \alpha - \beta x)^2 ]$. Letting $P$ be the distribution of $(x, y)$, we can write the problem as

$$
(\alpha^*, \beta^*)
    = \gamma(P) := \argmin_{\alpha, \beta \in \RR} \int [ (t - \alpha - \beta
    s)^2 ] P(\diff s, \diff t)
$$ (eq-lstmp)

There is a corresponding statistical problem when the underlying distribution $P$ is unknown. The problem is to produce the best linear predictor based only on a sample $\zdata = ((x_1, y_1), \ldots, (x_N, y_N))$ of observations from $P$. We proceed as follows: Given $\zdata$, we form the empirical distribution

```{math}
:enumerated: false

\hat P_N(B) = \frac{1}{N} \sum_{n=1}^N \1\{(x_n, y_n) \in B\}
```

Applying the sample analogue principle, we plug $\hat P_N$ into {eq}`eq-lstmp`, yielding the estimator

$$
(\hat \alpha_N, \hat \beta_N)
    = \gamma(\hat P_N) 
    = \argmin_{\alpha, \beta \in \RR} 
        \sum_{n=1}^N [ (y_n - \alpha - \beta x_n)^2 ]
$$ (eq-lstmpe)

(The term $\frac{1}{N}$ has been dropped from the right-hand side of {eq}`eq-lstmpe` because it doesn't affect the minimizers.) This is the simple linear bivariate least squares problem. The minimizers are

$$
\hat \beta_N = \frac{\sum_{n=1}^N (x_n - \bar x_N)(y_n - \bar y_N)}
                    {\sum_{n=1}^N (x_n - \bar x_N)^2}
    \quad \text{and} \quad
    \hat \alpha_N = \bar y_N 
        - \hat \beta_N \bar x_N
$$ (eq-solse)

(See [ex. %s](#ex-solse).) Comparing with {eq}`eq-slpp`, we see that these estimators are themselves the sample analogues of the underlying features $\alpha^*$ and $\beta^*$.

(sss-lim)=
#### Limitations
While the sample analogue principle produces sensible estimators in many situations, there are also instances where it fails completely. Let's look at one example.

Let $\pP$ be the set of absolutely continuous distributions on $\RR$, so that the distribution $P$ generating the data is assumed to have a density. The density of $P$ can be regarded as the feature

$$
\gamma(P) = DP
$$ (eq-denasf)

Here $DP :=$ the derivative of the CDF $F$ of $P$.

Let $\hat P_N$ be the empirical distribution from any given sample. What happens when we plug this empirical distribution into the right-hand side of {eq}`eq-denasf`? In other words, what information do we get when we differentiate the function $\hat F_N$? Since $\hat F_N$ is a step function, the derivative is zero everywhere, except at a finite number of jump points where the derivative is undefined. Thus the density estimate produced by the sample analogue principle is useless. This statement remains true even when $N$ is enormous.

One way to understand this outcome is that when we wish to estimate complex objects with finite samples, we need to combine the data with some kind of **regularization**. Loosely speaking, this means that we penalize complexity in our search for a solution, or impose some kind of "smoothing\" a priori. This idea has a long history in the field of numerical methods, where certain inverse problems are unstable, or "ill-posed,\" and their solution requires regularization.

Regularization is an essential element of statistics too. It says that we should not grant the empirical distribution equal status with the unknown true distribution. That is, we should regard the empirical distribution only as partial information, and seek to combine it with some form of prior information or external theory. Regularization is one application of prior information: It applies the idea that probability mass most likely falls in places other than just the sample points observed so far. We pick up the topic of regularization again in [chapter 14](#c-rms).

(ss-erm)=
### Empirical Risk Minimization
Let's now look at another inductive principle that combines the sample analogue principle with additional structure. The principle is called empirical risk minimization, or ERM. The terminology and main concepts come from the machine learning literature. Nonetheless, many standard estimators in econometrics are special cases of ERM, including maximum likelihood and least squares. The ideas discussed here help to frame these tools as part of a broader set of methods.

Consider a setting where we observe an input $\boldx \in \RR^K$ to a system, followed by a scalar output $y$. Both are random variables and the joint distribution of $\boldz := (\boldx, y)$ is $P$. Our aim is to predict new output values from observed input values. We'll do this by choosing a function $f$ such that $f(\boldx)$ is our prediction of $y$ once $\boldx$ is observed. In the machine learning literature, $f$ is called a **prediction rule**. In economics, $f$ is called a **strategy** or **policy function**.

Incorrect prediction incurs a loss. The size of this loss is written as $L(y, f(\boldx))$. The function $L$ is called the **loss function**. Common choices for the loss function include

- the **quadratic loss function** $L(y, f(\boldx)) = (y - f(\boldx))^2$,

- the **absolute deviation loss function** $L(y, f(\boldx)) = |y - f(\boldx)|$, and

- the **discrete loss function** $L(y, f(\boldx)) = \1\{y \neq f(\boldx)\}$.

Given loss function $L$, we consider choosing $f$ so as to minimize the **prediction risk**, which is defined as the expected loss

$$
R(f) := \EE L(y, f(\boldx))
$$ (eq-rf)

Sometimes $R(f)$ is also called the expected prediction error.[^2] The expectation in {eq}`eq-rf` is computed using the joint distribution $P$ of $(\boldx, y)$. Sometimes we'll write $\EEP L(y, f(\boldx))$ to emphasize this.

```{prf:example}
:label: eg-ebrf2

Let $L$ be the quadratic loss function. As shown in {ref}`ss-ce`, the minimizer of {eq}`eq-rf` is the regression function, defined at $\boldx$ by $f^*(\boldx) = \EE[y \given \boldx]$. In particular, from {eq}`eq-ebrf`, for any alternative policy $g$ we have

$$
R(g) =  R(f^*) + \EE[ (f^*(\boldx) - g(\boldx))^2 ]
$$ (eq-ebrf2)
In other words, with quadratic loss, good prediction equates to choosing $g$ to be close to the regression function, thereby minimizing the second term on the right-hand side of {eq}`eq-ebrf2`. In view of {eq}`eq-ebrf2`, the term $R(f^*)$ represents a lower bound for prediction risk.
```


In a statistical setting we are constrained in our ability to minimize risk by the fact that we cannot evaluate $\EEP$, as this requires knowledge of the joint distribution $P$. Suppose, however, that we have access to data $\boldz_1, \ldots, \boldz_N$, where each pair $\boldz_n = (\boldx_n, y_n)$ is an independent draw from $P$. To make use of the sample, we will apply the sample analogue principle, replacing $P$ in {eq}`eq-rf` with the empirical distribution $\hat P_N$. This produces the new objective function

$$
\Remp(f) 
    := \mathbb{E}_{\hat P_N} L(y, f(\boldx))
    = \frac{1}{N} \sum_{n=1}^N L(y_n, f(\boldx_n))
$$ (eq-erf)

This function is called the **empirical risk**. The empirical risk is the prediction risk evaluated under the empirical distribution of the sample rather than the true distribution $P$. Minimizing the empirical risk is called **empirical risk minimization** (ERM).

There is an additional step. When we choose a decision rule the problem that we actually solve is

$$
\hat f = \argmin_{f \in \hH} \Remp(f)
$$ (eq-minempr)

We have restricted the domain to a class of functions $\hH$. This set of functions is called the **hypothesis space**. Choice of $\hH$ is structure we impose on the estimation problem.

It might seem at first pass that we should set $\hH$ to be the set of all functions $f \colon \RR^K \to \RR$. After all, if the risk-minimizing function $f^* := \argmin_f R(f)$ is not in $\hH$, as visualized in [figure %s](#f-circles), then the solution to {eq}`eq-minempr` is not equal to $f^*$, and we are making a suboptimal choice.

```{figure} figures/fig_circles.svg
:name: f-circles

Choosing the hypothesis space
```

Actually this reasoning is false. In fact we usually want to be quite restrictive in our choice of $\hH$. This is because minimizing the empirical risk is not the same as minimizing the prediction risk. Occam's razor comes in to play here: we are solving a complex problem on the basis of limited information. It would be a serious mistake to act as if we had unlimited information. These ideas are explored further in {ref}`ss-ermls`.

```{prf:example}
:label: eg-lsqerm

Specializing to scalar $x$ and quadratic loss function $L(y, f(x)) = (y - f(x))^2$, and observing that the term $\frac{1}{N}$ makes no difference to the solution, the ERM problem becomes

$$
\min_{f \in \hH} \; \sum_{n=1}^N (y_n - f(x_n))^2
$$ (eq-ermlsq)
For obvious reasons, this optimization problem is called the **least squares problem**. If we specialize $\hH$ to be the set of affine functions

$$
\llL := \{ \text{ all functions of the form } \ell(x) = \alpha + \beta x \}
$$ (eq-allaff)
then the problem reduces to the **simple linear least squares problem**

$$
\min_{\ell \in \, \llL} \; \sum_{n=1}^N (y_n - \ell(x_n))^2
        = \min_{\alpha, \, \beta } \, \sum_{n=1}^N  (y_n - \alpha - \beta x_n)^2
$$ (eq-llsp)
This is exactly the same minimization problem we obtained in {eq}`eq-lstmpe`. The solutions are given in {eq}`eq-solse`.
```


While we ended up with the same minimization problem in two different ways, the reasoning here is more natural. When we first came to this problem in {eq}`eq-lstmpe` it was by applying the sample analogue principle to the problem of finding best linear predictors. Here we are simply trying to find the best predictor. At the same time we have restricted ourselves to linear approximations in recognition of the fact that we are basing our estimate on limited information. Linearity implements Occam's razor by bounding the complexity of the learning algorithm.

(sss-qreg)=
#### Quantile Regression
Quantile regression is used when we want to estimate a quantile of a given distribution based on a set of predictive variables. For example, quantile regression has been used to estimate quantiles of CEO compensation as a function of the market value of their firms ({cite:author}`koenker2001quantile` {cite:year}`koenker2001quantile`). In this section we show how quantile regression can be regarded as a special case of empirical risk minimization.

To describe the idea, let $F$ be a strictly increasing CDF on $\RR$ and let $\tau \in (0, 1)$ be given. Recall from {ref}`ss-quantile` that the $\tau$th quantile of $F$ is the $\xi$ that solves $F(\xi) = \tau$. Although it might not be immediately obvious, the $\tau$th quantile can also be defined as the solution to the optimization problem

$$
\min_{\xi \in \RR} \, \EE L_\tau (y, \xi)
$$ (eq-altxi)

where $y$ is a random variable with distribution $F$ and

```{math}
:enumerated: false

L_\tau (y, \xi)
    := | (y - \xi) (\tau - \1\{y < \xi\}) |
```

{prf:ref}`ex-quantreg` asks you to confirm that the solution to {eq}`eq-altxi` is precisely the $\xi$ that solves $F(\xi) = \tau$.

Suppose now that we want to estimate the $\tau$th quantile of $F$ using some input variable $x$ (e.g., a CEO compensation quantile as a function of firm size). Motivated by {eq}`eq-altxi`, we can frame the search for a suitable function as a problem of minimizing the prediction risk

$$
R(f) := \, \EE L_\tau (y, f(x))
$$ (eq-altxi2)

where $L_\tau$ is as defined above. If we employ the principle of empirical risk minimization with data set $(x_1, y_1), \ldots, (x_N, y_N)$, we get $\min_{f \in \hH} \; \sum_{n=1}^N L_\tau (y_n , f(x_n))$, where $\hH$ is the hypothesis space. When $\hH = \llL$ as defined in {eq}`eq-allaff`, the ERM problem is

```{math}
:enumerated: false

\min_{\alpha, \, \beta } \, \sum_{n=1}^N 
    | (y_n - \alpha - \beta x_n) (\tau - \1\{y_n < \alpha + \beta x_n\}) |
```

This is the standard expression for the quantile regression problem.

If $\tau = 0.5$, then the objective function is proportional to $\sum_{n=1}^N | y_n - \alpha - \beta x_n |$. This case is called **median regression** or **least absolute deviation regression**.

(ss-ermls)=
### The Choice of Hypothesis Space
Let's return to the issue of choosing the hypothesis space. We begin with a straightforward result:

````{prf:property}
:label: fa-emprd

Let $\hH_1$ and $\hH_2$ be two hypothesis spaces. If $\hat f_i$ is the empirical risk minimizer over $\hH_i$ as defined in {eq}`eq-minempr`, then

```{math}
:enumerated: false

\hH_1 \subset \hH_2 \implies \Remp(\hat f_1) \geq \Remp(\hat f_2)
```
````


In other words, we can always decrease empirical risk by increasing the hypothesis space. The reason that [fact %s](#fa-emprd) holds is that we are minimizing over a larger set.

At the same time, our actual objective is to minimize the prediction risk. The prediction risk of a given predictor $f$ measures the **out-of-sample fit** of $f$, which is the expected performance of $f$ when confronted with new data. A statistical procedure that takes a given data set and produces a predictor $f$ with low prediction risk (relative to the benchmark $f^*$—see {eq}`eq-ebrf2`) means that we have succeeded in meeting the goal we set at the start of the book: we have *generalized* from data.

Of course, prediction risk is unobservable. Sometimes it is tempting to use empirical risk $\Remp(f)$ as an estimator of $R(f)$. However, empirical risk, which measures **in-sample fit** of $f$, is a downward biased estimator of risk. In particular, as we vary the hypothesis space, prediction risk can rise even when empirical risk is falling.

Let's illustrate this by way of an example, where empirical risk is minimized over progressively larger hypothesis spaces. The model we will consider is one that generates input-output pairs via

$$
x \sim U[-1,1] \;\; \text{ and then } \;\; y = \cos(\pi x) + u
    \;\; \text{ where }\;\;  u \sim N(0,1)
$$ (eq-csmod)

Here $U[-1,1]$ is the uniform distribution on the interval $[-1,1]$. Our hypothesis spaces for predicting $y$ from $x$ will be sets of polynomial functions. To fix notation, let $\pP_d$ be the set of all polynomials of degree $d$. That is,

```{math}
:enumerated: false

\pP_d := \{ \text{ all functions } f_d(x) = c_0 x^0 + c_1 x^1 + \cdots c_d
    x^d \text{ where each } c_i \in \RR \}
```

This sequence of hypothesis spaces is increasing, in the sense that

```{math}
:enumerated: false

\pP_1 \subset \pP_2 \subset \pP_3 \subset \cdots
```

Indeed, if $f$ is a polynomial of degree $d$, then $f$ can be represented as a polynomial of degree $d+1$ just by setting the last coefficient $c_{d+1}$ to zero. The set of linear functions $\llL$ defined in {eq}`eq-allaff` is equal to $\pP_1$.

If we seek to predict $y$ from $x$ using quadratic loss and the set $\pP_d$ as our candidate functions, the risk minimization problem is

$$
\min_{f \in \pP_d} R(f) 
    \quad \text{where} \quad 
    R(f) = \EE[ (y - f(x))^2 ]
$$ (eq-rminpoly)

while the empirical risk minimization problem is

$$
\min_{f \in \pP_d} \Remp(f) 
    \quad \text{where} \quad 
    \Remp(f) = \frac{1}{N} \sum_{n=1}^N (y_n - f(x_n))^2
$$ (eq-erminpoly)

To illustrate the difference between risk and empirical risk, we first generate $N=25$ data points from the model {eq}`eq-csmod`. Taking this as our data set, we then solve {eq}`eq-erminpoly` repeatedly, once for each $d$ in $1,2,\ldots,20$. The solution to the $d$th minimization problem is denoted $\hat f_d$, and is, by construction, a polynomial of degree $d$. Finally, we compare the risk $R(\hat f_d)$ and empirical risk $\Remp(\hat f_d)$.[^3] The results are in [figure %s](#f-rvempr).

```{figure} figures/rvempr.svg
:name: f-rvempr

Risk and empirical risk as a function of $d$
```

As expected, empirical risk falls monotonically with $d$. But the risk decreases slightly and then increases rapidly. For large $d$, the minimizer $\hat f_d$ of the empirical risk is associated with very high risk in the sense of large expected loss.

We can get a sense for what is happening by plotting the data and the functions. In [figures %s](#f-ofit1)–[%s](#f-ofit14), the $N$ data points are plotted alongside the function $y = \cos(\pi x)$ from the true model {eq}`eq-csmod` in black, and fitted polynomial $\hat f_d$ in red. The function $y = \cos(\pi x)$ is the risk minimizer, and represents the ideal prediction function. In [figure %s](#f-ofit1) we have $d=1$, and the fitted polynomial $\hat f_1$ is the linear regression line. In [figures %s](#f-ofit3), [%s](#f-ofit11) and [%s](#f-ofit14) we have $d = 3$, $d=11$ and $d=14$ respectively. The fitted polynomials are $\hat f_3$, $\hat f_{11}$ and $\hat f_{14}$.

On the one hand, when $d=1$, the hypothesis space $\pP_d = \pP_1$ is small and no function in this class can effectively fit the underlying model. This is called **underfitting**, and is reflected in the poor fit of the red line to the black line in [figure %s](#f-ofit1).

When $d=3$, the class of functions $\pP_d = \pP_3$ is considerably larger. Given that the data are relatively noisy, and that we only have 25 observations, the fit of the function is good ([figure %s](#f-ofit3)). If we look at the risk for $d=3$ on the right-hand side of [figure %s](#f-rvempr), we see that it is lower than for $d=1$.

On the other hand, if we take $d$ larger, the fit to the underlying model becomes poor and the risk is high. Examining [figures %s](#f-ofit11) and [%s](#f-ofit14), which correspond to $d=11$ and $d=14$, we see that the fitted polynomial has been able to fit the observed data closely, passing near many of the data points. Too much emphasis has been given to this particular realization of the data. When a new input $x$ is drawn, the prediction $\hat f_d(x)$ is likely to be a poor predictor of $y$. This situation is called **overfitting**.

```{figure} figures/ofit1.svg
:name: f-ofit1

Fitted polynomial, $d=1$
```

```{figure} figures/ofit3.svg
:name: f-ofit3

Fitted polynomial, $d=3$
```

```{figure} figures/ofit11.svg
:name: f-ofit11

Fitted polynomial, $d=11$
```

```{figure} figures/ofit14.svg
:name: f-ofit14

Fitted polynomial, $d=14$
```

To summarize, the choice of $\hH$ is central to our ability to generalize from the data. When $\hH$ is too small, no function in $\hH$ provides a good fit to the regression function. If this is so, the prediction risk cannot be made small regardless of the data. This can be seen directly from {eq}`eq-ebrf2`. Any estimate $\hat f$ in $\hH$ must satisfy

```{math}
:enumerated: false

R(\hat f) \geq R(f^*) + \min_{f \in \hH} \EE[ (f^*(\boldx) - f(\boldx))^2 ]
```

Conversely, if $\hH$ is too large, then we can attain low empirical risk but the prediction risk itself is large.

Of course, in statistical applications we do not have the luxury of knowing the true data-generating process when we choose $\hH$. The best scenario is that we have firm theory that guides us to a suitable hypothesis space. The worst scenario is that we have no idea and choose blindly. Once again, the message is that statistical learning equals prior knowledge plus data.

## Some Parametric Methods

In this section we review some standard parametric estimation methods, including maximum likelihood, Bayesian estimation, and the generalized method of moments.

(ss-ti)=
### Maximum Likelihood
One standard approach to deriving estimators in a parametric setting is the **principle of maximum likelihood**. Let's begin with a simple case and then progress to more complex applications.

Suppose that we observe a single draw $x$ from distribution $\nN(\mu,1)$ where $\mu$ is unknown. The variance is known and set to $1$ for simplicity. Our task is to guess the value of $\mu$ given the observation $x$. A guess $\hat \mu$ of $\mu$ also pins down the distribution $\nN(\hat \mu, 1)$, so we could equivalently state that our job is to guess the distribution that generated $x$. In guessing this distribution, if we centered it around some number $\hat \mu$ much larger than $x$, say, then our observed data point $x$ would be an "unlikely" outcome for this distribution. See [figure %s](#f-norm1). The same logic would apply if we centered the density at a point much smaller than $x$.[^4]

```{figure} figures/fig_norm1.svg
:name: f-norm1

Maximizing the likelihood
```

In fact, in the absence of any additional information, the most obvious guess would be that the normal density is centered on $x$. To center the density on $x$, we just set $\hat \mu = x$.

Maximum likelihood formalizes these steps: The density of $x$ is the density $p(s ; \mu)$ of the distribution $\nN(\mu, 1)$. Plugging the observed value $x$ into this density gives

```{math}
:enumerated: false

p(x ; \mu) = (2 \pi )^{-1/2} 
   \exp \left\{ - \frac{(x - \mu)^2}{2} \right\}
```

Think of $p(x ; \mu)$ as representing the probability of realizing our sample point $x$. The principle of maximum likelihood suggests that we take as our guess $\hat \mu$ of $\mu$ the value that maximizes this probability. It's not hard to show that $\hat \mu = x$ is the maximizer. This coincides with our intuition from [figure %s](#f-norm1).

The same principle applies when the data are independent draws $x_1, \ldots, x_N$ from $\nN(\mu, 1)$, where $\mu$ is unknown. The joint density of the sample is the product of the marginals. Plugging the sample values into the joint density, one then maximizes the joint density with respect to $\mu$:

$$
\hat \mu := \argmax_{\mu \in \RR} 
    \;\;
    \frac{1}{(2 \pi)^{N/2}}
    \prod_{n=1}^N \exp \left\{ - \frac{(x_n - \mu)^2}{2} \right\}
$$ (eq-mnd)

The maximizer $\hat \mu$ is precisely the sample mean of $x_1, \ldots, x_N$ (see [ex. %s](#ex-maxjn)).

We can generalize these ideas in several ways. Let's suppose now that the data $x_1,\ldots,x_N$ has joint density $p$ in the sense of {eq}`eq-defjd00`. We will assume that $p = p(\cdot \,; \boldtheta)$ is a member of a parametric class $\pP$ indexed by parameter vector $\boldtheta \in \Theta$. Each choice of $\boldtheta$ pins down a particular density $p = p(\cdot \,; \boldtheta)$, but the value of $\boldtheta$ that generated the data is unknown.

In this setting, the **likelihood function** is $p$ evaluated at the sample $x_1,\ldots,x_N$, and regarded as a function of $\boldtheta$:

$$
L(\boldtheta) := p(x_1,\ldots,x_N ; \boldtheta)
    \qquad (\boldtheta \in \Theta)
$$ (eq-deflf0)

The principle of maximum likelihood tells us to estimate $\boldtheta$ by maximizing $L(\boldtheta)$ over $\boldtheta \in \Theta$. A statistic $\hat \boldtheta$ is called a **maximum likelihood estimate** (MLE) of $\boldtheta$ if

$$
\hat \boldtheta \in \argmax_{\boldtheta \in \Theta} L(\boldtheta)
$$ (eq-defmle)

This is equivalent to maximizing the **log likelihood function**

```{math}
:enumerated: false

\ell(\boldtheta) := \ln (L(\boldtheta)) 
    \qquad (\boldtheta \in \Theta)
```

The set of MLEs can in general be a singleton, contain multiple elements or be empty.

In the preceding discussion, $p$ was a density function, but it can be a probability mass function as well. Examples are given below.

To implement maximum likelihood estimation, we need to compute the joint distribution of the data. As we saw in {eq}`eq-mnd`, when the data points are independent this is easy because the joint density $p$ is the product of the marginals. More generally, if each $x_n$ is drawn independently from fixed arbitrary (marginal) density $p_n(\cdot \,; \boldtheta)$ on $\RR$, then

$$
L(\boldtheta) = \prod_{n=1}^N p_n(x_n ; \boldtheta)
    \quad \text{and} \quad
    \ell(\boldtheta) = \sum_{n=1}^N \ln p_n(x_n ; \boldtheta)
$$ (eq-deflf)

If each data point is multivariate, just replace $x_n$ with $\boldx_n$.

```{prf:example}
:label: eg-mlen

Suppose that $x_1, \ldots, x_N$ are IID draws from a normal distribution $\nN(\mu, v)$ with $\boldtheta = (\mu, v)$ unknown. The log likelihood function is

$$
\ell(\mu, v) 
            = -\frac{N}{2} \ln (2 \pi v)
            - \frac{1}{2} \sum_{n=1}^N \frac{(x_n - \mu)^2}{v}
$$ (eq-mnds)
Joint maximization over $(\mu, v)$ gives the maximum likelihood estimators

$$
\hat \mu = \frac{1}{N} \sum_{n=1}^N x_n
        \quad \text{and} \quad
        \hat v = \frac{1}{N} \sum_{n=1}^N (x_n - \bar x_N)^2
$$ (eq-mlen)
Thus, the MLEs of $\mu$ and $v$ are the sample mean and sample variance respectively.
```


Maximum likelihood estimation is a much celebrated theory of estimation. MLE estimators typically have excellent asymptotic properties. But good finite sample properties are not guaranteed (more about this soon) and the attractive asymptotic theory is dependent on correct specification of the underlying parametric class. Indeed we need to bring a lot of knowledge to the table just to form the estimator. (To specify the likelihood function, we must specify the entire joint distribution of the sample.)

````{prf:example}

Let $x_1, \ldots, x_N$ be IID observations from some unknown density on $\RR$. If the class of candidate distributions is too large, maximum likelihood estimation will fail. For example, there is no solution to the problem

```{math}
:enumerated: false

\max_p \prod_{n=1}^N p(x_n)
```

where the maximization is over all densities. We can make this expression arbitrarily large by choosing $p$ to concentrate its mass in small neighborhoods around each sample point. The limit of this process is not a density.
````


(ss-cmle)=
#### Conditional Maximum Likelihood
As in {ref}`ss-erm`, suppose that we observe inputs $\boldx_1, \ldots, \boldx_N$ to some system and corresponding outputs $y_1, \ldots, y_N$. The pairs $(\boldx_n, y_n)$ are assumed to be IID. Our aim is to estimate $\boldtheta$ in $p(y \given \boldx ; \boldtheta)$ in order to pin down the conditional density of $y$ given $\boldx$.

The principle of maximum likelihood tells us to maximize

```{math}
:enumerated: false

\ell(\boldtheta) = \sum_{n=1}^N \ln p(\boldx_n, y_n; \boldtheta)
    \quad \text{where} \quad
    p := \; \text{the joint density of } (\boldx_n, y_n)
```

Letting $\pi$ be the marginal density of $\boldx$, we can use the decomposition {eq}`eq-decompjd` to write

```{math}
:enumerated: false

p(\boldx, y; \boldtheta)
    = p( y \given \boldx; \boldtheta) \, \pi(\boldx)
```

The density $\pi$ is unknown but we have not parameterized it because we aren't trying to estimate it. We can now rewrite the log likelihood as

```{math}
:enumerated: false

\ell(\boldtheta) 
    = \sum_{n=1}^N \ln p(y_n \given \boldx_n; \boldtheta) 
    + \sum_{n=1}^N \ln \pi(\boldx_n)
```

The second term on the right-hand side is independent of $\boldtheta$ and as such it does not affect the maximizer. Hence the MLE is

```{math}
:enumerated: false

\argmax_{\boldtheta \in \Theta} \sum_{n=1}^N \ln p(y_n \given \boldx_n; \boldtheta)
```

The objective function here is called the **conditional log likelihood**. The preceding argument tells us that, when we want to estimate parameters in the conditional density of $y$ given $\boldx$, we can maximize the conditional log likelihood directly, which is simpler and more direct than maximizing the full log likelihood.

````{prf:example}

Consider a discrete response model with binary output $y_n$, where $y_n=1$ indicates that the $n$th individual in a sample of women participates in the labor force. This decision is influenced by a vector $\boldx_n$ measuring characteristics such as income from the rest of the household. Let

```{math}
:enumerated: false

q(\bolds) 
        := \PP\{ y = 1 \given \boldx = \bolds \}
        \qquad \left(\bolds \in \RR^K \right)
```

One modeling approach is to take $q(\bolds) = F(\boldbeta^\T \bolds)$, where $\boldbeta$ is a vector of parameters and $F$ is a specified CDF. We can then write

```{math}
:enumerated: false

\PP\{ y = i \given \boldx = \bolds\} 
        = F(\boldbeta^\T \bolds)^i 
            (1 - F(\boldbeta^\T \bolds))^{1-i}
        \quad \text{for } \bolds \in \RR^K \text{ and } i \in \{0, 1\}
```

This is the conditional PMF of $y$ given $\boldx$, so the conditional log likelihood of the sample is

```{math}
:enumerated: false

\begin{aligned}
\ell(\boldbeta) 
        & = \sum_{n=1}^N \ln 
            [ F(\boldbeta^\T \boldx_n)^{y_n} (1 - F(\boldbeta^\T \boldx_n))^{1-y_n} ]
        \\
        & = \sum_{n=1}^N 
            y_n \ln F(\boldbeta^\T \boldx_n) + \sum_{n=1}^N (1 - y_n) \ln (1 -
            F(\boldbeta^\T \boldx_n))
\end{aligned}
```

If $F$ is the standard normal CDF $\Phi$, then this model is called the **probit** model. If $F$ is the logistic CDF $F(\bolds) = 1/(1 + e^{-s})$, then it is called the **logit** model.
````


To find the MLE in this example, we can differentiate $\ell$ to obtain first order conditions, but this does not in general lead to an analytical solution. Instead, numerical optimization of the likelihood is required. Discussion of numerical optimization is given in {ref}`ss-nra`.

### Maximum Likelihood via ERM

Maximum likelihood is a special case of empirical risk minimization. To see this, suppose that we wish to learn unknown density $q$ on the basis of observations from $q$. We take our loss function to be

```{math}
:enumerated: false

L(p, x) := - \ln p(x)
```

The meaning is this: if our guess of $q$ is $p$ and the value $x$ is realized, then our loss is $- \ln p(x)$. Thus we suffer large loss when $p$ puts small probabilities near where $x$ is realized. The corresponding risk function is

```{math}
:enumerated: false

R(p) = \EE_q [L(p, x)] = - \int \ln [p(s)] q(s) \diff s
```

Now suppose that we observe IID draws $x_1, \ldots, x_N$ from $q$. To estimate $q$, the ERM principle indicates we should solve for

```{math}
:enumerated: false

\hat p 
    := \argmin_{p \in \pP} \; 
        \left\{ \frac{1}{N} \sum_{n=1}^N - \ln p(x_n) \right\}
    = \argmax_{p \in \pP} \; 
        \left\{ \sum_{n=1}^N \ln p(x_n) \right\}
```

where $\pP$ is a hypothesis space of densities. To clarify the connection with maximum likelihood, take $\pP$ to be a parametric class $\{p(\cdot \,; \boldtheta)\}_{\boldtheta \in \Theta}$. Choosing our estimate $\hat p$ of $q$ now reduces to choosing an estimate $\hat \boldtheta$ of $\boldtheta$. Rewriting our optimization problem for this case, we obtain

```{math}
:enumerated: false

\hat \boldtheta 
    = \argmax_{\boldtheta} \; 
        \left\{ \sum_{n=1}^N \ln p(x_n ; \boldtheta) \right\}
    = \argmax_{\boldtheta} \; \ell(\boldtheta)
```

Here $\ell$ is the log-likelihood. It follows from this expression that the ERM estimator is precisely the maximum likelihood estimator.

Incidentally, the minimizer of the true risk $R(p)$ is the unknown density $q$. To see this, let us first transform our expression for the risk function to

```{math}
:enumerated: false

R(p) =   \int \ln \left[ \frac{q(s)}{p(s)} \right] q(s) \diff s 
        - \int \ln[q(s)] q(s) \diff s
```

The term on the far right is called the **entropy** of $q$, and does not involve $p$. Hence minimization of the risk comes down to minimization of

$$
D(q, p) := \int \ln \left[ \frac{q(s)}{p(s)} \right] q(s) \diff s
$$ (eq-kldist)

This quantity is called the **Kullback–Leibler (KL) deviation** between $q$ and $p$. The KL deviation is possibly infinite, always nonnegative, and zero if and only if $p = q$.[^5] It follows that the unique minimizer of the risk is the true density $q$.

### The Method of Moments and GMM

Suppose that we wish to estimate a vector $\boldtheta$ that solves an equation of the form

$$
g(\boldtheta) = \EE h( \boldx )
$$ (eq-mm)

In this expression, $g$ and $h$ are observable vector-valued functions (taking values in $\RR^K$, say). We cannot solve this expression because the distribution $P$ of $\boldx$ is unknown, and hence the expectation cannot be evaluated. If, however, we have observations $\boldx_1, \ldots, \boldx_N$ from $P$ then we can apply the sample analogue principle and replace the population expectation with the expectation under the empirical distribution. This yields the **method of moments estimator**, which is the solution $\hat \boldtheta$, if it exists, to the equation

$$
g(\hat \boldtheta) = \frac{1}{N} \sum_{n=1}^N h(\boldx_n)
$$ (eq-mme)

````{prf:example}
:label: eg-estpareto

The mean of the Pareto distribution with scale parameter $s_0 = 1$ and shape parameter $\alpha$ is $\alpha / (\alpha - 1)$. (The condition for existence of the mean is $\alpha > 1$.) Letting $g(\alpha) := \alpha / (\alpha - 1)$, we can write the same statement as

```{math}
:enumerated: false

g(\alpha) = \EE x
        \quad \text{where} \quad
        \lL(x) = \text{Pareto}(\alpha, 1)
```

The first equation is a version of {eq}`eq-mm`. To estimate $\alpha$ with observations $x_1, \ldots, x_N$ from a $\text{Pareto}(\alpha, 1)$ distribution, we can apply the method of moments, which tells us to solve $g(\hat \alpha) = \frac{1}{N} \sum_{n=1}^N x_n$ for $\hat \alpha$. The result is $\hat \alpha := \bar x_N/(\bar x_N - 1)$.
````


**Generalized method of moments** (GMM) is a small step from method of moments. If we express {eq}`eq-mm` as $\EE[g(\boldtheta) - h(\boldx)] = \boldzero$, then it becomes natural to consider the more general expression

$$
\EE G(\boldtheta, \boldx) = \boldzero
$$ (eq-gmm)

This expression is called the **orthogonality condition**. The generalized method of moments estimator of $\boldtheta$ is the solution $\hat \boldtheta$ to the empirical counterpart, which is

$$
\frac{1}{N} \sum_{n=1}^N G(\hat \boldtheta, \boldx_n) = \boldzero
$$ (eq-gmme)

Of course, there is no guarantee that a solution will exist here, partly because the function $G$ can be nonlinear and partly because the number of equations can be greater than the number of unknowns. If the number of equations is greater, then the estimation problem is said to be **overidentified**.

Our study of overdetermined systems of equations in {ref}`ss-odse` suggests a logical way to handle the overidentified case: Minimize the norm of the left-hand side of {eq}`eq-gmme`. This leads to the expression

$$
\hat \boldtheta = \argmin_{\boldtheta \in \Theta} 
        \left \|
            \frac{1}{N} \sum_{n=1}^N G(\boldtheta, \boldx_n) 
        \right \|
$$ (eq-gmme2)

In practice, we usually make two adjustments to this expression. The first is trivial: we minimize the squared norm, instead of the norm, leaving the minimizer unchanged. The second is to replace the Euclidean norm $\| \cdot \|$ with a weighted norm $\| \cdot \|_W$ defined by $\| \boldx \|_W^2 = \boldx^\T \boldW \boldx$, where $\boldW$ is a positive definite **weighting matrix**. Positive definiteness is desired here because it shares with the Euclidean norm the property that led us to {eq}`eq-gmme2`. That is, $\| \boldx \|_W = 0$ if and only if $\boldx = \boldzero$. The estimation problem is then

$$
\hat \boldtheta 
    = \argmin_{\boldtheta \in \Theta} 
    \left[ \frac{1}{N} \sum_{n=1}^N G(\boldtheta, \boldx_n) \right]^\T
    \hat \boldW
    \left[ \frac{1}{N} \sum_{n=1}^N G(\boldtheta, \boldx_n) \right]
$$ (eq-gmme3)

The weighting matrix has been written as $\hat \boldW$ because it is allowed to depend on the sample. Evidently the choice of $\hat \boldW$ affects the minimizer, and the objective in choosing this matrix is to produce an estimator that has small variance asymptotically.

````{prf:example}

GMM is often used to estimate and test asset pricing models. Many asset pricing models lead to equations of the form

```{math}
:enumerated: false

\boldp_t =
        \EE \left[ M_{t+1} \boldx_{t+1} \given \gG_t \right]
```

where $\boldx_{t+1}$ is a vector of payoffs on $K$ assets at $t+1$, $\boldp_t$ is a corresponding vector of time $t$ asset prices, $M_{t+1}$ is a stochastic discount factor and $\gG_t$ is the time $t$ information set. See, for example, {cite}`hansen2014nobel`. If $\boldZ_t$ is a conformable matrix of observable variables adapted to the filtration $\gG_t$, then we can postmultiply by $\boldZ_t$, pass through the conditional expectation and rearrange to obtain

```{math}
:enumerated: false

\EE \left[ 
                M_{t+1} \boldx_{t+1} \boldZ_t 
                - \boldp_t \boldZ_t \given \gG_t 
            \right] = \boldzero
```

Taking the unconditional expectation and using the law of iterated expectations gives an expression in the form of {eq}`eq-gmm`. In {cite}`hansen1982generalized`, the variables in $\boldZ_t$ include lagged values of asset returns and aggregate consumption growth.
````


(ss-be)=
### Bayesian Estimation
Bayesian inference adopts a rather different strategy from the methods discussed so far. The main idea is to treat parameters as unknown quantities for which we hold subjective beliefs regarding their values. These subjective beliefs are called **priors**. The Bayesian approach to estimation suggests that we take both data and prior knowledge into account when forming an estimate or prediction.

```{prf:example}

Consider the expression "when you hear hooves, think horses, not zebras.\" Prior knowledge should be given some weighting when assessing evidence.
```


For the purposes of estimation, a prior can be thought of as a distribution over $\pP$, the set of distributions in play. The standard Bayesian approach is parametric, so we can specialize this further to a density over parameter space. Thus the primitives in our analysis are:

- $\boldtheta$, the parameter vector, which takes values in $\Theta \subset \RR^J$,

- $\pi$, the **prior distribution**, a density over $\Theta$,

- $\boldx$, the data, and

- $p(\cdot \given \boldtheta)$, the joint density of the data given $\boldtheta$.

Note that $L(\boldtheta) := p(\boldx \given \boldtheta)$ is the likelihood function.

Priors are reassessed based on evidence in the data. This process leads to an updated density over parameter space called the **posterior distribution**, which we represent by $\pi(\boldtheta \given \boldx)$. The posterior is obtained via an application of Bayes' law, which leads us to

$$
\pi(\boldtheta \given \boldx) 
    = \frac{p(\boldx \given \boldtheta) \pi(\boldtheta)}
    {p(\boldx)}
    = \frac{p(\boldx \given \boldtheta) \pi(\boldtheta)}
    {\int p(\boldx \given \boldtheta') \pi(\boldtheta') \diff \boldtheta'}
$$ (eq-brstat)

Here $p(\boldx)$ represents the unconditional density of $\boldx$ evaluated at the outcome. The term on the far right shows why this unconditional density is not listed as a primitive: we can recover it from the other primitives using the law of total probability.

The same method can be applied when the densities are replaced with probability mass functions and the integral is replaced with a sum. Here's a standard example that mixes a density over priors with a binomial PMF for the likelihood:

````{prf:example}
:label: eg-bandit

Consider a one-armed bandit (slot machine) with binary response $v$ indicating a fixed payout ($v=1$) or nothing ($v=0$). We would like to know the probability $\theta$ of $v=1$. Let $v_1, \ldots, v_N$ be a sequence of independent outcomes and let $x := \sum_{n=1}^N v_n$ be the total number of payouts. Recalling [example %s](#eg-binpmf), the likelihood for $x$ conditional on $\theta$ is

```{math}
:enumerated: false

p(x \given \theta) 
        = {N \choose x} \theta^{x} (1-\theta)^{N-x}
```

For our prior we take a Beta$(\alpha, \beta)$ distribution, so

$$
\pi(\theta) 
        = \frac{\theta^{\alpha-1} (1-\theta)^{\beta-1}}{B(\alpha, \beta)}
$$ (eq-betaprior)
for $0 < \theta < 1$. (Parameters like $\alpha$ and $\beta$ that are used to define the prior density are called **hyperparameters**.) Applying {eq}`eq-brstat` gives

$$
\pi(\theta \given x)
        = \frac{ \theta^{x + \alpha - 1} (1 - \theta)^{N - x + \beta - 1}}
                {c(x)}
$$ (eq-betapost)
where $c(x) := p(x) B(\alpha, \beta) / {N \choose x}$. We could try to calculate $c(x)$ directly but there's an easier way. We know that {eq}`eq-betapost` is a density in $\theta$ given $x$. Hence $c(x)$ must be the normalizing constant at $x$. Moreover, in comparing {eq}`eq-betaprior` with {eq}`eq-betapost`, it's clear that $\pi(\theta \given x)$ is a beta density. This leads us to the full form of the posterior, which is

$$
\pi(\theta \given x)
        = \frac{ \theta^{\alpha + x - 1} (1 - \theta)^{N - x + \beta - 1}}
        {B(x + \alpha, N - x + \beta)}
$$ (eq-betapost2)
````


{numref}`f-beta_bayes` shows evolution of the posterior density {eq}`eq-betapost2` in a simulation. The prior is set to Beta$(3, 5)$. The true payout probability is $\theta_0 = 0.7$. Despite the poor prior, the data shift probability mass towards $\theta_0$. See `johnstachurski.net/emet` for code.

```{figure} figures/beta_bayes.svg
:name: f-beta_bayes

Evolution of the posterior from Beta$(3, 5)$ prior
```

Point estimates are extracted from the posterior distribution based on some measure of central tendency, such as the mean, the median, or the **mode** of the posterior (i.e., the maximizer in the unimodal case).

````{prf:example}

The mean of the posterior in {eq}`eq-betapost2` yields the estimator

```{math}
:enumerated: false

\hat \theta := \frac{\alpha + x}{\alpha + \beta + N}
```

More payouts shift our estimator upwards. In the limit, $\hat \theta$ is near $\frac{x}{N}$, which is the MLE of $\theta$. This illustrates a common theme: The difference between maximum likelihood and Bayesian estimates typically concerns finite sample properties.
````


The posterior we came up with in {eq}`eq-betapost2` belongs to the same parametric class as the prior. Priors where the parametric class is preserved under Bayesian updating for a specific likelihood function are called **conjugate**. Conjugacy makes application of Bayes' law particularly simple.

In applications, the existence of a closed-form solution for the posterior is rare, and the integration over the parameter space in {eq}`eq-brstat` has to be carried out numerically. The standard technique is Markov chain Monte Carlo. An exposition of the Metropolis–Hastings algorithm was given in {ref}`ss-mcmc`. If you return to {eq}`eq-alpha`, you will notice how implementation of this algorithm only requires us to evaluate ratios of of the posterior, which means that the integral term on the right-hand side of {eq}`eq-brstat` cancels out. Because the integral can be very high dimensional, this often eliminates a huge amount of complexity.

Bayesian inference has experienced a surge of popularity in recent years. One reason is that in high dimensions, exploring the posterior via MCMC has proved to be more successful in practice than the numerical optimization required to obtain maximum likelihood estimates. Another is that Bayesian estimation provides a form of regularization that stabilizes and typically improves estimation of complex models. See {ref}`ss-abp`. A third is that Bayesian estimation comes with an elegant, unified decision-theoretic approach to inference. Some discussion is given in {ref}`ss-risk`.

## Further Reading

The ERM principle is a very general principle for solving statistical problems and producing estimators. For such a general method it is difficult to give a set of strong results showing that ERM produces good estimators. Indeed, there will be instances when ERM produces poor estimators, as discussed in {ref}`ss-ermls`. Having said that, some rather general consistency results have been obtained. Additional discussion can be found in {cite}`vapnik2013nature`.

There are many good treatments of Bayesian estimation available in the literature, including {cite}`geweke2005contemporary`, {cite}`geweke2011oxford`, and {cite}`kroese2014statistical`.

{cite}`chernozhukov2011inference` discuss quantile regression of the $\tau$th quantile in cases where $\tau$ is close to zero or one.

## Exercises

```{exercise}
:label: ex-minqrisk

Let $x$ be a random variable with $\mu := \EE[x]$ and finite second moment. Consider the risk function given by $R(\theta) = \EE[ (\theta  - x)^2 ]$. Show that $\mu$ is the minimizer of $R(\theta)$ over all $\theta \in \RR$.
```


```{exercise}
:label: ex-solse

Confirm the solutions in {eq}`eq-solse` by differentiating {eq}`eq-lstmpe`.
```


```{exercise}
:label: ex-ch08_estimators-auto-1

Show that the method of moments estimator $\hat \alpha$ in [example %s](#eg-estpareto) converges in probability to $\alpha$ whenever $\alpha > 1$. (Later we will call this property consistency.)
```


```{exercise}
:label: ex-maxjn

Confirm that the maximizer of {eq}`eq-mnd` is the sample mean of $x_1, \ldots, x_N$.
```


```{exercise}
:label: ex-ch08_estimators-auto-2

Confirm the results in {eq}`eq-mlen` by differentiating {eq}`eq-mnds`.
```


```{exercise}
:label: ex-quantreg

Let $F$ be a strictly increasing CDF on $\RR$ and let $\tau \in (0, 1)$ be given. Let $y$ be a random variable with $\lL(y) = F$. Adopting the notation of {ref}`sss-qreg`, show that the solution to the optimization problem $\min_{\xi \in \RR} \, \EE L_\tau (y, \xi)$ is the $\xi$ that solves $F(\xi) = \tau$.
```


```{exercise}
:label: ex-sm

Consider the one-armed bandit problem in [example %s](#eg-bandit). Directly obtain the log likelihood function of $\theta$ given $v_1, \ldots, v_N$. Use the principle of maximum likelihood to show that the sample mean $\bar v_N$ is the MLE for $\theta$.
```


````{exercise}
:label: ex-likerat

Let $f$ and $g$ be two fixed densities, let $x_1,\ldots,x_N$ be IID and consider the **likelihood ratio statistic**

```{math}
:enumerated: false

y_n := \prod_{i=1}^n \frac{g(x_i)}{f(x_i)}
```

Show that if $\lL(x_n) = f$ for all $n$, then $\{y_n\}$ is a martingale with respect to the filtration $\{\fF_n\}$ defined by $\fF_n := \{x_1,\ldots,x_n\}$.
````


````{exercise}
:label: ex-normkl

Let $D(p_1, p_2)$ be the KL deviation between normal densities $p_1$ and $p_2$ where $p_i = \nN(\mu_i, \sigma^2_i)$ for $i=1,2$. Show that

```{math}
:enumerated: false

D(p_1, p_2) 
        = \ln \frac{\sigma_2}{\sigma_1}
        + \frac{\sigma_1^2 + (\mu_1 - \mu_2)^2}{2 \sigma_2^2}
        - \frac{1}{2}
```
````


### Solutions to Selected Exercises

````{solution} ex-minqrisk

Adding and subtracting $\mu$, we can express $R(\theta)$ as

```{math}
:enumerated: false

R(\theta) = \EE \{ [ (\theta - \mu) + (\mu - x) ]^2 \}
```

Expanding this out and using $\EE[x] = \mu$, we obtain $R(\theta) = (\theta - \mu)^2 + \var x$. Evidently a minimum is obtained when $\theta = \mu$.
````


````{solution} ex-quantreg

We can write the minimization problem as

```{math}
:enumerated: false

\min_{\xi \in \RR}
        \left\{
            (\tau - 1)
            \int_{-\infty}^{\xi} (t - \xi) \, F(\diff t)
            + 
            \tau
            \int_\xi^{\infty} (t - \xi) \, F(\diff t)
        \right\}
```

The first-order condition is $(1 - \tau) F(\xi) + \tau (F(\xi) - 1) = 0$. Simplifying gives $F(\xi)=\tau$, as was to be shown.
````


````{solution} ex-sm

Each $v_n$ is binary, with PMF given by $p(s; \theta) := \theta^s (1 - \theta)^{1 - s}$ for $s \in \{0,1\}$. By independence, the joint distribution is the product of the marginals, and hence the log likelihood is

```{math}
:enumerated: false

\ell(\theta) 
        = \sum_{n=1}^N \log p(v_n; \theta)
        = \sum_{n=1}^N \left[ v_n \log \theta + (1 - v_n) \log (1 - \theta) \right]
```

Differentiating with respect to $\theta$ and setting the result equal to zero yields $\hat \theta = \bar v_N$ as claimed.
````


````{solution} ex-likerat

It is clear that $y_n$ is $\fF_n$-measurable, and hence $\{y_n\}$ is adapted to $\{\fF_n\}$. In addition, we have

```{math}
:enumerated: false

\EE [ y_{n+1} \given \fF_n ]
        = \EE \left[ \prod_{i=1}^{n+1} \frac{g(x_i)}{f(x_i)} \given \fF_n \right]
        = \prod_{i=1}^n \frac{g(x_i)}{f(x_i)}
            \EE \left[ \frac{g(x_{n+1})}{f(x_{n+1})} \given \fF_n \right]
```

Since $\{x_n\}$ is IID, $\lL(x_{n+1}) = f$ and $g$ is a density, we have

```{math}
:enumerated: false

\EE \left[ \frac{g(x_{n+1})}{f(x_{n+1})} \given \fF_n \right]
        = \EE \left[ \frac{g(x_{n+1})}{f(x_{n+1})} \right]
        = \int  \frac{g(s)}{f(s)} f(s) \, \diff s
        = \int  g(s) \, \diff s 
        = 1
```

```{math}
:enumerated: false

\fore
        \EE [ y_{n+1} \given \fF_n ]
        = \prod_{i=1}^n \frac{g(x_i)}{f(x_i)}
        = y_n
        
```
````


[^1]: Why not always try to infer $P$ itself? Although it's true that when we know $P$ we can in principle recover any feature $\gamma(P)$, we should also bear in mind the following general principle of inference with limited information: In solving a given problem, try to avoid first solving a more general problem as an intermediate step. Distributions are much more complicated objects than real numbers, so if we care only about the median of a distribution, say, it might be best to try to discover this single value directly rather than trying to infer the entire distribution first.

[^2]: See, for example, §2.4 of {cite}`friedman2001elements`.

[^3]: The risk $R(\hat f_d)$ is evaluated by substituting $\hat f_d$ into the expression for $R$ in {eq}`eq-rminpoly`. You can find the code at `johnstachurski.net/emet`.

[^4]: The scare quotes on "unlikely\" are just a reminder that our distribution is absolutely continuous and hence all individual outcomes $s \in \RR$ have probability zero. When we say that an outcome $x$ is unlikely, we mean that there is little probability mass in the neighborhood of $x$.

[^5]: More precisely, $D(q, p) = 0$ if and only if $p = q$ almost everywhere. Densities are equal almost everywhere when the set of points at which they fail to be equal has Lebesgue measure zero.
