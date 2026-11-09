---
title: "Confidence Intervals and Tests"
label: c-ciht
---

Until now we've focused primarily on point estimation. Since data are not infinite, point estimation involves uncertainty. We should not be overconfident in our claim to know the data-generating process. Rather, we should be only as confident as the data allows us to be.

We discussed uncertainty previously when we considered sampling distributions. In this chapter we further quantify estimation uncertainty using the language of confidence sets. We also cover hypothesis tests, which are a second key method of inference in classical statistics.

To tie in with standard presentations we'll switch to parametric language. The notation is as except that

(i) the marginal distribution $P$ of each observation $\boldz_n$ is written as $P_{\theta}$,

(ii) the universe $\pP$ of possible distributions is $\setntn{P_{\theta}}{\theta \in \Theta}$, and

(iii) the joint distribution of the data $\zdata = (\boldz_1, \ldots, \boldz_N)$ is written as $\Pdatatheta$.

However, we won't actually restrict ourselves to the case of finitely many parameters. Instead, $\theta$ can be thought of more generally as an index over the set of distributions.

(s-cons)=
## Confidence Sets
Confidence sets are random subsets of model space that contain the true data generating process with high probability. Let's start with the key definitions and then move on to examples.

(ss-pe)=
### Finite Sample Confidence Sets
Fix $\alpha \in (0,1)$. A random set $C(\zdata) \subset \Theta$ is called a **$1 - \alpha$ confidence set** for $\theta$ if $C(\zdata)$ is observable given the data and, for any $\theta \in \Theta$,

$$
\PP \{ \theta  \in C(\zdata) \} \geq 1 - \alpha
    \quad \text{whenever} \quad
    \lL(\zdata) = \Pdatatheta
$$ (eq-proci)

Note that it's the set that's random here, not $\theta$.

Confidence sets that are intervals are usually called **confidence intervals**.

Suppose, for example, that we want to estimate the mean of an unknown distribution based on a set of observations $\zdata$. If we can come up with a 95% confidence set $C(\zdata)$ in the sense of {eq}`eq-proci`, then we know that the set will contain the mean in about 95% of our experiments, regardless of the underlying distribution.

The last part of the sentence is essential: We don't know what $\theta$ is, and so the set *must be designed* such that the probability of containing $\theta$ exceeds $1 - \alpha$ regardless of which $\theta$ is generating the data. The difficulty is to actually construct a confidence set with this property for a reasonably large class of distributions. Here's a traditional approach that uses strong parametric assumptions.

````{prf:example}
:label: eg-cisk

Let $\zdata = \boldx = (x_1,\ldots,x_N)$ where each $x_n$ is an independent draw from $\nN(\mu, \sigma^2)$. Suppose for now that only $\mu \in \RR$ is unknown. We wish to form a confidence interval for $\mu$. By {eq}`eq-sdsmnc`, we have

$$
\lL \left[
            \sqrt N \frac{(\bar x_N - \mu)}{\sigma}  
        \right]
        = \nN(0, 1)
$$ (eq-cisk)
Notice that this statement is true regardless of the values of $\mu$ and $\sigma$. Applying {eq}`eq-symcdf` now gives

```{math}
:enumerated: false

\PP
            \left\{ 
            \frac{\sqrt{N}}{\sigma} | \bar x_N - \mu | 
            \leq z_{\alpha/2}
            \right\}
        = 1 - \alpha
        \quad \text{when} \quad
        z_{\alpha/2} := \Phi^{-1} \left(1 - \frac{\alpha}{2} \right)
```

Here $\Phi$ is the standard normal CDF. Some rearranging gives

```{math}
:enumerated: false

\PP
            \left\{ 
            \bar x_N - 
                \frac{\sigma}{\sqrt{N}} z_{\alpha/2}
            \leq \mu \leq 
            \bar x_N + 
                \frac{\sigma}{\sqrt{N}} z_{\alpha/2}
            \right\}
        = 1 - \alpha
```

Since this argument is true regardless of the value of $\mu$, we conclude that the interval

```{math}
:enumerated: false

C(\boldx) 
        := (\bar x_N - e_n, \, \bar x_N + e_n)
        \quad \text{with} \quad
        e_n := \frac{\sigma}{\sqrt{N}} z_{\alpha/2}
```

is a $1 - \alpha$ confidence interval for $\mu$. (Note that $C(\boldx)$ is indeed observable under our assumption that $\sigma$ is known.)
````


In [example %s](#eg-cisk), the assumption that $\sigma$ is known is clearly unrealistic. If $\sigma$ is unknown, a natural approach is to replace this term with the sample standard deviation $s_N$. In doing so, we will use the following fact, which [exercise %s](#ex-ctdr) asks you to prove.

```{prf:property}
:label: fa-ctdr

If $x_1,\ldots,x_N$ are IID draws from $\nN(\mu, \sigma^2)$, then

$$
\sqrt{N - 1} \frac{(\bar x_N - \mu)}{s_N}
$$ (eq-cisk2)
has a Student's $t$-distribution with $N-1$ degrees of freedom.
```


The random variable {eq}`eq-cisk2` is called a **pivot**, which means that its distribution does not depend on the unknown parameters.

````{prf:example}
:label: eg-cisk2

Consider then the setting of [example %s](#eg-cisk) but with $\sigma$ also unknown. Let $F_{N-1}$ be the CDF of the $t$-distribution with $N-1$ degrees of freedom. Using the same reasoning as [example %s](#eg-cisk) with $z_{\alpha/2}$ replaced by $t_{\alpha/2} := F^{-1}_{N-1}(1 - \alpha/2)$, we obtain

```{math}
:enumerated: false

\PP 
            \left\{ 
            \bar x_N - 
                \frac{s_N}{\sqrt{N-1}} t_{\alpha/2}
            \leq \mu \leq 
            \bar x_N + 
                \frac{s_N}{\sqrt{N-1}} t_{\alpha/2}
            \right\}
        = 1 - \alpha
```
````


Recall that the standard deviation of $\bar x_N$ is $\sigma/ \sqrt{N}$. The term $s_N/\sqrt{N-1}$ is a sample estimate of $\sigma/ \sqrt{N}$. Following the discussion in {ref}`sss-se`, we call $s_N/\sqrt{N-1}$ the standard error of $\bar x_N$ and write it as $\se(\bar x_N)$. We can then write the confidence interval in [example %s](#eg-cisk2) as

$$
C(\boldx) := (\bar x_N - \se(\bar x_N) t_{\alpha/2},\;
         \bar x_N + \se(\bar x_N) t_{\alpha/2})
$$ (eq-citc)

### Asymptotic Methods

While strong parametric assumptions of the form used in [examples %s](#eg-cisk) and [%s](#eg-cisk2) are difficult to justify in most econometric scenarios, there is another approach to forming confidence intervals that requires much less in terms of structure imposed on the underlying distributions, provided that the sample size is relatively large. Let's start with a definition: given $\alpha \in (0,1)$, a set $C_N(\zdata) \subset \Theta$ is called an **asymptotic $1 - \alpha$ confidence set** for $\theta$ if $C_N(\zdata)$ is observable given the data and

$$
\lim_{N \to \infty} \PP \{ \theta \in C_N(\zdata)  \} \geq 1 - \alpha
$$ (eq-procia)

for all $\theta \in \Theta$. Note that $C_N(\zdata)$ is really a sequence of sets and the definition concerns this sequence.

How can we create asymptotic confidence sets? Let's consider first the case of scalar $\theta$. Suppose that we have an estimator $\hat \theta_N$ of $\theta$ that is asymptotically normal for all $\theta \in \Theta$. Specializing the results in {ref}`sss-ad` to the scalar case, this means that, for all $\theta \in \Theta$, there exists a positive constant $v(\theta)$ such that

$$
\sqrt{N} (\hat \theta_N - \theta) \tod \nN(0, v(\theta)) 
        \quad \text{as} \quad
        N \to \infty
$$ (eq-asnorm2)

The constant $v(\theta)$ is called the asymptotic variance of $\hat \theta_N$. Suppose we have a sequence of statistics $\se(\hat \theta_N)$ such that

$$
\sqrt{N} \se(\hat \theta_N)
        \toprob \sqrt{v(\theta)} 
        \quad \text{as} \quad
        N \to \infty
$$ (eq-secfs)

As [exercise %s](#ex-nbsez) asks you to show, {eq}`eq-asnorm2` and {eq}`eq-secfs` imply that

$$
\frac{\hat \theta_N - \theta}{\se(\hat \theta_N)} 
    \, \tod \, \nN(0, 1)
    \quad \text{as} \quad
    N \to \infty
$$ (eq-nbsez)

We can now create our asymptotic confidence intervals as

$$
C_N := 
    (\hat \theta_N - \se(\hat \theta_N) z_{\alpha/2}, 
    \, \hat \theta_N + \se(\hat \theta_N) z_{\alpha/2})
$$ (eq-sacine)

To see this, take the limit and apply {eq}`eq-nbsez` to get

```{math}
:enumerated: false

\begin{aligned}
\lim_{N \to \infty} \PP \{ \theta \in C_N(\zdata) \} 
    & = \lim_{N \to \infty} 
        \PP
            \left\{ 
            \hat \theta_N - 
                \se(\hat \theta_N) z_{\alpha/2}
            \leq \theta \leq 
            \hat \theta_N + 
                \se(\hat \theta_N) z_{\alpha/2}
            \right\} \\
    & = \lim_{N \to \infty} 
        \PP 
            \left\{ 
                -z_{\alpha/2}
                \leq \frac{\hat \theta_N - \theta}{\se(\hat \theta_N)}
                \leq z_{\alpha/2}
            \right\} \\
    & = 1 - \alpha
\end{aligned}
```

````{prf:example}
:label: eg-ldss

Let $\bar x_N$ be the sample mean of IID data $\{x_n\}$. Suppose that these data come from some distribution with finite second moment. Let $\mu$ denote the common mean and let $\sigma^2$ denote the variance. Let $s_N$ be the sample standard deviation. Combining the central limit theorem and consistency of $s_N$ for $\sigma$, we have

$$
\sqrt{N}(\bar x_N - \mu) \tod \nN(0, \sigma^2)
        \;  \text{ and } \;
        \sqrt{N}\se(\bar x_N) \toprob \sigma
        \quad  \text{for} \quad
        \se(\bar x_N) := \frac{s_N}{\sqrt{N}}
$$ (eq-smnse)
It follows that, with this definition of $\se(\bar x_N)$, the set

```{math}
:enumerated: false

(\bar x_N - \se(\bar x_N) z_{\alpha/2}, 
            \, \bar x_N + \se(\bar x_N) z_{\alpha/2})
```

is an asymptotic $1-\alpha$ confidence interval for $\bar x_N$.
````


(ss-ane)=
### A Nonparametric Example
Here's an example of a nonparametric confidence set. In {ref}`ss-eds` we learned that if $\boldx = (x_1, \ldots, x_N)$ is a vector of IID draws from some CDF $F$ and $\hat F_N$ is the corresponding ECDF, then $\| F - \hat F_N \|_{\infty}$ converges in probability to $0$. In 1933, A. N. Kolmogorov used an extension of the central limit theorem to obtain an asymptotic distribution for this term. In particular, he showed that when $F$ is continuous,

$$
\sqrt{N} \sup_{s \in \RR} |\hat F_N(s) - F(s)| \tod K
$$ (eq-kr)

where $K$ is the **Kolmogorov** CDF

$$
K(s) := \frac{\sqrt{2 \pi}}{s} \sum_{i=1}^{\infty}
    \exp \left[ - \frac{(2 i - 1)^2 \pi^2}{8 s^2} \right]
    \qquad (s \geq 0)
$$ (eq-kdf)

As in the CLT, the limiting distribution $K$ is independent of the CDF $F$ that generates the data.

We can use {eq}`eq-kr` to produce an asymptotic $1-\alpha$ confidence set for $F$. To do so, let $\mathfrak F$ be the set of all CDFs on $\RR$, let $x_1,\ldots,x_N$ be IID draws from $F \in \mathfrak F$, let $k_{1-\alpha} := K^{-1}(1 - \alpha)$, and let

```{math}
:enumerated: false

C_N(\boldx) := \left\{ G \in \mathfrak F
        \; : \; \hat F_N(s) - \frac{k_{1-\alpha}}{\sqrt{N}} \leq G(s) \leq 
        \hat F_N(s) + \frac{k_{1-\alpha}}{\sqrt{N}} \; \text{ for all } \; s \in \RR \right\}
```

The set $C_N(\boldx) \subset \mathfrak F$ is an asymptotic $1-\alpha$ confidence set for $F$. Indeed, after rearranging the expression, we get

```{math}
:enumerated: false

\begin{aligned}
F \in C_N(\boldx) 
    & \iff  - k_{1-\alpha} \leq \sqrt{N}(\hat F_N(s) - F(s)) \leq k_{1-\alpha}
    \; \text{ for all } \; s  \\
    & \iff  \sqrt{N}|\hat F_N(s) - F(s)| \leq k_{1-\alpha}
    \; \text{ for all } \; s  \\
    & \iff  \sup_s \sqrt{N}|\hat F_N(s) - F(s)| \leq k_{1-\alpha}
\end{aligned}
```

Applying {eq}`eq-kr` now confirms our claim:

```{math}
:enumerated: false

\lim_{N \to \infty} \PP \{ F \in C_N(\boldx) \}
     = \lim_{N \to \infty} \PP 
     \left\{ \sup_s \sqrt{N}|\hat F_N(s) - F(s)| \leq k_{1-\alpha} \right\}
     = 1 - \alpha
```

Given our data $x_1,\ldots,x_N$ and the corresponding ECDF $\hat F_N$, we can present the confidence set $C_N(\boldx)$ visually by highlighting the area between the lower bound $\hat F_N(s) - k_{1-\alpha}/\sqrt{N}$ and the upper bound $\hat F_N(s) + k_{1-\alpha}/\sqrt{N}$ over $s$. This is done in [figures %s](#f-ks_sim1) and [%s](#f-ks_sim2). The data are drawn from a $t$-distribution $F$ with 10 degrees of freedom and $\alpha$ is set to $0.05$. The values of $N$ are 250 for [figure %s](#f-ks_sim1) and 1,000 for [figure %s](#f-ks_sim2).

```{figure} figures/ks_sim1.svg
:name: f-ks_sim1

Confidence set for the ECDF with 250 observations
```

```{figure} figures/ks_sim2.svg
:name: f-ks_sim2

Confidence set for the ECDF with 1,000 observations
```

(s-hyp)=
## Hypothesis Tests
This section considers problems of inference where we

(i) hold a belief or theory concerning the probabilities generating the data and

(ii) consider whether the data provides evidence for or against that theory.

```{prf:example}

While asset price returns data tend to have heavier tails than the normal distribution, for certain assets normality might still be a reasonable and convenient approximation. Let's look at whether or not normality provides a reasonable fit in the case of daily returns on the Nikkei 225. {numref}`f-nikkei_ecdf` shows a comparison of the empirical distribution of standardized returns from January 2014 until July 2015 with the CDF of the standard normal distribution.[^1] While the ECDF of the data and the CDF are not identical, we wouldn't expect this even if the hypothesis of normality is completely true. After all, we are only observing a sample, and a sample is subject to randomness. So the question becomes, is the sample "unlikely\" given the hypothesis of normality?
```


```{figure} figures/nikkei_ecdf.svg
:name: f-nikkei_ecdf

ECDF for standardized returns versus the standard normal CDF
```

To take the next steps, we must formalize the process of setting up a hypothesis and quantifying the extent to which the data points in this or some other direction.

Hypothesis testing begins with the specification of a **null hypothesis**, which is a statement that the observed data are being generated by a certain model, or one of a certain class of models. One might imagine that the standard procedure in statistics is to show the validity of the null hypothesis but it is not. Rather, a hypothesis test is an attempt to reject the null.

To see why, let's recall a famous example due to Karl Popper (1902–1994). Consider the claim that all swans are white. It is futile to attempt to prove correctness of this hypothesis: No quantity of white swan sightings can confirm the claim that *all* swans on planet earth are white. Yet a single black swan can show the claim to be false. In this sense, attempting to falsify a theory on the basis of observation is more constructive than attempting to confirm it.

Although we attempt to reject the null hypothesis, we do so only if strong evidence against it is observed. The rationale is this: Suppose that we have a collection of theories about how the economy works. The procedure would then be to step through the theories, at each stage taking correctness of the theory as the null hypothesis and attempting to reject. If the theory is rejected, then we can discard it. This is a useful process of elimination. However, we don't want to mistakenly discard a good theory. Hence we only reject when there is strong evidence against the null hypothesis.

### Constructing Tests

We adopt the setting and notation used in {ref}`ss-pe`. A null hypothesis is a specification of a set of models that we believe the data-generating process belongs to. This amounts to specifying a subset $\Theta_0$ of $\Theta$. The null hypothesis is often written as

```{math}
:enumerated: false

H_0 \colon \theta \in \Theta_0
```

If $\Theta_0$ is a singleton, then the null hypothesis is called a **simple hypothesis**. If not, then the null hypothesis is called a **composite hypothesis**. A test of the null hypothesis is a test of whether or not the observed data were generated by $P_{\theta}$ for some $\theta \in \Theta_0$.

```{prf:example}

The hypothesis of purchasing power parity is tested by considering models such as

$$
p = a + \beta e p^* + \sigma u
$$ (eq-ppp)
where $p$ is a price or index of domestic prices, $e$ is an exchange rate, $p^*$ is a corresponding foreign price, $u$ is a disturbance term and $a, \beta, \sigma$ are parameters. The null hypothesis for absolute purchasing power parity is that {eq}`eq-ppp` holds with $a = 0$ and $\beta=1$. Since $\sigma$ is not pinned down by the null, this is a composite null hypothesis.
```


Formally, a **test** of $H_0$ is a binary function $\phi$ mapping the observed data $\zdata$ into $\{0,1\}$. The decision rule is

```{math}
:enumerated: false

\begin{aligned}
\text{if } \phi(\zdata) & = 1, \text{ then reject } H_0   \\
    \text{if } \phi(\zdata) & = 0, \text{ then do not reject } H_0
\end{aligned}
```

(Failing to reject $H_0$ should not be confused with accepting $H_0$. More on this below.) Our aim will be to design $\phi$ such that it takes the value $1$ when the data present strong evidence against $H_0$.

```{prf:example}

Let $\zdata$ be a set of pairs $\boldz_n = (x_n, y_n)$ from some unknown bivariate distribution $P_{\theta}$. If $H_0$ is the hypothesis that the correlation between $x$ and $y$ is negative, then a large positive sample correlation would constitute evidence against $H_0$. Hence our test might take the form $\phi(\zdata) = \1\{\hat \rho > c\}$, where $\hat \rho$ is the sample correlation. An appropriate value of $c$ remains to be determined.
```


The outcome of our test depends on the random sample $\zdata$, and being random, its realization can be misleading. There are two different ways in which the realization can mislead us. First, we can mistakenly reject the null hypothesis when it is in fact true. This is called **type I error**. Second, we can fail to reject the null hypothesis when it is false. This is called **type II error**.

The **power function** associated with the test $\phi$ is the function

```{math}
:enumerated: false

\beta(\theta) := \PP_{\theta} \{\phi(\zdata) = 1\}
    \qquad (\theta \in \Theta)
```

Here and below, the symbol $\PP_{\theta}$ indicates that we are computing probabilities under the assumption that $\lL(\zdata) = \Pdatatheta$. Thus, $\beta(\theta)$ is the probability that the test rejects when the data are generated by the probabilistic model identified by $\theta$.

Ideally, we would like $\beta(\theta) = 0$ when $\theta \in \Theta_0$, and $\beta(\theta) = 1$ when $\theta \notin \Theta_0$. In practice, this is difficult to achieve.

As discussed above, we tend to be conservative in rejecting the null because we don't want to discard good theories. For this reason it is traditional to keep the probability of type I error small. Then, if our test tells us to reject the null, it's unlikely the null is true. Because of this, the standard procedure is to choose a small number $\alpha$ such as 0.05 or 0.01, and then adjust the test such that

$$
\beta(\theta) \leq \alpha
    \quad \text{for all } \theta \in \Theta_0
$$ (eq-defsize)

If {eq}`eq-defsize` holds, then the test is said to be of **size $\alpha$**. Similarly the **size** of any test with power function $\beta$ is

$$
\alpha := \sup_{\theta \in \Theta_0} \beta(\theta)
$$ (eq-defsize2)

This is the maximal rejection probability when the null hypothesis is true.

In constructing tests, a common (but not universal) setup is to define a real-valued **test statistic** $T$ and a **critical value** $c$, and then set

$$
\phi(\zdata) := \1\{T(\zdata) > c\}
$$ (eq-tfftc)

The pair $(T, c)$ then defines the test, and the rule becomes:

```{math}
:enumerated: false

\text{reject } H_0 \text{ if and only if } T(\zdata) > c
```

````{prf:example}
:label: eg-normslz

Suppose that $x_1,\ldots,x_N$ are independent draws from $\nN(\mu, 1)$ where the value of $\mu$ is unknown. The null hypothesis is $\mu \leq 0$, or $\Theta_0 = (-\infty, 0]$. Since we want to make inference about the mean, a natural choice for our test statistic is the sample mean. Thus

```{math}
:enumerated: false

T(\zdata) :=: T(x_1,\ldots,x_N) := \bar x_N
```

Each $c \in \RR$ now gives us a test via {eq}`eq-tfftc`, with power function $\beta(\mu) = \PP_{\mu} \{ \bar{x}_N > c\}$. To evaluate $\beta$, recall that $\lL(\bar{x}_N) = \nN(\mu, 1/N)$. As a result, if $\Phi$ is the CDF of the standard normal distribution and $\lL(z) = \Phi$, then

```{math}
:enumerated: false

\PP_{\mu} \{ \bar{x}_N \leq c \} 
        = \PP\{ \mu + N^{-1/2} z  \leq c \} 
        = \PP\{ z \leq N^{1/2}(c - \mu) \} 
        = \Phi[N^{1/2}(c - \mu)]
```

$$
\fore
        \beta(\mu) = 1 - \Phi[N^{1/2}(c - \mu)]
$$ (eq-pfct)
Given $c$, the power function is increasing in $\mu$ because higher $\mu$ pushes up the mean of $\bar{x}_N$, making the event $\{ \bar{x}_N > c\}$ more likely. Given $\mu$, the function is decreasing in $c$, because higher $c$ makes the event $\{ \bar{x}_N > c\}$ less likely.
````


Plots of the power function $\beta$ in {eq}`eq-pfct` are presented in [figure %s](#f-power) for two different values of $c$. Here $N$ is fixed at 10. Since $\Theta_0 = (-\infty, 0]$, the size of the test in each case is $\sup_{\mu \leq 0} \beta(\mu)$. Since the power curves are increasing, this is just $\beta(0)$.

{numref}`f-power` illustrates a typical trade-off between type I and type II error. If we increase $c$, we make rejection less likely for all values of $\mu$. This pushes down type I error, but it also increases the probability that we fail to reject a false null.

```{figure} figures/power.svg
:name: f-power

Power function $\beta$
```

### Choosing Critical Values

Let's think a bit more about the test in {eq}`eq-tfftc`. Typically, the choice of $T$ is suggested by the problem. For example, if our hypothesis is a statement about the second moment of a random variable, then we might take $T$ to be the sample second moment. Once $T$ is fixed, we need to adjust $c$ such that $(T,c)$ attains the appropriate size. Thus the standard procedure is to

(i) choose a desired size $\alpha$ according to our tolerance for type I error,

(ii) identify a suitable test statistic $T$, and

(iii) choose a critical value $c$ so that $(T,c)$ is of size $\alpha$.

In performing the last step, we can balance our desire to minimize type II error while maintaining a size $\alpha$ test by choosing $c$ such that {eq}`eq-defsize` holds with equality. In this case, the problem is to choose $c$ to solve

$$
\alpha = \sup_{\theta \in \Theta_0} \PP_{\theta}\{T(\zdata) > c\}
$$ (eq-solfc)

{numref}`f-c_alpha` gives an illustration. In the figure, we've taken $\Theta_0$ to be the two element set $\{\theta_a, \theta_b\}$. The blue line gives an imaginary distribution of $T(\zdata)$, represented as a density, when $\zdata$ is generated by $\theta_a$. The black line gives the same for $\theta_b$. Assuming that a value of $\alpha$ is prescribed, our next step is to determine the $c$ such that {eq}`eq-solfc` holds. Here we choose $c$ such that the largest of the two shaded areas is equal to $\alpha$.

```{figure} figures/c_alpha.svg
:name: f-c_alpha

Determining the critical value
```

````{prf:example}
:label: eg-normslz2

Let's look again at [example %s](#eg-normslz), where $x_1,\ldots,x_N \iidsim \nN(\mu, 1)$ for $\mu$ unknown, and our null hypothesis is $\mu \leq 0$. Given $\alpha$, our task is to find the appropriate critical value $c$ so that the test $(T,c)$ is of size $\alpha$. To solve for $c$ given $\alpha$, we use {eq}`eq-solfc`. Applying the expression for the power function in {eq}`eq-pfct`, this becomes

```{math}
:enumerated: false

\alpha = \sup_{\mu \leq 0} \{ 1 - \Phi[N^{1/2}(c - \mu)] \}
```

The right-hand side is increasing in $\mu$, so the supremum is obtained by setting $\mu = 0$. Setting $\mu = 0$ and solving for $c$, we obtain

```{math}
:enumerated: false

c(\alpha) := N^{-1/2} \Phi^{-1} (1 - \alpha)
```

where $\Phi^{-1}$ is the quantile function of the standard normal distribution. In R, this function can be evaluated using `qnorm`. For example,

``` r
> alpha = 0.05
> qnorm(1 - alpha)
[1] 1.644854
```

Since $\Phi^{-1}$ is increasing, we see that smaller $\alpha$ corresponds to larger $c(\alpha)$—we reduce the probability of type I error by increasing the critical value that the mean $\bar{x}_N$ must obtain for rejection to occur. Also higher $N$ brings down $c(\alpha)$: more data allows us to reduce the critical value without increasing the probability of rejecting a true null.
````


### Asymptotic Tests

In general, constructing a test of a given size requires some knowledge of the distribution of the test statistic under the null hypothesis. In many cases, however, we know relatively little about the distribution of the test statistic. Fortunately, we still have options. One approach is to bootstrap. Another is to try to determine the asymptotic distribution of the test statistic. As we've learned from the central limit theory, asymptotics hold out the possibility that we can pin down distributions without actually assuming a parametric structure *a priori*.

To proceed, we need to switch to a notion of asymptotic size for tests, rather than finite sample size. Writing $\beta_N$ instead of $\beta$ to emphasize the fact that the power function typically depends on sample size, a test is called **asymptotically of size $\alpha$** if

$$
\lim_{N \to \infty}
    \beta_N(\theta) \leq \alpha
    \quad \text{for all } \theta \in \Theta_0
$$ (eq-asmposa)

```{prf:example}
:label: eg-asymptest

Let $x_1,\ldots,x_N$ be an IID sample with mean $\theta$ and variance $\sigma^2$. Assume that both $\theta$ and $\sigma$ are unknown. We wish to test the hypothesis $H_0 \colon \theta = \theta_0$. Consider the statistic

$$
t_N 
        := \sqrt{N} \left\{ \frac{\bar x_N - \theta_0}{s_N}  \right\}
        =  \frac{\bar x_N - \theta_0}{\se(\bar x_N)}
$$ (eq-easn2)
where $\se(\bar x_N) := \frac{s_N}{\sqrt{N}}$ as in {eq}`eq-smnse`. As discussed in [example %s](#eg-ldss), this expression converges in distribution to a standard normal. Hence

$$
\phi_N(\boldx) := \1 \{ |t_N| > z_{\alpha/2} \}
$$ (eq-easn3)
is asymptotically of size $\alpha$ (see [ex. %s](#ex-nwald)).
```


#### A Test for the Empirical Distribution

As another example of how to construct tests which are asymptotically of size $\alpha$, let's return to the standardized daily returns data presented in [figure %s](#f-nikkei_ecdf). Let $\Phi$ be the standard normal CDF as before, and let $H_0$ be that standardized returns are IID draws from $\Phi$. Let $\alpha$ be given, and let $k_{1-\alpha} = K^{-1}(1-\alpha)$ be the $1-\alpha$ quantile of the Kolmogorov distribution $K$, as defined in {eq}`eq-kdf`. Finally, let $\hat F_N$ be the ECDF of the data. If the null hypothesis is true, then, by {eq}`eq-kr`, we have

$$
\sqrt{N} \sup_{s \in \RR} |\hat F_N(s) - \Phi(s)| \tod K
$$ (eq-kr2)

For the test

```{math}
:enumerated: false

\phi_N(\boldx) 
    := \1 \left\{ 
           \sqrt{N} \sup_{s \in \RR} |\hat F_N(s) - \Phi(s)| > k_{1-\alpha}
        \right\}
```

let $\beta_N(\Phi)$ be the value of the power function when the null hypothesis is true. By {eq}`eq-kr2`, we have

```{math}
:enumerated: false

\lim_{N \to \infty} \beta_N(\Phi) 
    = \lim_{N \to \infty} 
        \PP \left\{ 
            \sqrt{N}\sup_{s \in \RR} |\hat F_N(s) - \Phi(s)| > k_{1-\alpha}
        \right\}
    = \alpha
```

Hence {eq}`eq-asmposa` is verified, and the test is asymptotically of size $\alpha$. Notice how the construction of the test inverts the confidence set in {ref}`ss-ane`.

The Nikkei data set used here has 364 observations. The value of the test statistic $\sqrt{N}\sup_{s \in \RR} |\hat F_N(s) - \Phi(s)|$ is 1.548. If $\alpha = 0.05$, then, as shown in {ref}`ss-ane`, the critical value $k_{1-\alpha}$ is 1.36. Hence the test statistic exceeds the critical value, and the null hypothesis is rejected.

While the departure from normality is not easy to see in [figure %s](#f-nikkei_ecdf), it can be seen more clearly in a comparison of the histogram and density, as in [figure %s](#f-nikkei_hist). The data are noticeably more peaked and there is evidence of heavy tails.

This application also illustrates how tests often depend on a range of assumptions beyond those that are the focus of a particular hypothesis. Here, for example, a full description of $H_0$ is that standardized returns are IID draws from the standard normal density. The rejection suggests this null is false, but it may be the IID assumption rather than the normality assumption that makes our null a poor fit to the data.[^2]

```{figure} figures/nikkei_hist.svg
:name: f-nikkei_hist

Standardized daily returns for the Nikkei 225
```

#### $p$-Values

Typically, a test that rejects at size 0.05 will also reject at size 0.1, but may not reject at size 0.01. Lower $\alpha$ means less tolerance for type I error, and forces the critical value to become larger. Hence, for a fixed value of the test statistic, the result of the test may switch from reject to don't reject. A natural question then is: What is the smallest value of $\alpha$ at which we can still reject a given test statistic? This value is called the $p$-value.

To state this more precisely, consider again the setting of {ref}`s-hyp`. The null hypothesis is $H_0 : \theta \in \Theta_0$ and, for each $\alpha \in (0,1)$, a test $(T,c(\alpha))$ of size $\alpha$. We assume here that $c(\alpha)$ is determined via the relationship in {eq}`eq-solfc`. In this setting, the $p$-value of the test is defined as

```{math}
:enumerated: false

p(\zdata) := \inf \setntn{\alpha \in (0,1) \, }{\, c(\alpha) <
    T(\zdata) }
```

Roughly speaking, this is the $\alpha$ at which the test switches from reject to don't reject. Typically, $\alpha \mapsto c(\alpha)$ is continuous, and in this case the expression for $p(\zdata)$ reduces to

$$
p(\zdata) := \text{ the } \alpha \text{ such that } c(\alpha) =
    T(\zdata)
$$ (eq-edp)

```{prf:example}
:label: eg-pvsc

Recall the test {eq}`eq-easn3`. Here $c(\alpha) := \Phi^{-1}(1 - \alpha / 2)$, and $c(\alpha)$ is continuous in $\alpha$, so we can apply the definition of $p(\zdata)$ in {eq}`eq-edp`. To solve for $p(\zdata)$, then, we need to solve for $\alpha$ in the expression $\Phi^{-1}(1 - \alpha/2) = |t_N(\zdata)|$. With a bit of rearranging and an application of symmetry, we obtain

$$
p(\zdata) = 2 \Phi(-|t_N(\zdata)|)
$$ (eq-pvsc)
```


### Accepting the Null?

Sometimes failure to reject a null hypothesis comes to be seen as presenting evidence in favor of that hypothesis. This is a slippery slope. Here's one illustration.

In the 1980s econometricians found a new set of toys called unit root tests. First it was discovered that a newly developed test due to {cite}`dickey1979distribution` did not reject the null hypothesis of a unit root when applied to a variety of economic time series ({cite:author}`nelson1982trends` {cite:year}`nelson1982trends`). It was argued that this finding called into question a range of earlier time series methods that assumed trend stationarity, as well as a variety of business cycle studies that assumed stationary innovations around a trend. A raft of new unit root tests appeared. The idea that innovations to major economic time series contain a permanent component became something like a stylized fact.

Let's think about this carefully. For the sake of concreteness, suppose that we are considering unemployment rate data, where the unit root hypothesis was associated with the concept of hysteresis (see, e.g., {cite:author}`mitchell1993testing` {cite:year}`mitchell1993testing`). In this setting the unit root null can be expressed most simply as the hypothesis that $a = 1$ in the AR(1) process

$$
u_{t+1} = a u_t + b + \epsilon_{t+1}
$$ (eq-unemp)

Here $\{u_t\}$ is unemployment, $a$ and $b$ are parameters and $\{\epsilon_t\}$ is a zero-mean innovation (see eq. (5) of {cite:author}`mitchell1993testing` {cite:year}`mitchell1993testing`). Letting $T$ be the length of the sample, the standard one-sided Dickey–Fuller test of the unit root null is to

(i) estimate the parameters as $\hat a$ and $\hat b$ by least squares and

(ii) reject the null if $T (\hat a - 1) < c$, where $c$ is a critical value.

Intuitively, when $a < 1$, the test statistic will evaluate to a negative number in large samples. While we omit the derivation, {cite}`dickey1979distribution` showed that the test statistic converges in distribution under the null and tabulated critical values for different test sizes.

A number of studies have been unable to reject a unit root null for unemployment rate data (see {cite:author}`mitchell1993testing` {cite:year}`mitchell1993testing` for an earlier study or {cite:author}`cheng2012hysteresis` {cite:year}`cheng2012hysteresis` for a study using a newer set of tests). How should we interpret failure to reject the null in this setting? One possibility is that {eq}`eq-unemp` is a good model for the data with $a=1$. This is not the only possibility, however. Another possibility is that the data are being generated by some other process against which this test has little power.

The second possibility seems more plausible than the first. Indeed, if $a=1$ in {eq}`eq-unemp`, then $\PP\{u_t \in B\} \to 0$ as $t \to \infty$ for every bounded $\bB$-measurable $B \subset \RR$. But employment rates don't diverge. A more likely story is that while unemployment rates are highly persistent for some range of values, there are equilibrium forces that become increasingly stronger outside this band. For example, the further the unemployment rate falls, the more tightness in the labor market boosts wages, encouraging employers to substitute toward other factors of production and putting upward pressure on the unemployment rate.

We can easily reconcile these equilibrium forces with the idea that the unemployment rate is highly persistent at commonly observed values by adopting a nonlinear model. For example, consider the following simple variation on a random walk:

$$
u_{t+1} = h(u_t) + \epsilon_{t+1}
$$ (eq-unempnl)

where $h$ is the generalized logistic function depicted in [figure %s](#f-glu), along with the 45 degree line. For unemployment rates in a band between about 5 and 15, the function has a slope close to one, and the process exhibits strong persistence, just like a random walk. For values above and below this number the shape of $h$ causes drift back to the band.

```{figure} figures/glu.svg
:name: f-glu

A nonlinear process for unemployment dynamics
```

```{figure} figures/unempl_sim.svg
:name: f-unempl_sim

Simulated unemployment dynamics
```

{numref}`list-unit_root` gives code for a simulation when this model generates the data. The shocks $\{\epsilon_t\}$ are independent Gaussian. In each loop a time series of length 100 is drawn and the Dickey–Fuller test statistic is calculated. The size of the test is set to 0.05, in which case the critical value for the test is $-13.7$ (see, e.g., {cite:author}`hayashi2000econometrics` {cite:year}`hayashi2000econometrics`, p. 576). If you run the code you will see that for this random seed, the rejection frequency over 5,000 repetitions is around 0.05. This is the same rejection frequency as when the null is true. For this data size, the great majority of tests would fail to reject the unit root model, even though the data-generating process is {eq}`eq-unempnl`.

```{code-block} julia
:caption: DF test applied to a nonlinear time series (Julia)
:name: list-unit_root
using LinearAlgebra, Random, Statistics
N, alpha = 100, 0.05  # Length of time series, size of test
critical = -13.7      # Critical value for DF test of size 0.05
Random.seed!(123)     # Set random seed

function h(x, a=0.0, k=22, b=.16, q=0.6, m=10)
    return a + (k - a) / ((1 + q * exp(-b * (x - m)))^(1/q))
end

function dickey_fuller_stat(y)
    y1 = y[1:end-1]
    y2 = y[2:end]
    z2 = y2 .- mean(y2)
    z1 = y1 .- mean(y1)
    rho = dot(z1, z2) / dot(z1, z1) 
    return length(y1) * (rho - 1)
end

function sim_dgp(N) 
    y = Vector{Float64}(undef, N)
    y[1] = 8  # Initial value of unemployment rate
    for i in 1:N-1
        y[i+1] = h(y[i]) + 0.5 * randn()
    end
    return y
end

M = 5000
rejections = 0
for m in 1:M
    y = sim_dgp(N)
    global rejections += (dickey_fuller_stat(y) < critical)
end
println("Fraction under critical value: ", rejections / M)
```

The interesting thing about this experiment is that the nonlinear data-generating process has very different properties from the null. Not only is it stationary, it is uniformly ergodic, as can be seen by checking the conditions of theorem 16.0.2 of {cite}`meyn2009markov`. This is the strongest possible form of stationarity after the IID property. It means that the effect of initial conditions dies out at a geometric rate that is bounded away from one uniformly over the state space. There is no "permanent component.\"

### Statistical Tests in Economics

While the mechanics of hypothesis testing are relatively straightforward, the scientific paradigm based around hypothesis testing is the subject of continual debate. Macroeconomics is an interesting case in point. In the 1970s the rational expectations revolution ushered in a new class of macroeconomic models. One of the first things the proponents of rational expectations did was to test these models against data. Thomas Sargent's account runs as follows ({cite:author}`evans2005interview` {cite:year}`evans2005interview`):

> My recollection is that Bob Lucas and Ed Prescott were initially very enthusiastic about rational expectations econometrics. After all, it simply involved imposing on ourselves the same high standards we had criticized the Keynesians for failing to live up to. But after about five years of doing likelihood ratio tests on rational expectations models, I recall Bob Lucas and Ed Prescott telling me that those tests were rejecting too many good models.

In the following years, many macroeconomists moved away from formal hypothesis testing. A standard response to classical testing was: "I already know that my model is wrong. Rejection of my model using standard inference tells me nothing new." For example, {cite}`kydland1991econometrics` define econometrics as a process without tests:

> The main steps in econometric analysis are as follows: defining the question; setting up the model; calibrating the model; and reporting the findings.

How does this line of thinking compare with the standard statistical paradigm? The standard paradigm certainly recognizes that all models are wrong. All sensible models are wrong along some dimensions because they are, by definition, simplified representations of reality. Rejection of models is intended to facilitate survival of the fittest: Those models that explain observed phenomena effectively, and, at the same time, are most difficult to reject, will survive.[^3]

At the same time, it's reasonable to feel sympathy for the frustration with classical hypothesis testing expressed by some economists. While the hypothesis that water boils at 100 degrees C will never be refuted (at least at sea level), in economics even the best models capture only some phenomena of interest. Inevitably they fail to match with reality when viewed in other ways. With enough data almost all sensible economic models can be rejected. But should we discard them? Perhaps some of these models are still useful in guiding our thinking.

```{prf:example}

The consumption-based asset pricing model of {cite}`lucas1978asset` fails to provide a good fit for returns on most classes of assets. But the key pricing equation still provides an important benchmark, as well as a way for econometricians to organize their thoughts on asset prices in a quantitative manner.
```


How to resolve these different points of view? Some researchers argue that calibration and matching of moments are sufficient as a means of validating their models. But this form of validation imposes little discipline on model builders. Clever economists can find multiple ways to replicate most interesting phenomena observed in the data, each with their own distinct policy implications. More important, the objective of matching moments is largely an in-sample criterion. The dangers of such an approach were illustrated in [figure %s](#f-rvempr).

One middle ground is testing of partially specified models, where tests are only concerned with a specific restriction implied by theory. A second path forward is more emphasis on out-of-sample performance in economic studies, which provides a much stricter test of true understanding. A third is offered by the Bayesian approach to model selection, where the emphasis is on comparing and improving models rather than rejection. A fourth is improvements in economic theory, to better match our models to the data.

## Further Reading

For a technical discussion of the results surrounding the convergence to the Kolmogorov distribution in {eq}`eq-kdf`, see section 19.4 of {cite}`van2000asymptotic`.

For a discussions of statistical testing in econometrics see, for example, {cite}`mccloskey1996standard` and {cite}`freedman2009statistical`. Discussions of calibration in relation to estimation and hypothesis testing can be found in {cite}`kydland1991econometrics`, {cite}`kim1995econometric`, and many other places. {cite}`hansen2014nobel` presents the main arguments in favor of tests of restrictions in partially specified models. {cite}`silver2012signal` provides an interesting read on estimation and uncertainty with applications in forecasting.

An excellent summary of unit root tests is given in {cite}`hayashi2000econometrics`. For a critical overview of unit root testing, see {cite}`cochrane1991critique`.

## Exercises

```{exercise}
:label: ex-ctdr

Confirm [fact %s](#fa-ctdr).
```


```{exercise}
:label: ex-nbsez

Show that {eq}`eq-nbsez` is valid when {eq}`eq-asnorm2` and {eq}`eq-secfs` hold.
```


```{exercise}
:label: ex-nwald

Let $\alpha \in (0,1)$ and set $z_{\alpha/2} := \Phi^{-1}(1 - \alpha / 2)$. Show that if $t_N(\boldx) \tod \nN(0,1)$ whenever $H_0$ is true, and $T_N(\boldx) = |t_N(\boldx)|$, the sequence of tests $\phi_N(\boldx) := \1\{ T_N(\boldx) > z_{\alpha/2} \}$ is asymptotically of size $\alpha$.
```


````{exercise}
:label: ex-normslzcons

If a test is well-designed, then the rejection probability will converge to one as $N \to \infty$ whenever $H_0$ is false:

```{math}
:enumerated: false

\beta_N(\theta) \to 1 \text{ as } N \to \infty
        \quad \text{whenever} \quad
        \theta \notin \Theta_0
```

Such a test is said to be **consistent**. Show that the test in [example %s](#eg-asymptest) is consistent whenever $\alpha \in (0,1)$.[^4]
````


```{exercise}
:label: ex-chst

The chi-squared goodness of fit test is used to test whether or not a given data set was generated from a particular discrete distribution, described by a PMF $p_1,\ldots,p_J$ over values $1,\ldots,J$. More precisely, let $x_1,\ldots,x_N$ be $N$ random variables, each $x_n$ taking integer values between $1$ and $J$. The null hypothesis of the test is that the sample is IID with $\PP\{x_n = j\} = p_j$ for $n \in \{1,\ldots,N\}$ and $j \in \{1,\ldots,J\}$. The test statistic is given by

$$
X := N \sum_{j=1}^J \frac{(q_j - p_j)^2}{p_j}
$$ (eq-ts)
where $q_j$ is the fraction of the sample $x_1,\ldots,x_N$ taking the value $j$. Write a function called `chsqts` that takes arguments `observations` and `p`, where `observations` is a vector storing the sample $x_1,\ldots,x_N$ and `p` is a vector storing the values $p_1,\ldots,p_J$. The function call `chsqts(observations, p)` should return the value $X$ in equation {eq}`eq-ts`.[^5]
```


```{exercise}
:label: ex-ch10_intervals_tests-auto-1

This exercise continues on from [ex. %s](#ex-chst). Under the null hypothesis, $X$ in {eq}`eq-ts` is asymptotically chi-squared with $J - 1$ degrees of freedom. Let $J = 3$ with $p_1 = 0.2$, $p_2 = 0.2$ and $p_3 = 0.6$. Let $N = 20$. By repeatedly simulating 20 IID observations $x_1,\ldots,x_{20}$ from this PMF,[^6] generate 5,000 independent observations of the statistic $X$, and store them in a vector called `obsX`. Plot the ECDF of `obsX`. In the same figure, plot the chi-squared CDF with 2 degrees of freedom. The functions should be close. (You can improve the fit further by taking $N$ larger.)
```


### Solutions to Selected Exercises

````{solution} ex-ctdr

Assume the conditions of [fact %s](#fa-ctdr). Let $s_N^2$ be the sample variance. Combining [fact %s](#fa-astudt), [example %s](#eg-sdq), [fact %s](#fa-ismsv) and equation {eq}`eq-cisk`, the distribution of

```{math}
:enumerated: false

\sqrt N \frac{(\bar x_N - \mu)}{\sigma}  \,
        \sqrt{\frac{N-1}{N s^2_N /\sigma^2}}
```

has a $t$ distribution with $N-1$ degrees of freedom. Some algebra gives {eq}`eq-cisk2`.
````


````{solution} ex-nbsez

We aim to show that $(\hat \theta_N - \theta)/\se(\hat \theta_N)$ converges in distribution to a standard normal when $\sqrt{N} (\hat \theta_N - \theta) \tod \nN(0,v(\theta))$ and $\sqrt{N} \se(\hat \theta_N) \toprob \sqrt{v(\theta)}$. To see this, observe that

```{math}
:enumerated: false

\frac{\hat \theta_N - \theta}{\se(\hat \theta_N)} 
        = \xi_N \eta_N 
        \quad \text{where} \quad
        \xi_N :=  \frac{\sqrt{N} (\hat \theta_N - \theta)}{\sqrt{v(\theta)}}
        \quad \text{and} \quad
        \eta_N := \frac{\sqrt{v(\theta)}}{\sqrt{N} \se(\hat \theta_N)}
```

Using the various rules for convergence in probability and distribution (check them carefully—see [facts %s](#fa-reconpro) and [%s](#fa-cmtetc)) we obtain $\xi_N \tod Z \sim \nN(0,1)$ and $\eta_N \toprob 1$. Applying Slutsky's theorem now gives $\xi_N \eta_N \tod Z$, as was to be shown.
````


````{solution} ex-nwald

Fix $\alpha \in (0, 1)$, and let $z \sim \nN(0,1)$. In view of {eq}`eq-symcdf`, we have $\PP \{ |z| > z_{\alpha/2} \} = \alpha$. If $H_0$ is true, then $t_N(\boldx) \tod z$ by assumption. Since $g(s) := |s|$ is continuous, [fact %s](#fa-cmtetc) implies that $|t_N(\boldx)| \tod |z|$. As a result we have

```{math}
:enumerated: false

\lim_{N \to \infty} \beta_N(\theta)
        = \lim_{N \to \infty} \PP_{\theta} \{ |t_N(\boldx)| > z_{\alpha/2} \} 
        = \PP \{ |z| > z_{\alpha/2} \} 
        = \alpha
```

This confirms {eq}`eq-asmposa`, and the exercise is done.
````


[^1]: Returns are standardized by subtracting the sample mean and dividing by the sample standard deviation. Data are sourced from the Federal Reserve of St. Louis FRED data set.

[^2]: The test we've implemented can be modified to tackle the case of dependent data. See {cite}`negri2010review` for one such test and more references.

[^3]: Incidentally, given this description of model selection, it might appear there is an easy way to produce a model that survives forever in the pool of currently acceptable theories, without ever being culled via rejection: just make sure that the model has no testable implications. However, this strategy does not work because a model with no testable implications cannot be considered as a theory of anything. The *definition* of a scientific theory, as proposed by Karl Popper, is that the theory has one or more testable implications. In other words, the theory can be *falsified*.

[^4]: Hint: In essence, you need to show that the absolute value of the test statistic $|t_N|$ in {eq}`eq-easn2` diverges to infinity when $\theta \neq \theta_0$. Try the add and subtract strategy, replacing the expression $\bar x_N - \theta_0$ in {eq}`eq-easn2` with $(\bar x_N - \theta) + (\theta - \theta_0)$.

[^5]: `R` has a built-in function called `chisq.test` that you can use to test your results.

[^6]: In particular, draw each $x_n$ such that $\PP\{x_n = j\} = p_j$ for $j=1,2,3$.
