---
title: "Regularization"
label: c-rms
---

In several sections of this book we touched on the topic of regularization (see, e.g., {ref}`sss-lim` and {ref}`ss-ermls`). A variety of statistical procedures and machine learning algorithms employ regularization (under different names) to improve out-of-sample fit. Good out-of-sample fit means generalization from observed data, which, as we've stressed before, is the key problem of statistics. This chapter introduces a number of methods that use regularization and discusses their statistical properties.

(s-panp)=
## Nonparametric Density Estimation
Nonparametric density estimation is an application of regularization to the problem of recovering distributions from data. It combines the data with a prior belief that probability mass most likely falls in places other than just the sample points observed so far. As well as being of interest from a theoretical perspective, nonparametric density estimation is used in a great variety of applied studies. We begin our analysis with a review of parametric density estimation and then proceed to nonparametric methods.

(ss-ped)=
### Introduction
Suppose our data consist of IID observations $\boldx_1, \ldots, \boldx_N$ from unknown distribution $P$ on $\RR^d$. We assume throughout this section that $P$ is absolutely continuous. Our aim is to estimate the density of $P$, denoted below by $f$.

We know how to do this in a parametric setting. For example, let's add the assumption that $f$ belongs to the class of normal densities on $\RR$, so that $f = f(\cdot \, ; \mu, \sigma) =$ the normal density for distribution $\nN(\mu, \sigma^2)$. The MLEs of the parameters are $\hat \mu_N := \bar x_N$ and $\hat \sigma_N := s_N$ . Plugging these back into $f$ gives density estimate $f(\cdot \,; \bar x_N, s_N)$. Since $\bar x_N$ and $s_N$ are consistent (see {ref}`ss-ape`), the random density $f(\cdot \,; \bar x_N, s_N)$ will be close to $f(\cdot \,; \mu, \sigma)$ with high probability for large $N$.

We can clarify convergence of the densities themselves if we extend the notion of consistency from vectors to densities. To extend consistency to densities we need a notion of global deviation between densities. Perhaps the most important measures of global distance are the $L_p$ distances. Let's now define and discuss them.

For any $\bB$-measurable function $f \colon \RR^d \to \RR$ and $p \geq 1$, we set

$$
\| f \|_p 
    := \left\{ \int | f |^p \right\}^{1/p}
    := \left\{ \int | f(\bolds) |^p \diff  \bolds \right\}^{1/p}
$$ (eq-l1norm)

where integration is over all of $\RR^d$. If this expression is finite, then we write $f \in L_p$. For densities $f$ and $g$ on $\RR^d$ the $L_p$ distance between these densities is then defined as

$$
d_p(f, g) := \| f - g \|_p
$$ (eq-l1dist)

The norm $\| \cdot \|_p$ satisfies most of the properties of the norms we've met so far. For example, we have the triangle inequality

$$
\| f - g\|_p \leq \|f - h\|_p + \| h - g \|_p
$$ (eq-lptri)

for all $p \geq 1$ and $f, g, h \in L_p$. See, for example, theorem 5.1.5 of {cite}`dudley2002real`.

Specializing to $p=1$ gives the $L_1$ distance, which sums over absolute deviation. Specializing to $p=2$ gives the popular $L_2$ distance, which is a variation of the $L_2$ distance we used in {ref}`ss-tsl2`. See, in particular, {eq}`eq-ndrv`, where distance between random variables is evaluated based on expected deviation.

The $L_2$ distance is popular largely because it's so tractable, and because it's similarity to Euclidean vector distance means that we have $L_2$ analogies for a variety of fundamental results on vector space. Nevertheless, the $L_1$ distance is arguably a better choice for studying deviation between densities. For example, the $L_1$ distance between densities is always well-defined (see [ex. %s](#ex-l1wd)), which is essential if one hopes to provide universal consistency results.

````{prf:property}

**Scheffés lemma**: If $\{f_n\}$ and $f$ are densities on $\RR^d$, then

```{math}
:enumerated: false

f_n(\bolds) \to f(\bolds)
        \text{ for all } \bolds \text{ in } \RR^d
        \;\; \implies \;\;
        \| f_n - f \|_1 \to 0
```
````


```{prf:property}
:label: fa-pinsk

For any densities $f, g$ and $h$ on $\RR^d$ we have

(i) $\| f - g \|_1 \leq  \sqrt{ 2 D(f, g) }$, where $D(f, g)$ is the KL deviation defined in {eq}`eq-kldist`, and

(ii) $\| f - g \|_1 = 2 \sup_{B \in \bB(\RR^d)} \left| \int_B f - \int_B g \right|$.
```


The bound in (i) is called **Pinsker's inequality**, while (ii) is called **Scheffé's identity**. Scheffé's identity tells us that $L_1$ distance measures something that we directly care about: when $L_1$ deviation is small, so is the maximal deviation between probabilities assigned to events.[^1]

In what follows, we will say that a sequence $\{\hat f_N\}$ of random densities on $\RR^d$ is **$L_p$-consistent** for a density $f$ on $\RR^d$ if

```{math}
:enumerated: false

\| \hat f_N - f \|_p \toprob 0
    \quad \text{as} \quad
    N \to \infty
```

```{prf:example}
:label: eg-ndl1

Let $\hat f_N = f(\cdot \,; \bar x_N, s_N)$ be the $N$th element of the sequence of normal densities described above, where $x_1, \ldots, x_N$ are independent draws from a normal density $f = f(\cdot \,; \mu, \sigma)$ and $\bar x_N$ and $s_N$ the sample mean and standard deviation respectively. This sequence of densities is $L_1$-consistent for $f$. See [exercise %s](#ex-ndl1).
```


#### Failure of Consistency

The risk with the parametric approach is that the parametric assumption is incorrect, in the sense that the parametric class doesn't contain the density generating the data or any good approximation. If this is the case, a parametric approach is typically not consistent. More precisely, if we estimate $f$ with parametric class $\{f_{\boldtheta}\}_{\boldtheta \in \Theta}$, then the $L_p$ deviation between our estimate and $f$ is bounded below by

$$
\delta(f) := \inf_{\boldtheta \in \Theta} \| f - f_{\boldtheta} \|_p
$$ (eq-xi)

This value will be zero only when $f$ can be attained as the limit of elements of $\{f_{\boldtheta}\}_{\boldtheta \in \Theta}$.

```{prf:example}

Consider again the setting of {prf:ref}`eg-ndl1` but now suppose that the true density $f$ is not Gaussian. Then either the sequence $\hat f_N$ is not $L_1$-consistent for any density, or, if it is $L_1$-consistent for some density, then that density is not $f$. The reason is that $\delta(f)$ in {eq}`eq-xi` is always positive when the parametric class is Gaussian and $f$ is not, since the set of normal densities is closed under the taking of limits in $L_1$.
```


(ss-npkde)=
### Kernel Density Estimation
Sometimes we can make good choices for parametric classes by using descriptive statistics or by appealing to some theory with sharp quantitative implications. At other times this is difficult. In such settings it is best to use a nonparametric approach that is consistent under weaker assumptions.

Let's start with an intuitive discussion and then turn to theory. Suppose that we have IID data $\boldx_1,\ldots,\boldx_N$ generated from unknown density $f$ on $\RR^d$. To estimate $f$ using the data we will employ a **kernel density estimator** (KDE), which takes the form

$$
\hat f_N(\bolds) 
    := \frac{1}{N h^d} \sum_{n=1}^N 
        K \left( \frac{\bolds - \boldx_n}{h} \right)
$$ (eq-npkde)

Here $K$ is called the **kernel function** of the estimator, and $h$ is called the **bandwidth**. The kernel function $K$ is required to be a density on $\RR^d$ but is otherwise unrestricted. The bandwidth $h$ is any positive number. {prf:ref}`ex-vnpde` asks you to confirm that $\hat f_N$ is always a density.

To get a feeling for the estimate $\hat f_N$ in {eq}`eq-npkde`, let's look at a simple instance created from just three data points $x_1, x_2, x_3$ on $\RR$. For $K$ we take the standard normal density. Since $N=3$, the function $\hat f_N$ is the sum of three individual functions, the $n$th of which can be written as

```{math}
:enumerated: false

g_n(s) 
        = \frac{1}{N h} 
        K \left( \frac{s - x_n}{h} \right)
        = 
        \frac{1}{N h}
        \frac{1}{\sqrt{2 \pi} }
           \exp 
           \left\{ 
               - \frac{(s - x_n)^2}{2 h^2} 
           \right\}
```

This is $\frac{1}{N}$ times the density of a $\nN(x_n, h^2)$ random variable. So the three functions $g_n$ are smooth bumps, each centered on one of the data points and with the degree of concentration around $x_n$ governed by $h$. The top panel of [figure %s](#f-npkde) shows these three functions along with the data points for when $h=1$. The black line is the pointwise sum of these functions, which is the density estimate $\hat f_N$. By construction, it is large near the data points and small away from data points.

```{figure} figures/npkde.svg
:name: f-npkde

Nonparametric KDE, different bandwidths
```

The lower panel of [figure %s](#f-npkde) shows the same functions when the bandwidth is increased to $1.4$. Each $g_n$ becomes more spread out, and the sum $\hat f_N$ is smoother.

As discussed in more detail below, the role of the bandwidth is to add smoothing to the empirical distribution, so we can generalize from a finite sample. But what is the right amount of smoothing? {numref}`f-one_dim_kde` illustrates the trade-off associated with smoothing. The shaded distribution is the density $f$, which we imagine to be unknown. 40 observations are drawn from this distribution and represented as dots. The black line in each panel is the KDE $\hat f_N$ built from these observations. The kernel is the standard normal density again, while the bandwidth varies across the panels.

As the bandwidth goes to zero, the kernel density becomes similar to the empirical distribution, with all probability mass in very small regions around the sample points. Hence, when we use a very small bandwidth , we make the mistake of treating the empirical distribution as the true distribution. This is overfitting. At the same time, excessive smoothing adds too much bias, hiding the features of the true distribution.[^2]

```{figure} figures/one_dim_kde.svg
:name: f-one_dim_kde

Effect of changing the bandwidth
```

The optimal bandwidth in terms of minimizing $L_p$ deviation depends on the unknown density $f$. For example, if $f$ is smooth, then a relatively large bandwidth should be used. Of course $f$ is unknown. As such there are two standard approaches. One is to make assumptions on $f$ and choose the bandwidth accordingly. Another is to use cross-validation. For a review of both procedures see {cite}`scott2015multivariate`. We'll say more about cross-validation in the context of ridge regression below.

### Theory

A natural way to study kernel density estimates is via the notion of convolutions. The convolution of an arbitrary distribution $Q$ and a density $K$ on $\RR^d$ is the density on $\RR^d$ defined by

$$
(K \star Q) (\bolds') = \int K(\bolds' - \bolds) Q(\diff \bolds)
    \qquad (\bolds' \in \RR^d)
$$ (eq-gconv)

We already encountered this concept in the discrete case in [exercise %s](#ex-discreteconv), where you were asked to prove a version of

```{prf:property}
:label: fa-ioc

For any density $K$ and arbitrary distribution $Q$ on $\RR^d$, the density $K \star Q$ equals $\lL(\boldx + \boldy)$ when $\boldx$ and $\boldy$ are independent with $\lL(\boldx) = K$ and $\lL(\boldy) = Q$.
```


```{prf:example}

If $K = \nN(0, \sigma^2)$ for some $\sigma > 0$ and $Q$ is a distribution on $\RR$ that puts mass $q_n$ on points $s_1, \ldots, s_N$, then by {eq}`eq-gconv` and the rule for integrating over discrete distributions (see {eq}`eq-hasp2`),

$$
(K \star Q)(s') = \sum_{n=1}^N K(s' - s_n) q_n
$$ (eq-gconvd)
This distribution is a mixture of normals. {prf:ref}`ex-sioc` illustrates the connection between {eq}`eq-gconvd` and [fact %s](#fa-ioc).
```


We will be particularly interested in a certain class of convolutions, induced by densities of the form

$$
K_h(\bolds) := \frac{1}{h^d} K \left( \frac{\bolds}{h} \right)
$$ (eq-kh)

where $K$ is any density and $h > 0$ is a parameter. The density $K_h$ in {eq}`eq-kh` is the density of $h \boldx$ when $\boldx$ is a random vector on $\RR^d$ with density $K$. (See [ex. %s](#ex-skh).)

```{prf:example}
:label: eg-ctn

{numref}`f-convolve` shows the convolution of a tent-shaped distribution $Q$ and $K_h$ when $K$ is the standard normal density and $h$ takes different values.
```


```{figure} figures/convolve.svg
:name: f-convolve

Smoothing induced by convolutions
```

Note how in [figure %s](#f-convolve) the distribution of $K_h \star Q$ is close to $Q$ when $h$ is near zero. This makes sense, since $K_h \star Q$ is the distribution of $x + hy$ where $x$ is drawn from $Q$ and $y$ is standard normal. As we take $h \downarrow 0$, the distribution of $x + hy$ converges to the distribution of $x$. In fact this is always true:

````{prf:property}
:label: fa-stz

Let $f$ and $K$ be any densities on $\RR^d$, and let $K_h$ be defined from $K$ via {eq}`eq-kh`. If $f \in L_p$, then $K_h \star f \in L_p$, and

```{math}
:enumerated: false

\lim_{h \downarrow 0} \, \| K_h \star f - f \|_p = 0
```
````


For a proof, see theorem 9.1 of {cite}`devroye2012combinatorial`.

#### Convolution and KDEs

Let $h$ be a positive number and let $K_h$ be as defined in {eq}`eq-kh`. Using this notation, we can rewrite the KDE in {eq}`eq-npkde` as $\hat f_N(\bolds) = \frac{1}{N} \sum_{n=1}^N K_h(\bolds - \boldx_n)$. Letting $\hat P_N$ denote the empirical distribution of the sample and recalling the expression for integrals with respect to $\hat P_N$, we can also write this as $\hat f_N(\bolds') =  \int K_h(\bolds' - \bolds) \hat P_N(\diff \bolds)$, or, more simply,

$$
\hat f_N = K_h \star \hat P_N
$$ (eq-kdeac)

This equation says something interesting. In {ref}`sss-lim` we first tried to estimate the density by a naive application of the sample analogue principle. The method failed. In {eq}`eq-kdeac` we have an alternative method that employs smoothing, or regularization. We apply the "smoothing operator\" $K_h \star$ in addition to replacing the target distribution with the empirical distribution. As the next section shows, this method can be successful under very mild assumptions.

#### Consistency

It can be proved that, for *any* density $f$ and kernel $K$, the nonparametric kernel density estimator $\hat f_N$ is $L_1$-consistent for $f$. See theorem 9.2 of {cite}`devroye2012combinatorial`. The $L_1$ proof involves a large amount of measure-theoretic machinery. Here we'll state and prove the corresponding $L_2$ result instead.

```{prf:theorem}
:label: t-kdec

Let $f$ and $K$ be densities on $\RR^d$ and elements of $L_2$. If

(i) $\{\boldx_n\}_{n \geq 1}$ is an IID sequence of draws from $f$, and

(ii) the bandwidth sequence $\{h_N\}$ satisfies $h_N \to 0$ and $N h_N^d \to \infty$ as $N \to \infty$,

then the sequence of density estimates $\{\hat f_N\}$ defined in {eq}`eq-npkde` is $L_2$-consistent for $f$.
```


Let's break this claim down into smaller pieces. For the remainder of this section, $\| \cdot \|$ is the $L_2$ norm and $h := h_N$. Using the expression for $\hat f_N$ in {eq}`eq-kdeac` and the triangle inequality in {eq}`eq-lptri`, we have

$$
\| \hat f_N - f \|
    \leq \| K_h \star \hat P_N - K_h \star f \| + \| K_h \star f - f \|
$$ (eq-dnpk)

The first term is called the **estimation error**, and the second is called the **approximation error** or **bias**. The first error is caused by the fact that we only observe the empirical distribution $\hat P_N$ rather than the true distribution $f$. The second error is caused by the smoothing we deliberately added to control the estimation error. (Compare with the error decomposition.)

The fact that $h_N \to 0$ in condition (ii) of [theorem %s](#t-kdec) means that we shrink the amount of smoothing as the sample size increases, thereby reducing the approximation error. The requirement $N h_N^d \to \infty$ ensures that smoothing is not reduced too quickly, allowing us to control estimation error.

One way to think about the relationship between inference and the sample is that, when viewed as a "generalized density\" with all its mass on data points, the empirical distribution is too rough. To generalize we need to add smoothing. At the same time, the need for smoothing decreases as the sample size gets large.

````{prf:proof}
:nonumber:

**Proof of [theorem %s](#t-kdec).**

By [fact %s](#fa-stz), the approximation error converges to zero under the conditions of [theorem %s](#t-kdec). Hence it suffices to show that the estimation error converges in probability to zero. Fix $\delta > 0$. By Chebyshev's inequality, we have

```{math}
:enumerated: false

\PP \left\{ 
            \| K_h \star \hat P_N - K_h \star f \| \geq \delta 
            \right\} 
        = \PP \left\{ 
            \| K_h \star \hat P_N - K_h \star f \|^2 \geq \delta^2 
              \right\} 
        \leq \frac{\xi_N}{\delta^2}
```

where

```{math}
:enumerated: false

\xi_N := 
            \EE \left\{
                \| K_h \star \hat P_N - K_h \star f \|^2
            \right\}
```

To complete the proof, we need only show that $\xi_N$ converges to zero. Let

```{math}
:enumerated: false

\bar K_N (\bolds)
         := (K_h \star \hat P_N)(\bolds) - (K_h \star f)(\bolds) 
         = \frac{1}{N} \sum_{n=1}^N
             \left\{
             K_h(\bolds - \boldx_n) - \EE [ K_h(\bolds - \boldx_n) ]
             \right\}
```

We can then write

$$
\xi_N
        =
        \EE \left\{ 
            \int \left[ \bar K_N(\bolds) \right]^2 
            \diff \bolds
            \right\}
        =
        \int 
        \EE \left\{ \left[ \bar K_N(\bolds) \right]^2 \right\}
            \diff \bolds
$$ (eq-hktz)
(The interchange of order of expectation and integration is valid for nonnegative integrands—see theorem 4.4.5 of {cite}`dudley2002real`.) Since $\bar K_N (\bolds)$ is the sample mean of $N$ IID zero-mean random variables, we have

```{math}
:enumerated: false

\EE \left\{ \left[ \bar K_N(\bolds) \right]^2 \right\}
        = \var  \left[ \bar K_N(\bolds) \right]
        = \frac{1}{N} \var[ K_h(\bolds - \boldx_n) ]
```

Moreover

```{math}
:enumerated: false

\var[ K_h(\bolds - \boldx_n) ]
        = \EE \left\{ [ K_h(\bolds - \boldx_n) ]^2 \right\}
        - \left\{ \EE[ K_h(\bolds - \boldx_n)]  \right\}^2
        \leq 
        \EE \left\{ [ K_h(\bolds - \boldx_n) ]^2 \right\}
```

In summary,

$$
\xi_N \leq 
        \frac{1}{N}
        \int 
        \EE \left\{ [ K_h(\bolds - \boldx_n) ]^2 \right\}
            \diff \bolds
$$ (eq-xintz)
Now observe that, switching the order of integration once more,

```{math}
:enumerated: false

\int 
        \EE \left\{ [ K_h(\bolds - \boldx_n) ]^2 \right\}
            \diff \bolds
        =
        \int 
        \left\{
            \int [ K_h(\bolds - \bolds') ]^2 \diff \bolds
        \right\}
            f(\bolds') \diff \bolds'
```

From the definition of $K_h$ and a change of variable argument,

```{math}
:enumerated: false

\int [ K_h(\bolds - \bolds') ]^2 \diff \bolds
         = \frac{1}{h^{2d}}
            \int \left[ K \left( 
            \frac{\bolds - \bolds'}{h} 
                     \right) 
                 \right]^2 \diff \bolds
         = \frac{1}{h^d}
            \int [ K (\boldu) ]^2 \diff \boldu
```

Putting this together with {eq}`eq-xintz` gives the bound

```{math}
:enumerated: false

\xi_N 
        \leq 
        \int 
        \frac{1}{Nh^d} \, \| K \|^2  \, f(\bolds') \diff \bolds'
        = \frac{1}{Nh^d}  \,\| K \|^2
```

The term $\| K \|^2 := \| K \|^2_2$ is finite by assumption. Recalling that $Nh^d = Nh^d_N \to \infty$, we see that $\xi_N \to 0$ as required. ◻
````


### Commentary

In some fields of science, researchers have considerable knowledge about parametric classes and specific functional forms. For example, the theory of Brownian motion describes how the location of a tiny particle in liquid is approximately normally distributed. Thus the underlying theory provides an exact parametric class. In this kind of setting the parametric approach excels. It allows us to generalize by combining data with information we have on functional forms.

Unfortunately, the quantitative foundations of economics and other social sciences are looser and more prone to shifting over time. Econometricians usually come to the table with less knowledge of parametric classes and functional forms. These facts make nonparametric techniques attractive.

At the same time, nonparametric methods do not solve all our problems. The theoretical results presented above are purely asymptotic. Finite sample results are available, but it isn't possible to obtain strong results in this direction without correspondingly strict assumptions on the target density. There is no uniform rate of convergence for all target densities ({cite:author}`devroye2012combinatorial` {cite:year}`devroye2012combinatorial`, p. 85). This is because nonparametric methods have relatively little structure in the form of prior knowledge and hence require abundant data.

(s-modsel)=
## Controlling Complexity
In the last section we discussed flexible estimation methods that perform well asymptotically. In this section we turn our attention to finite sample properties. A common thread will be ridge regression, which is a popular method of estimation in both econometrics and machine learning (see, e.g., {cite:author}`peysakhovich2015machine` {cite:year}`peysakhovich2015machine`, {cite:author}`kim2014forecasting` {cite:year}`kim2014forecasting`, or {cite:author}`varian2014big` {cite:year}`varian2014big`). Ridge regression also connects to some deep ideas at the heart of finite sample theory, including complexity, prior knowledge and the bias–variance trade-off.

(ss-ridge)=
### Ridge Regression
Let's start off in the classical OLS setting of {ref}`ss-tolsa`, with [assumptions %s](#a-lols)–[%s](#a-bols2) all in force. The usual OLS estimator is $\hboldbeta$ as given in {eq}`eq-thelsqe2`. It can also be expressed as

```{math}
:enumerated: false

\hboldbeta = \argmin_{\boldb \in \RR^K} \sum_{n=1}^N (y_n - \boldx_n^\T \boldb)^2
```

As shown in {ref}`ss-tlse`, the OLS estimator minimizes the empirical risk under quadratic loss when the hypothesis space is the set of linear functions. Under our current assumptions, it is also unbiased for $\boldbeta$ ([theorem %s](#t-olsiu)) and, by the Gauss–Markov theorem, it has the lowest variance among all linear unbiased estimators of $\boldbeta$ ([theorem %s](#t-gm)).

While the Gauss–Markov theorem has historical importance, a more natural way to evaluate estimators is to consider their mean squared error, which tells us directly how much probability mass the estimator puts around the object it is trying to estimate (see {ref}`ss-vvb`). Recalling {eq}`eq-mse`, the MSE of an estimator $\hat \boldb$ of $\boldbeta$ is defined as

```{math}
:enumerated: false

\mse(\hat \boldb, \boldbeta) := \EE \left\{ \| \hat \boldb - \boldbeta \|^2  \right\}
```

As [exercise %s](#ex-dmse2) asks you to show, the following representation is also valid:

$$
\mse(\hat \boldb, \boldbeta) 
    = \EE \left\{ \|  \, \hat \boldb - \EE[\hat \boldb] \, \|^2 \right\}
        + \| \, \EE[\hat \boldb] - \boldbeta \, \|^2
$$ (eq-dmse2)

This equation extends {eq}`eq-dmse`, and tells us that MSE is the sum of a variance and a bias term. Minimization of MSE involves a trade-off between these two terms. Typically, the optimal choice is not at either extreme: MSE is minimized when some amount of bias is admitted.

Applying this idea to the OLS setting, {cite}`hoerl1970ridge` showed that there exists a biased linear estimator that has lower mean squared error than $\hboldbeta$. The estimator is defined as the solution to the modified least squares problem

$$
\min_{\boldb \in \RR_K}
    \left\{
     \sum_{n=1}^N (y_n - \boldx_n^\T \boldb)^2 + \lambda \| \boldb \|^2
     \right\}
$$ (eq-rrp)

where $\lambda \geq 0$ is called the **regularization parameter**. In solving {eq}`eq-rrp`, we are minimizing the empirical risk plus a term that penalizes large values of $\| \boldb \|$. The effect is to "shrink" the solution relative to the unpenalized solution $\hboldbeta$. Some calculus shows that the solution to {eq}`eq-rrp` is

$$
\hboldbeta_{\lambda} 
        := (\boldX^\T \boldX + \lambda \boldI)^{-1} \boldX^\T \boldy
$$ (eq-rrs)

The estimator $\hboldbeta_{\lambda}$ is called the **ridge regression estimator**. Note that

(i) $\hboldbeta_\lambda$ is the OLS estimator when $\lambda = 0$ and

(ii) $\hboldbeta_\lambda$ is biased whenever $\lambda > 0$ ([ex. %s](#ex-rrbiased)).

The following result is proved in {cite}`hoerl1970ridge`.

````{prf:theorem}
:label: t-rrg

Under the OLS [assumptions %s](#a-lols)–[%s](#a-bols2), there exists a $\lambda > 0$ such that

```{math}
:enumerated: false

\mse \left(\hboldbeta_{\lambda}, \boldbeta \right) < 
        \mse \left( \hboldbeta, \boldbeta \right)
```
````


The reduction in MSE over the least squares estimator occurs because, for some intermediate value of $\lambda$, the variance of $\hboldbeta_{\lambda}$ falls by more than enough to offset the bias induced by regularization.

Note that, for suitable choice of $\lambda$, the ridge regression estimator $\hboldbeta_{\lambda}$ outperforms $\hboldbeta$ even though all of the classical OLS assumptions are valid. At the same time, the right choice of $\lambda$ is a nontrivial problem. This problem falls under the heading of model selection, which is the topic treated in the next few sections.

#### Interpretation

The traditional view of ridge regression runs as follows: The standard OLS assumptions are treated as valid. There may, however, be instances where $\boldX^\T \boldX$ is almost singular due to strong correlation between regressors. In this case the process of inverting $\boldX^\T \boldX$ is numerically unstable. We can stabilize the inversion by adding some positive value of $\lambda$ in {eq}`eq-rrs`.

Here's another view: The standard OLS assumptions are implausible. Since our loss function is quadratic, we would ideally like to obtain the regression function $f^*$ but recovering this infinite dimensional object with a finite amount of data is ill-posed. We need to control the complexity of the candidate functions we use to approximate the regression function. The regularization term in ridge regression provides a means of managing complexity.

#### Tikhonov Regularization

Let's build some more intuition for the result in [theorem %s](#t-rrg) by running simulations. We'll frame the simulations in a general setting, of which ridge regression is a special case.

We know that the least squares estimator is the solution to an overdetermined system of equations. There is an existing theory for solving ill-posed linear systems in high dimensions. The basic idea is that any attempt to back out or infer a complex object by solving a system about which we have limited information requires a degree of regularization.

To illustrate, suppose that

(i) $\boldA \boldb = \boldc$ is an overdetermined system, where $\boldA$ is $N \times K$ with $N > K$.

(ii) Due to measurement error, we only observe an approximation $\boldc_0$ of $\boldc$.

(iii) $\boldb^*$ is the (unobservable) least squares solution $\argmin_{\boldb} \| \boldA \boldb - \boldc \|^2$.

In the absence of additional information, a natural approach to approximating $\boldb^*$ is to solve $\boldA \boldb = \boldc_0$ by least squares. An alternative, less obvious approach, is to minimize

$$
m(\lambda) := \| \boldA \boldb - \boldc_0 \|^2 + \lambda \| \boldb \|^2
$$ (eq-mla)

for some small but positive $\lambda$. This second approach is called **Tikhonov regularization**. We are minimizing least squares plus a penalty term.

Let's look at a simulation where $\boldA$ is chosen stochastically but with a tendency towards multicollinearity.[^3] To aid clarity, we first set $\boldb^* := (10,10,\ldots,10)^\T$, and then set $\boldc := \boldA \boldb^*$. By construction, $\boldb^*$ is a solution to the system $\boldA \boldb^* = \boldc$, and also the least squares solution.

```{figure} figures/tikreg.svg
:name: f-tikreg

Effect of Tikhonov regularization, $\lambda=1$
```

Measurement of $\boldc$ is corrupted with a Gaussian shock. In particular, $\boldc_0$ is drawn from $\nN(\boldc, \sigma^2 \boldI)$ where $\sigma$ is a small positive number. We then plot the OLS solution based on $\boldc_0$, which minimizes $\| \boldA \boldb - \boldc_0 \|^2$, and the regularized solution, which minimizes $m(\lambda)$ in {eq}`eq-mla`. The former are plotted against their index in green, while the latter are plotted in blue. The true solution $\boldb^*$ is plotted in black. The result is [figure %s](#f-tikreg).

The figure shows 10 solutions each for the ordinary and regularized solutions, corresponding to 10 draws of $\boldc_0$. The regularized solutions are much closer to the true solution on average.

Note that this result is dependent on a reasonable choice for $\lambda$. If you experiment with the code on the text website, you will see that for very small values of $\lambda$, the regularized solutions are almost the same as the unregularized solutions. Conversely, very large values of $\lambda$ pull the regularized solutions too close to the zero vector.

(ss-ssrr)=
### Subset Selection and Ridge Regression
One problem frequently faced in regression problems is which variables to include. For example, if we are comparing crime rates across different cites, we can think of any number of variables that might be relevant (median wage, unemployment, police density, etc.). The same is true if we are trying to model credit default rates for some group of individuals, or educational attainment across schools. A similar problem arises in time series models, where we want to know how many lags of the state variables to include. The general problem is known as **subset selection**, since we are trying to choose the right subset of all candidate regressors.

Further dimensions appear when we consider basis functions. Given a set of covariates $\boldx$, we have the option to map this into a larger vector $\boldphi(\boldx)$ using basis functions, as discussed in {ref}`ss-basis`. For example, given a single covariate $x$, we may consider mapping it into $\boldphi(x) = (1,x,x^2,\ldots,x^d)$ and regressing $y$ on $\boldphi(x)$. As we saw in [figures %s](#f-ofit1)–[%s](#f-ofit14), a good choice of $d$ is crucial. Choosing $d$ is another example of the subset selection problem because we are trying to decide whether to include the regressor $x^j$ for some given $j$.

Subset selection is a version of the empirical risk minimization problem. Suppose that we have output $y$ and inputs $\boldx \in \RR^K$, in the sense that $\boldx$ contains $K$ candidate regressors. If we want to include all regressors, then we can minimize empirical risk over $\llL$ in {eq}`eq-allaff2`, the hypothesis space of linear functions from $\RR^K$ to $\RR$. If we wish to exclude some subset of regressors, we can set $I \subset \{1,\ldots,K\}$ to be the set of indices of the regressors we want to exclude and regress $y$ on the remainder, which is equivalent to minimizing the empirical risk over the hypothesis space

```{math}
:enumerated: false

\hH_{-I} := \left\{
        \text{ all functions } 
            f(\boldx) = \boldb^\T \boldx
            \text{ with }
            b_k = 0 \text{ for all } k \in I
           \right\}
```

We are back to the problem of choosing a suitable hypothesis space over which to minimize empirical risk.

The subset selection problem has been tackled by many researchers. Well-known approaches include the Akaike Information Criterion (AIC), the Bayesian Information Criterion (BIC), and Mallow's $C_p$ statistic. For example, Mallow's $C_p$ statistic consists of two terms, one increasing in the size of the empirical risk, and the other increasing in $|I|$, the size of the subset selected. The objective is to minimize the statistic, which involves trading off poor fit (large empirical risk) against excess complexity of the hypothesis space (large $|I|$).

One of the problems with subset selection is that $K$ regressors means $2^K$ subsets to step through.[^4] To avoid this computational problem, one alternative is to use ridge regression. With ridge regression, the regularization term leads us to choose an estimate with smaller norm. What this means, in practice, is that the coefficients of less helpful regressors are driven towards zero, thereby "almost excluding" those regressors. While the model selection problem is not solved, it has been reduced to tuning a single parameter.

We can illustrate the idea by reconsidering the regression problem discussed in {ref}`ss-ermls`. {numref}`f-ofit1`–{numref}`f-ofit14` showed the fit we obtained by minimizing empirical risk over larger and larger hypothesis spaces. The hypothesis spaces were the sets $\pP_d$ of degree $d$ polynomials for different values of $d$. For each $d$ we minimized the empirical risk over $\pP_d$, which translates into solving

```{math}
:enumerated: false

\min_{\boldb} \sum_{n=1}^N 
        \left[ y_n - \boldb^\T\boldphi(x_n) \right]^2
        \quad \text{where} \quad
        \boldphi(x) = (x^0, x^1,\ldots, x^d)
```

As discussed above, choosing the right $d$ is isomorphic to subset selection, since we are deciding what powers of $x$ to include as regressors. {numref}`f-rvempr` showed that intermediate values of $d$ did best at minimizing risk.

We can do a similar thing using ridge regression. First, let's take $\pP_{14}$ as our hypothesis space. This space is certainly large enough to provide a good fit to the data, but with empirical risk minimization the result is overfitting (see [figure %s](#f-ofit14)). Here, instead of using empirical risk minimization, we solve the regularized problem

```{math}
:enumerated: false

\min_{\boldb} \sum_{n=1}^N 
        \left\{
            \left[ y_n - \boldb^\T\boldphi(x_n) \right]^2 + \lambda \| \boldb \|^2
        \right\}
```

for different values of $\lambda$. The data used here are exactly the same data used in the original [figures %s](#f-ofit1)–[%s](#f-ofit14) from {ref}`ss-ermls`. The solution for each $\lambda$ we denote by $\hboldbeta_{\lambda}$, which is the ridge regression estimator, and the resulting prediction function we denote by $\hat f_{\lambda}$, so that $\hat f_{\lambda}(x) = \hboldbeta_{\lambda}^\T \boldphi(x)$.

The function $\hat f_{\lambda}$ is plotted in red for increasingly larger values of $\lambda$ over [figures %s](#f-ridgeplot1)–[%s](#f-ridgeplot8). The black line is the risk-minimizing function. In [figure %s](#f-ridgeplot1), the value of $\lambda$ is too small to impose any real restriction, and the procedure overfits. In [figure %s](#f-ridgeplot5), the value of $\lambda$ is larger and the fit is good. In [figure %s](#f-ridgeplot8) the value of $\lambda$ is too large and the estimate is poor.

As in {ref}`ss-ermls`, we can compute the risk of each function $\hat f_{\lambda}$, since we know the underlying model (see {eq}`eq-csmod`). The risk is plotted against $\lambda$ in [figure %s](#f-ridge_risk). The $x$-axis is on log-scale. The risk is smallest for small but nonzero values of $\lambda$.

```{figure} figures/ridgeplot1.svg
:name: f-ridgeplot1

Fitted polynomial, $\lambda \approx 0$
```

```{figure} figures/ridgeplot5.svg
:name: f-ridgeplot5

Fitted polynomial, $\lambda \approx 0.0004$
```

```{figure} figures/ridgeplot8.svg
:name: f-ridgeplot8

Fitted polynomial, $\lambda \approx 18$
```

```{figure} figures/ridge_risk.svg
:name: f-ridge_risk

Risk of $\hat f_{\lambda}$ plotted against $\lambda$
```

(ss-abp)=
### Bayesian Methods and Regularization
The ideal case with model selection is that we have clear guidance from economic theory on which regressors to include, which functional forms to use, which values of our regularization parameter to choose, and so on. If theory or prior knowledge provides this information, then every effort should be made to exploit it. One technique for injecting prior information into statistical estimation is via Bayesian analysis. Let's now look at Bayesian linear regression and how it compares to ridge regression.

Suppose that our regression data take the linear form $\boldy = \boldX \boldbeta + \boldu$. To simplify the presentation we will assume that $\boldX$ is nonrandom. (Taking $\boldX$ to be random leads to the same conclusions but with a longer derivation). As before, $\boldu$ is random and unobservable. The new feature provided by the Bayesian perspective is that we take $\boldbeta$ to be random and unobservable as well. In addition we are assumed to have subjective prior beliefs regarding likely values for these variates, expressed in the form of probability distributions. Here we take the priors to be $\lL(\boldu) = \nN(\boldzero, \sigma^2 \boldI)$ and $\lL(\boldbeta) =\nN(\boldzero, \tau^2 \boldI)$.

Given our model $\boldy = \boldX \boldbeta + \boldu$, our prior on $\boldu$ implies that the density of $\boldy$ given $\boldbeta$ is $\nN(\boldX \boldbeta, \sigma^2 \boldI)$. In generic notation, we can write our distributions as

$$
p(\boldy \given \boldbeta) = \nN(\boldX \boldbeta, \sigma^2 \boldI)
    \quad \text{and} \quad
    p(\boldbeta) = \nN(\boldzero, \tau^2 \boldI)
$$ (eq-priors)

Applying Bayes' law (see {eq}`eq-brdc`) to the pair $(\boldy, \boldbeta)$, we obtain

$$
p(\boldbeta \given \boldy) 
        = \frac{p(\boldy \given \boldbeta) p(\boldbeta)}{p(\boldy)}
$$ (eq-bf2)

The left-hand side is the posterior density of $\boldbeta$ given the data $\boldy$, and represents our new beliefs updated from the prior on the basis of the data $\boldy$.

One way to summarize the information contained in the posterior is by examining its maximum value. The maximizer of the posterior is called the **maximum a posteriori** (MAP) probability estimate. Taking logs of {eq}`eq-bf2` and dropping the term that does not contain $\boldbeta$, it can be expressed as

$$
\hboldbeta_{M} := \argmax_{\boldbeta} 
    \left\{ \ln p(\boldy \given \boldbeta) + \ln p(\boldbeta) \right\}
$$ (eq-betamap0)

Inserting our functional forms into {eq}`eq-priors`, dropping constant terms and multiplying by $-1$, we obtain the expression

$$
\hboldbeta_{M} = \argmin_{\boldbeta} 
    \left\{ 
     \sum_{n=1}^N (y_n - \boldx_n^\T \boldbeta)^2 
    + \frac{\sigma^2}{\tau^2} \| \boldbeta \|^2 
    \right\}
$$ (eq-betamap)

This is precisely the penalized least squares problem {eq}`eq-rrp`, where the regularization parameter $\lambda$ is equal to $(\sigma/\tau)^2$. In view of {eq}`eq-rrs`, the solution is

```{math}
:enumerated: false

\hboldbeta_{M} 
        := (\boldX^\T \boldX + (\sigma/\tau)^2 \boldI)^{-1} \boldX^\T \boldy
```

Thus Bayesian estimation provides a principled derivation of the penalized objective function from ridge regression. Bayesian analysis provides the same effect as Tikhonov regularization, but now regularization arises out of combining prior knowledge with the data. Moreover, at least in principle, the value $(\sigma/\tau)^2$ is part of our prior knowledge, and hence there is no model selection problem.

In practice, one could of course question the assertion that we have so much prior knowledge that the regularization parameter $\lambda := (\sigma/\tau)^2$ is pinned down. If such knowledge is lacking, then we are back at the model selection problem. In the next section we forgo the assumption that this strong prior knowledge is available, and consider a more automated approach to choosing $\lambda$.

(ss-cval)=
### Cross-Validation
A natural way to think about model selection is to think about minimizing risk. Recall that, given loss function $L$ and a system producing input–output pairs $(\boldx, y) \in \RR^{K+1}$ with joint distribution $P$, the prediction risk of a function $f \colon \RR^K \to \RR$ is the expected loss

```{math}
:enumerated: false

R(f) := \EE [ L(y, f(\boldx)) ] = \int \int L(t, f(\bolds)) 
        P(\diff t, \diff \bolds)
```

that occurs when we use $f(\boldx)$ to predict $y$. Now suppose that we observe $N$ IID input–output pairs $\zdata := \{(\boldx_1, y_1), \ldots (\boldx_N, y_N) \}$. Given a selection of models, we would like to find the one that takes this data set and returns a predictor $\hat f$ such that $\hat f$ has lower prediction risk than the predictors returned by the other models.

We have to be precise in defining risk because, if we define it as $\EE [ L(y, \hat f(\boldx)) ]$, then we are taking expectation over all randomness, including that in $\hat f$, which depends on the data set $\zdata$. What we want to do now is take the data set as given, and see how well we can do in terms of predicting new values as evaluated by expected loss. Hence we define the prediction risk of $\hat f$ to be the expected loss taking $\zdata$ (and hence $\hat f$) as given:

```{math}
:enumerated: false

R(\hat f \given \dD) := \EE [ L(y, \hat f(\boldx)) \given \dD] 
    = \int \int L(t, \hat f(\bolds)) P(\diff t, \diff \bolds)
```

If we have a collection of models $M$ indexed by $m$, and $\hat f_m$ is the predictor produced by fitting model $m$ with data $\dD$, then we would like to find the model $m^*$ such that

```{math}
:enumerated: false

R(\hat f_{m^*} \given \dD) \leq R(\hat f_m \given \dD)
    \quad \text{for all } 
    m \in M
```

The obvious problem with this idea is that risk is unobservable. If we knew the joint distribution $P$ then we could calculate it; but then again, if we knew $P$ there would be no need to estimate anything in the first place.

Looking at this problem, you might have the following idea: Although we don't know $P$, we do have the data $\zdata$, which consists of IID draws from $P$. From the law of large numbers, we know that expectations can be approximated by averages over IID draws, so we could approximate $R(\hat f \given \dD)$ by

```{math}
:enumerated: false

\frac{1}{N} \sum_{n=1}^N L(y_n, \hat f(\boldx_n))
```

where the pairs $(\boldx_n, y_n)$ are from the data set $\zdata$.

However, this is just the empirical risk, and the empirical risk is a highly biased estimator of the risk. This point was discussed extensively in {ref}`ss-ermls`. See, in particular, [figure %s](#f-rvempr). The point that figure made was that complex models tend to overfit, producing low empirical risk, but high risk. In essence, the problem is that we are using the data $\zdata$ twice, for conflicting objectives. First, we are using it to fit the model, producing $\hat f$. Second, we are using it to evaluate the predictive ability of $\hat f$ on new observations.

So what we really need is fresh data. New data will tell us how $\hat f$ performs out-of-sample. If we had $J$ new observations $(y^v_j, \boldx^v_j)$, then we could estimate the risk by

```{math}
:enumerated: false

\frac{1}{J} \sum_{j=1}^J L(y^v_j, \hat f(\boldx^v_j))
```

Of course, this isn't a genuine solution because we don't have any new data in general. One way to work around this problem is to take $\zdata$ and split it into two disjoint subsets, called the **training set** and the **validation set**. The training set is used to fit $\hat f$ and the validation set is used to estimate the risk of $\hat f$. We then repeat this for all models and choose the one with lowest estimated risk.

Since data are scarce, a more common procedure is **cross-validation**, which attempts to use the whole data set for both fitting the model and estimating the risk. To illustrate the idea, suppose that we partition the data set into two subsets $\dD_1$ and $\dD_2$. First, we use $\dD_1$ as the training set and $\dD_2$ as the validation set. Next we use $\dD_2$ as the training set, and $\dD_1$ as the validation set. The estimate of the risk is the average of the estimates of the risk produced in these two steps.

Of course, we could divide the data into more than two sets. The extreme is to partition the data into $N$ subsets. This procedure is called **leave-one-out cross-validation**. Letting $\dD_{-n} := \zdata \setminus \{ (\boldx_n, y_n) \}$, the data set with just the $n$th data point $(\boldx_n, y_n)$ omitted, the leave-one-out cross validation algorithm can be expressed as follows:

- **for** $n = 1,\ldots,N$ **do**
  - fit $\hat f_{-n}$ using data $\dD_{-n}$
  - set $r_n := L(y_n, \hat f_{-n}(\boldx_n))$
- **end**
- return the risk estimate $r := \frac{1}{N} \sum_{n=1}^N r_n$

At each step inside the loop, we fit the model using all but the $n$th data point, and then predict the $n$th data point using the fitted model. The prediction quality is evaluated in terms of loss. Repeating this $n$ times, we produce an estimate of the risk using average loss. On an intuitive level, the procedure is attractive because we are using the available data intensively but still evaluating based on out-of-sample error.

In terms of model selection, the idea is to run each model through the cross-validation procedure, and then select the one that produces the lowest value of $r$, the estimated risk. Let's illustrate this idea, by considering again the ridge regression procedure used in {ref}`ss-ssrr`. In this problem the set of models is indexed by $\lambda$, the regularization parameter in the ridge regression. The data set $\zdata$ is the set of points shown in [figures %s](#f-ridgeplot1)–[%s](#f-ridgeplot8). For each $\lambda$, the fitted function $\hat f_{\lambda}$ is

```{math}
:enumerated: false

\hat f_{\lambda}(x) = \hboldbeta_{\lambda}^\T \boldphi(x)
    \quad \text{where} \quad
    \hboldbeta_{\lambda} :=
    \argmin_{\boldb} \sum_{n=1}^N 
        \left\{
        (y_n - \boldb^\T\boldphi(x_n))^2 + \lambda \| \boldb \|^2
        \right\}
```

Recall here that $\boldphi(x) = (x^0, x^1,\ldots, x^d)$ with $d$ fixed at 14, so we are fitting a polynomial of degree 14 to the data by minimizing regularized least squares error. The amount of regularization is increasing in $\lambda$. The resulting functions $\hat f_{\lambda}$ were shown for different values of $\lambda$ in [figures %s](#f-ridgeplot1)–[%s](#f-ridgeplot8). Intermediate values of $\lambda$ produced the best fit in terms of minimizing risk (see [figures %s](#f-ridgeplot5) and [%s](#f-ridge_risk)).

In that discussion, we used the fact that we knew the underlying model to evaluate the risk, and hence the values of $\lambda$ that produce low risk ([figure %s](#f-ridge_risk)). In real estimation problems, risk is unobservable, and we need to choose $\lambda$ on the basis of the data alone (assuming we don't have prior knowledge, as in the Bayesian case—see {ref}`ss-abp`). Let's see how a data-based procedure such as cross-validation performs in terms of selecting a good value of $\lambda$.

In this experiment, for each $\lambda$ in the grid `exp(seq(-22, 10, length=10))` we perform leave-one-out cross-validation. The fit at each step within the loop is via ridge regression, omitting the $n$th data point, and the resulting polynomial is used to predict $y_n$ from $x_n$. The prediction error is measured by squared loss. In other words, for each $\lambda$ in the grid, we use the following algorithm to estimate the risk:

- **for** $n = 1,\ldots,N$ **do**
  - set $\hboldbeta_{\lambda, -n} := \argmin_{\boldb} \sum_{i \neq n} \{ (y_i - \boldb^\T\boldphi(x_i))^2 + \lambda \| \boldb \|^2 \}$
  - set $r_{\lambda, n} := (y_n - \hboldbeta_{\lambda, -n}^\T\boldphi(x_n))^2$
- **end**
- return $r_{\lambda} := \frac{1}{N} \sum_{n=1}^N r_{\lambda, n}$

The value of $\lambda$ producing the smallest estimated risk $r_{\lambda}$ is around 0.015. This is in fact very close to the value that minimizes the actual risk (see [figure %s](#f-ridge_risk)). The associated function $\hat f_{\lambda}$ is plotted in red in [figure %s](#f-cvbest), and indeed the fit is close. In this instance, our fully automated procedure is successful.[^5]

```{figure} figures/cvbest.svg
:name: f-cvbest

Fitted polynomial, $\lambda \approx 0.015$
```

## Further Reading

Textbook treatments of nonparametric methods in econometrics can be found in {cite}`pagan1999nonparametric`, {cite}`li2007nonparametric`, {cite}`ullah2012semiparametric`, {cite}`henderson2015applied` and {cite}`hansen2010econometrics`. {cite}`bosq2012nonparametric` covers nonparametrics for time series models. Recent studies in the econometrics literature using nonparametric or semiparametric methods include {cite}`chen2012testing`, {cite}`jeong2012consistent`, {cite}`henderson2012empirical`, {cite}`christensen2013nonparametric`, {cite}`canay2013testability`, {cite}`su2013nonparametric`, {cite}`bajari2013econometrics`, {cite}`newey2013nonparametric`, {cite}`hansen2014nonparametric`, {cite}`mastromarco2014effect`, {cite}`hickman2014replacing`, {cite}`matzkin2010estimation`, {cite}`bhattacharya2015nonparametric`, and {cite}`du2015nonparametric`. A beautiful overview of nonparametrics and the method of sieves can be found in {cite}`geman1982nonparametric`.

Good discussions of ridge regression and cross-validation can be found in {cite}`friedman2001elements` and {cite}`abu2012learning`. One alternative to ridge regression is the LASSO method ({cite:author}`tibshirani1996regression` {cite:year}`tibshirani1996regression`), where the term $\lambda \sum_{k=1}^K b_k^2$ in {eq}`eq-rrp` is replaced with an absolute deviation penalty $\lambda \sum_{k=1}^K |b_k|$. Also related is penalized maximum likelihood estimation. A recent application can be found in {cite}`shapiro2015measuring`.

An important contribution to the literature on subset selection is the least angle regression of {cite}`efron2004least`. A discussion of this method in relation to other subset techniques can be found in {cite}`hesterberg2008least`.

Recent papers treating regularization and related topics within the econometrics literature include {cite}`florens2003inverse`, {cite}`carrasco2007linear`, {cite}`chen2011rate`, {cite}`florens2013regularizing`, {cite}`darolles2011nonparametric`, {cite}`hoderlein2011demand`, {cite}`horowitz2014adaptive`, {cite}`hautsch2012blocking`, {cite}`li2014estimating`, {cite}`lunde2015econometric`, {cite}`chernozhukov2015post`, and {cite}`ando2015selecting`.

Regarding implementations, one powerful and increasingly popular package for Bayesian learning is Stan. See `http://mc-stan.org/`. A popular and high-quality routine for cross validation in the context of ridge regression can be found in the `sklearn.linear_model.RidgeCV` method in Python's `scikit-learn` package.

## Exercises

```{exercise}
:label: ex-l1wd

Let $\| \cdot \|_1$ denote the $L_1$ norm as defined in {eq}`eq-l1dist`. Show that $\| f - g \|_1$ is bounded over the set of all densities.[^6]
```


```{exercise}
:label: ex-pti

Prove the triangle inequality in the case of $L_1$ (see {eq}`eq-lptri`).
```


```{exercise}
:label: ex-ndl1

Prove the claim in [example %s](#eg-ndl1).[^7]
```


```{exercise}
:label: ex-vnpde

For the case $d=1$, show that $\hat f$ in {eq}`eq-npkde` is a density for every $N$, every $h > 0$ and every realization of the sample.[^8]
```


```{exercise}
:label: ex-sioc

Let $x$ be a finite random variable taking value $s_n$ in $\{s_n\}_{n=1}^N$ with probability $p_n$. Let $y$ be independent of $x$ with density $g$ on $\RR$. Show that $z = x+y$ has density $f(s') = \sum_{n=1}^N g(s' -s_n) p_n$.
```


```{exercise}
:label: ex-skh

Let $\boldx$ be multivariate Gaussian, with distribution $\nN(\boldmu, \boldSigma)$ and let $h > 0$. Suppose that $\boldSigma$ is nonsingular. Let $K$ be the density of $\boldx$ and let $K_h$ be the density of $h \boldx$. Show that $K$ and $K_h$ satisfy the relationship shown in {eq}`eq-kh`. (The result is true more generally but the Gaussian case is a nice exercise.)
```


````{exercise}
:label: ex-ch14_regularization-auto-1

Let $x_1, \ldots, x_N$ be IID scalar random variables with common density $f$ and let $\hat f$ be as defined in {eq}`eq-npkde`. Show that if $K(t) = \1\{ -1/2 < t < 1/2 \}$, then, for any number $s$, we have

```{math}
:enumerated: false

\EE [ \hat f(s) ] 
            = \frac{1}{h} \int_{s - h/2}^{s + h/2} f(t)\diff t
```
````


```{exercise}
:label: ex-dmse2

Verify the claim in {eq}`eq-dmse2`.
```


```{exercise}
:label: ex-rrbiased

By taking derivatives in {eq}`eq-rrp` to find stationary points, derive the expression for the ridge regression estimator $\hboldbeta_{\lambda}$. Show that $\hboldbeta_{\lambda}$ is a biased estimator of $\boldbeta$ whenever $\lambda > 0$.
```


```{exercise}
:label: ex-dnfr

In deriving $\hboldbeta_\lambda$ in {eq}`eq-rrs`, we do not require the full rank assumption (see [assumption %s](#a-fr)) whenever $\lambda > 0$. Explain why.
```


```{exercise}
:label: ex-betamap

Verify {eq}`eq-betamap` using {eq}`eq-priors` and {eq}`eq-betamap0`.
```


### Solutions to Selected Exercises

````{solution} ex-l1wd

Let $f$ and $g$ be any two densities on $\RR^J$. By the triangle inequality, for any given $\bolds \in \RR^J$, we have

```{math}
:enumerated: false

| f(\bolds) - g(\bolds) | 
        \leq | f(\bolds) | + | g(\bolds) | 
        =  f(\bolds)  +  g(\bolds)
```

Integrating both sides of this inequality gives $\| f - g \|_1 \leq 2$.
````


```{solution} ex-pti

Fix densities $f, g$ and $h$ on $\RR^d$. Pick any $\bolds \in \RR^d$. By the scalar triangle inequality, we have $|f(\bolds) - g(\bolds) | \leq |f(\bolds) - h(\bolds)| + |h(\bolds) - g(\bolds) |$. As per the hint in [exercise %s](#ex-l1wd), integration preserves this bound. This gives the triangle inequality in $L_1$.
```


````{solution} ex-ndl1

By Pinsker's inequality and [exercise %s](#ex-normkl), we have

```{math}
:enumerated: false

\| f(\cdot \,; \mu, \sigma) - f(\cdot \,; \bar x_N, s_N)  \|_1
        \leq 
        \sqrt{ \frac{\delta_N}{2} }
        \; \text{ where } 
        \delta_N :=
        \ln \frac{s_N}{\sigma}
        + \frac{\sigma^2 + (\mu - \bar x_N)^2}{2 s_N^2}
        - \frac{1}{2}
```

Since $\bar x_N \toprob \mu$ and $s_N \toprob \sigma$ (see {ref}`ss-ape`), applying the rules in [fact %s](#fa-reconpro), we have $\delta_N \toprob 0$. The claim follows.
````


````{solution} ex-vnpde

The nonnegativity of $\hat f$ is obvious. To show that $\int \hat f(s) \diff s = 1$, it is enough to show that $\int K \left( \frac{s - a}{h} \right) \diff s = h$ for any given number $a$. This equality can be obtained by the change of variable $u := (s - a) / h$, which leads to

```{math}
:enumerated: false

\int K \left( \frac{s - a}{h} \right) \diff s 
        = \int K (u) h \diff u
        = h \int K (u) \diff u
```

Since $K$ is a density, the proof is done.
````


````{solution} ex-sioc

Fix $s' \in \RR$. By the law of total probability we have

```{math}
:enumerated: false

\PP\{z \leq s'\}
         = \sum_{n=1}^N \PP\{x + y \leq s' \given x = s_n\}
            \PP\{x = s_n\}
         = \sum_{n=1}^N \PP\{y \leq s' - s_n \} p_n
```

Differentiating with respect to $s'$ gives the stated form for the density of $z$.
````


````{solution} ex-skh

Adopting the setting and notation of [ex. %s](#ex-skh), our aim is to show that $K$ and $K_h$ satisfy {eq}`eq-kh`. By our rules for linear transforms of Gaussians, the distribution of $h\boldx$ is $\nN(h \boldmu, h^2 \boldSigma)$. That is,

```{math}
:enumerated: false

K_h(\bolds) = (2 \pi)^{-d/2} \det(h^2 \boldSigma)^{-1/2} 
        \exp \left\{ 
            - \frac{1}{2} (\bolds - h\boldmu)^\T 
                (h^2 \boldSigma)^{-1} (\bolds - h \boldmu) 
        \right\}
```

Applying scalar multiple rules for inverses and determinants (see page  and page ), we can also write this as

```{math}
:enumerated: false

K_h(\bolds) =  \frac{1}{h^d} 
            (2 \pi)^{-d/2} \det(\boldSigma)^{-1/2} 
        \exp \left\{ 
            - \frac{1}{2} (\bolds/h - \boldmu)^\T 
                \boldSigma^{-1} (\bolds/h - \boldmu) 
        \right\}
```

That is, $K_h(\bolds) := \frac{1}{h^d} K ( \frac{\bolds}{h} )$, where $K$ is the density of $\nN(\boldmu, \boldSigma)$.
````


```{solution} ex-dnfr

The full rank assumption is not necessary because the matrix $\boldZ := \boldX^\T \boldX + \lambda \boldI$ is always invertible. To show this it suffices to show that $\boldZ$ is positive definite. This is not difficult to verify using the definition.
```


[^1]: For a proof, see p. 39 of {cite}`devroye2012combinatorial`.

[^2]: {numref}`f-one_dim_kde` is produced using `scikit-learn`. See `johnstachurski.net/emet` for code. In R, KDEs can be constructed from the `density` function. Try, for example, `plot(density(runif(200)))`.

[^3]: Full details on the choice of $\boldA$ are given in the code for this chapter at `johnstachurski.net/emet`.

[^4]: Remember {cite}`sala1997just` and his two million growth regressions?

[^5]: The code to run this experiment can be found at `johnstachurski.net/emet`.

[^6]: Hint: It is a basic property of integration that $g \leq h$ pointwise implies $\int g \leq \int h$.

[^7]: Hint: Use [fact %s](#fa-pinsk) and [exercise %s](#ex-normkl).

[^8]: Hint: Try a change-of-variables argument.
