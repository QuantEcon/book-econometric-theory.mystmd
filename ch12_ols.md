---
title: "Ordinary Least Squares"
label: c-tcolse
---

In [chapter 11](#c-reg) we imposed only mild regularity conditions on the data. To further investigate linear least squares estimation and provide additional interpretation of the estimated coefficients, we need to make more assumptions. The purpose of this chapter is to describe the properties of linear least squares estimation under the classical OLS assumptions, where OLS stands for **ordinary least squares**.

## Estimation under OLS

In this section we introduce the OLS assumptions and study the finite sample properties of the OLS estimators under these assumptions.

(ss-tolsa)=
### Assumptions
Let $\boldX$ and $\boldy$ be as defined in {ref}`s-olsfp`. The OLS assumptions are as follows:

```{prf:assumption}
:label: a-fr2

$\boldX$ has full column rank with probability one.
```


```{prf:assumption}
:label: a-lols

The observations satisfy $\boldy = \boldX \boldbeta + \boldu$ for some unknown $K$-vector of parameters $\boldbeta$ and unobservable vector of shocks $\boldu$.
```


```{prf:assumption}
:label: a-bols

$\EE [ \boldu \given \boldX ] = \boldzero$.
```


```{prf:assumption}
:label: a-bols2

$\EE [ \boldu \boldu^\T \given \boldX ] = \sigma^2 \boldI$ for some unknown $\sigma > 0$.
```


Let's break these assumptions down and see what they imply.

{prf:ref}`a-fr2` is just a repetition of [assumption %s](#a-fr) from [chapter 11](#c-reg). {prf:ref}`a-lols` can be decomposed into the separate equations

$$
y_n = \boldx_n^\T \, \boldbeta + u_n,
    \qquad n = 1, \ldots, N
$$ (eq-lm)

Here are some examples of models that produce relationships in the form of {eq}`eq-lm`.

````{prf:example}
:label: eg-cdprod

The **Cobb–Douglas production function** relates capital and labor inputs with output via $y = A k^a \ell^{\delta}$, where $A$ is a random, firm-specific productivity term and $a$ and $\delta$ are parameters. Taking logs yields the linear regression model

```{math}
:enumerated: false

\ln y = \gamma + a \ln k + \delta \ln \ell + u
```

where the random term $\ln A$ is represented by $\gamma + u$.
````


```{prf:example}

The **gravity model** relates international trade flows between country $\ell$ and country $n$ via the equation $T_{\ell n} = \lambda \xi_{\ell n} G_\ell^{\alpha} G_n^\beta /D_{\ell n}^{\gamma}$. Here $T_{\ell n}$ is exports from country $\ell$ to country $n$, $\lambda$ is a constant term, $\xi_{\ell n}$ is a shock, $G_\ell$ and $G_n$ are GDP in country $\ell$ and $n$ respectively, and $D_{\ell n}$ is distance between them. Taking logs gives

$$
\ln T_{\ell n} 
        = \ln \lambda 
        +  \alpha \ln G_\ell + \beta \ln G_n - \gamma \ln D_{\ell n}
        + \ln \xi_{\ell n}
$$ (eq-gravity)
For motivation and discussion, see, for example, {cite}`van2010gravity`.
```


```{prf:property}
:label: fa-exfoss

If the linearity [assumption %s](#a-lols) holds, then

(i) $\boldM \boldy = \boldM \boldu$.

(ii) $\boldP \boldy = \boldX \boldbeta + \boldP \boldu$.

(iii) $\rss = \boldu^\T \boldM \boldu$.
```


The notation in [fact %s](#fa-exfoss) follows {ref}`ss-stnot`. {prf:ref}`ex-exfoss` requests a proof.

```{prf:property}
:label: fa-sruc

If [assumption %s](#a-bols) holds, then

(i) $\EE [\boldu] = \boldzero$,

(ii) $\EE [ u_m \given x_{nk} ] = 0$ for any $m, n, k$,

(iii) $\EE [ u_m x_{nk} ] = 0$ for any $m, n, k$ (orthogonality), and

(iv) $\cov [ u_m , x_{nk} ] = 0$ for any $m, n, k$.
```


{prf:ref}`ex-psruc` asks you to give proofs.

```{prf:property}
:label: fa-vciu

If [assumption %s](#a-bols2) holds, then

(i) $\var [\boldu] = \EE [ \boldu \boldu^\T ] = \sigma^2 \boldI$,

(ii) $\EE [ u_i^2 \given \boldX ] = \EE [ u_j^2 \given \boldX ] = \sigma^2$ for any $i, j$ in $1,\ldots,N$, and

(iii) $\EE [ u_i u_j \given \boldX ] = 0$ whenever $i \neq j$.
```


Parts (ii) and (iii) are called **homoskedasticity** and zero correlation respectively.

Combining [assumptions %s](#a-bols) and [%s](#a-bols2) gives

```{math}
:enumerated: false

\var [ \boldu  \given \boldX ] 
    := \EE [ \boldu \boldu^\T \given \boldX ]
            -  \EE [ \boldu \given \boldX ] \EE [ \boldu^\T \given \boldX ]
    = \EE [ \boldu \boldu^\T \given \boldX ]
```

Thus the conditional variance–covariance matrix is the diagonal matrix $\sigma^2 \boldI$.

(ss-toest)=
### The OLS Estimators
The standard estimator of $\boldbeta$ in {eq}`eq-lm` is the least squares estimator $\hboldbeta$ defined in {eq}`eq-thelsqe`. That is,

$$
\hboldbeta := (\boldX^\T \boldX)^{-1} \boldX^\T \boldy
$$ (eq-thelsqe2)

Unlike [chapter 11](#c-reg), it's now understood as an estimator of the unknown parameter vector $\boldbeta$. In this context, $\hboldbeta$ is also called the **OLS estimator** of $\boldbeta$. Substituting $\boldy = \boldX \boldbeta + \boldu$ into {eq}`eq-thelsqe2` and cancelling terms gives the useful comparison

$$
\hboldbeta = 
    \boldbeta + (\boldX^\T \boldX)^{-1} \boldX^\T \boldu
$$ (eq-hg)

The usual OLS estimator of the parameter $\sigma^2$ introduced in [assumption %s](#a-bols2) is

$$
\hat \sigma^2 := \frac{\rss}{N - K}
$$ (eq-hatsig)

While we can compute $\hboldbeta$ and $\hat \sigma^2$ directly using the expressions in {eq}`eq-thelsqe2` and {eq}`eq-hatsig`, in practice it's more common to pass them to a routine in a regression package. Here's an example that runs the gravity model regression {eq}`eq-gravity` on world trade data using Python.

First we import some Python libraries often used for statistics:

``` python
from numpy import log  # used in the formula below
import pandas as pd, statsmodels.formula.api as smf
```

Then we read in the data to a `pandas DataFrame` from a local CSV file:

``` python
data = pd.read_csv("./data/gravity_dataset_2013.csv")
```

Next we build a model using a formula to indicate the regression we want to run symbolically (similar to `R`).

``` python
formula = "log(value) ~ log(egdp) + log(igdp) + log(dist)"
model = smf.ols(formula, data)
```

The formula used above is analogous to {eq}`eq-gravity`. We omit the shock (the last term in {eq}`eq-gravity`) because it is not part of the data set, as well as the constant term, which is included automatically. The names `egdp` and `igdp` represent exporter and importer GDP respectively.

Now we can perform the estimation and print a table summarizing results:

``` python
result = model.fit(cov_type='HC1')
print(result.summary())
```

The output is as follows:


                                OLS Regression Results                            
    ==============================================================================
    Dep. Variable:             log(value)   R-squared:                       0.652
    Model:                            OLS   Adj. R-squared:                  0.652
    Method:                 Least Squares   F-statistic:                 1.203e+04
    Date:                Wed, 04 Nov 2015   Prob (F-statistic):               0.00
    Time:                        16:02:42   Log-Likelihood:                -47185.
    No. Observations:               19655   AIC:                         9.438e+04
    Df Residuals:                   19651   BIC:                         9.441e+04
    Df Model:                           3                                         
    Covariance Type:                  HC1                                         
    ==============================================================================
                     coef    std err          z      P>|z|      [95.0% Conf. Int.]
    ------------------------------------------------------------------------------
    Intercept    -30.2350      0.394    -76.773      0.000       -31.007   -29.463
    log(egdp)      1.2783      0.008    153.772      0.000         1.262     1.295
    log(igdp)      1.0287      0.009    118.885      0.000         1.012     1.046
    log(dist)     -1.3483      0.023    -58.113      0.000        -1.394    -1.303
    ==============================================================================
    Omnibus:                     1819.031   Durbin-Watson:                   1.798
    Prob(Omnibus):                  0.000   Jarque-Bera (JB):             3361.159
    Skew:                          -0.639   Prob(JB):                         0.00
    Kurtosis:                       4.571   Cond. No.                         708.
    ==============================================================================

The value of the OLS estimate $\hboldbeta$ is the column under `coef`. Objects such as standard errors are discussed below. The argument `cov_type=’HC1’` to the `model.fit` call that produced this table requests heteroskedasticity-consistent or "robust\" standard errors, which will be discussed in {ref}`ss-heter`.

Let's check that the values for $\hboldbeta$ reported under `coef` coincide with the expression for $\hboldbeta$ given in {eq}`eq-thelsqe2`. First we build the design matrix $\boldX$ (see the URL given above for details) and then compute as follows:

    betahat = inv(X.T @ X) @ X.T @ y

When we print the result the output agrees with the output in the table

    [-30.23498073   1.27825004   1.02865139  -1.34830012]    

The expression for the $R^2$ in the table above is for $R_c^2$, as given in {eq}`eq-cecrs`.

{numref}`f-partial_reg_plot` shows a **partial regression plot** for the linear model just described. It is produced from the `result` object generated above via the commands

``` python
import matplotlib.pyplot as plt
import statsmodels.api as sm
fig = plt.figure()
fig = sm.graphics.plot_partregress_grid(result, fig=fig)
```

To understand the figure, consider, for example, the bottom left subplot, which shows the relationship between `log(igdp)` and the dependent variable. The data points in the plot correspond to the residuals in {eq}`eq-resregy` and {eq}`eq-resregx`. In particular, the values on the horizontal axis are the residuals from regressing `log(igdp)` on all other columns in $\boldX$, while those on the vertical axis are the residuals from regressing the dependent variable `log(value)` on all other columns in $\boldX$.

The relationship between these residuals measures the change in `log(value)` associated with changing `log(igdp)` while holding all other variables fixed. The slope of the regression line running through the plot is exactly the coefficient estimate on `log(igdp)`. For clarification see the discussion surrounding {eq}`eq-resregy` and {eq}`eq-resregx`.

The partial regression plot can be useful both for interpretation and diagnostics. For now, however, let's look at some general theory concerning properties of the OLS estimators.

```{figure} figures/partial_reg_plot.png
:name: f-partial_reg_plot

Partial regression plot for the gravity model
```

### Finite Sample Properties

Let's start with classical results on finite sample properties under the OLS assumptions. These are elegant results, although you should also read {ref}`s-modsel` if you want to view them with proper perspective. Our first result shows that, under the OLS assumptions, the OLS estimators are unbiased.

```{prf:theorem}
:label: t-olsiu

If [assumptions %s](#a-fr2)–[%s](#a-bols) hold, then

$$
\EE \hboldbeta = \EE [\hboldbeta \given \boldX] = \boldbeta
$$ (eq-hbub)
If [assumption %s](#a-bols2) also holds, then

$$
\EE  \hat \sigma^2  = \EE [ \hat \sigma^2 \given \boldX ] = \sigma^2
$$ (eq-sbub)
```


````{prf:proof}
:nonumber:

By {eq}`eq-hg` and [assumption %s](#a-bols), we have $\EE [ \hboldbeta \given \boldX ] = \boldbeta + (\boldX^\T \boldX)^{-1} \boldX^\T \EE [ \boldu \given \boldX ] = \boldbeta$. Taking the unconditional expectation gives {eq}`eq-hbub`.

Regarding {eq}`eq-sbub`, if we use the expression for $\rss$ in [fact %s](#fa-exfoss), we get

```{math}
:enumerated: false

\EE [\rss \given \boldX]
        = \EE [\boldu^\T \boldM \boldu \given \boldX]
        = \sigma^2 \EE [\boldxi^\T \boldM \boldxi \given \boldX]
        \quad \text{where} \quad
        \boldxi := \sigma^{-1} \boldu
```

It follows now from [fact %s](#fa-eoqf) that $\EE [\rss \given \boldX] = \sigma^2 \trace \boldM$. Hence

```{math}
:enumerated: false

\EE [ \hat \sigma^2 \given \boldX ]
        = \frac{\sigma^2 \trace(\boldM)}{N - K}
```

By [fact %s](#fa-trop) we have $\trace \boldM = N - K$. Inserting this into the preceding equation and taking unconditional expectations gives {eq}`eq-sbub`. ◻
````


Now that we know $\hboldbeta$ is unbiased, let's consider the variance.

```{prf:theorem}
:label: t-olsev

If [assumptions %s](#a-fr2)–[%s](#a-bols2) hold, then

$$
\var [\hboldbeta \given \boldX ] = \sigma^2 (\boldX^\T \boldX)^{-1}
$$ (eq-vhbb)
```


For a proof, see [exercise %s](#ex-olsev) and its solution.

Is the value on the right-hand side of {eq}`eq-vhbb` small or large? A famous result in OLS theory is that under the stated conditions, the OLS estimator $\hboldbeta$ is best linear unbiased. That is, in the class of linear unbiased estimators of $\boldbeta$, the OLS estimator has the lowest variance.

````{prf:theorem}
:label: t-gm

Let [assumptions %s](#a-fr2)–[%s](#a-bols2) hold and let $\boldb$ be an estimator of $\boldbeta$. If $\boldb$ is linear and unbiased, then

```{math}
:enumerated: false

\var [\hboldbeta \given \boldX] \leq \var [\boldb \given \boldX]
```
````


Here $\var [\hboldbeta \given \boldX] \leq \var [\boldb \given \boldX]$ means that $\var [\boldb \given \boldX] - \var [\hboldbeta \given \boldX]$ is nonnegative definite. Since matrices have no standard ordering, this is one way to assert that $\var [\boldb \given \boldX]$ is "larger." One implication is that $\var [b_k \given \boldX] \geq \var [\hat \beta_k \given \boldX]$ for all $k$. (Why?)

Also, linearity of $\boldb$ means that $\boldb$ is linear as a function of $\boldy$, taking $\boldX$ as fixed. In view of [theorem %s](#t-lmaeq), this is equivalent to requiring that $\boldb = \boldC \boldy$ for some matrix $\boldC$. The matrix $\boldC$ is allowed to depend on $\boldX$ but not $\boldy$.

The statement that $\boldb$ is unbiased means $\EE [\boldb \given \boldX] = \EE [ \boldC \boldy \given \boldX] = \boldbeta$ for all $\boldbeta \in \RR^K$.

````{prf:proof}
:nonumber:

**Proof of [theorem %s](#t-gm).**

Let $\boldb = \boldC \boldy$, as described above, and let $\boldD := \boldC - \boldA$, where $\boldA := (\boldX^\T\boldX)^{-1} \boldX^\T$. Then

$$
\boldb = \boldC \boldy = \boldD \boldy + \boldA \boldy
        = \boldD (\boldX \boldbeta + \boldu) + \hboldbeta
        = \boldD \boldX \boldbeta + \boldD \boldu + \hboldbeta
$$ (eq-bdbb)
Taking conditional expectations and using the fact that $\boldD$ is a function of $\boldX$, we obtain

```{math}
:enumerated: false

\begin{aligned}
\EE [ \boldb \given \boldX] 
        & = \EE[ \boldD \boldX \boldbeta \given \boldX]  + 
            \EE[ \boldD \boldu \given \boldX ] 
            + \EE[ \hboldbeta \given \boldX] 
            \\
        & = \boldD \boldX \, \EE[ \boldbeta \given \boldX]  + 
             \boldD  \,\EE[\boldu \given \boldX ] 
            + \EE[ \hboldbeta \given \boldX] 
            \\
        & = \boldD \boldX \boldbeta + \boldzero + \boldbeta
\end{aligned}
```

Since $\boldb$ is unbiased and, in particular, $\EE [ \boldb \given \boldX] = \boldbeta$ for any given $\boldbeta$, we have

```{math}
:enumerated: false

\boldbeta = \boldD \boldX \boldbeta + \boldbeta
        \quad \text{ for all } \boldbeta \in \RR^K
```

```{math}
:enumerated: false

\fore
        \boldzero = \boldD \boldX \boldbeta 
        \quad \text{ for all } \boldbeta \in \RR^K
```

In light of [exercise %s](#ex-ictz), we conclude that $\boldD \boldX = \boldzero$. Combining this result with {eq}`eq-bdbb`, we obtain $\boldb = \boldD \boldu + \hboldbeta$. Hence $\boldb$ is equal to the OLS estimator plus zero-mean noise.

To complete the proof, observe that

```{math}
:enumerated: false

\var [ \boldb \given \boldX]
        = \var[ \boldD \boldu + \hboldbeta \given \boldX ]
        = \var[ (\boldD + \boldA ) \boldu \given \boldX ] 
        = (\boldD + \boldA ) \var[ \boldu \given \boldX ] (\boldD + \boldA )^\T
```

Using [assumption %s](#a-bols2) and [fact %s](#fa-owtra), the right-hand side of this expression becomes

```{math}
:enumerated: false

\sigma^2 (\boldD + \boldA )(\boldD^\T + \boldA^\T)
        = \sigma^2 (\boldD \boldD^\T + \boldD \boldA^\T + \boldA \boldD^\T + \boldA \boldA^\T)
```

Since

```{math}
:enumerated: false

\boldD \boldA^\T = \boldD \boldX (\boldX^\T \boldX)^{-1} 
        = \boldzero (\boldX^\T \boldX)^{-1}
        = \boldzero
```

and

```{math}
:enumerated: false

\boldA \boldA^\T
        =  (\boldX^\T \boldX)^{-1} \boldX^\T \boldX (\boldX^\T \boldX)^{-1} 
        = (\boldX^\T \boldX)^{-1}
```

we conclude that

```{math}
:enumerated: false

\var [ \boldb \given \boldX]
        = \sigma^2 [\boldD \boldD^\T + (\boldX^\T \boldX)^{-1}]
        = \sigma^2 \boldD \boldD^\T + \var [ \hboldbeta \given \boldX]
```

The matrix $\sigma^2 \boldD \boldD^\T$ is nonnegative definite, so the proof is now done. ◻
````


(ss-polse)=
#### Precision of the Estimates
{prf:ref}`t-olsev` shows that the variance–covariance matrix of $\hboldbeta$ given $\boldX$ is $\sigma^2 (\boldX^\T \boldX)^{-1}$. The scalar variances of the individual OLS coefficient estimates $\hat \beta_1, \ldots, \hat \beta_K$ are given by the principal diagonal of this matrix. Since each $\hat \beta_k$ is unbiased ([theorem %s](#t-olsiu)), small variance means probability mass concentrated around the true parameter $\beta_k$. In this case, we say that the estimator has high **precision**.

The Gauss–Markov theorem states that, as far as unbiased linear estimators go, the OLS estimates will have highest precision. But what if we hold the estimation technique fixed, as well as the sample size, and vary the application? Which problems will have high precision estimates, and which will have low precision estimates?

To answer this question, let's focus on the variance of a fixed coefficient $\beta_k$. We can write $\boldy = \boldX \boldbeta + \boldu$ as

$$
\boldy = \boldX_1 \boldbeta_1 + \col_k (\boldX) \beta_k + \boldu
$$ (eq-cfrm)

where $\col_k(\boldX)$ is the vector of observations of the $k$th regressor, $\boldX_1$ contains as its columns the observations of the other regressors, and $\hboldbeta_1$ is the OLS estimates of the corresponding coefficients. By the FWL theorem,

$$
\hat \beta_k 
    = (\col_k(\boldX)^\T \boldM_1 \col_k(\boldX))^{-1} \col_k(\boldX)^\T \boldM_1 \boldy
$$ (eq-exir)

where $\boldM_1$ is the residual projection $\boldM_1 := \boldI - \boldP_1$ and $\boldP_1 := \boldX_1 (\boldX_1^\T\boldX_1)^{-1} \boldX_1^\T$ projects onto $\colspace \boldX_1$. Applying $\boldM_1$ to both sides of {eq}`eq-cfrm` yields

```{math}
:enumerated: false

\boldM_1 \boldy =  \boldM_1 \col_k(\boldX) \beta_k + \boldM_1 \boldu
```

Substituting into {eq}`eq-exir` gives

$$
\hat \beta_k 
    = \beta_k + (\col_k(\boldX)^\T \boldM_1 \col_k(\boldX))^{-1}
    \col_k(\boldX)^\T \boldM_1 \boldu
$$ (eq-exir2)

Some calculations then show ([ex. %s](#ex-voir)) that

$$
\var[ \hat \beta_k \given \boldX ]
    = \sigma^2 
    (\col_k(\boldX)^\T \boldM_1 \col_k(\boldX))^{-1}
    = \sigma^2 \|\boldM_1 \col_k(\boldX)\|^{-2}
$$ (eq-voir)

Thus the variance of $\hat \beta_k$ depends on two components, the variance $\sigma^2$ of the shock $u$ and the norm of the vector $\boldM_1 \col_k \boldX$.

The variance in $\sigma^2$ is unavoidable: some data are noisier than others. The larger the variance in the unobservable shock, the harder it will be to estimate $\beta_k$ with good precession.

The other term is more interesting. $\boldM_1 \col_k \boldX$ is the residuals from regressing $\col_k \boldX$ on $\boldX_1$, and $\|\boldM_1 \col_k \boldX\|$ is the norm of this vector. If this norm is small, then the variance of $\hat \beta_k$ will be large.

When will this norm be small? It will be so when $\col_k \boldX$ is "almost" a linear combination of the other regressors and hence close to $\colspace \boldX_1$. Then

```{math}
:enumerated: false

\| \boldM_1 \col_k \boldX \| =
        \| \col_k \boldX  - \boldP_1 \col_k \boldX  \| \approx 0
```

This situation is sometimes referred to as **multicollinearity**.

{numref}`f-betahat_var` gives an illustration of the effect of multicollinearity on the variance of OLS estimate $\hat \beta_2$. In this example $\boldX$ has two columns. They are related by $\col_2 \boldX = \delta \col_1 \boldX + (1 - \delta) \boldz$, where $\boldz$ is a vector of $N$ independent draws from the standard normal distribution. Larger $\delta$ means more dependence between $\col_1 \boldX$ and $\col_2 \boldX$. The true parameter is $\beta_2 = 1$. The figure shows the distribution of $\hat \beta_2$ for different values of $\delta$. The estimator $\hat \beta_2$ is unbiased but its variance increases with $\delta$. (The normal shape of the distribution is explained in the next section.)

```{figure} figures/betahat_var.svg
:name: f-betahat_var

Distribution of $\hat \beta_2$ when $\col_2 \boldX = \delta \col_1 \boldX + (1 - \delta) \boldz$
```

(ss-ne)=
### Inference with Normal Errors
If we want to compute confidence intervals or test hypotheses about the coefficients in finite samples (i.e., without appealing to asymptotics), we need to strengthen our OLS assumptions even further, by specifying the parametric class of the error vector $\boldu$. Because of its many attractive properties, the normal distribution is the default.

```{prf:assumption}
:label: a-bolsg

$\boldX$ and $\boldu$ are independent and $\lL(\boldu) = \nN(\boldzero, \sigma^2 \boldI)$.
```


{prf:ref}`a-bolsg` implies both [assumption %s](#a-bols) and [assumption %s](#a-bols2). It also implies that the conditional distribution of $\hboldbeta$ given $\boldX$ is normal, as can be seen from {eq}`eq-hg` and [fact %s](#fa-nipre). The next theorem records this result.

````{prf:theorem}
:label: t-dools

If [assumptions %s](#a-fr2)–[%s](#a-lols) and [%s](#a-bolsg) hold, then

```{math}
:enumerated: false

\lL \left[ \, \hboldbeta \given \boldX \, \right] 
        = \nN(\boldbeta, \sigma^2 (\boldX^\T \boldX)^{-1})
```
````


It follows from [theorem %s](#t-dools) and {eq}`eq-mntm` that the distribution of individual coefficient $\hat \beta_k$ given $\boldX$ is also normal, with

$$
\lL \left[ \, \hat \beta_k \given \boldX \, \right] 
    = \lL \left[ \, \bolde_k^\T \, \hboldbeta  \given \boldX \, \right] 
    = \nN(\beta_k, \sigma^2 v_k(\boldX))
$$ (eq-dicof)

Here $v_k (\boldX) :=$ the $(k,k)$th element of $(\boldX^\T \boldX)^{-1}$.

{prf:ref}`t-dools` and the normality of the marginal distributions along the individual coefficients is illustrated in [figure %s](#f-dist_of_betahat). The central plot shows the contours of a density estimate that can be understood as a smoothed histogram based on simulated observations of $\hboldbeta$ (see {ref}`s-panp` for details on these density estimates). Each observation $\hboldbeta = (\hat \beta_1, \hat \beta_2)$ is obtained by fitting the regression model $y_n = \beta_1 + \beta_2 x_n + u_n$, where $\beta_1 = 0.5$, $\beta_2 = 1$ and $u_n$ is standard normal. The sample size for each fit is $N=250$. The univariate distributions on the top and right are univariate density estimates based on the observations of $\hat \beta_1$ and $\hat \beta_2$ respectively.

```{figure} figures/dist_of_betahat.svg
:name: f-dist_of_betahat

Simulated observations of $\hboldbeta$ and its marginals
```

We also want to know the distribution of $\hat \sigma^2$, as defined in {eq}`eq-hatsig`. It turns out to be more convenient to work with a transformation of $\hat \sigma^2$ defined as

$$
Q := \frac{\rss}{\sigma^2}
    = (N-K) \frac{\hat \sigma^2}{\sigma^2}
$$ (eq-Q)

````{prf:theorem}
:label: t-dools2

If [assumptions %s](#a-fr2)–[%s](#a-lols) and [%s](#a-bolsg) hold, then

```{math}
:enumerated: false

\lL \left[ \, Q \given \boldX \, \right] = \chi^2(N - K)
```
````


````{prf:proof}
:nonumber:

Recalling [fact %s](#fa-exfoss), we have

```{math}
:enumerated: false

Q 
        = \frac{\boldu^\T \boldM \boldu}{\sigma^2}
        = (\sigma^{-1} \boldu)^\T \boldM (\sigma^{-1} \boldu)
```

The claim now follows from [fact %s](#fa-mcfi) and the fact that $\trace \boldM = N - K$. (The last equality was discussed in the proof of [theorem %s](#t-olsiu).) ◻
````


```{prf:property}
:label: fa-ioes

If [assumptions %s](#a-fr2)–[%s](#a-lols) and [%s](#a-bolsg) hold, then the random elements $\hat \sigma^2$ and $\hboldbeta$ are independent given $\boldX$ (see [ex. %s](#ex-ioes)).
```


(ss-ttt)=
#### The $t$-Test
Let's consider the problem of testing a hypothesis about an individual coefficient $\beta_k$. Specifically, we consider the null hypothesis

```{math}
:enumerated: false

H_0 \colon \beta_k = \beta_k^0
```

where $\beta_k^0$ is any number. If we knew $\sigma^2$, we could construct a test of $H_0$ based around the following observation:

$$
z_k := \frac{\hat \beta_k - \beta_k} { \sigma \sqrt{v_k(\boldX)} }
    \; \implies \;
    \lL [ \, z_k \given \boldX \, ]
        = \nN(0, 1)
$$ (eq-istn)

(The implication follows from {eq}`eq-dicof`.) Since we don't, the standard methodology is to replace $\sigma^2$ with its estimator $\hat \sigma^2$, and determine the distribution of the resulting test statistic. Our next result implements this idea. In doing so, we will make use of the following notation:

```{math}
:enumerated: false

\se(\hat \beta_k) := \sqrt{ \hat \sigma^2 v_k(\boldX)}
```

The term $\se(\hat \beta_k)$ is called the **standard error** of $\hat \beta_k$. It can be regarded as the sample estimate of the standard deviation of $\hat \beta_k$. Replacing this standard deviation with its sample estimate $\se(\hat \beta_k)$ and $\beta_k$ with $\beta_k^0$ in {eq}`eq-istn`, we obtain the **$t$-statistic**

$$
t_k := \frac{\hat \beta_k - \beta_k^0}{\se(\hat \beta_k)}
$$ (eq-tstat)

associated with the hypothesis $H_0$. The distribution of this statistic under the null is described in the next theorem.

````{prf:theorem}
:label: t-ttest

Let [assumptions %s](#a-fr2)–[%s](#a-lols) and [%s](#a-bolsg) hold. If $H_0$ is true, then

```{math}
:enumerated: false

\lL \left[ \, t_k \given \boldX \, \right] = 
        \text{Student's t with } N-K \text{ degrees of freedom}
```
````


````{prf:proof}
:nonumber:

Taking $\boldX$ as fixed, we need to show that $t_k$ can be written as

$$
t_k = z_k \sqrt{ \frac{N-K}{Q} }
$$ (eq-nkqr)
where $z_k$ is standard normal, $Q$ is $\chi^2(N-K)$ and $z_k$ and $Q$ are independent. The desired result then follows from [fact %s](#fa-astudt).

To establish these properties, observe that, under $H_0$,

```{math}
:enumerated: false

t_k 
         = \frac{\hat \beta_k - \beta_k^0}
         {\sqrt{ \hat \sigma^2 v_k(\boldX)}}
         = \frac{\hat \beta_k - \beta_k}
         {\sqrt{ \sigma^2 v_k(\boldX) }}
            \sqrt{\frac{\sigma^2}{\hat \sigma^2}}
```

In view of the definition of $z_k$ in {eq}`eq-istn`, the definition of $Q$ in {eq}`eq-Q` and [theorem %s](#t-dools2), it remains only to show that $z_k$ and $Q$ are independent given $\boldX$. Since we can write $z_k$ as a function of $\hboldbeta$ and $Q$ as a function of $\hat \sigma^2$, this follows from [fact %s](#fa-ioes). ◻
````


Let $T := |t_k|$ and let a desired size $\alpha$ for our test of $H_0$ be given. In view of {eq}`eq-solfc`, to generate a test of size $\alpha$, we choose $c = c_{\alpha}$ to solve $\alpha = \PP_{\theta}\{T > c\}$, or

```{math}
:enumerated: false

1 - \alpha = \PP_{\theta}\{|t_k| \leq c\}
```

From {eq}`eq-symcdf`, the solution is $c_{\alpha} = F^{-1}(1 - \alpha/2)$, where $F$ is the Student's $t$ CDF with $N-K$ degrees of freedom. In view of [example %s](#eg-pvsc), the corresponding $p$-value is $2 F(- |t_k|)$.

````{prf:example}

A common implementation of the $t$-test is the test that a given coefficient is equal to zero. For the $k$th coefficient $\beta_k$, this leads to the statistic

```{math}
:enumerated: false

t_k := \frac{\hat \beta_k}{\se(\hat \beta_k)}
```

This statistic is sometimes called the **Z-score**. An application with simulated data is given in [listing %s](#list-tstats). Running the code shows that the Z-scores calculated by the `zscore` function agree with the "t value" column of the summary table produced by R's `summary` function (the last line of [listing %s](#list-tstats), output not shown). It is left as an exercise to check that the $p$-values in the same table agree with the formula $2 F(- |t_k|)$ given above.
````


```{code-block} r
:name: list-tstats
:caption: Calculating Z-scores (R)
:linenos:


set.seed(1234) 
N <- 50; K <- 3
beta <- rep(1, K)
X <- cbind(runif(N), runif(N), runif(N))
u <- rnorm(N) 
y <- X %*% beta + u

betahat <- solve(t(X) %*% X) %*% t(X) %*% y 
residuals <- y - X %*% betahat
sigmahat <- sqrt(sum(residuals^2) / (N - K))

# Compute t-stat (Z-score) for k-th regressor
zscore <- function(k) {
    se <- sigmahat * sqrt(solve(t(X) %*% X)[k, k])
    return(betahat[k] / se)
}
# Print t-stats
for (k in 1:3) {
    cat("t-stat, k =", k, ":", zscore(k), "\n")
}
# For comparison:
print(summary(lm(y ~ X - 1)))

```

#### The $F$-Test

The $t$-test is used to test hypotheses about individual regressors. For hypotheses concerning multiple regressors, the most common test is the $F$-test. The $F$-test can test quite general hypotheses. For simplicity we will focus on null hypotheses that restrict a subset of the coefficients to be zero.

Let $\boldX_1$, $\boldX_2$, $\hboldbeta_1$, $\hboldbeta_2$, $\boldP_1$, and $\boldM_1$ be as defined in {ref}`ss-sap`. As our null hypothesis we take

```{math}
:enumerated: false

H_0 \colon \boldbeta_2 = \boldzero
```

Since

$$
\boldy 
    = \boldX \boldbeta + \boldu 
    = \boldX_1 \boldbeta_1 + \boldX_2 \boldbeta_2 + \boldu
$$ (eq-altr)

it follows that under the null hypothesis we have

$$
\boldy = \boldX_1 \boldbeta_1 + \boldu
$$ (eq-nullr)

Letting

```{math}
:enumerated: false

\ussr := \| \boldM \boldy \|^2
    \quad \text{and} \quad
    \rssr := \| \boldM_1 \boldy \|^2
```

be the residual sums of squares for the unrestricted regression {eq}`eq-altr` and restricted regression {eq}`eq-nullr` respectively, the standard test statistic for our null hypothesis is

$$
F := \frac{(\rssr - \ussr) / K_2}{ \ussr / (N - K)}
$$ (eq-fstatr)

Large residuals in the restricted regression {eq}`eq-nullr` relative to those in {eq}`eq-altr` result in large values for $F$, which translates to evidence against the null hypothesis.

```{prf:theorem}
:label: t-ftest

Let [assumptions %s](#a-fr2)–[%s](#a-lols) and [%s](#a-bolsg) hold. If $H_0$ is true, then, given $\boldX$, the statistic $F$ defined in {eq}`eq-fstatr` has the $F$ distribution with parameters $(K_2, N - K)$.
```


````{prf:proof}
:nonumber:

Let $Q_1 := (\rssr - \ussr)/ \sigma^2$ and let $Q_2 := \ussr/\sigma^2$, so that

```{math}
:enumerated: false

F = \frac{Q_1 / K_2}{ Q_2 / (N - K)}
```

In view of [fact %s](#fa-rcif), it now suffices to show that, under the null hypothesis,

(a) $Q_1$ is chi-squared with $K_2$ degrees of freedom,

(b) $Q_2$ is chi-squared with $N - K$ degrees of freedom, and

(c) $Q_1$ and $Q_2$ are independent.

Part (b) was established in [theorem %s](#t-dools2). Regarding part (a), observe that, under the null hypothesis,

- $\ussr = \| \boldM \boldy \|^2 = \| \boldM (\boldX_1 \boldbeta_1 + \boldu) \|^2 = \| \boldM \boldu \|^2 = \boldu^\T \boldM \boldu$ and

- $\rssr = \| \boldM_1 \boldy \|^2 = \| \boldM_1 (\boldX_1 \boldbeta_1 + \boldu)\|^2 = \| \boldM_1 \boldu \|^2 = \boldu^\T \boldM_1 \boldu$.

It follows that

```{math}
:enumerated: false

\rssr - \ussr 
        = \boldu^\T \boldM_1 \boldu  - \boldu^\T \boldM \boldu 
        = \boldu^\T (\boldM_1 - \boldM) \boldu
```

Using the definitions of $\boldM$ and $\boldM_1$, we obtain

```{math}
:enumerated: false

Q_1 =
        \frac{\rssr - \ussr}{\sigma^2}
        = \frac{\boldu^\T (\boldI - \boldP_1 - \boldI + \boldP) \boldu}{\sigma^2}
        = (\sigma^{-1} \boldu)^\T (\boldP - \boldP_1) (\sigma^{-1} \boldu)
```

It is an exercise to show that $(\boldP - \boldP_1)$ is symmetric and idempotent.[^1] Hence

```{math}
:enumerated: false

\rank(\boldP - \boldP_1)
        = \trace(\boldP - \boldP_1)
        = \trace \boldP - \trace \boldP_1
        = K - K_1 = K_2
```

Via [fact %s](#fa-mcfi), we conclude that $\lL(Q_1) =\chi^2(K_2)$, as was to be shown.

To complete the proof, it remains to show that, under the null hypothesis and taking $\boldX$ as given, $Q_1$ and $Q_2$ are independent. To see this, observe that $Q_1$ is a function of $(\boldP - \boldP_1)\boldu$, while $Q_2$ is a function of $\boldM \boldu$. Hence it suffices to show that $(\boldP - \boldP_1)\boldu$ and $\boldM \boldu$ are independent given $\boldX$. By [fact %s](#fa-ucimin), it suffices to show that their covariance is zero. To verify this, observe that

```{math}
:enumerated: false

\cov [ (\boldP - \boldP_1) \boldu, \boldM \boldu \given \boldX ]
        = \EE [ (\boldP - \boldP_1) \boldu (\boldM \boldu)^\T \given \boldX ]
        = \EE [ (\boldP - \boldP_1) \boldu \boldu^\T \boldM  \given \boldX ]
```

Since $\boldP$, $\boldP_1$, and $\boldM$ are just functions of $\boldX$, this becomes

```{math}
:enumerated: false

(\boldP - \boldP_1) \EE [\boldu \boldu^\T\given \boldX ] \boldM  
        = \sigma^2 (\boldP - \boldP_1) \boldM  
        = \sigma^2 (\boldP - \boldP_1) (\boldI - \boldP)
```

Using idempotence and [fact %s](#fa-subsub), the matrix product on the right is

```{math}
:enumerated: false

(\boldP - \boldP_1) (\boldI - \boldP)
        = \boldP - \boldP^2 - \boldP_1 + \boldP_1 \boldP
        = \boldP - \boldP - \boldP_1 + \boldP_1 
        = \boldzero
```

This completes the proof of independence, and hence of [theorem %s](#t-ftest). ◻
````


The most common implementation of the $F$-test is the test that all coefficients of nonconstant regressors are zero. In this case {eq}`eq-altr` becomes

$$
\boldy = \boldone \beta_1 + \boldX_2 \boldbeta_2  + \boldu
$$ (eq-altrwi)

where $\boldbeta_2$ is the vector of coefficients corresponding to the nonconstant regressors. Since $\boldX_1 = \boldone$, we then have $\boldM_1 = \boldM_c$, where the latter is defined in {eq}`eq-pczm`, and hence $\rssr$ is the squared norm of $\boldy - \bar y \boldone$. An application with simulated data is given in [listing %s](#list-fstats). If you run the program, you will find that the F statistic calculated by the theoretical formula agrees with the F statistic produced by R's `summary` function.

```{code-block} r
:name: list-fstats
:caption: Calculating the F statistic (R)
:linenos:


set.seed(1234) 
N <- 50; K <- 3
beta <- rep(1, K)
x2 <- runif(N); x3 <- runif(N)
X <- cbind(rep(1, N), x2, x3)
u <- rnorm(N) 
y <- X %*% beta + u

betahat <- solve(t(X) %*% X) %*% t(X) %*% y 
residuals <- y - X %*% betahat
ussr <- sum(residuals^2)
rssr <- sum((y - mean(y))^2)

Fa <- (rssr - ussr) / 2
Fb <- ussr / (N - K)

cat("F =", Fa / Fb, "\n")        # Print result
print(summary(lm(y ~ x2 + x3)))  # For comparison

```

It is an exercise to show that in the case of {eq}`eq-altrwi`, the F statistic in {eq}`eq-fstatr` can be rewritten as

$$
F = \frac{R^2_c}{1 - R^2_c} \frac{N - K}{K_2}
$$ (eq-fstat2)

where $R^2_c$ is the centered R squared defined in {ref}`ss-cenob`. Can you provide some intuition as to why large $F$ is evidence against the null?

## Problems and Extensions

The standard OLS assumptions are strict and the results we have obtained are sensitive to their failure. Let's discuss some situations where the OLS assumptions fail and investigate the implications.

(ss-heter)=
### Nonspherical Errors
Heteroskedasticity occurs when the variance of the error term is not constant across observations, in contradiction to [assumption %s](#a-bols2). In particular, the errors are called **heteroskedastic** if the diagonal terms of $\EE [ \boldu \boldu^\T \given \boldX ]$ are not constant. If the off-diagonal terms are nonzero, then the errors are said to have **serial correlation**.

Failure of [assumption %s](#a-bols2) does not alone cause bias in the OLS estimator $\hboldbeta$. This is easily verified by reading the proof of the first part of [theorem %s](#t-olsiu). Indeed, the proof that $\hboldbeta$ is unbiased for $\boldbeta$ does not require [assumption %s](#a-bols2). At the same time, the claim that $\hat \sigma^2$ is unbiased for $\sigma$ in the second part of [theorem %s](#t-olsiu) does not make sense without [assumption %s](#a-bols2), and the expression $\var [\hboldbeta \given \boldX] = \sigma^2 (\boldX^\T \boldX)^{-1}$ from {eq}`eq-vhbb` is no longer valid. The Gauss–Markov theorem also breaks down. The most important practical consequence is that the results on inference in {ref}`ss-ne`, which depend on the expression for the variance given above, are no longer valid.

What happens if we replace the assumption $\EE [ \boldu \boldu^\T \given \boldX ] = \sigma^2 \boldI$ with the more general assumption

$$
\EE [ \boldu \boldu^\T \given \boldX ] = \boldOmega
$$ (eq-glsa)

for some positive definite matrix $\boldOmega$? ({prf:ref}`a-fr2`–{prf:ref}`a-bols` are maintained.) Using {eq}`eq-hg`, we can write the variance of the OLS estimator $\hboldbeta$ under these conditions as

$$
\var [\hboldbeta \given \boldX]
    = \EE [ (\hboldbeta - \boldbeta)(\hboldbeta - \boldbeta)^\T \given \boldX]
    = (\boldX^\T \boldX)^{-1} \boldX^\T \boldOmega  \boldX
        (\boldX^\T \boldX)^{-1}
$$ (eq-genvar)

If $\boldOmega$ is known, then we can estimate under [assumptions %s](#a-fr2)–[%s](#a-bols) and {eq}`eq-glsa` using **generalized least squares**. The idea is as follows: In view of [fact %s](#fa-cholesky), there exists a nonsingular matrix $\boldC$ such that $\boldC^\T \boldC = \boldOmega^{-1}$. Using $\boldC$, we transform the regression model by multiplying both sides of $\boldy = \boldX \boldbeta + \boldu$ by $\boldC$. This gives

```{math}
:enumerated: false

\boldy_c = \boldX_c \boldbeta + \boldu_c  
    \quad \text{where} \quad
    \boldy_c := \boldC \boldy, \,
    \boldX_c := \boldC \boldX
    \text{ and } \boldu_c := \boldC \boldu
```

{prf:ref}`ex-gls` asks you to show that

$$
\EE [ \boldu_c \given \boldX_c ] = \boldzero
    \quad \text{and} \quad
    \EE [ \boldu_c \boldu_c^\T \given \boldX_c ] = \boldI
$$ (eq-evug)

Hence the full set of OLS [assumptions %s](#a-fr2)–[%s](#a-bols2) is satisfied for the transformed data. If we use the transformed data to estimate $\boldbeta$ via $(\boldX_c^\T \boldX_c)^{-1} \boldX_c^\T \boldy_c$, we recover the best linear unbiased property. Expanding out the definitions of $\boldX_c$ and $\boldy_c$ and using $\boldC^\T \boldC = \boldOmega^{-1}$, we can show the new estimator to be

$$
\hboldbeta_{{\rm GLS}} 
    := (\boldX^\T \boldOmega^{-1} \boldX)^{-1} \boldX^\T \boldOmega^{-1} \boldy
$$ (eq-glsest)

This estimator is called the **generalized least squares (GLS) estimator** of $\boldbeta$.

In practice, $\boldOmega$ is rarely known. While we can estimate it, estimation is difficult unless additional structure is imposed. This is because $\boldOmega$ is $N \times N$, implying that, absent further restrictions, the number of parameters contained in $\boldOmega$ grows *faster* than the sample. Our estimation problem gets harder the more data that we have.

A common way of adding structure is to assume that $\boldOmega$ is diagonal. That is, $\boldOmega = \diag(\sigma_1^2, \ldots, \sigma_N^2)$, where each $\sigma_n$ is strictly positive and unknown. The diagonal assumption rules out correlation but maintains the possibility of heteroskedasticity. One estimator of $\boldOmega$ is

```{math}
:enumerated: false

\hat \boldOmega := \diag(\hat u^2_1, \ldots, \hat u^2_N)
    \quad \text{where} \quad
    \hboldu := \boldM \boldy
```

Replacing $\boldOmega$ with $\hat \boldOmega$ in {eq}`eq-genvar` leads to the estimate of $\var [\hboldbeta \given \boldX]$ given by

$$
\hat \boldV :=
    (\boldX^\T \boldX)^{-1} \boldX^\T \hat \boldOmega  \boldX
        (\boldX^\T \boldX)^{-1}
$$ (eq-robse)

If we use this expression in place of the usual variance estimate $\hat \sigma^2 (\boldX^\T \boldX)^{-1}$ when we form our standard errors, we get the **heteroskedasticity-consistent standard errors** proposed by {cite}`white1980heteroskedasticity`.

(ss-eb)=
### Bias
In {ref}`ss-heter` we considered what happens when [assumption %s](#a-bols2) (spherical errors) fails. What about the other OLS assumptions? {prf:ref}`a-fr2` (full rank) is fairly benign. {prf:ref}`a-lols` and {prf:ref}`a-bols` (linearity and exogeneity) are more problematic. Let's have a look at how they might fail and the consequences.

The linearity assumption $y_n = \boldx_n^\T \, \boldbeta + u_n$ is fundamental. If we replace it with a more general assumption such as $y_n = f(\boldx_n) + u_n$ for some arbitrary function $f$, then OLS theory falls apart. For example, we cannot even *ask* if the OLS estimator is unbiased because $\boldbeta$ has vanished. There is no parameter vector for $\hboldbeta$ to estimate.[^2]

The remaining assumption we need to consider is exogenous regressors. What happens if all other OLS assumptions hold but [assumption %s](#a-bols) fails?

Observe that, by {eq}`eq-hg`,

$$
\EE [ \hboldbeta \given \boldX ] - 
    \boldbeta = (\boldX^\T \boldX)^{-1} \boldX^\T \EE[ \boldu \given \boldX]
$$ (eq-bien)

Without $\EE [ \boldu \given \boldX ] = \boldzero$, we cannot assert that the right-hand side is zero. Hence $\hboldbeta$ is biased. Bias due to failure of [assumption %s](#a-bols) is called **endogeneity bias**.

As discussed in several parts of this text, striving for unbiasedness is somewhat misguided. Better estimators can usually be constructed by admitting some degree of bias in order to penalize complexity (thereby acknowledging that the empirical distribution is not the true distribution). However, there is an important difference with the bias in {eq}`eq-bien`. Good estimators reduce bias as the sample size increases, so that the object they aim to compute can be recovered asymptotically. This is not true for the OLS estimator when exogeneity fails.

To see this, recall from [theorem %s](#t-cwa) that $\hboldbeta \, \toprob  \, \EE[ \boldx \boldx^\T]^{-1} \, \EE[\boldx \, y]$ as $N \to \infty$. The right-hand side is the vector of coefficients in the best linear predictor. If $y = \boldx^\T \boldbeta + u$, then this becomes

```{math}
:enumerated: false

\hboldbeta \, 
     \toprob  \,
     \EE[ \boldx \boldx^\T]^{-1} \, \EE[\boldx \, (\boldx^\T \boldbeta + u)]
    =\boldbeta  + \EE[ \boldx \boldx^\T]^{-1} \, \EE[\boldx \,  u]
```

With exogeneity we get $\EE[\boldx \,  u] =  \EE[ \boldx \, \EE[ u \given \boldx ] ] = 0$ and the OLS estimator is consistent. Without exogeneity we cannot assert the same.

#### Examples of Endogeneity Bias

Let's look at some applications where [assumption %s](#a-bols) fails. First up, consider again the Cobb–Douglas example, which yields the regression model

```{math}
:enumerated: false

\ln y_n = \gamma + a \ln k_n + \delta \ln \ell_n + u_n
```

Here $y$ is output, $k$ is capital, $\ell$ is labor, and subscript $n$ indicates observation on the $n$th firm. The term $u_n$ is a firm specific productivity shock. A likely problem here is that the firm will choose higher levels of both capital and labor when it anticipates high productivity in the current period. This will lead to endogeneity bias.

To illustrate, suppose that $u_{n,-1}$ is the productivity shock received by firm $n$ last period and this value is observable to the firm. Suppose further that productivity follows a random walk, with $u_n = u_{n, -1} + \eta_n$, where $\eta_n$ is zero-mean white noise. As a result the firm forecasts period $n$ productivity as $\EE[u_n \given u_{n,-1}] = u_{n,-1}$. Finally, suppose that the firm increases labor input when productivity is anticipated to be high, with the specific relationship $\ell_n = a + b \EE[u_n \given u_{n,-1}]$ for $b > 0$. When all shocks are zero-mean, we then have

```{math}
:enumerated: false

\EE[ \ell_n u_n ] 
    = \EE[ (a + b u_{n,-1})(u_{n,-1} + \eta_n) ]
    = \EE[ b u_{n,-1}^2 ]
```

This term will be strictly positive whenever $u_{n,-1}$ has positive variance. Thus the conditions of [fact %s](#fa-sruc) fail, and therefore [assumption %s](#a-bols) does not hold (because [assumption %s](#a-bols) implies [fact %s](#fa-sruc)).

This source of endogeneity bias in estimating production functions has been discussed many times in the literature. The best solution is better modeling. For one illustration of careful modeling, see {cite}`olley1996dynamics`.

As a second example of endogeneity bias, suppose next that we have in hand data generated according to the simple AR(1) model

$$
x_0 = 0 
    \quad \text{and} \quad
    x_n = \beta x_{n-1} + u_n
    \quad \text{for } 
    n = 1,\ldots,N
$$ (eq-sar1)

Let $\{u_n\}_{n=1}^N$ be IID with distribution $\nN(0,\sigma^2)$. The unknown parameters are $\beta$ and $\sigma^2$. Letting $\boldy := (x_1,\ldots,x_N)$, $\boldx := (x_0,\ldots,x_{N-1})$, and $\boldu := (u_1,\ldots,u_N)$, we can write the $N$ equations in {eq}`eq-sar1` as $\boldy = \beta \boldx + \boldu$. The OLS estimate of $\beta$ is $\hat \beta := (\boldx^\T \boldx)^{-1} \boldx^\T\boldy$.

By [fact %s](#fa-sruc), if [assumption %s](#a-bols) holds, then $\EE [ u_m x_{n+1} ] = 0$ for any $m$ and $n$. In the current setup this equality fails. For example, if $n \geq m$, then we can ([ex. %s](#ex-exfy)) write $x_n$ as

$$
x_n = \sum_{j=0}^{n-1} \beta^j u_{n-j}
$$ (eq-exfy)

Hence $\EE[ x_n u_m] = \sum_{j=0}^{n-1} \beta^j \EE[ u_{n-j} u_m] = \beta^{n - m} \sigma^2$. In particular, [assumption %s](#a-bols) fails whenever $\beta \neq 0$.

To illustrate the bias in $\hat \beta$, [listing %s](#list-ghb) generates the data 10,000 times, computes $\hat \beta$ on each occasion and histograms the result. The parameters are $N = 25$ and $\beta = 0.9$. The output is [figure %s](#f-hghb). If you modify the code to print the sample mean, it gives around 0.82. Since the number of replications is large, this is close to $\EE[\hat \beta]$. An asymptotic 95% confidence interval for our estimate is around (0.818, 0.824).

```{code-block} python
:caption: Generates and histograms observations of $\hat \beta$ (Python)
:name: list-ghb
import numpy as np
import matplotlib.pyplot as plt

N = 25
x = np.zeros(N)
beta = 0.9                    
num_reps = 10000               
betahat_obs = np.empty(num_reps)
                        
for j in range(num_reps):
    u = np.random.randn(N)
    for t in range(N-1):
        x[t+1] = beta * x[t] + u[t+1]
    y = x[1:]       # x_1 ,..., x_N          
    x_vec = x[:-1]  # x_0, ..., x_{N-1}
    betahat_obs[j]  = np.sum(x_vec * y) / np.sum(x_vec**2)

plt.hist(betahat_obs, bins=50, alpha=0.6, density=True)
plt.show()
```

```{figure} figures/hghb.svg
:name: f-hghb

Observations of $\hat \beta$ when $\beta=0.9$
```

(ss-iv)=
### Instrumental Variables
The discussion in {ref}`ss-eb` shows that the OLS estimator is biased and inconsistent when regressors are endogenous. Here we show that a consistent estimator of $\boldbeta$ can be found when extra information is available in the form of "exogenous\" variables called **instruments**. These instruments will be collected into a $N \times J$ matrix $\boldZ$, where each column is observations on one exogenous variable. We connect $\boldZ$ and the other variables via the following assumptions:

(i) $\boldy = \boldX \boldbeta + \boldu$.

(ii) $\EE[ \boldu \given \boldZ] = \boldzero$.

(iii) $\EE[ \boldu \boldu^\T \given \boldZ] = \sigma^2 \boldI$ for some positive constant $\sigma$.

(iv) $\boldZ^\T \boldX$ has full column rank.

As [exercise %s](#ex-mbid) asks you to show, assumption (iv) implies that $J \geq K$, which is to say that we have at least as many instruments as regressors. If $J > K$, the model is said to be **overidentified**. If $J=K$, the model is called **exactly identified**. (When $J < K$, the model is called **underidentified**.)

The idea of instrumental variables least squares is to multiply $\boldy = \boldX \boldbeta + \boldu$ by $\boldZ^\T$ to produce

$$
\boldZ^\T \boldy = \boldZ^\T \boldX \boldbeta + \boldw
        \quad \text{where} \quad
        \boldw := \boldZ^\T \boldu
$$ (eq-olami)

Applying (ii) and (iii) from the IV assumptions gives $\EE[ \boldw \given \boldZ] = \boldzero$ and $\EE [ \boldw \boldw^\T \given \boldZ] = \sigma^2 \boldZ^\T \boldZ$. This suggests that we might apply generalized least squares estimation to the linear equation in {eq}`eq-olami`. We can use the expression for the GLS estimator from {eq}`eq-glsest` after substituting $\boldX$ for $\boldZ^\T \boldX$, $\boldy$ for $\boldZ^\T \boldy$, and $\boldOmega$ for $\sigma^2 \boldZ^\T \boldZ$. This leads to the **instrumental variable least squares (IVLS) estimator**

$$
\hboldbeta_{{\rm IVLS}} 
    := (\boldX^\T \boldZ (\boldZ^\T \boldZ)^{-1} \boldZ^\T \boldX)^{-1}
    \boldX^\T \boldZ (\boldZ^\T \boldZ)^{-1} \boldZ^\T \boldy
$$ (eq-ivls)

In the case where $J = K$, this reduces to the simpler expression

$$
\hboldbeta_{{\rm IVLS}} 
    := (\boldZ^\T \boldX)^{-1} \boldZ^\T \boldy
$$ (eq-ivls2)

In the preceding discussion, the analogy with GLS was not exact. (Direct application of GLS to {eq}`eq-olami` would require that $\EE[ \boldw \given \boldZ^\T \boldX] = \boldzero$ rather than $\EE[ \boldw \given \boldZ] = \boldzero$.) As a consequence we cannot claim unbiasedness and other properties by appealing to the GLS theory. Indeed $\hboldbeta_{{\rm IVLS}}$ is, in general, biased (see [ex. %s](#ex-ivib)).

At the same time, the bias is typically smaller than OLS and vanishes asymptotically (see, e.g., §8.3 of {cite:author}`davidson2004econometric` {cite:year}`davidson2004econometric`). {numref}`f-iv_example` gives an illustration. This figure corresponds to OLS and IVLS estimates of $\beta$ in $y = x \beta + u$ when $\beta=10$, $z$ is drawn independently of $u$, and $x$ is set equal to $\alpha z + (1 - \alpha) u$. Thus $z$ is exogenous and $x$ is a mixture of $z$ and $u$. The value $\alpha$ is a constant. In the simulation, $\alpha$ is set to $0.8$, so $z$ is a relatively strong instrument. The histograms are observations of the two estimators over 25,000 repetitions, each one having a sample size of $250$. IV provides a much better estimator of $\beta$ in this case.

```{figure} figures/iv_example.svg
:name: f-iv_example

Simulated draws of OLS and IVLS estimators of $\beta$
```

There is another way to derive the expression given for $\hboldbeta_{{\rm IVLS}}$ on the right-hand side of {eq}`eq-ivls`. It involves a method called **two-stage least squares**. The key idea is to project the endogenous regressors onto the column space of the exogenous instruments before running the regression. {prf:ref}`ex-2slq` walks you through the details.

### Causality

Causality is a slippery concept. Let's touch on some basic points, particularly concerning the connection between causality and regression. (References to more detailed treatments are given in {ref}`ss-frc`.)

You might have heard $R^2$ interpreted as measuring the "explanatory power" of the regressors in a particular econometric problem. The idea is that regression amounts to decomposing $\boldy$ into the sum of two orthogonal parts, the fitted values $\boldP \boldy$ and the residuals $\boldM \boldy$. As in {eq}`eq-bdtss`, the squared norm $\| \boldy\|^2 =: \tss$ can then be decomposed as $\| \boldP \boldy\|^2 =: \ess$ plus $\| \boldM \boldy\|^2 =: \rss$. This is sometimes paraphrased as "the total variation in $\boldy$ is the sum of explained and unexplained variation." The value of $R^2$ is then claimed to be the fraction of the variation in $\boldy$ "explained" by the regressors.

This terminology is problematic because the notion that the regressors "explain" variation in $\boldy$ suggests causation, and $R^2$ says nothing about causation *per se*. Instead, $R^2$ is better thought of as a measure of correlation (see {ref}`ss-cenob`). As has been observed on many occasions, correlation and causation are not the same thing. Two informal examples are as follows:

- We often see car crashes and ambulances together (correlation). This does not imply that ambulances cause crashes.

- Say we observe that people sleeping with their shoes on often wake up with sore heads. One possible explanation is that wearing shoes to bed causes headaches. Another is that both phenomena are caused by too many pints at the pub the night before.

Identifying causality in observational data can be very challenging. Here's another example, which, although it falls in the realm of traffic safety, is nonetheless interesting because it involves choices and human behavior: It has been observed that motorcycles fitted with ABS brakes have much lower crash rates than those without anti-lock brakes (see, e.g., {cite:author}`teoh2011effectiveness` {cite:year}`teoh2011effectiveness`). Does that mean fitting ABS to a given motorcycle will reduce the probability that the bike is involved in an accident by a similar amount? If so, should traffic authorities mandate the use of ABS brakes?

There's an obvious problem in jumping to this conclusion: selection bias in the sample. Cautious motorcyclists are more likely to choose bikes with safety features such as ABS. Lower crash rates for ABS bikes might therefore be caused primarily by more careful riding.

We can eliminate the selection bias by randomly assigning bikes of both types to a large group of riders and observing the results. One might expect that such randomization would reduce the observed effect, although perhaps not eliminate it entirely. In fact such a study has already been conducted for the case of passenger vehicles. In Munich, half of a large fleet of taxis were fitted with anti-lock brakes and half were fitted with conventional brakes ({cite:author}`aschenbrenner1994improved` {cite:year}`aschenbrenner1994improved`). The study lasted three years and drivers were randomly assigned across vehicles. A total of 747 accidents were observed. The authors found that accidents rates where slightly *higher* for taxis fitted with ABS.

The conclusion of the authors is that, while anti-lock brakes do improve braking performance, drivers compensate by driving faster, making shaper turns, poorer lane changes, and so forth. (This altered behavior is called "risk compensation.\") The net effect is roughly zero: mandatory ABS brakes will cause no appreciable drop in accident rates.

In general, understanding causality requires either good experiments or good theory. There is no magic econometric procedure that can reliably extract causality from observational data without a good model.

(ss-frc)=
## Further Reading
General econometrics texts treating OLS and its extensions include {cite}`hayashi2000econometrics`, {cite}`kennedy2008guide`, {cite}`greeneh`, {cite}`davidson2004econometric`, {cite}`wooldridge2010econometric`, and {cite}`cameron2005microeconometrics`. A comprehensive treatment of theory can be found in {cite}`ruud2000introduction`.

Regarding robust standard errors, asymptotics of the estimator in {eq}`eq-robse` are provided in {cite}`white1980heteroskedasticity`. Closely related is the covariance matrix estimator proposed in {cite}`newey1987simple`. See {cite}`hayashi2000econometrics`, {cite}`hansen2010econometrics`, and {cite}`cameron2005microeconometrics` for textbook treatments.

Measurement errors form another source of bias in OLS estimators. There is a standard approach to so-called errors-in-variables models based around instrumental variables. See, for example, chapter 8 of {cite}`davidson2004econometric`. For a treatment of nonlinear errors-in-variables, see {cite}`chen2011nonlinear`.

For more discussion of causality, see {cite}`freedman2009statistical`, {cite}`angrist2008mostly`, {cite}`pearl2009causality`, or {cite}`imbens2015causal`. An important notion of causality in time series models is the ability of one time series to predict another. A seminal contribution is {cite}`granger1969investigating`. See {cite}`dufour2010short` for a more recent discussion. {cite}`chernozhukov2013inference` consider how to perform inference on counterfactual distributions. {cite}`belloni2014high` discuss treatment effects in high dimensional settings.

## Exercises

```{exercise}
:label: ex-ambex

Argue that the sample mean of a random sample $y_1,\ldots,y_N$ from a given distribution $F$ can be viewed as a least squares estimator of the mean of $F$.
```


```{exercise}
:label: ex-exfoss

Prove the claims in [fact %s](#fa-exfoss).
```


```{exercise}
:label: ex-psruc

Prove the claims in [fact %s](#fa-sruc).
```


```{exercise}
:label: ex-iuc

Confirm part (i) of [fact %s](#fa-vciu).
```


```{exercise}
:label: ex-olsev

Prove [theorem %s](#t-olsev).
```


```{exercise}
:label: ex-ioes

Prove [fact %s](#fa-ioes).[^3]
```


```{exercise}
:label: ex-covpmz

Show by direct calculation that $\cov[\boldP \boldy, \boldM \boldy \given \boldX ]= \boldzero$.
```


```{exercise}
:label: ex-covbu

Show that under [assumptions %s](#a-fr2)–[%s](#a-bols2), we have $\cov [\hboldbeta, \boldM \boldu \given \boldX ] = \boldzero$.
```


```{exercise}
:label: ex-stato

Show that $\EE[\boldP \boldy \given \boldX ] = \boldX \boldbeta$ and $\var[\boldP \boldy \given \boldX ] = \sigma^2 \boldP$.
```


```{exercise}
:label: ex-statt

Show that $\EE[\boldM \boldy \given \boldX ] = \boldzero$ and $\var[\boldM \boldy \given \boldX ] = \sigma^2 \boldM$.
```


```{exercise}
:label: ex-voir

Confirm {eq}`eq-voir`.
```


```{exercise}
:label: ex-voir2

Using {eq}`eq-voir`, show that for the simple OLS model $\boldy = \beta_1 \boldone + \beta_2 \boldx + \boldu$, the variance of $\hat \beta_2$ given $\boldx$ is $\sigma^2 / \sum_{n=1}^N (x_n - \bar x)^2$.
```


```{exercise}
:label: ex-fstat2

Show that in the case of {eq}`eq-altrwi`, the F statistic in {eq}`eq-fstatr` can be rewritten as {eq}`eq-fstat2`.
```


```{exercise}
:label: ex-pymyun

Suppose that [assumption %s](#a-bolsg) holds, so that $\boldX$ and $\boldu$ are independent and $\lL(\boldu) = \nN(\boldzero, \sigma^2 \boldI)$. Show that, conditional on $\boldX$,

(i) $\boldP \boldy$ and $\boldM \boldy$ are normally distributed, and

(ii) $\boldP \boldy$ and $\boldM \boldy$ are independent.
```


```{exercise}
:label: ex-exfy

Verify the expression for $x_n$ in {eq}`eq-exfy`.
```


```{exercise}
:label: ex-ch12_ols-auto-1

Obtain [fact %s](#fa-ctdr) as a special case of [theorem %s](#t-ttest).[^4]
```


```{exercise}
:label: ex-gls

Prove the claims $\EE [ \boldu_c \given \boldX_c ] = \boldzero$ and $\EE [ \boldu_c \boldu_c^\T \given \boldX_c ] = \boldI$ from {eq}`eq-evug`.
```


```{exercise}
:label: ex-mbid

In connection with the theory of instrumental variables in {ref}`ss-iv`, show that if $\boldZ$ is $N \times J$, $\boldX$ is $N \times K$ and $\boldZ^\T \boldX$ has full column rank, then $J \geq K$.
```


```{exercise}
:label: ex-tdgp

Let all OLS assumptions hold except {prf:ref}`a-lols`. Suppose that the data-generating process is

$$
\boldy = \boldX \boldbeta + \boldZ \boldtheta + \boldu
$$ (eq-tdgp)
with $\boldtheta \neq \boldzero$ and $\EE[\boldu \given \boldX, \boldZ] = \boldzero$. Suppose further that we mistakenly ignore $\boldZ$, and simply regress $\boldy$ on $\boldX$. Show that the OLS estimator is biased whenever $\boldX^\T \boldZ = \boldzero$.
```


```{exercise}
:label: ex-ivib

Show that the assumptions of {ref}`ss-iv` are not sufficient for unbiasedness of the IVLS estimator of $\boldbeta$.
```


```{exercise}
:label: ex-2slq

This exercise asks you to construct $\hboldbeta_{{\rm IVLS}}$ in {eq}`eq-ivls` via an alternative approach called two-stage least squares. Your task is to show that the estimator produced by the following two steps is identical to the expression given in the right-hand side of {eq}`eq-ivls`.

(i) Construct $\hat \boldX$ by projecting each column of $\boldX$ onto $\colspace \boldZ$.

(ii) Regress $\boldy$ on $\hat \boldX$ using OLS.
```


### Solutions to Selected Exercises

````{solution} ex-ambex

Letting $\mu := \EE[y_n]$ we can write $y_n = \mu + u_n$ when $u_n := y_n - \mu$. In other words, $\boldy = \boldone \mu + \boldu$. The OLS estimate of $\mu$ is

```{math}
:enumerated: false

\hat \mu 
        := (\boldone^\T \boldone)^{-1} \boldone^\T \boldy
        = \frac{1}{N} \boldone^\T \boldy
        = \frac{1}{N} \sum_{n=1}^N y_n = \bar y_N
```

Reading right to left, the sample mean of $\boldy$ is the OLS estimate of the mean.
````


````{solution} ex-exfoss

From [fact %s](#fa-anox) we have $\boldM \boldX = \boldzero$ and $\boldP \boldX = \boldX$. From this (i) and (ii) follow easily. For example,

```{math}
:enumerated: false

\boldM \boldy = \boldM (\boldX \boldbeta + \boldu) = \boldM \boldu
```

Regarding (iii), since $\boldM$ is symmetric and idempotent , we have

```{math}
:enumerated: false

\rss 
        := \| \boldM \boldy \|^2
        = \| \boldM \boldu \|^2
        = (\boldM \boldu)^\T (\boldM \boldu) 
        = \boldu^\T \boldM^\T \boldM \boldu 
        = \boldu^\T \boldM\boldu 
        
```
````


```{solution} ex-psruc

The respective proofs for claims 1–4 are as follows:

(i) $\EE [\boldu] = \EE [ \EE[ \boldu \given \boldX ] ] = \EE [ \boldzero ] = \boldzero$

(ii) $\EE [ u_m \given x_{nk} ] = \EE [ \EE [u_m \given \boldX] \given x_{nk} ] = \EE [ 0 \given x_{nk} ] = 0$

(iii) $\EE [ u_m x_{nk} ] = \EE [ \EE [ u_m x_{nk} \given x_{nk} ] ] = \EE [ x_{nk} \EE [ u_m  \given x_{nk} ] ] = 0$

(iv) $\cov [ u_m, x_{nk} ] = \EE [ u_m x_{nk} ] - \EE[u_m] \EE[x_{nk}] = 0$
```


````{solution} ex-iuc

By definition, $\var \boldu = \EE[ \boldu \boldu^\T ]  - \EE[ \boldu ] \EE[ \boldu^\T ]$. Since $\EE[ \boldu ] = \EE[ \EE [ \boldu \given \boldX ]] = \boldzero$, this reduces to $\var \boldu = \EE[ \boldu \boldu^\T ]$. Moreover

```{math}
:enumerated: false

\EE[ \boldu \boldu^\T ] 
        = \EE[ \EE [ \boldu \boldu^\T \given \boldX ]]
        = \EE[ \sigma^2 \boldI ]
        = \sigma^2 \boldI 
        
```
````


````{solution} ex-olsev

If $\boldA := (\boldX^\T \boldX)^{-1} \boldX^\T$, then by {eq}`eq-hg` we have $\hboldbeta = \boldbeta + \boldA \boldu$, and

```{math}
:enumerated: false

\var [ \hboldbeta \given \boldX ]
        = \var [ \boldbeta + \boldA \boldu \given \boldX ]
        = \var [ \boldA \boldu \given \boldX ]
```

Treating $\boldA$ as nonrandom given $\boldX$ and applying [fact %s](#fa-ltavc), we have

```{math}
:enumerated: false

\var [ \boldA \boldu \given \boldX ]
        = \boldA \var [ \boldu \given \boldX ] \boldA^\T
        = \boldA ( \sigma^2 \boldI ) \boldA^\T
```

Moreover, $\boldA ( \sigma^2 \boldI ) \boldA^\T = \sigma^2 \boldA \boldA^\T = \sigma^2 (\boldX^\T \boldX)^{-1} \boldX^\T \boldX (\boldX^\T \boldX)^{-1} = \sigma^2 (\boldX^\T \boldX)^{-1}$. Hence

```{math}
:enumerated: false

\var [ \hboldbeta \given \boldX ]
        = \var [ \boldA \boldu \given \boldX ]
        = \sigma^2 (\boldX^\T \boldX)^{-1} 
        
```
````


````{solution} ex-ioes

Note first that $\boldP \boldu$ and $\boldM \boldu$ are independent by [fact %s](#fa-poci). By [fact %s](#fa-rviifi), functions of independent random vectors are likewise independent. Hence it suffices to show that $\hboldbeta$ and $\hat \sigma^2$ are functions of $\boldP \boldu$ and $\boldM \boldu$ respectively. The case of $\hat \sigma^2$ is straightforward. For $\hboldbeta$ we have

```{math}
:enumerated: false

\hboldbeta = \boldbeta + (\boldX^\T \boldX)^{-1} \boldX^\T \boldP \boldu
```

as can be seen by writing out $\boldP$ in terms of $\boldX$ and canceling to get {eq}`eq-hg`.
````


````{solution} ex-covpmz

Note that

```{math}
:enumerated: false

\EE[ \boldM \boldy \given \boldX] 
            = \EE[ \boldM \boldu \given \boldX] 
            = \boldM \EE[  \boldu \given \boldX] = \boldzero
```

Using these facts, we obtain

```{math}
:enumerated: false

\cov[\boldP \boldy, \boldM \boldy \given \boldX ]
        = \cov[\boldX \boldbeta + \boldP \boldu, \boldM \boldu \given \boldX ]
        = \EE[(\boldX \boldbeta + \boldP \boldu) (\boldM \boldu)^\T \given \boldX ]
```

By linearity of expectations and symmetry of $\boldM$, this becomes

```{math}
:enumerated: false

\cov[\boldP \boldy, \boldM \boldy \given \boldX ]
        = \EE[\boldX \boldbeta \boldu^\T \boldM \given \boldX ]
            + \EE[ \boldP \boldu \boldu^\T \boldM \given \boldX ]
```

Regarding the first term on the right-hand side, we have

```{math}
:enumerated: false

\EE[\boldX \boldbeta \boldu^\T \boldM \given \boldX ]
        = \boldX \boldbeta \EE[\boldu^\T \given \boldX ] \boldM 
        = \boldzero
```

Regarding the second term on the right-hand side, we have

```{math}
:enumerated: false

\EE[ \boldP \boldu \boldu^\T \boldM \given \boldX ]
        =  \boldP \EE[\boldu \boldu^\T  \given \boldX ]\boldM
        =  \boldP \sigma^2 \boldI \boldM
        =  \sigma^2 \boldP \boldM
        = \boldzero
```

```{math}
:enumerated: false

\fore
        \cov[\boldP \boldy, \boldM \boldy \given \boldX ]
        = \boldzero
        
```
````


````{solution} ex-covbu

Since $\boldM$ is a function of $\boldX$, we have $\EE [ \boldM \boldu \given \boldX ] = \boldM \EE [  \boldu \given \boldX ] = \boldzero$. As a result

```{math}
:enumerated: false

\begin{aligned}
\cov [\hboldbeta, \boldM \boldu \given \boldX ]
        & = \EE [\hboldbeta \, (\boldM \boldu)^\T \given \boldX ] \\
        & = \EE [ \boldbeta(\boldM \boldu)^\T + (\boldX^\T\boldX)^{-1} \boldX^\T \boldu 
            \, (\boldM \boldu)^\T \given \boldX ] \\
        & = \EE [ (\boldX^\T\boldX)^{-1} \boldX^\T \boldu 
            \, (\boldM \boldu)^\T \given \boldX ] \\
        & = \EE [ (\boldX^\T\boldX)^{-1} \boldX^\T \boldu \boldu^\T \boldM \given
            \boldX ] 
        = \sigma^2  (\boldX^\T\boldX)^{-1} \boldX^\T \boldM
\end{aligned}
```

The claim now follows from $\boldX^\T \boldM = (\boldM \boldX)^\T = \boldzero$.
````


````{solution} ex-stato

Regarding the claim that $\EE[\boldP \boldy \given \boldX ] = \boldX \boldbeta$, our previous results and linearity of expectations gives

```{math}
:enumerated: false

\EE[\boldP \boldy \given \boldX ] 
            = \EE[\boldX \boldbeta + \boldP \boldu \given \boldX ] 
            = \boldX \boldbeta + \boldP \EE[ \boldu \given \boldX ] 
            = \boldX \boldbeta
```

Regarding the claim that $\var[\boldP \boldy \given \boldX ] = \sigma^2 \boldP$, our rules for manipulating variances yield

```{math}
:enumerated: false

\var[\boldP \boldy \given \boldX ] 
            = \var[\boldX \boldbeta + \boldP \boldu \given \boldX ]
            = \var[\boldP \boldu \given \boldX ]
            = \boldP \var[ \boldu \given \boldX ] \boldP^\T
            = \boldP \sigma^2 \boldI \boldP^\T
```

Using symmetry and idempotence of $\boldP$, we obtain $\var[\boldP \boldy \given \boldX ] = \sigma^2 \boldP$.
````


```{solution} ex-statt

Similar to the solution of [exercise %s](#ex-stato).
```


````{solution} ex-voir

Repeating {eq}`eq-exir2` we have

$$
\hat \beta_k 
        = \beta_k + (\col_k(\boldX)^\T \boldM_1 \col_k(\boldX))^{-1}
        \col_k(\boldX)^\T \boldM_1 \boldu
$$ (eq-exir21)
Since $\beta_k$ is constant, taking the variance of {eq}`eq-exir21` conditional on $\boldX$, we obtain

```{math}
:enumerated: false

\var[ \hat \beta_k  \given \boldX ]
        = \var[ \boldA  \boldu \given \boldX ]
        \quad \text{where }
        \boldA := (\col_k(\boldX)^\T \boldM_1 \col_k(\boldX))^{-1} 
            \col_k(\boldX)^\T \boldM_1
```

Since $\boldA$ is a function of $\boldX$, we can treat it as constant given $\boldX$, and we obtain

```{math}
:enumerated: false

\var[ \hat \beta_k  \given \boldX ]
        = \boldA \var[ \boldu \given \boldX ] \boldA^\T
        = \boldA \sigma^2 \boldI \boldA^\T
        = \sigma^2 \boldA \boldA^\T
```

To complete the proof, we just observe that

```{math}
:enumerated: false

\begin{aligned}
\boldA \boldA^\T 
        & = (\col_k(\boldX)^\T \boldM_1 \col_k(\boldX))^{-1} 
            \col_k(\boldX)^\T \boldM_1
            \boldM_1 ^\T \col_k(\boldX) 
        (\col_k(\boldX)^\T \boldM_1 \col_k(\boldX))^{-1}  \\
        & = (\col_k(\boldX)^\T \boldM_1 \col_k(\boldX))^{-1}
\end{aligned}
```

where the last equality is due to symmetry and idempotence of $\boldM_1$. We conclude that {eq}`eq-voir` is valid.
````


````{solution} ex-voir2

We will apply {eq}`eq-voir`. For the model $\boldy = \beta_1 \boldone + \beta_2 \boldx + \boldu$, the matrix $\boldM_1$ is the residual projection $\boldM_c$ associated with $\boldone$, and $\col_k \boldX$ is just $\boldx$. Hence, from {eq}`eq-voir`, we have

```{math}
:enumerated: false

\var[ \hat \beta_2 \given \boldX ]
        = \sigma^2 (\boldx^\T \boldM_c \boldx)^{-1}
```

Using symmetry and idempotence of $\boldM_c$, this becomes

```{math}
:enumerated: false

\var[ \hat \beta_2 \given \boldX ]
        = \sigma^2 [(\boldM_c \boldx)^\T \boldM_c \boldx]^{-1}
        = \sigma^2 / \sum_{n=1}^N (x_n - \bar x)^2
```

Finally, since the only random variables in $\boldX$ are the random variables in $\boldx$, we obtain the solution

```{math}
:enumerated: false

\var[ \hat \beta_2 \given \boldx ]
        = \sigma^2 / \sum_{n=1}^N (x_n - \bar x)^2
        
```
````


````{solution} ex-fstat2

We need to show that in the special case {eq}`eq-altrwi` we have

```{math}
:enumerated: false

\frac{(\rssr - \ussr) / K_2}{ \ussr / (N - K)}
        = \frac{R^2_c}{1 - R^2_c} \frac{N - K}{K_2}
```

or, equivalently,

$$
\frac{\rssr - \ussr}{ \ussr } = \frac{R^2_c}{1 - R^2_c}
$$ (eq-fstatrp)
Consider first the left-hand side of {eq}`eq-fstatrp`. In the case of {eq}`eq-altrwi`, this becomes

```{math}
:enumerated: false

\frac{\rssr - \ussr}{ \ussr } 
        = \frac{ \| \boldM_c \boldy \|^2 - \| \boldM \boldy \|^2 }{\| \boldM \boldy \|^2}
```

On the other hand, regarding the right-hand side of {eq}`eq-fstatrp`, the definition of $R^2_c$ and some minor manipulation gives

```{math}
:enumerated: false

\frac{R^2_c}{1 - R^2_c} = 
        \frac{ \| \boldP \boldM_c \boldy \|^2 }
            {\| \boldM_c \boldy \|^2 - \| \boldP \boldM_c \boldy \|^2}
```

Hence, to establish {eq}`eq-fstatrp`, we need to show that

```{math}
:enumerated: false

\frac{ \| \boldM_c \boldy \|^2 - \| \boldM \boldy \|^2 }{\| \boldM \boldy \|^2}
        = \frac{ \| \boldP \boldM_c \boldy \|^2 }
            {\| \boldM_c \boldy \|^2 - \| \boldP \boldM_c \boldy \|^2}
```

This is can be established using {eq}`eq-rsamy`.
````


```{solution} ex-pymyun

We have shown previously that $\boldM \boldy = \boldM \boldu$ and $\boldP \boldy = \boldX \boldbeta + \boldP \boldu$. Since we are conditioning on $\boldX$, we can treat it as constant. When $\boldX$ is constant, $\boldP$, $\boldM$ and $\boldX \boldbeta$ are all constant. Since linear (or affine) transformations of normal random vectors are normal, $\boldM \boldu$ and $\boldX \boldbeta + \boldP \boldu$ are both normally distributed.

It remains to show that $\boldP \boldy$ and $\boldM \boldy$ are independent given $\boldX$. Since they are normally distributed given $\boldX$, we need only show that they are uncorrelated given $\boldX$. This was already proved in [exercise %s](#ex-covpmz).
```


````{solution} ex-gls

The proof of $\EE [ \boldu_c \given \boldX_c ] = \boldzero$ is easy. Regarding the second claim, recall that $\boldC$ is a nonsingular matrix chosen to satisfy $\boldC^\T \boldC = \boldOmega^{-1}$. (Existence is from the Cholesky decomposition; see page .) From $\boldu_c = \boldC \boldu$ we have

```{math}
:enumerated: false

\EE [ \boldu_c \boldu_c^\T \given \boldX_c ] 
        = \boldC\, \EE[ \boldu \boldu^\T  \given \boldX_c ] \boldC^\T 
        = \boldC\, \boldOmega \boldC^\T
```

Since $\boldC$ is nonsingular, we have

```{math}
:enumerated: false

\boldC \boldOmega \boldC^\T 
        = \boldC (\boldC^\T \boldC)^{-1} \boldC^\T 
        = \boldC \boldC^{-1} (\boldC^\T)^{-1} \boldC^\T 
        = \boldI
        
```
````


```{solution} ex-mbid

If not then $\boldZ^\T \boldX$ has fewer rows than columns, and hence its columns cannot be linearly independent. See [theorem %s](#t-exth).
```


````{solution} ex-tdgp

The OLS estimator is $\hboldbeta = (\boldX^\T \boldX)^{-1} \boldX^\T \boldy$. Substituting into {eq}`eq-tdgp` gives

```{math}
:enumerated: false

\hboldbeta 
         = (\boldX^\T \boldX)^{-1} \boldX^\T 
            (\boldX \boldbeta + \boldZ \boldtheta + \boldu)
         = \boldbeta +
            (\boldX^\T \boldX)^{-1} \boldX^\T \boldZ \boldtheta 
            + (\boldX^\T \boldX)^{-1} \boldX^\T \boldu
```

Therefore

```{math}
:enumerated: false

\EE[ \hboldbeta ]
        = \EE[ \EE[ \hboldbeta \given \boldX, \boldZ] ]
        = \boldbeta +
        \EE[ (\boldX^\T \boldX)^{-1} \boldX^\T \boldZ ] \boldtheta
```

If $\boldX^\T \boldZ = \boldzero$ then the last term on the right-hand side drops out, and $\hboldbeta$ is unbiased. Otherwise $\hboldbeta$ is biased.
````


````{solution} ex-ivib

To solve this exercise it suffices to give an example of a model that satisfies the assumptions of {ref}`ss-iv` and yet produces a biased IV estimator. To this end, suppose that we have only one observation, and let all variables be scalar. Let $z$ and $u$ be independent random variables with $\lL(u) = \nN(0, \sigma^2)$ and let $x = u$. Let $y = x\beta + u$. While the IV assumptions of {ref}`ss-iv` are satisfied, the IVLS estimator is biased. Indeed, the expression {eq}`eq-ivls2` combined with $y = x\beta + u$ gives us

```{math}
:enumerated: false

\hat \beta_{{\rm IVLS}} 
        = (z x)^{-1} z y
        = \beta + (z x)^{-1}{z u}
        = \beta + \frac{u}{x}
```

Using $x=u$ gives $\EE [ \hat \beta_{{\rm IVLS}} ]  =\hat \beta_{{\rm IVLS}} = \beta + 1$.
````


````{solution} ex-2slq

By [theorem %s](#t-projon), the matrix that implements orthogonal projection onto $\colspace \boldZ$ is $\boldP_Z := \boldZ (\boldZ^\T\boldZ)^{-1} \boldZ^\T$. Applying this operation to each column of $\boldX$ is equivalent to premultiplying $\boldX$ by $\boldP_Z$. Regressing $\boldy$ on $\boldP_Z \boldX$ using OLS produces the estimator

```{math}
:enumerated: false

\hboldbeta_{{\rm 2SLS}} 
        = [ (\boldP_Z \boldX)^\T (\boldP_Z\boldX) ]^{-1} 
            (\boldP_Z \boldX)^\T \boldy
        = [ \boldX^\T \boldP_Z\boldX ]^{-1} \boldX^\T \boldP_Z  \boldy
```

(The second equality uses symmetry and idempotence of $\boldP_Z$. See [fact %s](#fa-pmsi).) Using the definition of $\boldP_Z$ gives {eq}`eq-ivls`.
````


[^1]: Hint: See [fact %s](#fa-subsub).

[^2]: It is an interesting habit of economists to quibble endlessly about whether someone's estimators are biased, even when their models are wildly speculative in terms of specific functional forms. (For example, do we really believe that log GDP per capita is a linear function of various measures of institutional quality and an Africa dummy, say, or even that this is a close depiction of reality? If not, what exactly does it mean to say that $\hboldbeta$ is biased?)

[^3]: Hint: Use [fact %s](#fa-poci).

[^4]: Hint: Consider the regression model $\boldy = \boldone \beta + \boldu$.
