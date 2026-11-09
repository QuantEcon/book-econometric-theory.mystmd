---
title: "Regression"
label: c-reg
---

(s-olsfp)=
## Linear Regression
Linear regression is one of the core topics of statistics. It is even more central to econometrics, where a shortage of controlled experiments often leaves econometricians trying to factor heterogeneity out ex post. We begin with a nontraditional view of linear regression based on minimal assumptions.

(ss-mvls)=
### The Setup
Let's start with the kind of prediction problem discussed in {ref}`ss-erm`. We study a system with vector input $\boldx_n \in \RR^K$ followed by scalar output $y_n$. For example,

- $\boldx_n$ is a description of a lottery (probabilities, possible outcomes, etc.) in a controlled experiment and $y_n$ is willingness to pay in order to participate (see, e.g., {cite:author}`peysakhovich2015machine` {cite:year}`peysakhovich2015machine`).

- $\boldx_n$ is a set of household characteristics (ethnicity, age, location, etc.) and $y_n$ is household wealth at some later date (see, e.g., {cite:author}`mckernan2014disparities` {cite:year}`mckernan2014disparities`).

- $\boldx_n$ is price of electricity, prices of alternatives, temperature, household income, and measurements of the regional income distribution, while $y_n$ is regional electricity consumption (see, e.g., {cite:author}`auffhammer2014powering` {cite:year}`auffhammer2014powering`).

Although we don't exclude the possibility that $y$ is categorical (i.e., discrete), our loss function will be oriented toward regression (where $y$ takes values in $\RR$).

Suppose that we have $N$ observations $\boldz_n := (\boldx_n,y_n)$, all of which are draws from some fixed joint distribution $P$. Since $P$ is fixed, we are assuming that the system is stationary across the set of draws. Our aim is to predict new output values from input values on the basis of this data. In particular, our problem is to

$$
\text{choose a function } f \colon \RR^K \to \RR \text{ such that }
    f(\boldx) \text{ is a good predictor of } y
$$ (eq-prp)

To define "good predictor\" mathematically, we need a loss function. Throughout this chapter we will be using quadratic loss. Thus, in the language of {ref}`ss-erm`, our aim is to minimize the prediction risk

$$
R(f) := \EEP  (y - f(\boldx))^2
$$ (eq-rf2)

As we saw in {ref}`ss-ce`, the minimizer of {eq}`eq-rf2` over the set of all $\bB$-measurable functions is the regression function $f^*(\boldx) := \EEP  [y \given \boldx]$. If we could compute this, then all our problems would be solved. But we cannot compute it because $P$ is not known. Instead we apply the principle of empirical risk minimization (see {ref}`ss-erm`), which leads to the problem

$$
\min_{f \in \hH} \Remp(f)
    \quad \text{where} \quad
    \Remp(f) := \frac{1}{N} \sum_{n=1}^N (y_n - f(\boldx_n))^2
$$ (eq-ermm)

Here $\hH$ is the hypothesis space, a set of candidate functions mapping $\RR^K$ into $\RR$. For obvious reasons, the problem {eq}`eq-ermm` is called a **least squares** problem.

If we take $\hH$ to be the set of all functions from $\RR^K$ to $\RR$, then, provided the input vectors are all distinct, we can set the empirical risk $\Remp(f)$ to zero by choosing any function $f$ satisfying $y_n = f(\boldx_n)$ for all $n$. However, as discussed at length in {ref}`ss-ermls`, minimizing empirical risk is different from minimizing the prediction risk $R(f)$. The latter is what we actually want to minimize. Thus $\hH$ must be restricted.

In this chapter we consider the case $\hH = \llL$, where $\llL$ is all linear functions from $\RR^K$ to $\RR$. Recalling [theorem %s](#t-lmaeq), we can write

$$
\llL = 
    \left\{ \, \text{all } f \colon \RR^K \to \RR \, 
        \text{ such that }
        f(\boldx) = \boldx^\T \boldb \text{ for some } \boldb \in \RR^K 
    \right\}
$$ (eq-allaff2)

The problem {eq}`eq-ermm` then reduces to

$$
\min_{\boldb \in \RR^K} \;
        \sum_{n=1}^N (y_n - \boldx_n^\T \boldb)^2
$$ (eq-olslq)

The term $\frac{1}{N}$ has been dropped because it does not affect the minimizer. This is the multivariate version of {eq}`eq-llsp`.

The idea of choosing $\boldb$ to minimize {eq}`eq-olslq` is intuitive: We are choosing a "line of best fit\" to minimize in-sample prediction error. This optimization problem has a long tradition. It dates back at least as far as Carl Gauss's work on the orbital position of Ceres, published in 1801.

You might be wondering whether the choice $\hH = \llL$ is a suitable one. This is an excellent question. It might not be. However, there are good reasons to start with $\llL$, even in this setting where no linearity assumptions are imposed. First, $\llL$ is a natural starting point when seeking a class of simple, well-behaved functions. Second, as we'll see, setting $\hH = \llL$ allows us to obtain an analytical expression for the minimizer, which simplifies both analysis and computation. Third, the technique has an extension from $\llL$ to broader classes of functions, as described in {ref}`ss-basis`.

(ss-tlse)=
### The Least Squares Estimator
Now let's solve {eq}`eq-olslq`. With our knowledge of overdetermined systems (see {ref}`ss-odse`), we already have all the necessary tools. This will be more obvious after we switch to matrix notation. To do this, let

$$
\boldy := 
    \left(
    \begin{array}{c}
        y_1 \\
        y_2 \\
        \vdots \\
        y_N
    \end{array}
    \right),
    \quad
    \boldx_n := 
    \left(
    \begin{array}{c}
        x_{n1} \\ 
        x_{n2} \\
        \vdots \\
        x_{nK}
    \end{array}
    \right)
    = \text{ } n \text{th observation of all regressors}
$$ (eq-design0)

and

$$
\boldX := 
    \left(
    \begin{array}{c}
        \boldx_1^\T    \\
        \boldx_2^\T  \\
        \vdots     \\
        \boldx_N^\T  
    \end{array}
    \right)
    :=:
    \left(
    \begin{array}{cccc}
        x_{11} & x_{12} & \cdots & x_{1K} \\
        x_{21} & x_{22} & \cdots & x_{2K} \\
        \vdots & \vdots &  & \vdots \\
        x_{N1} & x_{N2} & \cdots & x_{NK} 
    \end{array}
    \right)
$$ (eq-design)

Sometimes $\boldX$ is called the **design matrix**. By construction, $\col_k \boldX =$ all observations on the $k$th regressor. Also, for any $\boldb \in \RR^K$, we have

```{math}
:enumerated: false

\boldX \boldb = 
    \left(
    \begin{array}{c}
        \boldx_1^\T \boldb  \\
        \boldx_2^\T \boldb\\
        \vdots \\
        \boldx_N^\T \boldb
    \end{array}
    \right)
```

It follows that the objective function in {eq}`eq-olslq` can be written as

```{math}
:enumerated: false

\sum_{n=1}^N (y_n - \boldx_n^\T \boldb)^2 
        = \| \boldy - \boldX \boldb \|^2
```

Since strictly increasing transforms preserve the set of minimizers (see {ref}`s-opt`),

$$
\argmin_{\boldb \in \RR^K} \| \boldy - \boldX \boldb \|^2
     = \argmin_{\boldb \in \RR^K} \| \boldy - \boldX \boldb \|
$$ (eq-olses)

We already know how to solve for the minimizer on the right-hand side of {eq}`eq-olses`. By [theorem %s](#t-lssol), the solution is

$$
\hboldbeta := (\boldX^\T \boldX)^{-1} \boldX^\T \boldy
$$ (eq-thelsqe)

Traditionally, this random vector $\hboldbeta$ is called the **least squares estimator**. Once we move to more classical assumptions it will be an estimator of a particular parameter vector. At this stage it just defines our answer to the problem posed in {eq}`eq-prp`. That is,

```{math}
:enumerated: false

\text{given } \boldx \in \RR^K \text{, our prediction of } y \text{ is }
    f(\boldx) = \boldx^\T \hboldbeta
```

In terms of geometric interpretation, since $\boldX \hboldbeta$ solves {eq}`eq-olses`, it is the closest point in $\colspace \boldX$ to $\boldy$. In particular,

```{math}
:enumerated: false

\boldP \boldy = \boldX \hboldbeta   
    \quad \text{when} \quad
    \boldP := \proj \, (\colspace \boldX )
```

(See {eq}`eq-ibprojon`.) In what follows, $\boldM$ is the residual projection, as defined in {eq}`eq-ann0`.

#### Assumptions

{prf:ref}`t-lssol` and our definition of $\hboldbeta$ in {eq}`eq-thelsqe` require that $\boldX$ has full column rank (or, equivalently, that the columns of $\boldX$ are linearly independent—see page ).

```{prf:assumption}
:label: a-fr

$\boldX$ has full column rank with probability one.
```


By [theorem %s](#t-exth), $N \geq K$ is a necessary condition for [assumption %s](#a-fr) to hold. (If $N < K$, then $\RR^N$, which is necessarily spanned by $N$ vectors, cannot contain $K$ linearly independent vectors.)

If [assumption %s](#a-fr) fails, then a minimizer of {eq}`eq-olses` still exists but is no longer unique (see [ex. %s](#ex-lsiosol)). While we can treat this case, it rarely occurs in well-designed regression problems.

Let's also put some mild regularity conditions on the common joint distribution $P$ of each data point $\boldz_n := (\boldx_n, y_n)$.

```{prf:assumption}
:label: a-poprank

$P$ is such that all elements of $\EE_P [ \boldz_n \boldz_n^\T]$ are finite. Moreover

$$
\Sigmax := \EE_P [ \boldx_n \boldx_n^\T ] \; \text{ is finite and positive definite}
$$ (eq-smx)
```


Finite second moments are imposed because we want to evaluate expected squared errors. This assumption cannot be weakened unless we are willing to work with a different loss function. Positive definiteness of $\Sigmax$ ensures that the asymptotic limit of our estimator is well defined.[^1]

(ss-stnot)=
#### Notation
There's a range of standard notation associated with linear least squares estimation. Let's collect it in one place. First, the projection

```{math}
:enumerated: false

\hboldy := \boldX \hboldbeta = \boldP \boldy
```

is called the **vector of fitted values**. The $n$th fitted value $\hat y_n$ is the prediction $\boldx_n^\T \, \hboldbeta$ associated with least squares estimate and the $n$th observation $\boldx_n$ of the input vector. The vector $\boldM \boldy$ is often denoted $\hboldu$, and called the **vector of residuals**:

```{math}
:enumerated: false

\hboldu := \boldM \boldy = \boldy - \hboldy
```

The vector of residuals corresponds to the error that occurs when $\boldy$ is approximated by $\boldP \boldy$. From [fact %s](#fa-opt3) we have

$$
\boldM \boldy \perp \boldP \boldy 
    \quad \text{and} \quad 
    \boldy =  \boldP \boldy + \boldM \boldy
$$ (eq-orr)

In other words, $\boldy$ can be decomposed into two orthogonal vectors $\boldP \boldy$ and $\boldM \boldy$, where the first represents the best approximation to $\boldy$ in $\colspace \boldX$, and the second represents the residual.

Related to the fitted values and residuals, we have some standard definitions:

- **Total sum of squares** $:= \tss := \| \boldy \|^2$.

- **Residual sum of squares** $:= \rss := \| \boldM \boldy \|^2$.

- **Explained sum of squares** $:= \ess := \| \boldP \boldy \|^2$.

By {eq}`eq-orr` and the Pythagorean law,

$$
\tss = \ess + \rss
$$ (eq-bdtss)

When running regressions it is conventional to report the **coefficient of determination**, or $R^2$. The plain vanilla definition of $R^2$ is

$$
R^2 := \frac{\ess}{\tss}
$$ (eq-rsd)

Many regression packages report an alternative definition of $R^2$. See {ref}`ss-cenob` below.

### Out-of-Sample Fit

We have stressed a number of times that learning from data (statistics) means generalization from current observations to new ones. As such, the most important measure of success for a statistical procedure is out-of-sample fit. So how does linear least squares perform out-of-sample? We start with a general observation about linear predictors.

````{prf:theorem}
:label: t-dclee

If $\ell$ is the linear function $\ell (\boldx) = \boldx^\T \boldb$, then

```{math}
:enumerated: false

R(\ell)
        = \EE (y - f^*(\boldx))^2  
        + \EE (f^*(\boldx) - \boldx^\T \boldb^*)^2 
        +  (\boldb^* - \boldb)^\T 
            \Sigmax (\boldb^* - \boldb)
```
````


Here $f^*$ is the regression function and $\boldb^* = \Sigmax^{-1} \, \EE[\boldx \, y]$ is the vector of coefficients in the best linear predictor . $R(f)$ is the prediction risk of $f$ and expectations are taken under the unknown joint distribution $P$ of the pairs $(\boldx, y)$.

The proof of [theorem %s](#t-dclee) is given in {ref}`ss-rpro`. For now let's look at interpretation. The general question is how well we can generalize (i.e., reduce prediction risk) using linear functions. {prf:ref}`t-dclee` decomposes the prediction risk of an arbitrary linear predictor $\ell (\boldx) = \boldx^\T \boldb$ into three terms:

(i) The **intrinsic risk** $\EE (y - f^*(\boldx))^2$.

(ii) The **approximation error** $\EE (f^*(\boldx) - \boldx^\T \boldb^*)^2$.

(iii) The **estimation error** $(\boldb^* - \boldb)^\T \Sigmax (\boldb^* - \boldb)$.

The intrinsic risk is also called Bayes risk (see [example %s](#eg-ebrf2)). It is the residual error after $y$ is approximated with the best possible predictor (i.e., the regression function). It is large to the extent that $y$ is hard to predict using $\boldx$.

The approximation error or *bias* is the deviation between the best predictor and the best linear predictor. It reflects the cost of our decision to approximate the regression function using a linear architecture. If this architecture is held fixed, the approximation error is also fixed and cannot be reduced in the estimation process.

The estimation error is caused by the deviation of our estimator from the best linear predictor $\boldb^*$. This deviation occurs because we are predicting using finite sample information on the joint distribution of $(\boldx, y)$.

{prf:ref}`t-dclee` tells us that once $\hH$ is set to the class of linear functions, the best we can do is find an estimation method (a learning algorithm) that produces an estimate that is close to $\boldb^*$ on average when the sample size is sufficiently large. The next result states that the least squares estimator $\hboldbeta$ has this property.

```{prf:theorem}
:label: t-cwa

Let [assumptions %s](#a-poprank)–[%s](#a-fr) hold and let $\hboldbeta_N$ be the least squares estimator given sample size $N$. If the observations $\{\boldz_n\}$ are independent, then

$$
\hboldbeta_N \; \toprob \; \boldb^*
          \quad \text{as} \quad
        N \to \infty
$$ (eq-cwa)
```


The proof is below. Independence is required only for the LLN to function. We could weaken this to ergodicity (see {ref}`ss-spdef`) and obtain the same conclusion.[^2]

On one level, [theorems %s](#t-dclee) and [%s](#t-cwa) are reassuring. They tell us that if the underlying process is relatively linear, then we will attain small risk asymptotically. On the other hand, the Glivenko–Cantelli theorem tells us that we can learn everything about the underlying distribution in the limit. Here we are bounded away from this kind of consistency whenever the approximation error is positive.

As will be discussed in [chapter 14](#c-rms), a general principle of induction is that we should introduce bias in finite samples to avoid overfitting, while at the same time reducing bias asymptotically, as the empirical distribution converges to the true distribution. By comparison, in standard linear regression the bias is held fixed by the linearity assumption.

(ss-rpro)=
#### Proofs
````{prf:proof}
:nonumber:

**Proof of [theorem %s](#t-dclee).**

Fix $\boldb \in \RR^K$ and let $\ell (\boldx) = \boldx^\T \boldb$. In view of {eq}`eq-ebrf2`, we have the prediction risk

```{math}
:enumerated: false

R(\ell)
        = \EE[ (y - f^*(\boldx))^2 ] 
        + \EE[ (f^*(\boldx) - \boldx^\T \boldb)^2 ]
```

Hence the result will be established if we can show that

$$
\EE[ (f^*(\boldx) - \boldx^\T \boldb)^2 ]
        = \EE[ (f^*(\boldx) - \boldx^\T \boldb^*)^2 ]
        + \EE[ (\boldb^* - \boldb)^\T 
            \boldx \boldx^\T (\boldb^* - \boldb) ]
$$ (eq-odbb)
To see that {eq}`eq-odbb` holds, observe that

$$
f^*(\boldx) - \boldx^\T \boldb
        = f^*(\boldx) - \boldx^\T \boldb^* 
        + \boldx^\T ( \boldb^* - \boldb)
$$ (eq-odbb2)
The terms $f^*(\boldx) - \boldx^\T \boldb^*$ and $\boldx^\T ( \boldb^* - \boldb)$ are orthogonal. The reason is that $\boldx^\T \boldb^*$ is the orthogonal projection of $f^*(\boldx)$ onto $S = \Span\{\boldx\}$, the linear subspace of $L_2$ spanned by all linear combinations of the form $\bolda^\T \boldx$. (See [ex. %s](#ex-rtblp).) As such, $f^*(\boldx) - \boldx^\T \boldb^*$ is orthogonal to every element of the target subspace $\Span\{\boldx\}$. This includes $\boldx^\T ( \boldb^* - \boldb)$.

For any orthogonal elements $u$ and $v$ of $L_2$ we have $\EE[(u + v)^2] = \EE[u^2] + \EE[v^2]$. (This is the Pythagorean law in $L_2$.) Squaring both sides of {eq}`eq-odbb2`, taking expectations and applying this law gives {eq}`eq-odbb`. The proof of [theorem %s](#t-dclee) is done. ◻
````


````{prf:proof}
:nonumber:

**Proof of [theorem %s](#t-cwa).**

The proof of [theorem %s](#t-cwa) is not hard if we express $\hboldbeta_N$ in a slightly different way. Multiplying and dividing by $N$ in the definition of $\hboldbeta_N$ and then expanding out the matrix products (see [ex. %s](#ex-veseeo)) gives

$$
\hboldbeta_N 
        = \left[ \frac{1}{N} \, \boldX^\T\boldX \right]^{-1} 
                \cdot \; \frac{1}{N} \, \boldX^\T \boldy
        = \left[ \frac{1}{N} \sum_{n=1}^N \boldx_n \boldx_n^\T \right]^{-1} 
            \cdot \; \frac{1}{N} \sum_{n=1}^N \boldx_n y_n
$$ (eq-seeo0)
By the matrix LLN in [fact %s](#fa-llnmat), we have

```{math}
:enumerated: false

\frac{1}{N} \sum_{n=1}^N \boldx_n \boldx_n^\T 
             \toprob \Sigmax
         \quad \text{and} \quad
        \frac{1}{N} \sum_{n=1}^N \boldx_n y_n
             \toprob \EE[\boldx y]
         \quad \text{as} \quad
         N \to \infty
```

By [fact %s](#fa-cmtetcv1), convergence in probability is preserved over the taking of inverses and products. Hence $\hboldbeta_N \, \toprob \, \Sigmax^{-1} \EE[\boldx \, y] = \boldb^*$, as was to be shown. ◻
````


(s-gof)=
### In-Sample Fit
In-sample fit measures how well a given model fits the same data set that it was estimated on. The difference between in-sample fit (empirical risk) and out-of-sample fit (risk) was discussed in {ref}`ss-ermls`. In-sample fit of a regression is often measured with $R^2$ (see {eq}`eq-rsd`). Let's make some further comments on $R^2$ and then discuss how $R^2$ relates to in-sample fit.

```{prf:property}

$0 \leq R^2 \leq 1$ with $R^2 = 1$ if and only if $\boldy \in \colspace \boldX$.
```


That $R^2 \leq 1$ is immediate from $\| \boldP \boldy \| \leq \| \boldy \|$ (cf. [theorem %s](#t-opt2)). {prf:ref}`ex-yibx` asks you to prove the second claim. More generally, a high $R^2$ indicates $\boldy$ is relatively close to $\colspace \boldX$. This fact suggests that we can increase $R^2$ at least weakly by adding regressors. As we do so the column space of $\boldX$ expands, pushing out towards $\boldy$. Here's a formal statement:

````{prf:property}
:label: fa-rsqe

Let $\boldX_a$ and $\boldX_b$ be two design matrices. If $R^2_a$ and $R^2_b$ are the respective coefficients of determination, then

```{math}
:enumerated: false

\colspace \boldX_a \subset \colspace \boldX_b 
         \; \implies \;
         R^2_a \leq R^2_b
```
````


For a proof, see [exercise %s](#ex-rsqe) and its solution.

High $R^2$ is sometimes equated with successful regression. This is a misunderstanding of the aim of statistics. The correct definition of statistical learning is effective generalization from existing data. In the present context this means that the linear predictor produced by regression attains low risk.

So how does $R^2$ relate to risk? What $R^2$ actually measures is the degree to which empirical risk is minimized. To see this, note that

```{math}
:enumerated: false

R^2 
    = 1 - \frac{\rss}{ \tss }
    = 1 - N \frac{\Remp(\hat f) }{\tss}
```

where $\Remp$ is as defined in {eq}`eq-ermm` and $\hat f$ is our linear predictor $\hat f(\boldx) = \boldx^\T \hboldbeta$. Thus high $R^2$ means low empirical risk and good in-sample fit. But low empirical risk is no guarantee of low prediction risk, as was emphasized in {ref}`ss-ermls`.

Here's a simulation that shows how we can produce high $R^2$ without estimating anything meaningful. We take $x_n$ and $y_n$ as independent draws from a uniform distribution on $[0,1]$. By construction, there is no relationship between these two variables. For the regressors we take the powers $1, x, x^2, \ldots, x^K$, where $K$ is a positive integer. The `R` code below runs these regressions for different values of $K$. At $K=25$ the value of $R^2$ is around 0.95. This is despite the fact that no relationship exists between $x$ and $y$.

``` r

set.seed(1234)
N <- 25
y <- runif(N)
x <- runif(N)
X <- rep(1, N)

Kmax <- 25
for (K in 1:Kmax) {
    X <- cbind(X, x^K)
    results <- lm(y ~ 0 + X)
    Py2 <- sum(results$fitted.values^2)
    y2 <- sum(y^2)
    cat("K =", K, "R^2 =", Py2 / y2, "\n")
}

```

(You can obtain all code from the text at `johnstachurski.net/emet`.)

To finish this section, let's draw a connection between [fact %s](#fa-rsqe), which says that the value of $R^2$ is at least weakly increasing in the number of right-hand side variables, and [fact %s](#fa-emprd). Suppose that $\boldx$ lists a large number of possible regressors. Let the hypothesis space be

$$
\hH_j := 
    \left\{ \text{all } f \colon \RR^j \to \RR \st
            f(\boldx) = \boldx^\T \boldb \text{ for some }
        \boldb \in \RR^j
    \right\}
$$

Here $1 \leq j \leq K$. Empirical risk minimization over $\hH_j$ is equivalent to linear regression over the first $j$ regressors. Empirical risk falls as $j$ increases by [fact %s](#fa-emprd). Hence $R^2$ increases. This is the same conclusion as [fact %s](#fa-rsqe).

## The Geometry of Least Squares

In this section we cover transformations of the data and an important theorem on subsets of the least squares estimator.

(ss-basis)=
### Transformations and Basis Functions
In discussing the decision to set $\hH = \llL$ in {ref}`ss-mvls`, we mentioned that we can use many of the same ideas when extending $\hH$ to a broader class of functions. The idea is to first transform the data using some arbitrary function $\boldphi \colon \RR^K \to \RR^J$. The action of $\boldphi$ on $\boldx \in \RR^K$ is

```{math}
:enumerated: false

\boldx \mapsto \boldphi(\boldx)
    =
    \left(
    \begin{array}{c}
        \phi_1(\boldx) \\ 
        \phi_2(\boldx) \\
        \vdots \\
        \phi_J(\boldx)
    \end{array}
    \right) \in \RR^J
```

The individual functions $\phi_1, \ldots, \phi_J$ mapping $\RR^K$ into $\RR$ are sometimes called **basis functions**. In machine learning texts, the range of $\boldphi$ is called **feature space**. Linear least squares is now applied in feature space. That is, we solve the empirical risk minimization problem when the hypothesis space is

```{math}
:enumerated: false

\hH_{\boldphi} := \{ \text{all functions } \ell \circ \boldphi,  
    \text{ where } \ell \text{  is a linear function from } \RR^J \text{ to } \RR \}
```

The empirical risk minimization problem is then

$$
\min_{\ell } \;
        \sum_{n=1}^N \{ y_n - \ell(\boldphi(\boldx_n)) \}^2
        =
        \min_{\boldgamma \in \RR^J} \;
        \sum_{n=1}^N (y_n - \boldgamma^\T \boldphi(\boldx_n))^2
$$ (eq-olslqphi)

Switching to matrix notation, if

$$
\boldPhi := 
    \left(
    \begin{array}{ccc}
        \phi_1(\boldx_1) & \cdots & \phi_J(\boldx_1)   \\
        \phi_1(\boldx_2) & \cdots & \phi_J(\boldx_2)   \\
        \vdots & \cdots & \vdots   \\
        \phi_1(\boldx_N) & \cdots & \phi_J(\boldx_N)   \\
    \end{array}
    \right) \in \RR^{N \times J}
$$ (eq-phidesign)

then the objective in {eq}`eq-olslqphi` can be expressed as $\| \boldy - \boldPhi \boldgamma \|^2$. Since increasing functions don't affect minimizers, the problem becomes

$$
\argmin_{\boldgamma \in \RR^J} \| \boldy - \boldPhi \boldgamma \|
$$ (eq-olsesphi)

Assuming that $\boldPhi$ is full column rank, the solution is

```{math}
:enumerated: false

\hboldgamma := (\boldPhi^\T \boldPhi)^{-1} \boldPhi^\T \boldy
```

````{prf:example}

Adding an intercept to a regression can be regarded as a transformation of the data. Indeed adding an intercept is equivalent to applying the transformation

```{math}
:enumerated: false

\boldphi(\boldx)
        =
        \left(
        \begin{array}{c}
            1 \\ 
            \boldx
        \end{array}
        \right) 
        =
        \left(
        \begin{array}{c}
            1 \\ 
            x_1 \\
            \vdots \\
            x_K
        \end{array}
        \right)
```

In practice, adding an intercept means fitting an extra parameter, and this extra degree of freedom allows a more flexible fit in our regression.
````


```{prf:example}
:label: eg-polymod

Let $K = 1$, so that $x_n \in \RR$. Consider the mononomial basis functions $\phi_j(x) := x^{j-1}$, so that

$$
\boldgamma^\T \boldphi(x_n)
        = \boldgamma^\T 
        \left(
        \begin{array}{c}
            x_n^0   \\
            x_n^1   \\
            \vdots  \\
            x_n^{J-1}  
        \end{array}
        \right)
        = \sum_{j = 1}^J \gamma_j x_n^{j-1}
$$ (eq-polymod)
The mononomial basis transformation applied to scalar $x$ corresponds to univariate polynomial regression, as discussed in {ref}`ss-ermls`. Under this transformation, the matrix $\boldPhi$ in {eq}`eq-phidesign` is called the **Vandermonde matrix**. By the Weierstrass approximation theorem, polynomials of sufficiently high order can effectively approximate any one-dimensional continuous nonlinear relationship.
```


```{prf:example}

{prf:ref}`eg-polymod` uses mononomials as the basis functions. A common alternative is to use orthogonal polynomials such as Chebychev polynomials or Hermite polynomials. Other alternatives include wavelets and splines. In econometrics this procedure is often referred to as nonparametric series regression. A key topic is the optimal number of basis functions.[^3]
```


{numref}`f-lap` and {numref}`f-lap3d` help illustrate how transformations can reduce approximation error. In [figure %s](#f-lap) it is clear that no linear function mapping $x$ to $y$ can produce small approximation error. {numref}`f-lap3d` shows the data after applying the transformation $\RR \ni x \mapsto \boldphi(x) := (x, x^3)^\T \in \RR^2$. The plane drawn in [figure %s](#f-lap3d) represents a linear function $\ell \colon \RR^2 \to \RR$. The composition $\ell \circ \boldphi$ has low approximation error. The two figures illustrate how nonlinear data can become linear when projected into higher dimensions.

```{figure} figures/linear_after_proj.svg
:name: f-lap

Nonlinear relationship between $x$ and $y$
```

```{figure} figures/linear_after_proj3d.svg
:name: f-lap3d

Approximate linearity after projecting the data to $\RR^2$
```

Below, when $y$ is regressed on $\boldx$, we can imagine that the data have already been transformed, and $\boldx$ is the result. Hence we use $\boldX$ to denote the design matrix instead of $\Phi$ without loss of generality.

(ss-sap)=
### The Frisch–Waugh–Lovell Theorem
The Frisch–Waugh–Lovell (FWL) theorem yields an expression for an arbitrary sub-vector of the least squares estimator $\hboldbeta$ obtained by regressing $\boldy$ on $\boldX$. In stating the theorem, we continue with the assumptions of {ref}`s-olsfp`. Let $\boldy$ and $\boldX$ be given and let $\hboldbeta$ be as in {eq}`eq-thelsqe`. In addition, let $K_1$ be an integer with $1 \leq K_1 < K$, and let

- $\boldX_1$ be a matrix consisting of the first $K_1$ columns of $\boldX$,

- $\boldX_2$ be a matrix consisting of the remaining $K_2 := K - K_1$ columns,

- $\hboldbeta_1$ be the $K_1 \times 1$ vector consisting of the first $K_1$ elements of $\hboldbeta$.

- $\hboldbeta_2$ be the $K_2 \times 1$ vector consisting of the remaining $K_2$ elements of $\hboldbeta$,

- $\boldP_1 := \proj ( \colspace \boldX_1)$, and

- $\boldM_1 := \boldI - \boldP_1 =$ the corresponding residual projection.

```{prf:theorem}
:label: t-fwl

The vector $\hboldbeta_2$ satisfies

$$
\hboldbeta_2 
        = (\boldX_2^\T \boldM_1 \boldX_2)^{-1} \boldX_2^\T \boldM_1 \boldy
$$ (eq-fwle)
```


````{prf:proof}
:nonumber:

From {eq}`eq-orr` and the definitions above we have

```{math}
:enumerated: false

\boldy 
        = \boldX \hboldbeta + \boldM \boldy
        = \boldX_1 \hboldbeta_1 + \boldX_2 \hboldbeta_2 + \boldM \boldy
```

Premultiplying both sides of this expression by $\boldX_2^\T \boldM_1$, we obtain

$$
\boldX_2^\T \boldM_1 \boldy 
        = \boldX_2^\T \boldM_1 \boldX_1 \hboldbeta_1 +
            \boldX_2^\T \boldM_1 \boldX_2 \hboldbeta_2 + 
            \boldX_2^\T \boldM_1 \boldM \boldy
$$ (eq-mbmx)
The first term on the right-hand side is zero by [fact %s](#fa-anox). The last term is also zero because

```{math}
:enumerated: false

(\boldX_2^\T \boldM_1 \boldM \boldy)^\T 
        = \boldy^\T  \boldM^\T \boldM_1^\T \boldX_2
        = \boldy^\T  \boldM \boldM_1 \boldX_2
        = \boldy^\T  \boldM \boldX_2
        = \boldzero
```

In the first equality, we used the usual property of transposes ([fact %s](#fa-owtra)); in the second, we used symmetry of $\boldM$ and $\boldM_1$; in the third, we used [fact %s](#fa-subsub2); and in the fourth, we used [fact %s](#fa-anox) again (which tells us that $\boldM$ maps all columns of $\boldX$, and hence all columns of $\boldX_2$, to the zero vector).

In light of the above, {eq}`eq-mbmx` becomes $\boldX_2^\T \boldM_1 \boldy = \boldX_2 ^\T\boldM_1 \boldX_2 \hboldbeta_2$. To go from this equation to {eq}`eq-fwle`, we just need to check that the matrix premultiplying $\hboldbeta_2$ is nonsingular. The proof is left as an exercise ([ex. %s](#ex-fpfwl)). ◻
````


As [exercise %s](#ex-cnxfwl) asks you to show, the expression for $\hboldbeta_2$ in [theorem %s](#t-fwl) can be rewritten as

$$
\hboldbeta_2 
        = [(\boldM_1\boldX_2)^\T \boldM_1 \boldX_2]^{-1} (\boldM_1 \boldX_2)^\T \boldM_1 \boldy
$$ (eq-sxb2)

Close inspection of this formula confirms the following claim: There is another way to obtain $\hboldbeta_2$ besides just regressing $\boldy$ on $\boldX$ and then extracting the last $K_2$ elements: we can also regress $\boldM_1 \boldy$ on $\boldM_1 \boldX_2$ to produce the same result.

To get some feeling for what this means, let's look at a special case, where $\boldX_2$ is the single column $\col_K \boldX$, containing the observations on the $K$th regressor. To tie into this notation let's write $\boldX_1$ as $\boldX_{-K}$ to remind us that it stands for all columns of $\boldX$ expect the $K$th one, and similarly for $\boldM_1$. In view of the preceding discussion, the least squares estimate $\hat \beta_K$ can be found by regressing

$$
\tilde \boldy 
    := \boldM_{-K} \boldy = \text{residuals of regressing } \boldy \text{ on
    } \boldX_{-K}
$$ (eq-resregy)

on

$$
\tilde \boldx_K 
    := \boldM_{-K} \col_K \boldX 
        = \text{residuals of regressing
            } \col_K \boldX \text{ on } \boldX_{-K}
$$ (eq-resregx)

Loosely speaking, these two residual terms $\tilde \boldy$ and $\tilde \boldx_K$ can be thought of as the parts of $\boldy$ and $\col_K \boldX$ that are "not explained by" $\boldX_{-K}$. Thus, on an intuitive level, the process for obtaining the least squares estimate $\hat \beta_K$ is as follows:

(i) Remove effects of all other regressors from $\boldy$ and $\col_K \boldX$, producing $\tilde \boldy$ and $\tilde \boldx_K$.

(ii) Regress $\tilde \boldy$ on $\tilde \boldx_K$.

This is obviously different from the process for obtaining the coefficient of the vector $\col_K \boldX$ in a simple univariate regression, the latter being just

(i) Regress $\boldy$ on $\col_K \boldX$.

In words, the difference between the univariate least squares estimated coefficient of the $K$th regressor and the multiple regression least squares coefficient is that the multiple regression coefficient $\hat \beta_K$ measures the *isolated relationship* between $x_K$ and $y$, without taking into account indirect channels involving other variables.

We can illustrate this idea further with a small simulation. Suppose that

```{math}
:enumerated: false

y = x_1 + x_2 + u  
    \quad \text{where} \quad
    u \iidsim \nN(0,1)
```

If we generate $N$ independent observations from this model and regress $y$ on the observations of $(x_1, x_2)$, then, provided that $N$ is sufficiently large, the coefficients for $x_1$ and $x_2$ will both be close to unity (see [theorem %s](#t-cwa)). However, if we regress $y$ on $x_1$ alone, then the coefficient for $x_1$ will depend on the relationship between $x_1$ and $x_2$. For example:

``` r
> N <- 1000
> x1 <- runif(N)
> x2 = 10 * exp(x1) + rnorm(N)
> y <- x1 + x2 + rnorm(N)
> results <- lm(y ~ 0 + x1)
> results$coefficients
      x1 
      30.83076 
```

Here the coefficient for $x_1$ is much larger than unity, because an increase in $x_1$ tends to have a large positive effect on $x_2$, which in turn increases $y$. The coefficient in the univariate regression reflects this total effect.

(ss-sr)=
#### Application: Simple Regression
Here's an easy application of the FWL theorem: deriving the familiar expression for the slope coefficient in simple linear regression (see {ref}`sss-blpsa`) from the multivariate expression. Simple linear regression is a special case of multivariate regression when $\boldone$ is the first column of $\boldX$ and $K = 2$. In this section, the second column of $\boldX$ will be denoted by $\boldx$. As we saw in {eq}`eq-solse`, the least squares estimates are

```{math}
:enumerated: false

\hat \beta_2 = \frac{\sum_{n=1}^N (x_n - \bar x)(y_n - \bar y)}
                    {\sum_{n=1}^N (x_n - \bar x)^2}
    \quad \text{and} \quad
    \hat \beta_1 = \bar y - \hat \beta_2 \bar x
```

where $\bar x$ is the sample mean of $\boldx$ and $\bar y$ is the sample mean of $\boldy$. We can rewrite the slope coefficient $\hat \beta_2$ more succinctly as

$$
\hat \beta_2 = 
    [(\boldx - \bar x \boldone)^\T(\boldx - \bar x \boldone)]^{-1}
    (\boldx - \bar x \boldone)^\T(\boldy - \bar y \boldone)
$$ (eq-sre1)

By the FWL theorem (equation {eq}`eq-sxb2`), we also have

$$
\hat \beta_2 
        = [(\boldM_c\boldx)^\T \boldM_c \boldx]^{-1} (\boldM_c \boldx)^\T \boldM_c \boldy
$$ (eq-sre2)

where $\boldM_c$ is the residual projection associated with the linear subspace $S = \Span\{\boldone\}$, as defined in {eq}`eq-pczm`. For this residual projection $\boldM_c$ and any $\boldz$, we have $\boldM_c \boldz = \boldz - \bar z \boldone$. Hence the right-hand sides of {eq}`eq-sre1` and {eq}`eq-sre2` coincide.

(ss-cenob)=
### Centered Observations
Let's generalize the preceding discussion to the case where there are multiple nonconstant regressors. The only difference is that instead of one column $\boldx$ of observations on a single nonconstant regressor, we have a matrix $\boldX_2$ containing multiple columns, each a vector of observations on a nonconstant regressor. If the least squares estimate $\hboldbeta$ is partitioned into $(\hat \beta_1, \hboldbeta_2)$, then we can write

```{math}
:enumerated: false

\boldX \hboldbeta = \boldone \beta_1 + \boldX_2 \hboldbeta_2
```

Applying the FWL theorem (equation {eq}`eq-sxb2`) once more, we can write $\hboldbeta_2$ as

```{math}
:enumerated: false

\hboldbeta_2 
        = [(\boldM_c\boldX_2)^\T \boldM_c \boldX_2]^{-1} (\boldM_c \boldX_2)^\T \boldM_c \boldy
```

where $\boldM_c$ is the residual projection in {eq}`eq-pczm`. As we saw in the last section, $\boldM_c \boldy$ is $\boldy$ centered around its mean. Similarly, $\boldM_c\boldX_2$ is a matrix formed by taking each column of $\boldX_2$ and centering it around its mean. It follows that in a least squares regression with an intercept, the estimated coefficients of the nonconstant regressors are equal to the estimated coefficients of a zero-intercept regression performed after all variables have been centered around their mean.

Let's use some related ideas to discuss an alternative to the coefficient of determination introduced in {eq}`eq-rsd`. There are several versions of $R^2$ reported in common regression packages. One of these is so called centered $R^2$. The version in {eq}`eq-rsd` will henceforth be called the uncentered $R^2$ for clarity.

One motivation for introducing an alternative to uncentered $R^2$ is that it fails to be invariant to certain changes of units. While it is invariant to changes of units that involve rescaling of the regressand $\boldy$ (see [ex. %s](#ex-riscy)), it is not invariant to changes of units that involve addition or subtraction (actual inflation versus inflation in excess of a certain level, income versus income over a certain threshold, etc.) whenever $\boldX$ contains an intercept. {prf:ref}`ex-snir` asks you to prove this.

This is one reason many econometricians use the **centered $R^2$** rather than $R^2$, at least when the regression contains an intercept. For the purposes of this section, let's assume that this is the case (or, more generally, that $\boldone \in \colspace \boldX$). Centered $R^2$ is defined as

$$
R_c^2 
     := \frac{\| \boldP \boldM_c \boldy \|^2}{\| \boldM_c \boldy \|^2}
     = \frac{\|  \boldM_c \boldP \boldy \|^2}{\| \boldM_c \boldy \|^2}
$$ (eq-crs)

$\boldM_c$ is as defined in {eq}`eq-pczm`. The equality of the two expressions for $R^2_c$ is left as an exercise ([ex. %s](#ex-crs)). Adding a constant to each element of $\boldy$ will have no effect on $R^2_c$ because $\boldM_c$ maps constant vectors to $\boldzero$ (see [example %s](#eg-pvonepm)).

Centered $R^2$ can be rewritten ([ex. %s](#ex-cecrs)) as

$$
R_c^2 = \frac{ \sum_{n=1}^N (\hat y_n - \bar y )^2 }
            { \sum_{n=1}^N (y_n - \bar y )^2 }
$$ (eq-cecrs)

It is a further exercise ([ex. %s](#ex-crsicc)) to show that, in the case of the simple regression, the $R_c^2$ is equal to the square of the sample correlation between the regressor and regressand, as defined in {eq}`eq-sampcorr`. Thus $R_c^2$ is a measure of correlation.

## Further Reading

For additional references on the material covered in this chapter see, for example, {cite}`friedman2001elements`, {cite}`ruud2000introduction`, {cite}`cameron2005microeconometrics`, or {cite}`davidson2004econometric`.

## Exercises

```{exercise}
:label: ex-fposd

Let $\boldx := (z, az)^\T$, where $a \in \RR$ and $z$ is any scalar random variable. Show that $\EE \boldx \boldx^\T$ is nonnegative definite but fails to be positive definite.
```


```{exercise}
:label: ex-riscy

Show that uncentered $R^2$ is invariant to changes of units that involve rescaling of the regressand $\boldy$ (dollars versus cents, kilometers versus miles, etc.)
```


```{exercise}
:label: ex-snir

Fix $\boldX, \boldy$ and consider regressing $\boldy$ on $\boldX$. Suppose that $\boldX$ contains the intercept, in the sense that $\boldone \in \colspace \boldX$. Let $R^2$ represent the uncentered coefficient of determination. Let $R^2_\alpha$ represent the same when $\boldy$ is replaced by $\boldy + \alpha \boldone$. Show that $R^2_\alpha \to 1$ as $\alpha \to \infty$.
```


```{exercise}
:label: ex-cdirr

Show that $R^2$ is invariant to a rescaling of the regressors.
```


```{exercise}
:label: ex-crsicc

Show that, in the case of the simple regression model in {ref}`ss-sr`, $R_c^2$ is equal to the square of the sample correlation between $\boldx$ and $\boldy$.
```


```{exercise}
:label: ex-crs

Confirm the equality of the two alternative expressions for $R_c^2$ in {eq}`eq-crs`.
```


```{exercise}
:label: ex-cecrs

Verify the expression for $R_c^2$ in {eq}`eq-cecrs`.
```


```{exercise}
:label: ex-rsqe

Prove [fact %s](#fa-rsqe).
```


```{exercise}
:label: ex-veseeo

Verify expression {eq}`eq-seeo0`.
```


```{exercise}
:label: ex-obolsd

Let's show that $\hboldbeta$ solves the least squares problem in a slightly different way: Let $\boldb$ be any $K \times 1$ vector, and let $\hboldbeta := (\boldX^\T \boldX)^{-1} \boldX^\T \boldy$.

(i) Show that $\| \boldy - \boldX \boldb \|^2 = \| \boldy - \boldX \hboldbeta \|^2 + \| \boldX(\hboldbeta - \boldb) \|^2$.

(ii) Using (i), argue that $\hboldbeta$ is the minimizer of $\| \boldy - \boldX \boldb \|^2$ over all $K \times 1$ vectors $\boldb$.
```


```{exercise}
:label: ex-ch11_regression-auto-1

Verify that $\sum_{n=1}^N (y_n - \boldb^\T \boldx_n)^2 = \| \boldy - \boldX \boldb \|^2$.
```


```{exercise}
:label: ex-soavv

Show carefully that any solution to $\min_{\boldb \in \RR^K} \| \boldy - \boldX \boldb \|^2$ is also a solution to $\min_{\boldb \in \RR^K} \| \boldy - \boldX \boldb \|$, and vice versa.
```


```{exercise}
:label: ex-wpoo

Confirm that $\boldP \boldone = \boldone$ whenever $\boldone \in \colspace \boldX$.
```


```{exercise}
:label: ex-icime0

Show that, for any regression containing the intercept, the vector of residuals must sum to zero.
```


```{exercise}
:label: ex-icime

Show that, for any regression containing the intercept, the mean of the fitted values $\hboldy = \boldP \boldy$ is equal to the mean of $\boldy$.
```


```{exercise}
:label: ex-ch11_regression-auto-2

Show that $\boldP \boldM = \boldM \boldP = \boldzero$. Using this fact (instead of the orthogonal projection theorem), show that the vector of fitted values and the vector of residuals are orthogonal.
```


```{exercise}
:label: ex-yibx

Show that if $R^2 = 1$, then the vector of residuals is identically zero, $\boldP \boldy = \boldy$, and $\boldy \in \colspace \boldX$.
```


```{exercise}
:label: ex-ch11_regression-auto-3

Suppose that the regression contains an intercept. Let $\bar y$ be the sample mean of $\boldy$, and let $\bar \boldx$ be a $1 \times K$ row vector such that the $k$th element of $\bar \boldx$ is the sample mean of the $k$th column of $\boldX$. Show that $\bar y = \bar \boldx \hboldbeta$.
```


```{exercise}
:label: ex-myae

Suppose the regression contains an intercept. Let $\boldM_c$ be as defined in {eq}`eq-pczm`. Show that the following identity always holds:

$$
\| \boldM \boldy \|^2 = \| \boldM_c \boldy \|^2 - \| \boldP \boldM_c \boldy \|^2
$$ (eq-rsamy)
```


```{exercise}
:label: ex-scsp

Let $\boldX_a$ and $\boldX_b$ be $N \times K_a$ and $N \times K_b$ respectively. Suppose that every column of $\boldX_a$ is also a column of $\boldX_b$. Show that $\colspace \boldX_a \subset \colspace\boldX_b$.
```


````{exercise}
:label: ex-nesc

Let $\boldx := (x_1,\ldots,x_N)$ and $\boldy := (y_1,\ldots,y_N)$ be sequences of scalar random variables. Show that the sample correlation $\hat \rho$ between $\boldx$ and $\boldy$ (defined in {eq}`eq-sampcorr`) can be written as

```{math}
:enumerated: false

\hat \rho = \frac{  (\boldM_c \boldx)^\T (\boldM_c \boldy)  }
            { \| \boldM_c \boldx \| \| \boldM_c \boldy \| }
```
````


```{exercise}
:label: ex-cnxfwl

Show that the two expressions for $\hboldbeta_2$ in {eq}`eq-fwle` and {eq}`eq-sxb2` are equal.[^4]
```


```{exercise}
:label: ex-fpfwl

At the end of the proof of [theorem %s](#t-fwl), it was claimed that the matrix $\boldX_2^\T \boldM_1 \boldX_2$ is nonsingular. Verify this claim.
```


```{exercise}
:label: ex-ch11_regression-auto-4

Build an arbitrary data set $\boldX$, $\boldy$ by simulation. Run a regression with the intercept, and record the values of the estimated coefficients of the nonconstant (i.e., $k \geq 2$) regressors. Confirm that these values are equal to the estimated coefficients of the no-intercept regression after all variables have been centered around their mean.
```


### Solutions to Selected Exercises

````{solution} ex-fposd

The expression $\EE \boldx \boldx^\T$ is nonnegative definite for any $\boldx$ because, given $\boldc \in \RR^K$, we have

```{math}
:enumerated: false

\boldc^\T \EE \boldx \boldx^\T \boldc
        = \EE (\boldc^\T \boldx) (\boldx^\T \boldc)
        = \EE (\boldx^\T \boldc)^2 \geq 0
```

Here we used the fact that the transpose of a scalar is equal to the scalar. However, if $\boldx := (z, az)^\T$, then

```{math}
:enumerated: false

\EE \boldx \boldx^\T
        = \EE z^2 
        \begin{pmatrix}
            1  & a \\
            a  & a^2
        \end{pmatrix}
```

The second column is a multiple of the first, so the matrix is singular and the determinant is zero. As such it cannot be positive definite (see [fact %s](#fa-ipde)).
````


````{solution} ex-riscy

If $\boldy$ is scaled by $\alpha \in \RR$, then

```{math}
:enumerated: false

\frac{\| \boldP \alpha \boldy \|^2}{\| \alpha \boldy \|^2}
         = \frac{\| \alpha \boldP \boldy \|^2}{\| \alpha \boldy \|^2}
         = \frac{\alpha^2 \| \boldP \boldy \|^2}{\alpha^2 \| \boldy \|^2}
         = \frac{\| \boldP \boldy \|^2}{\| \boldy \|^2}
```

Hence the uncentered $R^2$ is unchanged.
````


````{solution} ex-snir

Fix $\alpha \in \RR$. By definition, $R^2_\alpha$ is

```{math}
:enumerated: false

\frac{\| \boldP (\boldy + \alpha \boldone) \|^2}
                         {\| \boldy + \alpha \boldone \|^2}
         = \frac{\| \boldP \boldy + \alpha \boldP \boldone \|^2}
                         {\| \boldy + \alpha \boldone \|^2}
         = \frac{\| \boldP \boldy + \alpha \boldone \|^2}
                         {\| \boldy + \alpha \boldone \|^2}
         = \frac{\alpha^2 \| \boldP \boldy / \alpha + \boldone \|^2}
                         {\alpha^2 \| \boldy / \alpha + \boldone \|^2}
         = \frac{ \| \boldP \boldy / \alpha + \boldone \|^2}
                         { \| \boldy / \alpha + \boldone \|^2}
```

where the second inequality follows from the fact that $\boldone \in \colspace \boldX$. Taking the limit as $\alpha \to \infty$, we find that $R^2_\alpha \to 1$ as $\alpha \to \infty$.
````


````{solution} ex-cdirr

This follows immediately from the definition of $R^2$, and the fact that, for any $\alpha \neq 0$,

```{math}
:enumerated: false

\boldP 
        = \boldX (\boldX^\T\boldX)^{-1} \boldX^\T
        = \frac{\alpha^2}{\alpha^2} \boldX (\boldX^\T\boldX)^{-1} \boldX^\T
        = (\alpha \boldX) ( (\alpha \boldX)^\T(\alpha \boldX))^{-1} (\alpha \boldX^\T)
        
```
````


````{solution} ex-crsicc

From [exercise %s](#ex-nesc), the squared sample correlation between $\boldx$ and $\boldy$ can be written as

```{math}
:enumerated: false

\hat \rho^2 
        = \frac{ [ (\boldM_c \boldx)^\T (\boldM_c \boldy) ]^2 }
            { \| \boldM_c \boldx \|^2 \| \boldM_c \boldy \|^2 }
```

Also, $R_c^2 = \|  \boldM_c \boldP \boldy \|^2/\| \boldM_c \boldy \|^2$. Hence it suffices to show that, for the simple linear regression model in {ref}`ss-sr`, we have

$$
\|  \boldM_c \boldP \boldy \| 
        = \frac{ |(\boldM_c \boldx)^\T (\boldM_c \boldy)|  }
            { \| \boldM_c \boldx \| }
$$ (eq-seos)
Let $\boldX = (\boldone, \boldx)$ be the design matrix, where the first column is $\boldone$ and the second column is $\boldx$. Let

```{math}
:enumerated: false

\hat \beta_1 := \bar y - \hat \beta_2 \bar x
        \quad \text{and} \quad
        \hat \beta_2 := \frac{ (\boldM_c \boldx)^\T (\boldM_c \boldy)  }
            { \| \boldM_c \boldx \|^2 }
```

be the least squares estimators of $\beta_1$ and $\beta_2$ respectively (see {ref}`ss-sr`). We then have

```{math}
:enumerated: false

\boldP \boldy 
        = \boldX \hboldbeta 
        = \boldone \hat \beta_1 + \boldx \hat \beta_2
```

```{math}
:enumerated: false

\fore
        \boldM_c \boldP \boldy 
        = \boldM_c \boldX \hboldbeta 
        = \boldM_c \boldx \hat \beta_2
```

```{math}
:enumerated: false

\fore
        \| \boldM_c \boldP \boldy \|
        = \| \boldM_c \boldx \hat \beta_2 \|
        = |\hat \beta_2| \| \boldM_c \boldx  \|
        = \frac{ |(\boldM_c \boldx)^\T (\boldM_c \boldy)|  }
            { \| \boldM_c \boldx \|^2 } \| \boldM_c \boldx  \|
```

Canceling $\|\boldM_c \boldx  \|$ we get {eq}`eq-seos`. This completes the proof.
````


````{solution} ex-crs

It is sufficient to show that $\boldP \boldM_c = \boldM_c \boldP$. Since $\boldone \in \colspace \boldX$ by assumption, we have $\boldP \boldone = \boldone$, and $\boldone^\T \boldP = (\boldP^\T \boldone)^\T = (\boldP \boldone)^\T = \boldone^\T$. Therefore $\boldP \boldone \boldone^\T = \boldone \boldone^\T \boldP$, and

```{math}
:enumerated: false

\boldP \boldM_c 
        = \boldP - \frac{1}{N} \boldP \boldone \boldone^\T 
        = \boldP - \frac{1}{N} \boldone \boldone^\T \boldP 
        = \boldM_c \boldP 
        
```
````


````{solution} ex-cecrs

It suffices to show that $\| \boldP \boldM_c \boldy \|^2 = \sum_{n=1}^N (\hat y_n - \bar y )^2$. This holds because

```{math}
:enumerated: false

\sum_{n=1}^N (\hat y_n - \bar y )^2 
        = \| \boldP \boldy - \boldone \bar y \|^2
        = \| \boldP \boldy - \boldP \boldone \bar y \|^2
        = \| \boldP (\boldy - \boldone \bar y) \|^2
        = \| \boldP \boldM_c \boldy \|^2
        
```
````


````{solution} ex-rsqe

Adopt the notation and assumptions of [fact %s](#fa-rsqe). Let $\boldy$ be given. Let $\boldP_a$ and $\boldP_b$ be the projections onto the column spaces of $\boldX_a$ and $\boldX_b$ respectively, so that $R^2_i = \| \boldP_i \boldy \|^2/\| \boldy \|^2$ for $i \in \{a, b\}$. From [fact %s](#fa-subsub) we have $\boldP_a \boldP_b \boldy = \boldP_a \boldy$. Using this fact, and setting $\boldy_b := \boldP_b \boldy$, gives

```{math}
:enumerated: false

\frac{R^2_a}{R^2_b} 
        = \left( \frac{\| \boldP_a \boldy \|}{\| \boldP_b \boldy \|} \right)^2
        = \left( \frac{\| \boldP_a \boldP_b \boldy \|}{\| \boldP_b \boldy \|} \right)^2
        = \left( \frac{\| \boldP_a \boldy_b \|}{\| \boldy_b \|} \right)^2
        \leq 1
```

where the final inequality follows from [theorem %s](#t-opt2). Hence $R^2_b \geq R^2_a$, and regressing with $\boldX_b$ produces (weakly) larger $R^2$.
````


````{solution} ex-obolsd

Part (ii) follows immediately from part (i). Regarding part (i), observe that

```{math}
:enumerated: false

\| \boldy - \boldX \boldb \|^2 
            = \| \boldy - \boldX \hboldbeta + \boldX(\hboldbeta - \boldb) \|^2
```

By the Pythagorean law, the claim

```{math}
:enumerated: false

\| \boldy - \boldX \boldb \|^2 
            = \| \boldy - \boldX \hboldbeta \|^2 
            + \| \boldX(\hboldbeta - \boldb) \|^2
```

will be confirmed if $\boldy - \boldX \hboldbeta \; \perp\; \boldX(\hboldbeta - \boldb)$. This follows from the definition of $\hboldbeta$, because for arbitrary $\bolda \in \RR^K$ we have

```{math}
:enumerated: false

(\boldX \bolda)^\T (\boldy - \boldX \hboldbeta)
        = \bolda^\T (\boldX^\T \boldy - \boldX^\T \boldX (\boldX^\T\boldX)^{-1}
        \boldX^\T \boldy)
        = \bolda^\T (\boldX^\T \boldy - \boldX^\T \boldy)
        = 0
        
```
````


````{solution} ex-soavv

Let $\hboldbeta$ be a solution to $\min_{\boldb \in \RR^K} \| \boldy - \boldX \boldb \|^2$, which is to say that

```{math}
:enumerated: false

\| \boldy - \boldX \hboldbeta \|^2 \leq \| \boldy - \boldX \boldb \|^2
        \quad \text{for any } \boldb \in \RR^K
```

If $a$ and $b$ are nonnegative constants with $a \leq b$, then $\sqrt{a} \leq \sqrt{b}$, and hence

```{math}
:enumerated: false

\| \boldy - \boldX \hboldbeta \| \leq \| \boldy - \boldX \boldb \|
        \quad \text{for any } \boldb \in \RR^K
```

In other words, $\hboldbeta$ is a solution to $\min_{\boldb \in \RR^K} \| \boldy - \boldX \boldb \|$. The "vice versa" argument follows along similar lines.
````


````{solution} ex-icime0

To see why the vector of residuals sums to zero, observe that

```{math}
:enumerated: false

\boldone^\T \boldM \boldy
        = \boldone^\T \boldy - \boldone^\T\boldP \boldy
        = \boldone^\T \boldy - (\boldP^\T\boldone)^\T \boldy
        = \boldone^\T \boldy - (\boldP\boldone)^\T \boldy
```

where the last equality uses the fact the $\boldP$ is symmetric. Moreover, by [exercise %s](#ex-wpoo), we have $\boldP\boldone = \boldone$ whenever $\boldone$ is a column of $\boldX$. Therefore

```{math}
:enumerated: false

\boldone^\T \boldM \boldy
        = \boldone^\T \boldy - (\boldP\boldone)^\T \boldy
        = \boldone^\T \boldy - \boldone^\T \boldy
        = 0
        
```
````


````{solution} ex-icime

Since $\boldone \in \colspace \boldX$, we have $\boldP\boldone = \boldone$. It follows that

```{math}
:enumerated: false

\frac{1}{N} \sum_{n=1}^N \hat y_n
        = \frac{1}{N} \boldone^\T \hboldy
        = \frac{1}{N} \boldone^\T \boldP \boldy
        = \frac{1}{N} \boldone^\T \boldy 
        = \frac{1}{N} \sum_{n=1}^N y_n
        
```
````


```{solution} ex-yibx

If $R^2 = 1$, then, by {eq}`eq-bdtss`, we have $\| \boldM \boldy \|^2 = 0$, and hence have $\hat \boldu = \boldM \boldy  = \boldzero$. Since $\boldM \boldy + \boldP \boldy = \boldy$, this implies that $\boldP \boldy = \boldy$. But then $\boldy \in \colspace \boldX$ by (vi) of [theorem %s](#t-opt2).
```


````{solution} ex-myae

[Fact %s](#fa-opt3) tells us that for any conformable vector $\boldz$ we have $\boldz = \boldP \boldz + \boldM \boldz$, where the two vectors on the right-hand side are orthogonal. Letting $\boldz = \boldM_c \boldy$, we obtain

```{math}
:enumerated: false

\boldM_c \boldy = \boldP \boldM_c \boldy + \boldM \boldM_c \boldy
```

From [fact %s](#fa-subsub2) we have $\boldM \boldM_c \boldy = \boldM \boldy$. Using this result, orthogonality, and the Pythagorean law, we obtain

```{math}
:enumerated: false

\| \boldM_c \boldy \|^2 
        = \| \boldP \boldM_c \boldy \|^2 + \| \boldM \boldy \|^2
```

Rearranging gives {eq}`eq-rsamy`
````


```{solution} ex-scsp

This is just [fact %s](#fa-xsys).
```


````{solution} ex-fpfwl

By [fact %s](#fa-ipde), to show that $\boldX_2^\T \boldM_1 \boldX_2$ is nonsingular, it suffices to show that the matrix is positive definite. By idempotence and symmetry of $\boldM_1$,

```{math}
:enumerated: false

\boldX_2^\T \boldM_1 \boldX_2  
        = \boldX_2^\T \boldM_1\boldM_1 \boldX_2  
        = \boldX_2^\T \boldM_1^\T\boldM_1 \boldX_2  
        = (\boldM_1\boldX_2)^\T \boldM_1 \boldX_2
```

Take any $\bolda \neq \boldzero$. We need to show that

```{math}
:enumerated: false

\bolda^\T (\boldM_1\boldX_2)^\T \boldM_1 \boldX_2 \bolda
        = (\boldM_1\boldX_2 \bolda)^\T \boldM_1 \boldX_2 \bolda 
        = \| \boldM_1\boldX_2 \bolda \|^2 > 0
```

Since the only vector with zero norm is the zero vector, it now suffices to show that $\boldM_1\boldX_2 \bolda$ is nonzero. From [fact %s](#fa-opt3), we see that $\boldM_1\boldX_2 \bolda = \boldzero$ only when $\boldX_2 \bolda$ is in the column span of $\boldX_1$. Thus, the proof will be complete if we can show that $\boldX_2 \bolda$ is not in the column span of $\boldX_1$.

Indeed, $\boldX_2 \bolda$ is not in the column span of $\boldX_1$. If it were, then we could write $\boldX_1 \boldb = \boldX_2 \bolda$ for some $\boldb \in \RR^{K_1}$. Rearranging, we get $\boldX \boldc = \boldzero$ for some nonzero $\boldc$ (recall $\bolda \neq \boldzero$). This contradicts linear independence of the columns of $\boldX$.
````


[^1]: In essence, positive definiteness of $\Sigmax$ requires that no random variable in $\boldx$ can be written as a linear combination of other variables in $\boldx$. See [exercise %s](#ex-fposd).

[^2]: Later, in {ref}`ss-lslse`, we'll do something similar (i.e., weaken independence to ergodicity) in a setting with some additional structure.

[^3]: See, for example, {cite}`hong1995consistent`, {cite}`sun2011robust`, or {cite}`chen2015optimal`.

[^4]: Hint: Use the symmetry and idempotence of the matrix $\boldM_1$.
