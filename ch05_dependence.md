---
title: "Modeling Dependence"
---

To study relationships between random variables we need to add more concepts to our toolkit. In this chapter we review multivariate distributions and give a rigorous construction of conditional expectation.

## Random Vectors and Matrices

Let's begin with random vectors and matrices.

(ss-rvandm)=
### Random Vectors
In {ref}`ss-rvs` we defined a random variable to be a function $x \colon \Omega \to \RR$ with the property that $\{x \in B\}$ is an event for any Borel set $B \in \bB(\RR)$. The definition of a random vector is analogous. A **random vector** $\boldx$ in $\RR^N$ is a function from $\Omega$ to $\RR^N$ with the property that

$$
\setntn{\omega \in \Omega}{\boldx(\omega) \in B} \in \fF
    \quad \text{for all } B \in \bB(\RR^N)
$$ (eq-dorvv)

This is just {eq}`eq-dorv` with $\bB(\RR^N)$ replacing $\bB(\RR)$. In statistical applications we can almost always cook up a sample space $\Omega$ and a set of events $\fF$ on $\Omega$ such that all random vectors we care about satisfy {eq}`eq-dorvv`.[^1]

It's also possible to define a **random vector** $\boldx$ in $\RR^N$ as a tuple of $N$ random variables $(x_1, \ldots, x_N)$. The two definitions are equivalent. As with ordinary vectors, random vectors have no intrinsic shape, and we'll write them in rows or columns according to convenience. However, in expressions involving matrix multiplication, random vectors will default to column vectors unless otherwise stated.

```{prf:example}

Recall the blindfolded monkey experiment. The sample space is the unit disk $\Omega := \setntn{(h, v) \in \RR^2}{\| (h, v) \| \leq 1}$ and the event space is the Borel sets in $\Omega$. If $\boldx$ is the identity on $\Omega$, then it simply reports the outcome $(h,v)$. This is one instance of a random vector.
```


```{prf:example}

Consider a random sample listing the income $y_n$ of $n=1,\ldots,N$ individuals from a given population. The vector $(y_1, \ldots, y_N)$ that reports the outcome of this sampling can be regarded as a random vector in $\RR^N$.
```


The definition of a random vector ensures that $\{\boldx \in B\}$ is a well-defined event for every $B \in \bB(\RR^N)$. If we want to ensure that some transformation $\boldy = f(\boldx)$ is likewise a random vector, we need a restriction on $f$ so that $\{\boldy \in B\}$ is also a well-defined event. The appropriate restriction is $\bB$-measurability, which, analogous to the scalar case , means that $f \colon \RR^N \to \RR^M$ satisfies $f^{-1}(B) \in \bB(\RR^N)$ for all $B \in \bB(\RR^M)$. As in the scalar case, every kind of transformation we use in this text is $\bB$-measurable, so the restriction is never binding.

Expectations of random vectors are defined element-by-element. In particular, if $\boldx = (x_1, \ldots, x_N)$ is a random vector in $\RR^N$, then its **expectation** $\EE \boldx$ is defined as the vector of expectations

```{math}
:enumerated: false

\EE \boldx
    = 
    \EE
    \left(
    \begin{array}{c}
        x_1 \\
        x_2 \\
        \vdots \\
        x_N
    \end{array}
    \right)
    :=
    \left(
    \begin{array}{c}
        \EE x_1 \\
        \EE x_2 \\
        \vdots \\
        \EE x_N
    \end{array}
    \right)
```

An $M \times N$ **random matrix** $\boldX$ is an $M \times N$ array of random variables. Its expectation is defined as

```{math}
:enumerated: false

\EE \boldX
    := 
    \begin{pmatrix}
        \EE x_{11} & \cdots & \EE x_{1N} \\
        \vdots &  & \vdots \\
        \EE x_{M1} & \cdots & \EE x_{MN}
    \end{pmatrix}
```

````{prf:property}
:label: fa-lemc

If $\boldX$ and $\boldY$ are random matrices or vectors and $\boldA$ and $\boldB$ are constant and conformable, then

```{math}
:enumerated: false

\EE [ \boldA \boldX + \boldB \boldY ]
            =  \boldA \EE [\boldX ] + \boldB \EE[\boldY]
```
````


[Fact %s](#fa-lemc) is obtained by repeated application of the linearity of scalar expectations stated in [fact %s](#fa-exprop).

The **variance–covariance matrix** of a random vector $\boldx$ in $\RR^N$ with $\boldmu := \EE \boldx$ is the $N \times N$ matrix

$$
\var [\boldx] := \EE [ (\boldx - \boldmu) (\boldx - \boldmu)^\T]
$$ (eq-vcd1)

Expanding this out gives

```{math}
:enumerated: false

\var [\boldx]
    = 
    \left(
    \begin{array}{ccc}
        \EE [(x_1 - \mu_1)(x_1 - \mu_1)]
            & \cdots & \EE [(x_1 - \mu_1)(x_N - \mu_N)] \\
        \vdots & & \vdots \\
        \EE [(x_N - \mu_N)(x_1 - \mu_1)]  
            & \cdots & \EE [(x_N - \mu_N)(x_N - \mu_N)] \\
    \end{array}
    \right)
```

The $j,k$th term is the scalar covariance between $x_j$ and $x_k$ and the principal diagonal contains the variance of each $x_n$.

```{prf:property}
:label: fa-ltavc

For any random vector $\boldx$ with $\EE[ \boldx^\T\boldx ] < \infty$,

(i) $\var [\boldx]$ exists and is nonnegative definite,

(ii) $\var[\boldx] = \EE [\boldx \boldx^\T] - \boldmu \boldmu^\T$, and

(iii) $\var[\boldA \boldx + \boldb] = \boldA \var[\boldx] \boldA^\T$ (any $\boldA, \boldb$ constant and conformable).
```


The **cross-covariance** between random vectors $\boldx$ and $\boldy$ is defined as

```{math}
:enumerated: false

\cov [\boldx, \boldy]
    := \EE [ (\boldx - \EE [\boldx]) (\boldy - \EE [\boldy])^\T]
```

Evidently $\var[\boldx] = \cov[\boldx, \boldx]$.

Here's a useful fact concerning expectations of quadratic forms. The proof is a solved exercise (see [ex. %s](#ex-eoqf)).

````{prf:property}
:label: fa-eoqf

If $\boldz$ is a random vector in $\RR^N$ satisfying $\EE[\boldz \boldz^\T] = \boldI$ and $\boldA$ is any constant $N \times N$ matrix, then

```{math}
:enumerated: false

\EE[\boldz^\T \boldA \boldz] = \trace \boldA
```
````


(ss-muvd)=
### Multivariate Distributions
The definition of a distribution over $\RR^N$ generalizes that of a distribution over $\RR$ given in {ref}`ss-repdist`. In particular, a **distribution** or **law** $P$ on $\RR^N$ is a probability measure over the Borel sets $\bB(\RR^N)$. By definition, it satisfies $P(\RR^N) = 1$ and $P(\cup_{n=1}^{\infty} B_n) = \sum_{n=1}^\infty P(B_n)$ for any disjoint sequence $\{B_n\}$ in $\bB(\RR^N)$. See also the properties of probability measures listed in [fact %s](#fa-bpro).

{numref}`f-gaussian_example` shows a visualization of a distribution $P$ on $\RR^2$, with hot colors indicating regions of high probability. (The distribution is multivariate Gaussian, to be defined below.) The sets $A$ and $B$ are typical Borel sets. In the figure, $B$ has the same area as $A$ but is assigned greater probability by $P$. By additivity, we know that $P(A \cup B) = P(A) + P(B)$ must hold.

```{figure} figures/gaussian_example.svg
:name: f-gaussian_example

Example distribution and events $A$ and $B$
```

As was the case for $\RR$, distributions on $\RR^N$ are fully characterized by their values on certain more elementary families of subsets of $\RR^N$. For example, for any two distributions $Q, P$ on $\RR^N$,

```{math}
:enumerated: false

P(I) = Q(I) 
    \text{ for all } I \in \rR
    \quad \implies \quad
    P(B) = Q(B) 
    \text{ for all } B \in \bB(\RR^N)
```

where $\rR$ denotes the family of all **cylinder sets**

```{math}
:enumerated: false

C = B_1 \times \cdots \times B_N, \qquad B_n \in \bB(\RR), \; n=1,\ldots, N
```

or, specializing further, the family of all rectangles[^2]

```{math}
:enumerated: false

I = (-\infty, s_1] \times \cdots \times (-\infty, s_N], 
    \qquad s_n \in \RR, \; n=1,\ldots, N
```

(See prop 3.7 of {cite:author}`cinlar` {cite:year}`cinlar`.) {numref}`f-gaussian_example_2` shows an example set $I = (-\infty, s_1] \times (-\infty, s_2]$ in the case of $\RR^2$.

```{figure} figures/gaussian_example_2.svg
:name: f-gaussian_example_2

Distributions are characterized by their values on elements of $\rR$
```

Thus, analogous to {eq}`eq-ntf`, any distribution $P$ on $\RR^N$ is characterized by the function

$$
F(\bolds) :=
    F(s_1, \ldots, s_N) 
    := P \left( \times_{n=1}^N (-\infty, s_n] \right)
    \qquad (\bolds \in \RR^N)
$$ (eq-ntf2)

This is an example of a **multivariate cumulative distribution function**, which is a function $F \colon \RR^N \to [0, 1]$ that is

(i) right-continuous in each of its arguments,

(ii) increasing in each of its arguments, and

(iii) satisfies

```{math}
:enumerated: false

F(\bolds_j) \to 1 \text{ as }
                \bolds_j \to \infty 
                \quad \text{and} \quad
                F(s_1, \ldots, s_{nj}, \ldots, s_N) \to 0
                \text{ as }
                s_{nj} \to -\infty
```

As in the one-dimensional case, we will call a distribution $P$ on $\RR^N$ **discrete** if $P$ is supported on a countable subset of $\RR^N$, and **absolutely continuous** if $P(B) = 0$ whenever $B$ has zero Lebesgue measure. The latter condition is again necessary and sufficient for existence of a density representation, by which we mean existence of a density $p$ on $\RR^N$ such that

$$
P(B) = \int_B p(\bolds) \, \diff  \bolds
    \qquad \text{for all } B \in \bB(\RR^N)
$$ (eq-drbd3)

The right-hand side of {eq}`eq-drbd3` is a multivariate integral that we can write out more explicitly as

$$
\int_{-\infty}^{\infty}
        \cdots
        \int_{-\infty}^{\infty} 
        \1_B (s_1,\ldots,s_N)
        p(s_1,\ldots,s_N) 
        \diff s_1 \cdots \diff s_N
$$ (eq-defjd0)

Conversely, if $p$ is any density on $\RR^N$, then {eq}`eq-drbd3` defines a distribution.

```{prf:example}
:label: eg-mnden

The **multivariate normal density** or **multivariate Gaussian density** on $\RR^N$ is a function $p$ of the form

$$
p(\bolds) = (2 \pi)^{-N/2} \det(\boldSigma)^{-1/2} 
        \exp \left\{ 
            - \frac{1}{2} (\bolds - \boldmu)^\T \boldSigma^{-1} (\bolds - \boldmu) 
        \right\}
$$ (eq-mnormden)
where $\boldmu$ is any $N \times 1$ vector and $\boldSigma$ is a positive definite $N \times N$ matrix. In symbols, we represent this distribution by $\nN(\boldmu, \boldSigma)$. The case $\nN(\boldzero, \boldI)$ is called the **multivariate standard normal distribution**. {numref}`f-bivar_gaussian_3d` shows the standard normal density when $N=2$.
```


```{figure} figures/bivar_gaussian_3d.svg
:name: f-bivar_gaussian_3d

Bivariate standard normal density
```

Let's describe one important way to construct multivariate distributions.

````{prf:property}
:label: fa-prodist

Given distributions $P_1, \ldots, P_N$ on $\RR$, there exists a unique and well-defined distribution $\mathring{P}$ on $\RR^N$ such that

```{math}
:enumerated: false

\mathring{P}(B_1 \times \cdots \times B_N)
        = \prod_{n=1}^N P_n(B_n)
        \qquad \text{for all } B_n \in \bB(\RR), \; n=1,\ldots, N
```
````


This distribution is called the **product distribution** of $P_1, \ldots, P_N$. If $P_n = P$ for all $n$, then we'll call $\mathring P$ the **$N$th product of $P$**. The uniqueness part of [fact %s](#fa-prodist) follows from the claim above that distributions on $\RR^N$ are uniquely pinned down by their values on cylinder subsets of $\RR^N$. As we'll see, product distributions correspond to joint distributions of independent random variables.

(ss-marg)=
#### Marginals of Multivariate Distributions
Given any distribution $P$ on $\RR^N$, the $n$th **marginal distribution** of $P$ is the distribution on $\RR$ defined by

```{math}
:enumerated: false

P_n(B)
    = P (\RR \times \cdots \times \RR \times B \times \RR \times \cdots
    \times \RR)
```

Here $B$ is the $n$th element of the Cartesian product. We can also write this as

```{math}
:enumerated: false

P_n(B) 
    = P \setntn{\bolds \in \RR^N}{\bolds^\T \bolde_n \in B}
```

The $n$th marginal tells us the probabilistic implications of $P$ when we only concern ourselves with the outcomes for the $n$th coordinate. From $P_n$ we can also extract the **marginal CDF $F_n$** via {eq}`eq-ntf`.

If $P_n$ is absolutely continuous, it has density $p_n$. When the joint distribution $P$ has a density $p$, the marginal distribution $P_n$ has a density, and $p_n$ can be obtained directly from $p$ by integrating out the other variables. For example, in the bivariate case,

```{math}
:enumerated: false

p_1(s_1)
    = \int_{-\infty}^{\infty} p(s_1, s_2) \diff s_2
```

{numref}`f-jointplot` shows a bivariate density represented by contours and the two marginals obtained from projecting onto each axis.

```{figure} figures/jointplot.svg
:name: f-jointplot

Bivariate joint density and its two marginals
```

In general, the joint distribution cannot be determined from the marginals alone, since the marginals do not tell us about the interactions across coordinates. The exception is when there is no interaction. This is the case for product distributions, as described in [fact %s](#fa-prodist) and further discussed in {ref}`ss-indep`.

(ss-drvecs)=
### Distributions of Random Vectors
Let $\boldx$ be a random vector in $\RR^N$. The **distribution** of $\boldx$ is the probability measure $P$ on $\bB(\RR^N)$ defined by

$$
P(B) = \PP\{\boldx \in B\}
    \qquad (B \in \bB(\RR^N))
$$ (eq-dmnu)

In this context, $P$ is also called the **joint distribution** of $x_1,\ldots, x_N$, and we write $\lL(\boldx) = P$. The joint distribution can also be represented by the multivariate CDF $F \colon \RR^N \to [0, 1]$ defined by {eq}`eq-ntf2`, which in this case can be written as $F(s_1,\ldots,s_N) = \PP\{ x_1 \leq s_1, \ldots, x_N \leq s_N\}$, or, in vector notation,

```{math}
:enumerated: false

F(\bolds) = \PP\{\boldx \leq \bolds\}
    \qquad (\bolds \in \RR^N)
```

It follows from {eq}`eq-drbd3` that when the distribution $P$ of $\boldx$ is absolutely continuous, there exists a nonnegative function $p$ on $\RR^N$ satisfying

$$
\int_B p(\bolds) \diff  \bolds = \PP\{\boldx \in B\}
    \qquad (B \in \bB(\RR^N))
$$ (eq-defjd00)

Here $\int_B$ is a multivariate integral over the region $B$. The function $p$ is called either the density or **joint density** of $\boldx$. For {eq}`eq-defjd00` to hold, it suffices that

$$
\int_{-\infty}^{s_N} 
    \cdots
     \int_{-\infty}^{s_1} 
    p(t_1,\ldots,t_N) 
    \diff t_1 \cdots \diff t_N
    = F(s_1,\ldots,s_N)
$$ (eq-defjd)

for all $s_n  \in \RR$, $n=1,\ldots, N$.

#### Joint versus Marginal

If $\boldx = (x_1,\ldots, x_N)$ is a random vector in $\RR^N$, then each $x_n$ is a random variable on $\RR$. Let $P_n = \lL(x_n)$, so that

$$
P_n(B) = \PP \{x_n \in B \}
    \qquad (B \in \bB(\RR), \; n = 1, \ldots, N)
$$ (eq-deffx2)

This is of course the marginal distribution relative to the joint distribution of $\boldx$, as introduced in {ref}`ss-marg`. In this context it is called the **marginal distribution of $x_n$** in order to distinguish it from the (joint) distribution of $\boldx$. From a mathematical perspective, the marginal distribution of $x_n$ and the distribution of $x_n$ are identical.

If $P_1 = P_2 = \cdots = P_N$, then we say that $x_1, \ldots, x_N$ are **identically distributed**.

(sss-grv)=
#### Gaussian Random Vectors
Let $z$ be standard normal on $\RR$. A random variable $x$ is said to be **normally distributed** if $x = \mu + \sigma z$ for some $\sigma \geq 0$, and we write $\lL(x) = \nN(\mu, \sigma)$. Here we explicitly allow the possibility that $\sigma = 0$, in which case $\var[x] =0$ and $\PP\{x = \mu\} = 1$.

A random vector $\boldx$ in $\RR^N$ is called **multivariate normal** or **multivariate Gaussian** if there exists a $K \times 1$ standard normal random vector $\boldz$, an $N \times K$ matrix $\boldC$ and an $N \times 1$ vector $\boldmu$ such that

$$
\boldx = \boldmu + \boldC \boldz
$$ (eq-imn)

Once again we explicitly allow degenerate cases, since they occur frequently in statistics. If $\boldx$ is multivariate normal, then we write $\lL(\boldx ) = \nN(\boldmu, \boldSigma)$, where

```{math}
:enumerated: false

\boldmu := \EE \boldx 
    \quad \text{and} \quad
    \boldSigma := \var \boldx
```

In the case of {eq}`eq-imn` we have $\boldSigma = \boldC\boldC^\T$. See [fact %s](#fa-ltavc).

Be aware that $\lL(\boldx) = \nN(\boldmu, \boldSigma)$ does not imply that $\boldx$ has the multivariate normal density in {eq}`eq-mnormden`. The reason is that the distribution of $\boldx$ can fail to be absolutely continuous. For example, if $\boldC = \boldzero$ in {eq}`eq-imn`, then $\boldx$ puts all mass on a single point in $\RR^N$. Without absolute continuity, no density representation can exist.

Absolute continuity of the distribution of $\boldx$ coincides with the setting where $\boldSigma := \var \boldx$ is nonsingular. (Perhaps because of this, the scenario where a random vector has a density is sometimes called the "nonsingular case\" by econometricians.) For $\boldx$ defined in {eq}`eq-imn`, nonsingularity of $\boldSigma$ will be true if and only if $\boldC^\T$ has full column rank (see [ex. %s](#ex-fr)).

````{prf:property}
:label: fa-nipre

Let $\boldx$ be a random vector in $\RR^N$. The following statements are true:

(i) The vector $\boldx$ is multivariate normal if and only if $\bolda^\T \boldx$ is normally distributed in $\RR$ for every constant $N \times 1$ vector $\bolda$.

(ii) If $\lL(\boldx) = \nN(\boldmu, \boldSigma)$, then

```{math}
:enumerated: false

\lL(\boldA \boldx + \boldb) = \nN(\boldA \boldmu + \boldb, \boldA
                        \boldSigma \boldA^\T)
```

for all constant conformable $\boldA, \boldb$.
````


One corollary of [fact %s](#fa-nipre) is that if $\boldx = (x_1, \ldots, x_N)$ is multivariate normal, then the marginal distribution of $x_n$ is univariate normal. This can be seen by taking $\bolda = \bolde_n$ in part (i). In fact the rule in part (ii) implies that

$$
\lL(x_n)
    = \lL(\bolde_n^\T \boldx)
    = \nN(\bolde_n^\T \, \boldmu, \bolde_n^\T \boldSigma \bolde_n)
    = \nN(\mu_n, \sigma^2_n)
$$ (eq-mntm)

where $\sigma^2_n$ is the $n,n$th element of $\boldSigma$.

Here's a common exam question: Is it true that the joint distribution of $N$ univariate normal random variables is always multivariate normal? The answer is no. A counterexample is given in {ref}`ss-cop`. See also [fact %s](#fa-lcinorm) below.

(sss-efdv)=
#### Expectations from Distributions
Let $h \colon \RR^N \to \RR$ be any $\bB$-measurable function and let $P$ be a distribution on $\RR^N$. In {ref}`ss-igwd` we saw that one can assign meaning to the expression $\int h(s) P(\diff s)$ by treating $h$ as a random variable on $(\RR, \bB(\RR), P)$ and taking its expectation. We can do exactly the same thing in $\RR^N$, with $h$ now regarded as a random variable on $(\RR^N, \bB(\RR^N), P)$. Its expectation will be written as

$$
\EEP h :=: \int h(\bolds) P(\diff \bolds)
$$ (eq-eph)

Here's a multivariate version of [fact %s](#fa-hasp):

```{prf:property}
:label: fa-haspvec

Let $h \colon \RR^N \to \RR$ be $\bB$-measurable and let $P$ be a distribution on $\RR^N$. If $P$ is discrete, with PMF $\{p_j\}_{j \geq 1}$ and support $\{\bolds_j\}_{j \geq 1}$, then

$$
\int h(\bolds) P(\diff \bolds) = \sum_{j\geq 1} h(\bolds_j) p_j
$$ (eq-hasp2)
If $P$ is absolutely continuous with density $p$, then

$$
\int h(\bolds)P(\diff \bolds) 
        = \int h(\bolds) p(\bolds) \diff \bolds
$$ (eq-hasd2)
```


The right-hand side of {eq}`eq-hasd2` should be understood as

```{math}
:enumerated: false

\int_{-\infty}^\infty
        \cdots
        \int_{-\infty}^\infty
        h(s_1, \ldots, s_N) \,
        p(s_1,\ldots,s_N)  \,
        \diff s_1 \cdots \diff s_N
```

In {ref}`ss-expfd` we saw that expectations of functions of random variables can be computed from their distributions. The same is true for functions of random vectors. In particular, suppose that we want to compute $\EE h(\boldx)$, where $h \colon \RR^N \to \RR$ is $\bB$-measurable. We can do this either by evaluating the expectation of $h(\boldx)$ on the original probability space $(\Omega, \fF, \PP)$, or we can take the distribution $P$ of $\boldx$ and compute the expectation of $h$ on $(\RR^N, \bB(\RR^N), P)$. This is formalized in the following extension to [fact %s](#fa-efd).

````{prf:property}
:label: fa-efdvec

If $\lL(\boldx) = P$ and $h \colon \RR^N \to \RR$ is a $\bB$-measurable function such that $h(\boldx)$ is integrable, then

```{math}
:enumerated: false

\EE h(\boldx) = \int h(\bolds) P(\diff \bolds)
```

where $\int h(\bolds) P(\diff \bolds)$ is the expectation of $h$ on $(\RR^N, \bB(\RR^N), P)$.
````


As in the univariate case, this means that objects like moments are properties of the distribution. For example, let $\boldx$ be a random vector in $\RR^K$ with $\lL(\boldx) = P$. The variance–covariance matrix $\var [\boldx]$ of $\boldx$ has $i,j$th element $\EE [x_i x_j] - \EE[x_i]\EE[x_j]$. We can write this in terms of $P$. In particular, if

$$
\boldSigma_P = (\sigma_{ij}) 
    \quad \text{where} \quad
    \sigma_{ij} := \int (s_i s_j) P(\diff \bolds) 
        - \int s_i P(\diff \bolds) \cdot \int s_j P(\diff \bolds)
$$ (eq-vcd2)

then $\boldSigma_P = \var[\boldx]$.

(ss-indep)=
### Independence
A collection of $N$ random variables $x_1, \ldots, x_N$ is said to be **independent** if

$$
\PP \bigcap_{n=1}^N \{x_n \in B_n\}
    = \prod_{n=1}^N \PP\{x_n \in B_n\}
$$ (eq-bcap)

for any $B_1, \ldots, B_N$, where each $B_n$ is a Borel subset of $\RR$. In view of the definition of independence of events, the random variables $x_1, \ldots, x_N$ are independent precisely when sets of the form $\{x_1 \in B_1\}, \ldots, \{x_N \in B_N\}$ are independent events.

An infinite set of random variables $\{x_n\}_{n=1}^{\infty}$ is called independent if any finite subset of $\{x_n\}_{n=1}^{\infty}$ is independent in the sense of the preceding definition. We use the abbreviation **IID** for collections of random variables that are both independent and identically distributed. If we state that $\{x_n\}_{n=1}^{\infty}$ are **IID copies of $x$**, we mean that the random variables $\{x_n\}_{n=1}^{\infty}$ are IID with the same distribution as $x$.

There are a number of equivalent ways that we can rephrase the definition of independence. Since $\cap_{n=1}^N \{x_n \in B_n\} = \{ (x_1, \ldots, x_N) \in B_1 \times \cdots \times B_N \}$, a direct translation of {eq}`eq-bcap` using {eq}`eq-dmnu` and {eq}`eq-deffx2` gives

$$
P(B_1 \times \cdots \times B_N) 
    = \prod_{n=1}^N P_n(B_n)
$$ (eq-bcap2)

Here $P$ is the joint distribution of $\boldx = (x_1, \ldots, x_N)$ and $P_n$ is its $n$th marginal. Thus, the elements of a random vector are independent if and only if their joint distribution equals the product distribution formed from their marginals (see [fact %s](#fa-prodist)).

The left- and right-hand sides of {eq}`eq-bcap2` agree for all $B_1, \ldots, B_N$ if they agree on certain families of sets, like the rectangles discussed in {ref}`ss-muvd`. For this reason a necessary and sufficient condition for independence of $x_1, \ldots, x_N$ is

$$
F(s_1,\ldots,s_N) 
    = \prod_{n=1}^N F_n(s_n)
$$ (eq-pdind)

for all $(s_1, \ldots, s_N) \in \RR^N$, where $F$ is the CDF of $\boldx$ and $F_1, \ldots, F_N$ are the marginal CDFs . If the distribution of $\boldx$ is absolutely continuous we can also test independence via its density:

````{prf:property}
:label: fa-immd

If $\boldx = (x_1, \ldots, x_N)$ has joint density $p$ and marginals $p_1, \ldots, p_N$, then $x_1, \ldots, x_N$ are independent if and only if

```{math}
:enumerated: false

p (s_1,\ldots, s_N) = \prod_{n=1}^N p_n(s_n)
        \qquad \text{for all} \quad
        (s_1, \ldots, s_N) \in \RR^N
```
````


````{prf:example}
:label: eg-fmg

Let $\lL(\boldx) = \lL(x_1, \ldots, x_N) = N(\boldmu, \boldSigma)$. Suppose in addition that $\boldSigma$ is diagonal, with $n$th diagonal component $\sigma_n > 0$. Then $x_1, \ldots, x_N$ are independent. To see this observe that for any $\bolds = (s_1, \ldots, s_N) \in \RR^N$ we have

```{math}
:enumerated: false

\begin{aligned}
p(\bolds) 
        & = (2 \pi)^{-N/2} \determinant (\boldSigma)^{-1/2} 
        \exp \left\{ 
            - \frac{1}{2} (\bolds - \boldmu)^\T \boldSigma^{-1} (\bolds - \boldmu) 
        \right\}
        \\
        & = \frac{1}{(2 \pi)^{N/2} \prod_{n=1}^N \sigma_n}
        \exp \left\{ 
            - \frac{1}{2} \sum_{n=1}^N (s_n - \mu_n)^2 \sigma_n^{-2}
        \right\}
\end{aligned}
```

where the computation of the determinant and inverse of $\boldSigma$ used [facts %s](#fa-dtmat) and [%s](#fa-powdm). The last expression can be factored further, leading to

```{math}
:enumerated: false

p(\bolds) 
        =
        \prod_{n=1}^N \frac{1}{(2 \pi)^{1/2} \sigma_n}
            \exp \left\{ 
                \frac{-(s_n - \mu_n)^2}{2 \sigma_n^2}
            \right\}
         =
        \prod_{n=1}^N p_n(s_n)
```

where $p_n$ is the density of $\nN(\mu_n, \sigma^2_n)$. Independence now follows from [fact %s](#fa-immd).
````


The next fact is true in general but most easily verified in the density case, where you can apply {eq}`eq-hasd2` and [fact %s](#fa-immd).

````{prf:property}
:label: fa-imme

If $x_1, \ldots, x_N$ are independent and each $x_n$ is integrable, then

```{math}
:enumerated: false

\EE \left[ \prod_{n=1}^N x_n \right]
       = \prod_{n=1}^N \EE[ x_n ]
```
````


The definition of independence of random vectors is almost identical to the scalar case. Random vectors $\boldx_1, \ldots, \boldx_N$ in $\RR^K$ are called **independent** if

$$
\PP \bigcap_{n=1}^N \{\boldx_n \in B_n\}
    = \prod_{n=1}^N \PP\{\boldx_n \in B_n\}
$$ (eq-bcapv)

for any $B_1, \ldots, B_N$, where each $B_n$ is a Borel subset of $\RR^K$.

```{prf:property}
:label: fa-rviifi

If $\boldx_1, \ldots, \boldx_N$ are independent random vectors in $\RR^K$ and $f_1, \ldots, f_N$ are any $\bB$-measurable functions, then $f_1(\boldx_1), \ldots, f_N(\boldx_N)$ are also independent.
```


To see why, observe that $f_n(\boldx_n) \in B_n$ if and only if $\boldx_n \in f^{-1}(B_n)$, as follows from the definition of preimage of a set (see {ref}`ss-functions`). This leads to

```{math}
:enumerated: false

\bigcap_{n=1}^N \{f_n(\boldx_n) \in B_n\}
    = \bigcap_{n=1}^N \{\boldx_n \in f^{-1}(B_n)\}
```

Applying independence of $\boldx_1, \ldots, \boldx_N$ now gives

```{math}
:enumerated: false

\PP
    \bigcap_{n=1}^N \{f_n(\boldx_n) \in B_n\}
    = \prod_{n=1}^N \PP \{\boldx_n \in f^{-1}(B_n)\}
    = \prod_{n=1}^N \PP \{f_n(\boldx_n) \in B_n\}
```

This completes the proof of [fact %s](#fa-rviifi).

```{prf:property}
:label: fa-indzcov

If $\boldx$ and $\boldy$ are independent, then $\cov(\boldx,\boldy) = 0$.
```


Note that the converse is not true: one can construct examples of dependent random variables with zero covariance. However,

```{prf:property}
:label: fa-ucimin

If $\boldx$ is multivariate Gaussian and $\boldA$ and $\boldB$ are conformable constant matrices, then $\boldA \boldx$ and $\boldB \boldx$ are independent if and only if $\cov(\boldA \boldx, \boldB \boldx) = \boldzero$.
```


We already proved a simple version of this in [example %s](#eg-fmg). A full proof can be put together using the theory of characteristic functions (see, e.g., ch. 5 of {cite:author}`gutic` {cite:year}`gutic`). Below are two useful related results. The first one links orthogonal projection and multivariate Gaussians ([ex. %s](#ex-poci) requests a proof).

```{prf:property}
:label: fa-poci

Let $S$ be any linear subspace of $\RR^N$, let $\boldP := \proj S$ and let $\boldM$ be the residual projection. If $\lL(\boldz) = \nN(\boldzero, \sigma^2 \boldI)$ in $\RR^N$ for some $\sigma^2 > 0$, then $\boldP \boldz$ and $\boldM \boldz$ are independent.
```


````{prf:property}
:label: fa-lcinorm

If $w_1, \ldots, w_N$ are independent with $\lL(w_n) = \nN(\mu_n, \sigma_n^2)$ for all $n$, then

```{math}
:enumerated: false

\lL \left[ \alpha_0 + \sum_{n=1}^N \alpha_n w_n \right]
        = \nN \left( 
            \alpha_0 + \sum_{n=1}^N \alpha_n \mu_n,\;
            \sum_{n=1}^N \alpha_n^2 \sigma_n^2
            \right)
```
````


This last fact is important because sums of arbitrary normals are not always normal. What is required is a multivariate normal distribution. [Fact %s](#fa-lcinorm) works because independence pins down the joint distribution as multivariate normal. In particular, $\lL(w_1, \ldots, w_N) = \nN(\boldmu, \boldSigma)$ where $\bolde_n^\T \, \boldmu = \mu_n$ and $\boldSigma = \diag(\sigma_1^2, \ldots, \sigma_N^2)$. The expression for the variance in [fact %s](#fa-lcinorm) follows from [fact %s](#fa-vacov).

(ss-cop)=
### Copulas
A **copula** $C$ on $\RR^N$ is a multivariate CDF supported on the unit hypercube $[0, 1]^N$ with the property that all its marginals are uniform on $[0, 1]$. In other words, $C$ is a function of the form

$$
C(s_1, \ldots, s_N) = \PP \{u_1 \leq s_1, \ldots, u_N \leq s_N\}
$$ (eq-icif)

where $0 \leq s_n \leq 1$ and $\lL(u_n) = U[0, 1]$ for all $n$. There are many possible copulas because, while each $u_n$ has its marginal distribution pinned down, there are infinitely many ways to specify the joint distribution.

```{prf:example}

The function $C(s_1, s_2) = s_1 s_2$ on $[0, 1]^2$ is called the **independence copula**. The marginal distributions are $C(s_1, 1) = s_1$ and $C(1, s_2) = s_2$ as required. (These are CDFs for the $U[0, 1]$ distribution.)
```


````{prf:example}

The **Gumbel copulas** are the class of functions on $[0, 1]^2$ defined by

```{math}
:enumerated: false

C(s_1, s_2) = 
        \exp 
        \left\{
            - \left[ 
                (-\ln s_1)^\theta + (-\ln s_2)^\theta
            \right]^{1/\theta}
        \right\}
        \qquad \quad (\theta \geq 1)
```

The **Clayton copulas** are given by

```{math}
:enumerated: false

C(s_1, s_2) = 
        \left\{
            \max \left[
                s_1^{-\theta} +  s_2^{-\theta} -1, \, 0
                \right]
        \right\}^{-1/\theta}
        \qquad \quad (\theta \geq -1, \, \theta \neq 0)
```

Both of these belong to a general class called the **Archimedean copulas**.
````


One reason that copulas are useful is that we can take univariate CDFs $F_1, \ldots, F_N$ and a copula $C$ and create a multivariate CDF on $\RR^N$ via

$$
F(s_1, \ldots, s_N) = C(F_1(s_1), \ldots, F_N(s_N))
    \qquad \quad (s_n \in \RR, \, n=1, \ldots, N)
$$ (eq-ffcop)

The benefit of creating a multivariate CDF using this two step procedure is that in specifying a multivariate distribution we can separate out specification of the marginals, which give probabilities for isolated coordinates, and specification of the joint distribution, which defines the interactions.

```{prf:example}

{cite}`bonhomme2009assessing` use copulas to model one component of earnings dynamics in a study based on three-year panels from the French Labor Force Survey. The cross sections are relatively large (around 30,000), allowing for flexible modeling of the marginal distributions via a mixture of normals. However, the time series dimension is short, so a one-parameter family of copulas is used to bind the marginals across time in a parsimonious way.
```


How rich is the set of multivariate CDFs we can create via {eq}`eq-ffcop`? It is universal. This is one part of a well-known theorem due to {cite}`sklar1959fonctions`:

```{prf:theorem}
:label: t-sklar

If $F$ is any CDF on $\RR^N$ with marginals $F_1, \ldots, F_N$, then there exists a copula $C$ such that {eq}`eq-ffcop` holds. If each $F_n$ is continuous, then this representation is unique.
```


Copulas can be used to clarify a point raised in {ref}`sss-grv`: normal marginals are not sufficient for the multivariate normal property. If $F_1, \ldots, F_N$ are univariate normal, then $C(F_1(s_1), \ldots, F_N(s_N))$ will equal the multivariate normal CDF for one choice of copula, called the Gaussian copula. Other choices lead to different distributions. {numref}`f-copula` gives an illustration where the top two distributions are multivariate normal and the bottom two are created by pairing the Clayton copula with univariate normal marginals. The distributions are represented as densities.

```{figure} figures/copula.svg
:name: f-copula

Bivariate Gaussian (top) and non-Gaussian (bottom)
```

A related point is that while linear combinations of normal random variables are often thought to be normal, this isn't true unless the original random variables are multivariate normal (cf., [fact %s](#fa-nipre)) or independent (which is a special case of multivariate normal). As an example, [figure %s](#f-sum_not_normal) shows a histogram of draws from $w = x-y$ when $(x,y)$ is distributed according to the density in the lower right panel of [figure %s](#f-copula). For comparison, a normal density with the same mean and standard deviation is plotted.

(ss-pnd)=
### Properties of Named Distributions
Here are some useful facts relating independence and certain common distributions.

```{prf:property}
:label: fa-schis

If $x_1,\ldots,x_N$ are independent and $\lL(x_n) = \chi^2(k_n)$, then $\lL(\sum_n x_n) = \chi^2(\sum_n k_n)$.
```


```{figure} figures/sum_not_normal.svg
:name: f-sum_not_normal

Normality is not necessarily preserved under linear transformations
```

````{prf:property}
:label: fa-astudt

If $z$ and $x$ are independent with $\lL(z) = \nN(0,1)$ and $\lL(x) = \chi^2(k)$, then

```{math}
:enumerated: false

z \sqrt{\frac{k}{x} } 
        \;\;
        \text{ is } t \text{ distributed with } k \text{ degrees of freedom}
```
````


````{prf:property}
:label: fa-rcif

If $x_1$ and $x_2$ are independent with $\lL(x_i) = \chi^2(k_i)$ for $i=1,2$, then

```{math}
:enumerated: false

\frac{x_1 / k_1}{x_2 / k_2} 
        \;\;
        \text{ is } F \text{ distributed with parameters } k_1, k_2
```
````


Let's note some facts connecting chi-squared and normal random variables. First, if $\lL(x) = \nN(0, 1)$, then $\lL(x^2) = \chi^2(1)$. {prf:ref}`ex-chis1` asks you to prove this. Second, if a collection of standard normals are independent, then so are their squares ([fact %s](#fa-rviifi)). Combining this result with [fact %s](#fa-schis) gives

```{prf:property}
:label: fa-ssnac

If $\lL(z_1, \ldots, z_N) = \nN(\boldzero, \boldI)$, then $\lL(\sum_{n=1}^N z_n^2) = \chi^2(N)$.
```


The next result generalizes [fact %s](#fa-ssnac) and crops up frequently in statistics.

````{prf:property}
:label: fa-mcfi

If $\lL(\boldz) = \nN(\boldzero, \boldI)$ and $\boldA$ is symmetric and idempotent, then

```{math}
:enumerated: false

\lL \left( \boldz^\T \boldA \boldz \right) = \chi^2(K)  
        \quad \text{where} \quad
        K := \trace \boldA
```
````


We can obtain [fact %s](#fa-mcfi) from [fact %s](#fa-ssnac) as follows. Since $\boldA$ is symmetric, we use the spectral decomposition theorem to write it as $\boldA = \boldQ \boldLambda \boldQ^\T$, where $\boldQ$ is an orthogonal matrix and $\boldLambda = \diag(\lambda_1, \ldots, \lambda_N)$ with $\lambda_n :=$ the $n$th eigenvalue of $\boldA$. Setting $\boldy := \boldQ^\T \boldz$, we have

```{math}
:enumerated: false

\boldz^\T \boldA \boldz
    = \boldz^\T \boldQ \boldLambda \boldQ^\T \boldz
    = \boldy^\T \boldLambda \boldy
    = \sum_{n=1}^N \lambda_n y_n^2
```

By [exercise %s](#ex-eigid), idempotence of $\boldA$ implies all eigenvalues lie in $\{0, 1\}$. By [fact %s](#fa-eiprop), the number of nonzero eigenvalues is equal to $\trace \boldA$. It follows that $\boldz^\T \boldA \boldz$ is the sum of $y^2_n$ over $K := \trace \boldA$ nonzero terms. By [fact %s](#fa-ssnac), a random variable with this property is $\chi^2(K)$ whenever $\boldy$ is standard normal. Combining $\lL(\boldz) = \nN(\boldzero, \boldI)$ and [fact %s](#fa-nipre), we have $\lL(\boldy) = \nN(\boldzero, \boldQ^\T\boldQ)$. But $\boldQ$ is an orthogonal matrix, so $\boldQ^\T\boldQ = \boldI$.

## Conditioning and Expectation

Conditional expectation is one of the most important concepts in both economic theory and econometrics. This section gives a construction of expectation based around projection, which, in my view, is more intuitive and rewarding than both the elementary treatments using conditional distributions and the advanced treatments using Radon–Nykodym derivatives. The construction begins by formalizing the concepts of information and conditioning. We then build on our knowledge of orthogonal projection from {ref}`s-pd` to frame conditional expectation as optimal prediction given limited information.

### Conditional Distributions

Let's briefly discuss conditional distributions, focusing on the density case. Let $x_1$ and $x_2$ be random variables. The **conditional density** of $x_2$ given $x_1 = s_1$ is defined as

$$
p(s_2 \given s_1) 
        := \frac{p(s_1, s_2)}{p(s_2)}
$$

whenever the denominator is nonzero. Here $p$ stands in for either joint, marginal, or conditional density, with the type determined by the argument. (For example, $p(s_i)$ is the marginal density of $x_i$ evaluated at $s_i$.)

The law of total probability extends to the density case as follows: If $(x_1, x_2)$ is a random vector on $\RR^2$, then

$$
p(s_2)
    = \int_{-\infty}^{\infty} p(s_2 \given s_1) p(s_1) \diff s_1
    \qquad \qquad 
    (s_2 \in \RR)
$$ (eq-dltp)

To see this, fix $s_2 \in \RR$ and integrate the joint density to get the marginal, giving

```{math}
:enumerated: false

p(s_2) = \int_{-\infty}^{\infty} p(s_1, s_2) \diff s_1
```

Combining this result with $p(s_2 \given s_1) = p(s_1, s_2)/p(s_1)$ yields {eq}`eq-dltp`. Bayes' law also extends to the density case, giving

$$
p(s_2 \given s_1) = \frac{p(s_1 \given s_2) p(s_2)}{p(s_1)}
$$ (eq-brdc)

The same ideas extend to additional arguments in a natural way. For example, the conditional density of $x_{k+1},\ldots, x_N$ given $x_1 = s_1,\ldots, x_k = s_k$ is defined by

$$
p(s_{k+1},\ldots,s_N \given s_1,\ldots,s_k) 
        = \frac{p(s_1,\ldots,s_N)}{p(s_1,\ldots,s_k)}
$$ (eq-condden)

Rearranging this expression, we obtain a useful decomposition of the joint density:

$$
p(s_1,\ldots,s_N) 
    = p(s_{k+1},\ldots,s_N \given s_1,\ldots,s_k) p(s_1,\ldots,s_k)
$$ (eq-decompjd)

(ss-tsl2)=
### The Space $L_2$
Suppose we want to predict random variable $y$ using another variable $x$. In this case we should choose $x$ such that $x$ and $y$ are expected to be close under most realizations of uncertainty. But what does "expected to be close\" mean?

Perhaps the most common measure of distance between random variables is **mean squared error** (MSE), which is defined in this case as $\EE [ (x - y)^2 ]$. We met this concept previously in the linear prediction problem of {ref}`ss-blp`. A common variation on MSE is **root mean squared error**, which, as the name suggests, is defined as

$$
\| x - y \| := \sqrt{ \EE [ (x - y)^2 ] }
$$ (eq-ndrv)

On the right-hand side of {eq}`eq-ndrv` we have a measure of distance between random variables. On the left-hand side we have used the same symbol as we used for Euclidean norm distance between two vectors. This is both deliberate and conventional, as there are many parallels between ordinary vector space with the Euclidean norm and the set of random variables combined with the "norm\" defined in {eq}`eq-ndrv`. Let's now clarify and formalize these ideas.

The first geometric concept we defined for vectors was inner product. Analogously, we define here the **inner product between two random variables** $x$ and $y$ as

$$
\inner{x, y} := \EE[xy]
$$ (eq-ipl2)

Of course we want this definition to make sense, and the Cauchy–Schwarz inequality for random variables tells us that $\EE[xy]$ will be finite and well-defined whenever $x$ and $y$ both have finite second moments. Hence we restrict attention to this set, which is commonly denoted $L_2$. That is,

```{math}
:enumerated: false

L_2 
    := \{ \text{ all random variables } x \text{ on } (\Omega, \fF, \PP)
    \text{ with } \EE [x^2] < \infty \}
```

Here's a list of some of the ways that $L_2$ inner product is similar to the Euclidean inner product on $\RR^N$.

```{prf:property}
:label: fa-innpp2

For any $\alpha, \beta \in \RR$ and any $x, y, z \in L_2$, the following statements are true:

(i) $\inner{x, y} = \inner{y, x}$.

(ii) $\inner{\alpha x, \beta y} =  \alpha \beta \inner{x, y}$.

(iii) $\inner{x, \alpha y + \beta z} =  \alpha \inner{x, y} + \beta \inner{x, z}$.
```


If you compare [fact %s](#fa-innpp2) with [fact %s](#fa-innpp) you'll see that the correspondence is essentially one-to-one. The properties in [fact %s](#fa-innpp2) follow from the definition of the inner product and linearity of $\EE$ (see [fact %s](#fa-exprop)).

Analogous to the Euclidean vector norm, we define the **$L_2$ norm** by

$$
\| x \| 
    := \sqrt{ \inner{x, x} } 
    := \sqrt{ \EE[x^2] } 
    \qquad (x \in L_2)
$$ (eq-rvnorm)

From this norm we get a notion of distance $\|x - y\|$ between random variables that agrees with the notion of root MSE in {eq}`eq-ndrv`. As with inner product, the properties of the $L_2$ norm and vector norm are largely identical. (Compare the following with [fact %s](#fa-profno).)

```{prf:property}
:label: fa-l2se

For any $\alpha \in \RR$ and any $x, y \in L_2$, the following statements are true:

(i) $\| x \| \geq 0$ and $\| x \| = 0$ if and only if $x = 0$

(ii) $\| \alpha x \| = |\alpha| \| x \|$

(iii) $\| x + y \| \leq  \| x \| + \| y \|$

(iv) $| \! \inner{x, y}\!  | \leq  \| x \| \| y \|$
```


Point (i) requires a clarification, which is given below. Property (ii) is immediate from {eq}`eq-rvnorm` and the linearity of $\EE$. Property (iii) is called the **triangle inequality**, as in the vector case. Property (iv) is just the **Cauchy–Schwarz inequality**, restated in new notation. As in the vector case, the triangle inequality can be proved from the Cauchy–Schwarz inequality. For a proof, see [exercise %s](#ex-csrv).

Returning to (i), it isn't true that $\| x \| = 0$ implies $x(\omega) = 0$ for all $\omega \in \Omega$. What we can say is that if $\| x \| = 0$ then $\PP\{x =0\} =1$ (see, e.g., {cite:author}`williams1991probability` {cite:year}`williams1991probability`, §5.2). If we agree not to distinguish between random variables that only differ with zero probability, then (i) is valid. This is the standard convention when dealing with $L_2$.

#### Linear Subspaces in $L_2$

From (ii) and (iii) of [fact %s](#fa-l2se), we see that if $\alpha, \beta \in \RR$ and $x, y \in L_2$, then $\alpha x + \beta y$ also has finite norm, and is therefore in $L_2$. More generally, any **linear combination**

```{math}
:enumerated: false

\alpha_1 x_1 + \cdots \alpha_K x_K,
    \qquad \alpha_k \in \RR, \; x_k \in L_2
```

is again in $L_2$. As in the vector case (see {ref}`ss-lcs`), when $X$ is a subset of $L_2$, the set of finite linear combinations that can be formed from elements of $X$ is called the **span** of $X$, and denoted by $\Span X$.

```{prf:example}
:label: eg-llasp

If $x \in L_2$ and $\1 := \1_\Omega$ is the constant random variable always equal to $1$, then $\Span\{\1, x\}$ is the set of random variables

$$
\alpha + \beta x 
         := \alpha \1 + \beta x
         \quad \text{for scalars} \quad
         \alpha, \beta
$$ (eq-llasp)
This is the set $\lL$ introduced when we discussed best linear predictors.
```


A subset $S$ of $L_2$ is called a **linear subspace** of $L_2$ if it is closed under addition and scalar multiplication. That is, for each $x, y \in S$ and $\alpha, \beta \in \RR$, we have $\alpha x + \beta y \in S$.

```{prf:example}

The span of any set of elements of $L_2$ is a linear subspace in $L_2$. (Compare with [example %s](#eg-silsub).)
```


````{prf:example}

The set $Z := \setntn{x \in L_2}{\EE x = 0}$ is a linear subspace of $L_2$ because

```{math}
:enumerated: false

x, y \in Z
        \text{ and }
        \alpha, \beta \in \RR
        \quad \implies \quad
        \EE[\alpha x + \beta y] = \alpha \EE[x] + \beta \EE[y] = 0
```
````


As in $\RR^N$, an **orthonormal basis** of a linear subspace $S$ of $L_2$ is a set $\{u_1, \ldots, u_K\} \subset S$ with the property that

```{math}
:enumerated: false

\inner{u_j, u_k} = \1\{j = k\}
    \quad \text{and} \quad
    \Span\{u_1, \ldots, u_K\} = S
```

````{prf:example}
:label: eg-orbrv

Let $x \in L_2$ so that $S := \Span\{\1, x\}$ is as in {eq}`eq-llasp`. If we define

```{math}
:enumerated: false

u_1 := \1 
        \quad \text{and} \quad
        u_2 := \frac{x - \mu}{\sigma_x}
```

Then

```{math}
:enumerated: false

\inner{ u_1, u_2 } 
        = \EE[u_1 u_2] 
        = \EE \left[ \frac{x - \mu}{\sigma_x} \right]
        = 0
```

Clearly, $\|u_1\| = \|u_2\| = 1$, so this pair is orthonormal. It's also straightforward to show that $\Span\{u_1, u_2\} = \Span\{\1, x\}$, so $\{u_1, u_2\}$ is an orthonormal basis for $S$.
````


(ss-projl2)=
### Projections in $L_2$
Given that inner products and norms in $L_2$ share so many properties with their Euclidean cousins, it's not surprising that many results about $\RR^N$ carry over to $L_2$. Often we can just lift the proofs across with barely a change to notation. We know a great deal about vectors and Euclidean space, so this process yields many insights about random variables. Let's illustrate this in the context of orthogonal projection.

As in the Euclidean case, if $\inner{x,y}=0$, then we say that $x$ and $y$ are **orthogonal**, and write **$x \perp y$**. This terminology is used frequently in econometrics. Its importance stems partly from the next fact, the proof of which is a short exercise.

```{prf:property}

If $x, y \in L_2$ and $\EE x=0$ or $\EE y=0$ then $x \perp y \iff  \cov[x,y] = 0$.
```


Now consider the following problem. Given $y \in L_2$ and linear subspace $S \subset L_2$, we seek the closest element $\hat y$ of $S$ to $y$. Closeness is in terms of $L_2$ norm, so $\hat y$ is the minimizer of $\| y - z \|$ over all $z \in S$. That is, we seek

$$
\hat y 
    = \argmin_{z \in S} \|y - z\|
    = \argmin_{z \in S} \sqrt{\EE [ (y - z)^2 ]}
$$ (eq-mproj2)

The next theorem mimics [theorem %s](#t-opt). A nice treatment that covers the current setting can be found in chapter 2 of {cite}`cheney2013analysis`.

```{prf:theorem}
:label: t-opt4

Let $y \in L_2$ and let $S$ be any nonempty closed linear subspace of $L_2$. The following statements are true:

(i) The optimization problem {eq}`eq-mproj2` has exactly one solution.

(ii) $\hat y \in L_2$ is the unique solution[^3] to {eq}`eq-mproj2` if and only if $\hat y \in S$ and $y - \hat y \perp S$.
```


The statement that $S$ is closed means that $\{x_n \} \subset S$ and $x \in L_2$ with $\| x_n - x \| \to 0$ implies $x \in S$. This condition is true for all the linear subspaces we want to work with.

Analogous with the case of $\RR^N$, the random variable $\hat y$ in [theorem %s](#t-opt4) is called the **orthogonal projection of $y$ onto $S$**. Holding $S$ fixed, we can think of the operation

```{math}
:enumerated: false

y \; \mapsto \text{ the orthogonal projection of } y \text{ onto } S
```

as a function from $L_2$ to $L_2$. The function will be denoted by $\boldP$, so that, for each $y \in L_2$, the symbol $\boldP y$ represents the image of $y$ under $\boldP$, which is the orthogonal projection $\hat y$. As before, $\boldP$ is called the **orthogonal projection onto $S$**, and we write $\boldP = \proj S$. Since $\boldP y$ solves {eq}`eq-mproj2`, we can and do interpret $\boldP y$ as the *best predictor of $y$ from within the collection of random variables contained in $S$*.

```{prf:property}
:label: fa-l2opt

Let $S$ be any linear subspace of $L_2$, and let $\boldP = \proj S$. The following statements are true:

(i) $\boldP$ is a linear function.

Moreover, for any $y \in L_2$, we have

(ii) $\boldP y \in S$,

(iii) $y - \boldP y \perp S$,

(iv) $\| y \|^2 = \| \boldP y \|^2 + \| y - \boldP y \|^2$,

(v) $\| \boldP y \| \leq \| y \|$, and

(vi) $\boldP y = y$ if and only if $y \in S$.
```


Compare with [theorem %s](#t-opt2). In (i), the statement that $\boldP$ is linear means that $\boldP(\alpha x + \beta y) = \alpha \boldP x + \beta \boldP y$ for all $x, y \in L_2$ and $\alpha,\beta \in \RR$.

The next two results are $L_2$ versions of [fact %s](#fa-subsub) and [fact %s](#fa-projon) respectively. The proofs are essentially the same.

```{prf:property}
:label: fa-subsubl2

Let $S_i$ be a linear subspace of $L_2$ for $i=1,2$ and let $\boldP_i = \proj S_i$. If $S_1 \subset S_2$, then $\boldP_1 \boldP_2 y = \boldP_1 y$ for all $y \in L_2$.
```


```{prf:property}
:label: fa-l2opt2

If $\{u_1, \ldots, u_K\}$ is an orthonormal basis of $S$, then, for all $y \in L_2$,

$$
\boldP y = \sum_{k=1}^K \inner{y, u_k} \, u_k
$$ (eq-projon2)
```


````{prf:example}
:label: eg-projam

The mean of a random variable $x$ can be thought of as the "best predictor of $x$ within the set of constants.\" To see this from a new perspective, let $S := \Span\{\1\}$, where $\1 := \1_\Omega$, and let $\boldP := \proj S$. The object $\boldP x$ is precisely the best predictor of $x$ within the class of constant random variables. Not surprisingly, $\boldP x = \mu \1$, where $\mu := \EE x$. The easiest way to check this is to observe that $\{\1\}$ is an orthonormal set spanning $S$, and hence, by {eq}`eq-projon2`,

```{math}
:enumerated: false

\boldP x 
        = \inner{x, \1} \, \1
        = \EE[x \1] \1 = \EE[x] \1 = \mu \, \1
```

You can also check the claim that $\mu \1$ is the projection of $x$ onto $S$ by verifying the conditions in (ii) of [theorem %s](#t-opt4).
````


````{prf:example}
:label: eg-rvblp

Let's revisit the problem of best linear predictors from {ref}`ss-blp`, establishing the same result in a different way. Fix $x, y \in L_2$ and consider projecting $y$ onto $S := \Span\{\1, x\}$. The set $S$ is as given in {eq}`eq-llasp`. The problem of projecting $y$ onto $S$ is equivalent to the best linear prediction problem in {eq}`eq-llsrf`. To implement the projection recall from [example %s](#eg-orbrv) that

```{math}
:enumerated: false

u_1 := \1 
        \quad \text{and} \quad
        u_2 := \frac{x - \mu}{\sigma_x}
```

form an orthonormal basis for $S$. Hence, letting $\boldP = \proj S$ and applying {eq}`eq-projon2`,

```{math}
:enumerated: false

\boldP y 
         = \la y, u_1 \ra u_1 + \la y, u_2 \ra u_2
         = \EE[ y ] + 
            \frac{\cov[x, y]}{\var[x]}(x - \EE[x])
```

Alternatively, we can write

```{math}
:enumerated: false

\boldP y 
        = \alpha^* + \beta^* x
        \quad \text{where} \quad
        \beta^* := \frac{\cov[x, y]}{\var[x]}
        \quad \text{and} \quad
        \alpha^* := \EE[y] - \beta^* \EE[x]
```

These are the same solutions we obtained using calculus in {ref}`ss-blp`.
````


(sss-pr)=
#### Population Regression
Let's consider an extension of the best linear prediction problem of [example %s](#eg-rvblp) to a more general setting, where the information for predicting $y$ is a random vector $\boldx$ in $\RR^K$. As before, prediction will be linear, with the aim of minimizing squared error. Hence the object we seek is the $L_2$ orthogonal projection of $y$ onto the linear subspace

```{math}
:enumerated: false

\Span\{\boldx\} 
    := \text{ all random variables of the form
        } \boldx^\T \boldb \text{ for some } \boldb \in \RR^K
```

To make sure that $\Span\{\boldx\} \subset L_2$, we will assume that $\EE[ \boldx^\T \boldx] < \infty$.

```{prf:property}
:label: fa-mpr

If $\EE[\boldx \boldx^\T]$ is positive definite, then the projection $\boldP y$ of any $y \in L_2$ onto $\Span\{\boldx\}$ is given by

$$
\hat y = \boldx^\T \boldb^*      
      \quad \text{where} \quad
      \boldb^* := \EE[\boldx \boldx^\T]^{-1} \EE[\boldx y]
$$ (eq-mpr)
```


{prf:ref}`ex-mpr` asks you to prove this. Positive definiteness of $\EE[\boldx \boldx^\T]$ ensures invertibility, without which the vector $\boldb^*$ is not uniquely defined (see [ex. %s](#ex-lsiosol)). Note that, by the definition of orthogonal projections, $\boldb^*$ necessarily satisfies

```{math}
:enumerated: false

\boldb^* = \argmin_{\bolda \in \RR^K} \EE[ (y - \boldx^\T \bolda)^2 ]
```

The linear prediction problem considered here is also called **population linear regression**. The "population\" qualifier is because we are using the true joint distribution of $(\boldx, y)$ when we compute expectations. Population regression has a sample counterpart called multivariate linear regression, based on observations of $(\boldx, y)$. We will discuss it at length in [chapter 11](#c-reg).

(ss-m)=
### Measurability
While the use of projections to obtain best linear predictors in [example %s](#eg-rvblp) and [fact %s](#fa-mpr) is elegant and useful, we don't always want to constrain ourselves to linear predictions. To drop the linearity requirement, we change the linear subspaces used for projection from the set of linear functions of $\boldx$ to the set of arbitrary functions of $\boldx$. As we now show, the best predictor that results is exactly the conditional expectation with respect to $\boldx$.

The key to this process is defining the subspace to project onto. In the language of probability, the subspace of arbitrary real-valued functions of $\boldx$ is called the $\boldx$-measurable functions. Let's take these ideas one step at a time, starting with the notion of measurability with respect to information sets.

Let $\gG := \{x_1, \ldots, x_D\}$ be any set of random variables and let $z$ be any other random variable. The variable $z$ is said to be **$\gG$-measurable** if there exists a $\bB$-measurable function $g \colon \RR^D \to \RR$ such that

$$
z = g(x_1, \ldots, x_D)
$$ (eq-gmeas)

As usual, this equality between random variables should be interpreted pointwise, so the more explicit version is

$$
z(\omega) 
    = g(x_1(\omega), \ldots, x_D(\omega)) \qquad \text{for all } \omega \in \Omega
$$ (eq-gmeas2)

In this context, $\gG$ is sometimes referred to as the **information set**. When convenient we'll also write $\boldx = (x_1, \ldots, x_D)$ and say that $z$ is $\boldx$-measurable. Similar terminology will be used for scalars and matrices. For example, if $\boldX$ is a random matrix, then $\boldX$-measurability means $\gG$-measurability when $\gG$ lists all elements of $\boldX$.

The intuitive meaning of $\gG$-measurability of $z$ is that $z$ is completely determined by the elements in $\gG$. Imagine it like this: When uncertainty is realized, some $\omega$ is selected from $\Omega$. We don't get to view $\omega$ itself, but we do get to observe the realized values $x_1(\omega), \ldots, x_D(\omega)$ on the right-hand side of {eq}`eq-gmeas2`. If $z$ is $\gG$-measurable, we can now calculate the realized value $z(\omega)$ of $z$, even without knowing $\omega$.

```{prf:example}

Let $x, y$ and $z$ be random variables and let $\alpha$ and $\beta$ be scalars. If $z = \alpha x + \beta y$, then $z$ is $\{x, y\}$-measurable (take $g(s, t) := \alpha s + \beta t$).
```


```{prf:example}

If $x_1, \ldots, x_N$ are random variables and $\gG := \{x_1, \ldots, x_N\}$, then the sample mean $\bar x_N := \frac{1}{N} \sum_{n=1}^N x_n$ is $\gG$-measurable.
```


```{prf:example}

Let $\boldx$ and $y$ be independent and nondegenerate. Then $y$ is not $\boldx$-measurable, for if it were, then we would have $y = g(\boldx)$ for some function $g$, contradicting independence of $\boldx$ and $y$.[^4]
```


```{prf:example}

Let $y$ be a noisy signal of $x$, in the sense that $y = x + u$ where $x$ and $u$ are independent and $u$ is nondegenerate. Then $y$ is not $x$-measurable. Intuitively, even if we know the value of $x$, the realization of $y$ cannot be computed until we know the value of the disturbance term $u$. More formally, if $y$ is $x$-measurable, then $y = g(x)$ for some $g$, in which case $u = g(x) - x$. But then $u$ is $x$-measurable, contradicting independence.
```


```{prf:example}
:label: eg-cmag

Let $y = \alpha$, where $\alpha$ is a constant. This degenerate random variable is $\gG$-measurable for any information set $\gG$, because $y$ is already deterministic. For example, if $\gG = \{x_1, \ldots, x_p\}$, then we can take $y = g(x_1, \ldots, x_p) = \alpha + \sum_{i=1}^p 0 x_i$.
```


If $x$ and $y$ are known given the information in $\gG$, then a third random variable that depends only on $x$ and $y$ is likewise known given $\gG$. Hence $\gG$-measurability is preserved under the taking of sums, products, and so on. In particular,

```{prf:property}
:label: fa-lcapom

Let $\alpha, \beta$ be any scalars, and let $x$ and $y$ be random variables. If $x$ and $y$ are both $\gG$-measurable, then $u := xy$ and $v := \alpha x + \beta y$ are also $\gG$-measurable.
```


Suppose now that $\gG \subset L_2$, which just means that all the random variables in $\gG$ have finite second moment. Consider the set

```{math}
:enumerated: false

L_2(\gG) 
    := \{ \text{all } \gG \text{-measurable random variables in } L_2 \}
```

In view of [fact %s](#fa-lcapom), we have the following result:

```{prf:property}

For any $\gG \subset L_2$, the set $L_2(\gG)$ is a linear subspace of $L_2$.
```


This furnishes us with a subspace to project onto, allowing us to define conditional expectations. But before we do so let's cover ordering of information sets by set inclusion.

```{prf:property}
:label: fa-mwrtgh

If $\gG \subset \hH$ and $z$ is $\gG$-measurable, then $z$ is $\hH$-measurable.
```


This is intuitive: If $z$ is known once the variables in $\gG$ are known, then it is certainly known when the extra information provided by $\hH$ is available.

````{prf:example}

Let $x_1$, $x_2$ and $y$ be random variables and let

```{math}
:enumerated: false

\gG := \{x_1\} \subset \{x_1, x_2\} =: \hH
```

If $y$ is $\gG$-measurable, then $y = g(x_1)$ for some $\bB$-measurable $g$. But then $y$ will also be $\hH$-measurable. For example, we can write $y = h(x_1,x_2)$ where $h(x_1,x_2) = g(x_1) + 0 x_2$.
````


The next result is immediate from [fact %s](#fa-mwrtgh).

```{prf:property}

If $\gG \subset \hH$, then $L_2(\gG) \subset L_2(\hH)$.
```


Incidentally, most advanced probability texts give a definition of $\gG$-measurability framed in terms of $\sigma$-algebras. This is necessary if you wish to work through all related proofs. But the definition given here is equivalent and intuitive. See theorem 4.4 of {cite}`cinlar` if you wish to learn more.

(ss-ce)=
### Conditional Expectation
Now it's time to define conditional expectations. Let $\gG \subset L_2$ and $y$ be some random variable in $L_2$. The **conditional expectation** of $y$ given $\gG$ is written as $\EE[ y \given \gG]$ or $\EE^{\gG} [y]$, and defined as the closest $\gG$-measurable random variable to $y$. More formally,

$$
\EE [y \given \gG] := \argmin_{z \in L_2(\gG)} \| y - z \|
$$ (eq-defcd)

This definition makes sense. Our understanding of the conditional expectation $\EE[y \given \gG]$ is that it is the best predictor of $y$ given the information contained in $\gG$. The definition in {eq}`eq-defcd` says the same thing. It simultaneously restricts $\EE [y \given \gG]$ to be $\gG$-measurable, so we can actually compute it once the variables in $\gG$ are realized, and selects $\EE [y \given \gG]$ as the closest such variable to $y$ in terms of RMSE.

But does the minimizer generally exist? And is it unique? The answer to both questions is yes. Comparing {eq}`eq-defcd` and {eq}`eq-mproj2`, we see that

```{math}
:enumerated: false

\EE[ y \given \gG] = \boldP y 
    \quad \text{when }  
    \boldP := \proj L_2(\gG)
```

By [theorem %s](#t-opt4), the projection exists and is unique. {prf:ref}`t-opt4` also characterizes $\EE[y \given \gG]$ as the unique point $\hat y$ in $L_2$ such that

```{math}
:enumerated: false

\hat y \in L_2(\gG)
    \quad \text{and} \quad
    y - \hat y  \perp z
    \; \text{ for all } z \in L_2(\gG)
```

Rewriting these conditions in a slightly different way, we can give an alternative (and equivalent) definition of conditional expectation: $\hat y  \in L_2$ is the **conditional expectation** of $y$ given $\gG$ if

(i) $\hat y$ is $\gG$-measurable and

(ii) $\EE[\hat y  \, z] = \EE[yz]$ for all $\gG$-measurable $z \in L_2$.

When convenient we'll also use symbols like $\EE[ y \given x_1, \ldots, x_D]$ or $\EE[y \given \boldx]$. These are the same as $\EE[y \given \gG]$ when $\gG$ is defined as the information set containing the variables we condition on.

````{prf:example}

If $x$ and $u$ are independent, $\EE u=0$ and $y = x + u$, then $\EE[y \given x] = x$. To prove this we need to show that $x$ satisfies (i)–(ii) above. Clearly, $x$ is $x$-measurable. For (ii) we need to show that $\EE[x \, z] = \EE[y\, z]$ for all $x$-measurable $z$. This translates to the claim that

```{math}
:enumerated: false

\EE[ x g(x)] = \EE[ (x + u) g(x)]
```

for any $\bB$-measurable $g$, which is true from independence and $\EE u = 0$.
````


```{prf:property}
:label: fa-ceisf

Given $\boldx \in \RR^D$ and $y$ in $L_2$, there exists a $\bB$-measurable function $f^* \colon \RR^D \to \RR$ such that $\EE[y \given \boldx] = f^*(\boldx)$.
```


This is obvious because $\EE[y \given \boldx]$ is by definition $\boldx$-measurable. A prediction of $y$ given $\boldx$ has to be a deterministic function of $\boldx$. The particular function $f^*$ satisfying $f^*(\boldx) = \EE[y \given \boldx]$ is called the **regression function** of $y$ given $\boldx$. It is, by construction, the best predictor of $y$ given $\boldx$ in terms of mean squared error.

The next example shows that when $x$ and $y$ are linked by a conditional density, then our definition of conditional expectation reduces to the one seen in elementary probability texts. The proof of the claim in the example is the topic of [exercise %s](#ex-cexpdc).

````{prf:example}
:label: eg-cexpdc

If $x$ and $y$ are random variables and $p(y \given x)$ is the conditional density of $y$ given $x$, then

```{math}
:enumerated: false

\EE[ y \given x] = \int t p(t \given x) \diff t
```
````


Here are some goodies we can harvest using the fact that $\EE[\cdot \given \gG] = \proj L_2(\gG)$.

```{prf:property}
:label: fa-coce

Let $x$ and $y$ be random variables in $L_2$, let $\alpha$ and $\beta$ be scalars, and let $\gG$ and $\hH$ be subsets of $L_2$. The following properties hold:

(i) Linearity: $\EE[ \alpha x + \beta y \given \gG] = \alpha \EE[ x \given \gG ] + \beta \EE[ y \given \gG]$.

(ii) If $\gG \subset \hH$, then $\EE[ \EE[ y \given \hH] \given \gG] = \EE[ y \given \gG]$ and $\EE[ \EE[ y \given \gG] ] = \EE[y]$.

(iii) If $y$ is independent of the variables in $\gG$, then $\EE[ y \given \gG] = \EE[y]$.

(iv) If $y$ is $\gG$-measurable, then $\EE[ y \given \gG] = y$.

(v) If $x$ is $\gG$-measurable, then $\EE[ x y \given \gG] = x \EE[y \given \gG]$.
```


Some of these claims are given as solved exercises. Others follow directly from the definition and from the properties of projections. For example, in (iv) we are saying that if $y \in L_2(\gG)$, then $y$ is projected into itself. This is immediate from (vi) in [fact %s](#fa-l2opt).

Part (ii) is called the **the law of iterated expectations**. The law follows from the property of orthogonal projections given in [fact %s](#fa-subsub): Projecting onto the bigger subspace $L_2(\hH)$ and from there onto $L_2(\gG)$ is the same as projecting directly onto the smaller subspace $L_2(\gG)$.

Part (v) is called **conditional determinism**. It says that if $x$ is known given $\gG$ then we can treat it as a constant in the conditional expectation.

#### Conditional Expectations and Prediction

To recap, given $y \in L_2$ and random vector $\boldx$ in $\RR^D$, the conditional expectation $\EE[y \given \boldx]$ is a function $f^*$ of $\boldx$, called the regression function of $y$ given $\boldx$, such that $f^*(\boldx)$ is the best predictor of $y$ in terms of MSE. In other words,

$$
f^* = \argmin_{g \in G} \EE[ (y - g(\boldx))^2 ]
$$ (eq-mmse)

where $G$ is the functions from $\RR^D$ to $\RR$ with $g(\boldx) \in L_2$. In fact a more explicit result can be obtained: For any $g \in G$, we have

$$
\EE [ (y - g(\boldx))^2 ]
    = \EE [ (y - f^*(\boldx) )^2 ] +  \EE[ (f^*(\boldx) - g(\boldx))^2 ]
$$ (eq-ebrf)

This implies {eq}`eq-mmse` because $(f^*(\boldx) - g(\boldx))^2 \geq 0$.

To prove {eq}`eq-ebrf`, let $f^*$ be the regression function, pick any $g \in G$ and observe that

```{math}
:enumerated: false

\begin{aligned}
(y - g(\boldx))^2 
    & = (y - f^*(\boldx) + f^*(\boldx) - g(\boldx))^2  \\
    & = (y - f^*(\boldx) )^2 + 2(y - f^*(\boldx))(f^*(\boldx) - g(\boldx)) +
    (f^*(\boldx) - g(\boldx))^2
\end{aligned}
```

Consider the expectation of the cross-product term. From the law of iterated expectations,

$$
\EE \{ (y - f^*(\boldx)) (f^*(\boldx) - g(\boldx)) \}
     = \EE \{ \, \EE[ (y - f^*(\boldx)) (f^*(\boldx) - g(\boldx)) \given \boldx ]\,  \}
$$ (eq-crop)

Using conditional determinism we can rewrite the term inside the curly brackets on the right-hand side of {eq}`eq-crop` as $(f^*(\boldx) - g(\boldx)) \EE[ (y - f^*(\boldx))  \given \boldx ]$. For the second term in this product we have

```{math}
:enumerated: false

\EE[ y - f^*(\boldx) \given \boldx] 
     = \EE[ y \given \boldx] - \EE[ f^*(\boldx) \given \boldx] 
     = \EE[ y \given \boldx] -  f^*(\boldx)
     = 0
```

Hence the expectation in {eq}`eq-crop` is zero. Equation {eq}`eq-ebrf` now follows.

### The Vector Case

Conditional expectations of random matrices are defined using the notion of conditional expectations for scalar random variables. For example, given random matrices $\boldX$ and $\boldY$, we set

```{math}
:enumerated: false

\EE [\boldY \given \boldX]
    := 
    \left(
    \begin{array}{cccc}
        \EE [y_{11} \given \boldX]  & \cdots & \EE [y_{1K} \given \boldX] \\
        \vdots &  & \vdots \\
        \EE [y_{N1} \given \boldX]  & \cdots & \EE [y_{NK} \given \boldX]
    \end{array}
    \right)
```

We also define

(i) $\cov [\boldx, \boldy \given \boldZ] := \EE[ \boldx \boldy^\T \given \boldZ ] - \EE[\boldx\given \boldZ ] \EE[\boldy\given \boldZ ]^\T$

(ii) $\var[\boldx \given \boldZ] := \EE [\boldx \boldx^\T\given \boldZ ] - \EE [\boldx\given \boldZ ] \EE[\boldx\given \boldZ ]^\T$

Using the definitions, one can show that the properties of scalar conditional expectations in [fact %s](#fa-coce) carry over to the matrix setting. Here's a partial list:

```{prf:property}
:label: fa-liconex2

If $\boldX$, $\boldY$ and $\boldZ$ are random matrices and $\boldA$ and $\boldB$ are constant and conformable, then

(i) $\EE [ \boldY \given \boldZ ]^\T = \EE [ \boldY^\T \given \boldZ ]$.

(ii) $\EE [ \boldA \boldX + \boldB \boldY \given \boldZ ] = \boldA \EE[\boldX \given \boldZ ] + \boldB \EE[ \boldY \given \boldZ  ]$.

(iii) $\EE [\EE [ \boldY \given \boldX ] ] = \EE[\boldY]$ and $\EE [\EE [ \boldY \given \boldX, \boldZ ] \given \boldX ] = \EE[\boldY \given \boldX]$.

(iv) If $\boldX$ and $\boldY$ are independent, then $\EE [ \boldY \given \boldX ]  = \EE[\boldY]$.

(v) If $g(\boldX)$ is a matrix depending only on $\boldX$, then

    1.  $\EE [ g(\boldX) \given \boldX ] = g(\boldX)$

    2.  $\EE [ g(\boldX) \, \boldY \given \boldX ] = g(\boldX) \EE [ \boldY \given \boldX ]$ and $\EE [ \boldY \, g(\boldX) \given \boldX ] = \EE [ \boldY \given \boldX] \, g(\boldX)$
```


Part (v) shows examples of conditional determinism in the matrix case.

## Further Reading

Comprehensive treatments of multivariate distributions can be found in monographs such as {cite}`izenman2008modern` and {cite}`joe1997multivariate`. Beautiful treatments of conditional expectation via orthogonal projection in $L_2$ can be found in {cite}`williams1991probability`, {cite}`taylor2012introduction`, {cite}`schilling2005measures`, and {cite}`cinlar`.

## Exercises

```{exercise}
:label: ex-trden

Confirm the expression for the density $q$ of $y = \psi(x)$ in {eq}`eq-trden` in the case where $\psi$ is strictly increasing.[^5]
```


```{exercise}
:label: ex-ch05_dependence-auto-1

Let $x_1$ and $x_2$ be random variables with joint density $p$ and marginals $p_1$ and $p_2$. Show that $x_1$ and $x_2$ are independent in the sense that {eq}`eq-pdind` holds whenever $p(s_1,s_2) = p_1(s_1) p_2(s_2)$ for every $(s_1, s_2) \in \RR^2$.[^6]
```


````{exercise}
:label: ex-siex

Let $p \colon \RR^2 \to \RR$ be defined by

```{math}
:enumerated: false

p(s_1, s_2) = 
    \begin{cases}
        \exp(s_1 + s_2) & \text{if } s_1 \leq 0 \text{ and } s_2 \leq 0 \\
        0 & \text{otherwise}
    \end{cases}
```

Let $x_1, x_2$ have joint density $p$. Show that $x_1$ and $x_2$ are independent.
````


```{exercise}
:label: ex-chis1

Show that if $\lL(x) = \nN(0, 1)$, then $\lL(x^2) = \chi^2(1)$.[^7]
```


```{exercise}
:label: ex-ch05_dependence-auto-2

Confirm [fact %s](#fa-indzcov): If $\boldx$ and $\boldy$ are independent, then $\cov[\boldx,\boldy] = 0$.
```


```{exercise}
:label: ex-ostatu

Let $x$ and $y$ be independent uniform random variables on $[0,1]$. Let $z := \max\{x, y\}$. Compute the CDF, density and mean of $z$.[^8] In addition, compute the CDF of $w := \min\{x, y\}$.
```


```{exercise}
:label: ex-eoqf

Confirm [fact %s](#fa-eoqf).
```


```{exercise}
:label: ex-ch05_dependence-auto-3

Let $\boldx$ be a random vector. Confirm the claim in [fact %s](#fa-ltavc) that $\var[\boldx]$ is nonnegative definite.
```


```{exercise}
:label: ex-fnwg

Let $\lL(\boldx) = \lL(x_1,\ldots, x_N) = \nN(\boldzero, \boldI)$.

(i) What is the distribution of $x_1^2/x_2^2$?

(ii) What is the distribution of $x_1  [2 / (x_2^2 + x_3^2)]^{1/2}$?

(iii) What is the distribution of $\| \boldx \|^2$?

(iv) If $\bolda$ is an $N \times 1$ constant vector, what is the distribution of $\bolda^\T \boldx$?
```


```{exercise}
:label: ex-discreteconv

Let $x$ and $y$ be independent random variables taking values in $\{0\} \cup \NN$ with PMFs $p = \{p(j)\}$ and $q = \{q(j)\}$ respectively. (I'm using function rather than subscript style notation for this exercise to make the arguments a bit easier to follow.) The **convolution** of $p$ and $q$ is the PMF $p \star q$ defined by

$$
(p \star q)(j) = \sum_{k \geq 0} p(j - k) q(k)
$$ (eq-discreteconv)
Show that $p \star q$ is the PMF of $z := x + y$ when $x$ and $y$ are independent with $\lL(x) = p$ and $\lL(y) = q$.
```


```{exercise}
:label: ex-poci

Prove [fact %s](#fa-poci).
```


```{exercise}
:label: ex-lssl2eg

Fix $z \in L_2$ with $\EE z=0$ and let $S := \setntn{x \in L_2}{ \cov[ x, z] = 0}$. Show that $S$ is a linear subspace of $L_2$.
```


```{exercise}
:label: ex-pcs

Let $x \in L_2$ with $\EE x = \mu$ and $\var[x] = \sigma^2$. Let $S := \Span\{\1\}$, as in [example %s](#eg-projam). Let $\boldP := \proj S$ and let $\boldM$ be the residual projection defined by $\boldM x = x - \boldP x$. Show that $\| \boldM x \| = \sigma$.
```


```{exercise}
:label: ex-precsrv

Let $x, y \in L_2$, let $S := \Span\{x\}$, and let $\boldP := \proj S$. Show that

$$
\boldP y = \EE[x^2]^{-1} \EE[xy]\, x
$$ (eq-pryonxa)
```


```{exercise}
:label: ex-csrv

Prove the Cauchy–Schwarz inequality for random variables, which was first stated as [fact %s](#fa-csrv) as $| \EE [xy] | \leq ( \EE[x^2] \EE[y^2] )^{1/2}$ and repeated as part (iv) of [fact %s](#fa-l2se). Use the results on orthogonal projections in $L_2$ as found in {ref}`ss-projl2`.[^9]
```


```{exercise}
:label: ex-mpr

Prove that $\hat y$ in {eq}`eq-mpr` is the orthogonal projection of $y$ onto $\Span\{\boldx\}$ under the conditions stated in [fact %s](#fa-mpr).
```


```{exercise}
:label: ex-rtblp0

Let $\boldx \in \RR^K$ be a random vector and let $\| \cdot \|$ be the ordinary vector norm. Show that if $\| \boldx \|$ is in $L_2$, then the random variable $z = \ell(\boldx)$ is in $L_2$ whenever $\ell \colon \RR^K \to \RR$ is linear.
```


```{exercise}
:label: ex-rtblp

Consider the setting of the population regression discussion in {ref}`sss-pr`. In particular, let $y$ be an element of $L_2$, let $\boldb^*$ be the coefficients of the best linear predictor, as defined in {eq}`eq-mpr`, let $S = \Span\{\boldx\}$, and let $f^*(\boldx) = \EE[y \given \boldx]$. Show that $\ell^*(\boldx) := \boldx^\T \boldb^*$ is the orthogonal projection of $f^*(\boldx)$ onto $S$.
```


```{exercise}
:label: ex-ivce

Prove part (iii) of [fact %s](#fa-coce) when $\gG = \{\boldx\}$.
```


```{exercise}
:label: ex-lccoce

Prove part (v) of [fact %s](#fa-coce).
```


````{exercise}
:label: ex-ch05_dependence-auto-4

Let $\var[y \given x] := \EE[y^2 \given x] - (\EE[y \given x])^2$. Show that

```{math}
:enumerated: false

\var y = \EE [\var[y \given x]] + \var[\EE[y \given x]]
```
````


```{exercise}
:label: ex-ceoac

Show that the conditional expectation of a constant $\alpha$ is $\alpha$. In particular, using the results in [fact %s](#fa-coce) as appropriate, show that if $\alpha$ is a constant and $\gG$ is any information set, then $\EE[\alpha \given \gG] = \alpha$.
```


```{exercise}
:label: ex-cexpdc

Prove the claim in [example %s](#eg-cexpdc).
```


### Solutions to Selected Exercises

````{solution} ex-trden

Let $x$, $p$ and $y = \psi(x)$ be as stated in [fact %s](#fa-trden). Assume in addition that $\psi$ is strictly increasing. Let $F$ be the CDF of $x$ and let $G$ be the CDF of $y$. By [fact %s](#fa-cdftr) we have $G(s) = F(\psi^{-1}(s))$. Differentiating gives

```{math}
:enumerated: false

q(s) = G'(s) 
        = p(\psi^{-1}(s)) \, \frac{\diff \psi^{-1}(s)}{\diff s}
```

Since $\psi^{-1}$ is increasing, this agrees with {eq}`eq-trden`.
````


````{solution} ex-siex

Let $p_1$ be the marginal density of $x_1$ and $p_2$ be the marginal density of $x_2$. It suffices to show that

$$
p(s_1, s_2) = p_1(s_1) p_2(s_2)
        \quad
        \text{for all } (s_1, s_2) \in \RR^2
$$ (eq-m)
To see that this is the case, first we obtain the marginal distributions by integrating out the other variable. To start, let $s_1 \leq 0$ and observe that

```{math}
:enumerated: false

p_1(s_1) 
    = \int_{-\infty}^0 p(s_1, s_2) \diff s_2
    = \int_{-\infty}^0 e^{s_1} e^{s_2} \diff s_2
    = e^{s_1} \int_{-\infty}^0 e^{s_2} \diff s_2
    = e^{s_1}
```

Conversely, if $s_1 > 0$, then $p(s_1, s_2) = 0$ for any $s_2$, and hence

```{math}
:enumerated: false

p_1(s_1) 
    = \int_{-\infty}^0 p(s_1, s_2) \diff s_2 = 0
```

From this and a symmetric argument in the case of $p_2$, we get

```{math}
:enumerated: false

p_1(s_1) = 
    \begin{cases}
        e^{s_1} & \text{if } s_1 \leq 0 \\
        0 & \text{otherwise}
    \end{cases}
    \qquad \text{and} \qquad
    p_2(s_2) = 
    \begin{cases}
        e^{s_2} & \text{if } s_2 \leq 0 \\
        0 & \text{otherwise}
    \end{cases}
```

It follows that {eq}`eq-m` is valid.
````


````{solution} ex-chis1

Let $F(s)$ be the CDF of $x^2$ and let $f$ be the density. Let $\Phi$ and $\phi$ be the $\nN(0, 1)$ cdf and density respectively. Following arguments similar to the solution to [exercise %s](#ex-avsc) we have

```{math}
:enumerated: false

F(s) = \PP\{x^2 \leq s\} = \PP\{-\sqrt{s} \leq x \leq \sqrt{s}\} 
            = 2 \Phi(\sqrt{s}) - 1
```

Differentiating with respect to $s$ gives $f(s) = \phi(\sqrt{s}) s^{-1/2}$. Plugging in the definition of the $\nN(0,1)$ density and using $\Gamma(1/2) = \sqrt{\pi}$ gives the density in [example %s](#eg-chisq) when $k=1$.
````


````{solution} ex-ostatu

As in the statement of the exercise, $x$ and $y$ are independent uniform random variables on $[0,1]$, $z := \max\{x, y\}$ and $w := \min\{x, y\}$. As a first step to the proofs, note that if $a$, $b$ and $c$ are three numbers, then

- $\max\{a, b\} \leq c$ if and only if $a \leq c$ and $b \leq c$.

- $\min\{a, b\} \leq c$ if and only if $a \leq c$ or $b \leq c$.

From these facts we can see that, for any $s \in \RR$,

```{math}
:enumerated: false

\{z \leq s\} = \{x \leq s\} \cap \{y \leq s\}
        \quad \text{and} \quad
        \{w \leq s\} = \{x \leq s\} \cup \{y \leq s\}
```

Now, for $s \in [0,1]$, we have

```{math}
:enumerated: false

\PP\{z \leq s\} 
        = \PP [ \{x \leq s\} \cap \{y \leq s\} ]
        = \PP\{x \leq s\} \PP\{y \leq s\} 
        = s^2
```

By differentiating, we get the density $p(s) = 2s$, and by integrating $\int_0^1 s p(s) \diff s$, we get $\EE[z] = 2/3$. Finally, regarding the CDF of $w$, for $s \in [0,1]$, we have

```{math}
:enumerated: false

\begin{aligned}
\PP\{w \leq s\} 
        &= \PP [ \{x \leq s\} \cup \{y \leq s\} ] \\
        &= \PP  \{x \leq s\} + \PP  \{y \leq s\} - \PP [ \{x \leq s\} \cap \{y \leq s\} ]
\end{aligned}
```

Hence $\PP\{w \leq s\} = 2s - s^2$.
````


````{solution} ex-eoqf

Since $\boldz^\T \boldA \boldz$ has the form $\sum_{i=1}^N \sum_{j=1}^N  a_{ij} z_i z_j$, under the stated assumptions we have

```{math}
:enumerated: false

\EE [\boldz^\T \boldA \boldz]
        = \sum_{i=1}^N \sum_{j=1}^N  a_{ij} \EE[ z_i z_j ]
        = \trace \boldA
        
```
````


```{solution} ex-fnwg

The respective solutions are as follows:

(i) $\lL(x_1^2/x_2^2) = F(1,1)$ by [fact %s](#fa-rcif).

(ii) $\lL[x_1  [2 / (x_2^2 + x_3^2)]^{1/2}] = t(2)$ by [fact %s](#fa-astudt).

(iii) $\lL(\| \boldx \|^2) = \chi^2(N)$ by [fact %s](#fa-ssnac).

(iv) $\lL(\bolda^\T \boldx) = \nN(0, \|\bolda\|^2)$ by [fact %s](#fa-nipre).
```


````{solution} ex-discreteconv

Let $p \star q$ be as defined in {eq}`eq-discreteconv`, let $z := x + y$ where $x$ and $y$ are independent, and let $\lL(x) = p$ and $\lL(y) = q$. The claim is that $(p \star q)(j) = \PP\{z = j\}$ for any $j \geq 0$. To see this, fix any such $j$ and apply the law of total probability to get

```{math}
:enumerated: false

\PP\{z = j\}
         = \sum_{k \geq 0} \PP\{ z = j \given y = k\} \PP\{y = k \}
         = \sum_{k \geq 0} \PP\{ x + k = j \} q(k)
```

Using $\PP\{ x + k = j \} =\PP\{ x = j -k \} = p(j - k)$ completes the proof.
````


````{solution} ex-poci

Let $S$, $\boldP$, $\boldM$, and $\boldz$ be as in the statement of the question. In view of [fact %s](#fa-ucimin), to prove independence of $\boldP \boldz$ and $\boldM \boldz$ it suffices to show that $\cov[\boldP \boldz, \boldM \boldz] = \boldzero$. Since these vectors are zero-mean and $\boldM$ is symmetric , the covariance is equal to

```{math}
:enumerated: false

\EE[ (\boldP \boldz) (\boldM \boldz)^\T] 
        =\EE[ \boldP \boldz \boldz^\T \boldM] 
        =\boldP \EE[ \boldz \boldz^\T]  \boldM 
        = \sigma^2 \boldP \boldM 
        = \boldzero
```

The final equality follows from [fact %s](#fa-opt3).
````


````{solution} ex-lssl2eg

Let $z$ and $S$ be as described in the question. $S$ is a linear subspace of $L_2$ because if $x, y \in S$ and $\alpha, \beta \in \RR$, then

```{math}
:enumerated: false

\cov [ \alpha x + \beta y, z]
        = \alpha \EE [ xz ]  + \beta \EE[yz]
        = \alpha \cov [ x,z ]  + \beta \cov[y,z]
        = 0
        
```
````


```{solution} ex-pcs

We saw in [example %s](#eg-projam) that $\boldP x = \mu \1$, so $\boldM x = x - \mu \1$, or, written more conventionally, $\boldM x = x - \mu$. The norm of $\boldM x$ is $\| \boldM x \| = \sqrt{\EE[ (x - \mu)^2 ]} = \sigma$, as was to be shown.
```


```{solution} ex-precsrv

Let $x, y, S$, and $\boldP$ be as in the statement of the exercise, and let $\boldP y$ be as defined in {eq}`eq-pryonxa`. As a scalar multiple of $x$ it's clear that $\boldP y \in S$, so we need only show that $y - \boldP y \perp x$, or, equivalently, $\EE[xy - x\boldP y] = 0$. Substituting in the expression for $\boldP y$ in {eq}`eq-pryonxa` and using linearity of $\EE$ confirms this equality.
```


````{solution} ex-mpr

Assume the conditions of [fact %s](#fa-mpr) and adopt the notation stated there. Let $y$ be any element of $L_2$. In view of [theorem %s](#t-opt4), to show that $\hat y = \boldP y$, we need to show that $\hat y \in S := \Span\{\boldx\}$ and that $y - \hat y \perp S$. The first part is trivial, since $\hat y = \boldx^\T \boldb^*$ where $\boldb^*$ is the vector $\EE[\boldx \boldx^\T]^{-1} \EE[\boldx y]$. To prove the second part, it suffices to show that

```{math}
:enumerated: false

\EE[ \bolda^\T \boldx y] 
        = \EE[\bolda^\T \boldx \boldx^\T \boldb^* ]
```

Substituting in the definition of $\boldb^*$ and some easy manipulations confirm that this equality holds.
````


```{solution} ex-rtblp0

Let $z = \ell(\boldx)$ where $\ell$ is a linear function. The claim is that $z$ is square integrable. By [theorem %s](#t-lmaeq), we can write $z$ as $z = \boldx^\T \bolda$ for some $\bolda \in \RR^K$. Applying the Cauchy–Schwarz inequality for vectors, we have $z^2 \leq \| \boldx\|^2 \, \| \bolda \|^2$. Since $\| \boldx \|$ is in $L_2$, it follows that $\EE[z^2] < \infty$.
```


```{solution} ex-rtblp

Adopt the setting and notation of [exercise %s](#ex-rtblp). Let $\gG$ be the set of all $\boldx$-measurable functions in $L_2$. Let $S = \Span\{\boldx\} =$ the linear $\boldx$-measurable functions in $L_2$. (As in {ref}`sss-pr`, we are assuming that $\| \boldx \|$ is in $L_2$, which assures us that elements of $S$ are in fact in $L_2$. See [ex. %s](#ex-rtblp0).) Clearly $S \subset \gG$. The result now follows from [fact %s](#fa-subsubl2).
```


```{solution} ex-ivce

Let $y$ be independent of $\boldx$. From the definition of conditional expectation, to show that $\EE[ y \given \boldx] = \EE[y]$ we need to show that

(i) $\EE[y]$ is $\gG$-measurable, and

(ii) $\EE[ \EE[y] g(\boldx)] = \EE[ y g(\boldx)]$ for any $\bB$-measurable $g$.

Part (i) is trivial because $\EE[y]$ is constant (see [example %s](#eg-cmag)). Regarding (ii), if $g$ is any such function, then by [facts %s](#fa-imme) and [%s](#fa-rviifi) we have $\EE[ y g(\boldx)] = \EE[y] \EE[g(\boldx)]$. The result follows.
```


````{solution} ex-lccoce

We need to show that if $x$ is $\gG$-measurable, then $\EE[ x y \given \gG] = x \EE[y \given \gG]$. To confirm this, we must show that

(i) $x \EE[y \given \gG]$ is $\gG$-measurable, and

(ii) $\EE[x \EE[y \given \gG] z] = \EE[ xy z]$ for any $z \in L_2(\gG)$.

Regarding (i), $\EE[y \given \gG]$ is $\gG$-measurable by definition, and $x$ is $\gG$-measurable by assumption, so $x\EE[y \given \gG]$ is $\gG$-measurable by [fact %s](#fa-lcapom). Regarding (ii), fix $z \in L_2(\gG)$, and let $u:=xz$. Since $x \in L_2(\gG)$, we have $u \in L_2(\gG)$. We need to show that

```{math}
:enumerated: false

\EE[ \EE[y \given \gG] u ] = \EE[ y u ]
```

Since $u \in L_2(\gG)$, this is immediate from the definition of $\EE[y \given \gG]$.
````


```{solution} ex-ceoac

By [fact %s](#fa-coce), we know that if $\alpha$ is $\gG$-measurable, then $\EE[ \alpha \given \gG] = \alpha$. {prf:ref}`eg-cmag` tells us that $\alpha$ is indeed $\gG$-measurable.
```


````{solution} ex-cexpdc

As in [example %s](#eg-cexpdc), let $x$ and $y$ be random variables where $p(y \given x)$ is the conditional density of $y$ given $x$. Let $g(x) := \int t p(t \given x) \diff t$. The claim is that $\EE[ y \given x]  = g(x)$. To prove this, we need to show that $g(x)$ is $x$-measurable, and that

$$
\EE[ g(x) h(x) ] = \EE[ y h(x) ]  
        \quad \text{for any } \bB \text{-measurable } h
$$ (eq-cexpdcp)
The first claim is obvious. Regarding {eq}`eq-cexpdcp`, let $h$ be any such function. Using the notation in {eq}`eq-condden`, we can write

```{math}
:enumerated: false

\begin{aligned}
\EE[ g(x) h(x) ] 
        & = \EE \left[ \int t p(t \given x) \diff t \; h(x) \right] \\
        & = \int \int t p(t \given s) \diff t \; h(s)  p(s) \diff s \\
        & = \int \int t \frac{p(s, t)}{p(s)} \diff t \; h(s)  p(s) \diff s 
         = \int \int t \, h(s) p(s, t) \diff t \diff s
\end{aligned}
```

This is equal to the right-hand side of {eq}`eq-cexpdcp`, and the proof is done.
````


[^1]: For technical details on these constructions see, for example, chapter 4 of {cite}`cinlar`.

[^2]: See page  for more discussion of rectangles.

[^3]: Uniqueness is a probability one statement: Any two solutions $\hat y_1$ and $\hat y_2$ satisfy $\PP\{\hat y_1 = \hat y_2\}=1$. We've already agreed not to distinguish between such elements. See the discussion after [fact %s](#fa-l2se).

[^4]: A tighter proof can be given using the language of $\sigma$-algebras, although the intuition is less immediate.

[^5]: Hint: Use [fact %s](#fa-cdftr).

[^6]: Hint: Apply {eq}`eq-defjd`.

[^7]: Hints: You need to match the density in [example %s](#eg-chisq) when $k=1$. Use the fact that $\Gamma(1/2) = \sqrt{\pi}$. Try to find the CDF of $z^2$ and differentiate it. In finding the CDF, the solution to [exercise %s](#ex-avsc) might be helpful.

[^8]: Hint: Fix $s \in \RR$ and compare the sets $\{z \leq s\}$ and $\{x \leq s\} \cap \{y \leq s\}$. What is the relationship between these two sets?

[^9]: Hint: Use {eq}`eq-pryonxa` and review the solution to [exercise %s](#ex-pcsin).
