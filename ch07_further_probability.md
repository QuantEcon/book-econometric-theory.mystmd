---
title: "Further Topics in Probability"
label: c-tsm
---

There are two important topics in probability that we haven't yet managed to cover. One is stochastic process theory and the other is simulation. Let's go over some of the key results.

## Stochastic Processes

A **stochastic process** on $\RR^K$ is a sequence of random vectors $\{\boldx_t\}_{t\geq 0}$ defined on a common probability space $(\Omega, \fF, \PP)$. The definition has relatively little content, but the fact that all these random vectors live on the same probability space should not be overlooked. It means that probabilities such as $\PP\{\boldx_t \geq 0 \text{ for all } t\}$ or $\PP\{\boldx_t \leq -10 \text{ for some } t\}$ are well defined. Putting the elements of the process on a common probability space pins down the joint distribution across $t$.

To do any kind of analysis with stochastic processes, we need to put some structure on the joint distribution, or, equivalently, on the random elements. In practice, this is almost always done recursively. Examples are given below. For the moment, however, let's look at stochastic processes in an abstract setting and state some definitions.

One kind of stochastic process is an IID sequence. Such sequences are ideal to work with in a statistical setting. For example, the law of large numbers tells us that in an IID environment, finite expectations can always be approximated with sample averages. If we want to know how good the approximation is, we can look to the central limit theorem.

Once we move out of the IID world, however, these helpful properties can break down. Without restrictions on our stochastic process, we can say very little about estimation. The first task, then is to describe some properties that make our stochastic process resemble an IID process in *some* way.

(ss-spdef)=
### Stationarity and Ergodicity
As a first step, let's define a type of process that preserves one of the properties of IID processes (the "ID\" part). A stochastic process $\{ \boldx_t \}_{t \geq 0}$ is called **stationary** if it has the property that the distribution of any subset of these random vectors is unaffected by shifting them forward in time. That is,

$$
\lL ( \boldx_{t_1}, \ldots, \boldx_{t_k} )
    =
    \lL ( \boldx_{t_1 + m}, \ldots, \boldx_{t_k + m} )
$$ (eq-statproc)

for any $m \in \NN$ and any sequence of integers $t_1, \ldots, t_k$.

Clearly, any stationary sequence is identically distributed. In other words, for any stationary sequence, the sequence of marginals $\{\lL(\boldx_t)\}$ is a constant sequence.

```{prf:example}

An IID process $\{\boldx_t\}$ is stationary, since independence implies that the laws on both sides of {eq}`eq-statproc` are just $k$ products of the marginal $\lL(\boldx_1)$.
```


````{prf:example}
:label: eg-llnfail0

For a simple example of a stationary process, let $x_1$ be a draw from some arbitrary law $P$ on $\RR$ and let $x_{t+1} = x_t$ for all $t$. For any collection of integers $t_1, \ldots, t_k$ and any Borel sets $B_1, \ldots, B_k$, we have

```{math}
:enumerated: false

\PP\{x_{t_1} \in B_1, \ldots, x_{t_k} \in B_k\}
        = \PP \{x_1 \in \cap_{i=1}^k B_k \}
        = P ( \cap_{i=1}^k B_k )
```

If we replace each $x_{t_i}$ with $x_{t_i + m}$, we get the same number. Hence $\{x_t\}$ is stationary.
````


```{prf:example}

A **random walk** on $\RR$ is a stochastic process $\{x_t\}$ where $x_t = \sum_{j=1}^t w_j$ for some IID zero-mean process $\{w_j\}$. Let $\sigma^2 := \var w_t$. If $\sigma > 0$, then $\{x_t\}$ is not stationary, since $\var x_t = t \sigma^2$. In particular, $\lL(x_t)$ depends on $t$.
```


Stationarity alone isn't strong enough to give us properties like the LLN and CLT.

```{prf:example}
:label: eg-llnfail

For the process $\{x_t\}$ in [example %s](#eg-llnfail0), the LLN fails whenever $P$ is nondegenerate. Indeed $\bar x_T := \frac{1}{T} \sum_{t=1}^T x_t  = x_1$, and hence $\lL(\bar x_T) = \lL(x_1) = P$ for all $t$. In particular, $\bar x_T$ does not converge in probability to any constant.
```


The problem with the process in [example %s](#eg-llnfail) is excessive dependence. New realizations provide insufficient (in fact no) new information. On the flip side, dependent processes can satisfy versions of the LLN if dependence is not too excessive (e.g., if $x_t$ and $x_{t+k}$ are "almost independent\" for $k$ sufficiently large.) These are the so-called ergodic processes.

The formal definition of ergodicity varies across different literatures. Here we define a stationary stochastic process $\{\boldx_t\}$ on $\RR^K$ to be **ergodic** if the LLN holds; that is, if

$$
\frac{1}{T} \sum_{t=1}^T h(\boldx_t) \toprob \mu_h := \EE h(\boldx_t) 
         \quad \text{ as } \quad T \to \infty
$$ (eq-ergodic)

for any $\bB$-measurable $h \colon \RR^N \to \RR$ such that $\EE |h(\boldx_t)|  < \infty$. (Since $\{\boldx_t\}$ is stationary, this expectation does not depend on $t$.) One of the main tasks of this chapter will be deriving sufficient conditions for ergodicity in the kinds of processes economists often work with.

Incidentally, a casual definition of ergodicity one often hears is that "cross sectional and time series averages coincide.\" What does that mean in the present context?

To answer the question, imagine that $\boldx_t = x_t$ is a binary variable with $x_t = 1$ if a given individual is employed at time $t$ and zero otherwise. In this case, $\frac{1}{T} \sum_{t=1}^T x_t$ is the fraction of time spent employed over the period $1, \ldots, T$. With $h(x)=x$, the right-hand side of {eq}`eq-ergodic` is $\EE x_t = \PP\{x_t = 1\}$, or the stationary probability of employment. If we observe time $t$ employment outcomes $x_t^1, \ldots, x_t^N$ of individuals $1, \ldots, N$ and, moreover, the individuals all follow the same model and are sufficiently independent, the cross sectional average $\frac{1}{N} \sum_{n=1}^N x_t^n$ will also be close to $\PP\{x_t = 1\}$. Hence, with ergodicity and large samples, cross-sectional and time series averages (almost) coincide.

(ss-srsq)=
### Stochastic Recursive Sequences
In this section let's turn to a more concrete view of stochastic processes. Afterward we'll circle back and consider issues like ergodicity again.

{numref}`f-gdp` shows the evolution of US GDP per capita since the early postwar era. The top panel is per capita GDP itself, while the second panel is the growth rate. How could we go about modeling such time series in a parsimonious way?

```{figure} figures/gdp.svg
:name: f-gdp

Postwar US GDP per capita. Source: Penn World Tables
```

No reasonable model can be deterministic, since shocks play an important role in year-on-year change in GDP, and an IID process is too simple. We can, however, get close to both of these time series with one of the most common nontrivial stochastic processes: the scalar Gaussian AR(1) model

$$
x_{t+1} = b + a x_t + c w_{t+1}
    \quad \text{with } \{ w_t \} \iidsim \nN(0,1)
    \text{ and } x_0 \text{ given}
$$ (eq-sgar1)

Here $a, b$, and $c$ are parameters, while $x_t$ is called the **state variable**. Equation {eq}`eq-sgar1` is an example of a **stochastic difference equation**. The process $\{x_t\}$ it defines is called a stochastic process or **stochastic recursive sequence**. A realization of the process is called a **time series**.

The dynamics of $\{x_t\}$ depends on the parameters. We can get some sense of this from simulation. Six individual time series $\{x_t\}$ are shown in [figure %s](#f-ar1_dynam), each generated using a different value of $a$. The values of $b$ and $c$ are both held fixed at 1, while $x_0 = 0$. `R` code for generating the figure is given in [listing %s](#list-ar1_dynam).

As suggested by the figure, the simulated time paths are sensitive to the value of $a$. If $a$ is outside the interval $(-1,1)$, the series tend to diverge. If $|a| < 1$, then the opposite is true. For example, if you look at the time series for $a = 0.9$ in [figure %s](#f-ar1_dynam), you will see that, after an initial burn in period where the series is affected by the initial condition $x_0$, the process settles down to random motion within a band (between about 5 and 15 in this case).

```{figure} figures/ar1_dynam.svg
:name: f-ar1_dynam

Dynamics of the linear AR(1) model
```

```{code-block} r
:caption: Code for figure {numref}`f-ar1_dynam` (R)
:name: list-ar1_dynam
# Generates an AR(1) time series starting from x = init
ar1ts <- function(init, n, b, a) {
    x <- numeric(n)
    x[1] <- init 
    w <- rnorm(n-1)
    for (t in 1:(n-1)) {
        x[t+1] <- b + a * x[t] + w[t]
    }   
    return(x)
}   

avec <- c(0.1, -0.1, 0.9, -0.9, 1.1, -1.1)
N <- 200
par(mfrow=c(3,2))  # Arrangement of figures 
for (a in avec) {
    plot(ar1ts(0, N, 1, a), type="l", xlab=paste("a = ", a), ylab="")
}   
```

If the scalar Gaussian AR(1) model in {eq}`eq-sgar1` is unable to capture important aspects of a time series that we wish to model, we can generalize it in several different directions. For example, we can remove the assumption that the shocks are normally distributed, in which case it is simply called the scalar AR(1) model. This might be desirable, if, say, the state process has heavy tails.

We can also generalize the scalar AR(1) model to $\RR^K$, yielding the vector AR(1) model, or VAR(1) model

$$
\boldx_{t+1} = \boldb + \boldA \boldx_t + \boldC \boldw_{t+1}
$$ (eq-var1)

Here $\{ \boldw_t \}$ is assumed to be IID and to satisfy $\EE \boldw_t = \boldzero$ and $\EE[\boldw_t\boldw_t'] = \boldI$. If $\boldw_t$ is multivariate normal, then {eq}`eq-var1` is called the Gaussian VAR(1). The vector $\boldx_t$ is called the **state vector**. A huge variety of models can be represented as a VAR(1) process or contain VAR processes as subcomponents (e.g., linear state space models).

An alternative generalization of the scalar AR(1) model is the scalar AR($p$) model, where the next state $x_{t+1}$ is a linear function not just of the current state $x_t$, but of the last $p$ previous states. For example, the AR(2) process has dynamics

```{math}
:enumerated: false

x_{t+1} = b + a x_t + \gamma x_{t-1} + w_{t+1}
```

Although $x_{t+1}$ is a function of two lagged states, $x_t$ and $x_{t-1}$, we can reformulate it as a first-order model. To begin, we define an additional state variable $y_t$ via $y_t = x_{t-1}$. The dynamics can then be expressed as

```{math}
:enumerated: false

\begin{aligned}
x_{t+1}  &= b + a x_t + \gamma y_t + w_{t+1}  \\
    y_{t+1}  &= x_t
\end{aligned}
```

We can write this in matrix form as

```{math}
:enumerated: false

\left( 
    \begin{array}{c}
        x_{t+1} \\
        y_{t+1} 
    \end{array}
    \right)
    = 
    b
    \left( 
    \begin{array}{c}
        1 \\
        0
    \end{array}
    \right)
    +
    \left( 
    \begin{array}{cc}
        a & \gamma \\
        1    & 0
    \end{array}
    \right)
    \left( 
    \begin{array}{c}
        x_t \\
        y_t 
    \end{array}
    \right)
    +
    \left( 
    \begin{array}{c}
        1 \\
        0
    \end{array}
    \right)
    w_{t+1}
```

This is a special case of the VAR(1) model in {eq}`eq-var1`. This example shows how higher order processes can be reduced to first-order processes by increasing the number of state variables. For this reason most theoretical discussion below takes place in the context of first-order models.

(ss-nltsm)=
#### Nonlinear Models
Linear models are simple but that isn't always an advantage. Sometimes they are unable to capture the kinds of dynamics we observe in data. For example, the bursts of volatility in the graph of daily returns on the Hang Seng Index from 1990 to 2014 shown in [figure %s](#f-hsi) are difficult to replicate without introducing nonlinearities into our model. It's also true that many theoretical modeling exercises produce models that are nonlinear.

```{figure} figures/hsi.svg
:name: f-hsi

Daily returns on the Hang Seng Index. Source: Yahoo! Finance
```

One well-known nonlinear model is the $p$th order autoregressive conditional heteroskedasticity model (ARCH($p$) model). The model arose from the need to attach a time varying volatility component $\sigma_t$ to asset return models. Suppose, for example, that $x_t = \sigma_t w_t$ where $\sigma_{t+1}^2 = \alpha_0 + \alpha_1 x_t^2$. Combining these equations gives the ARCH(1) model

$$
x_{t+1} = (\alpha_0 + \alpha_1 x_t^2)^{1/2} w_{t+1}
    \quad \text{with} \quad \{ w_t \} \iidsim \nN(0,1)
$$ (eq-arch)

and $\alpha_0 > 0$, $\alpha_1 \geq 0$. Still better fits to asset returns data can be obtained via Generalized ARCH models, the simplest of which is the GARCH(1,1) process

```{math}
:enumerated: false

\begin{aligned}
x_t &= \sigma_t w_t \\
    \sigma_{t+1}^2 &= \alpha_0 + \alpha_1 x_t^2 + \alpha_2 \sigma_t^2
\end{aligned}
```

Now next period volatility depends on its own lagged state as well as $x_t$.

Another popular nonlinear model is the **smooth transition autoregression** (STAR) model

$$
x_{t+1} = g(x_t) + w_{t+1}
$$ (eq-starrs)

where $g$ is of the form

```{math}
:enumerated: false

g(s) := (b_0 + a_0 s)(1-\tau(s)) + (b_1 + a_1 s)\tau(s)
```

Here $\tau \colon \RR \to [0,1]$ is increasing and satisfies $\lim_{s \to -\infty} \tau(s)=0$ and $\lim_{s \to \infty}\tau(s) =1$. For $s$ near $-\infty$, $g(s)$ is similar to $b_0 + a_0 s$. For $s$ near $+\infty$, $g(s)$ is close to $b_1 + a_1 s$. Thus the dynamics transition between two different linear models, with the smoothness of the transition depending on the shape of $\tau$.

## Markov Processes

All of the examples in {ref}`ss-srsq` are special cases of a general class of process called Markov processes. Markov processes are very rewarding to study because the concepts are broad enough to include a vast array of different stochastic processes, while at the same time having enough structure to support elegant and powerful theory. Here we cover the foundations, from definitions and basic properties to stability and asymptotics.

(ss-mm)=
### The Markov Assumption
Let $S$ be a subset of $\RR^K$ and let $\bB(S)$ be the Borel subsets of $S$. The fundamental primitive of a discrete time Markov process on $S$ is a **stochastic kernel** or **transition probability function**, which is a function $Q \colon S \times \bB(S) \to [0, 1]$ such that

(i) $Q(\bolds, \cdot)$ is a probability measure over $\bB(S)$ for all $\bolds \in S$ and

(ii) $g(\bolds) := Q(\bolds, B)$ is $\bB$-measurable for each $B \in \bB(S)$.

$S$ is called the **state space** of the model. Part (ii) is a mild regularity condition that ensures certain integrals are well defined. Part (i) is the major defining characteristic, yielding a distribution supported on $S$ for every point $\bolds$ in $S$. We think of $Q$ and some initial distribution $P_0$ as generating a corresponding Markov process $\{\boldx_t\}$ taking values in $S$ via

- draw $\boldx_0$ from $P_0$
- **for** $t$ in $0, 1, 2, \ldots$ **do**
  - draw $\boldx_{t+1}$ from the distribution $Q(\boldx_t, \cdot)$
- **end**

The value $Q(\bolds, B)$ is then the probability that $\boldx_{t+1}$ lands in $B$ given $\boldx_t = \bolds$. The sequence $\{\boldx_t\}$ so generated is also called a **sample path** for $Q$.

There is a measure-theoretic procedure for building this process $\{\boldx_t\}$ more formally but the construction is usually unnecessary because most Markov processes are introduced via a stochastic difference equation, which specifies both the process $\{\boldx_t\}$ and, implicitly, the associated stochastic kernel $Q$. The canonical stochastic difference equation for a **first-order Markov process** takes the form

$$
\boldx_{t+1} = G(\boldx_t, \boldw_{t+1})
    \quad \text{with} \quad
    \lL(\boldx_0) = P_0
$$ (eq-amp)

Here $\{\boldw_t\}_{t \geq 1}$ is an IID sequence of $\RR^M$-valued shocks with common distribution $\Psi$ and $G$ is a given $\bB$-measurable function mapping the current state $\boldx_t \in S$ and shock $\boldw_{t+1} \in \RR^M$ into the new state $\boldx_{t+1} \in S$. The initial condition $\boldx_0$ and the shocks $\{\boldw_t\}_{t \geq 1}$ are also assumed to be independent of each other. The corresponding stochastic kernel is

```{math}
:enumerated: false

Q(\bolds, B) = \PP\{G(\bolds, \boldw_{t+1}) \in B\}
    \qquad (\bolds \in S, \; B \in \bB(S))
```

The choice of date $t+1$ is irrelevant here because of the IID assumption, and is used only to tie in with {eq}`eq-amp`. Recalling that $\Psi$ is the distribution of $\boldw_{t+1}$, the stochastic kernel can also be written as

```{math}
:enumerated: false

Q(\bolds, B) 
    = \Psi \left\{\boldw \in \RR^M \; : \; G(\bolds, \boldw) \in B \right\}
```

This expression makes it clear that $Q$ is defined by the two primitives $G$ and $\Psi$.

````{prf:example}
:label: eg-archsk

For the ARCH(1) model of {eq}`eq-arch`, the stochastic kernel is

```{math}
:enumerated: false

Q(s, B) = \PP\{ (\alpha_0 + \alpha_1 s^2)^{1/2} w_{t+1} \in B \}
```

when $w_{t+1}$ is standard normal.
````


Returning to the general process {eq}`eq-amp`, repeated substitution gives

```{math}
:enumerated: false

\begin{aligned}
\boldx_1 & = G(\boldx_0, \boldw_1) \\
    \boldx_2 & = G(G(\boldx_0, \boldw_1), \boldw_2) \\
    \boldx_3 & = G(G(G(\boldx_0, \boldw_1), \boldw_2),\boldw_3)
\end{aligned}
```

and so on. Continuing in this fashion, we see that the state vector $\boldx_t$ can be written as a function of $\boldx_0$ and the shocks $\boldw_1,\ldots,\boldw_t$ for any $t$. In other words, for each $t$, there exists a function $H_t$ such that

$$
\boldx_t = H_t(\boldx_0, \boldw_1, \boldw_2,\ldots, \boldw_t)
$$ (eq-xaft)

```{prf:example}
:label: eg-sar1

For the scalar linear AR(1) process {eq}`eq-sgar1` we can write the right-hand side of {eq}`eq-xaft` explicitly as

$$
x_t = b \sum_{k=0}^{t-1} a^k 
            + \sum_{k=0}^{t-1} a^k c w_{t-k} + a^t x_0
$$ (eq-ar1ft)
This is called the **moving average representation** of $x_t$. The proof is an exercise.
```


Often there is no neat expression for the function $H_t$. Nonetheless, {eq}`eq-xaft` clarifies the fact that {eq}`eq-amp` pins down each $\boldx_t$ as a well-defined random vector, depending on the initial condition and the shocks up until date $t$. Since $\boldx_t$ is a function of $\boldx_0, \boldw_1, \boldw_2,\ldots, \boldw_t$, and these random vectors are all, by the IID assumption, independent of $\boldw_{t+k}$ for any $k \geq 1$, [fact %s](#fa-rviifi) tells us that the current state $\boldx_t$ and future shocks are independent. Let's record this:

```{prf:property}
:label: fa-aamp

For the process {eq}`eq-amp`, the pair $\boldx_t$ and $\boldw_{t+k}$ are independent for all $k \geq 1$.
```


Although it's rarely made explicit, behind all the random vectors $\{\boldx_t\}$ lies a single probability space $(\Omega, \fF, \PP)$. The idea is that an element $\omega$ of $\Omega$ is selected by "nature" at the start of time. This determines the initial condition $\boldx_0$ and the shocks $\boldw_t$ as

```{math}
:enumerated: false

\boldx_0(\omega), \boldw_1(\omega), \boldw_2(\omega),
    \boldw_3(\omega),\ldots
```

From these, each state vector $\boldx_t$ is determined via

```{math}
:enumerated: false

\boldx_t(\omega) = H_t(\boldx_0(\omega), \boldw_1(\omega),
        \boldw_2(\omega),\ldots, \boldw_t(\omega))
```

where $H_t$ is the function in {eq}`eq-xaft`.

#### The Density Case

If $Q(\bolds, \cdot)$ is absolutely continuous for all $\bolds$, then we can define the corresponding **stochastic density kernel** or the **transition density** via

$$
q(\bolds, \cdot) := \text{ the density of } Q(\bolds, \cdot)
    \text{ for all } \bolds \in S
$$ (eq-td)

In other words, $q(\bolds, \cdot)$ is the conditional density of $\boldx_{t+1}$ given $\boldx_t = \bolds$. Heuristically, $q(\bolds, \bolds') \diff \bolds'$ is the probability of transitioning from $\bolds$ to $\bolds'$ in one step.

```{prf:example}
:label: eg-archskd

Recall the ARCH(1) process, with stochastic kernel $Q$ as given in [example %s](#eg-archsk). The transition density $q(s, \cdot)$ is the density of $y = (\alpha_0 + \alpha_1 s^2)^{1/2} w_{t+1}$ when $\lL(w_{t+1}) = \nN(0, 1)$. Hence

$$
q(s, s') 
        = 
        \frac{1}{\sqrt{2 \pi \sigma_s^2}}
           \exp \left\{ - 
           \frac{(s')^2}{2\sigma^2_s} \right\} 
        \quad \text{where} \quad
        \sigma^2_s := \alpha_0 + \alpha_1 s^2
$$ (eq-arch_td)
```


A visualization of the ARCH kernel {eq}`eq-arch_td` is given in [figure %s](#f-arch_sk) when $\alpha_0 = 0.5$ and $\alpha_1 = 0.8$. The diagonal line is the 45 degree line. Each vertical fiber of the graph corresponds to one density $s' \mapsto q(s, s')$. For large $s$ (to the right on the horizontal axis), the corresponding density puts most mass below the 45 degree line, indicating a tendency to decrease. Conversely, when $-s$ is large, the process tends to increase.

```{figure} figures/arch_sk.svg
:name: f-arch_sk

Transition density $q(s, s')$ for the ARCH kernel {eq}`eq-arch_td`
```

```{prf:example}

Consider a general model with additive shock, so that

$$
\boldx_{t+1} = G(\boldx_t, \boldw_{t+1}) = g(\boldx_t) + \boldw_{t+1}
        \quad \text{with} \quad
        \{\boldw_t\} \iidsim \psi
$$ (eq-ampadd)
for some function $g$ and density $\psi$. In this case the transition density has the form

$$
q( \bolds, \bolds') = \psi(\bolds' - g(\bolds))
$$ (eq-cdncg)
In the scalar case, {eq}`eq-cdncg` follows from [example %s](#eg-ndat), since $q(s, \cdot)$ is the density of $y = g(s) + w_{t+1}$ with $\lL(w_{t+1}) =\psi$. This gives $q(s, s') = \psi(s' - g(s))$. The multivariate proof uses a similar argument.
```


### Marginal and Joint Distributions

Let's note some features of the marginal and joint distributions of a Markov process. Consider our generic process $\boldx_{t+1} = G(\boldx_t, \boldw_{t+1})$ from {eq}`eq-amp`, where $\lL(\boldx_0) = P_0$ and $\{\boldw_t\}$ is IID on $\RR^M$ with common distribution $\Psi$. The **marginal distribution of the state vector** is just $\lL(\boldx_t)$. Since $\boldx_t$ is a well-defined random vector (see in particular {eq}`eq-xaft`), $\lL(\boldx_t)$ is also well-defined, and we denote it in what follows by $P_t$. That is,

```{math}
:enumerated: false

P_t(B) = \PP\{\boldx_t \in B\}
    \qquad (B \in \bB(S))
```

````{prf:example}

Consider the Gaussian AR(1) process in {eq}`eq-sgar1`. Suppose that $\lL(x_0)= \nN(\mu_0, \sigma_0^2)$, where $\mu_0$ and $\sigma_0$ are given constants, and that $x_0$ is independent of the shock process $\{w_t\}$. Combining the moving average representation of $x_t$ in {eq}`eq-ar1ft` with [fact %s](#fa-lcinorm), we have

```{math}
:enumerated: false

P_t 
        = \nN 
            \left(
                \mu_t, \, \sigma_t^2
            \right)
        \quad \text{where} \quad
        \mu_t := b \sum_{k=0}^{t-1} a^k + a^t \mu_0
        \; \text{ and } \;
        \sigma^2_t := \sum_{k=0}^{t-1} a^{2k} c^2  + a^{2t} \sigma_0^2
```
````


{numref}`f-norm_den_seq` shows a visualization of this sequence when $\mu_0=0$, $\sigma_0 = b = c =1$ and $a = 0.8$. In fact we don't actually need to recompute the sums in the definition of $\mu_t$ and $\sigma_t^2$ at each step because these sequences satisfy the recursions

$$
\mu_{t+1} = b + a \mu_t
    \quad \text{and} \quad
    \sigma_{t+1}^2 = a^2 \sigma_t^2 + c^2
$$ (eq-msse)

You can verify this by taking expectations and variances on both side of {eq}`eq-sgar1`.

```{figure} figures/norm_den_seq.svg
:name: f-norm_den_seq

Sequence of marginal densities for the Gaussian AR(1) process
```

For a general Markov process, tracking evolution of the first two moments is usually more complicated, and even then fails to provide full information on the distributions. We do, however, have the following general result:

```{prf:property}
:label: fa-rdli

The marginal distributions of a Markov process with stochastic kernel $Q$ obey the recursion

$$
P_{t+1}(B) = \int Q(\bolds, B) P_t(\diff  \bolds)
        \qquad (B \in \bB(S),\; t \geq 0)
$$ (eq-rdli)
```


{prf:ref}`ex-pfr` asks you to prove this. [Fact %s](#fa-rdli) is a version of the law of total probability . To calculate the probability of reaching $B$ next period, we sum the probability of that event given outcome $\bolds$ today, weighted by the probability that $\bolds$ occurs today. Here's a parallel result for the density case.

```{prf:property}
:label: fa-detsmp

If $Q(\bolds, \cdot)$ is absolutely continuous for all $\bolds \in S$, then $P_t$ is absolutely continuous for all $t \geq 1$. Letting $p_t$ be the corresponding density, the sequence $\{p_t\}$ satisfies

$$
p_{t+1}(\bolds') = \int q( \bolds, \bolds') p_t(\bolds) \diff \bolds
        \qquad (B \in \bB(S),\; t \geq 1)
$$ (eq-rdli2)
where $q$ is the transition density corresponding to $Q$.
```


Let's now look at joint distributions. We'll cover these only in the density case, which is the most important case in econometrics.

```{prf:property}

Let $\{\boldx_t\}$ be a Markov process with transition density $q$ and initial density $p_0$. Under the conditions of [fact %s](#fa-detsmp), the joint density $p_T$ of $\boldx_0, \ldots, \boldx_T$ is

$$
p_T(\bolds_0, \ldots, \bolds_T) 
        = p_0(\bolds_0) \prod_{t=0}^{T-1} q(\bolds_t, \bolds_{t+1} )
$$ (eq-jftd)
```


This result can be obtained in a few different ways. One is to repeatedly apply {eq}`eq-decompjd` to factor the joint density into

```{math}
:enumerated: false

p_T(\bolds_0, \ldots, \bolds_T) 
    = p_0(\bolds_0) \prod_{t=0}^{T-1} p(\bolds_{t+1} \given \bolds_t, \ldots, \bolds_0)
```

This is the expression for the joint density of a general stochastic process, where the current state depends on lagged states of all orders. With the first-order Markov property, an expression like $p(\bolds_{t+1} \given \bolds_t, \ldots, \bolds_0)$ reduces to $p(\bolds_{t+1} \given \bolds_t)$, which is just $q(\bolds_t, \bolds_{t+1})$.

(ss-snry)=
### Stationarity of Markov Processes
Consider the notion of stationarity introduced in {ref}`ss-spdef`. Not every Markov process is stationary. As stated above, a necessary condition is that the marginals are constant, and in [figure %s](#f-norm_den_seq) we saw an example of how the marginals $\{P_t\}$ can change over time. {numref}`f-long_norm_den_seq` shows the same sequence but now longer, with identical parameters. Probability mass moves from left to right.

```{figure} figures/long_norm_den_seq.svg
:name: f-long_norm_den_seq

Convergence of marginal densities, Gaussian AR(1) case
```

While $\{P_t\}$ is not typically constant, for the particular sequence shown in [figure %s](#f-long_norm_den_seq), the difference between successive distributions is diminishing over time. We can see this analytically as well by recalling that $P_t = \nN(\mu_t, \sigma^2_t)$, where $\mu_t$ and $\sigma^2_t$ are as given in {eq}`eq-msse`. Provided that $|a| < 1$, these sequences converge, with

```{math}
:enumerated: false

\mu_t \to \mu_{\infty} := \frac{b}{1 - a}
    \quad \text{and} \quad
    \sigma^2_t \to \sigma^2_{\infty} := \frac{c^2}{1 - a^2}
```

Let $P_\infty$ denote the distribution associated with these limits. That is,

$$
P_\infty 
    = \nN(\mu_\infty, \sigma_\infty^2)
    = \nN \left[ \frac{b}{1 - a}, \frac{c^2}{1 - a^2} \right]
$$ (eq-sdsg)

What happens if we actually attain this limiting distribution, or if we take our initial draw $x_0$ from this distribution? On an intuitive level, it seems that the distribution sequence will be constant, remaining at $P_\infty$ over time. Indeed this is the case, as [exercise %s](#ex-inspi) asks you to confirm.

The limiting distribution is an example of what is called a stationary distribution. More generally, a **stationary distribution** for a Markov process on $S$ with stochastic kernel $Q$ is any distribution $P_\infty$ on $S$ satisfying

```{math}
:enumerated: false

P_\infty(B) = \int Q(\bolds, B) P_\infty(\diff \bolds)   
    \qquad \text{for all } B \in \bB(S)
```

Comparing with the recursion $P_{t+1}(B) = \int Q(\bolds, B) P_t(\diff \bolds)$ from [fact %s](#fa-rdli), we see that the definition of stationarity of $P_\infty$ means precisely that

```{math}
:enumerated: false

\lL(\boldx_t) = P_\infty
    \; \implies \; 
    \lL(\boldx_{t+1}) = P_\infty
```

From this we see that any process starting at $P_\infty$ will be identically distributed. In fact it is stationary, in the sense that {eq}`eq-statproc` holds. Let's record this:

```{prf:property}
:label: fa-idmp

Let $\{\boldx_t\}$ be a Markov process with stationary distribution $P_{\infty}$. If $\lL(\boldx_0) =  P_{\infty}$, then $\{\boldx_t\}$ is a stationary stochastic process, with $\lL(\boldx_t) = P_{\infty}$ for all $t$.
```


```{prf:example}
:label: eg-sgid

Consider the scalar Gaussian linear AR(1) process from {eq}`eq-sgar1`. If $|a| < 1$ and $\lL(x_0) = P_\infty$, where $P_\infty$ is as in {eq}`eq-sdsg`, then $\{x_t\}$ is stationary.
```


Moments of the stationary distribution are called the **stationary moments** of the Markov process. Other names for the stationary distribution include the "invariant\" or "unconditional\" distribution, and the moments are then called the invariant or unconditional moments.

```{prf:property}

If $Q(\bolds, \cdot)$ is absolutely continuous for all $\bolds \in S$ with transition density $q(\bolds, \cdot)$, then every stationary distribution $P_\infty$ is absolutely continuous, and its density $p_\infty$ satisfies

$$
p_{\infty}(\bolds') 
        = \int q(\bolds, \bolds' ) p_{\infty}(\bolds) \, \diff \bolds
        \quad \text{for all } \bolds' \in S
$$ (eq-dsdist)
```


Such a density is called a **stationary density**.

#### Existence

Not every Markov process has a stationary distribution. For example, in the linear scalar AR(1) model $x_{t+1} = b + a x_t + c w_{t+1}$, we know that the marginal variance evolves according to $\sigma^2_{t+1} = a^2 \sigma_t^2 + c^2$. If $|a| \geq 1$, then this sequence diverges. Because the variance is not constant, the distribution sequence cannot be constant either.

Let's state one sufficient condition for existence. The condition we give slightly specializes a famous result known as the Krylov–Bogolyubov theorem. It uses the notion of a **coercive** function, which is any nonnegative function $V$ on the state space $S$ such that

```{math}
:enumerated: false

C_\gamma := \setntn{\bolds \in S}{V(\bolds) \leq \gamma}
```

is a compact subset of $S$ for all $\gamma \in \RR$. Intuitively, coercive functions are functions that go to $+\infty$ as we move "away from the center\" of the state space.

```{prf:example}

If $S = \RR^K$ and $V(\bolds) = \| \bolds \|$, then $V$ is coercive on $S$. Indeed, for this function, each $C_\gamma$ is a closed sphere centered on the origin. Closed and bounded subsets of $\RR^K$ are compact.
```


```{prf:example}
:label: eg-sicoer

If $S = \RR$ and $V(s) = s^2$, then $V$ is coercive on $S$. To see this, just observe that $C_\gamma = [-\sqrt{\gamma}, \sqrt{\gamma}]$. This is a closed, bounded subset of $S$.
```


Now we can state the following:

```{prf:theorem}
:label: t-esd

Let $\boldx_{t+1} = G(\boldx_t, \boldw_{t+1})$ be the process on $S$ defined in {eq}`eq-amp`. If

(i) $\bolds \mapsto G(\bolds, \boldw)$ is continuous for each fixed $\boldw \in \RR^M$, and

(ii) there exists a coercive function $V$ on $S$ and constants $\lambda$ and $L$ such that $0 \leq \lambda < 1$ and

$$
\EE  V[ G(\bolds, \boldw_{t+1} ) ] 
                    \leq \lambda V(\bolds) + L
                    \qquad (\bolds \in S)
$$ (eq-cwrtc)

then the Markov process has at least one stationary distribution.
```


See §11.2 of {cite}`stachurski2009economic` for proofs and references. Condition {eq}`eq-cwrtc` is an example of a **drift condition**. One way to understand it is to write it as

```{math}
:enumerated: false

\frac{\EE  V[ G(\bolds, \boldw_{t+1} )  ] }{V(\bolds)}
        \leq \lambda + \frac{L}{V(\bolds) }
```

As we move "away from the center\" of the state space, $V(\bolds)$ becomes large, and the term on the right is eventually $< 1$. Looking at the left hand side, this tells us that the value attached to the state by $V$ is expected to decrease, which in turn suggests movement away from the "edges\" of the state space and back towards the center.

````{prf:example}
:label: eg-archas

Consider again the ARCH(1) process $x_{t+1} = (\alpha_0 + \alpha_1 x_t^2)^{1/2} w_{t+1}$ from {eq}`eq-arch`, where, as usual $\alpha_0 > 0$ and $\alpha_1 \geq 0$. Let $S = \RR$. If $\alpha_1 < 1$, then this process has a stationary distribution on $S$. To see this, observe that

```{math}
:enumerated: false

G(s, w) = (\alpha_0 + \alpha_1 s^2)^{1/2} w
        \qquad (s, w \in \RR)
```

is continuous in $s$ for each fixed $w$. Moreover $V(s) = s^2$ is coercive on $S$ by [example %s](#eg-sicoer), and

$$
\EE G(s, w_{t+1})^2  
        = (\alpha_0 + \alpha_1 s^2) \EE w_{t+1}^2
        = \alpha_0 + \alpha_1 s^2
$$ (eq-archas)
Hence {eq}`eq-cwrtc` holds with equality when $V(s) :=s^2$, $\lambda :=\alpha_1$ and $L := \alpha_0$.
````


(ss-larlln)=
### Asymptotics of Markov Processes
When do Markov processes exhibit nice properties such as ergodicity? To gain some intuition, let's first consider the process $x_{t+1} = a x_t + w_{t+1}$, where $\{w_t\}$ is IID $\nN(0, 1)$ and $|a| < 1$. We consider a stationary process, where each $x_t$ has marginal

$$
P_{\infty} := \nN \left( 0, \, (1- a^2)^{-1} \right)
$$ (eq-ar1statd)

A law of large numbers for this process takes the form

$$
\bar x_T \toprob \int s\, P_{\infty}(\diff s) = 0
    \quad \text{as} \quad T \to \infty
$$ (eq-ar1llnw)

[Fact %s](#fa-ffci) tells us that for {eq}`eq-ar1llnw` to hold it suffices that

$$
\EE \bar x_T \to 0
    \quad \text{and} \quad
    \var \bar x_T \to 0
    \quad \text{as} \quad T \to \infty
$$ (eq-sclln)

The first part is obvious, since $\EE \bar x_T = \frac{1}{T} \sum_{t=1}^T \EE x_t$ and $\EE x_t = 0$. For the second part, since $\EE \bar x_T = 0$,

```{math}
:enumerated: false

\var \left[ \frac{1}{T} \sum_{t=1}^T x_t \right]
    =  \EE \left[ \frac{1}{T} \sum_{t=1}^T x_t \right]^2
    =  \frac{1}{T^2} \sum_{n=1}^T \sum_{m=1}^T \EE[x_n x_m]
```

Using $\EE[x_n x_m] = \EE[x_m x_n]$ we obtain the bound

$$
\var \bar x_T  
    \leq \frac{2}{T^2} \sum_{1 \leq n \leq m \leq T} \EE[x_n x_m] 
     = \frac{2}{T^2} \sum_{j=0}^{T-1} \sum_{t=1}^{T-j} \EE[x_{t+j} x_t]
$$ (eq-bvcp)

Thus, to obtain $\var \bar x_T \to 0$ as $T \to \infty$, we require the covariance terms $\EE[x_{t+j} x_t]$ to converge to zero in $j$ sufficiently quickly. That they do in this particular instance is the subject of [exercise %s](#ex-cbv).

To summarize, if correlations die out sufficiently quickly, then we can regain the LLN. The same is true for the CLT when second moment conditions hold. We now state a much more general result for sample paths of Markov processes. The proofs are omitted but the essential ideas are the same.

````{prf:theorem}
:label: t-umspt

Let $\{\boldx_t\}$ be a Markov process on $S$ with stochastic kernel $Q$. Suppose that $Q(\bolds, \cdot)$ is absolutely continuous for all $\bolds \in S$ and let $q(\bolds, \cdot )$ be the corresponding transition density. If

(a) $q$ is strictly positive and continuous on $S \times S$, and

(b) there exists a coercive function $V$ on $S$ and constants $\lambda$ and $L$ such that $0 \leq \lambda < 1$ and

```{math}
:enumerated: false

\int V(\bolds') q(\bolds, \bolds' ) \diff  \bolds'
                \leq \lambda V(\bolds) + L
                \qquad (\bolds \in S)
```

then the following statements are true:

(i) $Q$ has a unique stationary distribution $P_\infty$ with density $p_\infty$.

(ii) $P_t \tow P_\infty$ as $t \to \infty$ for all initial $P_0$.

(iii) If $h \colon S \to \RR$ is $\bB$-measurable and $\int |h(\bolds)| p_{\infty}(\bolds)\diff \bolds < \infty$, then

$$
\frac{1}{T} \sum_{t=1}^T h(\boldx_t) 
                    \toprob 
                    \int h(\bolds) p_{\infty}(\bolds) \diff \bolds
                    \quad \text{as} \quad
                    T \to \infty
$$ (eq-mllnr)

(iv) If, in addition, $h^2 \leq V$, then there exists a $\sigma^2_h \geq 0$ such that

$$
\sqrt{T}
                    \left\{
                        \frac{1}{T} \sum_{t=1}^T h(\boldx_t) -
                        \int h(\bolds) p_{\infty}(\bolds) \diff \bolds
                    \right\}
                    \tod \nN(0, \sigma^2_h)
                    \quad \text{as} \quad
                    T \to \infty
$$ (eq-mkclt)
````


Conditions (a) and (b) are not dissimilar to the conditions for existence in [theorem %s](#t-esd). Condition (b) is the drift condition in {eq}`eq-cwrtc` specialized to the density setting. The major difference is that $q$ is required to be strictly positive. This gives enough mixing to generate uniqueness and the sample path properties. {prf:ref}`t-umspt` can be proved from results in {cite}`meyn2009markov`.[^1]

{prf:ref}`t-umspt` is but one of a vast range of results on asymptotic properties of Markov processes. Discussion of alternatives is given in {ref}`s-frsp`. In fact the conditions are rather strict and can certainly be weakened at the expense of a more complex theorem statement. Nonetheless, they suffice for many applications of interest.

The constant $\sigma^2_h$ in {eq}`eq-mkclt` is given by the expression

$$
\sigma^2_h 
    := \var h(\boldx_0)
    + 2 \sum_{t=1}^{\infty} \cov [h(\boldx_0), h(\boldx_t)]
$$ (eq-mavr)

with the understanding that $\{\boldx_t\}$ is a stationary process. If $\sigma^2_h = 0$, then the limiting distribution is degenerate.

```{prf:property}
:label: fa-gsmp

Let $\boldx_{t+1} = g(\boldx_t) + \boldw_{t+1}$ be the additive shock process on $\RR^K$ from {eq}`eq-ampadd`. If $\EE \|\boldw_t\|$ is finite, $g$ is continuous, the density $\psi$ of $\boldw_t$ is continuous and positive everywhere on $\RR^K$, and there exist constants $\lambda$ and $L$ such that $0 \leq \lambda < 1$ and

$$
\| g(\bolds) \| \leq \lambda \| \bolds \| + L
        \quad \text{ for all } \bolds \in \RR^K
$$ (eq-drift)
then the conditions of [theorem %s](#t-umspt) are satisfied.
```


To see why, recall from {eq}`eq-cdncg` that the transition density satisfies $q(\bolds, \bolds' ) = \psi(\bolds' - g(\bolds))$. Since $\psi$ is assumed continuous and positive on $\RR^K$, the transition density has these properties on $\RR^K \times \RR^K$. Hence condition (a) in [theorem %s](#t-umspt) is satisfied. Condition (b) is also satisfied because, taking $\lambda$ and $L$ as in {eq}`eq-drift`, we have

```{math}
:enumerated: false

\int \|\bolds'\| q(\bolds, \bolds' ) \diff  \bolds'
    = \EE \|g(\bolds) + \boldw_{t+1} \| 
    \leq \lambda \| \bolds \| + L + \int \|\boldw\| \psi(\boldw) \diff \boldw
```

```{prf:example}
:label: eg-ar1ssc

Let $\{x_t\}$ be a scalar AR(1) process $x_{t+1} = a x_t + b + c w_{t+1}$, where $\{w_t\}$ is IID and standard normal. If $|a| < 1$ and $c > 0$, then the conditions of [fact %s](#fa-gsmp) are satisfied. The density of $cw_t$ is $\nN(0, c^2)$, which is continuous and strictly positive on $\RR$. To see that {eq}`eq-drift` holds, rewrite the model as $x_{t+1} = g(x_t) + cw_{t+1}$, where $g(s) := a s + b$. By the triangle inequality, $|g(s)| \leq \lambda |s| + L$, where $\lambda := |a|$ and $L := |b|$.
```


````{prf:example}
:label: eg-ztar

Here's a variation on the scalar "threshold autoregressive" model:

```{math}
:enumerated: false

x_{t+1} = a |x_t| + (1 - a^2)^{1/2} w_{t+1}  
        \quad \text{with} \quad
        |a| < 1
        \quad \text{and} \quad
        \{w_t\} \iidsim \nN(0, 1)
```

The conditions of [fact %s](#fa-gsmp) are satisfied. To see this, we can rewrite the model as

```{math}
:enumerated: false

x_{t+1} = g(x_t) + v_{t+1}, \quad
        g(s) = a |s| 
        \quad \text{and} \quad
        \{v_t\} \iidsim \nN(0, 1 - a^2)
```

Clearly {eq}`eq-drift` is satisfied with $\lambda = |a|$ and $L = 0$. Moreover $\lL(v_t) = \nN(0, 1 - a^2)$, which has finite mean and a bounded density everywhere positive on $\RR$.
````


For most Markov processes the stationary density has no known closed-form. However, in the case of [example %s](#eg-ztar), the stationary density is known to have the form $p_{\infty}(s) = 2 \phi(s) \Phi(q s)$, where $q := a (1 - a^2)^{-1/2}$, $\phi$ is the standard normal density and $\Phi$ is the standard normal CDF.

````{prf:example}

Consider the STAR model from {eq}`eq-starrs`, where the function $g$ is given by $g(s) := (b_0 + a_0 s)(1-\tau(s)) + (b_1 + a_1 s)\tau(s)$ for some $\tau \colon \RR \to [0,1]$. Applying the triangle inequality and using $0 \leq \tau(s) \leq 1$, we have

```{math}
:enumerated: false

|g(s)| 
      \leq |b_0| + |b_1| + \max \{ |a_0|, |a_1|\} |s|
```

Letting $L := |b_0| + |b_1|$ and $\lambda := \max\{|a_0|, |a_1|\}$, we then have $|g(s)| \leq  \lambda |s| + L$. If both $|a_0|$ and $|a_1|$ are strictly less than one, then condition {eq}`eq-drift` is satisfied.
````


```{prf:example}

Recall the ARCH process $x_{t+1} = (\alpha_0 + \alpha_1 x_t^2)^{1/2} w_{t+1}$ from {eq}`eq-arch`, where $\{ w_t \}$ is IID and standard normal. In [example %s](#eg-archas) we saw that this model has a stationary distribution on $S=\RR$ whenever $\alpha_0 > 0$ and $0 \leq \alpha_1 < 1$. In particular, we showed in {eq}`eq-archas` that it satisfies drift with respect to the coercive function $V(s) = s^2$. This gives condition (b) of [theorem %s](#t-umspt). Condition (a) also holds because the transition density is continuous and everywhere positive.
```


### The Linear Case

Recall the VAR process in {eq}`eq-var1`, which, to repeat, takes the form

$$
\boldx_{t+1} = \boldb + \boldA \boldx_t + \boldC \boldw_{t+1}
$$ (eq-var12)

where $\{ \boldw_t \}$ is assumed to be IID and to satisfy $\EE \boldw_t = \boldzero$ and $\EE[\boldw_t\boldw_t'] = \boldI$. Regarding the asymptotics of this process, if the shock vector has a continuous and positive density, then we can apply [fact %s](#fa-gsmp) to obtain asymptotic stability and LLN and CLT results whenever the matrix norm of $\boldA$ satisfies $\|\boldA \| < 1$. To see this, let $g(\bolds) = \boldb + \boldA \bolds$. By the norm triangle inequality and the definition of the matrix norm, we have

```{math}
:enumerated: false

\|g(\bolds)\| 
    = \| \boldb + \boldA \bolds\|
    \leq \| \boldb \|  + \| \boldA \bolds\|
    \leq \| \boldb \|  + \| \boldA \| \| \bolds\|
```

Hence {eq}`eq-drift` holds.

As mentioned above, the stability and sample path conditions considered above are relatively strict. In the linear case, we can certainly do better. Here's one such result:

```{prf:property}

For the VAR process {eq}`eq-var12`, the conclusions of [theorem %s](#t-umspt) hold with $V = \| \cdot \|$ whenever $\rho(\boldA) < 1$.
```


Here $\rho(\boldA)$ is the spectral radius of $\boldA$, as defined. Under this condition, a unique stationary distribution exists. Taking expectations across {eq}`eq-var12`, we see that its mean $\boldmu_{\infty}$ must satisfy $\boldmu_{\infty} = \boldb + \boldA \boldmu_{\infty}$. Applying the Neumann series lemma, this equation has the unique solution

```{math}
:enumerated: false

\boldmu_\infty 
    = (\boldI - \boldA)^{-1} \boldb 
    = \sum_{i=0}^{\infty} \boldA^i  \boldb
```

Let $\boldSigma_\infty$ be the asymptotic variance (i.e., the variance–covariance matrix of the stationary distribution). From {eq}`eq-var12` it must satisfy

```{math}
:enumerated: false

\boldSigma_\infty = \boldA \boldSigma_\infty \boldA^\T + \boldC \boldC^\T
```

This is an example of a **discrete time Lyapunov equation**. Under the condition $\rho(\boldA) < 1$ it can be solved by iteration or by using existing numerical solvers. Here's an example in Julia that uses a solver from the QuantEcon package:

``` julia
using QuantEcon

A = [0.8 -0.2; 
    -0.1 0.7]
C = [0.5 0.4;
     0.4 0.6]

solve_discrete_lyapunov(A, C * C')
```

## Martingales

Martingales generalize the notion of a random walk, where increments to the current state are draws from a white noise process. The concept of martingales has been central to many advances in probability since the mid twentieth century. Econometricians also like martingales because they neatly capture the idea of "unpredictable\" change.

### Definitions

To define martingales, we need the notion of a filtration, which is an increasing sequence of information sets. Recall from {ref}`ss-m` that an information set is just a set of random variables. Let $\{\fF_t\}$ be a sequence of information sets. That is, $\fF_t$ is an information set for each $t$. The sequence $\{\fF_t\}$ is called a **filtration** if, in addition, $\fF_t \subset \fF_{t+1}$ for all $t$. Intuitively, expanding information sets reflect the idea that more information is revealed over time.

```{prf:example}

Let $\{\boldx_t\}$ be a sequence of random vectors, and let

$$
\fF_0 := \emptyset
        \quad \text{and} \quad
        \fF_t := \{\boldx_1, \ldots, \boldx_t\}
        \quad \text{for all } t \geq 1
$$ (eq-canfil)
The sequence $\{\fF_t\}$ is a filtration. It is called the **filtration generated by** $\{\boldx_t\}$.
```


Now let $\{m_t\}$ be a scalar stochastic process and let $\{\fF_t\}$ be a filtration. We say that $\{m_t\}$ is **adapted** to $\{\fF_t\}$ if $m_t$ is $\fF_t$-measurable for every $t$. (For the definition of measurability, see {ref}`ss-m`.)

The intuitive meaning is this: $\fF_t$ represents the variables we know as of time $t$. If $\{m_t\}$ is adapted to $\fF_t$, then we can compute $m_t$ at time $t$ as well.

```{prf:example}
:label: eg-smiad

Let $\{x_t\}$ be any stochastic process, let $\{\fF_t\}$ be the filtration generated by $\{x_t\}$ and let $\bar x_t := \frac{1}{t} \sum_{j=1}^t x_j$ be the sample mean up until date $t$. Then $\{\bar x_t\}$ is adapted to $\{\fF_t\}$.
```


```{prf:example}
:label: eg-notad

Let $\{w_t\}$ be IID standard normal, let $\{\fF_t\}$ be the filtration generated by $\{w_t\}$ and let $x_t := \sup_{j \geq t} w_j$. This process is not adapted to $\{\fF_t\}$ because the value $x_t$ cannot be determined from $\{w_1, \ldots, w_t\}$.
```


```{prf:property}
:label: fa-emfv

If $\{m_t\}$ is adapted to $\{\fF_t\}$, then $\EE[m_t \given \fF_{t+j}] = m_t$ for any $j \geq 0$.
```


The proof goes like this: By adaptedness, $m_t$ is $\fF_t$-measurable. From the definition of a filtration we have $\fF_t \subset \fF_{t+j}$. Applying [fact %s](#fa-mwrtgh), it follows that $m_t$ is $\fF_{t+j}$-measurable. [Fact %s](#fa-emfv) now follows from [fact %s](#fa-coce).

A sequence of integrable random variables $\{m_t\}$ adapted to a filtration $\{\fF_t\}$ is called a **martingale** with respect to $\{\fF_t\}$ if

```{math}
:enumerated: false

\EE [ m_{t+1} \given \fF_t] = m_t 
    \quad \text{for all} \quad t
```

The expression $\EE [ m_{t+1} \given \fF_t]$ is sometimes simplified to $\EE_t [ m_{t+1}]$ when the filtration is understood.

````{prf:example}
:label: eg-rw

Let $\{\eta_t\}$ be an IID sequence of zero-mean random variables and let $m_t := \sum_{j=1}^t \eta_j$. For example, $\eta_t$ might be the payoff on the $t$th round of a fair game and $m_t$ is the wealth of a gambler after the $t$th round.[^2] We can also express $\{m_t\}$ recursively via $m_{t+1} = m_t + \eta_{t+1}$ with $m_0 = 0$. The process $\{m_t\}$ is a martingale with respect to the filtration $\{\fF_t\}$ generated by $\{\eta_t\}$. That $\{m_t\}$ is adapted to $\{\fF_t\}$ follows immediately from the definitions. Moreover

```{math}
:enumerated: false

\EE[ m_{t+1} \given \fF_t]  
         =  \EE[  m_t + \eta_{t+1} \given \fF_t] 
         =  m_t + \EE[ \eta_{t+1} \given \fF_t]
```

where the last equality follows from linearity and [fact %s](#fa-emfv). Since $\{\eta_t\}$ is IID, we know that $\eta_{t+1}$ is independent of the variables in $\fF_t$, and hence $\EE[ \eta_{t+1} \given \fF_t] = \EE[ \eta_{t+1}] = 0$. Therefore $\{m_t\}$ is a martingale with respect to $\{\fF_t\}$ as claimed.
````


````{prf:example}

A well-known example of a martingale from economic theory is the consumption process in {cite}`hall1978stochastic`. To understand the idea, consider the Euler equation

```{math}
:enumerated: false

u'(c_t) = \EE_t \left[ \frac{1+r_{t+1}}{1+\rho} u'(c_{t+1}) \right]
```

where $u'$ is the derivative of a utility function $u$, $r_t$ is an interest rate and $\rho$ is a discount factor. The expectation $\EE_t[\cdot]$ can be thought of as $\EE[\cdot \given \fF_t]$, where $\fF_t$ contains all variables observable at time $t$. Specializing to the case $r_{t+1} = \rho$ and $u(c) = c - a c^2 / 2$, the Euler equation reduces to

```{math}
:enumerated: false

c_t 
        = \EE_t [c_{t+1} ]
        =: \EE [c_{t+1} \given \fF_t ]
```

Thus consumption is a martingale with respect to $\{\fF_t\}$.
````


#### Martingale Difference Sequences

A stochastic process $\{d_t\}$ is called a **martingale difference sequence** (or **MDS**) with respect to some filtration $\{\fF_t\}$ if

```{math}
:enumerated: false

\EE [ d_{t+1} \given \fF_t] = 0 \quad \text{for all} \quad t
```

An MDS is so named because if $\{m_t\}$ is a martingale with respect to $\{\fF_t\}$, then $d_t = m_t - m_{t-1}$ is an MDS with respect to $\{\fF_t\}$. See [exercise %s](#ex-dmds).

```{prf:property}
:label: fa-ummd

If $\{d_t\}$ is an MDS, then $\EE d_t=0$ for all $t$.
```


Indeed, by the law of iterated expectations, $\EE d_t = \EE[ \EE[ d_t \given \fF_{t-1} ]] = \EE 0 = 0$.

### Martingale Difference LLN and CLT

Next we consider asymptotics for martingale difference sequences. Such sequences are good candidates for the LLN and CLT. To see this, recall that integrable IID sequences always satisfy the LLN. With finite second moment they also satisfy the CLT. Suppose now that $\{m_t\}$ is a stationary MDS with respect to some filtration $\{\fF_t\}$. We don't have independence in this case but we do have the following:

```{prf:property}
:label: fa-lcmds

If $\{m_t\}$ is an MDS, then $\cov[m_{t+j}, m_t] = 0$ for all $j \geq 1$.
```


To see this, fix $t \geq 0$ and $j \geq 1$. We have

```{math}
:enumerated: false

\cov[m_{t+j}, m_t] 
    = \EE[m_{t+j} m_t] 
    = \EE[ \EE[ m_{t+j} m_t \given \fF_{t+j-1}]]
```

Since $t + j - 1 \geq t$ and $\{\fF_t\}$ is a filtration, we know that $m_t$ is $\fF_{t+j-1}$-measurable. Hence

```{math}
:enumerated: false

\EE[ \EE[ m_{t+j} m_t \given \fF_{t+j-1}]]  
    = \EE[ m_t \EE[ m_{t+j} \given \fF_{t+j-1}]] 
    = \EE[ m_t \cdot 0 ]
    = 0
```

The lack of correlation in [fact %s](#fa-lcmds) is a key part of the next result.

```{prf:theorem}
:label: t-mdclt

If $\{m_t\}$ is a stationary MDS with respect to some filtration $\{\fF_t\}$, then

$$
\frac{1}{T} \sum_{t=1}^T m_t \toprob 0
        \quad \text{as} \quad
        T \to \infty
$$ (eq-mdlln)
If, in addition, $\gamma^2 := \EE[m_t^2] < \infty$ and $\frac{1}{T} \sum_{t=1}^T \EE[m_t^2 \given \fF_{t-1} ] \toprob \gamma^2$ as $T \to \infty$, then

$$
T^{-1/2}  \sum_{t=1}^T m_t
        \tod \nN(0, \gamma^2)
        \quad \text{as} \quad
        T \to \infty
$$ (eq-mdclt)
```


The proof of {eq}`eq-mdlln` is identical to that of the classical LLN in [theorem %s](#t-lln0). {prf:ref}`t-mdclt` can be established from a martingale CLT proved in {cite}`durrett2010probability`.

(s-sim)=
## Simulation
Simulation is invaluable as a tool for exploring statistical concepts and solving applied problems. For example, the current popularity of Bayesian methods for econometrics and machine learning has been driven to a large extent by clever simulation algorithms that have proved successful at solving hard estimation problems.

Modern scientific programming environments all have good routines for simulating independent observations from a range of common named distributions. We've seen examples of this in the displayed computer code. But it's easy to run into situations where the distribution we are working with is nonstandard, or has no closed-form solution. We then need to roll our own routines for producing random variates. Let's review two of the most important algorithms.

(ss-it)=
### Inverse Transforms
Suppose that we wish to draw from a distribution on $\RR$ with CDF $F$. Assuming that we can generate uniform random variables, we can do so via the **inverse transform** method:

- draw $u$ from the uniform distribution on $[0, 1]$
- return $F^{-1}(u)$, where $F^{-1}$ is the quantile function of $F$ (see {eq}`eq-quantwi2`)

We wish to show that this technique is valid, or, equivalently, that

$$
\lL(u) = U[0, 1] 
    \quad \implies \quad
    \lL[F^{-1}(u)] = F
$$ (eq-itmp)

In the proof we'll assume that $F$ is strictly increasing. (This isn't necessary, but the proof of the general case provides no additional intuition.) Let $\lL(u) = U[0, 1]$ and let $z := F^{-1}(u)$. To see that $\lL(z) = F$, observe that, by strict monotonicity, $F(a) \leq F(b)$ if and only if $a \leq b$. As a result, for any $s \in \RR$,

```{math}
:enumerated: false

\PP\{z \leq s\}
    = \PP\{F^{-1}(u) \leq s\}
    = \PP\{u \leq F(s)\}
```

The defining property of a $U[0, 1]$ random variable $u$ is that $\PP\{u \leq t\} = t$ whenever $0 \leq t \leq 1$. Using this fact in the preceding calculation gives $\PP\{z \leq s\} = F(s)$, proving {eq}`eq-itmp`.

Here's one extension for the basic inverse transform method. Suppose that we want to simulate from a joint distribution on $\RR^N$ that is defined in terms of a copula (see {ref}`ss-cop`). In particular, suppose that $F_n$ is a continuous CDF on $\RR$ for $n=1,\ldots, N$, that $C$ is a copula on $\RR^N$ and that the joint distribution of interest is

```{math}
:enumerated: false

F(s_1, \ldots, s_N) = C(F_1(s_1), \ldots, F_N(s_N))
```

(This is just {eq}`eq-ffcop`.) If we can simulate from $C$, then we can simulate from $F$. The method is

- draw $u_1, \ldots, u_N$ from $C$
- return $(x_1, \ldots, x_N) := (F^{-1}_1 (u_1), \ldots, F^{-1}_N(u_N))$

{prf:ref}`ex-simcop` asks you to prove that $(x_1, \ldots, x_N)$ is a draw from $F$.

(ss-mcmc)=
### Markov Chain Monte Carlo
Markov chain Monte Carlo (MCMC) is a way of simulating from a given density $\pi$ on $S \subset \RR^K$. The idea is to construct a stochastic kernel $P$ on $S$ such that

(i) $\pi$ is a stationary distribution for $P$, and

(ii) $P$ is sufficiently ergodic that its sample path averages converge to expectations under $\pi$

See {eq}`eq-mllnr` for a clearer statement of (ii). If we specialize to indicator functions we can approximate probabilities under $\pi$ too.

There are two major strands of MCMC algorithms: Metropolis–Hastings and the Gibbs sampler. Here we'll cover the former. While boldface symbols are omitted in what follows, the ideas are valid in any dimension. In the section on Bayesian estimation in {ref}`ss-be`, I will explain why Metropolis–Hastings fits so beautifully with Bayesian updating.

The **Metropolis–Hastings algorithm** starts with a Markov process in the form of a stochastic density kernel $q = q(s, s')$ called the **proposal density**. Draws from the proposal density are called **proposals**. Each time we draw a proposal, we either (i) accept it and move to that new state, or (ii) reject it and stay put. The probability of accepting is structured so that the chain tends to stay in regions where $\pi$ puts most probability mass. As a consequence, for a sequence $\{x_t\}$ generated by this process, we get

```{math}
:enumerated: false

\text{fraction of time spent in } B = 
    \;\;
    \frac{1}{T} \sum_{t=1}^T \1\{x_t \in B\} \approx \pi(B)
    \quad \text{for large } T
```

This is a version of {eq}`eq-mllnr` with $h = \1_B$.

The ideas are subtle, so let's build up in steps. Consider first a general accept-reject style process of the kind discussed above, where the probability of acceptance is determined by an arbitrary function $\alpha = \alpha(x_t, y)$ of the current state $x_t$ and proposal $y$. The function takes values in $[0, 1]$. The algorithm for generating $x_{t+1}$ from $x_t$ given $\alpha$ and proposal density $q$ is given in [algorithm %s](#alg-qa).

```{prf:algorithm} A Markov process with proposal density $q$ and acceptance function $\alpha$
:label: alg-qa

- draw $y$ from $q(x_t, \cdot)$
- draw $u$ independently from the uniform $[0, 1]$ distribution
- **if** $u \leq \alpha(x_t, y)$:
  - (– with probability $\alpha(x_t, y)$)
  - set $x_{t+1} = y$
- **else**:
  - (– with probability $1 - \alpha(x_t, y)$)
  - set $x_{t+1} = x_t$
```

The stochastic kernel $P$ associated with $\{x_t\}$ in [algorithm %s](#alg-qa) has the form

$$
P(s, B) = \int_B p(s, s') \diff s' + (1 - \lambda(s)) \1\{s \in B\}
$$ (eq-mhk)

for all $s \in S$ and $B \in \bB(S)$, where

```{math}
:enumerated: false

p(s, s') := q(s, s') \alpha(s, s')
    \quad \text{and} \quad
    \lambda(s) := \int p(s, s') \diff s'
```

To show this, we can break the probabilities up as follows. Fix any Borel set $B$ and let $x_t$ be given at $s$. Let $x_{t+1}$ be drawn according to [algorithm %s](#alg-qa) when $x_t = s$. By additivity of $\PP$ we can write $\PP\{x_{t+1} \in B\}$ as

$$
\PP\{x_{t+1} \in B\} \cap \{u \leq \alpha(s, y)\}
        + \PP\{x_{t+1} \in B\} \cap \{u > \alpha(s, y)\}
$$ (eq-bprob)

Our aim is to show that the two terms in {eq}`eq-bprob` match the two terms on the right-hand side of {eq}`eq-mhk`.

Using the fact that intersections become products when we switch to indicators ([fact %s](#fa-cbif), page ), we can applying the law of iterated expectations to obtain

```{math}
:enumerated: false

\begin{aligned}
\PP\{x_{t+1} \in B\} \cap \{u \leq \alpha(s, y)\}
    & = \PP\{y \in B\} \cap \{u \leq \alpha(s, y)\}
    \\
    & = \EE \left[ \1\{y \in B\} \1 \{u \leq \alpha(s, y)\} \right]
    \\
    & = \EE \left\{
        \EE \left[ \1\{y \in B\} \1 \{u \leq \alpha(s, y)\} \given y \right]
          \right\}
    \\
    & = \EE \left\{  
            \1\{y \in B\}
            \EE \left[  \1 \{u \leq \alpha(s, y)\} \given y \right]
          \right\}
\end{aligned}
```

The term $\EE \left[  \1 \{u \leq \alpha(s, y)\} \given y \right]$ is just $\PP\{u \leq \alpha(s, y)\}$ with $y$ given, which, since $u$ is uniform, evaluates to $\alpha(s, y)$. Hence the final expression in the last display evaluates to

```{math}
:enumerated: false

\int \1\{s' \in B\} \alpha(s, s') q(s, s') \diff s'
        = \int_B \alpha(s, s') q(s, s') \diff s'
```

$$
\fore
    \PP\{x_{t+1} \in B\} \cap \{u \leq \alpha(s, y)\}
         = \int_B p(s, s') \diff s'
$$ (eq-ppxb)

Equation {eq}`eq-ppxb` says that the first term in {eq}`eq-bprob` is equal to the first term on the right-hand side of {eq}`eq-mhk`. For the second terms, note that $x_{t+1} = x_t = s$ when $u > \alpha(s, y)$, so

```{math}
:enumerated: false

\begin{aligned}
\PP \{x_{t+1} \in B\} \cap \{u > \alpha(s, y)\}
    & = \EE \left[ \1 \{s \in B\} \1\{u > \alpha(s, y)\} \right]
    \\
    & = \1\{s \in B\} \EE [ \1 \{u > \alpha(s, y)\} ]
    \\
    & = \1\{s \in B\} [1 - \PP \{u \leq \alpha(s, y)\} ]
\end{aligned}
```

To match this expression with the second term in {eq}`eq-mhk`, we just need to show that $\PP \{u \leq \alpha(s, y)\} = \int p(s, s') \diff s'$. But this is immediate from {eq}`eq-ppxb` if we take $B = S$.

So far we have verified that the kernel $P$ in {eq}`eq-mhk` is the stochastic kernel corresponding to [algorithm %s](#alg-qa). Now let's double back and design the selection function $\alpha$ so that the stationary distribution of $P$ is the target density $\pi$. In the Metropolis–Hastings algorithm, this is done by setting

$$
\alpha(s, s') 
    := 
        \min \left\{ 
            \frac{\pi(s') q(s', s)}{\pi(s) q(s,s')}
            , 1
             \right\}
$$ (eq-alpha)

If $\pi(s) q(s,s') = 0$, then $\alpha(s, s') := 1$. The min is in place so $\alpha(s, s')$ can be interpreted as a probability. We don't actually need to impose it in [algorithm %s](#alg-qa), since for $u \in [0, 1]$ we have $u \leq \min\{a, 1\}$ iff $u \leq a$. The important element is the ratio inside the min, which is large when the new proposal $s'$ has high probability under $\pi$ relative to the current location $s$ and small when the reverse is true. This keeps the chain in regions where $\pi$ has high mass, and leads to our next result.

```{prf:theorem}
:label: t-mh

Let $P$ be the stochastic kernel in {eq}`eq-mhk`. If $\alpha$ is as defined {eq}`eq-alpha`, then $\pi$ is a stationary distribution for $P$.
```


````{prf:proof}
:nonumber:

The proof turns on the fact that for this choice of $\alpha$ we have $p(s, s') \pi(s) = p(s', s) \pi(s')$ for all $s, s'$ (see [ex. %s](#ex-sta)). From that observation and the definition of $P$ we obtain

```{math}
:enumerated: false

\begin{aligned}
\int P(s, B) \pi(s) \diff s
        & = \int \left\{ 
            \int_B p(s, s') \diff s'
            \right\}
                 \pi(s) \diff s
                 + \int (1 - \lambda(s)) \1\{s \in B\} \pi(s) \diff s
            \\
            & = \int_B \left\{ 
                 \int p(s, s') \pi(s) \diff s
             \right\}
                 \diff s'
                 + \int_B (1 - \lambda(s)) \pi(s) \diff s
            \\
            & = \int_B \left\{ 
                \int p(s', s) \pi(s') \diff s
            \right\}
                 \diff s'
                 + \int_B (1 - \lambda(s)) \pi(s) \diff s
            \\
        & = \int_B \lambda(s') \pi(s') \diff s'
                 + \int_B (1 - \lambda(s)) \pi(s) \diff s
            \\
        & = \int_B \pi(s) \diff s
\end{aligned}
```

Thus $\pi$ is a stationary distribution of $P$, as was to be shown. ◻
````


It remains to establish that $P$ has the kind of stability properties discussed in [theorem %s](#t-umspt) (and in particular the LLN result in {eq}`eq-mllnr` of that theorem). This isn't universally true, and requires at least some conditions on the primitives $q$ and $\pi$. One well-known result is given in corollary 3 of {cite}`tierney1994markov` for the case where $q$ and $\pi$ are continuous and strictly positive on $S$ and $S$ itself is compact. In essence, the continuity and positivity corresponds to condition (a) of [theorem %s](#t-umspt), while the compactness of $S$ is a strict version of the drift condition in (b) of the same theorem.[^3] For further results and references, see {cite}`roberts2004general`.

One disadvantage of MCMC methods is that the variates they generate are not in general independent. Sample averages from dependent samples are typically less informative than IID samples. (See [example %s](#eg-llnfail) for the extreme case.) It is desirable, then, that dependence between $x_t$ and $x_{t+k}$ dies out quickly as $k$ increases. Whether it does or not is determined by the choice of $q$. Extensive discussion of such issues can be found in {cite}`brooks2011handbook`.

A function implementing the Metropolis–Hastings algorithm in Julia with a random walk proposal density is given in [listing %s](#list-rw_metropolis). The proposal density takes the form $q(s, s') = \psi(s' - s)$, where $\psi$ is a given density. This is the system {eq}`eq-ampadd`–{eq}`eq-cdncg` with $g$ set to the identity. If $\psi$ is symmetric, then $q(s, s') = q(s', s)$, and $\alpha(s, s')$ simplifies to $\min \{ \pi(s')/ \pi(s) , 1 \}$ if $\pi(s) > 0$ and $1$ otherwise. This is the case shown in [listing %s](#list-rw_metropolis), with $\psi$ set to the $\nN(0, \sigma^2)$ distribution. Look up Julia's ternary operator if you find it hard to follow the code.

```{code-block} julia
:caption: The random walk Metropolis–Hastings algorithm (Julia)
:name: list-rw_metropolis
function rw_metropolis(pi_density, T, init=0, sigma=1)
    xvec = Vector{Float64}(undef, T)
    xvec[1] = init
    for t in 1:(T - 1)
        x = xvec[t]
        y = x + sigma * randn()
        alpha = pi_density(x) > 0 ? pi_density(y) / pi_density(x) : 1 
        xvec[t+1] = rand() < alpha  ? y : x
    end
    return xvec
end
```

{numref}`f-rw_metropolis` shows a histogram of $10^6$ draws generated from this algorithm when $\pi$ is the beta distribution with parameters $(2, 5)$. The black line is the target density. This is just a toy example, since routines for sampling from beta distributions are readily available.

```{figure} figures/rw_metropolis.svg
:name: f-rw_metropolis

Observations from the random walk Metropolis–Hastings algorithm
```

(s-frsp)=
## Further Reading
The modern theory of Markov processes owes much to a brilliant mathematician called Wolfgang Doeblin. His short life and shining contributions to probability theory are discussed in {cite}`bru2002comments`.

An excellent overview of Markov processes and their stability can be found in {cite}`meyn2009markov`. For more on Markov processes from a density perspective see {cite}`lasota1994chaos` or chapter 8 of {cite}`stachurski2009economic`. For a nice survey of sample path properties with a focus on the central limit theorem, see {cite}`jones2004markov`.

{cite}`diaconis2009markov` provides a fascinating overview of MCMC and its impact on different fields. {cite}`richey2010evolution` is a more general survey. Good references on the theoretical side include {cite}`roberts2004general`, {cite}`tierney1994markov`, and {cite}`geweke1996monte`. See {cite}`kroese2014statistical` for a more practical discussion with a focus on applications.

## Exercises

```{exercise}
:label: ex-inspi

Show that the process in [example %s](#eg-sgid) is identically distributed.
```


```{exercise}
:label: ex-covsm

Let $x_{t+1} = a x_t + w_{t+1}$, where $\{w_t\}$ is IID standard normal and $\lL(x_0) = \nN(0, 1/(1-a^2))$. Show that $\cov[x_t, x_{t+j}] = a^j (1 - a^2)^{-1}$ for all $j \geq 0$.
```


```{exercise}
:label: ex-cbv

Show that the right-hand side of {eq}`eq-bvcp` converges to zero in $T$.
```


```{exercise}
:label: ex-pfr

Prove [fact %s](#fa-rdli).[^4]
```


```{exercise}
:label: ex-ch07_further_probability-auto-1

Consider the additive shock process in [fact %s](#fa-gsmp). Assume the conditions stated in the fact and let $\EE \| \boldx_0 \|$ be finite. Show that $\EE \| \boldx_t \|$ is bounded in $t$.
```


```{exercise}
:label: ex-ch07_further_probability-auto-2

Consider again the scalar Markov sequence $x_{t+1} = a x_t + w_{t+1}$. Assume that $\{w_t\}$ is IID, having Student's $t$-distribution with 2 degrees of freedom, and that $|a| < 1$. Show that the conditions of [fact %s](#fa-gsmp) are satisfied.
```


```{exercise}
:label: ex-dmds

Let $\{\fF_t\}$ be a filtration. Show that if $\{m_t\}$ is a martingale with respect to $\{\fF_t\}$, then $d_t = m_t - m_{t-1}$ is a martingale difference sequence with respect to $\{\fF_t\}$.
```


```{exercise}
:label: ex-ekg

Let $\{\fF_t\}$ be a filtration, and let $\{m_t\}$ be a martingale with respect to $\{\fF_t\}$. Let $\alpha_t := m_t + 1$, let $\kappa_t := 2 m_t$, and let $\gamma_t := m_t^2$. Which of these is a martingale with respect to $\{\fF_t\}$?
```


```{exercise}
:label: ex-ch07_further_probability-auto-3

Let $\{\eta_t\}$ be an IID sequence of scalar random variables with $\EE \eta_1 = 0$ and $\var \eta_1 = \sigma^2 > 0$. Let $\{\fF_t\}$ be the filtration defined by $\fF_t := \{\eta_1, \ldots, \eta_t\}$, and let $z_t := t \cdot \eta_t$ for each $t$.

(i) Is $\{z_t\}$ IID? Why or why not?

(ii) Is $\{z_t\}$ a martingale difference sequence with respect to $\{\fF_t\}$? Why or why not?
```


````{exercise}
:label: ex-ch07_further_probability-auto-4

Let $\{\eta_t\}$ be an IID sequence of scalar random variables with $\PP\{\eta_1 = 1\} = \PP\{\eta_1 = -1 \} = 0.5$ and let $\{\fF_t\}$ be the filtration defined by $\fF_t := \{\eta_1, \ldots, \eta_t\}$. Let

```{math}
:enumerated: false

m_t := \sum_{j=1}^t \eta_j
        \quad \text{and} \quad
        \kappa_t := m_t^2 - t
```

Show that $\{\kappa_t\}$ is a martingale with respect to $\{\fF_t\}$.
````


```{exercise}
:label: ex-cuwmm

Consider the scalar sequence $x_{t+1} = a x_t + w_{t+1}$, where $\{w_t\} \iidsim \nN(0,1)$ and $x_0 = 0$. Let $\fF_t := \{w_1, \ldots, w_t\}$. Give conditions on $a$ such that

(i) $\{x_t\}$ is a martingale with respect to $\{\fF_t\}$ and

(ii) $\{x_t\}$ is a martingale difference sequence with respect to $\{\fF_t\}$.
```


````{exercise}
:label: ex-ch07_further_probability-auto-5

Let $\{\fF_t\}$ be a filtration. A sequence of random variables $\{\gamma_t\}$ is called **previsible** with respect to $\{\fF_t\}$ if $\gamma_t$ is $\fF_{t-1}$-measurable for every $t$. Show that if $\{x_t\}$ is a martingale with respect to $\{\fF_t\}$, and $\{\gamma_t\}$ is previsible with respect to $\{\fF_t\}$, then the process $\{y_t\}$ defined by

```{math}
:enumerated: false

y_t := \sum_{i=1}^t \gamma_i (x_i - x_{i-1})
```

is a martingale with respect to $\{\fF_t\}$.
````


````{exercise}
:label: ex-simcop

Let $F(s_1, \ldots, s_N) = C(F_1(s_1), \ldots, F_N(s_N))$, where $C$ is a copula and each $F_n$ is a continuous CDF on $\RR$. Let $\lL(u_1, \ldots, u_N) = C$ and let

```{math}
:enumerated: false

(x_1, \ldots, x_N) := (F^{-1}_1(u_1), \ldots, F^{-1}_N(u_N))
```

Prove that $\lL(x_1, \ldots, x_N) = F$.
````


```{exercise}
:label: ex-sta

Show that if $a \geq 0$ and $b$ and $c$ are any real numbers, then $a \min\{b, c\} = \min \{ ab, ac\}$. Next let $\alpha$ be as defined in {eq}`eq-alpha` and let $p(s, s') := q(s, s') \alpha(s, s')$. Show that $\pi(s) p(s, s') = p(s', s) \pi(s')$ for all $s, s'$.
```


### Solutions to Selected Exercises

````{solution} ex-inspi

Here's a proof for the case $b = 0$ and $c=1$. For the induction step, suppose that $\lL(x_t) = P_\infty$, where $P_\infty$ is as defined in {eq}`eq-sdsg`. With $b = 0$ and $c=1$, this becomes $P_\infty = \nN(0, 1/(1-a^2))$. Since $x_{t+1} = a x_t + w_{t+1}$ and linear combinations of independent normals are normal, we know that $x_{t+1}$ is normal. Thus it remains only to show that $\EE x_{t+1} = 0$ and $\var x_{t+1} = 1 / (1 - a^2)$. The first claim is true because, by linearity of expectations, $\EE x_{t+1} = \EE[ a x_t + w_{t+1}] = a \EE[x_t] + \EE[w] = 0$. The second claim also holds because

```{math}
:enumerated: false

\var x_{t+1} 
        = \var[ a x_t + w_{t+1}] 
        = a^2 \var[x_t] + \var[w_{t+1}] + 2 a \cov[x_t, w_{t+1}]
```

Since $x_t$ and $w_{t+1}$ are independent the final term is zero, and we get

```{math}
:enumerated: false

\var x_{t+1} 
        = a^2 \frac{1}{1 - a^2} + 1 
        = \frac{1}{1 - a^2} 
        
```
````


````{solution} ex-covsm

We aim to show that $\cov[x_t, x_{t+j}] = a^j / (1 - a^2)$. Fix $j \geq 1$ and iterate with $x_{t+1} = a x_t + w_{t+1}$ to get $x_{t+j} = a^j x_t + \sum_{k=1}^j a^{j-k} w_{t+k}$. Using this expression and independence of $x_t$ and $w_{t+k}$ gives

```{math}
:enumerated: false

\EE[x_{t+j} x_t] 
         = \EE 
             \left[ 
                 a^j x_t^2 + \sum_{k=1}^j a^{j-k} w_{t+k} x_t
             \right]  
         =  a^j \EE [x_t^2] 
            + \sum_{k=1}^j a^{j-k} \EE[ w_{t+k} ] \EE[ x_t ]
```

Since $\EE[w_{t+k}] = 0$, this gives $\cov[x_{t+j}, x_t] = \EE[x_{t+j} x_t] = a^j \EE [x_t^2]$. The process $\{x_t\}$ is identically distributed ([fact %s](#fa-idmp)) with common variance $1/(1 - a^2)$. The result follows.
````


````{solution} ex-cbv

As shown in [exercise %s](#ex-covsm), we have $\EE[x_t x_{t+j}] = a^j / (1 - a^2)$. Hence

```{math}
:enumerated: false

\frac{2}{T^2} \sum_{j=0}^{T-1} \sum_{t=1}^{T-j} \frac{a^j}{1 - a^2}
        = \frac{2}{T^2}  \sum_{j=0}^{T-1} (T-j) \frac{a^j}{1 - a^2}
        \leq \frac{2}{T}  \frac{1}{(1-a)(1 - a^2)}
```

The term on the right converges to $0$ as $T \to \infty$.
````


````{solution} ex-pfr

Fix $B \in \bB(S)$ and let $P_t$ be the $t$th marginal. By the law of iterated expectations, we have

```{math}
:enumerated: false

\PP\{\boldx_{t+1} \in B\}
        = \EE[ \1\{\boldx_{t+1} \in B\} ]
        = \EE[ \EE[ \1\{\boldx_{t+1} \in B\} \given \boldx_t ] ]
        = \EE[ Q(\boldx_t, B) ]
```

In other words, $P_{t+1}(B) = \int Q(\bolds, B) P_t(\diff  \bolds)$, as was to be shown.
````


````{solution} ex-dmds

To prove that $\{d_t\}$ is adapted to $\{\fF_t\}$, we need to show that $d_t$ is a function of variables in $\fF_t$. To see this, observe that $m_t$ is a function of variables in $\fF_t$ and $m_{t-1}$ is a function of variables in $\fF_{t-1}$. By the definition of filtrations, every variable in $\fF_{t-1}$ is also in $\fF_t$. Hence $m_{t-1}$ is also a function of variables in $\fF_t$. Since both $m_t$ and $m_{t-1}$ are functions of variables in $\fF_t$, so is $d_t = m_t - m_{t-1}$.

Next we show that $\EE[ d_{t+1} \given \fF_t ] = 0$. Since $\EE[ m_{t+1} \given \fF_t ] = m_t$ and $\EE[ m_t \given \fF_t ] = m_t$ (by [fact %s](#fa-emfv)), we have

```{math}
:enumerated: false

\EE[ d_{t+1} \given \fF_t ] 
        = \EE[ m_{t+1} \given \fF_t]  - \EE[ m_t \given \fF_t ]
        = m_t - m_t = 0
        
```
````


````{solution} ex-ekg

First, $\{\alpha_t\}$ is a martingale with respect to $\{\fF_t\}$. To see that $\{\alpha_t\}$ is adapted to $\{\fF_t\}$, note that $\{m_t\}$ is adapted to $\{\fF_t\}$ by assumption. Therefore $m_t$ is a function of variables in $\fF_t$ and hence $\alpha_t = m_t + 1$ is a function of variables in $\fF_t$. Moreover, using the fact that $\{m_t\}$ is a martingale with respect to $\{\fF_t\}$, we have

```{math}
:enumerated: false

\EE[ \alpha_{t+1} \given \fF_t ] 
        = \EE[ m_{t+1} + 1 \given \fF_t ] 
        = \EE[ m_{t+1} \given \fF_t ] + 1
        =  m_t + 1
        = \alpha_t
```

Second, $\{\kappa_t\}$ is a martingale with respect to $\{\fF_t\}$. The proof is similar and hence omitted.

Third, $\{\gamma_t\}$ is not generally a martingale with respect to $\{\fF_t\}$. For example, let $\{m_t\}$ be the random walk martingale in [example %s](#eg-rw). The process $\{\gamma_t\}$ given by $\gamma_t = m_t^2$ is not a martingale whenever $\sigma^2 := \EE[\eta_1^2]$ is strictly positive. To see this, observe that

```{math}
:enumerated: false

\gamma_{t+1} 
        = m_{t+1}^2 
        = \left( \sum_{j=1}^t \eta_j + \eta_{t+1} \right)^2
        = ( m_t + \eta_{t+1} )^2
        =  m_t^2 + 2 m_t \eta_{t+1} + \eta_{t+1}^2
```

```{math}
:enumerated: false

\fore
        \EE[ \gamma_{t+1} \given \fF_t] 
        = \EE[ m_t^2 \given \fF_t ] 
            + 2  \EE[ m_t \eta_{t+1} \given \fF_t ] + \EE[ \eta_{t+1}^2 \given \fF_t ]
```

Since

```{math}
:enumerated: false

\EE[ m_t \eta_{t+1} \given \fF_t ]
        = m_t \EE[ \eta_{t+1} \given \fF_t ]
        = m_t \EE[ \eta_{t+1} ]
        = 0
```

and $\EE[ \eta_{t+1}^2 \given \fF_t ] = \EE[ \eta_{t+1}^2 ] = \sigma^2$ we now have

```{math}
:enumerated: false

\EE[ \gamma_{t+1} \given \fF_t] 
        = \EE[ m_t^2 \given \fF_t ]  + \sigma^2
        = m_t^2  + \sigma^2
        = \gamma_t + \sigma^2 > \gamma_t
```

Hence $\{\gamma_t\}$ is not a martingale with respect to $\{\fF_t\}$.
````


````{solution} ex-cuwmm

It is clear from {eq}`eq-xaft` that $\{x_t\}$ is adapted to the filtration for all values of $a$. Regarding part (i), if $a = 1$, then $\{x_t\}$ is the random walk in [example %s](#eg-rw). Hence $\{x_t\}$ is a martingale with respect to $\{\fF_t\}$. Regarding part (ii), if $a = 0$, then $x_t = w_t$, and

```{math}
:enumerated: false

\EE[ x_{t+1} \given \fF_t ] 
        = \EE[ w_{t+1} \given \fF_t ]
        = \EE[ w_{t+1} ]
        = 0
```

Hence $\{x_t\}$ is a martingale difference sequence with respect to $\{\fF_t\}$.
````


````{solution} ex-simcop

We aim to show that $\lL(x_1, \ldots, x_N) = F$. We have

```{math}
:enumerated: false

\PP \cap_n \{x_n \leq s_n\}
        = \PP \cap_n \{F^{-1}_n (u_n) \leq s_n\}
        = \PP \cap_n \{u_n \leq F_n(s_n)\}
```

Since $C$ is the CDF of $u_1, \ldots, u_N$, the right-hand side is equal to $C(F_1(s_1), \ldots, F_N(s_N))$. This is $F(s_1, \ldots, s_N)$ by the definition of $F$.
````


```{solution} ex-sta

The first claim is immediate, since $a \geq 0$ implies that $f(x) = ax$ is an order preserving transform. The second claim is just algebra. You've made it this far through the book, so I guess you don't need my help with that.
```


[^1]: The stated results are a consequence of theorem 17.0.1 in {cite}`meyn2009markov`, which in turn requires that the process is $V$-uniformly ergodic. In our case, the latter property follows from aperiodicity and irreducibility, which are immediate consequences of strict positivity of $q$ on $S \times S$, and from the drift condition (b). In the language of {cite}`meyn2009markov`, the function $V$ is "unbounded off petite sets.\" This holds for any coercive function $V$ when $q$ is positive and continuous, since compact sets are then petite.

[^2]: We are assuming that wealth starts at zero and may take arbitrarily large negative values without the gambler getting knee-capped by the mafia.

[^3]: When $S$ is compact, every constant nonnegative function on $S$ is coercive, and if $V$ is constant then the drift condition in {eq}`eq-cwrtc` is trivial.

[^4]: Hint: Try converting probabilities into expectations of indicators and using the law of iterated expectations.
