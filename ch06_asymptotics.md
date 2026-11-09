---
title: "Asymptotics"
label: c-asy
---

In this chapter we cover some fundamental results from probability theory concerning asymptotics of random sequences.

(s-asymp)=
## LLN and CLT
The law of large numbers and central limit theorem are two of the pillars of econometrics and statistics. Below we review both theorems, starting with the necessary concepts of convergence in probability and distribution.

(ss-moc)=
### Convergence of Random Vectors
A sequence of random vectors $\{\boldx_n\}$ is said to **converge in probability** to a random vector $\boldx$ if,

$$
\text{for all } \delta > 0,
    \quad
    \PP \{ \|\boldx_n - \boldx\| > \delta \} \to 0 
    \quad \text{as} \quad
    n \to \infty
$$ (eq-coninp)

In symbols, we write **$\boldx_n \toprob \boldx$**. In the scalar case $\|\boldx_n - \boldx\|$ reduces to $|x_n - x|$.

```{prf:example}
:label: eg-ilcp

If $\lL(\boldx_n) =  \nN(\boldzero,\sigma_n\boldI)$ and $\sigma_n \to 0$, then $\boldx_n \toprob \boldzero$ as $n\to \infty$. {numref}`f-conv_prob` illustrates this visually for the scalar case. The variance is $\sigma_n = 1/n$. With fixed $\delta > 0$, the probability $\PP \{ |x_n | > \delta \}$ is shown for different values of $n$. This probability collapses to zero as $n \to \infty$. If we now fix $\delta$ at a smaller positive value, $\PP \{ |x_n | > \delta \}$ can again be made arbitrarily small by increasing $n$. In other words, {eq}`eq-coninp` holds.
```


```{figure} figures/conv_in_prob.svg
:name: f-conv_prob

$\PP \{ |x_n | > \delta \} \to 0$ when $\lL(x_n) = \nN(0, 1/n)$
```

```{prf:property}
:label: fa-reconpro

The following statements are true:

(i) $\boldx_n \toprob \boldx \iff \| \boldx_n - \boldx \| \toprob 0$.

(ii) $\boldx_n \toprob \boldx \implies g(\boldx_n) \toprob g(\boldx)$ whenever $g$ is continuous at $\boldx$.

(iii) $\boldx_n \toprob \boldx$ and $\boldy_n \toprob \boldy$ $\implies$ $\boldx_n + \boldy_n \toprob \boldx + \boldy$ and $\boldx_n^\T \boldy_n \toprob \boldx^\T \boldy$.

(iv) $\boldx_n \toprob \boldx$ and $\bolda_n \to \bolda$ $\implies$ $\boldx_n + \bolda_n \toprob \boldx + \bolda$ and $\boldx_n^\T \bolda_n \toprob \boldx^\T \bolda$.

(v) $\boldx_n \toprob \boldx$ $\iff$ $\bolda^\T \boldx_n \toprob \bolda^\T \boldx$ for any $\bolda \in \RR^K$.
```


Part (i) is an immediate consequence of the definition. Other proofs are touched on in the exercises. For more details, see, for example, {cite}`cinlar`.

The scalar sequence $\{x_n\}$ is said to converge to $x$ **in mean square** if

$$
\EE (x_n - x)^2   \to 0 
    \quad \text{as} \quad n \to \infty
$$ (eq-cims)

and we write **$x_n \toms x$**. Unlike convergence in probability, for convergence in mean square to be defined we require that our variables have finite second moments.

```{prf:property}
:label: fa-ffci

Let $\{x_n\}$ and $x$ have finite second moments and let $\alpha$ be any constant. The following statements are true:

(i) $x_n \toms x \implies x_n \toprob x$.

(ii) $x_n \toms \alpha$ $\iff$ $\EE x_n \to \alpha$ and $\var[x_n] \to 0$.
```


Part (i) follows from Chebyshev's inequality . In particular, from monotonicity of $\PP$ and {eq}`eq-chein` we get

```{math}
:enumerated: false

\PP\{|x_n - x| > \delta \} 
    \leq \PP\{|x_n - x| \geq \delta \} 
    \leq \frac{\EE (x_n - x)^2}{\delta^2}
```

Part (ii) is implied by the next fact.

```{prf:property}
:label: fa-dmse0

For any $x \in L_2$ and any constant $\alpha$ we have

$$
\EE [ (x - \alpha)^2 ] = \var[x] + (\EE[x] - \alpha)^2
$$ (eq-dmse0)
```


As we'll see later, {eq}`eq-dmse0` plays an important role in statistics. Its verification is an exercise. (The proof is straightforward. If you want more of a challenge, try to obtain it as a consequence of (iv) of [fact %s](#fa-l2opt).)

```{prf:example}
:label: eg-ilcp2

In [figure %s](#f-conv_prob) we showed visually how $\lL(x_n) = \nN(0,1/n)$ implies $x_n \toprob 0$. This also follows from (i)–(ii) of [fact %s](#fa-ffci), since $\EE x_n = 0$ and $\var[x_n] = 1/n \to 0$.
```


(ss-lln)=
### The Law of Large Numbers
The law of large numbers (LLN) is one of the most fundamental theorems of both probability and statistics. As a prelude to the LLN, let's consider the effects of averaging over independent random quantities. We'll take the perspective of investors, who often diversify their portfolios in order to decrease risk. To formalize this idea let

- $x_n$ be the payoff from holding one dollar of asset $n$,

- $\EE x_n = \mu$ and $\var[x_n] = \sigma^2$ for all $n$, and

- $\cov[x_j, x_k] = 0$ when $j \neq k$.

If we hold just asset 1, then the payoff is $x_1$, the expected payoff is $\mu$ and the variance is $\sigma^2$. If we diversify by spreading one dollar evenly over $N$ of these assets, our payoff is

```{math}
:enumerated: false

\bar x_N := \frac{1}{N} \sum_{n=1}^N x_n
```

The expected payoff is unchanged at

$$
\EE  \bar x_N 
        = \EE \left[ \frac{1}{N} \sum_{n=1}^N x_n \right]
        = \frac{1}{N} \sum_{n=1}^N \EE  x_n 
        = \mu
$$ (eq-meanub)

But the variance declines at rate $\frac{1}{N}$ because

```{math}
:enumerated: false

\begin{aligned}
\EE [(\bar x_N - \mu)^2 ]  
    & = \EE \left\{ \left[ \frac{1}{N} \sum_{i=1}^N (x_i - \mu) \right]^2 \right\}
    \\
    & = \frac{1}{N^2} \sum_{i=1}^N \sum_{j=1}^N \EE (x_i - \mu)(x_j - \mu) 
    = \frac{1}{N^2} \sum_{i=1}^N \EE (x_i - \mu)^2 
    = \frac{\sigma^2}{N}
\end{aligned}
```

The important equality here is the third one, which holds because of the zero covariance between assets. In this step $N^2$ terms become $N$ terms. To summarize,

$$
\EE \bar x_N = \mu 
    \quad \text{and} \quad
    \var[ \bar x_N ] = \frac{\sigma^2}{N}
    \quad \text{for all } N
$$ (eq-sosm)

Now imagine someone is kind enough to sell us an unlimited quantity of such assets, so we can diversify without limit. By taking $N \to \infty$ and combining {eq}`eq-sosm` with [fact %s](#fa-ffci) we obtain a proof of the **law of large numbers**:

```{prf:theorem}
:label: t-lln0

Let $\{x_n\}$ be IID copies of $x$. If $x$ is integrable, then

$$
\frac{1}{N} \sum_{n=1}^N x_n \toprob \EE x 
         \quad \text{ as } \quad N \to \infty
$$ (eq-lln0)
```


Actually we haven't fully proved the LLN because we assumed existence of second moments above, while [theorem %s](#t-lln0) assumes only existence of the first moment. However, the proof above captures the spirit of the LLN most cleanly. (See theorem 8.3.5 of {cite:author}`dudley2002real` {cite:year}`dudley2002real` for the general case.)

At first glance, the law of large numbers {eq}`eq-lln0` appears to only be a statement about the sample mean but in fact it's self-improving: We can go from {eq}`eq-lln0` to arbitrary functions of random variables, or even random vectors. In particular, if $\boldx$ is any random vector, $\{\boldx_n\}$ are IID copies and $h \colon \RR^N \to \RR$ is any $\bB$-measurable function such that $h(\boldx)$ is integrable, then

$$
\frac{1}{N} \sum_{n=1}^N h(\boldx_n) \toprob \EE  h(\boldx) 
         \quad \text{ as } \quad N \to \infty
$$ (eq-lln0g)

This can be confirmed by letting $y_n := h(\boldx_n)$ and then applying [theorem %s](#t-lln0). The independence of $\{y_n\}$ is confirmed by [fact %s](#fa-rviifi). That each $y_n$ has the same distribution as $y := h(\boldx)$ follows from

```{math}
:enumerated: false

\PP\{y_n \in B\} 
    = \PP\{h(\boldx_n) \in B\} 
    = \PP\{\boldx_n \in h^{-1}(B)\} 
    = \PP\{\boldx \in h^{-1}(B)\}
```

The law of large numbers applies to probabilities as well as expectations. To see this, fix $B \subset \bB(\RR^N)$, and consider the probability $\PP\{\boldx \in B\}$. Letting $h$ be the indicator function $h(\bolds) = \1_B(\bolds) = \1\{\bolds \in B\}$ and using the principle that expectations of indicators equal probabilities of events, we have

```{math}
:enumerated: false

\EE h(\boldx)  = \EE \1\{\boldx \in B\}  = \PP \{\boldx \in B\}
```

Combining this equality with {eq}`eq-lln0g`, we see that if $\{\boldx_n\}$ is IID with distribution $P$, then

$$
\frac{1}{N} \sum_{n=1}^N \1\{\boldx_n \in B\} \toprob P(B)
$$ (eq-lln0lg)

In words, the fraction of the sample that falls in $B$ converges to the probability that the distribution assigns to $B$.

#### An Illustration

To illustrate the law of large numbers, consider flipping a coin until 10 heads have occurred. The coin is not fair: The probability of heads is 0.4. Let $x$ be the number of tails observed in the process. This random variable is known to have the **negative binomial distribution**, with $\EE x =15$. The LLN predicts that if we simulate a large number of observations of $x$ and take the average, we get a value close to 15. Running the Julia code in [listing %s](#list-nebisi) gives outcomes consistent with this claim.[^1]

```{code-block} julia
:caption: Illustrates the LLN (Julia)
:name: list-nebisi
num_reps = 10^6
outcomes = Vector{Float64}(undef, num_reps)

for i in 1:num_reps
    num_tails = num_heads = 0
    while num_heads < 10
        b = rand()
        num_heads = num_heads + (b < 0.4)   # +1 with prob 0.4
        num_tails = num_tails + (b >= 0.4)  # +1 with prob 0.6
    end
    outcomes[i] = num_tails
end

println(sum(outcomes) / num_reps)
```

{numref}`f-cauchy_samples` shows what can happen when the finite first moment condition in the LLN is not enforced. The draws are from the standard Cauchy distribution, which fails to have a finite first moment. The top panel shows observations and the bottom shows a running sample mean of observations $\bar x_N$, with $N$ on the horizontal axis. A few extreme observations dominate the average and the sample mean does not converge.[^2]

```{figure} figures/cauchy_samples.svg
:name: f-cauchy_samples

Samples from the Cauchy distribution and sample mean
```

### Convergence in Distribution

As a prelude to the central limit theorem we need to review convergence in distribution. This is conceptually quite different to convergence of random variables, whether in probability or mean square. For example, the distributions $\nN(0, 1)$ and $\nN(10^{-10}, 1)$ are very close, but independent draws from these respective distributions need not be.

Let $\{P_n\}$ and $P$ be distributions on $\RR^K$ in the sense of {ref}`ss-muvd`. To represent the notion that $P_n$ converges to $P$ we could require that $P_n(B) \to P(B)$ for all $B \in \bB(\RR^K)$. However, this turns out to be rather strict and also indifferent to closeness in the sense of Euclidean norm. For example, suppose that $P_n$ is a distribution on $\RR$ that puts all mass on $\frac{1}{n}$, while $P$ puts all mass on $0$. Intuition suggests that $P_n \to P$, but if we take $B = (0, \infty)$, we have $P_n(B) = 1$ for all $n$, while $P(B) = 0$.

For this reason the common notion of convergence of distributions, which is called weak convergence, requires only that $P_n(B) \to P(B)$ for all "continuity sets\" in $\RR^K$. There's an equivalent way to phrase this that's easier to work with in applications: We say that $\{P_n\}$ **converges weakly** to $P$ if

```{math}
:enumerated: false

\int h(\bolds) P_n(\diff \bolds) \to \int h(\bolds) P(\diff \bolds)
    \quad 
    \text{for all continuous bounded } h \colon \RR^K \to \RR
```

and write **$P_n \tow P$**. For the definition of these integrals see {ref}`sss-efdv`. Continuity of $h$ enforces respect for the topology on $\RR^K$, while boundedness is imposed only to ensure that the integrals are finite.

````{prf:property}

Let $F_n$ be the CDF of $P_n$ and let $F$ be the CDF of $P$. In the univariate case ($K=1$) we have

```{math}
:enumerated: false

P_n \tow P 
        \quad \iff \quad
        F_n(s) \to F(s) 
        \; \text{ for all } s \text{ at which } F \text{ is continuous}
```
````


```{prf:example}

It can be shown that the $t$-distribution with $k$ degrees of freedom converges weakly to the standard normal distribution as $k \to \infty$; see [figure %s](#f-t_to_norm).
```


```{figure} figures/t_to_norm.svg
:name: f-t_to_norm

$t$-Distribution with $k$ df converges to $\nN(0,1)$ as $k \to \infty$
```

Sometimes densities are easier to work with than CDFs. In this connection, the following sufficient condition can be useful.

```{prf:property}
:label: fa-cdicd

Let $\{P_n\}$ and $P$ be absolutely continuous probability measures on $\RR^K$, with densities $p_n$ and $p$. If $p_n(\bolds) \to p(\bolds)$ for all $\bolds \in \RR^K$, then $P_n \tow P$.
```


Let $\{\boldx_n\}$ and $\boldx$ be random vectors. We say that $\boldx_n \to \boldx$ **in distribution** if their respective distributions converge weakly. This convergence is symbolized by $\boldx_n \tod \boldx$. Thus

```{math}
:enumerated: false

\boldx_n \tod \boldx
    \; \iff \;
    \lL(\boldx_n) \tow \lL(\boldx)
```

In view of [fact %s](#fa-efdvec), this is equivalent to

```{math}
:enumerated: false

\EE[ h(\boldx_n)]  \to \EE[ h(\boldx) ]
    \quad 
    \text{for all continuous bounded } h \colon \RR^K \to \RR
```

```{prf:property}
:label: fa-cmtetc

The following statements are true:

(i) If $g \colon \RR^K \to \RR^J$ is continuous and $\boldx_n \tod \boldx$, then $g(\boldx_n) \tod g(\boldx)$.

(ii) If $\bolda^\T \boldx_n \tod \bolda^\T \boldx$ for any $\bolda \in \RR^K$, then $\boldx_n \tod \boldx$.

(iii) $\boldx_n \toprob \boldx \implies \boldx_n \tod \boldx$.

(iv) If $\bolda$ is a constant vector and $\boldx_n \tod \bolda$, then $\boldx_n \toprob \bolda$.
```


Part (i) of [fact %s](#fa-cmtetc) is called the **continuous mapping theorem**. Part (ii) is called the Cramér–Wold theorem, or the **Cramér–Wold device**. The next result is sometimes known as **Slutsky's theorem**.

```{prf:property}
:label: fa-slut

If $\alpha$ is constant, $x_n \toprob \alpha$ and $y_n \tod y$, then $x_n + y_n \tod \alpha + y$ and $x_n y_n \tod \alpha y$.
```


An immediate but useful consequence is that

```{prf:property}
:label: fa-sluti

$x_n \toprob 0$ and $y_n \tod y$ $\implies$ $x_n y_n \toprob 0$.
```


(ss-clt)=
### The Central Limit Theorem
The **central limit theorem** is among the most striking and important results in all of mathematics. Relative to the LLN, it requires an additional second moment condition. Here's a statement of the theorem in the scalar IID case.

```{prf:theorem}
:label: t-clt0

Let $x$ have finite second moment and let $\{x_n\}$ be IID copies of $x$. If $\mu := \EE x$ and $\sigma^2 := \var x$, then

$$
\sqrt{N} (\bar x_N - \mu) 
        \tod \nN(0, \sigma^2)
         \quad \text{ as } \quad
         N \to \infty
$$ (eq-clt0)
```


The meaning of $\tod$ here is that the distribution of the left-hand side $\tow \nN(0, \sigma^2)$. We can interpret {eq}`eq-clt0` as follows: On one hand, $(\bar x_N - \mu) \toprob 0$ by the LLN; on the other, $\sqrt{N} \to \infty$. If we take the product, these two competing terms just balance. Remarkably, the distribution of the product approaches a zero-mean Gaussian as $N \to \infty$, regardless of the distribution of $x$.

We'll skip the proof of the CLT (which provides little intuition) and turn to simulation. For starters, see [figure %s](#f-clt3d), where $Q_N :=$ the distribution of the left-hand side of {eq}`eq-clt0` for $N=1,\ldots,5$. The initial distribution $Q=Q_1$ is multi-modal, and was constructed as a convex combination of three beta distributions. The CLT predicts that $Q_N \tow \nN(0, \sigma^2)$, where $\sigma^2 = \int s^2 Q_1(\diff s)$. Note how quickly the bell shape emerges, despite the irregularity of the initial distribution.

```{figure} figures/clt3d_2.svg
:name: f-clt3d

CLT in action, starting from a beta mixture
```

Another common statement of the central limit theorem is as follows: If all the conditions of [theorem %s](#t-clt0) are satisfied, then

$$
z_N := \sqrt{N} \left\{ \frac{\bar x_N - \mu}{\sigma}  \right\}
    \tod \nN(0,1)
    \quad \text{ as } \quad
    N \to \infty
$$ (eq-clt1)

{prf:ref}`ex-acmt` asks you to confirm this.

The convergence in {eq}`eq-clt1` is illustrated by [listing %s](#list-illus_clt), the output of which is given in [figure %s](#f-illus_clt). The listing generates 5,000 observations of the random variable $z_N$ defined in {eq}`eq-clt1`, where each $x_n$ is $\chi^2(5)$. (The mean of this distribution is 5, and the variance is $2 \times 5 = 10$.) The observations of $z_N$ are stored in the vector `outcomes`. The figure histograms this vector and then superimposes the density of the $\nN(0, 1)$ distribution.

```{code-block} python
:caption: Illustrates the CLT (Python)
:name: list-illus_clt
import numpy as np
import scipy.stats as st

num_reps = 5000
outcomes = np.empty(num_reps)
N, k = 1000, 5     # k = degrees of freedom
chi = st.chi2(k)

for i in range(num_reps):
    xvec = chi.rvs(N)
    outcomes[i] = np.sqrt(N / (2 * k)) * (xvec.mean() - k)

# code to histogram outcomes and plot st.norm.pdf
# for comparison is omitted here
```

```{figure} figures/illus_clt.svg
:name: f-illus_clt

Observations of $z_N$ in {eq}`eq-clt1` when the underlying distribution is $\chi^2(5)$
```

(s-asympe)=
## Extensions
In this section we cover some extensions to the classical LLN and CLT.

(ss-cid)=
### Convergence of Random Matrices
The notion of convergence in probability discussed in {ref}`ss-moc` extends naturally to matrices. Let $\{\boldX_n\}_{n=1}^{\infty}$ be a sequence of random $N \times K$ matrices. We say that $\boldX_n$ converges to a random $N \times K$ matrix $\boldX$ **in probability** and write $\boldX_n \toprob \boldX$ if

```{math}
:enumerated: false

\| \boldX_n - \boldX \| \toprob 0
    \quad \text{as} \quad n \to \infty
```

where $\| \cdot \|$ is the matrix norm defined in {ref}`ss-mn`.

````{prf:property}
:label: fa-cmtetcv1

Assuming conformability, the following statements are true:

(i) If $\boldX_n \toprob \boldX$ and $\boldX_n$ and $\boldX$ are nonsingular, then $\boldX_n^{-1} \toprob \boldX^{-1}$.

(ii) If $\boldX_n \toprob \boldX$ and $\boldY_n \toprob \boldY$, then

```{math}
:enumerated: false

\boldX_n + \boldY_n \toprob \boldX + \boldY,
                    \quad
                    \boldX_n \boldY_n \toprob \boldX \boldY,
                    \quad \text{and} \quad
                    \boldY_n \boldX_n \toprob \boldY \boldX
```

(iii) If $\boldX_n \toprob \boldX$ and $\boldA_n \to \boldA$, then

```{math}
:enumerated: false

\boldX_n + \boldA_n \toprob \boldX + \boldA,
                    \quad
                    \boldX_n \boldA_n \toprob \boldX \boldA,
                    \quad \text{and} \quad
                    \boldA_n \boldX_n \toprob \boldA \boldX
```

(iv) $\boldX_n \toprob \boldX$ if and only if $\boldX_n \bolda \toprob \boldX \bolda$ for any conformable vector $\bolda$.

(v) $\bolda^\T \boldX_n \bolda \toprob \bolda^\T \boldX \bolda$ whenever $\bolda$ is a conformable constant vector and $\boldX_n \toprob \boldX$.
````


In part (iii) of [fact %s](#fa-cmtetcv1), the matrices $\boldA_n$ and $\boldA$ are nonrandom and convergence is in matrix norm. Note that part (v) can be obtained from two applications of part (iii).

In econometrics we often use the vector version of Slutsky's theorem:

````{prf:property}
:label: fa-cmtetcv2

Let $\boldx_n$ and $\boldx$ be random vectors in $\RR^K$, let $\boldY_n$ be random matrices, and let $\boldC$ be a constant matrix. Assuming conformability, we have

```{math}
:enumerated: false

\boldY_n \toprob \boldC \text{ and } \boldx_n \tod \boldx
        \quad \implies \quad
            \boldY_n \boldx_n \tod \boldC \boldx
            \quad \text{and} \quad
            \boldY_n + \boldx_n \tod \boldC + \boldx
```
````


### Vector-Valued LLNs and CLTs

The scalar LLN and CLT that we discussed above extend to the vector case in a natural way. For example,

```{prf:theorem}
:label: t-vllnclt

Let $\boldx$ be a random vector in $\RR^K$ and let $\{\boldx_n\}$ be IID copies of $\boldx$. If $\boldmu := \EE \boldx$ is finite, then

$$
\bar \boldx_N :=
        \frac{1}{N} \sum_{n=1}^N \boldx_n \toprob \boldmu
         \quad \text{ as } \quad N \to \infty
$$ (eq-vlln)
If, in addition, $\EE \|\boldx\|^2 < \infty$, then

$$
\sqrt{N} \left( \bar \boldx_N - \boldmu \right)
        \tod \nN(\boldzero, \Sigma)
        \quad \text{where } \;
        \Sigma := \var \boldx
$$ (eq-vclt)
```


Here $\frac{1}{N} \sum_{n=1}^N \boldx_n$ should be understood in terms of vector addition and scalar multiplication. Thus $\bar \boldx_N$ is a convex combination of the sample points in $\RR^K$. {numref}`f-vector_mean` illustrates this in two dimensions. The green dot is $\boldmu$. The black dots are realizations of $\boldx_1,\ldots,\boldx_N$. The red dot is the sample mean $\frac{1}{N} \sum_{n=1}^N \boldx_n$. By {eq}`eq-vlln`, the red dot converges to the green dot.

```{figure} figures/vector_mean.svg
:name: f-vector_mean

LLN, vector case
```

The vector LLN in [theorem %s](#t-vllnclt) follows from the scalar LLN. To see this, let $\boldx_n$ be as in [theorem %s](#t-vllnclt), let $\bolda$ be any constant vector in $\RR^K$, let $y_n := \bolda^\T \boldx_n$ and let $y := \bolda^\T \boldx$. The sequence $\{y_n\}$ is IID (see [fact %s](#fa-rviifi)) with the same distribution as $y$. By the scalar LLN ([theorem %s](#t-lln0)), we have

```{math}
:enumerated: false

\frac{1}{N} \sum_{n=1}^N y_n \toprob \EE y 
    = \EE[\bolda^\T \boldx]
    = \bolda^\T\EE[ \boldx]
    = \bolda^\T \boldmu
```

At the same time,

```{math}
:enumerated: false

\frac{1}{N} \sum_{n=1}^N y_n 
    = \frac{1}{N} \sum_{n=1}^N \bolda^\T \boldx_n 
    = \bolda^\T \left[ \frac{1}{N} \sum_{n=1}^N  \boldx_n \right]
    = \bolda^\T \bar \boldx_N
```

We have shown that

```{math}
:enumerated: false

\bolda^\T \bar \boldx_N \toprob \bolda^\T \boldmu \; \text{ for any }\;
   \bolda \in \RR^K
```

The claim $\bar \boldx_N \toprob \boldmu$ now follows from [fact %s](#fa-reconpro).

The vector CLT in [theorem %s](#t-vllnclt) also follows from the scalar case and the proof is quite similar. See [exercise %s](#ex-vclt).

These results can be further extended to random matrices. In fact we will need a matrix LLN in our investigation of linear regression. Here's a suitable result:

```{prf:property}
:label: fa-llnmat

Let $\boldX$ be a random matrix and let $\{\boldX_n\}$ be IID copies of $\boldX$. If $\EE \| \boldX \|< \infty$, then

$$
\frac{1}{N} \sum_{n=1}^N \boldX_n \toprob \EE  \boldX 
         \quad \text{ as } \quad N \to \infty
$$ (eq-llnmat)
```


The proof is straightforward: By [fact %s](#fa-cmtetcv1), to prove {eq}`eq-llnmat`, it suffices to show that $\frac{1}{N} \sum_{n=1}^N \boldX_n \bolda \toprob \EE [ \boldX] \bolda$ for any conformable vector $\bolda$. Since $\boldX_n \bolda$ is a vector with expectation $\EE [ \boldX ] \bolda$, this is immediate from the vector LLN in [theorem %s](#t-vllnclt).

### The Delta Method

The asymptotic normality result in the central limit theorem is preserved under linear transformations, as follows from [fact %s](#fa-cmtetcv2). It turns out that the same thing is true even for functions that are only locally almost linear—that is, for differentiable functions. This greatly extends the power and reach of the central limit theorem.

```{prf:theorem}
:label: t-dm

Let $g \colon \RR^K \to \RR$, let $\boldtheta$ be a point in the domain of $g$, and let $\{\boldt_n\}$ be a sequence of random vectors in $\RR^K$. If

(i) $\sqrt{n} (\boldt_n - \boldtheta) \tod \nN(0, \boldSigma)$ for some positive definite $\boldSigma$ and

(ii) $\nabla g(\boldtheta)$ exists, is continuous, and each element is nonzero,

then

$$
\sqrt{n} \{ g(\boldt_n) - g(\boldtheta) \}
        \tod \nN(0, \nabla g(\boldtheta)^\T \boldSigma \nabla g(\boldtheta))
         \quad \text{ as } \quad
         n \to \infty
$$ (eq-dmmv)
```


Here $\nabla g(\boldtheta)$ is the **gradient vector** of $g$ at $\boldtheta$, which is defined as

```{math}
:enumerated: false

\nabla g(\boldtheta)
    := 
    \begin{pmatrix}
         g'_1(\boldtheta)
         \\
         \vdots
         \\
         g'_K(\boldtheta)
    \end{pmatrix}
    \quad \text{where} \quad
    g'_k (\boldtheta) := \frac{\partial g(\boldtheta)} {\partial \theta_k}
```

In the scalar case, {eq}`eq-dmmv` translates to

$$
\sqrt{n} \{ g(t_n) - g(\theta) \}
    \tod \nN(0, g'(\theta)^2 \sigma^2)
     \quad \text{ as } \quad
     n \to \infty
$$ (eq-dm0)

In statistical settings, the technique of obtaining asymptotic distributions by appealing to [theorem %s](#t-dm) is called the **delta method**. We will see its importance later on. The proof of [theorem %s](#t-dm) is based on a Taylor expansion of $g$ around the point $\theta$. See theorem 3.1 of {cite}`van2000asymptotic`. {prf:ref}`ex-dm` walks you through some of the key ideas.

## Further Reading

For high quality treatments of the law of large numbers and the central limit theorem see, for example, {cite}`dudley2002real`, {cite}`taylor2012introduction`, {cite}`durrett2010probability`, {cite}`pollard2002user`, {cite}`dasguptaasymptotic`, or {cite}`cinlar`. An advanced but thorough treatment of the delta method can be found in {cite}`van2000asymptotic`.

## Exercises

```{exercise}
:label: ex-orvnc

Let $\{x_n\}$ be a sequence of random variables satisfying $x_n = y$ for all $n$, where $y$ is a single random variable. Show that if $\PP\{y = -1\} = \PP\{y = 1\} = 0.5$, then $x_n \toprob 0$ fails. Show that if $\PP\{y = 0\} = 1$, then $x_n \toprob 0$ holds.
```


```{exercise}
:label: ex-cipos

Prove the following part of [fact %s](#fa-reconpro): $\boldx_n \toprob \boldx$ and $\boldy_n \toprob \boldy$ $\implies$ $\boldx_n + \boldy_n \toprob \boldx + \boldy$. (The proof is nontrivial if you are still getting used to these kinds of arguments.)
```


```{exercise}
:label: ex-pcwd

Show that $\boldx_n \toprob \boldx$ implies $\bolda^\T \boldx_n \toprob \bolda^\T \boldx$ for every $\bolda \in \RR^K$.
```


```{exercise}
:label: ex-cdnicp

We saw in [fact %s](#fa-cmtetc) that if $x_n \toprob x$, then $x_n \tod x$. Show that the converse is not in general true.
```


```{exercise}
:label: ex-ch06_asymptotics-auto-1

Show that $\boldX_n \toprob \boldX$ if and only if every element of $\boldX_n$ converges in probability to the corresponding element of $\boldX$.[^3]
```


```{exercise}
:label: ex-acmt

Confirm {eq}`eq-clt1`.[^4]
```


```{exercise}
:label: ex-ch06_asymptotics-auto-2

Let $\lL(u) = U[0,1]$ and $x_n := n \1\{0 \leq u \leq 1/n\}$ for $n \in \NN$.

(i) Calculate the expectation of $x_n^2$.

(ii) Show that $x_n \toprob 0$ as $n \to \infty$.

Conclude that $x_n \toprob x$ does not imply $x_n \toms x$.
```


```{exercise}
:label: ex-ch06_asymptotics-auto-3

Let $x$ be any random variable with $\EE x = \mu$ and $\var x = \sigma^2 < \infty$. Show that $x_n := x/n$ converges to zero in probability as $n \to \infty$.
```


````{exercise}
:label: ex-dm

Suppose that $\{t_n\}$ is a sequence of random variables, $\theta$ is a constant, and

```{math}
:enumerated: false

\sqrt{n} (t_n - \theta) \tod \nN(0, \sigma^2) 
        \quad \text{as} \quad
        n \to \infty
```

Let $g \colon \RR \to \RR$ be differentiable at $\theta$ with $g'(\theta) \neq 0$. Taking a first-order Taylor expansion of $g$ around $\theta$, we can write $g(t_n) = g(\theta) + g'(\theta) (t_n - \theta) + R(t_n - \theta)$, where $R(t_n - \theta)$ is a remainder term. As shown in theorem 3.1 of {cite}`van2000asymptotic`, under these conditions we have $\sqrt{n} R(t_n - \theta) \toprob 0$. Using this fact, prove that $\sqrt{n} \{g(t_n) - g(\theta)\} \tod \nN(0, g'(\theta)^2 \sigma^2)$.
````


```{exercise}
:label: ex-vclt

Confirm the claim $\sqrt{N} \left( \bar \boldx_N - \boldmu \right) \tod \nN(\boldzero, \Sigma)$ in [theorem %s](#t-vllnclt).
```


````{exercise}
:label: ex-cmvclt

Let $\{\boldx_n\}$ be an IID sequence of random vectors in $\RR^K$ with $\EE  \boldx_n  = \boldzero$ and $\var \boldx_n = \boldI$. Let

```{math}
:enumerated: false

\bar \boldx_N := \frac{1}{N} \sum_{n=1}^N \boldx_n 
        \quad \text{and} \quad
        y_N := N \cdot \| \bar \boldx_N \|^2
```

What is the asymptotic distribution of $\{y_N\}$?
````


### Solutions to Selected Exercises

```{solution} ex-orvnc

From the definition, $x_n \toprob 0$ means that, given any $\delta > 0$, we have $\PP \{ |x_n| > \delta \} \to 0$. If $x_n = y$ for all $n$ and $\PP\{y = -1\} = \PP\{y = 1\} = 0.5$, then this fails for $\delta := 1/2$, say, since $\PP \{ |x_n| > 1/2 \} = 1$ for all $n$.

However, if $\PP\{y = 0\} = 1$, then for any $\delta > 0$ we have $\PP \{ |x_n| > \delta \} = \PP \{ |y| > \delta \} = 0$. This sequence is constant at zero and hence $x_n \toprob 0$ holds.
```


````{solution} ex-cipos

Let $\boldx_n \toprob \boldx$ and $\boldy_n \toprob \boldy$. Fix $\delta > 0$. Using the triangle inequality and elementary reasoning we can establish the following ordering of events:

```{math}
:enumerated: false

\begin{aligned}
\{ \|\boldx_n + \boldy_n - (\boldx + \boldy)\| > \delta \}
        & \subset \{ \|\boldx_n - \boldx \| + \| \boldy_n - \boldy\| > \delta \}
        \\
        & \subset \{ \|\boldx_n - \boldx \| > \delta/2\}
        \cup \{ \| \boldy_n - \boldy\| > \delta/2 \}
\end{aligned}
```

By monotonicity and subadditivity of $\PP$, we then have

```{math}
:enumerated: false

\PP \{ \|\boldx_n + \boldy_n - (\boldx + \boldy) \| > \delta \}
        \leq \PP \{ \|\boldx_n - \boldx \| > \delta/2\}
        + \PP\{ \| \boldy_n - \boldy\| > \delta/2 \}
```

By hypothesis, both terms on the right-hand side converge to zero in $n$. Hence $\boldx_n + \boldy_n \toprob \boldx + \boldy$ as claimed.
````


````{solution} ex-pcwd

Let $\boldx_n \toprob \boldx$ in $\RR^K$ and fix any $\bolda \in \RR^K$. If $\bolda = \boldzero$ the claim is trivial, so suppose instead that $\bolda$ is nonzero. Let $\delta > 0$. By the Cauchy–Schwarz inequality,

```{math}
:enumerated: false

\{ |\bolda^\T \boldx_n - \bolda^\T \boldx | > \delta \}
        \subset \{ \| \bolda \| \|\boldx_n - \boldx \| > \delta\}
```

Using monotonicity of $\PP$ now gives

```{math}
:enumerated: false

\PP\{ |\bolda^\T \boldx_n - \bolda^\T \boldx | > \delta \}
        \leq \PP\{ \|\boldx_n - \boldx \| > \delta/\| \bolda \|\}
```

Since $\boldx_n \toprob \boldx$, the right-hand side converges to zero.
````


```{solution} ex-cdnicp

On one hand, if $\{x_n\}$ and $x$ are IID and $\nN(0, 1)$, then $\lL(x_n) = \lL(x)$ for all $n$, and hence $x_n \tod x$. On the other hand, $\lL(x_n - x) = \nN(0,2)$ for all $n$. Hence $\PP\{ |x_n - x| \geq \delta\}$ is a positive constant for any $\delta > 0$, and $x_n \toprob x$ fails.
```


````{solution} ex-vclt

Define

```{math}
:enumerated: false

\boldz_N := \sqrt{N} \left( \bar \boldx_N - \boldmu \right)
        \quad \text{and} \quad
        \lL(\boldz) = \nN(\boldzero, \Sigma)
```

We need to show that $\boldz_N \tod \boldz$. To do this, we apply the Cramér–Wold device ([fact %s](#fa-cmtetc), page ) and the scalar CLT ([theorem %s](#t-clt0), page ). To begin, fix $\bolda \in \RR^K$. Observe that

```{math}
:enumerated: false

\bolda^\T \boldz_n := \sqrt{N} \left( \bar y_n - \EE[y_n] \right)
```

where $y_n := \bolda^\T \boldx_n$. Since $y_n$ is IID and

```{math}
:enumerated: false

\var y_n 
        = \var [ \bolda^\T \boldx_n]
        = \bolda^\T \var [ \boldx_n] \bolda
        = \bolda^\T \Sigma \bolda
```

the scalar CLT yields

```{math}
:enumerated: false

\bolda^\T \boldz_N \tod \nN(0, \bolda^\T \Sigma \bolda)
```

Since $\lL(\bolda^\T \boldz) = \nN(0, \bolda^\T \Sigma \bolda)$, we have shown that $\bolda^\T \boldz_N \tod \bolda^\T \boldz$. Since $\bolda$ was arbitrary, the Cramér–Wold device tells us that $\boldz_N$ converges in distribution to $\boldz$.
````


````{solution} ex-cmvclt

By assumption, $\{\boldx_n\}$ is an IID sequence in $\RR^K$ with $\EE  \boldx_n  = \boldzero$ and $\var \boldx_n = \boldI$. It follows from the vector central limit theorem that

```{math}
:enumerated: false

\sqrt{N} \bar \boldx_N \tod \boldz
        \quad \text{where} \quad \lL(\boldz) = \nN(\boldzero, \boldI)
```

Letting $g(\bolds) := \| \bolds \|^2$ and applying the continuous mapping theorem ([fact %s](#fa-cmtetc)), we obtain

```{math}
:enumerated: false

y_N = \| \sqrt{N} \bar \boldx_N \|^2 \tod \| \boldz \|^2 =
        \sum_{k=1}^K z_k^2
```

From [fact %s](#fa-ssnac), we conclude that $y_N \tod \chi^2(K)$.
````


[^1]: A routine to draw from the negative binomial distribution is available from the *Distributions* package, but we've rolled our own for extra transparency.

[^2]: In fact it can be shown that the sample mean is itself standard Cauchy for all $N$.

[^3]: Hint: To prove the "if\" part use [fact %s](#fa-mnbd).

[^4]: Hint: Combine [theorem %s](#t-clt0) and [fact %s](#fa-cmtetc).
