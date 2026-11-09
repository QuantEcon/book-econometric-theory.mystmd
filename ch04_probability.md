---
title: "Foundations of Probability"
label: c-p
---

Probability theory is one of the foundation stones of statistics and econometrics. We'll start at the basics but move relatively quickly.

(s-pm)=
## Probabilistic Models
Probability is a technically demanding field. The main reason is that we want to assign probabilities to all manner of events. The set of possible events can be very large, and we need ways to manage this complexity. In this and the next few sections, we'll review the basic machinery involved and why this machinery is necessary. At the same time, we'll skip several details that would take us too far away from econometric theory. This will cost little in terms of understanding.

One concept worth mentioning before we start is countable and uncountable sets. In essence, a nonempty set $S$ is called **countable** if it is either finite or can be represented as a sequence, and **uncountable** if is not countable. In probability theory this distinction is important. A quick review can be found in {ref}`ss-card`.

(ss-ssev)=
### Sample Spaces and Events
In setting up a probabilistic model, we start with the notion of a **sample space**, which we think of as being a "list\" of all possible outcomes in a given random experiment. In general, the sample space can be any nonempty set, and is usually denoted $\Omega$. A typical element of $\Omega$ is denoted $\omega$. The general idea is that a realization of uncertainty will lead to the selection of a particular $\omega \in \Omega$.

```{prf:example}
:label: eg-ram00

In a random experiment that involves rolling a die once, the set of possible outcomes is naturally represented by $\Omega := \{1,\ldots, 6\}$.
```


```{prf:example}
:label: eg-dart

Imagine that Burton Malkiel's blindfolded monkey throws darts at a dartboard of radius $1$. To record the outcome we impose ordinary Cartesian coordinates with origin at the center of the board. Let $(h, v)$ be a typical location measured on the horizontal and vertical coordinates respectively. A natural sample space is $\Omega := \setntn{(h, v) \in \RR^2}{\| (h, v) \| \leq 1}$. This set is also called the **unit disk** in $\RR^2$.
```


```{prf:example}
:label: eg-roditwi0

Consider an experiment where we roll a die twice. A suitable sample space is the set of integer pairs $(i, j)$ with $1 \leq i,j \leq 6$. The first element $i$ represents the outcome of the first roll, while the second element $j$ represents the outcome of the second roll. Thus $\Omega := \setntn{(i, j)}{ i, j \in \{1,\ldots,6\}} = \{1,\ldots,6\} \times \{1,\ldots,6\}$.
```


Once we have a sample space in place we can start to think about events. In the language of probability theory, an **event** is just a subset of $\Omega$ (with some caveats that will be addressed below). An event $A$ occurs whenever the individual $\omega \in \Omega$ selected in the random experiment happens to lie in $A$. {numref}`f-eao` illustrates.

```{figure} figures/xfig_event_occurs.svg
:name: f-eao

Outcomes and events
```

#### Probabilities and Events

The next thing we need to do is to assign probabilities to outcomes. One idea is to start by assigning appropriate probabilities to every $\omega$ in $\Omega$. We could then derive the probability of events by looking at the probabilities attached to points contained in those events. It turns out, however, that this is not the best approach, at least if we want to consider uncountable $\Omega$. Rather, the standard approach is to *directly assign probabilities to subsets of $\Omega$ instead*.

To see why the first approach is problematic, consider the dart throwing model in [example %s](#eg-dart), where $\Omega$ is the unit disk in $\RR^2$. Imagine that hitting any point in $\Omega$ is "equally likely.\" We can model this by assuming that the probability the dart lands in some region $A$ of the dartboard (some $A \subset \Omega$) is simply proportional to the area of $A$.

Consider now the probability of hitting some individual point $\omega \in \Omega$. Intuitively, this should be less than the probability of hitting any region $A$ containing $\omega$. By choosing $A$ to be sufficiently small in area, we can show that the probability of hitting $\omega$ is less than $\epsilon$ for any given $\epsilon > 0$. We are forced to conclude that the probability of hitting the individual point $\omega$ is zero.[^1]

This example shows that in many sensible models the probability assigned to individual points in the sample space must be zero. If you reflect on this for even a moment you will agree that there is no way to build up a useful probabilistic model with such information alone. Hence we must assign probabilities directly to the events themselves.

(ss-soev)=
#### Choosing the Set of Events
Should we assign probabilities to every possible event—that is, to every possible subset of $\Omega$? In the dart model above, we could say that the probability the dart lands in any region $A \subset \Omega$ is equal to the area of $A$ divided by the area of $\Omega$ itself, which is $\pi$. In symbols,

$$
\PP(A) = \frac{\lambda(A)}{\pi}
    \quad \text{where } \, 
    \lambda(A) := \text{ the area of set } A \text{ }
$$ (eq-apar)

This seems sensible enough. For example, if $A = \Omega$, then $\PP(A) = \pi/\pi = 1$, as it should be. If $A$ is the top half of the dartboard then $\PP(A) = 0.5$, and so on.

Unfortunately, if we think carefully about this approach, some cracks start to appear. What exactly do we mean by "area of $A$\"? If $A$ is a circle, then its area is $\pi$ times the square of its radius. If $A$ is a rectangle, then its area is the product of its two sides. But how about the area of more general subsets of $\Omega$?

It turns out that assigning area to all subsets of $\Omega$ in a consistent way is problematic. The reason is that $\Omega$ contains an awful lot of subsets, and some of them can be manipulated to exhibit very strange phenomena. (If you are curious, look up the Banach–Tarski paradox for a hint of the problems that can occur when we try to assign notions like area or volume to arbitrary sets.) We can dream up all sorts of strange events that confound our ability to construct a sensible probabilistic model.

Because of this, when $\Omega$ is uncountable, we usually resist taking the set of events to be all subsets of $\Omega$. Instead, we take them to be some family of "well-behaved\" subsets of $\Omega$, denoted below by $\fF$, and only assign probabilities to elements of $\fF$.

There is of course a trade-off here. If we restrict $\fF$ to be a small fraction of the subsets of $\Omega$, then assigning probabilities becomes easy and we can eliminate strange behavior. But if we make $\fF$ too small, then we might not be able to assign probabilities to events that we care about. For example, in the dartboard model, if we restrict $\fF$ to be the set of rectangles, then we can easily measure area and hence assign probabilities via {eq}`eq-apar`. But this excludes interesting events such as the dart landing in the top half of the board.

How can we ensure that $\fF$ is large enough? To think about this in a systematic way, suppose that $A \in \fF$, and that $\PP(A)$ is well-defined and represents the "probability of event $A$." Now, given that we can assign a probability to the event $A$, it would be unfortunate if we couldn't assign a probability to the event "not $A$", which corresponds to $A^c$. So normally we require that if $A \in \fF$, then $A^c \in \fF$. When this is true, we say that $\fF$ is "closed under the taking of complements."

Further, suppose that $A$ and $B$ are both in $\fF$, and we can assign probabilities to these events. In this case, it would be natural to think about the probability of the event "$A$ or $B$", which corresponds to $A \cup B$. So we also require that if $A$ and $B$ are in $\fF$, then $A \cup B$ is likewise in $\fF$. We say that $\fF$ is "closed under the taking of unions." In fact for technical reasons it turns out that we prefer $\fF$ to be closed under the taking of countable unions, which means that if $A_1, A_2, \ldots$ is a sequence of sets in $\fF$, then its union is likewise in $\fF$.

A set $\fF$ of subsets of a given set $\Omega$ that satisfies all these properties is called a $\sigma$-algebra. Formally, $\fF$ is a **$\sigma$-algebra** on $\Omega$ if

(i) $A \in \fF \implies A^c \in \fF$,

(ii) $A_1, A_2, \ldots \in \fF \implies \cup_{n=1}^{\infty} A_n \in \fF$, and

(iii) $\Omega \in \fF$.

Note that (i) and (iii) together imply that $\emptyset \in \fF$. In this context, the empty set $\emptyset$ is called the **impossible event**. In contrast, $\Omega$ is called the **certain event** because it always occurs (regardless of which outcome $\omega$ is selected, $\omega \in \Omega$ is true by definition).

The $\sigma$-algebra of events $\fF$ we choose to pair with a given $\Omega$ can vary from problem to problem. There is however one particular $\sigma$-algebra that you should know at least something about, and that is the standard $\sigma$-algebra we use to pair with ordinary vector space $\RR^N$. This $\sigma$-algebra is called the Borel sets.

To understand the idea, consider the following: The property of $\fF$ being a $\sigma$-algebra seems to imply that $\fF$ is a relatively rich class of sets. After all, we can take complements and unions (and intersections too—see [ex. %s](#ex-dma)) without ever leaving $\fF$. However, there are trivial examples of $\sigma$-algebras, such as $\fF := \{ \emptyset, \Omega\}$. So if we want to pin down a large and useful class of events, the best way to proceed is to assume not only that $\fF$ is a $\sigma$-algebra, but also that $\fF$ contains some fundamental sets that we are likely to want to attach probabilities to, such as intervals and rectangles.

In this connection, the **Borel subsets of $\RR^N$** are defined as the smallest $\sigma$-algebra of subset of $\RR^N$ that contains all rectangles in $\RR^N$. (See {ref}`s-setfun` if you wish to learn more about rectangles.) This definition turns out to be highly successful: Restricting attention to Borel sets does a great job of excluding the kind of strange behavior found in the Banach-Tarski paradox, and, moreover, every subset of $\RR^N$ we need to deal with in day to day analysis is a Borel set. These include planes and hyperplanes, circles, spheres, polygons, finite sets, and sequences of points. In what follows we use $\bB(\RR^N)$ to represent the Borel sets of $\RR^N$.

(ss-probs)=
### Probabilities
Now let's take a $\sigma$-algebra of events as given and consider how to assign probabilities to elements of $\fF$. For given event $B$, the symbol $\PP(B)$ will be interpreted as representing "the probability that event $B$ occurs.\" The way you should think about it is this:

> $\PP(B)$ represents the probability that when uncertainty is resolved and some $\omega \in \Omega$ is selected by "nature," the statement $\omega \in B$ is true.

In order to make sure our model of probability is well behaved, we need to put certain restrictions on $\PP$. For example, it wouldn't do to have a $B$ with $\PP(B) = -93$, since negative probabilities don't make sense. The standard restrictions are imposed in the next definition:

Let $\Omega$ be a nonempty set and let $\fF$ be a $\sigma$-algebra of subsets of $\Omega$. A **probability** $\PP$ on $(\Omega, \fF)$ is a function from $\fF$ to $[0,1]$ that satisfies

(i) $\PP(\Omega) = 1$ and

(ii) $\PP(\cup_{n=1}^{\infty} A_n) = \sum_{n=1}^{\infty} \PP(A_n)$ for any disjoint sequence of sets $A_1, A_2, \ldots \in \fF$.

In this context $\PP$ is also called a **probability measure**. Together, the triple $(\Omega, \fF, \PP)$ is called a **probability space**.

In the definition, axiom (i) is straightforward. We require $\PP(\Omega) = 1$ because, by construction, every possible $\omega$ lies in the set $\Omega$. Axiom (ii) is called **countable additivity**. The concept of disjointness in the statement of the axiom is pairwise: Any distinct pair $A_i, A_j$ share no points in common. By taking some of the sets as empty, countable additivity can be shown to imply finite **additivity**, which requires that

$$
\PP(A_1 \cup \cdots \cup A_k)
    = \PP (A_1) + \cdots + \PP (A_k)
$$ (eq-add)

whenever $A_1, \ldots, A_k$ are disjoint.

The general notion of additivity is a natural requirement for probabilistic models. It says that to find the probability of a given event, we can determine all the different (think "disjoint\") ways that the event could occur, and then sum their probabilities. Probabilists boost this up to countable additivity because it helps certain limiting arguments work smoothly.

{numref}`f-additivity` gives an illustration of additivity. Suppose that the $N$ dots in the figure are the only possible outcomes. Suppose that each of these occurs with probability $1/N$, and hence the probability of an event is the number of dots it contains divided by $N$. Hence $\PP(A)=2/N$, $\PP(B)=4/N$ and $\PP(C)=3/N$. In addition there are $9$ dots in $A\cup B \cup C$, so

```{math}
:enumerated: false

\PP(A \cup B \cup C) = \frac{9}{N}  
    = \PP(A) + \PP(B) + \PP(C)
```

Notice how disjointness is important in making additivity work here: If the sets were not disjoint we would run the risk of double counting on the right hand side of this sum, breaking the equality.

```{figure} figures/xfig_additivity.svg
:name: f-additivity

Each of the $N$ dots occurs with probability $1/N$
```

````{prf:example}
:label: eg-ram0

Let $\Omega := \{1,\ldots, 6\}$ represent the six different faces of a die, as in [example %s](#eg-ram00). Since $\Omega$ is finite, we can and do take $\fF$ to be the set of all subsets of $\Omega$. We define a probability $\PP \colon \fF \to [0,1]$ by

$$
\PP(A) := \frac{|A|}{6}
        \quad \text{where } \, 
        |A| := \text{ number of elements in set } A
$$ (eq-defpd)
For example, $\PP \{2, 4, 6\} = 3/6 = 1/2$. It's easy to see that $0 \leq \PP(A) \leq 1$ for any $A \in \fF$, and that $\PP(\Omega) = 1$. Regarding additivity, suppose that $A$ and $B$ are two disjoint subsets of $\{1,\ldots, 6\}$. Then $|A \cup B| = |A| + |B|$, since, by disjointness, the number of elements in the union is just the number contributed by $A$ plus the number contributed by $B$. Hence

```{math}
:enumerated: false

\PP(A \cup B) 
         = \frac{|A \cup B|}{6} 
         = \frac{|A| + |B|} {6}
         = \frac{|A|}{6} + \frac{|B|}{6}
         = \PP(A) + \PP(B)
```

This proves additivity for pairs of sets. An analogous argument confirms additivity for any finite collection, in the sense that {eq}`eq-add` holds. Finite additivity is in this case equivalent to countable additivity, since the total number of distinct events is finite.
````


````{prf:example}
:label: eg-ram

Continuing [example %s](#eg-ram0), if we roll the die, the probability of getting an even number should be equal to the probability of getting a 2 plus that of getting a 4 plus that of getting a 6. This works because of additivity. Formally,

```{math}
:enumerated: false

\PP\{2, 4, 6\} 
        = \PP[\{2\} \cup \{ 4\} \cup \{6\}] 
        = \PP\{2\} + \PP\{ 4\} + \PP \{6\} 
        = 1/2
```
````


````{prf:example}
:label: eg-mchip

Imagine flipping a coin $N$ times and recording the outcome of each flip. One sample space for this experiment is

```{math}
:enumerated: false

\Omega_0 := \setntn{(b_1, \ldots, b_N)}
            {b_n \in \{\text{heads, tails}\} \text{ for each } n}
```

Letting zero represent tails and one represent heads yields the more practical space

```{math}
:enumerated: false

\Omega 
        := \{0, 1\}^N 
        :=
        \setntn{ (b_1,\ldots, b_N) }
            { b_n \in \{0, 1\} \text{ for each } n}
```

Thus $\Omega$ is the set of all binary sequences of length $N$. Since $\Omega$ is again finite, we take $\fF$ to be the set of all subsets of $\Omega$. As our probability, we define

```{math}
:enumerated: false

\PP(A) := 2^{-N} |A|
```

where $|A|$ is again the number of distinct elements in $A$. This is a kind of uniform probability, where all outcomes are equally likely. Indeed, if $(b_1,\ldots, b_N)$ is any complete sequence of flips and $A \in \fF$ is the singleton set containing just this draw, then $\PP(A) = 2^{-N}$. Hence all such sequences have the same probability under $\PP$.

Regarding the axioms, it is not hard to check that $0 \leq \PP(A) \leq 1$ for all $A \subset \Omega$ and $\PP(\Omega) = 1$ (note that the number of distinct binary sequences of length $N$ is $2^N$). In addition $\PP$ is additive. {prf:ref}`ex-mchip` asks you to confirm this.
````


```{prf:example}
:label: eg-itlm

Consider again the dartboard model from [example %s](#eg-dart), where $\Omega$ is the unit disk in $\RR^2$. For the event space, the set of all subsets of $\Omega$ is problematic (see {ref}`ss-ssev`), so we take instead $\fF$ to be the set of Borel subsets of $\RR^2$ that lie inside $\Omega$. For $\PP$ we follow the "uniform\" probability assignment given in {eq}`eq-apar`. That is, $\PP(B) = \lambda(B) / \pi$ for every $B \in \fF$. The function $\lambda$ that assigns area to Borel sets is known to be countably additive, in that $\lambda(\cup_n A_n) = \sum_{n=1}^\infty \lambda(A_n)$ whenever these sets are disjoint. From this it follows that $\PP$ is likewise countably additive. Evidently $\PP(\Omega) = 1$.
```


The function $\lambda$ from [example %s](#eg-itlm) that maps Borel sets to their "area\" is formally known as the **Lebesgue measure**. {ref}`ss-measure` provides a brief introduction to this concept.

#### Implications of the Axioms

Now let's go back to the general case, where $(\Omega, \fF, \PP)$ is an arbitrary probability space. From the axioms, we can derive a surprising number of properties. Let's list the key ones, starting with the next fact.

```{prf:property}
:label: fa-bpro

Let $(\Omega, \fF, \PP)$ be a probability space, and let $A, B \in \fF$. If $A \subset B$, then

(i) $\PP(B \setminus A) = \PP(B) - \PP(A)$,

(ii) $\PP(A) \leq \PP(B)$,

(iii) $\PP(A^c) = 1 - \PP(A)$, and

(iv) $\PP(\emptyset) = 0$.
```


(Here $A^c = \Omega \setminus A$. See {ref}`s-setfun` for a review of set-theoretic operations.)

To prove [fact %s](#fa-bpro), note that when $A \subset B$, we have $B = (B \setminus A) \cup A$. Since $B \setminus A$ and $A$ are disjoint, additivity of $\PP$ now gives

```{math}
:enumerated: false

\PP(B) = \PP(B \setminus A) + \PP(A)
    \quad \qquad (\text{whenever } A \subset B)
```

This equality implies (i)–(iv) of [fact %s](#fa-bpro). Rearranging gives (i), while nonnegativity of $\PP$ gives (ii). Specializing to $B = \Omega$ gives (iii), and setting $B = A$ gives (iv).

The property $A \subset B \implies \PP(A) \leq \PP(B)$ is called **monotonicity**, and is both fundamental and intuitive: If $A \subset B$, then we know that $B$ occurs whenever $A$ occurs, since $\omega$ lands in $B$ whenever it lands in $A$. Hence the probability of $B$ should be no smaller. Many results in probability turn on this idea.

````{prf:property}
:label: fa-bpro2

If $A$ and $B$ are any (not necessarily disjoint) events, then

```{math}
:enumerated: false

\PP(A \cup B) = \PP(A) + \PP(B) - \PP(A \cap B)
```
````


{prf:ref}`ex-bpro2` asks you to prove this. It follows from [fact %s](#fa-bpro2) that for any $A, B \in \fF$, we have $\PP(A \cup B) \leq \PP(A) + \PP(B)$. This inequality is called **subadditivity**. Thus probabilities are subadditive over arbitrary pairs of events and additive over disjoint pairs.

#### Dependence and Independence

If $A$ and $B$ are events, then the **conditional probability of $A$ given $B$** is

$$
\PP(A \,|\, B) := 
    \frac{ \PP(A \cap B) }{  \PP(B) }
$$ (eq-efcp)

It represents the probability that $A$ will occur, given the information that $B$ has occurred. For the definition to make sense, we need $\PP(B) > 0$. Events $A$ and $B$ are called independent if $\PP(A \cap B) = \PP(A) \PP(B)$. More generally, events $A_1, \ldots, A_N \in \fF$ are called **independent** if

$$
\PP \left(\cap_{n=1}^N A_n \right) = \prod_{n=1}^N \PP(A_n)
$$ (eq-indeps)

Note that if $A$ and $B$ are independent, then the conditional probability of $A$ given $B$ is just the probability of $A$.

````{prf:example}
:label: eg-roditwi

Recall [example %s](#eg-roditwi0), where we roll a die twice. The sample space is $\Omega := \setntn{(i, j)}{ i, j \in \{1,\ldots,6\}}$. Since $\Omega$ is finite we can again take $\fF$ to be the set of all subsets of $\Omega$. For our probability, we define $\PP(E) := |E| / 36$. (Here elements of $E$ are pairs, so $|E|$ is the number of pairs in $E$.) Now consider the events

```{math}
:enumerated: false

A := \setntn{(i, j) \in \Omega}{i \text{ is even}}
        \quad \text{and} \quad
        B := \setntn{(i, j) \in \Omega}{j \text{ is even}}
```

In this case we have

```{math}
:enumerated: false

A \cap B = \setntn{(i, j) \in \Omega}{i \text{ and } j  \text{ are even}}
```

The events $A$ and $B$ are independent under $\PP$. To check this, recall the basic principle for counting ordered tuples: The total number of possible tuples is the product of the number of possibilities for each element. For example, the number of distinct tuples

```{math}
:enumerated: false

(i, j, k) \text{ where } i \in I, j \in J \text{ and } k \in K
```

is $|I| \times |J| \times |K|$. Hence, the number of elements in $A$ is $3 \times 6 = 18$, the number of elements in $B$ is $6 \times 3 = 18$, and the number of elements in $A \cap B$ is $3 \times 3 = 9$. As a result

```{math}
:enumerated: false

\PP(A \cap B) = 9 / 36 = 1 /4 
            = (18 / 36) \times (18 / 36) = \PP(A) \PP(B)
```

This confirms that $A$ and $B$ are independent.
````


To state the next result, we need the concept of a **partition** of $\Omega$, by which we mean a countable set of events $\{B_m\}_{m \geq 1}$ with the property that distinct pairs are disjoint and $\cup_{m \geq 1} B_m = \Omega$.

```{prf:example}

Let $\Omega = (0, \infty)$ and $B_m = (m-1, m]$ for all $m \in \NN$. The events $\{B_m\}_{m \geq 1}$ form a partition of $\Omega$.
```


The following fact is known as the **law of total probability**.

````{prf:property}
:label: fa-ltp

Let $\{B_m\}_{m \geq 1}$ be a partition of $\Omega$ and let $A$ be any event. If $\PP(B_m) > 0$ for all $m$, then

```{math}
:enumerated: false

\PP(A) = \sum_{m \geq 1} \PP(A \, | \, B_m) \cdot \PP(B_m)
```
````


The proof is as follows: It is a basic result from set theory that $A \cap (\cup_m B_m) = \cup_m (A \cap B_m)$. From this equality, the definition of a partition and countable additivity of $\PP$, we get

```{math}
:enumerated: false

\PP(A) = \PP[A \cap (\cup_{m \geq 1} B_m) ]
     = \PP[\cup_{m \geq 1} (A \cap B_m)] 
     = \sum_{m \geq 1} \PP(A \cap B_m)
```

Applying {eq}`eq-efcp` completes the proof.

````{prf:example}

Consider a bank that lends to two types of customers, type A and type B. Types are not observable to the bank. Type A customers have a default probability of 5%. Type B customers have a default probability of 10%. Of potential customers, 90% are of type B and 10% are of type A. Given a randomly selected customer, we wish to calculate on behalf of the bank the default probability associated with this customer.

To this end, let $D$ be the event that the customer defaults, let $A$ be the event that the customer is of type $A$, and let $B$ be the event that the customer is of type $B$. Since the customer must be one of these two types, we have, by the law of total probability,

```{math}
:enumerated: false

\PP(D) = \PP(D \mid A)\PP(A) + \PP(D \mid B)\PP(B)
```

The right-hand side evaluates to $0.05 \times 0.1 + 0.1 \times 0.9 = 0.005 + 0.09 = 0.095$.
````


Another useful result concerning conditional probability is **Bayes' law**, which states that for any events $A$ and $B$ with positive probability,

$$
\PP(A \given B) = \frac{\PP(B \given A) \PP(A)}{\PP(B)}
$$ (eq-bt)

To prove {eq}`eq-bt`, observe that, from the definition of conditional probability,

```{math}
:enumerated: false

\PP(A \given B)=\frac{\PP(A \cap B)}{\PP(B)}
    \quad \text{and} \quad
    \PP(B \given A) = \frac{\PP(A \cap B)}{\PP(A)}
```

Hence $\PP(A \cap B) = \PP(A \given B)\, \PP(B) = \PP(B \given A)\, \PP(A)$. Rearranging yields {eq}`eq-bt`.

````{prf:example}

Banks use automated systems to detect fraudulent or illegal transactions. Consider a test that responds to each transaction with $P$ or $N$, where $P$ means "positive\" (transaction flagged as fraudulent) and $N$ means "negative\" (transaction flagged as normal). Letting $F$ mean fraudulent, suppose that

- $\PP( P \given F) = 0.99$ (the test flags 99% of fraudulent transactions),

- $\PP( P \given F^c ) = 0.01$ (rate of false positives), and

- $\PP( F ) = 0.001$ (prevalence of fraud).

We want to know the probability of fraud given a positive test. To answer this question, we use Bayes' law

```{math}
:enumerated: false

\PP(F \,|\, P) = \frac{\PP(P \given F)\PP(F)}{\PP(P)}
```

and the law of total probability $\PP(P) = \PP(P \given F) \PP(F) + \PP(P \given F^c) \PP(F^c)$ to get

```{math}
:enumerated: false

\PP(F \,|\, P) 
        = \frac{0.99 \times 0.001}{0.99 \times 0.001 + 0.01 \times 0.999}
        = \frac{11}{122}
        \approx \frac{1}{11}
```

Less than one in ten positives are actually fraudulent.
````


(ss-rvs)=
### Random Variables
In some probability courses, a random variable is defined as a "value that changes randomly," or something to that effect. This is a fairly useless definition. A better definition of a random variable is a function $x$ from a sample space $\Omega$ into $\RR$. Think of it this way:

(i) "Nature" picks out an element $\omega$ in $\Omega$ according to some probability.

(ii) The random variable $x$ sends this $\omega$ into numerical value $x(\omega)$.

Thus a random variable is a map $x$ that converts outcomes in sample space—which can be any kind of object—*into numerical outcomes.* This is valuable because numerical outcomes are easy to order, add, subtract and so forth. In other words, random variables "report\" the outcome of an experiment in a format amenable to further analysis.

````{prf:example}
:label: eg-fufh

Consider an experiment where we flip a coin until we get heads. Let $0$ represent tails and $1$ represent heads. As our sample space we take the set of infinite binary sequences

```{math}
:enumerated: false

\Omega := \setntn{(b_1, b_2, \ldots)}{b_n \in \{0, 1\} \text{ for each } n}
```

Here's a random variable giving the number of flips until the first heads:

```{math}
:enumerated: false

x(\omega) = x(b_1, b_2, \ldots) = \min \setntn{n \in \NN}{b_n = 1}
```

Suppose that we also wish to know whether any heads occur in the first 10 flips. If so we can introduce a new random variable

$$
y(\omega) = y(b_1, b_2, \ldots) 
        := \min \left\{ \sum_{n=1}^{10} b_n , 1 \right\}
$$ (eq-egbin)
that returns a 1 when this event occurs and zero otherwise. As per the definition, $x$ and $y$ are functions from $\Omega$ into $\RR$.[^2]
````


The random variable $y$ in {eq}`eq-egbin` is an example of a **binary** (or **Bernoulli**) random variable, which is a random variable taking values in the set $\{0, 1\}$. There is a generic way to create binary random variables that often comes in handy, using **indicator functions**. If $Q$ is a statement, such as "on the planet Uranus, there exists a tribe of three-headed monkeys," then $\1\{Q\}$ is considered as equal to 1 when the statement $Q$ is true, and zero when the statement $Q$ is false. Another common variation on the notation is, for arbitrary $A \in \fF$,

```{math}
:enumerated: false

\1_A(\omega) 
    := \1\{\omega \in A\} 
    := 
    \begin{cases}
        1 & \text{if } \omega \in A
        \\
        0 & \text{otherwise}
    \end{cases}
```

So defined, $\1_A$ is a binary random variable. In fact any binary random variable $x$ can be expressed as $\1_A$ by setting $A := \setntn{\omega \in \Omega}{x(\omega) = 1}$.

```{prf:property}
:label: fa-cbif

If $A_1, \ldots, A_N$ are subsets of $\Omega$, then

(i) $\1_{\cap_{n=1}^N A_n} = \prod_{n=1}^N \1_{A_n}$ and

(ii) $\1_{\cup_{n=1}^N A_n } = \sum_{n=1}^N \1_{A_n}$ whenever the sets are disjoint.
```


Here equality means when evaluated at any $\omega \in \Omega$. See [exercise %s](#ex-cbif).

#### Notational Conventions

Given random variable $x$, you will often see notation such as

```{math}
:enumerated: false

\{ x \text{ has some property} \}
```

This is a shorthand for the *event*

```{math}
:enumerated: false

\setntn{\omega \in \Omega}{x(\omega) \text{ has some property}}
```

````{prf:example}

Consider the claim that, for any random variable $x$,

$$
\PP\{x \leq a\} \leq \PP\{ x \leq b\}
        \quad \text{whenever}
        \quad
        a \leq b
$$ (eq-exno)
The proof goes as follows: Pick any $a, b \in \RR$ with $a \leq b$. The symbol $\{x \leq a\}$ is shorthand for the event $\setntn{\omega \in \Omega}{x(\omega) \leq a}$. Likewise, $\{x \leq b\}$ is shorthand for $\setntn{\omega \in \Omega}{ x(\omega) \leq b}$. Since $a \leq b$, we have

```{math}
:enumerated: false

\setntn{\omega \in \Omega}{x(\omega) \leq a}
        \subset \setntn{\omega \in \Omega}{ x(\omega) \leq b}
```

(Any $\omega$ giving $x(\omega) \leq a$ must also give $x(\omega) \leq b$.) The result in {eq}`eq-exno` now follows from monotonicity of $\PP$ ([fact %s](#fa-bpro)).
````


````{prf:example}
:label: eg-mchip3

Recall [example %s](#eg-mchip), with sample space $\Omega := \{0, 1\}^N$, events $\fF :=$ all subsets of $\Omega$ and $\PP(A) := 2^{-N} |A|$. Consider a random variable $x$ on $\Omega$ that returns the first element of any given sequence. That is, $x(\omega) = x(b_1, \ldots, b_N) = b_1$. The probability that $x = 1$ is $1/2$. Indeed

```{math}
:enumerated: false

\PP\{x = 1\} 
         := \PP \setntn{\omega \in \Omega}{x(\omega) = 1} 
         = \PP \setntn{ (b_1, \ldots, b_N) \in \Omega }{b_1 = 1}
```

Because the number of length $N$ binary sequences with $b_1 = 1$ is $2^{N-1}$, we have $\PP\{x = 1\} = 2^{-N} 2^{N-1} = 1/2$.
````


Another notational convention you should be aware of is that, for random variables, equalities, inequalities and arithmetic operations should be interpreted *pointwise*. For example,

- $x \leq y \; \iff \; x(\omega) \leq y(\omega)$ for all $\omega \in \Omega$,

- $x = y \; \iff \; x(\omega) = y(\omega)$ for all $\omega \in \Omega$, and

- $z = \alpha x + \beta y \; \iff \; z(\omega) = \alpha x(\omega) + \beta y(\omega)$ for all $\omega \in \Omega$.

#### Random Variables Are Measurable Functions

We've omitted a technical detail in our construction of random variables. It arises because we cannot always take the set of events $\fF$ to be the set of all subsets of $\Omega$ (since assigning probabilities to arbitrarily complex sets causes problems—see the discussion in {ref}`ss-ssev`).

To explain, let $(\Omega, \fF, \PP)$ be any probability space. Let $B$ be any subset of $\RR$ and consider evaluating the probability

```{math}
:enumerated: false

\PP\{x \in B\}
    := \PP\setntn{\omega \in \Omega}{x(\omega) \in B}
```

Here $x$ is some function from $\Omega$ to $\RR$. Now observe that if we put no restrictions on $x$, we might end up in a mess. The reason is that we have no way of being sure that the set $\setntn{\omega \in \Omega}{x(\omega) \in B}$ is in fact an element of $\fF$. Since $\PP$ is only defined over elements of $\fF$ (i.e., for *bona fide* events), it is possible that $\PP\{x \in B\}$ is undefined.

To prevent such embarrassment, it's standard to put some restrictions both on $x$ and $B$ that ensure $\{x \in B\} \in \fF$. For $B$, we naturally restrict attention to $\bB(\RR)$, the Borel subsets of $\RR$, since we know these sets to be both relatively well behaved and common enough to satisfy our purposes (see {ref}`ss-ssev`). For $x$, we simply require that $\{x \in B\} \in \fF$ whenever $B$ is a Borel set. Thus, our formal definition is that a **random variable** on $(\Omega, \fF)$ is a function $x \colon \Omega \to \RR$ satisfying

$$
\setntn{\omega \in \Omega}{x(\omega) \in B} \in \fF
    \quad \text{for all } B \in \bB(\RR)
$$ (eq-dorv)

These kinds of functions are also sometimes called $\fF$-measurable functions. In other texts you will see {eq}`eq-dorv` expressed using preimage notation. In the preimage notation, $x^{-1}(B)$ is all $\omega \in \Omega$ such that $x(\omega) \in B$. This is just the left-hand side of {eq}`eq-dorv`, so we can rewrite {eq}`eq-dorv` as

```{math}
:enumerated: false

x^{-1}(B) \in \fF \quad \text{for all } B \in \bB(\RR)
```

Thus $x$ "pulls back\" Borel sets to events.

(sss-mt)=
#### Measurable Transformations
So far so good. If $x$ is a random variable, then $\{x \in B\}$ is a well-defined event for every $B \in \bB(\RR)$. But suppose now that we want to discuss some transformation of $x$. For example, we might be interested in the properties of $y := e^x$. Is $y$ also a random variable? That is, is $\{y \in B\}$ always a well-defined event?

The answer to this question is affirmative provided that our transformation satisfies a regularity condition called Borel measurability. Formally, $f \colon \RR \to \RR$ is called **Borel measurable**, or **$\bB$-measurable**, if

$$
f^{-1}(B) \in \bB(\RR)
    \quad \text{for all } B \in \bB(\RR)
$$ (eq-bm)

Thus $f$ "pulls back\" Borel sets to Borel sets. Although we omit the details, the class of $\bB$-measurable functions is vast (any continuous function, any increasing function, etc., etc.), and all transformations we work with below are $\bB$-measurable.

So why is $\bB$-measurability of $f$ the key to ensuring that $y = f(x)$ is a random variable? Suppose that $f$ is $\bB$-measurable and that $x$ is a random variable. We claim that $\{y \in B\} \in \fF$ for all $B \in \bB(\RR)$. This holds because

$$
\{y \in B\} = \{f(x) \in B\} = \{x \in f^{-1}(B)\}
$$ (eq-byim)

for all $B \in \bB(\RR)$. By {eq}`eq-bm` we know that $f^{-1}(B)$ is a Borel set. Since $x$ is a random variable, it follows that $\{x \in f^{-1}(B)\}$ is an element of $\fF$.

(ss-exp)=
### Expectations
Our next task is to define expectations of random variables and state their basic properties. Elementary probability courses define expectations by summing or integrating functions over densities. The view of a random variable as a function on a sample space plays no part in these definitions. However, when we move to more sophisticated arguments, particularly those involving a large number of random variables, we need a better definition. This section outlines the modern treatment of expectations found in probability and statistics. As we'll see later, the definition in this section is a generalization of the elementary definition. The two agree whenever the latter is well defined.

For finite random variables (i.e., random variables with finite range) the expectation is just the sum of all possible values of the variable, weighted by the corresponding probabilities. More formally, given probability space $(\Omega, \fF, \PP)$ and random variable $x$ taking only finitely many distinct values $s_1, \ldots, s_J$, the **expectation** of $x$ is defined as

$$
\EE x = \sum_{j=1}^J s_j \, \PP\{x = s_j\}
$$ (eq-expec)

````{prf:example}
:label: eg-expecc

Let's apply this definition to the simplest possible case, which is a random variable $x$ satisfying $x(\omega) = \alpha$ for all $\omega \in \Omega$, where $\alpha$ is some constant scalar value. In this case the sum in {eq}`eq-expec` has only one term, and

```{math}
:enumerated: false

\EE x
        = \alpha \PP\{x = \alpha\} 
        = \alpha \PP\setntn{\omega \in \Omega}{x(\omega) = \alpha} 
        = \alpha \PP(\Omega) 
        = \alpha
```
````


```{prf:example}
:label: eg-brv

To evaluate the expectation of a binary random variable $x$, we apply {eq}`eq-expec` to obtain

$$
\EE x  = 1 \times \PP\{x=1\} + 0 \times \PP\{x = 0\} = \PP\{x=1\}
$$ (eq-expbin)
```


````{prf:example}
:label: eg-mchip4

Recall [example %s](#eg-mchip3), which models $N$ flips of a fair coin. The sample space is $\Omega := \{0, 1\}^N$, the events are $\fF :=$ all subsets of $\Omega$, and $\PP(A) := 2^{-N} |A|$ for all $A \in \fF$. Let $x(\omega) = x(b_1, \ldots, b_N) = \sum_{n=1}^N b_n$. This value corresponds to the number of "heads\" observed in one possible outcome $(b_1, \ldots, b_N)$. The probability $\PP$ we've assigned makes all draws equally likely. Hence a good guess is that half the flips will come up heads, or $\EE x  = N/2$. Let's try to obtain this from {eq}`eq-expec`.

Observe first that $0 \leq x \leq N$. Also, by the definition of $\PP$, for any $k$ we have $\PP\{x = k\} = 2^{-N} |A_k|$, where

```{math}
:enumerated: false

A_k :=
         \{x = k\}
         = \left\{
               (b_1, \ldots, b_N) \in \Omega
                \; : \;
                \sum_{n=1}^N b_n = k
           \right\}
```

It is a standard result from combinatorics that $|A_k| = \binom{N}{k}$, where the right-hand side is the so-called **binomial coefficient** for $N,k$. One identity it satisfies is $\sum_{k=0}^N k \binom{N}{k} = N 2^{N-1}$ for all $N$. Combining this identity with {eq}`eq-expec`, the expectation of $x$ is

```{math}
:enumerated: false

\EE x 
        = \sum_{k=0}^N k \, 2^{-N} |A_k|
        = 2^{-N} \sum_{k=0}^N k \, \binom{N}{k} 
        = \frac{N}{2}
```
````


(sss-earv)=
#### Expectations for Arbitrary Random Variables
Definition {eq}`eq-expec` only applies for random variables with finite range. We still need to define expectation for general $x$. This is done by bootstrapping up from the finite case. The procedure is based on approximating arbitrary random variables with finite random variables. One finite approximation $x_n$ to an arbitrary random variable $x$ is shown in [figure %s](#f-finite_rv_approx). The approximation can be improved without limit if we allow $x_n$ to take a larger and larger number of distinct values. This process yields a sequence of finite random variables $x_n$ converging to $x$ (e.g., see {cite}`cinlar`, lemma 2.16). The expectation of $x$ is then defined as

$$
\EE x := \lim_{n \to \infty} \EE x_n
$$ (eq-extar)

The value $\EE x$ is also referred to as the **Lebesgue integral** of $x$ with respect to $\PP$, with the alternative notation $\EE x = \int x(\omega) \PP(\diff \omega)$.

```{figure} figures/tikz_finite_rv_approx.svg
:name: f-finite_rv_approx

Finite approximation to a general random variable
```

Can we be sure that the limit in {eq}`eq-extar` makes sense? Suppose first that $x$ is nonnegative (i.e., $x(\omega) \geq 0$ for all $\omega \in \Omega$). It can be shown that the sequence $\{x_n\}$ approximating $x$ can then be chosen so that $\EE x_n$ is monotone increasing and hence convergent, and that the limiting value does not depend on the particular approximating sequence $\{x_n\}$. See, for example, {cite}`dudley2002real`, proposition 4.1.5. Hence $\EE x$ is well defined, although we cannot rule out divergence of the sequence, which gives $\EE x = \infty$.

We've assumed that $x$ is nonnegative. If we drop this assumption, we can still treat $x$ with similar methods by writing it as the difference between two nonnegative random variables. Typically this is done by writing $x = x^+ - x^{-}$, where $x^+ := \max\{x, 0\}$ and $x^- := -\min\{x, 0\}$. Both $x^+$ and $x^-$ are nonnegative random variables. The expectation can then be defined as

$$
\EE x := \EE x^+ - \EE x^-
$$ (eq-expec3)

This is fine except that {eq}`eq-expec3` might have the form $\infty - \infty$, which isn't well-defined. For this reason we often restrict attention to **integrable** random variables, which are all random variables $x$ such that $\EE |x| < \infty$. Since $x^+ \leq |x|$ and $x^- \leq |x|$, this is enough to ensure that {eq}`eq-expec3` is well-defined.

The definition of $\EE$ given in this section doesn't really tell us how to calculate $\EE x$ in specific instances, but we'll see how to do so in {ref}`ss-expfd`.

#### Properties of Expectation

Let's summarize the main results concerning expectations and list some properties of this operation.

```{prf:property}
:label: fa-exprop

Given any probability space $(\Omega, \fF, \PP)$, there exists a uniquely defined function $\EE$ that maps each integrable random variable $x$ on $(\Omega, \fF, \PP)$ into a value

$$
\EE x = \int x(\omega) \PP(\diff \omega)
$$ (eq-ealm)
in $\RR$, called the **expectation of $x$** under $\PP$. The function has the following properties:

(i) $\EE \alpha = \alpha$ for all $\alpha \in \RR$.

(ii) $\EE \1_A = \PP(A)$ for all $A \in \fF$.

(iii) $x \leq y \implies \EE x \leq \EE y$.

(iv) $\EE [ \alpha x + \beta y ] = \alpha \EE x + \beta \EE y$ for all integrable $x, y$ and constants $\alpha, \beta$.

If we need to remind ourselves of the underlying probability measure $\PP$ we will write $\mathbb{E}_{\PP} x$ instead of $\EE x$.
```


Since we only know how to take expectations of random variables, the expression $\EE \alpha$ in part (i) should be understood as the expectation of a constant random variable everywhere equal to $\alpha$. (We could write this object more precisely as $\alpha \1_\Omega$.) That the expectation equals $\alpha$ was already established in [example %s](#eg-expecc). It also follows from (ii) and (iv), by taking $x = \1_{\Omega}$ and $\beta = 0$.

Part (ii) follows from [example %s](#eg-brv). It is restated here to emphasize this fundamental link between probabilities and expectations. Part (iii) is called **monotonicity of expectations**. The statement $x \leq y$ means that $x(\omega) \leq y(\omega)$ for all $\omega \in \Omega$. {prf:ref}`ex-moncon` asks you to prove this monotonicity in a special case. For a general proof see lemma 4.1.9 of {cite}`dudley2002real`. Part (iv) is called **linearity of expectations**. Here $\alpha x + \beta y$ should be understood as the random variable $(\alpha x + \beta y)(\omega) := \alpha x(\omega) + \beta y(\omega)$. {prf:ref}`ex-adfrv0` and {prf:ref}`ex-adfrv` ask you to check linearity in some special cases. For the general case see theorem 4.1.10 of {cite}`dudley2002real`.

As an exercise, let's use these results to prove that if $x$ is a finite random variable with range $\{s_j\}_{j=1}^J$ and $h$ is any $\bB$-measurable function, then $h(x)$ has expectation

$$
\EE h(x) = \sum_{j=1}^J h(s_j) \PP\{x = s_j\}
$$ (eq-exh0)

First observe that $\sum_{j=1}^J \1\{x = s_j\} = 1$, and hence we can write $h(x)$ as

```{math}
:enumerated: false

h(x) 
    = h(x) \sum_{j=1}^J \1\{x = s_j\}  
    = \sum_{j=1}^J h(s_j) \1\{x = s_j\}
```

Using linearity of expectations gives $\EE h(x) =  \sum_{j=1}^J h(s_j) \EE \1\{x = s_j\}$. Applying part (ii) of [fact %s](#fa-exprop) leads to {eq}`eq-exh0`.

The next result is another useful link between probabilities and expectations, typically called **Chebyshev's inequality**.

```{prf:property}
:label: fa-chein

For any nonnegative random variable $x$ and any $\delta > 0$, we have

$$
\PP\{x \geq \delta \} \leq \frac{\EE x}{\delta}
$$ (eq-chein0)
```


A common variation of Chebyshev's inequality is the bound

$$
\PP\{|x| \geq \delta \} \leq \frac{\EE x^2}{\delta^2}
$$ (eq-chein)

{prf:ref}`ex-chein` asks you to check both {eq}`eq-chein0` and {eq}`eq-chein`.

(ss-blp)=
### Moments and Co-Moments
Let $x$ be a random variable and let $k \in \NN$. If $x^k$ is integrable, then

- $\EE [x^k]$ is called the **$k$th moment** of $x$.

- $\EE [ (x - \EE x)^k]$ is called the **$k$th central moment** of $x$.

If $\EE [ |x|^k ] = \infty$, then the $k$th moment is said not to exist. For some random variables even the first moment does not exist. For others every moment exists.

```{prf:property}
:label: fa-ofinm

If the $k$th moment of $x$ exists, then so does the $j$th for all $j \leq k$.
```


{prf:ref}`ex-ofinm` asks you to confirm this. The next fact is called the **Cauchy–Schwarz inequality for random variables**. Its proof will be given later, once we've covered orthogonal projection ([ex. %s](#ex-csrv)).

```{prf:property}
:label: fa-csrv

If $x$ and $y$ are random variables with finite second moment, then

$$
\left| \, \EE [xy] \, \right| \leq \sqrt{ \EE[x^2] \EE[y^2] }
$$ (eq-csrv)
```


The second central moment of $x$ is called the **variance** of $x$:

```{math}
:enumerated: false

\var x := \EE [ (x - \EE x)^2 ]
```

This gives a measure of the dispersion of $x$. The **standard deviation** of $x$ is

```{math}
:enumerated: false

\sigma_x := \sqrt{\var x}
```

The **covariance** of random variables $x$ and $y$ is defined as

```{math}
:enumerated: false

\cov[x, y] := \EE [ (x - \EE x)(y - \EE y) ]
```

```{prf:property}
:label: eq-efav

If $x$ and $y$ have finite second moments, then

(i) $\var x$ and $\cov[x, y]$ are finite,

(ii) $\var x = \EE [x^2] - [\EE x]^2$, and

(iii) $\cov[x, y] = \EE [ xy] -\EE[x] \EE[y]$.
```


Part (i) follows from (ii)–(iii), the Cauchy–Schwarz inequality and [fact %s](#fa-ofinm). Parts (ii)–(iii) follow from linearity of $\EE$ and some simple manipulations.

````{prf:property}
:label: fa-vacov

If $x_1,\ldots,x_N$ are random variables and $\alpha_0, \alpha_1,\ldots,\alpha_N$ are constant scalars, then

```{math}
:enumerated: false

\var \left[ \alpha_0 + \sum_{n=1}^N \alpha_n x_n \right]
        =  \sum_{n=1}^N \alpha_n^2 \var[ x_n ]   
                + 2 \sum_{n < m} \alpha_n \alpha_m \cov[x_n, x_m]
```
````


Some simple implications of this fact are that

(i) $\var[\alpha + \beta x] = \beta^2 \var[x]$ and

(ii) $\var[\alpha x + \beta y] = \alpha^2 \var[x] + \beta^2 \var[y] + 2 \alpha \beta \cov[x,y]$.

The **correlation** of $x$ and $y$ is defined as

```{math}
:enumerated: false

\corr[x,y] := \frac{\cov[x,y]}{\sigma_x \, \sigma_y}
```

If $\corr[x,y] = 0$, we say that $x$ and $y$ are **uncorrelated**. Positive correlation means that $\corr[x,y]$ is positive and negative correlation means that $\corr[x,y]$ is negative. The first part of the next fact follows immediately from [fact %s](#fa-csrv); the second is just algebra.

````{prf:property}
:label: fa-acorr

Given any two random variables $x$, $y$ and positive constants $\alpha$, $\beta$, we have

```{math}
:enumerated: false

-1 \leq \corr[x,y] \leq 1 
        \quad \text{ and }\quad  
        \corr[\alpha x, \beta y] = \corr[x,y]
```
````


(sss-blp)=
#### Best Linear Predictors
As an application of covariance and a prelude to our discussion of prediction, let's consider the problem of predicting the value of a random variable $y$ given knowledge of the value of a second random variable $x$ (and also knowledge of the underlying probability distributions, which makes this a problem in probability rather than statistics). Thus we seek a function $f$ such that $f(x)$ is close to $y$ on average. To measure the latter, we will use **mean squared error**, which amounts in this case to

```{math}
:enumerated: false

\EE [ (y - f(x))^2 ]
```

As we'll see in {ref}`ss-ce`, the minimizer of the mean squared deviation over all functions of $x$ is obtained by choosing $f(x) = \EE[y \given x]$, where the right-hand side is the conditional expectation of $y$ given $x$. Here we'll consider the simpler problem of finding a good predictor of $y$ within the class of "linear" functions

```{math}
:label: eq-llrv
:enumerated: false

\llL := \{ \text{ all functions of the form } \ell(x) = \alpha + \beta x \}
```

(We defer to common usage and call such functions linear in this context, despite the fact that they are, in general, affine; see [example %s](#eg-affun).) Thus we consider the problem

$$
\min_{\ell \in \, \llL} \EE[ (y - \ell(x))^2 ] 
    = \min_{\alpha, \beta \in \RR} \EE[ (y - \alpha - \beta x)^2 ]
$$ (eq-llsrf)

If $\alpha$ and $\beta$ solve {eq}`eq-llsrf`, then the function

$$
\ell^*(x) := \alpha^* + \beta^* x
$$ (eq-blpre)

is called the **best linear predictor** of $y$ given $x$.

```{prf:example}
:label: eg-mktbeta

In finance, the relationship between returns on a given asset $R_a$ and returns on a market benchmark $R_m$ is called the **beta** of the asset. It measures exposure to systemic risk, as opposed to idiosyncratic risk specific to the asset. The beta of $R_a$ is often defined as the coefficient $\beta^*$ in the best linear prediction {eq}`eq-blpre` when $x$ is returns on the market benchmark and $y=R_a$.
```


To solve {eq}`eq-llsrf` we expand the square on the right-hand side and use linearity of $\EE$ to write the objective function as

```{math}
:enumerated: false

\psi(\alpha, \beta) 
    := \EE [y^2] - 2 \alpha \EE[y] - 2 \beta \EE[xy] + 2\alpha\beta \EE[x]
        + \alpha^2 + \beta^2 \EE[x^2]
```

Computing the derivatives and solving the first-order conditions yields ([ex. %s](#ex-ctmi)) the minimizers

$$
\beta^* := \frac{\cov[x,y]}{\var[x]}
    \quad \text{and} \quad
    \alpha^* := \EE[y] - \beta^* \EE[x]
$$ (eq-slpp)

If you've studied simple linear least squares regression before, you might realize that $\alpha^*$ and $\beta^*$ are the "population" counterparts of coefficient estimates in an empirical setting. In [example %s](#eg-mktbeta), this allows us to take the coefficient $\beta^*$ as the definition of an asset's beta, independent of any particular data set. If we wish to evaluate $\beta^*$, however, we need to use data, since the underlying probabilities, and hence the expectation operator $\EE$, are not observable. It is at this stage that the population regression is replaced by simple linear least squares using data. More on this below.

(s-dist)=
## Distributions
Any random variable $x$ on a given probability space $(\Omega, \fF, \PP)$ assigns to each Borel set $B \subset \RR$ a number $\PP\{x \in B\}$ indicating the probability that $x$ lies in $B$. While this is an important construction, as a means of specifying probabilities for outcomes in $\RR$, it is somewhat indirect. When building statistical models, it's often more convenient to assign probabilities to subsets of $\RR$ directly. This is typically done by specifying a *distribution* over $\RR$. In this section we look at how one might go about specifying such distributions.

(ss-repdist)=
### Defining Distributions on $\RR$
How can we specify probabilities over subsets of $\RR$? In {ref}`ss-soev`–{ref}`ss-probs` we learned how to construct a set of events $\fF$ on an arbitrary sample space $\Omega$ and define a probability measure $\PP$ on $\fF$. If we want to describe probabilities over $\RR$, then the natural thing to do is to repeat the process with $\Omega$ specialized to $\RR$. For the set of events over $\RR$ we take the Borel sets $\bB(\RR)$, as discussed in {ref}`ss-soev`. A probability measure defined over $\bB(\RR)$ is called a **law** or **distribution**.

Distributions will be represented by symbols such as $P$ and $Q$. Referring back to the definition, a distribution $P$ is a map from $\bB(\RR)$ to $[0, 1]$ such that

(i) $P(\RR) = 1$ and

(ii) $P(\cup_{n=1}^{\infty} B_n) = \sum_{n=1}^\infty P(B_n)$ for any disjoint sequence $\{B_n\}$.

Being a probability measure, $P$ satisfies all the properties in [fact %s](#fa-bpro). If there exists a Borel set $S$ with $P(S) = 1$, then we say that $P$ is **supported** on $S$.

Probability measures are both general and comprehensive. Each one pins down probabilistic information for every event in $\RR$ that we might care about in econometrics or statistics. At the same time, this kind of generality can also be a hindrance. For example, it's not a trivial task to assign a specific value $P(B)$ to every single Borel set.

Fortunately, we can characterize distributions more succinctly. To begin, we recall the notion of a **cumulative distributed function**, or CDF, which is any function $F \colon \RR \to [0, 1]$ satisfying

(i) monotonicity: $s \leq s'$ implies $F(s) \leq F(s')$,

(ii) right-continuity: $F(s_n) \downarrow F(s)$ whenever $s_n \downarrow s$, and

(iii) $\lim_{s \to -\infty} F(s) = 0$ and $\lim_{s \to \infty} F(s) = 1$.

CDFs and distributions on $\RR$ can be put in one-to-one correspondence as follows: First, it's true by definition that any two distributions $Q$ and $P$ satisfying $Q(B) = P(B)$ for all $B \in \bB(\RR)$ are the same distribution, but it's also true that if $Q(B) = P(B)$ for all sets $B$ of the form $(-\infty, s]$, then $Q$ and $P$ are again the same distribution. In other words, a distribution $P$ is entirely characterized by the values of the function

$$
F(s) := P( (-\infty, s] )
    \qquad (s \in \RR)
$$ (eq-ntf)

Moreover this function turns out to have properties (i)–(iii) above. The next fact summarizes the relationship between $F$ and $P$.

```{prf:property}
:label: fa-cdf

The following statements are true:

(i) If $P$ is any distribution on $\RR$, then the function $F$ in {eq}`eq-ntf` is a CDF.

(ii) Given any CDF $F$ on $\RR$, there exists exactly one distribution $P$ satisfying {eq}`eq-ntf`.
```


For a full proof, see {cite}`williams1991probability`, lemma 1.6, or {cite}`dudley2002real`, theorem 9.1.1. Here let's restrict ourselves to showing that the function $F$ in {eq}`eq-ntf` satisfies part (i) of the definition of a CDF. To see why (i) holds for this $F$, observe that $s \leq s'$ implies $(-\infty, s] \subset (-\infty, s']$. Part (ii) of [fact %s](#fa-bpro) then gives $P((-\infty, s]) \leq P((-\infty, s'])$. Hence $F(s) \leq F(s')$ as claimed.

CDFs are relatively easy to specify and work with. Here are some examples.

````{prf:example}
:label: eg-normcdf

The **univariate normal distributions** or **Gaussian distributions** refer to the class of distributions identified by CDFs of the form

```{math}
:enumerated: false

F(s) =
        \frac{1}{\sqrt{2 \pi} \sigma}
        \int_{-\infty}^s
           \exp \left\{ - 
               \frac{(t - \mu)^2}{2\sigma^2} \right\} \diff t
       \qquad (s \in \RR)
```

where $\mu \in \RR$ and $\sigma > 0$. We represent the distribution associated with $(\mu, \sigma)$ by $\nN(\mu, \sigma^2)$. The distribution $\nN(0, 1)$ is called the **standard normal distribution**. We use the symbol $\Phi$ for its CDF. Examples of normal CDFs are shown in [figure %s](#f-normal_cdfs).
````


````{prf:example}
:label: eg-pareto

The **Pareto distributions** are the univariate distributions with CDFs of the form

```{math}
:enumerated: false

F(s) = 
        \begin{cases}
            0   
                & \text{ if } s < s_0
            \\
            1 - \left(
                    \frac{s_0}{s}
                \right)^{\alpha} 
                & \text{ if } s_0 \leq s
        \end{cases}
        \qquad (s \in \RR, \; s_0, \, \alpha > 0)
```

Pareto distributions are often used to model phenomena with heavy right-hand tails, such as the distribution of wealth or income.
````


```{figure} figures/normal_cdfs.svg
:name: f-normal_cdfs

Normal CDFs
```

````{prf:example}
:label: eg-betacdf

The class of **beta CDFs** is given by

```{math}
:enumerated: false

F(s) =
        \begin{cases}
            0 & \text{ if } s \leq 0
            \\
            \frac{1}{B(\alpha, \beta)}
                \int_0^s u^{\alpha-1} (1-u)^{\beta-1} \diff u
                & \text{ if } 0 < s < 1
            \\
            1   & \text{ if } 1 \leq s
        \end{cases}
```

where $\alpha, \beta > 0$. See [figure %s](#f-beta_cdfs). In this example $B(\alpha, \beta)$ is the **beta function**

```{math}
:enumerated: false

B(\alpha, \beta) 
            := \frac{\Gamma(\alpha)\Gamma(\beta)}{\Gamma(\alpha + \beta)}
        \qquad \text{where} \qquad
        \Gamma(a) := \int_0^{\infty} u^{a-1} e^{-u} \diff u
```

The function $\Gamma$ is called the **gamma function**.
````


```{figure} figures/beta_cdfs.svg
:name: f-beta_cdfs

Beta CDFs
```

````{prf:example}
:label: eg-cauchy_cdf

The class of **Cauchy CDFs** is given by

```{math}
:enumerated: false

F(s) 
        = \frac{1}{\pi} \arctan\left(  \frac{s - \tau}{\gamma} \right) 
            + \frac{1}{2}
       \qquad (s \in \RR)
```

The parameters $\tau \in \RR$ and $\gamma > 0$ are the location and scale parameters respectively. If $\tau=0$ and $\gamma=1$, then $F$ is called **standard Cauchy**. See [figure %s](#f-cauchy_cdfs).
````


```{figure} figures/cauchy_cdfs.svg
:name: f-cauchy_cdfs

Cauchy CDFs
```

````{prf:example}
:label: eg-udist0

Given $a < b$, the **uniform CDF** on $[a, b]$ is the CDF

```{math}
:enumerated: false

F(s) = 
        \begin{cases}
            0 & \text{ if } s \leq a
            \\
            \frac{s-a}{b-a} & \text{ if } a < s < b
            \\
            1   & \text{ if }  b \leq s
        \end{cases}
```

We represent this distribution symbolically by $U[a, b]$.
````


(ss-dandpmf)=
### Densities and PMFs
There are two important sub-cases for distributions: the discrete case and the absolutely continuous case. As we'll see, in these two sub-cases distributions have nice visual representations and expectations are easy to compute. Almost all "named\" distributions fall into one of these two categories.

#### The Discrete Case

A distribution $P$ is called **discrete** if it is supported on a countable set; that is, if there exists a countable set $\{s_j\}_{j \geq 1}$ with $P(\{s_j\}_{j \geq 1}) = 1$. For such a $P$ let

```{math}
:enumerated: false

p_j 
    := P\{s_j\} 
    := P(\{s_j\})
    = \text{probability mass on the single point } s_j
```

Property (ii) of $P$ and $P(\{s_j\}_{j \geq 1}) = 1$ imply that $\sum_{j \geq 1} p_j = 1$ ([ex. %s](#ex-psto)). Hence $\{p_j\}_{j \geq 1}$ is an example of a **probability mass function**, or **PMF**, which is any nonnegative sequence (finite or infinite) that sums to unity.

Together, $\{s_j\}_{j \geq 1}$ and $\{p_j\}_{j \geq 1}$ characterize the underlying distribution. For example, we can express the CDF corresponding to $P$ via {eq}`eq-ntf` as

$$
F(s) = \sum_{j \geq 1} \1\{s_j \leq s\} p_j
$$ (eq-cdfdrv)

This is intuitive: $F(s)$ represents the probability of an outcome in $(-\infty, s]$ under $P$. Since $P$ puts all its probability mass on $\{s_j\}_{j\geq 1}$, this should equal the sum of the probabilities of hitting $s_j$ over all $s_j \leq s$.

```{figure} figures/tikz_discrete_cdf.svg
:name: f-discrete_cdf

Discrete CDF
```

{numref}`f-discrete_cdf` gives a visualization when $\{s_j\}_{j \geq 1} =\{s_1, s_2\}$. $F$ is a step function, with a jump up of size $p_j$ at point $s_j$. The CDF is right-continuous but not continuous.

````{prf:example}
:label: eg-binpmf

Given $N \in \NN$ and $\pi \in (0, 1)$, the sequence $\{p_0, \ldots, p_N\}$ defined by

```{math}
:enumerated: false

p_j = {N\choose j}\pi^j(1-\pi)^{N-j}
```

is called the **binomial PMF**. The value $p_j$ represents the probability of $j$ successes in $N$ independent trials, each having success probability $\pi$.
````


```{figure} figures/binomial_pmfs.svg
:name: f-pmfs_binom

Binomial PMFs
```

#### The Absolutely Continuous Case

Another case of interest is what elementary books on probability call the "continuous\" case, where distributions can be described by densities. A **density** is a nonnegative function $p$ on $\RR$ that integrates to 1. We say that distribution $P$ is **represented by density $p$** (or **has density $p$**) if $p$ is a density and

$$
P(B) = \int_B p(s) \diff s
    \qquad \text{for all } B \in \bB(\RR)
$$ (eq-drbd)

Here the integral over $B$ can be understood as

```{math}
:enumerated: false

\int_B p(s) \diff s := \int_{-\infty}^\infty \1_B(s) p(s) \diff s
```

Recall from our discussion in {ref}`ss-repdist` that distributions are characterized by their values on simple sets, like half lines of the form $(-\infty, s]$. Thus, if $p$ is a density and

$$
P((-\infty, s]) 
    = \int_{-\infty}^s p(t) \diff t
    \qquad \text{for all } s \in \RR
$$ (eq-drbd2)

then {eq}`eq-drbd` also holds.

Not every distribution has a density representation. An exact necessary and sufficient condition for existence is absolute continuity. A distribution $P$ on the Borel subsets of $\RR$ is called **absolutely continuous** if $P(B) = 0$ whenever $B$ has Lebesgue measure zero. A definition of the latter property is given in {ref}`ss-measure`, although for us it is enough to know that any countable subset of $\RR$ has Lebesgue measure zero. This leads to the following fact:

```{prf:property}
:label: fa-denc

If $P$ is absolutely continuous, then $P(C) = 0$ whenever $C$ is countable.
```


Compare with the discrete distributions, which put all mass on a countable set. One important implication of [fact %s](#fa-denc) is that individual points receive no probability mass under absolute continuity, and the corresponding CDF contains no jumps. In fact, when {eq}`eq-drbd2` holds, the fundamental theorem of calculus tells us that $F(s) = P((-\infty, s])$ is differentiable at all continuity points of $p$, and, in particular,

```{math}
:enumerated: false

F'(s) = p(s)
    \qquad \text{for all } s \in \RR \text{ such that } p \text{  is continuous at } s
```

````{prf:example}
:label: eg-normden

The normal CDFs from [example %s](#eg-normcdf) are differentiable for all $\mu$, $\sigma$, with density

```{math}
:enumerated: false

p(s) = F'(s) = 
        \frac{1}{\sqrt{2 \pi} \sigma}
           \exp \left\{ - 
               \frac{(s - \mu)^2}{2\sigma^2} \right\}
```

See [figure %s](#f-normal_densities). We reserve the symbol $\phi$ for the standard normal density.
````


```{figure} figures/normal_densities.svg
:name: f-normal_densities

Normal densities
```

````{prf:example}
:label: eg-cauchyden

The Cauchy CDF from [example %s](#eg-cauchy_cdf) has density

```{math}
:enumerated: false

p(s) = 
        \frac{1}{\pi \gamma}
            \left[
                1 + \left( \frac{s - \tau}{\gamma} \right)^2
            \right]^{-1}
            \qquad (s \in \RR, \; \gamma > 0, \, \tau \in \RR)
```

The Cauchy densities are more peaked around their modes and have greater mass in their tails than normal densities. See [figure %s](#f-cauchy_densities).
````


````{prf:example}
:label: eg-betaden

The beta CDFs from [example %s](#eg-betacdf) have densities given by

```{math}
:enumerated: false

p(s) = \frac{s^{\alpha-1} (1-s)^{\beta-1}}{B(\alpha, \beta)}
            \qquad (\alpha, \beta > 0)
```

when $0 < s < 1$ and zero elsewhere.
````


```{figure} figures/cauchy_densities.svg
:name: f-cauchy_densities

Cauchy densities
```

````{prf:example}
:label: eg-udist

The $U[a, b]$ distribution is represented by the density

```{math}
:enumerated: false

p(s) = \frac{1}{b - a} \1\{a \leq s \leq b\}
        \qquad (s \in \RR, \; a, \, b \in \RR, \; a < b)
```
````


````{prf:example}
:label: eg-gammad

The **gamma distribution** with shape parameter $\alpha$ and scale parameter $\beta$ is the distribution with density

```{math}
:enumerated: false

p(s) = 
        \frac{s^{\alpha-1} e^{-s/\beta}}{\beta^{\alpha} \Gamma(\alpha)} 
            \qquad (\alpha, \beta > 0)
```

when $0 < s < 1$ and zero elsewhere.
````


````{prf:example}
:label: eg-chisq

The **chi-squared distribution with $k$ degrees of freedom** is the distribution with density

```{math}
:enumerated: false

p(s) := \frac{1}{2^{k/2} \Gamma(k/2)} s^{k/2 - 1} e^{-s/2} 
        \qquad (s > 0,\; k \in \NN)
```

This distribution is represented by the symbol $\chi^2(k)$. {numref}`f-chisq_densities` illustrates.
````


````{prf:example}
:label: eg-studt

**Student's $t$-distribution with $k$ degrees of freedom**, or, more simply, the $t$-distribution with $k$ degrees of freedom, is the distribution on $\RR$ with density

```{math}
:enumerated: false

p(s) := \frac{\Gamma(\frac{k+1}{2})}{(k\pi)^{1/2} \Gamma(\frac{k}{2})}
        \left( 1 + \frac{s^2}{k} \right)^{-(k+1)/2}
        \qquad (s \in \RR,\; k > 0)
```
````


```{figure} figures/chisq_densities.svg
:name: f-chisq_densities

Chi-squared densities
```

````{prf:example}

The **$F$-distribution** with parameters $k_1, k_2$ is the distribution with the unlikely looking density

```{math}
:enumerated: false

p(s) 
        := \frac{\sqrt{(k_1 s)^{k_1}k_2^{k_2}/[k_1 s + k_2]^{k_1 + k_2}}}
            {s B(k_1/2, k_2/2)}
        \qquad (s \geq 0, \; k_1, \, k_2 > 0)
```

The $F$-distribution arises in a number of hypothesis tests, as discussed below.
````


(ss-igwd)=
### Integrating with Distributions
We all understand what it means to take the ordinary integral $\int_a^b h(s)\diff s$ of a well-behaved function $h$ on some interval $[a, b]$. Sometimes we want to weight this integral, assigning more mass to different regions of $[a, b]$. This leads to integrals such as $\int_a^b h(s) p(s) \diff s$. For example, $h$ could be a welfare function and $p$ might record the density of agents of a certain type. In a probabilistic context, $p$ might be a density indicating probabilities of outcomes, $h$ might be a payoff function and $\int_a^b h(s) p(s) \diff s$ would be expected payoff.

We now wish to give meaning to expressions of the form $\int h(s) P(\diff s)$ that can be viewed as generalizations of these integrals. The nature of the generalization is that the probability $P$ might not have a density, and yet we still want to compute a weighted integral using $P$. As motivation we'll continue with the previous example, where $h(s)$ is the payoff from outcome $s \in \RR$, and the likelihood of outcomes is determined by probability measure $P$. To provide some basic regularity, we assume that $h$ is a $\bB$-measurable function on $\RR$.

Actually our task is straightforward given the foundations already laid. In particular, we already have a general notion of integration suitable for understanding the expression $\int h(s) P(\diff s)$. In {ref}`ss-exp` we saw how to take an arbitrary probability space $(\Omega, \fF, \PP)$ and construct an expectation operator $\EE$ that assigns to any integrable random variable a notion of expectation. Expectations can be thought of as integrals, as suggested by the alternative notation $\int x(\omega) \PP(\diff \omega)$ for $\EE x$ used in {eq}`eq-ealm`.

So now take our distribution $P$ on $\RR$ and consider $(\RR, \bB(\RR), P)$ as a probability space. This new probability space has its own expectation operator $\EEP$. (See [fact %s](#fa-exprop) if you wish to refresh your memory on this point.) Moreover we can use it to take the expectation of the payoff function $h$, which, being a $\bB$-measurable function on $\RR$, is by definition a random variable on $(\RR, \bB(\RR), P)$. This gives us a value $\EEP h$, which is exactly the integral we're after:

```{math}
:enumerated: false

\EEP h :=: \int h(s) P(\diff s) 
    := \text{the expectation of } h \text{ under } P
```

While the preceding discussion gives a neat definition of $\int h(s) P(\diff s)$, it would be disappointing if this definition didn't lead us back to familiar expressions like $\int h(s) p(s) \diff s$ when $P$ has a density. Fortunately, it does:

```{prf:property}
:label: fa-hasp

Let $h \colon \RR \to \RR$ be $\bB$-measurable and let $P$ be a distribution on $\RR$. If $P$ is discrete, with PMF $\{p_j\}_{j \geq 1}$ and support $\{s_j\}_{j \geq 1}$, then

$$
\int h(s) P(\diff s) = \sum_{j\geq 1} h(s_j) p_j
$$ (eq-hasp)
If $P$ is absolutely continuous with density $p$, then

$$
\int h(s)P(\diff s) = \int_{-\infty}^{\infty} h(s) p(s) \diff s
$$ (eq-hasd)
```


For a proof of [fact %s](#fa-hasp), see §6.12 of {cite}`williams1991probability`.

(ss-drvs)=
### Distributions of Random Variables
In {ref}`ss-rvs` we defined random variables. Over the last few sections we've defined and discussed distributions over the Borel subsets of $\RR$. Our next job is to connect them. We'll see that the connection runs in both directions: Every random variable defines a distribution on $\RR$. Conversely, given any distribution, we can reverse engineer a random variable having that distribution.

To begin, let $x$ be a random variable on some probability space $(\Omega, \fF, \PP)$. By construction, the probability $\PP\{x \in B\}$ is well-defined for every $B \in \bB(\RR)$ (see {eq}`eq-dorv`). The values obtained from this expression when $x$ is held fixed and $B$ is varied is called the **distribution of $x$**. In particular, the set function $P$ defined by

$$
P(B) = \PP\{x \in B\} 
    \qquad (B \in \bB(\RR))
$$ (eq-pidox)

is both the distribution of $x$ and also a distribution in the sense of {ref}`ss-rvs` (i.e., a probability measure over the Borel sets). It records the probabilities of all outcomes of interest for the variable $x$ when regarded as a single entity.[^3] In what follows we write $\lL(x)$ for the distribution (or "law\") of $x$. The relationship in {eq}`eq-pidox` is summarized by writing $\lL(x) = P$.

By {eq}`eq-ntf`, the CDF corresponding to the distribution $P$ of $x$ satisfies

$$
F(s) = \PP \{x \leq s\}
    \qquad (s \in \RR)
$$ (eq-deffx)

Sometimes $F$ is called the **distribution function** or just distribution of $x$. Since $F$ characterizes $P$, we also write $\lL(x) = F$ to indicate that $F$ represents the distribution of $x$.

```{prf:property}
:label: fa-axb

If $\lL(x)= F$, then $\PP\{a < x \leq b\} = F(b) - F(a)$ for any $a \leq b$.
```


Indeed, if $a \leq b$, then $\{a < x \leq b\} = \{x \leq b \} \setminus \{x \leq a\}$ and $\{x \leq a \} \subset \{x \leq b\}$. Applying [fact %s](#fa-bpro) gives

```{math}
:enumerated: false

\PP\{a < x \leq b\} 
    = \PP\{x \leq b\} - \PP\{x \leq a\}
    = F(b) - F(a)
```

From {eq}`eq-deffx` we get for each random variable a uniquely defined CDF. It's also true that for every CDF $F$, there exists a probability space $(\Omega, \fF, \PP)$ and a random variable $x \colon \Omega \to \RR$ such that $\lL(x) = F$; {ref}`ss-it` outlines one construction.

In what follows, if $\lL(x) = P$ and $P$ has density $p$, we'll say that **$x$ has density $p$**. If the distribution of $x$ is discrete, we'll call $x$ a **discrete random variable**. In the former case, where the distribution of $x$ is absolutely continuous, no probability mass is placed on finite sets, leading to the next fact.

````{prf:property}

If $x$ has a density, then $\PP\{x=s\}=0$ for all $s \in \RR$, and for any $a < b$,

```{math}
:enumerated: false

\PP\{a < x < b\}
        = \PP\{a < x \leq b\}
        = \PP\{a \leq x < b\}
        = \PP\{a \leq x \leq b\}
```
````


#### Distributions of Transformations

One often needs to obtain the distribution of the transform of a random variable. This is easy in the monotone case.

```{prf:property}
:label: fa-cdftr

If $\lL(x) = F$ and $y := \psi(x)$, where $\psi \colon \RR \to \RR$ is strictly increasing, then $\lL(y)= G$ where $G(s) := F(\psi^{-1}(s))$.
```


To see why, observe that under these hypotheses $\psi^{-1}$ exists and is (weakly) increasing. Hence

```{math}
:enumerated: false

\PP\{y \leq s\} 
    = \PP\{\psi(x) \leq s\} 
    = \PP\{x \leq \psi^{-1}(s) \} 
    = F( \psi^{-1}(s) )
```

Note how monotonicity is used in the second equality. (A related result can be derived for the strictly decreasing case, where the direction of the inequality flips.)

```{prf:example}

If $\lL(x) = F$ and $y := \exp(x)$, then the CDF of $y$ is $G(s) := F(\ln(s))$.
```


An important related question is to determine when $y = \psi(x)$ has a density, and what the form of this density is.

```{prf:property}
:label: fa-trden

If $x$ has density $p$ on $\RR$ and $y := \psi(x)$ where $\psi$ is a diffeomorphism on $\RR$, then the distribution of $y$ is absolutely continuous, with density

$$
q(s) = p(\psi^{-1}(s)) \left| \frac{\diff \psi^{-1}(s)}{\diff s} \right|
                 \qquad (s \in \RR)
$$ (eq-trden)
```


{prf:ref}`ex-trden` asks you to confirm {eq}`eq-trden`. The term **diffeomorphism** just means that $\psi$ is a bijection on $\RR$ and both $\psi$ and its inverse are differentiable.

````{prf:example}
:label: eg-ndat

If $x$ has density $p$ on $\RR$ and $\mu$ and $\sigma$ are constants with $\sigma > 0$, then the density of $y := \mu + \sigma x$ is

```{math}
:enumerated: false

q(s) = p \left( \frac{s - \mu}{\sigma} \right) 
                 \frac{1}{\sigma} 
                 \qquad (s \in \RR)
```
````


When $x$ is standard normal, this implies that $y = \mu + \sigma x$ is $\nN(\mu, \sigma^2)$, as can be seen by taking $p$ to be the standard normal density $\phi$ and comparing with the expression for the $\nN(\mu, \sigma^2)$ density.

[Fact %s](#fa-trden) can be generalized. See, for example, theorem 15.5 of {cite}`schilling2005measures`.

(ss-expfd)=
### Expectations from Distributions
Let $x$ be a random variable on probability space $(\Omega, \fF, \PP)$. The distribution of $x$ encodes all the information necessary to calculate expectations of $x$ or of any $\bB$-measurable transformation $h(x)$. This is clearest when $x$ is finite. For example, suppose that $\lL(x) = P$, that $h \colon \RR \to \RR$ is any $\bB$-measurable function, and that $P$ puts all mass on finite set $\{s_j\}_{j=1}^J$. Let $\{p_j\}$ be the PMF of $P$, so that $p_j = P\{s_j\}$. Applying the definition of expectations (see in particular {eq}`eq-exh0`) and using $\PP\{x = s_j\} = P\{s_j\}$ gives

$$
\EE h(x) 
    = \sum_{j=1}^J h(s_j) \PP\{x = s_j\}
    = \sum_{j=1}^J h(s_j) P\{s_j\}
    = \sum_{j=1}^J h(s_j) p_j
$$ (eq-ehx)

Comparing this with {eq}`eq-hasp`, we can also write $\EE h(x) = \int h(s) P(\diff s)$. Put differently, the expectation of $h(x)$ on $(\Omega, \fF, \PP)$ is equal to the expectation of $h$ on $(\RR, \bB(\RR), P)$. This is also true in the infinite case (see, e.g., theorem 5.2 of {cite}`cinlar`), as stated in the next fact.

````{prf:property}
:label: fa-efd

Let $x$ be a random variable on some probability space $(\Omega, \fF, \PP)$, let $\lL(x) = P$ and let $h$ be a $\bB$-measurable function such that $h(x)$ is integrable. The expectation $\EE h(x)$ is entirely determined by $h$ and $P$. In particular,

```{math}
:enumerated: false

\EE h(x) = \int h(s) P(\diff s)
```

where $\int h(s) P(\diff s)$ is the expectation of $h$ on $(\RR, \bB(\RR), P)$.
````


````{prf:example}

Let $x$ be a random variable whose distribution $P$ is the uniform distribution on $[a, b]$. Applying the definition of the uniform density in [example %s](#eg-udist), we obtain

```{math}
:enumerated: false

\EE x 
        = \int s P(\diff s) 
        = \int s p(s) \diff s 
        = \int_{-\infty}^{\infty} 
         \frac{s}{b - a} \1\{a \leq s \leq b\} \, \diff s
```

Solving the integral gives $\EE x = \mu := (a + b)/2$. The variance is

```{math}
:enumerated: false

\var [x]
        = \int (s - \mu)^2 P(\diff s)
        = \int_a^b
            \left(s -\frac{a + b}{2} \right)^2
            \frac{1}{b - a} \diff s
        = \frac{1}{12} (b-a)^2
```
````


````{prf:example}

Suppose that $\lL(x) = \nN(\mu, \sigma)$. In the case where $\sigma > 0$, the mean can be computed via

```{math}
:enumerated: false

\EE x 
        = \int_{-\infty}^{\infty} s 
        \frac{1}{\sqrt{2 \pi} \sigma}
           \exp \left\{ - 
               \frac{(s - \mu)^2}{2\sigma^2} \right\}  \diff s
        = \mu
```

The variance is given by

```{math}
:enumerated: false

\var [x]
        = \int_{-\infty}^{\infty}
            (s - \mu)^2
            \frac{1}{\sqrt{2 \pi} \sigma}
               \exp \left\{ - 
               \frac{(s - \mu)^2}{2\sigma^2} \right\}  \diff s
        = \sigma^2
```

For a discussion of how these integrals are solved, see §9.4 of {cite}`dudley2002real`.
````


#### Moments from Distributions

In {ref}`ss-blp` we defined moments of random variables. For example, $\EE x^k$ is the $k$th moment of $x$, assuming it exists. It is clear from [fact %s](#fa-efd) that any two random variables with the same distribution share the same moments. Hence moments are best thought of as a property of the distribution, not the random variable. Thus we define

- the **mean** of $P$ as $\mu = \int s P(\diff s)$,

- the **$k$th moment** of $P$ as $\int s^k P(\diff s)$,

- the **variance** of $P$ as $\int (s - \mu)^2 P(\diff s)$,

and so on.

(ss-quantile)=
### Quantile Functions
In this section it will be convenient to describe distributions in terms of CDFs. So let $F$ be a CDF on $\RR$, which, for the moment, we assume to be strictly increasing. Given $\tau \in (0, 1)$, the **$\tau$th quantile** of $F$ is the $\xi \in \RR$ that solves $F(\xi) = \tau$. Under our assumptions on $F$, such a $\xi$ exists and is uniquely defined. The $0.5$th quantile is called the **median** of $F$. It provides a measure of central tendency for the distribution.

We can of course express the $\tau$th quantile of $F$ using the inverse function:

$$
F^{-1}(\tau) := \text{ the unique } \xi \text{ such that } F(\xi) = \tau
    \qquad (0 < \tau < 1)
$$ (eq-quantwi)

Not surprisingly, the inverse of $F$ is also called the **quantile function**.

```{prf:example}

The quantile function associated with the standard Cauchy distribution is $F^{-1}(\tau) = \tan [\pi(\tau - 1/2)]$. {numref}`f-cauchy_quant` illustrates. The horizontal axis is $\tau \in (0, 1)$.
```


```{figure} figures/cauchy_quant.svg
:name: f-cauchy_quant

Cauchy quantile function
```

If $F$ is not strictly increasing, then there exist at least two distinct points $s, s'$ such that $F(s) = F(s')$, and hence $F^{-1}$ in {eq}`eq-quantwi` is not well-defined. This problem is negotiated by setting

$$
F^{-1}(\tau) := \inf \setntn{s \in \RR}{F(s) \geq \tau}
    \qquad (0 < \tau < 1)
$$ (eq-quantwi2)

(See {ref}`sss-sai` for a definition of inf.) In the case where $F$ is strictly increasing, {eq}`eq-quantwi2` reduces to {eq}`eq-quantwi`.

The quantile function is used in hypothesis testing to determine critical values. Let's outline a common scenario.

````{prf:property}
:label: fa-avsc

Let $x$ be a random variable with density $p$. If $p$ is symmetric, then the CDF $G$ of $y := |x|$ is given by

```{math}
:enumerated: false

G (s) 
        := \PP\{y\leq s\} 
        = 
        \begin{cases}
            2 F(s) - 1 & \text{ if } s \geq 0
            \\
            0 & \text{ otherwise}
        \end{cases}
```
````


Here **symmetry** of $p$ means that $p(s) = p(-s)$ for all $s \in \RR$. This is equivalent to $F(s) = 1 - F(-s)$, as can be seen by differentiating or integrating. {prf:ref}`ex-avsc` asks you to prove [fact %s](#fa-avsc).

Now take a random variable $x$ with $\lL(x) = F$ and constant $\alpha \in (0,1)$ as given, and consider the $c$ that solves $\PP\{-c \leq x \leq c \}  = 1 - \alpha$.

```{prf:property}

If $\lL(x) = F$, $x$ has a symmetric density and $F$ is strictly increasing, then

$$
c = F^{-1}(1 - \alpha/2) 
        \quad \implies \quad
        \PP\{-c \leq x \leq c\} = 1 - \alpha
$$ (eq-symcdf)
```


To see this, fix $\alpha \in (0,1)$ and let $c$ be as in {eq}`eq-symcdf`. From [fact %s](#fa-avsc), we have

```{math}
:enumerated: false

\PP\{|x| \leq c\} = 2 F(c) - 1 
    = 2 F[F^{-1}(1 - \alpha/2)] - 1 = 1 - \alpha
```

In the case where $F$ is the standard normal CDF $\Phi$, this value $c$ is usually denoted by $z_{\alpha/2}$. We will adopt the same notation:

$$
z_{\alpha/2} := \Phi^{-1}(1 - \alpha/2)
$$ (eq-zalpha2)

See [figure %s](#f-fcca).

```{figure} figures/fcca.svg
:name: f-fcca

Critical values for the standard normal density
```

## Further Reading

The modern axiomatic treatment of probability based on measure theory was formulated in the 1930s by Andrei N. Kolmogorov. A discussion of his life and work can be found in {cite}`kendall1990andrei`. Nice introductory texts on measure-theoretic probability include {cite}`williams1991probability` and {cite}`rosenthal2006first`. Beautiful expositions at a slightly higher level can be found in {cite}`dudley2002real`, {cite}`cinlar`, and {cite}`pollard2002user`.

## Exercises

```{exercise}
:label: ex-dma

Show that if $\fF$ is a $\sigma$-algebra and $A, B \in \fF$, then $A \cap B \in \fF$ also holds.[^4]
```


```{exercise}
:label: ex-bpro2

Prove [fact %s](#fa-bpro2): $\PP(A \cup B) = \PP(A) + \PP(B) - \PP(A \cap B)$ for any $A,B$.[^5]
```


```{exercise}
:label: ex-abc

Given sample space $\Omega := \{1, 2, 3\}$, let $A := \{1\}$, $B := \{2\}$, and $C := \{3\}$. Let $\PP(A) = \PP(B) = 1/3$. Compute $\PP(C)$, $\PP(A \cup B)$, $\PP(A \cap B)$, $\PP(A^c)$, $\PP(A^c \cup B^c)$, and $\PP(A \given B)$. Are $A$ and $C$ independent?
```


```{exercise}
:label: ex-lodi

A die is designed so that the probability of getting face $m$ is $q m$, where $m \in \{1,\ldots,6\}$ and $q$ is a constant. Compute $q$.
```


```{exercise}
:label: ex-cbif

Prove the claims in [fact %s](#fa-cbif).
```


```{exercise}
:label: ex-degprob

Let $\Omega$ be a nonempty finite set, and let $\omega_0$ be a fixed element of $\Omega$. For each $A \subset \Omega$, define $\PP(A) := \1\{\omega_0 \in A\}$. Is $\PP$ a probability on $\Omega$? Why or why not?
```


```{exercise}
:label: ex-aiiself

Let $\Omega$ be any sample space, and let $\PP$ be a probability on the subsets $\fF$. Let $A \in \fF$. Show that if $\PP(A) = 0$ or $\PP(A) = 1$, then $A$ is independent of every other event in $\fF$. Show that if $A$ is independent of itself, then either $\PP(A) = 0$ or $\PP(A) = 1$. Show that if $A$ and $B$ are independent, then $A^c$ and $B^c$ are also independent.
```


```{exercise}
:label: ex-mchip

Let $\PP$ and $\Omega$ be defined as in [example %s](#eg-mchip). Show that $\PP$ is additive, in the sense that if $A$ and $B$ are disjoint events, then $\PP(A \cup B) = \PP(A) + \PP(B)$.
```


```{exercise}
:label: ex-mchip2

Let $\PP$ and $\Omega$ be defined as in [example %s](#eg-mchip). Let $A$ be the event that the first flip is heads, and let $B$ be the event that the second flip is heads. Show that $A$ and $B$ are independent under $\PP$.
```


```{exercise}
:label: ex-firav

Show that when $\Omega$ is finite, each random variable on $\Omega$ is finite.[^6]
```


```{exercise}
:label: ex-fdoex0

Show that any finite random variable can be expressed as

$$
x = \sum_{j=1}^J s_j \1_{A_j}
$$ (eq-simple)
where the scalars $\{s_j\}$ are distinct and the sets $\{A_j\}$ are events that partition $\Omega$.
```


```{exercise}
:label: ex-fdoex

For the generic finite random variable $x$ defined in [exercise %s](#ex-fdoex0), show that $\EE x := \sum_{j=1}^J s_j \PP(A_j)$.
```


```{exercise}
:label: ex-adfrv0

[Fact %s](#fa-exprop) states, among other things, that if $x$ is a random variable and $\alpha$ is a constant, then $\EE [\alpha x] = \alpha \EE[x]$. Prove this for the case where $x$ is a finite random variable using the definition {eq}`eq-expec`.
```


```{exercise}
:label: ex-adfrv

[Fact %s](#fa-exprop) states that, among other things, if $x$ and $y$ are two random variables, then $\EE [x + y] = \EE[x] + \EE[y]$. Prove this for the binary random variables

$$
x(\omega) =  \1_A(\omega)
        \quad \text{and} \quad
        y(\omega) =  \1_B(\omega)
$$ (eq-twvsf)
using the definition {eq}`eq-expec`.
```


```{exercise}
:label: ex-psto

Let $P$ be a discrete distribution with $P(\{s_j\}_{j \geq 1}) = 1$, and let $p_j := P\{s_j\}$. Show that $\sum_{j \geq 1} p_j = 1$.
```


```{exercise}
:label: ex-cdfdrv

Let $P$, $s_j$ and $p_j$ be as in [exercise %s](#ex-psto). Show that if $F$ is the CDF of $P$, then $F(s) = \sum_{j \geq 1} \1\{s_j \leq s\} p_j$ for all $s \in \RR$.
```


```{exercise}
:label: ex-crp0

Let $P$, $s_j$ and $p_j$ be as in [exercise %s](#ex-psto). Show that $\lim_{s \to \infty} F(s) = 1$.
```


```{exercise}
:label: ex-avsc

Prove the claim in [fact %s](#fa-avsc).
```


```{exercise}
:label: ex-ispmf

Let $x$ be a discrete random variable taking values $s_1,\ldots,s_J$, and let $p_j := \PP\{x = s_j\}$. Show that $0 \leq p_j \leq 1$ for each $j$, and $\sum_{j=1}^J p_j = 1$.
```


```{exercise}
:label: ex-cdfsq

Let $\lL(x) = F$ where $F$ is the uniform CDF on $[0,1]$. Give an expression for the CDF $G$ of the random variable $y = x^2$.
```


```{exercise}
:label: ex-ch04_probability-auto-1

Let $F$ be the CDF on $\RR$ defined by $F(s) = e^s/(1 + e^s)$ for $s \in \RR$.

(i) Obtain the quantile function corresponding to $F$.

(ii) Let $\lL(x) = F$ and compute $\PP\{0 \leq x \leq \ln 2\}$.
```


```{exercise}
:label: ex-ch04_probability-auto-2

Let $\lL(y) = F$, where $F$ is a CDF. Show that $F(s) = \EE  \1\{y \leq s\}$ for any $s$.
```


```{exercise}
:label: ex-moncon

Confirm monotonicity of expectations (see [fact %s](#fa-exprop)) for the special case where $x$ and $y$ are the random variables in {eq}`eq-twvsf`.
```


```{exercise}
:label: ex-ofinm

Prove [fact %s](#fa-ofinm) (existence of $k$th moment implies existence of $j$th moment for all $j \leq k$).
```


```{exercise}
:label: ex-vacov

Confirm the expression for variance of linear combinations in [fact %s](#fa-vacov).
```


```{exercise}
:label: ex-ch04_probability-auto-3

Let $x$ and $y$ be scalar random variables. With reference to [fact %s](#fa-acorr), is it true that $\corr[\alpha x, \beta y] = \corr[x, y]$ for *any* constant scalars $\alpha$ and $\beta$? Why or why not?
```


```{exercise}
:label: ex-ctmi

Confirm the solutions in {eq}`eq-slpp`.
```


```{exercise}
:label: ex-vicor

Consider the setting of {ref}`ss-blp`. Let $\alpha^*$, $\beta^*$ and $\ell^*$ be as defined there. Let the prediction error $u$ be defined as $u := y - \ell^*(x)$. Show that

(i) $\EE \ell^*(x) = \EE y$

(ii) $\var[\ell^*(x)] = \corr[x, y]^2 \var[y]$

(iii) $\cov[\ell^*(x), u] = 0$

(iv) $\var[u] = (1 - \corr[x, y]^2) \var[y]$
```


```{exercise}
:label: ex-chein

Prove the Chebyshev inequalities {eq}`eq-chein0` and {eq}`eq-chein`.
```


### Solutions to Selected Exercises

```{solution} ex-dma

Let $\fF$ be a $\sigma$-algebra and let $A, B \in \fF$. The claim is that $A \cap B \in \fF$. To see this, note that by the definition of a $\sigma$-algebra we must have $A^c \in \fF$ and $B^c \in \fF$. Since $\fF$ is closed under the taking of unions, this means that $A^c \cup B^c \in \fF$. Applying [fact %s](#fa-bpos) gives $(A \cap B)^c \in \fF$. But then $A \cap B \in \fF$, as was to be shown.
```


````{solution} ex-bpro2

Pick any sets $A, B \in \fF$. To show that

```{math}
:enumerated: false

\PP(A \cup B) = \PP(A) + \PP(B) - \PP(A \cap B)
```

we start by decomposing $A$ into the union of two disjoint sets: $A = [(A \cup B) \setminus B] \cup (A \cap B)$. Using additivity of $\PP$, we then have

```{math}
:enumerated: false

\PP(A) = \PP[(A \cup B) \setminus B] + \PP(A \cap B)
```

Since $B \subset (A \cup B)$, we can apply part 1 of [fact %s](#fa-bpro) to obtain

```{math}
:enumerated: false

\PP(A) = \PP(A \cup B) - \PP(B) + \PP(A \cap B)
```

Rearranging this expression gives the result that we are seeking.
````


```{solution} ex-abc

First $\PP(C) = 1/3$ as $1=\PP(\Omega)=\PP(A \cup B \cup C)=\PP(A)+\PP(B)+\PP(C)=1/3+1/3+\PP(C)$, and hence $\PP(C)=1/3$. In addition $\PP(A \cup B)=2/3$, $\PP(A\cap B)=0$, $\PP(A^c)=2/3$, $\PP(A^c \cup B^c)=\PP( (A \cap B)^c) = \PP(\Omega) =1$, and $\PP(A \cap C)=0 \neq 1/9 = \PP(A) \PP(C)$. Therefore $A$ is not independent of $C$.
```


````{solution} ex-lodi

When the die is rolled one face must come up, so the sum of the probabilities is one. More formally, letting $\Omega = \{1,\ldots,6\}$ be the sample space, we have

```{math}
:enumerated: false

\PP\{1,\ldots,6\} 
        = \PP \cup_{m=1}^6 \{m\}
        = \sum_{m=1}^6 \PP \{m\}
        = \sum_{m=1}^6 qm = 1
```

Solving the last equality for $q$, we get $q = 1/21$.
````


```{solution} ex-cbif

Let $A_1, \ldots, A_N$ be subsets of $\Omega$. Pick any $\omega \in \Omega$. Regarding part (i), suppose first that $\1_{\cap_{n=1}^N A_n}(\omega) = 1$. Then $\omega \in \cap_{n=1}^N A_n$, and hence $\omega \in A_n$ for all $n$. But then $\prod_{n=1}^N \1_{A_n}(\omega) = 1$. Conversely, if $\1_{\cap_{n=1}^N A_n}(\omega) = 0$, then $\omega$ is not in at least one $A_n$, so $\prod_{n=1}^N \1_{A_n}(\omega) = 0$. This verifies part (i) of [fact %s](#fa-cbif). The proof of (ii) is similar.
```


````{solution} ex-degprob

To show that $\PP$ is a probability on $\Omega$, we need to check that

(i) $\1\{\omega_0 \in A\} \in [0,1]$ for every $A \subset \Omega$.

(ii) $\1\{\omega_0 \in \Omega\} = 1$.

(iii) If $A \cap B = \emptyset$, then $\1\{\omega_0 \in A \cup B\} = \1\{\omega_0 \in A \} + \1\{\omega_0 \in B \}$.

Part (i) is immediate from the definition of an indicator function. Part (ii) holds because $\omega_0 \in \Omega$. Regarding (iii), pick any disjoint $A$ and $B$. If $\omega_0 \in A$, then $\omega_0 \notin B$, and we have

```{math}
:enumerated: false

\1\{\omega_0 \in A \cup B \} 
        = 1
        = \1\{\omega_0 \in A \} + \1\{\omega_0 \in B \}
```

If $\omega_0 \in B$, then $\omega_0 \notin A$, and once again we have

```{math}
:enumerated: false

\1\{\omega_0 \in A \cup B \} 
        = 1
        = \1\{\omega_0 \in A \} + \1\{\omega_0 \in B \}
```

Finally, if $\omega_0$ is in neither $A$ nor $B$, then

```{math}
:enumerated: false

\1\{\omega_0 \in A \cup B \} 
        = 0
        = \1\{\omega_0 \in A \} + \1\{\omega_0 \in B \}
```

We have shown that (i)–(iii) hold, and hence $\PP$ is a probability on $\Omega$.
````


````{solution} ex-aiiself

Suppose that $\PP(A) = 0$ and that $B \in \fF$. We claim that $\PP(A \cap B) = \PP(A)\PP(B)$, or, in this case, $\PP(A \cap B) = 0$. Using nonnegativity and monotonicity of $\PP$ ([fact %s](#fa-bpro)), we obtain

```{math}
:enumerated: false

0 \leq \PP(A \cap B) \leq P(A) = 0
```

Therefore $\PP(A \cap B) = 0$ as claimed.

Now suppose that $\PP(A) = 1$. We claim that $\PP(A \cap B) = \PP(A)\PP(B)$, or, in this case, $\PP(A \cap B) = \PP(B)$. In view of [fact %s](#fa-bpro2), we have

```{math}
:enumerated: false

\PP(A \cap B) = \PP(A) + \PP(B) - \PP(A \cup B)
```

Since $\PP(A) = 1$, it suffices to show that $\PP(A \cup B) = 1$. This last equality is implied by monotonicity of $\PP$, because $1 = \PP(A) \leq \PP(A \cup B) \leq 1$.

Next, suppose that $A$ is independent of itself. In that case we have $\PP(A) = \PP(A \cap A) = \PP(A) \PP(A) = \PP(A)^2$. If $a = a^2$, then either $a = 0$ or $a=1$.

Finally, let $A$ and $B$ be independent. We have

```{math}
:enumerated: false

\PP(A^c \cap B^c) = \PP( (A \cup B)^c) = 1 - \PP(A \cup B)
```

Applying [fact %s](#fa-bpro2) and independence, we can transform the right-hand side to obtain

```{math}
:enumerated: false

\PP(A^c \cap B^c) 
            = (1 - \PP(A)) (1 - \PP(B))
            = \PP(A^c) \PP(B^c)
```

In other words, $A^c$ and $B^c$ are independent.
````


```{solution} ex-mchip

The proof is almost identical to the proof of additivity in [example %s](#eg-ram).
```


```{solution} ex-mchip2

The proof of independence is essentially the same as the proof of independence of $A$ and $B$ in [example %s](#eg-roditwi).
```


```{solution} ex-fdoex0

Let $x$ be finite-valued, with range $\{s_1, \ldots, s_J\}$. Let $A_j := \{x = s_j\}$. By the definition of a function (see {ref}`ss-functions`), each $\omega \in \Omega$ lies in exactly one $A_j$. Hence $A_1, \ldots, A_J$ is a partition of $\Omega$, and $x = \sum_{j=1}^J s_j \1_{A_j}$.
```


```{solution} ex-fdoex

Clearly $x$ is finite-valued, so, by {eq}`eq-expec`, the expectation is given by $\EE x = \sum_{j=1}^J s_j \, \PP\{x = s_j\}$. Hence, to confirm $\EE x := \sum_{j=1}^J s_j \PP(A_j)$, we need only show that $\{x = s_j\} = A_j$. That the last statement is true follows from {eq}`eq-simple` and the fact that $A_1, \ldots, A_J$ is a partition of $\Omega$.
```


````{solution} ex-adfrv0

Let $x$ be finite with range $\{s_1, \ldots, s_J\}$ and let $\alpha$ be any constant. We claim that $\EE[\alpha x] = \alpha \EE[x]$. If $\alpha=0$ the result is trivial, as both sides are clearly zero. If not then $\alpha x = \alpha s_j$ if and only if $x = s_j$. Hence, applying {eq}`eq-expec`,

```{math}
:enumerated: false

\EE[ \alpha x ]  
        = \sum_{j=1}^J \alpha s_j \PP\{x = s_j\}
        = \alpha \left[ \sum_{j=1}^J s_j \PP\{x = s_j\} \right]
        = \alpha \EE [x]
        
```
````


````{solution} ex-adfrv

Consider the sum $x + y$. By this, we mean the random variable $(x+y)(\omega) := x(\omega) + y(\omega)$. We claim that $\EE [ x + y] = \EE[x] + \EE[y]$. To see this, note first that

```{math}
:enumerated: false

(x + y)(\omega) 
        =  \1_{A \setminus B}(\omega) 
            +  \1_{B \setminus A}(\omega)
            + 2 \1_{ A \cap B }(\omega)
```

(To check this, just go through the different cases for $\omega$, and verify that the right-hand side of this expression agrees with $x(\omega) + y(\omega)$. Sketching a Venn diagram will help.) Therefore, by the definition of expectation,

$$
\EE[x + y]
        =  \PP (A \setminus B) +  \PP(B \setminus A)
            + 2 \PP(A \cap B)
$$ (eq-adexf)
Now observe that $A = (A \setminus B) \cup (A \cap B)$, and hence, by disjointness,

```{math}
:enumerated: false

\EE[x] :=  \PP(A) =  \PP(A \setminus B) + \PP(A \cap B)
```

Performing a similar calculation with $y$ produces

```{math}
:enumerated: false

\EE[y] := \PP(B) = \PP(B \setminus A) + \PP(A \cap B)
```

Summing these equations produces the value on the right-hand side of {eq}`eq-adexf`, thereby confirming $\EE[x + y] = \EE[x] + \EE[y]$.
````


```{solution} ex-psto

Let $P(\{s_j\}_{j \geq 1}) = 1$ and $p_j := P\{s_j\}$ as stated in the exercise. The set $\{s_j\}_{j \geq 1}$ can be written as $\cup_{j=1}^{\infty} \{s_j\}$. Since $P$ is a distribution it is countably additive and hence $1= P(\cup_{j=1}^{\infty} \{s_j\}) = \sum_{j=1}^\infty P\{s_j\} = \sum_{j=1}^\infty p_j$.
```


````{solution} ex-cdfdrv

Using additivity and the fact that $P\{s_j\}_{j \geq 1} = 1$, we have

```{math}
:enumerated: false

F(s)
        = P(-\infty, s]
        = P \left( 
            \bigcup_{j \st s_j \leq s} \{s_j\}
              \right)
        = \sum_{j \st s_j \leq s} P \{s_j\} 
        = \sum_{j=1}^J \1\{s_j \leq s\} p_j
        
```
````


````{solution} ex-crp0

We are assuming that $x$ has finite range, and hence takes only finitely many different values. Let $m$ be the largest such value. For this $m$, we have

```{math}
:enumerated: false

\lim_{s \to \infty} F_x(s) 
        \geq F_x(m) 
        = \PP\setntn{\omega \in \Omega}{x(\omega) \leq m}
        = \PP(\Omega) = 1
```

(The inequality is due to the fact that $F_x$ is increasing.) On the other hand,

```{math}
:enumerated: false

\lim_{s \to \infty} F_x(s) 
        = \lim_{s \to \infty} \PP\{x \leq s\}
        \leq \lim_{s \to \infty} \PP(\Omega) = 1
```

From these two inequalities we get $1 \leq \lim_{s \to \infty} F_x(s) \leq 1$, which is equivalent to $\lim_{s \to \infty} F_x(s) = 1$.
````


````{solution} ex-avsc

Fix $s \geq 0$. Using additivity over disjoint sets, we have

```{math}
:enumerated: false

F_{|x|} (s) 
        := \PP\{|x| \leq s\}
        = \PP\{-s \leq x \leq s\} 
        = \PP\{x = -s\} + \PP\{-s < x \leq s\}
```

By assumption, $\PP\{x = -s\} = 0$. Applying [fact %s](#fa-axb) then yields

```{math}
:enumerated: false

F_{|x|} (s) = \PP\{-s < x \leq s\} = F(s) - F(-s)
```

As discussed after [fact %s](#fa-avsc), symmetry is equivalent to $F(-s) = 1 - F(s)$. Plugging this in to the right-hand side of the previous expression gives $F_{|x|} (s) = 2 F(s) - 1$.
````


```{solution} ex-ispmf

That $0 \leq p_j \leq 1$ for each $j$ follows immediately from the definition of $\PP$. In addition, using additivity of $\PP$, we have

$$
\sum_{j=1}^J p_j
        = \sum_{j=1}^J \PP\{x = s_j\}
        = \PP \cup_{j=1}^J \{x = s_j\}
        = \PP (\Omega) = 1
$$ (eq-ispmf2)
(We are using the fact that the sets $\{x = s_j\}$ are disjoint. Why is this always true? Look at the definition of a function given in {ref}`ss-functions`.)
```


````{solution} ex-cdfsq

Evidently $G(s) = 0$ when $s < 0$. For $s \geq 0$ we have

```{math}
:enumerated: false

\PP\{x^2 \leq s\}
        = \PP\{|x| \leq \sqrt{s}\}
        = \PP\{x \leq \sqrt{s}\}
        = F(\sqrt{s})
```

Thus, $G(s) = F(\sqrt{s}) \1\{s \geq 0\}$.
````


````{solution} ex-moncon

If $x(\omega) := \1\{\omega \in A\} \leq \1\{\omega \in B\} =: y(\omega)$ for any $\omega \in \Omega$, then $A \subset B$. (If $\omega \in A$, then $x(\omega) = 1$. Since $x(\omega) \leq y(\omega) \leq 1$, we have $y(\omega) = 1$, and hence $\omega \in B$.) Using [fact %s](#fa-exprop) and monotonicity of $\PP$, we then have

```{math}
:enumerated: false

\EE x 
        = \EE \1\{\omega \in A\} 
        = \PP(A)
        \leq \PP(B)
        = \EE \1\{\omega \in B\} 
        = \EE y
        
```
````


```{solution} ex-ofinm

Let $a$ be any nonnegative number, and let $j \leq k$. If $a \geq 1$, then $a^j \leq a^k$. If $a < 1$, then $a^j \leq 1$. Thus, for any $a \geq 0$, we have $a^j \leq a^k + 1$, and for any random variable $x$ we have $|x|^j \leq |x|^k + 1$. Using linearity and monotonicity of expectations and $\EE 1 = 1$, we then have $\EE |x|^j \leq \EE |x^k| + 1$. Hence the $j$th moment exists whenever the $k$th moment exists.
```


````{solution} ex-vacov

We have

```{math}
:enumerated: false

\begin{aligned}
\var \left[ \sum_{n=1}^N \alpha_n x_n \right] 
    & =
      \EE \left[ \left( \sum_{n=1}^N \alpha_n x_n 
              - \EE \left[
      \sum_{n=1}^N \alpha_n x_n 
                  \right] \right)^2 \right]  
      \\
    &= \EE \left[ 
        \left( 
            \sum_{n=1}^N \alpha_n (x_n - \e{x_n})
        \right)^2 
          \right]  
     \\
    &= \EE \left[ 
        \sum_{n=1}^N \alpha_n^2 (x_n - \e{x_n})^2 + 2
        \sum_{n < m} \alpha_n \alpha_m (x_n - \e{x_n})(x_m - \e{x_m}) 
          \right]
       \\
    &= \sum_{n=1}^N \alpha_n^2 \var [x_n] + \sum_{n < m} \alpha_n
      \alpha_m \cov [x_n, x_m]
\end{aligned}
```

as required.
````


````{solution} ex-chein

Pick any nonnegative random variable $x$ and $\delta > 0$. By considering what happens at an arbitrary $\omega \in \Omega$, you should be able to convince yourself that

```{math}
:enumerated: false

x 
        = \1\{x \geq \delta\} x + \1\{x < \delta\} x
        \geq \1\{x \geq \delta\} \delta
```

Using the results in [fact %s](#fa-exprop) and rearranging gives {eq}`eq-chein0`. Regarding {eq}`eq-chein`, observe that

```{math}
:enumerated: false

x^2 
        = \1\{|x| \geq \delta\} x^2 + \1\{|x| < \delta\} x^2
                    \geq \1\{|x| \geq \delta\} \delta^2
```

Proceeding as before leads to {eq}`eq-chein`.
````


[^1]: If $p$ is a nonnegative number satisfying $p \geq 0$ and $p \leq \epsilon$ for all $\epsilon > 0$, then $p=0$ must hold.

[^2]: Is this true? If $\omega = \omega_0$ is an infinite sequence containing only zeros, then $\setntn{n \in \NN}{b_n = 1} = \emptyset$. What is $x(\omega_0)$ in this case? The convention is $\min \emptyset = \infty$, which is reasonable, but now $x$ can take the value $\infty$. However, in most applications this event has probability zero. We could set $x(\omega_0) = 0$ without changing anything significant. Now we're back to a well-defined function from $\Omega$ to $\RR$.

[^3]: When we want to combine $x$ with other random variables we'll need their joint distribution. See {ref}`ss-drvecs`.

[^4]: Hint: See the rules concerning set-theoretic operations in [fact %s](#fa-bpos).

[^5]: Hint: Sketching the Venn diagram, convince yourself that $A = [(A \cup B) \setminus B] \cup (A \cap B)$. Finish the proof using the definition of a probability and [fact %s](#fa-bpro).

[^6]: Hint: Review the definition of a function in {ref}`ss-functions`.
