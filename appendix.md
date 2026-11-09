---
title: "Appendix"
---

(s-setfun)=
## Sets
In the text we often refer to the **real numbers**. This set is denoted by $\RR$, and we understand it to contain "all the numbers." $\RR$ can be visualized as the "continuous\" real line:

![](figures/xfig_real_line.svg)

The term "real" sounds odd but it means "not imaginary.\"

$\RR$ is an example of a **set**. A set is a collection of distinct objects viewed as a whole. In the case of $\RR$, the objects are numbers. Other examples of sets are the set of all rectangles in the plane, or the set of all monkeys in Japan.

If $A$ is a set, then the statement $x \in A$ means that $x$ is contained in (alternatively, is an element of) $A$. If $B$ is another set, then $A \subset B$ means that any element of $A$ is also an element of $B$, and we say that $A$ is a **subset** of $B$. The statement $A = B$ means that $A$ and $B$ contain the same elements (or, equivalently, $A \subset B$ and $B \subset A$). For example, if $\mathbbm I$ is the irrational numbers,[^1] then $\mathbbm I \subset \RR$. Also, $0 \in \RR$, $\pi \in \RR$, $-3 \in \RR$, $e \in \RR$, and so on.

```{figure} figures/sets1.svg
:name: f-sets1

Sets $A$ and $B$ in $S$
```

Commonly used subsets of $\RR$ include the intervals. For arbitrary $a$ and $b$ in $\RR$, the **open interval** $(a, b)$ is defined as

```{math}
:enumerated: false

(a, b) := \setntn{x \in \RR}{a < x < b}
```

while the **closed interval** $[a, b]$ is defined as

```{math}
:enumerated: false

[a, b] := \setntn{x \in \RR}{a \leq x \leq b}
```

We also use half open intervals such as $[a, b) := \setntn{x \in \RR}{a \leq x < b}$, half lines such as $(-\infty, b) = \setntn{x \in \RR}{x < b}$, and so on.

Let $S$ be a set and let $A$ and $B$ be two subsets of $S$, as illustrated in [figure %s](#f-sets1).

The **union** of $A$ and $B$ is the set of elements of $S$ that are in $A$ or $B$ or both:

```{math}
:enumerated: false

A \cup B := 
        \setntn{x \in S}{x \in A \text{ or } x \in B}
```

Here and elsewhere, "or" is used in the mathematical sense. It means "and/or".

The **intersection** of $A$ and $B$ is the set of all elements of $S$ that are in both $A$ and $B$:

```{math}
:enumerated: false

A \cap B := 
        \setntn{x \in S}{x \in A \text{ and } x \in B}
```

The set $A \setminus B$ is all points in $A$ that are not points in $B$:

```{math}
:enumerated: false

A \setminus B := 
        \setntn{x \in S}{x \in A \text{ and } x \notin B}
```

The **complement** of $A$ is the set of elements of $S$ that are not contained in $A$:

```{math}
:enumerated: false

A^c := S \setminus A := \setntn{x \in S}{x \notin A}
```

{numref}`f-allsets` illustrate these definitions.

```{figure} figures/allsets.svg
:name: f-allsets

Union, intersection, and complement
```

For example, since $\RR$ consists of the irrationals $\mathbbm I$ and the rationals $\QQ$, we have

```{math}
:enumerated: false

\QQ \subset \RR, \;
    \mathbbm I \subset \RR, \;
    \QQ \cup \mathbbm I = \RR, \;
    \QQ^c = \mathbbm I, \;
    \text{etc.}
```

Also

```{math}
:enumerated: false

\NN := \{1, 2, 3,\ldots\} \subset \QQ \subset \RR
```

The **empty set** is the set containing no elements. It is denoted by $\emptyset$. If the intersection of $A$ and $B$ equals $\emptyset$, then $A$ and $B$ are said to be **disjoint**.

The next fact lists some well-known rules for set-theoretic operations.

```{prf:property}
:label: fa-bpos

Let $A$ and $B$ be subsets of $S$. The following statements are true:

(i) $A \cup B = B \cup A$ and $A \cap B = B \cap A$.

(ii) $(A \cup B)^c = B^c \cap A^c$ and $(A \cap B)^c = B^c \cup A^c$.

(iii) $A \setminus B = A \cap B^c$.

(iv) $(A^c)^c = A$.
```


If $\{A_{\alpha}\}_{\alpha \in \Lambda}$ is an infinite collection of subsets of $S$, the union $\cup_{\alpha \in \Lambda} A_\alpha$ is defined as all $x \in S$ such that $x \in A_\alpha$ for some $\alpha \in \Lambda$. The intersection $\cap_{\alpha \in \Lambda} A_\alpha$ is all $x \in S$ that lie in $A_\alpha$ for all $\alpha \in \Lambda$.

```{prf:example}

If $A_n = (-n, n)$ for each $n \in \NN$, then $\cup_{n \in \NN} A_n = \RR$.
```


### Cartesian Products

Collections with natural order are frequently written as **tuples**. These are finite sequences of terms, usually denoted using notation such as $(a_1, a_2)$ or $(x_1, x_2, x_3)$. Collections of tuples arise naturally in Cartesian products. The **Cartesian product** of set $A_1, \ldots, A_N$ is the set

```{math}
:enumerated: false

A_1 \times \cdots \times A_N
    := \setntn{(a_1, \ldots, a_N)}{a_n \in A_n
    \text{ for } n =1, \ldots, N}
```

The **vector space $\RR^N$** is one example of a Cartesian product, defined as

```{math}
:enumerated: false

\RR^N = \RR \times \cdots \times \RR \quad (N \text{ times})
```

Thus $\RR^N$ is all tuples $(x_1, \ldots, x_N)$ with $x_n \in \RR$ for each $n$.

A useful class of subsets of $\RR^N$ is those that are formed as Cartesian products of intervals. These are called **rectangles**. An example is

```{math}
:label: eq-rectang
:enumerated: false

I  = \times_{n=1}^N \, (a_n, b_n]
      := \setntn{(x_1, \ldots, x_N) \in \RR^N}{a_n < x_n \leq b_n \text{ for } n=1,\ldots, N}
```

where $a_1, b_1, \ldots, a_N, b_N$ are given numbers. The rectangle $(2, 5] \times (1, 2]$ is shown in [figure %s](#f-cart_prod).

```{figure} figures/cart_prod.svg
:name: f-cart_prod

Rectangle $(2, 5] \times (1, 2]$ in $\RR^2$
```

(ss-functions)=
## Functions
There are two fundamental primitives in mathematics: sets and functions.[^2] A **function** $f$ from set $A$ to set $B$ is a rule that associates to each element of $A$ a unique element of $B$. The symbols $f \colon A \to B$ means that $f$ is a function from $A$ to $B$. The set $A$ is called the **domain** of $f$, while $B$ is called the **codomain**. {numref}`f-function` illustrates.

```{figure} figures/tikz_function.svg
:name: f-function

A function from $A$ to $B$
```

If $f \colon A \to B$ and $g \colon B \to C$ then the function $h \colon A \to C$ defined by $h(a) = g(f(a))$ is called the **composition** of $g$ and $f$, and written as $g \circ f$.

```{figure} figures/tikz_composition.svg
:name: f-composition

The composition $g \circ f$
```

````{prf:example}

If $f \colon \RR \to \RR$ is defined by $f(x) = e^x :=: \exp(x)$ and $g \colon \RR \to \RR$ is defined by $g(x) = x^2$, then

```{math}
:enumerated: false

(g \circ f)(x) = g(f(x)) = \exp(2x)
```
````


If $a$ and $b$ are points such that $f(a) = b$, then $a$ is called a **preimage** of $b$. The set of preimages of $b$ can be empty, a singleton or contain multiple values. The set of points

```{math}
:enumerated: false

\setntn{b \in B}{f(a) = b \text{ for some } a \in A}
```

is called the **range** of $f$, and written as $\range f$. Thus $b$ is in the range of $f$ if it has at least one preimage. {numref}`f-range` depicts the range of a given function $f \colon A \to B$.

```{figure} figures/xfig_range.svg
:name: f-range

Range of $f \colon A \to B$
```

A function $f \colon A \to B$ is called **one-to-one** if distinct elements of $A$ are always mapped into distinct elements of $B$. That is, if

```{math}
:enumerated: false

a \in A, \; a' \in A,  \text{ and } a \neq a' 
   \implies f(a) \neq f(a')
```

or, equivalently, if $f(a) = f(a') \implies a = a'$. Evidently $f \colon A \to B$ is one-to-one if and only if each element of $B$ has at most one preimage under $f$. If a function $f \colon A \to B$ and $f$ is one-to-one, then the equation $f(x) = b$ has at most one solution.

A function $f \colon A \to B$ is called **onto** if $\range f = B$; this is if every $b \in B$ has at least one preimage. If $f \colon A \to B$ and $f$ is onto, then the equation $f(x) = b$ has at least one solution.

Finally, $f \colon A \to B$ is called a **bijection** if it is both one-to-one and onto. Bijections are sometimes called one-to-one correspondences.

```{prf:example}

The function $f(x) = 2x$ is a bijection from $\RR$ to $\RR$.
```


If $f$ is a bijection, then the equation $f(x) = b$ always has exactly one solution. We can therefore define the **inverse function** to $f$, denoted by $f^{-1}$, as the map from $B \to A$ such that $f^{-1}(b)$ is the unique $a$ with $f(a) = b$.

```{prf:example}

The function $f \colon A \to B$ in [figure %s](#f-function) is a bijection, with $f^{-1} \colon B \to A$ mapping $b_1$ to $a_1$ and $b_2$ to $a_2$.
```


```{prf:example}

If $f \colon \RR \to (0, \infty)$ is defined by $f(x) = \exp(x) := e^x$ and $\phi \colon (0, \infty) \to \RR$ is defined by $\phi(x) = \log(x)$, then $\phi = f^{-1}$ because, for any $a \in \RR$, we have $\phi(f(a)) = \log(\exp(a)) = a$.
```


```{prf:property}
:label: fa-bij

Let $f \colon A \to B$ and $g \colon B \to C$ be bijections.

(i) $f^{-1}$ is a bijection and its inverse is $f$.

(ii) $f^{-1}(f(a)) = a$ for all $a \in A$.

(iii) $f(f^{-1}(b)) = b$ for all $b \in B$.

(iv) $g \circ f$ is a bijection from $A$ to $C$ and $(g \circ f)^{-1} = f^{-1} \circ g^{-1}$.
```


Part (iv) of [fact %s](#fa-bij) is illustrated in [figure %s](#f-bij_inv).

```{figure} figures/tikz_bij_inv.svg
:name: f-bij_inv

The inverse of $g \circ f$ is $f^{-1} \circ g^{-1}$
```

### Preimage of Sets

Given any function $f$ from $A$ to $B$ and any subset $E$ of $B$, the symbol $f^{-1}(E)$ denotes the set of points that are mapped into $E$ under $f$. That is,

```{math}
:enumerated: false

f^{-1}(E) := \setntn{x \in A}{f(x) \in E}
```

The set $f^{-1}(E)$ is called the **preimage** of $E$ under $f$.

```{prf:example}

If $f(x) = x^2$ and $E = [0, 1]$, then $f^{-1}(E) = [-1, 1]$.
```


Preimages of sets interact nicely with set-theoretic operations. For example, given any $E, F \subset B$, we have

(i) $f^{-1}(E \cup F) = f^{-1}(E) \cup f^{-1}(F)$,

(ii) $f^{-1}(E \cap F) = f^{-1}(E) \cap f^{-1}(F)$, and

(iii) $f^{-1}(E^c) = [f^{-1}(E)]^c$.

The equalities in (i) and (ii) also extend to infinite intersections and unions.

(ss-card)=
## Cardinality and Measure
For finite sets we have a fairly clear idea of what "the same size\" means. For infinite sets things get tricker. It is nonetheless valuable to distinguish infinite sets according to size. To this end, when a bijection exists between sets $A$ and $B$, they are said to have the **same cardinality**, and we write $|A| = |B|$.

```{prf:example}

If $A = \{a_1, \ldots, a_n\}$ and $B = \{b_1, \ldots, b_n\}$ then $|A| = |B|$ because $f(a_n) = b_n$ is a bijection. (Remember that elements of a set are by definition distinct.)
```


More generally, we have the following result:

```{prf:property}

If $A$ and $B$ are finite then $|A| = |B|$ if and only if $A$ and $B$ have the same number of elements.
```


This fact tells us that classifying sets by cardinality generalizes the notion of classifying sets by number of elements.

```{prf:property}

For any sets $A, B, C$ we have $|A| = |B|$ and $|B| = |C|$ $\implies$ $|A| = |C|$.
```


In other words, cardinality is transitive. This follows from the fact that compositions of bijections are bijections ([fact %s](#fa-bij)).

A nonempty set $A$ is called **finite** if $|A| = |\{1, 2, \ldots, n\}|$ for some $n \in \NN$. Otherwise, $A$ is called **infinite**.

Sets that either are finite or have the same cardinality as $\NN$ are called **countable**.

```{prf:example}

$E := \{2, 4, \ldots\}$ is countable because it can be put in one-to-one correspondence with $\NN$ via $f \colon E \ni m \mapsto m/2 \in \NN$.
```


In essence, sets formed by sequences are countable.

**Uncountable** sets are sets that are not countable. For example, $\RR$ is uncountable, as is any set that contains a nonempty open interval $(a, b)$. Put simply, you cannot write down a sequence that will exhaust the elements of these sets. The story is similar for $\RR^N$. (Such sets are said to have **the power of the continuum**.) For a beautiful proof that $\RR$ is uncountable, look up the diagonalization argument of Georg Cantor.

(ss-measure)=
### Lebesgue Measure and Sets of Measure Zero
A set $B \subset \RR$ is said to have **Lebesgue measure zero** in $\RR$ if, given any $\epsilon > 0$, there exists a sequence of intervals $I_n = (a_n, b_n)$ such that $B$ is contained in the union of these intervals and $\sum_{n=1}^{\infty} (b_n - a_n) \leq \epsilon$. The definition extends to $\RR^N$ by replacing intervals with rectangles (see page  for a definition of the latter).

For example, any countable set has measure zero. To see this in $\RR$, let $C := \cup_{n=1}^{\infty} \{x_n\}$ with $x_n \in \RR$. Fix $\epsilon >0$. To each $x_n$ we assign an interval $(a_n, b_n)$ containing $x_n$ and satisfying $b_n -a_n < \epsilon 2^{-n}$. The union of intervals contains $C$, and, in addition, the sum of their lengths is $\epsilon \sum_{n=1}^\infty 2^{-n} = \epsilon$.

Uncountable sets can have measure zero too. For example, every linear subspace $S$ with $\dim S < N$ has measure zero in $\RR^N$. A more classical example of an uncountable set with measure zero is the Cantor set. Details can be found in many texts on analysis.

The concept of measure zero sets forms part of a branch of mathematics called **measure theory**. We have used aspects of measure theory in this text—for example, in our discussion of $\sigma$-algebras in {ref}`ss-soev` and probabilities in {ref}`ss-probs`. The basic problem of measure theory is how to provide a well-defined notion of "measure\" or "size\" to sets. This is easy if we just want to assign area to a circle (area $= \pi r^2$) or rectangle, or volume to a sphere. But what about an arbitrary subset of $\RR^N$? Is there a well-defined and widely accepted notion of size for such a set?

The answer to this question is affirmative. The so-called **Lebesgue outer measure** is a map, typically denoted $\lambda^*$, that gives a natural notion of size to any set. It works by approximating arbitrary sets by simpler sets, like rectangles, and agrees with standard notions of area or volume for common geometric objects. The only problem is that it isn't always additive, in the sense that

$$
\lambda^*(A \cup B) = \lambda^*(A) + \lambda^*(B)
$$ (eq-addlo)

can fail. To get around this problem, we exclude any sets that cause {eq}`eq-addlo` to fail. That is, we define the **Lebesgue measurable sets** to be those set on which additivity holds. This causes no difficulty since the discarded sets never come up in day-to-day analysis. In fact we usually focus on a smaller class of well-behaved sets, called the **Borel sets**. These were introduced in {ref}`ss-soev`.

(s-opt)=
## Real-Valued Functions
For any set $A$, if $f \colon A \to \RR$ then $f$ is called a **real-valued function**. If $f$ and $g$ are real-valued functions, then $f + g$ is defined by $(f + g)(x) = f(x) + g(x)$, while $\alpha f$ is defined by $(\alpha f)(x) = \alpha f(x)$. A **maximizer** of $f$ on $A$ is a point $a^* \in A$ such that $f(a^*) \geq f(a)$ for all $a \in A$. The value $f(a^*)$ is called the **maximum** of $f$ on $A$. A **minimizer** of $f$ on $A$ is a point $b^* \in A$ such that $f(b^*) \leq f(a)$ for all $a \in A$. The value $f(b^*)$ is called the **minimum** of $f$ on $A$.

Monotone increasing transformations of functions do not affect maximizers. To see this, let $f \colon \RR \to \RR$ and let $h$ be **monotone increasing**, in the sense that if $x \leq x'$, then $h(x) \leq h(x')$, and let $g$ be the function defined by $g(a) = h(f(a))$. Any maximizer of $f$ on $A$ is also a maximizer of $g$ on $A$. To see why, let $a \in A$. Since $a^*$ is a maximizer of $f$, it must be the case that $f(a) \leq f(a^*)$. Since $h$ is monotone increasing, this implies that $h(f(a)) \leq h(f(a^*))$. Given that $a$ was chosen arbitrarily, we have now shown that $g(a^*) \geq g(a)$ for all $a \in A$. In other words, $a^*$ is a maximizer of $g$ on $A$.

```{figure} figures/max_preserved.svg
Increasing transform $h(x) = \exp(x/2)$ preserves the maximizer
```

(sss-sai)=
### Sup and Inf
Maximizers and minimizers don't always exist. To illustrate, consider $f \colon (0, 1) \to (0, 1)$ defined by $f(x) = x$. The function $f$ has neither maximizer nor minimizer on $(0, 1)$. For example, no $a \in (0,1)$ is a minimizer because $b := a/2$ is in $(0,1)$ and $f(b) < f(a)$.

To negotiate this kind of problem, we can use the notion of supremum instead. If $A$ is a set, then the **supremum** $s := \sup A$ is the unique $s \in \RR$ such that (i) $a \leq s$ for every $a \in A$ and (ii) there exists a sequence $\{x_n\} \subset A$ such that $x_n \to s$. For example, $1$ is the supremum of both $(0,1)$ and $[0,1]$. The **infimum** $i := \inf A$ is the unique $i \in \RR$ such that (i) $a \geq i$ for every $a \in A$ and (ii) there exists a sequence $\{x_n\} \subset A$ such that $x_n \to i$. For example, $0$ is the infimum of both $(0,1)$ and $[0,1]$.

One can show that the supremum and infimum of any bounded set $A$ exist. The boundedness requirement can be dropped when the values $-\infty$ and $\infty$ are admitted as a possible infima and suprema.

Returning to our original example with $f(x) = x$, while $\max_{x \in (0,1)} f(x)$ is not well-defined, $\sup_{x \in (0,1)} f(x) := \sup \setntn{f(x)}{x \in (0,1)} = \sup (0, 1) = 1$.

[^1]: The **irrationals** are those numbers such as $\pi$ and $\sqrt 2$ that cannot be expressed as fractions of whole numbers.

[^2]: Actually functions can be represented as a special type of set containing ordered pairs, and hence in pure mathematics cannot claim to be as foundational as sets. However, for our purposes, it will be fine to think of a function as a primitive in its own right.
