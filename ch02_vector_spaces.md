---
title: "Vector Spaces"
label: c-la
---

Our first technical topic for this book is linear algebra, which is one of the foundation stones of applied mathematics in general, and econometrics and statistics in particular. Data ordered by observation are naturally stored in vectors. Related vectors are naturally grouped into matrices. Once our data are organized this way, we need to perform basic arithmetic operations or solve equations or quadratic minimization problems. In this chapter we cover the foundations of vector operations used in linear algebra. As we'll see, the conceptual aspects of linear algebra are clearest if we begin by stripping away details such as matrices and look instead at linear operations from a more abstract perspective.

(s-vm)=
## Vectors and Vector Space
Let's begin with vector space and basic vector operations.

(ss-v)=
### Vectors
For arbitrary $N \in \NN$, the symbol $\RR^N$ represents the set of all $N$-vectors, or vectors of length $N$. A typical element has the form

```{math}
:enumerated: false

\boldx = 
    \left(
    \begin{array}{c}
        x_1 \\
        x_2 \\
        \vdots \\
        x_N
    \end{array}
    \right)
    \quad
    \text{ where } x_n \in \RR \text{ for each } n
```

($\RR = \RR^1$ represents the set of all real numbers, which is the union of the rational and irrational numbers.) While $\boldx$ has been written vertically, as a column of numbers, we could also write it horizontally as $\boldx = (x_1, \ldots, x_N)$. For now vectors are just a sequences of numbers, and it makes no difference whether we write them vertically or horizontally. (Only when we work with matrix multiplication does it become necessary to distinguish between column and row vectors.)

The vector of ones will be denoted $\boldone$ and the vector of zeros will be denoted $\boldzero$:

```{math}
:enumerated: false

\boldone := 
    \left(
    \begin{array}{c}
        1 \\
        \vdots \\
        1
    \end{array}
    \right),
    \quad \quad
    \boldzero := 
    \left(
    \begin{array}{c}
        0 \\
        \vdots \\
        0
    \end{array}
    \right)
```

Although vectors are infinitesimal points in $\RR^N$, they are often represented visually as arrows from the origin to the point itself. {numref}`f-vecs` gives an illustration for the case $N=2$.

```{figure} figures/vecs.svg
:name: f-vecs

Three vectors in $\RR^2$
```

In vector space theory there are two fundamental algebraic operations: vector addition and scalar multiplication. For $\boldx, \boldy \in \RR^N$, their **vector sum** is

```{math}
:enumerated: false

\boldx + \boldy 
    = 
    \left(
    \begin{array}{c}
        x_1 \\
        \vdots \\
        x_N
    \end{array}
    \right)
    +
    \left(
    \begin{array}{c}
         y_1 \\
        \vdots \\
         y_N
    \end{array}
    \right)
    :=
    \left(
    \begin{array}{c}
        x_1 + y_1 \\
        \vdots \\
        x_N + y_N
    \end{array}
    \right)
```

If $\alpha \in \RR$, then the **scalar product** of $\alpha$ and $\boldx$ is defined to be

```{math}
:enumerated: false

\alpha \boldx 
    =
    \alpha 
    \left(
    \begin{array}{c}
        x_1 \\
        \vdots \\
        x_N
    \end{array}
    \right)
    :=
    \left(
    \begin{array}{c}
        \alpha x_1 \\
        \vdots \\
        \alpha x_N
    \end{array}
    \right)
```

Thus addition and scalar multiplication are defined in terms of ordinary addition and multiplication in $\RR$, computed element by element, adding and multiplying respectively.[^1] {numref}`f-vec_add` and {numref}`f-vec_scalar` show examples of vector addition and scalar multiplication in the case $N=2$.

```{figure} figures/tikz_vec_add.svg
:name: f-vec_add

Vector addition
```

```{figure} figures/tikz_vec_scalar.svg
:name: f-vec_scalar

Scalar multiplication
```

Subtraction of two vectors is performed element by element, just like addition. Subtraction is not a primitive operation because the definition can be given in terms of addition and scalar multiplication: $\boldx - \boldy := \boldx + (-1) \boldy$. An illustration of subtraction is given in [figure %s](#f-vec_minus). One way to remember this operation is to draw a line from $\boldy$ to $\boldx$ and then shift it to the origin.

```{figure} figures/tikz_vec_minus.svg
:name: f-vec_minus

Difference between vectors
```

The **inner product** of two vectors $\boldx$ and $\boldy$ in $\RR^N$ is denoted by $\inner{\boldx,  \boldy}$, and defined as the sum of the products of their elements:

$$
\inner{\boldx,  \boldy} := \sum_{n=1}^N x_n y_n
$$ (eq-dipv)

```{prf:property}
:label: fa-innpp

For any $\alpha, \beta \in \RR$ and any $\boldx, \boldy, \boldz \in \RR^N$, the following statements are true:

(i) $\inner{\boldx, \boldy} = \inner{\boldy, \boldx}$,

(ii) $\inner{\alpha \boldx, \beta \boldy} =  \alpha \beta \inner{\boldx, \boldy}$, and

(iii) $\inner{\boldx, \alpha \boldy + \beta \boldz} =  \alpha \inner{\boldx, \boldy} + \beta \inner{\boldx, \boldz}$.
```


These properties are easy to check from {eq}`eq-dipv`. For example, regarding the second equality, pick any $\alpha, \beta \in \RR$ and any $\boldx, \boldy \in \RR^N$. By the definitions of scalar multiplication and inner product respectively, we have

```{math}
:enumerated: false

\inner{\alpha \boldx, \beta \boldy}
    = \sum_{n=1}^N \alpha x_n \beta y_n
    = \alpha \beta \sum_{n=1}^N x_n y_n 
    = \alpha \beta \inner{\boldx, \boldy}
```

The (Euclidean) **norm** of a vector $\boldx \in \RR^N$ is defined as

$$
\| \boldx \| := \sqrt{\inner{\boldx, \boldx}}
$$ (eq-eucnorm)

and represents the length of the vector $\boldx$. (In the arrow representation of vectors in [figures %s](#f-vec_add)–[%s](#f-vec_minus), the norm of the vector is equal to the length of the arrow.)

```{prf:property}
:label: fa-profno

For any $\alpha \in \RR$ and any $\boldx, \boldy \in \RR^N$, the following statements are true:

(i) $\| \boldx \| \geq 0$ and $\| \boldx \| = 0$ if and only if $\boldx = \boldzero$,

(ii) $\| \alpha \boldx \| = |\alpha| \| \boldx \|$,

(iii) $\| \boldx + \boldy \| \leq  \| \boldx \| + \| \boldy \|$, and

(iv) $| \! \inner{\boldx, \boldy} \! | \leq  \| \boldx \| \| \boldy \|$.
```


Properties (i) and (ii) you can verify yourself without difficulty. Proofs for (iii) and (iv) are a bit harder. Property (iii) is called the **triangle inequality**, while (iv) is called the **Cauchy–Schwarz inequality**. The proof of the Cauchy–Schwarz inequality is given as a solved exercise after we've built up some more tools (see [ex. %s](#ex-pcsin)). If you're prepared to accept the Cauchy–Schwarz inequality for now, then the triangle inequality follows, because, by the properties of the inner product given in [fact %s](#fa-innpp),

```{math}
:label: eq-tifcs
:enumerated: false

\| \boldx + \boldy \|^2
    = \inner{\boldx + \boldy, \boldx + \boldy}
    = \inner{\boldx, \boldx} + 2 \inner{\boldx, \boldy} +  \inner{\boldy,
    \boldy}
    \leq \inner{\boldx, \boldx} + 2 | \! \inner{\boldx, \boldy} \! | +  \inner{\boldy,
    \boldy}
```

Applying the Cauchy–Schwarz inequality leads to $\| \boldx + \boldy \|^2 \leq ( \| \boldx \| + \| \boldy\|)^2$. Taking the square root gives the triangle inequality.

Given two vectors $\boldx$ and $\boldy$, the value $\| \boldx - \boldy \|$ has the interpretation of being the "distance" between these points. To see why, consult [figure %s](#f-vec_minus) again.

(ss-lcs)=
### Linear Combinations and Span
One of the most elementary ways to work with vectors is to combine them using linear operations. Given vectors $\boldx_1,\ldots, \boldx_K$ in $\RR^N$, a **linear combination** of these vectors is a new vector of the form

$$
\boldy 
        = \sum_{k=1}^K \alpha_k \boldx_k 
        =  \alpha_1 \boldx_1 + \cdots +  \alpha_K \boldx_K
$$ (eq-lc)

for some collection of scalars $\alpha_1,\ldots, \alpha_K$ (i.e., with $\alpha_k \in \RR$ for all $k$). {numref}`f-lin_comb` shows four different linear combinations $\boldy = \alpha_1 \boldx_1 + \alpha_2 \boldx_2$ where $\boldx_1, \boldx_2$ are fixed vectors in $\RR^2$ and the scalars $\alpha_1$ and $\alpha_2$ are varied.

```{figure} figures/lin_comb.svg
:name: f-lin_comb

Linear combinations of $\boldx_1, \boldx_2$
```

Given any nonempty $X \subset \RR^N$, the set of all vectors that can be made by (finite) linear combinations of elements of $X$ is called the **span** of $X$, and denoted by $\Span X$. For example, the set of all linear combinations of $X := \{\boldx_1,\ldots, \boldx_K\}$ is

```{math}
:enumerated: false

\Span X:= \left\{ \text{ all vectors } \sum_{k=1}^K \alpha_k \boldx_k 
    \text{ such that }
    \boldalpha := (\alpha_1,\ldots, \alpha_K) \in \RR^K \right\}
```

As will be discussed below, the span of certain collections of vectors turns out to have an intimate connection with existence of solutions to linear equations.

```{prf:example}

By construction, the four vectors labeled $\boldy$ in [figure %s](#f-lin_comb) lie in the span of $X = \{\boldx_1, \boldx_2\}$. Looking at this picture might lead you to wonder whether *any* vector in $\RR^2$ could be created as a linear combination of $\boldx_1, \boldx_2$. The answer is affirmative. We'll prove this in {ref}`ss-bad`.
```


```{prf:example}
:label: eg-lineinplane

Let $X = \{ \boldone \} = \{ (1,1) \} \subset \RR^2$. The span of $X$ is all vectors of the form $\alpha \boldone =  (\alpha, \alpha)$ with $\alpha \in \RR$. This constitutes a line in the plane. Since we can take $\alpha = 0$, it follows that the origin $\boldzero$ is in $\Span X$. In fact $\Span X$ is the unique line in the plane that passes through both $\boldzero$ and the vector $\boldone = (1,1)$.
```


```{prf:example}

Let $\boldx_1 = (3, 4, 2)$ and let $\boldx_2 = (3, -4, 0.4)$. The span of $\{\boldx_1, \boldx_2\}$ is a plane in $\RR^3$ that passes through both of these vectors and the origin, as shown in [figure %s](#f-span_plane).
```


```{figure} figures/span_plane.svg
:name: f-span_plane

Span of $X = \{\boldx_1, \boldx_2\}$
```

````{prf:example}
:label: eg-canonbv0

Consider the vectors $\{\bolde_1, \ldots, \bolde_N\} \subset \RR^N$, where $\bolde_n$ has all zeros except for a $1$ as the $n$th element:

```{math}
:enumerated: false

\bolde_1 := 
        \left(
        \begin{array}{c}
            1 \\
            0 \\
            \vdots \\
            0
        \end{array}
        \right),
        \quad 
        \bolde_2 := 
        \left(
        \begin{array}{c}
            0 \\
            1 \\
            \vdots \\
            0
        \end{array}
        \right),
        \quad 
        \cdots,
        \quad
        \bolde_N := 
        \left(
        \begin{array}{c}
            0 \\
            0 \\
            \vdots \\
            1
        \end{array}
        \right)
```

The case of $\RR^2$ is illustrated in [figure %s](#f-vec_canon). The vectors $\bolde_1, \ldots, \bolde_N$ are called the **canonical basis vectors** of $\RR^N$. One reason is that $\{\bolde_1, \ldots, \bolde_N\}$ spans all of $\RR^N$. Here's a proof for $N=2$: Observe that for any $\boldy \in \RR^2$, we have

```{math}
:enumerated: false

\boldy 
        :=
        \left(
        \begin{array}{c}
            y_1 \\
            y_2
        \end{array}
        \right)
        =
        \left(
        \begin{array}{c}
            y_1 \\
            0
        \end{array}
        \right)
        +
        \left(
        \begin{array}{c}
            0 \\
            y_1
        \end{array}
        \right)
        =
        y_1
        \left(
        \begin{array}{c}
            1 \\
            0
        \end{array}
        \right)
        +
        y_2
        \left(
        \begin{array}{c}
            0 \\
            1
        \end{array}
        \right)
        = y_1 \bolde_1 + y_2 \bolde_2
```

Thus $\boldy \in \Span \{\bolde_1, \bolde_2\}$. Since $\boldy$ is an arbitrary vector in $\RR^2$, we have shown that $\{\bolde_1, \bolde_2\}$ spans $\RR^2$.
````


```{figure} figures/tikz_vec_canon_x.svg
:name: f-vec_canon

Canonical basis vectors in $\RR^2$
```

```{prf:example}
:label: eg-planein3d

Consider the set $P := \setntn{(x_1, x_2, 0) \in \RR^3}{ x_1, x_2 \in \RR}$. Graphically, $P$ corresponds to the flat plane in $\RR^3$, where the height coordinate is always zero. If we take $\bolde_1 = (1, 0, 0)$ and $\bolde_2 = (0, 1, 0)$, then given $\boldy = (y_1, y_2, 0) \in P$, we have $\boldy = y_1 \bolde_1 + y_2 \bolde_2$. In other words, any $\boldy \in P$ can be expressed as a linear combination of $\bolde_1$ and $\bolde_2$. Equivalently, $P \subset \Span\{\bolde_1, \bolde_2\}$.
```


The next fact follows directly from the definition of span.

```{prf:property}
:label: fa-xsys

If $X, Y$ are nonempty subsets of $\RR^N$ and $X \subset Y$, then $\Span X \subset \Span Y$.
```


### Linear Independence

Linear independence is an apparently simple concept with implications that stretch deep into many aspects of analysis. If you wish to understand when a matrix is invertible, or when a system of linear equations has a unique solution, or when a least squares estimate is uniquely defined, the most important foundational idea is linear independence of vectors.

Let's begin with the definition. Consider a set of vectors $X := \{\boldx_1,\ldots, \boldx_K\}$. We can surely realize the origin $\boldzero$ as a linear combination of these vectors, just by setting all of the scalars $\alpha_k$ in $\sum_{k=1}^K \alpha_k \boldx_k$ to zero. The set $X$ is called linearly independent when this is the only possibility. That is, $X \subset \RR^N$ is called **linearly independent** if

$$
\alpha_1 \boldx_1
    + \cdots +
    \alpha_K \boldx_K
     = \boldzero 
    \quad \implies \quad
    \alpha_1 = \cdots = \alpha_K = 0
$$ (eq-lif)

We call $X$ **linearly dependent** if it is not linearly independent.

````{prf:example}
:label: eg-lcind

In [figure %s](#f-lin_comb), the two vectors are $\boldx_1 = (1.2, 1.1)$ and $\boldx_2 = (-2.2, 1.4)$. Suppose that $\alpha_1$ and $\alpha_2$ are scalars with

```{math}
:enumerated: false

\alpha_1
        \left(
        \begin{array}{c}
            1.2 \\
            1.1
        \end{array}
        \right)
        + 
        \alpha_2
        \left(
        \begin{array}{c}
            -2.2 \\
            1.4
        \end{array}
        \right)
        =
        \boldzero
```

This translates to the linear, two-equation system shown in [figure %s](#f-alpha_eq), where the unknowns are $\alpha_1$ and $\alpha_2$. The only solution is $\alpha_1 = \alpha_2 = 0$. Hence $\{\boldx_1, \boldx_2\}$ is linearly independent.
````


```{figure} figures/alpha_eq.svg
:name: f-alpha_eq

The only solution is $\alpha_1 = \alpha_2 = 0$
```

````{prf:example}
:label: eg-canonbv1

The set of $N$ canonical basis vectors $\{\bolde_1, \ldots, \bolde_N\}$ is linearly independent in $\RR^N$. To see this, let $\alpha_1, \ldots, \alpha_N$ be coefficients such that $\sum_{n=1}^N \alpha_n \bolde_n = \boldzero$. Equivalently,

```{math}
:enumerated: false

\alpha_1
        \left(
        \begin{array}{c}
            1 \\
            0 \\
            \vdots \\
            0
        \end{array}
        \right)
        +
        \alpha_2
        \left(
        \begin{array}{c}
            0 \\
            1 \\
            \vdots \\
            0
        \end{array}
        \right)
        + \cdots +
        \alpha_N
        \left(
        \begin{array}{c}
            0 \\
            0 \\
            \vdots \\
            1
        \end{array}
        \right)
        =
        \left(
        \begin{array}{c}
            \alpha_1 \\
            \alpha_2 \\
            \vdots \\
            \alpha_N
        \end{array}
        \right)
        =
        \left(
        \begin{array}{c}
            0 \\
            0 \\
            \vdots \\
            0
        \end{array}
        \right)
```

In particular, $\alpha_n = 0$ for all $n$.
````


```{prf:example}
:label: eg-vec_noncanon

Consider the vectors $\{\boldx_1, \boldx_2\}$ as given in

![](figures/tikz_vec_noncanon.svg)

This pair fails to be linearly independent, since $\boldx_2 = -2 \boldx_1$, and hence $2 \boldx_1 + \boldx_2 = \boldzero$.
```


How can we interpret linear independence? One way to understand it is as an indicator of the algebraic diversity of a given collection of vectors. In particular, in a linearly independent set, the span is relatively large, in the sense that every vector contributes to the span. Here's a formal statement of this idea.

```{prf:theorem}
:label: t-fdlinin

Let $X := \{\boldx_1,\ldots, \boldx_K\} \subset \RR^N$. For $K > 1$, the following statements are equivalent:

(i) $X$ is linearly independent.

(ii) $X_0$ is a proper subset[^2] of $X$ $\,\implies\,$ $\Span X_0$ is a proper subset of $\Span X$.

(iii) No vector in $X$ can be written as a linear combination of the others.
```


{prf:ref}`ex-ott` asks you to check these equivalences. For now let's just step through them in the context of two examples. First consider the pair of canonical basis vectors $\{\bolde_1, \bolde_2\}$ in $\RR^2$, as depicted in [figure %s](#f-vec_canon). As we saw in [examples %s](#eg-canonbv0) and [%s](#eg-canonbv1), this pair is linearly independent, and its span is all of $\RR^2$. Both vectors contribute to the span, since removing either one reduces the span to just a line in $\RR^2$. (For example, the span of $\{\bolde_1\}$ is just the horizontal axis in $\RR^2$.) Neither one of this pair can be written as a linear combination of the other.

Next, consider instead the pair $\{\boldx_1, \boldx_2\}$ in [example %s](#eg-vec_noncanon). These vectors fail to be linearly independent, as shown in that example. It is also clear that dropping either one does not change the span—it remains the horizontal axis in any case. Finally, we saw in [example %s](#eg-vec_noncanon) that $\boldx_2 = -2 \boldx_1$, which means that each vector can be written as a linear combination of the other.

```{prf:property}
:label: fa-slili

If $X := \{\boldx_1,\ldots, \boldx_K\}$ is linearly independent, then

(i) every subset of $X$ is linearly independent,

(ii) $X$ does not contain $\boldzero$, and

(iii) $X \cup \{ \boldx \}$ is linearly independent for all $\boldx \in \RR^N$ such that $\boldx \notin \Span X$.
```


The proof is a solved exercise ([ex. %s](#ex-slili)).

#### Linear Independence and Uniqueness

As we'll see below, the problem of existence of solutions to systems of linear equations comes down to whether or not a given point is contained in the span of a collection of vectors (which typically correspond to the columns of a matrix). This depends in general on the size of the span, and the size of the span in turn depends on whether or not the vectors are linearly independent.

Given that linear independence is the key condition for existence of solutions, it's surprising at first to learn that linear independence is the key condition for uniqueness as well. As we'll see, the connection between linear independence and uniqueness stems from the following result.

```{prf:theorem}
:label: t-unor

Let $X := \{\boldx_1,\ldots,\boldx_K\}$ be any collection of vectors in $\RR^N$. The following statements are equivalent:

(i) $X$ is linearly independent.

(ii) For each $\boldy \in \RR^N$ there exists at most one set of scalars $\alpha_1,\ldots,\alpha_K$ such that

$$
\boldy =  
                    \alpha_1 \boldx_1
                    + \cdots +
                    \alpha_K \boldx_K
$$ (eq-yilcx)
```


```{prf:proof}
:nonumber:

((i)$\implies$(ii)) Let $X$ be linearly independent and pick any $\boldy \in \RR^N$. Suppose that there are two sets of scalars such that {eq}`eq-yilcx` holds. In particular, suppose that $\boldy = \sum_{k=1}^K \alpha_k \boldx_k = \sum_{k=1}^K \beta_k \boldx_k$. It follows from the second equality that $\sum_{k=1}^K (\alpha_k - \beta_k) \boldx_k = \boldzero$. By linear independence, we then have $\alpha_k = \beta_k$ for all $k$. In other words, the representation is unique.

((ii)$\implies$(i)) If (ii) holds, then there exists at most one set of scalars such that $\boldzero = \sum_{k=1}^K \alpha_k \boldx_k$. Because $\alpha_1 = \cdots = \alpha_k = 0$ has this property, we conclude that no nonzero scalars yield $\boldzero = \sum_{k=1}^K \alpha_k \boldx_k$. In other words, $X$ is linearly independent. ◻
```


(ss-linsub)=
### Linear Subspaces
We will find it rewarding to study more closely the structure of spans generated by a collection of vectors. One of the defining features of the span of a set $X$ is that it is "closed" under the linear operations of vector addition and scalar multiplication, in the sense that

(i) $\boldx, \boldy \in \Span X \; \implies \; \boldx + \boldy \in \Span X$, and

(ii) $\boldy \in \Span X$ and $\gamma \in \RR \; \implies \; \gamma \boldy \in \Span X$.

For example, (i) holds because the sum of two linear combinations of elements of $X$ is another linear combination of elements of $X$.

The notion of a set being closed under scalar multiplication and vector addition is important enough to have its own name: A nonempty subset $S$ of $\RR^N$ is called a **linear subspace** (or just **subspace**) of $\RR^N$ if

$$
\boldx, \, \boldy \in S
    \text{ and } \alpha, \, \beta \in \RR
    \; \implies \; 
    \alpha \boldx + \beta \boldy \in S
$$ (eq-lsubd)

```{prf:example}
:label: eg-silsub

It follows from the preceding discussion that if $X$ is any nonempty subset of $\RR^N$, then $\Span X$ is a linear subspace of $\RR^N$. For this reason, $\Span X$ is often called the **linear subspace spanned by $X$.**
```


````{prf:example}
:label: eg-ipwass

Given any $\bolda \in \RR^N$, the set $A := \setntn{\boldx \in \RR^N}{\inner{\bolda, \boldx}  = 0}$ is a linear subspace of $\RR^N$. To see this let $\boldx, \boldy \in A$ and let $\alpha, \beta \in \RR$. We claim that $\boldz := \alpha \boldx + \beta \boldy \in A$, or, equivalently, that $\inner{\bolda,\boldz} = 0$. This is true because

```{math}
:enumerated: false

\inner{\bolda, \boldz} =
        \inner{\bolda, \alpha \boldx + \beta \boldy} 
        = \alpha \inner{ \bolda, \boldx}
        + \beta \inner{\bolda, \boldy} = 0 + 0 = 0
```
````


```{prf:example}

The entire space $\RR^N$ is a linear subspace of $\RR^N$ because any linear combination of $N$-vectors is an $N$-vector.
```


To visualize subspaces in $\RR^3$, think of lines and planes that pass through the origin. Here are some elementary facts about linear subspaces:

```{prf:property}
:label: fa-eplsub

If $S$ is a linear subspace of $\RR^N$, then

(i) $\boldzero \in S$,

(ii) $X \subset S \implies \Span X \subset S$, and

(iii) $\Span S = S$.
```


There's also one deep result about linear subspaces we need to cover, which forms a cornerstone of many foundational results:

```{prf:theorem}
:label: t-exth

Let $S$ be a linear subspace of $\RR^N$. If $S$ is spanned by $K$ vectors, then any linearly independent subset of $S$ has at most $K$ vectors.
```


In other words, if there exists a set $X = \{\boldx_1, \ldots, \boldx_K\}$ with $S \subset \Span X$, then any subset of $S$ with more than $K$ vectors will be linearly dependent. The proof can be found in most texts on linear algebra (e.g., §3.5 of {cite:author}`jan` {cite:year}`jan`).

```{prf:example}

We saw in [example %s](#eg-canonbv0) that $\RR^2$ is spanned by the pair $\{\bolde_1, \bolde_2\}$, where $\bolde_i$ is the $i$th canonical basis vector in $\RR^2$. (See also [figure %s](#f-vec_canon).) It follows immediately from this fact and [theorem %s](#t-exth) that the three vectors in $\RR^2$ shown in [figure %s](#f-vecs) are linearly dependent.
```


```{prf:example}

Consider the plane $P := \setntn{(x_1, x_2, 0) \in \RR^3}{ x_1, x_2 \in \RR}$ from [example %s](#eg-planein3d). We saw in that example that $P$ can be spanned by two vectors. As a consequence of [theorem %s](#t-exth), we now know that any three vectors in this plane—such as the three shown in [figure %s](#f-flat_plane)—are linearly dependent.
```


```{figure} figures/flat_plane.svg
:name: f-flat_plane

Any three vectors in $P$ are linearly dependent
```

(ss-bad)=
### Bases and Dimension
Consider again the pair $\boldx_1, \boldx_2$ shown in [figure %s](#f-lin_comb), and the four different vectors labeled $\boldy$ that we created from $\boldx_1, \boldx_2$ by way of linear combinations. Eyeballing the figure gives the impression that any $\boldy \in \RR^2$ could be constructed as a linear combination of $\boldx_1, \boldx_2$ with suitable choice of the scalars $\alpha_1, \alpha_2$. Indeed this is true. The reason is that the pair $\{\boldx_1, \boldx_2\}$ is linearly independent (see [example %s](#eg-lcind)), and any pair of linearly independent vectors in $\RR^2$ spans $\RR^2$. Here's a statement of this result for the general case:

```{prf:theorem}
:label: t-bsili

Let $X := \{ \boldx_1, \ldots, \boldx_N \}$ be any $N$ vectors in $\RR^N$. The following statements are equivalent:

(i) $\Span X = \RR^N$.

(ii) $X$ is linearly independent.
```


```{prf:proof}
:nonumber:

((i)$\implies$(ii)) Suppose that $\Span X = \RR^N$ but $X$ is not linearly independent. Then, by [theorem %s](#t-fdlinin), there exists a proper subset $X_0$ of $X$ with $\Span X_0 = \Span X$. Since $X_0$ is a proper subset of $X$ it contains $K < N$ elements. We now have $K$ vectors spanning $\RR^N$. In particular, the span of $K$ vectors contains the $N>K$ linearly independent vectors $\bolde_1, \ldots, \bolde_N$. This contradicts [theorem %s](#t-exth).

((ii)$\implies$(i)) Suppose that $X$ is linearly independent and yet there exists an $\boldx \in \RR^N$ with $\boldx \notin \Span X$. By [fact %s](#fa-slili), it follows that the $N+1$ element set $X \cup \{\boldx\} \subset \RR^N$ is linearly independent. Since $\RR^N$ is spanned by the $N$ canonical basis vectors, this statement also contradicts [theorem %s](#t-exth). ◻
```


We now come to a key definition. Let $S$ be a linear subspace of $\RR^N$ and let $B \subset S$. The set $B$ is called a **basis** of $S$ if

(i) $B$ spans $S$ and

(ii) $B$ is linearly independent.

The plural of basis is **bases**. In view of [theorem %s](#t-unor), when $B$ is a basis of $S$, each point in $S$ has exactly one representation as a linear combination of elements of $B$.

It follows immediately from [theorem %s](#t-bsili) that any $N$ linearly independent vectors in $\RR^N$ form a basis of $\RR^N$.

```{prf:example}
:label: eg-canonbv

The set of canonical basis vectors $\{\bolde_1, \ldots, \bolde_N\} \subset \RR^N$ described in [example %s](#eg-canonbv0) is linearly independent and spans all of $\RR^N$. As a result, it provides a basis for $\RR^N$—as anticipated by the name.
```


```{prf:example}

The pair $\{\boldx_1, \boldx_2\}$ from [figure %s](#f-lin_comb) forms a basis of $\RR^2$.
```


```{prf:example}

The pair $\{\bolde_1, \bolde_2\}$ is a basis for the set $P$ defined in [example %s](#eg-planein3d).
```


Here are the two most fundamental results about bases:

```{prf:theorem}
:label: t-basis

If $S$ is a linear subspace of $\RR^N$ distinct from $\{\boldzero\}$, then

(i) $S$ has at least one basis and

(ii) every basis of $S$ has the same number of elements.
```


The proof of part (i) is not particularly hard. See, for example, section 3.2 of {cite}`jan`. Part (ii) follows from [theorem %s](#t-exth) and is left as an exercise ([ex. %s](#ex-ebse)).

If $S$ is a linear subspace of $\RR^N$, then the common number identified in [theorem %s](#t-basis) is called the **dimension** of $S$, and written as $\dim S$.

```{prf:example}

For $P := \setntn{(x_1, x_2, 0) \in \RR^3}{ x_1, x_2 \in \RR}$ we have $\dim P = 2$ because $\{\bolde_1, \bolde_2\} \subset \RR^3$ is a basis (see [example %s](#eg-planein3d)) and $\{\bolde_1, \bolde_2\}$ has two elements.
```


```{prf:example}

$\dim \RR^N = N$ because $\{\bolde_1, \ldots, \bolde_N\} \subset \RR^N$ is a basis.
```


```{prf:example}

A line $\setntn{\alpha \boldx \in \RR^N}{\alpha \in \RR}$ through the origin is one dimensional.
```


In $\RR^N$ the singleton subspace $\{\boldzero\}$ is said to have zero dimension.

If we take a set of $K$ vectors, then how large will its span be in terms of dimension? The next theorem answers this question.

```{prf:theorem}
:label: t-dimspan

Let $X := \{\boldx_1,\ldots,\boldx_K\} \subset \RR^N$. Then

(i) $\dim \Span X \leq K$ and

(ii) $\dim \Span X = K$ if and only if $X$ is linearly independent.
```


{prf:ref}`ex-dimspan` asks you to prove these results.

Let's finish this section with facts that can be deduced from the preceding results.

```{prf:property}
:label: fa-ondsrn

The following statements are true:

(i) Let $S$ and $S'$ be $K$-dimensional linear subspaces of $\RR^N$. If $S \subset S'$, then $S = S'$.

(ii) If $S$ is an $M$-dimensional linear subspace of $\RR^N$ and $M < N$, then $S \neq \RR^N$.
```


Part (i) of [fact %s](#fa-ondsrn) implies that the only $N$-dimensional linear subspace of $\RR^N$ is $\RR^N$ itself.

(ss-lf)=
### Linear Maps
The single most important class of functions in applied mathematics is the linear functions. In high school we are told that linear functions as those whose graph is a straight line. Here's a better definition: A function $T \colon \RR^K \to \RR^N$ is **linear** if

$$
T(\alpha \boldx + \beta \boldy) = \alpha T\boldx + \beta T\boldy
        \; \text{ for any } \boldx, \boldy \in \RR^K \text{ and }
        \alpha, \beta \in \RR
$$ (eq-dlinf)

(Following a common convention, we'll write linear functions with uppercase letters and omit the parenthesis around the argument where no confusion arises. This convention has come about because the action of linear maps is essentially isomorphic to multiplication of vectors by matrices. More on that topic soon.)

```{prf:example}

The function $T \colon \RR \to \RR$ defined by $Tx = 2x$ is linear because, for any $\alpha, \beta, x, y$ in $\RR$, we have $T(\alpha x + \beta y) = 2(\alpha x + \beta y) = \alpha 2 x + \beta 2 y = \alpha Tx + \beta Ty$.
```


````{prf:example}

Given $\bolda \in \RR^K$, the function $T \colon \RR^K \to \RR$ defined by $T \boldx = \inner{\bolda, \boldx}$ is linear. Indeed, by the rules for inner products, for any $\alpha, \beta$ in $\RR$ and $\boldx, \boldy$ in $\RR^K$ we have

```{math}
:enumerated: false

T(\alpha \boldx + \beta \boldy)
      = \inner{\bolda, \alpha \boldx + \beta \boldy}
      = \alpha \inner{\bolda, \boldx}  + \beta \inner{\bolda, \boldy}
      = \alpha T \boldx + \beta T \boldy
```
````


```{prf:example}

The function $f \colon \RR \to \RR$ defined by $f(x) = x^2$ fails to be linear because, if $\alpha = \beta = x = y = 1$, then $f(\alpha x + \beta y)  = 4$, while $\alpha f(x) + \beta f(y) = 2$.
```


```{prf:example}
:label: eg-affun

The function $f \colon \RR \to \RR$ defined by $f(x) = 1 + 2x$ is not linear because if $\alpha = \beta = x = y = 1$, then $f(\alpha x + \beta y) = f(2) = 5$, while $\alpha f(x) + \beta f(y) = 3 + 3 = 6$. This kind of function is called an **affine** function. We see that identifying linear functions with functions whose graph is a straight line is not correct.
```


The definition in {eq}`eq-dlinf` tells us directly that if $T$ is linear then the exchange of order in $T [ \sum_{k=1}^K \alpha_k \boldx_k ] = \sum_{k=1}^K \alpha_k T \boldx_k$ will be valid whenever $K=2$. A simple inductive argument extends this to arbitrary $K$. As an application of this fact, consider the following: As discussed in [example %s](#eg-canonbv0), any $\boldx \in \RR^K$ can be expressed in terms of the basis vectors as $\sum_{k=1}^K \alpha_k \bolde_k$, for some suitable choice of scalars. Hence, for a linear function $T$, its range, denoted $\range T$, is the set of all points of the form

```{math}
:enumerated: false

T\boldx
    = T \left[ \sum_{k=1}^K \alpha_k \bolde_k \right]
    = \sum_{k=1}^K \alpha_k T \bolde_k
```

as we vary $\alpha_1, \ldots, \alpha_K$ over all scalar combinations. (See {ref}`ss-functions` for the definition of range.) In other words, the range of a linear map is the span of the image of the canonical basis functions. This will prove to be important later on. The next fact summarizes.

````{prf:property}
:label: fa-res

If $T \colon \RR^K \to \RR^N$ is a linear map, then

```{math}
:enumerated: false

\range T = \Span V,
        \quad \text{where } 
        V := \{T\bolde_1, \ldots, T\bolde_K\}
```
````


Soon we'll turn to the topic of determining when linear functions are bijections, an issue that is intimately related to invertibility of matrices. To this end it's useful to note that, for linear functions, the property of being one-to-one can be determined by examining the set of points it maps to the origin. To express this idea, for any $T \colon \RR^K \to \RR^N$, we define the **null space** or **kernel** of $T$ as

```{math}
:enumerated: false

\kernel T := \setntn{\boldx \in \RR^K}{T\boldx = \boldzero}
```

```{prf:property}
:label: fa-kernel

If $T \colon \RR^K \to \RR^N$ is linear, then

(i) $\kernel T$ is a linear subspace of $\RR^K$ and

(ii) $\kernel T = \{\boldzero\}$ if and only if $T$ is one-to-one.
```


The proofs are straightforward. For example, if $T\boldx = T\boldy$ for some $\boldx, \boldy \in \RR^K$, then $T(\boldx - \boldy) = \boldzero$, and hence $\boldx - \boldy \in \kernel T$. So if $\kernel T = \{\boldzero\}$, then $T\boldx = T\boldy$ implies $\boldx = \boldy$, which means that $T$ is one-to-one.

(ss-lfb)=
### Linear Independence and Bijections
Many statistical problems are "inverse\" problems, in the sense that we observe outcomes and wish to determine what generated them. For example, we might want to know what consumer preferences led to observed market behavior, or what kinds of expectations led to a given shift in exchange rates.

```{figure} figures/xfig_inverse_prob.svg
:name: f-inverse_prob

Inverse problem
```

Consider the generic inverse problem in [figure %s](#f-inverse_prob), where $F$ and $y$ are given, and we seek to obtain the unknown object $x$. Two immediate questions are: Does this problem have a solution? and Is it unique? To provide general answers to these questions, we need to know whether $F$ is one-to-one, onto, etc. (see {ref}`ss-functions` for definitions and further discussion). The best case is when $F$ is a bijection, for then we know that a unique solution $x$ exists for every possible $y$.

In general, functions can be onto, one-to-one, bijections, or none of the above. However, for linear functions from $\RR^N$ to $\RR^N$, the first three properties are all equivalent! The next theorem gives details.

```{prf:theorem}
:label: t-linfunc

If $T$ is a linear function from $\RR^N$ to $\RR^N$, then all of the following are equivalent:

(i) $T$ is a bijection.

(ii) $T$ is onto.

(iii) $T$ is one-to-one.

(iv) $\kernel T = \{ \boldzero \}$.

(v) $V := \{T\bolde_1, \ldots, T\bolde_N\}$ is linearly independent.

(vi) $V := \{T\bolde_1, \ldots, T\bolde_N\}$ forms a basis of $\RR^N$.
```


If any one of these conditions is true, then $T$ is called **nonsingular**. (Equivalently, a nonsingular function is a linear bijection.) Otherwise $T$ is called **singular**. The proof of [theorem %s](#t-linfunc) is a solved exercise ([ex. %s](#ex-linfunc)). {numref}`f-linbijec` provides intuition for the case of $N=1$. In the top panel all conditions in [theorem %s](#t-linfunc) are satisfied. In the lower panel none are. In particular, we can see that the condition for $T$ to be one-to-one and onto is exactly the same (i.e., $\alpha \neq 0$).

```{figure} figures/linbijec.svg
:name: f-linbijec

The case of $N=1$, nonsingular and singular
```

If $T$ is nonsingular, then, being a bijection, it must have an inverse function $T^{-1}$ that is also a bijection ([fact %s](#fa-bij)). It turns out that this inverse function inherits the linearity of $T$ (see [ex. %s](#ex-nsins)). In summary,

```{prf:property}
:label: fa-nsins

If $T \colon \RR^N \to \RR^N$ is nonsingular, then so is $T^{-1}$.
```


{prf:ref}`t-linfunc` only applies to linear maps between spaces of the *same* dimension. When linear functions map across distinct dimensions the situation changes:

```{prf:theorem}
:label: t-lfoc

For a linear map $T$ from $\RR^K \to \RR^N$, the following statements are true:

(i) If $K < N$, then $T$ is not onto.

(ii) If $K > N$, then $T$ is not one-to-one.
```


The most important implication is that if $N \neq K$, then we can forget about bijections. The proof of [theorem %s](#t-lfoc) is a solved exercise ([ex. %s](#ex-lfoc)).

(s-orthog)=
## Orthogonality
One of the core concepts in this book is orthogonality, not just of vectors but also of more complex objects such as random variables. Let's begin with the vector definition and some key implications.

(ss-odbp)=
### Definition and Basic Properties
Let $\boldx$ and $\boldz$ be vectors in $\RR^N$. If $\inner{\boldx,  \boldz} = 0$, then we write **$\boldx \perp \boldz$** and call $\boldx$ and $\boldz$ **orthogonal**. In $\RR^2$, $\boldx$ and $\boldz$ are orthogonal when they are perpendicular to one another, as in [figure %s](#f-xpz). For $\boldx \in \RR^N$ and $S \subset \RR^N$, we say that **$\boldx$ is orthogonal to $S$** if $\boldx \perp \boldz$ for all $\boldz \in S$ ([figure %s](#f-xpS)), and we write **$\boldx \perp S$**. A set of vectors $\{\boldz_1, \ldots, \boldz_K\} \subset \RR^N$ is called an **orthogonal set** if its elements are mutually orthogonal, that is, if $\boldz_j \perp \boldz_k$ whenever $j$ and $k$ are distinct.

````{prf:property}
:label: fa-pythag

If $\{\boldz_1,\ldots, \boldz_K\}$ is an orthogonal set, then

```{math}
:enumerated: false

\| \boldz_1 + \cdots + \boldz_K \|^2 
        = \| \boldz_1 \|^2 + \cdots + \| \boldz_K \|^2
```
````


```{figure} figures/tikz_xpz.svg
:name: f-xpz

$\boldx \perp \boldz$
```

```{figure} figures/tikz_xpS.svg
:name: f-xpS

$\boldx \perp S$
```

Orthogonal sets and linear independence are closely related. For example,

```{prf:property}
:label: fa-orthii

If $O \subset \RR^N$ is an orthogonal set and $\boldzero \notin O$, then $O$ is linearly independent.
```


While not every linearly independent set is orthogonal, an important partial converse to [fact %s](#fa-orthii) is given in {ref}`ss-odo`.

An orthogonal set $O \subset \RR^N$ is called an **orthonormal set** if $\| \boldu \| = 1$ for all $\boldu \in O$. An orthonormal set that lies in and spans a linear subspace $S$ of $\RR^N$ is called an **orthonormal basis** of $S$. It is, necessarily, a basis of $S$. (Why?) The standard example of an orthonormal basis for all of $\RR^N$ is the canonical basis $\{\bolde_1, \ldots, \bolde_N\}$.

By definition, if $O = \{\boldu_1, \ldots, \boldu_K\}$ is any basis of $S$, then, for any $\boldx \in S$, we can find unique scalars $\alpha_1, \ldots, \alpha_K$ such that $\boldx = \sum_{k=1}^K \alpha_k \boldu_k$. While the values of these scalars are not always transparent, for an orthonormal basis they are easy to compute:

```{prf:property}
:label: fa-bdons

If $\{\boldu_1, \ldots, \boldu_K\}$ is an orthonormal set and $\boldx \in \Span\{\boldu_1, \ldots, \boldu_K\}$, then

$$
\boldx = \sum_{k=1}^K \inner{\boldx, \boldu_k} \boldu_k
$$ (eq-exbyob)
```


The proof is an exercise. Given $S \subset \RR^N$, the **orthogonal complement** of $S$ is defined as

```{math}
:enumerated: false

S^{\perp} := \setntn{\boldx \in \RR^N}{\boldx \perp S}
```

{numref}`f-orth_comp` gives an example in $\RR^2$.

```{figure} figures/tikz_orth_comp.svg
:name: f-orth_comp

Orthogonal complement of $S$
```

```{prf:property}
:label: fa-ocls

For any nonempty $S \subset \RR^N$, the set $S^{\perp}$ is a linear subspace of $\RR^N$.
```


Indeed, if $\boldx, \boldy \in S^{\perp}$ and $\alpha, \beta \in \RR$, then $\alpha \boldx + \beta \boldy \in S^{\perp}$ because, for any $\boldz \in S$,

```{math}
:enumerated: false

\inner{\alpha \boldx + \beta \boldy, \boldz}
     = \alpha \inner{ \boldx, \boldz} + \beta \inner{\boldy, \boldz}
     = \alpha \times 0  + \beta \times 0 
     = 0
```

```{prf:property}
:label: fa-scscez

For $S \subset \RR^N$, we have $S \cap S^{\perp} = \{\boldzero\}$.
```


### The Orthogonal Projection Theorem

A central problem in linear regression and many other applications is approximation of some $\boldy$ of $\RR^N$ by an element of a given subspace $S$ of $\RR^N$. Stated more precisely, the problem is, given $\boldy$ and $S$, to find the closest element $\hat \boldy$ of $S$ to $\boldy$. Closeness is in terms of Euclidean norm, so $\hat \boldy$ is the minimizer of $\| \boldy - \boldz \|$ over all $\boldz \in S$:

$$
\hboldy = \argmin_{\boldz \in S} \|\boldy - \boldz\|
$$ (eq-mproj)

The next theorem tells us that a solution $\hboldy$ to this minimization problem always exists, as well as providing a means to identify it.

```{prf:theorem}
:label: t-opt

Let $\boldy \in \RR^N$ and let $S$ be any nonempty linear subspace of $\RR^N$. The following statements are true:

(i) The optimization problem {eq}`eq-mproj` has exactly one solution.

(ii) $\hboldy \in \RR^N$ solves {eq}`eq-mproj` if and only if $\hat \boldy \in S$ and $\boldy - \hat \boldy \perp S$.

The unique solution $\hboldy$ is called the **orthogonal projection of $\boldy$ onto $S$**.
```


The intuition is easy to grasp from a graphical presentation. {numref}`f-orth_proj2D0` illustrates. Looking at the figure, we can see that the closest point to $\boldy$ in $S$ is the one point $\hboldy \in S$ such that $\boldy - \hboldy$ is orthogonal to $S$.

For a full proof see, for example, theorem 5.16 of {cite}`erhan2013real`. Let's just cover sufficiency of the conditions in part (ii): Let $\boldy \in \RR^N$ and let $S$ be a linear subspace of $\RR^N$. Let $\hat \boldy$ be a vector in $S$ satisfying $\boldy - \hat \boldy \perp S$. Let $\boldz$ be any other point in $S$. We have

```{math}
:enumerated: false

\| \boldy - \boldz \|^2
    = \| (\boldy - \hat \boldy) + (\hat \boldy - \boldz) \|^2
    = \| \boldy - \hat \boldy \|^2  + \| \hat \boldy - \boldz  \|^2
```

The second equality follows from $\boldy - \hat \boldy \perp S$ (why?) and the Pythagorian law. Since $\boldz$ was an arbitrary point in $S$, we have $\| \boldy - \boldz \| \geq \| \boldy - \hat \boldy \|$ for all $\boldz \in S$. Hence {eq}`eq-mproj` holds.

```{figure} figures/tikz_orth_proj2D0.svg
:name: f-orth_proj2D0

Orthogonal projection
```

````{prf:example}
:label: eg-pvones

Let $\boldy \in \RR^N$ and let $\boldone \in \RR^N$ be the vector of ones. Let $S$ be the set of constant vectors in $\RR^N$, meaning that all elements are equal. Evidently $S$ is the span of $\{\boldone\}$. The orthogonal projection of $\boldy$ onto $S$ is $\hboldy := \bar y \boldone$, where $\bar y := \frac{1}{N} \sum_{n=1}^N y_n$. To see this, note that $\hboldy \in S$ clearly holds. Hence we only need to check that $\boldy - \hboldy$ is orthogonal to $S$, for which it suffices to show that $\inner{\boldy - \hboldy, \boldone} = 0$ (see [ex. %s](#ex-pecb)). This is true because

```{math}
:enumerated: false

\inner{\boldy - \hboldy, \boldone }
        = \inner{\boldy, \boldone} - \inner{\hboldy, \boldone}
        = \sum_{n=1}^N y_n - \bar y \inner{\boldone, \boldone}
        = 0
```
````


### Projection as a Mapping

In view of [theorem %s](#t-opt), for each fixed linear subspace $S$ in $\RR^N$, the operation

```{math}
:enumerated: false

\boldy \; \mapsto \text{ the orthogonal projection of } \boldy \text{ onto } S
```

is a well-defined function from $\RR^N$ to $\RR^N$. The function is typically denoted by $\boldP$. For each $\boldy \in \RR^N$, the symbol $\boldP \boldy$ represents the image of $\boldy$ under $\boldP$, which is the orthogonal projection $\hboldy$. $\boldP$ is called the **orthogonal projection onto $S$**, and we write

```{math}
:enumerated: false

\boldP = \proj S
```

{numref}`f-orth_proj2Dp` illustrates the action of $\boldP$ on two different vectors.

```{figure} figures/tikz_orth_proj2Dp.svg
:name: f-orth_proj2Dp

Orthogonal projection under $\boldP$
```

Using this notation, we can restate the orthogonal projection theorem, as well as adding some properties of $\boldP$:

```{prf:theorem}
:label: t-opt2

Let $S$ be any linear subspace of $\RR^N$, and let $\boldP = \proj S$. The following statements are true:

(i) $\boldP$ is a linear function.

Moreover, for any $\boldy \in \RR^N$, we have

(ii) $\boldP \boldy \in S$,

(iii) $\boldy - \boldP \boldy \perp S$,

(iv) $\| \boldy \|^2 = \| \boldP \boldy \|^2 + \| \boldy - \boldP \boldy \|^2$,

(v) $\| \boldP \boldy \| \leq \| \boldy \|$,

(vi) $\boldP \boldy = \boldy$ if and only if $\boldy \in S$, and

(vii) $\boldP \boldy = \boldzero$ if and only if $\boldy \in S^{\perp}$.
```


These results are not difficult to prove, given [theorem %s](#t-opt). Linearity of $\boldP$ is left as an exercise ([ex. %s](#ex-prpil)). Parts (ii)–(iii) follow directly from [theorem %s](#t-opt). To see (iv), observe that $\boldy$ can be decomposed as $\boldy  =  \boldP \boldy + \boldy - \boldP \boldy$. Now apply the Pythagorean law. Part (v) follows from part (iv). Part (vi) follows from the definition of $\boldP \boldy$ as the closest point to $\boldy$ in $S$. Part (vii) is an exercise.

```{prf:property}
:label: fa-projon

If $\{\boldu_1, \ldots, \boldu_K\}$ is an orthonormal basis for $S$, then, for each $\boldy \in \RR^N$,

$$
\boldP \boldy = \sum_{k=1}^K \inner{\boldy, \boldu_k} \boldu_k
$$ (eq-projon)
```


[Fact %s](#fa-projon) is a fundamental result. It's true because the right-hand side of {eq}`eq-projon` clearly lies in $S$ (being a linear combination of basis functions) and, for any $\boldu_j$ in the basis set

```{math}
:enumerated: false

\inner{\boldy - \boldP \boldy, \boldu_j}
    = \inner{\boldy, \boldu_j }
    -  \sum_{k=1}^K \inner{\boldy, \boldu_k} \inner{\boldu_k, \boldu_j}
    = \inner{\boldy, \boldu_j} - \inner{\boldy, \boldu_j }
    = 0
```

This is enough to confirm that $\boldy - \boldP \boldy \perp S$ (see [ex. %s](#ex-pecb)).

```{prf:example}
:label: eg-pvones2

Recall [example %s](#eg-pvones), where we showed that the projection of $\boldy \in \RR^N$ onto $\Span\{\boldone\}$ is $\bar y \boldone$, where $\bar y$ is the "sample mean\" $\bar y := \frac{1}{N} \sum_{n=1}^N y_n$. We can see this from {eq}`eq-projon` too. To apply {eq}`eq-projon`, we just need to find an orthonormal basis for $\Span\{\boldone\}$. The obvious candidate is $\{N^{-1/2} \boldone\}$. Applying {eq}`eq-projon` now gives $\boldP \boldy = \inner{N^{-1/2} \boldone, \boldy} N^{-1/2} \boldone$. As before, this leads us to $\bar y \boldone$.
```


There's one more essential property of $\boldP$ that we need to make note of: Suppose that we have two linear subspaces $S_1$ and $S_2$ of $\RR^N$, where $S_1 \subset S_2$. What then is the difference between (1) first projecting a point onto the bigger subspace $S_2$, and then projecting the result onto the smaller subspace $S_1$, and (2) projecting directly to the smaller subspace $S_1$? The answer is none—we get the same result:

````{prf:property}
:label: fa-subsub

Let $S_i$ be a linear subspace of $\RR^N$ for $i=1,2$ and let $\boldP_i = \proj S_i$. If $S_1 \subset S_2$, then

```{math}
:enumerated: false

\boldP_1 \boldP_2 \boldy = \boldP_2 \boldP_1 \boldy = \boldP_1 \boldy
        \quad \text{for all } \boldy \in \RR^N
```
````


(ss-odo)=
### The Residual Projection
Consider the setting of the orthogonal projection theorem. Our interest is in projecting $\boldy$ onto $S$, where $S$ is a linear subspace of $\RR^N$. The closest point to $\boldy$ in $S$ is $\hboldy := \boldP \boldy$ where $\boldP = \proj S$. Unless $\boldy$ was already in $S$, some error $\boldy - \boldP \boldy$ remains. Tracking and managing this residual will be important to us, so let's introduce an operator $\boldM$ that takes $\boldy \in \RR^N$ and returns the residual. We can define it as

$$
\boldM := \boldI - \boldP
$$ (eq-ann0)

where $\boldI$ is the identity mapping on $\RR^N$. For any $\boldy$ we have $\boldM \boldy = \boldI \boldy - \boldP \boldy = \boldy - \boldP \boldy$ as required. In regression analysis $\boldM$ shows up as a matrix called the "annihilator." This is a pretty cool name, but it's also not a great description of its function. In what follows we will refer to $\boldM$ as the **residual projection**.

```{prf:example}
:label: eg-pvonesr

Recall [example %s](#eg-pvones), where we found that the projection of $\boldy \in \RR^N$ onto $\Span\{\boldone\}$ is $\bar y \boldone$. The residual projection is $\boldM_c \boldy := \boldy - \bar y \boldone$. In econometric applications, we'll view this as a vector of errors obtained when the elements of a vector are predicted by its sample mean. The subscript reminds us that $\boldM_c$ centers vectors around their mean.
```


```{prf:property}
:label: fa-opt3

Let $S$ be a linear subspace of $\RR^N$, let $\boldP = \proj S$, and let $\boldM$ be the residual projection as defined in {eq}`eq-ann0`. The following statements are true:

(i) $\boldM = \proj S^{\perp}$.

(ii) $\boldy = \boldP \boldy + \boldM \boldy$ for any $\boldy \in \RR^N$.

(iii) $\boldP \boldy \perp \boldM \boldy$ for any $\boldy \in \RR^N$.

(iv) $\boldM \boldy = \boldzero$ if and only if $\boldy \in S$.

(v) $\boldP \circ \boldM = \boldM \circ \boldP = \boldzero$.
```


Part (v) means that $\boldP \boldM \boldy = \boldM \boldP \boldy = \boldzero$ for all $\boldy \in \RR^N$. {numref}`f-orth_proj2D` illustrates the action of $\boldM$. The results in [fact %s](#fa-opt3) can be seen in the figure.

```{figure} figures/tikz_orth_proj2D.svg
:name: f-orth_proj2D

The residual projection
```

If $S_1$ and $S_2$ are two subspaces of $\RR^N$ with $S_1 \subset S_2$, then $S_2^{\perp} \subset S_1^{\perp}$. This means that the result in [fact %s](#fa-subsub) is reversed for $\boldM$.

````{prf:property}
:label: fa-subsub2

Let $S_1$ and $S_2$ be two subspaces of $\RR^N$ and let $\boldy \in \RR^N$. Let $\boldM_1$ and $\boldM_2$ be the projections onto $S_1^{\perp}$ and $S_2^{\perp}$ respectively. If $S_1 \subset S_2$, then

```{math}
:enumerated: false

\boldM_1 \boldM_2 \boldy = \boldM_2 \boldM_1 \boldy = \boldM_2 \boldy
```
````


As an application of the ideas above, let's now discuss a procedure called Gram–Schmidt orthogonalization, which provides a fundamental link between the two major concepts discussed in this chapter: linear independence and orthogonality. It can be considered as a partial converse to [fact %s](#fa-orthii).

````{prf:theorem}
:label: t-gso

For each linearly independent set $\{\boldb_1, \ldots, \boldb_K\} \subset \RR^N$, there exists an orthonormal set $\{\boldu_1, \ldots, \boldu_K\}$ with

```{math}
:enumerated: false

\Span \{\boldb_1, \ldots, \boldb_k\}
        = 
        \Span \{\boldu_1, \ldots, \boldu_k\}
        \quad \text{for } \;
        k = 1, \ldots, K
```
````


The proof of [theorem %s](#t-gso) provides an important algorithm for generating the orthonormal set $\{\boldu_1, \ldots, \boldu_K\}$. The first step is to construct orthogonal sets $\{\boldv_1, \ldots, \boldv_k\}$ with span identical to $\{\boldb_1, \ldots, \boldb_k\}$ for each $k$. The construction of $\{\boldv_1, \ldots, \boldv_K\}$ uses the so called **Gram–Schmidt orthogonalization** procedure. First, for each $k = 1, \ldots, K$, let

(i) $B_k := \Span\{\boldb_1, \ldots, \boldb_k\}$,

(ii) $\boldP_k := \proj B_k$ and $\boldM_k := \proj B_k^{\perp}$,

(iii) $\boldv_k := \boldM_{k-1} \boldb_k$ where $\boldM_0$ is the identity mapping, and

(iv) $V_k := \Span\{\boldv_1, \ldots, \boldv_k\}$.

In step (iii) we map each successive element $\boldb_k$ into a subspace orthogonal to the subspace generated by $\boldb_1, \ldots, \boldb_{k-1}$. In the exercises you are asked to show that $\{\boldv_1, \ldots, \boldv_K\}$ is an orthogonal set and $V_k = B_k$ for all $k$ ([ex. %s](#ex-gsso)). To complete the argument, we introduce the vectors $\boldu_k$ by $\boldu_k := \boldv_k / \| \boldv_k \|$ and confirm that this set of vectors $\{\boldu_1, \ldots, \boldu_k\}$ is orthonormal with span equal to $V_k$. These results are also (solved) exercises.

## Further Reading

Good texts on vector spaces include {cite}`marcus1988introduction` and {cite}`jan`.

## Exercises

````{exercise}
:label: ex-ruinpr

Show that inner products of linear combinations satisfy the following rule:

```{math}
:enumerated: false

\inner{
        \sum_{k=1}^K \alpha_k \boldx_k,
        \sum_{j=1}^J \beta_j \boldy_j 
        }
          =  \sum_{k=1}^K \sum_{j=1}^J \alpha_k \beta_j 
          \inner{\boldx_k, \boldy_j }
```
````


```{exercise}
:label: ex-ch02_vector_spaces-auto-1

Show that the vectors $(1,1)$ and $(-1,2)$ are linearly independent.
```


```{exercise}
:label: ex-ch02_vector_spaces-auto-2

Use [fact %s](#fa-profno) to show that if $\boldy \in \RR^N$ is such that $\inner{\boldy, \boldx} = 0$ for every $\boldx \in \RR^N$, then $\boldy = \boldzero$.
```


````{exercise}
:label: ex-csap

Fix nonzero $\boldx \in \RR^N$. Consider the optimization problem

```{math}
:enumerated: false

\max_{\boldy} \inner{\boldx,\boldy}
        \quad \text{ subject to }  \quad
        \boldy \in \RR^N \text{ and } \| \boldy \| = 1
```

Show that the maximizer is $\hboldx := (1/ \| \boldx \|) \boldx$.[^3]
````


```{exercise}
:label: ex-2isls3

Is $\RR^2$ a linear subspace of $\RR^3$? Why or why not?
```


```{exercise}
:label: ex-zik

Show that if $T \colon \RR^K \to \RR^N$ is a linear function, then $\boldzero \in \ker T$.
```


```{exercise}
:label: ex-gmsi

Let $\{\boldx_1, \boldx_2\}$ be a linearly independent set in $\RR^2$ and let $\gamma$ be a nonzero scalar. Is it true that $\{\gamma \boldx_1, \gamma \boldx_2\}$ is also linearly independent?
```


````{exercise}
:label: ex-irvis

Is it true that

```{math}
:enumerated: false

\boldz := 
    \begin{pmatrix}
        -3.9  \\
        12.4 \\
        -6.8
    \end{pmatrix}
    \in
    \Span\{\boldx_1, \boldx_2, \boldx_3\}
    \quad \text{when} \quad
    \boldx_1 = 
    \begin{pmatrix}
        -4  \\
        0 \\
        0
    \end{pmatrix},
    \;
    \boldx_2 = 
    \begin{pmatrix}
        0  \\
        2 \\
        0
    \end{pmatrix},
    \;
    \boldx_3 = 
    \begin{pmatrix}
        0  \\
        0 \\
        -1
    \end{pmatrix}
```

Why or why not?
````


```{exercise}
:label: ex-inls

Show that if $S$ and $S'$ are two linear subspaces of $\RR^N$, then $S \cap S'$ is also a linear subspace of $\RR^N$.
```


```{exercise}
:label: ex-ch02_vector_spaces-auto-3

Prove [fact %s](#fa-eplsub).
```


```{exercise}
:label: ex-qalss

Let $Q := \setntn{(x_1, x_2, x_3) \in \RR^3}{x_2 = x_1 + x_3}$. Is $Q$ a linear subspace of $\RR^3$?
```


```{exercise}
:label: ex-ch02_vector_spaces-auto-4

Let $Q := \setntn{(x_1, x_2, x_3) \in \RR^3}{x_2 = 1}$. Is $Q$ a linear subspace of $\RR^3$?
```


```{exercise}
:label: ex-ch02_vector_spaces-auto-5

Show that if $T \colon \RR^N \to \RR^N$ is a linear function and $\lambda$ is any scalar, then $E := \setntn{\boldx \in \RR^N}{T \boldx = \lambda \boldx}$ is a linear subspace of $\RR^N$.
```


```{exercise}
:label: ex-pecb

Show that if $B \subset S$ with $\Span B = S$, then $\boldx \perp S$ if and only if $\boldx \perp \boldb$ for all $\boldb \in B$.
```


```{exercise}
:label: ex-ott

Prove the equivalences in [theorem %s](#t-fdlinin).
```


```{exercise}
:label: ex-slili

Prove [fact %s](#fa-slili).
```


```{exercise}
:label: ex-ebse

Show that if $S$ is a linear subspace of $\RR^N$ then every basis of $S$ has the same number of elements.
```


```{exercise}
:label: ex-ondsrn

Prove [fact %s](#fa-ondsrn).
```


```{exercise}
:label: ex-dimspan

Prove [theorem %s](#t-dimspan).
```


```{exercise}
:label: ex-nsins

Show that if $T \colon \RR^N \to \RR^N$ is nonsingular (i.e., a linear bijection), then $T^{-1}$ is also linear.
```


```{exercise}
:label: ex-linfunc

Prove [theorem %s](#t-linfunc).
```


```{exercise}
:label: ex-lfoc

Prove [theorem %s](#t-lfoc).
```


```{exercise}
:label: ex-ch02_vector_spaces-auto-6

Find two unit vectors (i.e., vectors with norm equal to one) that are orthogonal to $(1, -2)$.
```


```{exercise}
:label: ex-ch02_vector_spaces-auto-7

Prove the Pythagorean law ([fact %s](#fa-pythag)). See [ex. %s](#ex-ruinpr) if you need a hint.
```


```{exercise}
:label: ex-orthii

Prove [fact %s](#fa-orthii).
```


```{exercise}
:label: ex-opt3

Prove [fact %s](#fa-opt3) using [theorems %s](#t-opt) and [%s](#t-opt2).
```


```{exercise}
:label: ex-scscez

Prove [fact %s](#fa-scscez): If $S \subset \RR^N$, then $S \cap S^{\perp} = \{\boldzero\}$.
```


```{exercise}
:label: ex-ch02_vector_spaces-auto-8

Prove [fact %s](#fa-subsub).
```


```{exercise}
:label: ex-prpil

Let $\boldP$ be the orthogonal projection described in [theorem %s](#t-opt2). Confirm that $\boldP$ is a linear function from $\RR^N$ to $\RR^N$.
```


```{exercise}
:label: ex-fcpop

Let $S := \setntn{(x_1, x_2, x_3) \in \RR^3}{x_3 = 0}$ and let $\boldy := \boldone := (1,1,1)$. Using the orthogonal projection theorem, find the closest point in $S$ to $\boldy$.
```


```{exercise}
:label: ex-prpil2

Let $S$ be any linear subspace of $\RR^N$ and let $\boldP = \proj S$ (see [theorem %s](#t-opt2)). Is $\boldP$ one-to-one as a function on $\RR^N$?
```


```{exercise}
:label: ex-rtine

Prove the reverse triangle inequality. That is, given two vectors $\boldx$ and $\boldy$, show that $| \| \boldx \| - \| \boldy \| | \leq \| \boldx - \boldy \|$.[^4]
```


```{exercise}
:label: ex-pmiz

Show that $\boldP \boldy = \boldzero$ implies $\boldy \in S^{\perp}$.
```


In the next three exercises, the notation is as given in [theorem %s](#t-gso) and the discussion immediately afterwards.

```{exercise}
:label: ex-gsso

Show that $V_k = B_k$ for all $k$.
```


```{exercise}
:label: ex-gsst

Show that $\{\boldv_1, \ldots, \boldv_K\}$ is an orthogonal set.
```


```{exercise}
:label: ex-gssth

Show that $\{\boldu_1, \ldots, \boldu_k\}$ is an orthonormal set with span equal to $V_k$ for all $k$.
```


### Solutions to Selected Exercises

````{solution} ex-csap

Fix nonzero $\boldx \in \RR^N$. Let $\hboldx := \boldx / \| \boldx \|$. Comparing this point with any other $\boldy \in \RR^N$ satisfying $\| \boldy \| = 1$, the Cauchy–Schwarz inequality yields

```{math}
:enumerated: false

\inner{\boldy, \boldx}
        \leq |\inner{\boldy, \boldx}|
        \leq \| \boldy \| \| \boldx \|
        = \| \boldx \|
        = \frac{\inner{\boldx, \boldx}}{ \| \boldx \| }
        = \inner{\hboldx, \boldx}
```

Hence $\hboldx$ is the maximizer, as claimed.
````


```{solution} ex-2isls3

This is a bit of a trick question. To solve it, you need to look carefully at the definitions (as always). A linear subspace of $\RR^3$ is a subset of $\RR^3$ with certain properties. $\RR^3$ is a collection of 3-tuples $(x_1, x_2, x_3)$ where each $x_i$ is a real number. Elements of $\RR^2$ are 2-tuples (pairs), and hence not elements of $\RR^3$. Therefore $\RR^2$ is not a subset of $\RR^3$, and, in particular, not a linear subspace of $\RR^3$.
```


````{solution} ex-zik

Let $T$ be as in the question. We need to show that $T \boldzero = \boldzero$. Here's one proof. We know from the definition of scalar multiplication that $0 \boldx = \boldzero$ for any vector $\boldx$. Let $\boldx$ and $\boldy$ be any vectors in $\RR^K$ and apply the definition of linearity to obtain

```{math}
:enumerated: false

T\boldzero 
        = T(0 \boldx + 0 \boldy) = 0 T\boldx + 0 T \boldy 
        = \boldzero + \boldzero 
        = \boldzero
        
```
````


```{solution} ex-gmsi

The answer is yes. Suppose, to the contrary, that $\{\gamma \boldx_1, \gamma \boldx_2\}$ is linearly dependent. Then one element can be written as a linear combination of the others. In our setting with only two vectors, this translates to $\gamma \boldx_1 = \alpha \gamma \boldx_2$ for some $\alpha$. Since $\gamma \neq 0$, we can multiply each side by $1/\gamma$ to get $\boldx_1 = \alpha \boldx_2$. This contradicts linear independence of $\{\boldx_1, \boldx_2\}$.
```


````{solution} ex-irvis

There is an easy way to do this: We know that any linearly independent set of 3 vectors in $\RR^3$ will span $\RR^3$. Since $\boldz \in \RR^3$, this will include $\boldz$. So all we need to do is show that the 3 vectors are linearly independent. To this end, take any scalars $\alpha_1, \alpha_2, \alpha_3$ with

```{math}
:enumerated: false

\alpha_1
    \begin{pmatrix}
        -4  \\
        0 \\
        0
    \end{pmatrix}
    +
    \alpha_2
    \begin{pmatrix}
        0  \\
        2 \\
        0
    \end{pmatrix}
    +
    \alpha_3
    \begin{pmatrix}
        0  \\
        0 \\
        -1
    \end{pmatrix}
    =
    \boldzero
    :=
    \begin{pmatrix}
        0  \\
        0 \\
        0
    \end{pmatrix}
```

Written as three equations, this says that $-4 \alpha_1 = 0$, $2 \alpha_2 = 0$, and $-1 \alpha_3 = 0$. Hence $\alpha_1 = \alpha_2 = \alpha_3 = 0$, and therefore the set is linearly independent.
````


```{solution} ex-inls

Let $S$ and $S'$ be two linear subspaces of $\RR^N$. Fix $\boldx, \boldy \in S \cap S'$ and $\alpha, \beta \in \RR$. We claim that $\boldz := \alpha \boldx + \beta \boldy \in S \cap S'$. To see this, note that since $\boldx, \boldy \in S$ and $S$ is a linear subspace, we have $\boldz \in S$; and since $\boldx, \boldy \in S'$ and $S'$ is a linear subspace, we have $\boldz \in S'$. It follows that $\boldz \in S \cap S'$, as was to be shown.
```


```{solution} ex-qalss

If $\bolda := (1, -1, 1)$, then $Q$ is all $\boldx$ with $\inner{\bolda,\boldx} = 0$. This set is a linear subspace of $\RR^3$, as shown in [example %s](#eg-ipwass).
```


```{solution} ex-ott

We are asked to verify the equivalences in [theorem %s](#t-fdlinin) for the set $X := \{\boldx_1,\ldots, \boldx_K\}$. We will prove the cycle (i) $\implies$ (ii) $\implies$ (iii) $\implies$ (i).

((i)$\implies$(ii)) We aim to show that if (i) holds and $X_0$ is a proper subset of $X$, then $\Span X_0$ is a proper subset of $\Span X$. To simplify notation let's take $X_0 := \{\boldx_2,\ldots, \boldx_K\}$. Suppose, to the contrary, that $\Span X_0 = \Span X$. Since $\boldx_1 \in \Span X$, we must then have $\boldx_1 \in \Span X_0$, from which we deduce the existence of scalars $\alpha_2, \ldots, \alpha_K$ such that $\boldzero = -\boldx_1 + \alpha_2 \boldx_2 + \cdots + \alpha_K \boldx_K$. Since $-1 \neq 0$, this contradicts part (i).

((ii)$\implies$(iii)) The claim is that when (ii) holds, no vector in $X$ can be written as a linear combination of the others. Suppose, to the contrary, that $\boldx_1 = \alpha_2 \boldx_2 + \cdots + \alpha_K \boldx_K$, say. Let $\boldy \in \Span X$, so that $\boldy = \beta_1 \boldx_1 + \cdots + \beta_K \boldx_K$. If we use the preceding equality to substitute out $\boldx_1$, we get $\boldy$ as a linear combination of $\{\boldx_2, \ldots, \boldx_K\}$ alone. In other words, any element of $\Span X$ is in the span of the proper subset $\{\boldx_2, \ldots, \boldx_K\}$. Contradiction.

((iii)$\implies$(i)) The final claim is that $\alpha_1 = \cdots = \alpha_K = 0$ whenever $\alpha_1 \boldx_1 + \cdots  + \alpha_K \boldx_K = \boldzero$. Suppose, to the contrary, that there exist scalars with $\alpha_1 \boldx_1 + \cdots  + \alpha_K \boldx_K = \boldzero$ and yet $\alpha_k \neq 0$ for at least one $k$. It follows immediately that $\boldx_k = (1/\alpha_k) \sum_{j \neq k} \alpha_j \boldx_j$. This contradicts (iii).
```


```{solution} ex-slili

The aim is to prove [fact %s](#fa-slili). Regarding the part (i), let's take $X$ as linearly independent and show that the subset $X_0 := \{\boldx_1,\ldots, \boldx_{K-1}\}$ is linearly independent. (The argument for more general subsets is similar.) Suppose, to the contrary, that $X_0$ is linearly dependent. Then, by the definition, we can take $\alpha_1, \ldots, \alpha_{K-1}$ not all zero with $\sum_{k=1}^{K-1} \alpha_k \boldx_k = \boldzero$. Letting $\alpha_K =0$, we can write this as $\sum_{k=1}^K \alpha_k \boldx_k = \boldzero$. Since not all coefficients are zero, we have contradicted linear independence of $X$.

Regarding (ii), let $X := \{\boldx_1,\ldots, \boldx_K\}$ be linearly independent and suppose that $\boldx_j = \boldzero$. Then by setting $\alpha_k = 0$ for $k\neq j$ and $\alpha_j = 1$, we can form scalars not all equal to zero with $\sum_{k=1}^K \alpha_k \boldx_k = \boldzero$.

Regarding (iii), let $X := \{\boldx_1,\ldots, \boldx_K\} \subset \RR^N$ be linearly independent and let $\boldx_{K+1}$ be any point in $\RR^N$ such that $\boldx_{K+1} \notin \Span X$. The claim is that $X \cup \{\boldx_{K+1}\}$ is linearly independent. Suppose, to the contrary, that there exist $\alpha_1, \ldots, \alpha_K, \alpha_{K+1}$ not all zero such that $\sum_{k=1}^{K+1} \alpha_k \boldx_k = \boldzero$. There are two possibilities for $\alpha_{K+1}$, both of which lead to a contradiction: First, if $\alpha_{K+1} = 0$, then, since $\alpha_1, \ldots, \alpha_K, \alpha_{K+1}$ are not all zero, at least one of $\alpha_1, \ldots, \alpha_K$ are nonzero, and moreover $\sum_{k=1}^K \alpha_k \boldx_k = \sum_{k=1}^{K+1} \alpha_k \boldx_k = \boldzero$. This contradicts our assumption of independence on $X$. Second, if $\alpha_{K+1} \neq 0$, then from $\sum_{k=1}^{K+1} \alpha_k \boldx_k = \boldzero$ we can express $\boldx_{K+1}$ as a linear combination of elements of $X$. This contradicts the hypothesis that $\boldx_{K+1} \notin \Span X$.
```


```{solution} ex-ebse

Let $B_1$ and $B_2$ be two bases of $S$, with $K_1$ and $K_2$ elements respectively. By definition, $B_2$ is a linearly independent subset of $S$. Moreover, $S$ is spanned by the set $B_1$, which has $K_1$ elements. Applying [theorem %s](#t-exth), we see that $B_2$ has at most $K_1$ elements. That is, $K_2 \leq K_1$. Reversing the roles of $B_1$ and $B_2$ gives $K_1 \leq K_2$.
```


```{solution} ex-ondsrn

The aim is to prove [fact %s](#fa-ondsrn). Suppose that $S$ and $S'$ are $K$-dimensional linear subspaces of $\RR^N$ with $S \subset S'$. We claim that $S = S'$. To see this, observe that by the definition of dimension, $S$ is equal to $\Span B$ where $B$ is a set of $K$ linearly independent basis vectors $\{\boldb_1, \ldots, \boldb_K\}$. If $S \neq S'$, then there exists a vector $\boldx \in S'$ such that $\boldx \notin \Span B$. In view of [theorem %s](#t-fdlinin), the set $\{\boldx, \boldb_1, \ldots, \boldb_K\}$ is linearly independent. Moreover, since $\boldx \in S'$ and since $B \subset S \subset S'$, we now have $K+1$ linearly independent vectors inside $S'$. At the same time, being $K$-dimensional, we know that $S'$ is spanned by $K$ vectors. This contradicts [theorem %s](#t-exth).

Regarding part (ii), suppose that $S$ is an $M$-dimensional linear subspace of $\RR^N$ where $M < N$ and yet $S = \RR^N$. Then we have a space $S$ spanned by $M < N$ vectors that also contains the $N$ linearly independent canonical basis vectors. We are led to another contradiction of [theorem %s](#t-exth). Hence $S = \RR^N$ cannot hold.
```


```{solution} ex-dimspan

Regarding part (i), let $B$ be a basis of $\Span X$. By definition, $B$ is a linearly independent subset of $\Span X$. Since $\Span X$ is spanned by $K$ vectors, [theorem %s](#t-exth) implies that $B$ has no more than $K$ elements. Hence, $\dim \Span X \leq K$.

Regarding part (ii), suppose first that $X$ is linearly independent. Then $X$ is a basis for $\Span X$. Since $X$ has $K$ elements, we conclude that $\dim \Span X = K$. Conversely, if $\dim \Span X = K$, then $X$ must be linearly independent. If $X$ were not linearly independent, then there would exist a proper subset $X_0$ of $X$ such that $\Span X_0 = \Span X$. By part (i) of this theorem, we then have $\dim \Span X_0 \leq \# X_0 \leq K - 1$. Therefore $\dim \Span X \leq K - 1$, a contradiction.
```


````{solution} ex-nsins

Let $T \colon \RR^N \to \RR^N$ be nonsingular, and let $T^{-1}$ be its inverse. To see that $T^{-1}$ is linear, we need to show that for any pair $\boldx, \boldy$ in $\RR^N$ (which is the domain of $T^{-1}$) and any scalars $\alpha$ and $\beta$, the following equality holds:

$$
T^{-1}(\alpha \boldx + \beta \boldy) = \alpha T^{-1}\boldx + \beta T^{-1} \boldy
$$ (eq-il)
In the proof we will exploit the fact that $T$ is by assumption a linear bijection.

Pick any vectors $\boldx, \boldy \in \RR^N$ and any two scalars $\alpha, \beta$. Since $T$ is a bijection, we know that $\boldx$ and $\boldy$ have unique preimages under $T$. In particular, there exist unique vectors $\boldu$ and $\boldv$ such that $T\boldu = \boldx$ and $T\boldv = \boldy$. Using these definitions, linearity of $T$ and the fact that $T^{-1}$ is the inverse of $T$, we have

```{math}
:enumerated: false

T^{-1}(\alpha \boldx + \beta \boldy)
         = T^{-1}(\alpha T\boldu + \beta T \boldv) 
         = T^{-1}(T(\alpha \boldu + \beta  \boldv)) 
         = \alpha \boldu + \beta  \boldv 
         = \alpha T^{-1} \boldx + \beta T^{-1} \boldy
```

This chain of equalities confirms {eq}`eq-il`.
````


````{solution} ex-linfunc

A collection of equivalent statements such as this is usually proved via a cycle of implications, with the form (i) $\implies$ (ii) $\implies \cdots \implies$ (vi) $\implies$ (i). However, in this case the logic is clearer if we directly show that all statements are equivalent to linear independence of $V$.

First observe equivalence of the onto property and linear independence of $V$ via

```{math}
:enumerated: false

T \text{ onto } 
        \; \iff \;
        \range T  = \RR^N
        \; \iff \;
        \Span V  = \RR^N
```

by [fact %s](#fa-res). The last statement is equivalent to linear independence of $V$ by [theorem %s](#t-bsili).

Next let's show that $\kernel T = \{ \boldzero \}$ implies linear independence of $V$. To this end, suppose that $\kernel T = \{ \boldzero \}$ and let $\alpha_1, \ldots, \alpha_N$ be such that $\sum_{n=1}^N \alpha_n T\bolde_n = \boldzero$. By linearity of $T$, we then have $T(\sum_{n=1}^N \alpha_n \bolde_n) = \boldzero$. Since $\kernel T = \{ \boldzero \}$, this means that $\sum_{n=1}^N \alpha_n \bolde_n = \boldzero$, which in view of independence of $\{\bolde_1, \ldots, \bolde_N\}$, implies $\alpha_1 = \cdots = \alpha_N = 0$. This establishes that $V$ is linearly independent.

Now let's check that linear independence of $V$ implies $\kernel T = \{ \boldzero \}$. To this end, let $\boldx$ be a vector in $\RR^N$ such that $T\boldx = \boldzero$. We can represent $\boldx$ in the form $\sum_{n=1}^N \alpha_n \bolde_n$ for suitable scalars $\{\alpha_n\}$. From linearity and $T\boldx = \boldzero$, we get $\sum_{n=1}^N \alpha_n T\bolde_n = \boldzero$. By linear independence of $V$, this implies that each $\alpha_n = 0$, whence $\boldx = \boldzero$. Thus $\kernel T = \{ \boldzero \}$ as claimed.

From [fact %s](#fa-kernel) we have $\kernel T = \{ \boldzero \}$ iff $T$ is one-to-one, so we can now state the following equivalences

$$
T \text{ onto } 
        \iff
        V \text{ linearly independent }
        \iff
        T \text{ one-to-one }
$$ (eq-lmequiv)
If $T$ is a bijection, then $T$ is onto and hence $V$ is linearly independent by {eq}`eq-lmequiv`. Conversely, if $V$ is linearly independent then $T$ is both onto and one-to-one by {eq}`eq-lmequiv`. Hence $T$ is a bijection.

Finally, equivalence of linear independence of $V$ and the statement that $V$ forms a basis of $\RR^N$ is immediate from the definition of bases and [theorem %s](#t-bsili).
````


```{solution} ex-lfoc

Regarding part (i), let $K < N$ and let $T \colon \RR^K \to \RR^N$ be linear. $T$ cannot be onto because, if $T$ were onto, then we would have $\range T = \RR^N$, in which case the vectors in $V = \{T\bolde_1, \ldots, T\bolde_K\}$ in [fact %s](#fa-res) would span $\RR^N$, despite having only $K < N$ elements. This is impossible. (Why?)

Regarding part (ii), let $T \colon \RR^K \to \RR^N$ be linear and let $K > N$. Seeking a contradiction, suppose in addition that $T$ is one-to-one. Let $\{\alpha_k\}_{k=1}^K$ be such that $\sum_{k=1}^K \alpha_k T \bolde_k = \boldzero$. By linearity, $T(\sum_{k=1}^K \alpha_k \bolde_k) = \boldzero$, and since $T$ is one-to-one and $T \boldzero = \boldzero$, this in turn implies $\sum_{k=1}^K \alpha_k \bolde_k = \boldzero$. Since the canonical basis vectors are linearly independent, it must be that $\alpha_1 = \cdots = \alpha_K = 0$. From this we conclude that $\{T \bolde_1, \ldots, T\bolde_K\}$ is linearly independent. Thus $\RR^N$ contains $K$ linearly independent vectors, despite the fact that $N < K$. This is impossible by [theorem %s](#t-exth).
```


```{solution} ex-orthii

Let $O = \{\boldx_1, \ldots, \boldx_K\} \subset \RR^N$ be an orthogonal set that does not contain $\boldzero$. Let $\alpha_1, \ldots, \alpha_K$ be such that $\sum_{k=1}^K \alpha_k \boldx_k = \boldzero$. We claim that $\alpha_j = 0$ for any $j$. To see that this is so, fix $j$ and take the inner product of both sides of $\sum_{k=1}^K \alpha_k \boldx_k = \boldzero$ with respect to $\boldx_j$ to obtain $\alpha_j \| \boldx_j \|^2 = 0$. Since $\boldx_j \neq \boldzero$, we conclude that $\alpha_j = 0$. The proof is done.
```


```{solution} ex-scscez

Let $S \subset \RR^N$. We aim to show that $S \cap S^{\perp} = \{\boldzero\}$. Fix $\bolda \in S \cap S^{\perp}$. Since $\bolda \in S^{\perp}$, we know that $\inner{\bolda, \bolds} = 0$ for any $\bolds \in S$. Since $\bolda \in S$, we have in particular, $\inner{\bolda, \bolda} = \| \bolda \|^2 = 0$. As we saw in [fact %s](#fa-profno), the only such vector is $\boldzero$.
```


````{solution} ex-prpil

Fix $\alpha, \beta \in \RR$ and $\boldx, \boldy \in \RR^N$. The claim is that

```{math}
:enumerated: false

\boldP (\alpha \boldx + \beta \boldy) 
        = \alpha \boldP \boldx + \beta \boldP \boldy
```

To verify this equality, we need to show that the right-hand side is the orthogonal projection of $\alpha \boldx + \beta \boldy$ onto $S$. Going back to [theorem %s](#t-opt), we need to show that (i) $\alpha \boldP \boldx + \beta \boldP \boldy \in S$ and (ii) for any $\boldz \in S$, we have

```{math}
:enumerated: false

\inner{\alpha \boldx + \beta \boldy - 
        (\alpha \boldP \boldx + \beta \boldP \boldy), \boldz} = 0
```

Here (i) is immediate, because $\boldP \boldx$ and $\boldP \boldy$ are in $S$ by definition, and, moreover $S$ is a linear subspace. To see that (ii) holds, just note that

```{math}
:enumerated: false

\inner{\alpha \boldx + \beta \boldy - 
        (\alpha \boldP \boldx + \beta \boldP \boldy), \boldz }
        = \alpha \inner{\boldx - \boldP \boldx, \boldz} + 
        \beta \inner{\boldy - \boldP \boldy, \boldz}
```

By definition, the projections of $\boldx$ and $\boldy$ are orthogonal to $S$, so we have $\inner{\boldx - \boldP \boldx, \boldz} = \inner{\boldy - \boldP \boldy, \boldz} = 0$. We are done.
````


````{solution} ex-fcpop

Let $\boldx = (x_1, x_2, x_3)$ be the closest point in $S$ to $\boldy$. Note that $\bolde_1 \in S$ and $\bolde_2 \in S$. By the orthogonal projection theorem, we have (i) $\boldx \in S$ and (ii) $\boldy - \boldx \perp S$. From (i) we have $x_3 = 0$. From (ii) we have

```{math}
:enumerated: false

\la \boldy - \boldx, \bolde_1 \ra = 0
        \quad \text{and} \quad
        \la \boldy - \boldx, \bolde_2 \ra = 0
```

These equations can be expressed more simply as $1 - x_1 = 0$ and $1 - x_2 = 0$. We conclude that $\boldx = (1, 1, 0)$.
````


```{solution} ex-prpil2

If $S= \RR^N$, then $\boldP$ is the identity mapping (why?), which is one-to-one. If $S \neq \RR^N$, then take any $\boldx \notin S$. By definition, $\boldy := \boldP \boldx$ is in $S$, and hence $\boldy$ and $\boldx$ are distinct. But $\boldP$ maps elements of $S$ to themselves, so $\boldP\boldy = \boldy = \boldP \boldx$. Hence $\boldP$ is not one-to-one.
```


````{solution} ex-rtine

From the triangle inequality we have

```{math}
:enumerated: false

\| \boldx \| 
        = \| \boldx - \boldy + \boldy \|
        \leq \| \boldx - \boldy \| + \| \boldy \|
```

It follows that $\| \boldx \| - \| \boldy \| \leq \| \boldx - \boldy \|$. A similar argument reversing the roles of $\boldx$ and $\boldy$ gives $\| \boldy \| - \| \boldx \| \leq \| \boldx - \boldy \|$. Combining the last two inequalities gives

```{math}
:enumerated: false

- \| \boldx - \boldy \|
        \leq 
        \| \boldx \| - \| \boldy \| 
        \leq 
        \| \boldx - \boldy \|
```

This is equivalent to $| \| \boldx \| - \| \boldy \| | \leq \| \boldx - \boldy \|$.
````


```{solution} ex-pmiz

If $\boldP \boldy = \boldzero$, then $\boldy = \boldP \boldy + \boldM \boldy = \boldM \boldy$. Hence $\boldM$ does not shift $\boldy$. If an orthogonal projection onto a subspace doesn't shift a point, that's because the point is already in that subspace (see [theorem %s](#t-opt2)). In this case the subspace is $S^{\perp}$, and we conclude that $\boldy \in S^{\perp}$.
```


```{solution} ex-gsso

To see that $V_k = B_k$ for all $k$, fix $k$ and consider the claim that $V_k \subset B_k$. By definition, $\boldv_k = \boldb_k - \boldP_{k-1} \boldb_k$, and the two terms on the right-hand side lie in $B_k$. Hence $\boldv_k \in B_k$. Since spans increase as we add more elements, it follows that $\boldv_j \in B_k$ for $j \leq k$. In other words, $\{\boldv_1, \ldots, \boldv_k\} \subset B_k$. Since $B_k$ is a linear subspace, we have $V_k \subset B_k$.

An induction argument shows that $B_k \subset V_k$ also holds. Clearly it holds for $k=1$. Suppose that it also holds for $k-1$. From the definition of $\boldv_k$, we have $\boldb_k = \boldP_{k-1} \boldb_k + \boldv_k$. The first term on the right-hand side lies in $B_{k-1}$, which, by our induction hypothesis, satisfies $B_{k-1} \subset V_{k-1} \subset V_k$. The second term on the right-hand side is $\boldv_k$, which obviously lies in $V_k$. Hence both terms are in $V_k$, and therefore $\boldb_k \in V_k$. Using arguments analogous to the end of the last paragraph leads us to $B_k \subset V_k$.
```


```{solution} ex-gsst

To show that $\{\boldv_1, \ldots, \boldv_K\}$ is an orthogonal set, it suffices to check that $\boldv_k \perp \boldv_j$ whenever $j < k$. To see this, fix any pair $j < k$. By construction, $\boldv_k \in B_{k-1}^{\perp}$. But, as shown in the solution to [ex. %s](#ex-gsso), $\boldv_j \in B_{k-1}$ must hold. Hence $\boldv_k \perp \boldv_j$.
```


```{solution} ex-gssth

Since $\{\boldv_1, \ldots, \boldv_k\}$ is orthogonal, the family $\{\boldu_1, \ldots, \boldu_k\}$ will be orthonormal provided that the norm of each element is 1. This is true by construction, since $\boldu_k := \boldv_k / \| \boldv_k \|$. The only concern is that $\| \boldv_k \| = 0$ might hold for some $k$. But $\boldv_k = \boldzero$ is impossible because, if it was to hold, then by (vii) of [theorem %s](#t-opt2) we would have $\boldb_k \in B_{k-1}$, contradicting linear independence of $\{\boldb_1, \ldots, \boldb_K\}$.

The proof that $\Span\{\boldu_1, \ldots, \boldu_k\} = V_k$ is straightforward and left for you.
```


[^1]: In some instances, the notion of scalar multiplication includes multiplication of vectors by complex numbers. In what follows we will work almost entirely with real scalars. Hence scalar multiplication means real scalar multiplication unless otherwise stated.

[^2]: $A$ is a proper subset of $B$ if $A \subset B$ and $A \neq B$.

[^3]: Hint: There's no need to go taking derivatives and setting them equal to zero. An easier proof exists. If you're stuck, consider the Cauchy–Schwarz inequality.

[^4]: Hint: Use the triangle inequality.
