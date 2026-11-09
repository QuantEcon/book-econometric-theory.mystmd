---
title: "Linear Algebra and Matrices"
label: c-lam
---

The previous chapter was relatively abstract but it sets us up nicely to study the kinds of practical problems from linear algebra that we encounter in statistics and econometrics. This chapter treats many of the core problems in linear algebra, frequently relating them back to our analysis of vectors and linear maps in [chapter 2](#c-la).

## Matrices and Linear Equations

Matrices provide a convenient way to organize data and algebraic operations. Let's start with definitions.

(ss-mbd)=
### Basic Definitions
A **$N \times K$ matrix** is a rectangular array $\boldA$ of real numbers with $N$ rows and $K$ columns, written in the following way:

```{math}
:enumerated: false

\boldA = 
    \left(
    \begin{array}{cccc}
        a_{11} & a_{12} & \cdots & a_{1K} \\
        a_{21} & a_{22} & \cdots & a_{2K} \\
        \vdots & \vdots &  & \vdots \\
        a_{N1} & a_{N2} & \cdots & a_{NK} 
    \end{array}
    \right)
```

The symbol $a_{nk}$ stands for the element in the $n$th row of the $k$th column. Often these elements represent coefficients in a system of linear equations, such as

$$
\begin{array}{c}
        a_{11} x_1 + a_{12} x_2 + \cdots + a_{1K} x_K = b_1 \\
        \vdots  \\
        a_{N1} x_1 + a_{N2} x_2 + \cdots + a_{NK} x_K = b_N 
    \end{array}
$$ (eq-sleq)

For obvious reasons, $\boldA$ is also called a vector if either $N = 1$ or $K = 1$. In the former case, $\boldA$ is called a **row vector**, while in the latter case it is called a **column vector**. When convenient, we will use the notation $\row_n \boldA$ to refer to the $n$th row of $\boldA$, and $\col_k \boldA$ to refer to its $k$th column.

We extend the notation $\boldzero$ and $\boldone$ from vectors to matrices, in the sense that these symbols will also represent matrices with all elements equal to zero and one respectively. Dimensions will be stated explicitly or clear from context.

If $\boldA$ is $N \times K$ and $N = K$, then $\boldA$ is called **square**. For an $N \times N$ matrix $\boldA$, the $N$ elements of the form $a_{nn}$ for $n=1,\ldots,N$ are called the **principal diagonal**:

```{math}
:enumerated: false

\left(
    \begin{array}{cccc}
        \red{a_{11}} & a_{12} & \cdots & a_{1N} \\
        a_{21} & \red{a_{22}} & \cdots & a_{2N} \\
        \vdots & \vdots &  & \vdots \\
        a_{N1} & a_{N2} & \cdots & \red{a_{NN}} \\
    \end{array}
    \right)
```

The unique $N \times N$ matrix with ones on the principle diagonal and zeros elsewhere is called the **identity matrix**, and denoted by $\boldI$:

```{math}
:enumerated: false

\boldI := 
    \left(
    \begin{array}{cccc}
        1 & 0 & \cdots & 0 \\
        0 & 1 & \cdots & 0 \\
        \vdots & \vdots &  & \vdots \\
        0 & 0 & \cdots & 1 \\
    \end{array}
    \right)
```

Note that $\col_n \boldI = \bolde_n$, the $n$th canonical basis vector in $\RR^N$.

Just as was the case for vectors, a number of algebraic operations are defined for matrices. The first two, scalar multiplication and addition, are immediate generalizations of the vector case: For $\gamma \in \RR$, we let

```{math}
:enumerated: false

\gamma 
    \left(
    \begin{array}{cccc}
        a_{11} & a_{12} & \cdots & a_{1K} \\
        a_{21} & a_{22} & \cdots & a_{2K} \\
        \vdots & \vdots &  & \vdots \\
        a_{N1} & a_{N2} & \cdots & a_{NK} \\
    \end{array}
    \right)
    :=
    \left(
    \begin{array}{cccc}
        \gamma a_{11} & \gamma a_{12} & \cdots & \gamma a_{1K} \\
        \gamma a_{21} & \gamma a_{22} & \cdots & \gamma a_{2K} \\
        \vdots & \vdots &  & \vdots \\
        \gamma a_{N1} & \gamma a_{N2} & \cdots & \gamma a_{NK} \\
    \end{array}
    \right)
```

while

```{math}
:enumerated: false

\left(
    \begin{array}{ccc}
        a_{11} & \cdots & a_{1K} \\
        a_{21} & \cdots & a_{2K} \\
        \vdots & \vdots & \vdots \\
        a_{N1} & \cdots & a_{NK} \\
    \end{array}
    \right)
    +
    \left(
    \begin{array}{ccc}
        b_{11} & \cdots & b_{1K} \\
        b_{21} & \cdots & b_{2K} \\
        \vdots & \vdots & \vdots \\
        b_{N1} & \cdots & b_{NK} \\
    \end{array}
    \right)
    :=
    \left(
    \begin{array}{ccc}
        a_{11} + b_{11} &  \cdots & a_{1K} + b_{1K} \\
        a_{21} + b_{21} &  \cdots & a_{2K} + b_{2K} \\
        \vdots & \vdots & \vdots \\
        a_{N1} + b_{N1} &  \cdots & a_{NK} + b_{NK} \\
    \end{array}
    \right)
```

Addition is only defined for matrices that have identical shape.

The **matrix product** $\boldC = \boldA \boldB$ of matrices $\boldA$ and $\boldB$ is formed by taking as its $i,j$th element the inner product of the $i$th row of $\boldA$ and the $j$th column of $\boldB$. That is,

```{math}
:enumerated: false

c_{ij} = \inner{\row_i \boldA,\,  \col_j  \boldB} = \sum_{k=1}^K a_{ik} b_{kj}
```

Here's the picture for $i=j=1$:

```{math}
:enumerated: false

\left(
    \begin{array}{ccc}
        \red{a_{11}} & \red{\cdots} & \red{a_{1K}} \\
        a_{21} & \cdots & a_{2K} \\
        \vdots & \vdots & \vdots \\
        a_{N1} & \cdots & a_{NK} \\
    \end{array}
    \right)
    \left(
    \begin{array}{ccc}
        \red{b_{11}} & \cdots & b_{1J} \\
        \red{b_{21}} & \cdots & b_{2J} \\
        \red{\vdots} & \vdots & \vdots \\
        \red{b_{K1}} & \cdots & b_{KJ} \\
    \end{array}
    \right)
    =
    \left(
    \begin{array}{ccc}
        \red{c_{11}} & \cdots & c_{1J} \\
        c_{21} & \cdots & c_{2J} \\
        \vdots & \vdots & \vdots \\
        c_{N1} & \cdots & c_{NJ} \\
    \end{array}
    \right)
```

Since inner products are only defined for vectors of equal length, we need the length of the rows of $\boldA$ to be equal to the length of the columns of $\boldB$. In other words, if $\boldA$ is $N \times K$ and $\boldB$ is $J \times M$, then we require $K = J$. The resulting matrix $\boldA \boldB$ is $N \times M$. Here's the rule to remember:

```{math}
:enumerated: false

\boldA \text{ is } N \times K 
    \text{ and } 
    \boldB \text{ is } K \times M
    \; \implies \;
    \boldA \boldB
    \text{ is }  N \times M
```

Matrix multiplication is not commutative: $\boldA \boldB$ and $\boldB \boldA$ are not in general equal. In most other ways it behaves like ordinary multiplication:

```{prf:property}
:label: fa-bmari

For conformable matrices $\boldA$, $\boldB$, $\boldC$ and scalar $\alpha$, we have

(i) $\boldA (\boldB \boldC) = (\boldA \boldB) \boldC$,

(ii) $\boldA (\boldB + \boldC) = \boldA \boldB + \boldA \boldC$,

(iii) $(\boldA + \boldB) \boldC = \boldA \boldC + \boldB \boldC$,

(iv) $\boldA \alpha \boldB = \alpha \boldA \boldB$, and

(v) $\boldA \boldI = \boldA$ and $\boldI \boldA = \boldA$, where $\boldI$ is the identity matrix.
```


Here and below we use the word **conformable** to indicate dimensions are such that the operation in question makes sense. For example, two matrices are conformable for addition if they have the same number of rows and columns.

The **$k$th power** of a square matrix $\boldA$ is defined as

```{math}
:enumerated: false

\boldA^k := \underbrace{\boldA \cdots \boldA}_{k \text{ terms}}
```

If $\boldB$ is such that $\boldB^2 = \boldA$, then $\boldB$ is called the **square root** of $\boldA$ and written $\sqrt{\boldA}$.

Before going further, let's state the fundamental connection between matrix multiplication and the more elementary notion of linear combinations of vectors, as introduced in {ref}`ss-lcs`: Given $N \times K$ matrix $\boldA$ and $K \times 1$ column vector $\boldx$, the product $\boldA \boldx$ is an $N \times 1$ column vector formed as a linear combination of the columns of $\boldA$, with scalars $x_1, \ldots, x_K$. In symbols,

```{math}
:enumerated: false

\begin{aligned}
\boldA \boldx
    & = 
    \left(
    \begin{array}{cccc}
        a_{11} & a_{12} & \cdots & a_{1K} \\
        a_{21} & a_{22} & \cdots & a_{2K} \\
        \vdots & \vdots &  & \vdots \\
        a_{N1} & a_{N2} & \cdots & a_{NK} 
    \end{array}
    \right)
    \left(
    \begin{array}{c}
        x_1 \\
        x_2 \\
        \vdots \\
        x_K
    \end{array}
    \right)
    \\
    & =
    x_1 \left(
    \begin{array}{c}
        a_{11} \\
        a_{21} \\
        \vdots \\
        a_{N1} 
    \end{array}
    \right)
    +
    x_2 \left(
    \begin{array}{c}
        a_{12} \\
        a_{22} \\
        \vdots \\
        a_{N2} 
    \end{array}
    \right)
    + \cdots + 
    x_K \left(
    \begin{array}{c}
        a_{1K} \\
        a_{2K} \\
        \vdots \\
        a_{NK} 
    \end{array}
    \right)
    \\
    & = 
    \sum_{k=1}^K x_k \col_k \boldA
\end{aligned}
```

(ss-rank)=
### Matrices as Maps
One of the most useful ways to think about matrices is as maps from one vector space to another. In particular, an $N \times K$ matrix $\boldA$ can be thought of as a map sending a vector $\boldx \in \RR^K$ into a new vector $\boldy = \boldA \boldx$ in $\RR^N$. As the next theorem shows, these maps are always linear. In fact, they are the *only* linear functions. In other words, the set of linear functions from $\RR^K$ to $\RR^N$ and the set of $N \times K$ matrices are in one-to-one correspondence:

```{prf:theorem}
:label: t-lmaeq

Let $T$ be a function from $\RR^K$ to $\RR^N$. The following are equivalent:

(i) $T$ is linear.

(ii) There exists an $N \times K$ matrix $\boldA$ such that $T\boldx = \boldA \boldx$ for all $\boldx \in \RR^K$.
```


````{prf:proof}
:nonumber:

((i) $\implies$ (ii)) Let $T \colon \RR^K \to \RR^N$ be linear. We aim to construct a matrix $\boldA$ such that $T\boldx = \boldA \boldx$ for all $\boldx \in \RR^K$. As usual, let $\bolde_k$ be the $k$th canonical basis vector in $\RR^K$. Define an $N \times K$ matrix $\boldA$ by $\col_k \boldA = T\bolde_k$. Pick any $\boldx \in \RR^K$. We can also write $\boldx = \sum_{k=1}^K x_k \bolde_k$, where $x_k$ is the $k$th element of $\boldx$. By linearity, we have

```{math}
:enumerated: false

T\boldx 
        = \sum_{k=1}^K x_k T \bolde_k
        = \sum_{k=1}^K x_k \col_k \boldA
```

This is just $\boldA \boldx$, as shown in {ref}`ss-mbd`.

((ii) $\implies$ (i)) Fix $N \times K$ matrix $\boldA$ and consider the function $T \colon \RR^K \to \RR^N$ defined by $T \boldx = \boldA \boldx$. Pick any $\boldx$, $\boldy$ in $\RR^K$, and any scalars $\alpha$ and $\beta$. The rules of matrix arithmetic (see [fact %s](#fa-bmari)) tell us that

```{math}
:enumerated: false

T(\alpha \boldx + \beta \boldy) 
        := \boldA (\alpha \boldx + \beta \boldy)
        = \boldA \alpha \boldx  + \boldA \beta \boldy
        = \alpha \boldA \boldx  + \beta \boldA \boldy
        =: \alpha T\boldx + \beta T\boldy
```

Hence $T$ is linear. ◻
````


When we consider the problem of solving a system of linear equations such as $\boldA \boldx = \boldb$, the first issue we need to concern ourselves with is existence. Can we find an $\boldx$ that satisfies this equation, for any given $\boldb$? A bit of thought will convince you that this is the same question as: Is the corresponding linear map $T\boldx = \boldA \boldx$ an onto function? (See {ref}`ss-functions` for the definition.) Equivalently, is $\range T$ equal to all of $\RR^N$?

The range of $T$ is all vectors of the form $T\boldx = \boldA \boldx$ where $\boldx$ varies over $\RR^K$. We just saw in {ref}`ss-mbd` that, for any $\boldx \in \RR^K$, we have $\boldA \boldx = \sum_{k=1}^K x_k \col_k \boldA$. It follows that $\range T$ is equal to the **column space** of $\boldA$, which is by definition the span of the columns of $\boldA$. We represent it by the symbol

$$
\colspace \boldA := \Span \{ \col_1 \boldA, \ldots, \col_K \boldA \}
$$ (eq-secolsp)

To summarize the preceding discussion, we have

```{math}
:enumerated: false

\colspace \boldA = \range T = \setntn{\boldA \boldx}{ \boldx \in \RR^K}
```

How large is the column space of a given matrix? To answer that question, we have to say what "large\" means. In the context of linear algebra, size of subspaces is usually measured by dimension. The dimension of $\colspace \boldA$ is known as the **rank** of $\boldA$. That is,

```{math}
:enumerated: false

\rank \boldA := \dim \colspace \boldA
```

$\boldA$ is said to have **full column rank** if $\rank \boldA$ is equal to $K$, the number of its columns. The reason we say "full" rank here is that, by definition, $\colspace \boldA$ is the span of $K$ vectors. Hence, by part (i) of [theorem %s](#t-dimspan), we must have $\dim \colspace \boldA \leq K$. In other words, the rank of $\boldA$ is less than or equal to $K$. $\boldA$ is said to have full column rank when this maximum is achieved.

When is this maximum achieved? By part (ii) of [theorem %s](#t-dimspan), this will be the case precisely when the columns of $\boldA$ are linearly independent. Let's state this as a fact:

```{prf:property}
:label: fa-frie

Let $\boldA$ be an $N \times K$ matrix. The following statements are equivalent:

(i) $\boldA$ is of full column rank.

(ii) The columns of $\boldA$ are linearly independent.

(iii) $\boldA \boldx = \boldzero \; \implies \; \boldx = \boldzero$.
```


The last equivalence follows from [theorem %s](#t-fdlinin).

(ss-smi)=
### Square Matrices and Invertibility
Perhaps the most common problem in linear algebra is solving systems of linear equations. A generic representation is $\boldA \boldx = \boldb$, where $\boldx$ contains the unknowns and $\boldA$ and $\boldb$ are given. There are a variety of scenarios depending on the properties of $\boldA$. For now let's consider the case where $\boldA$ is $N \times N$, and seek conditions on $\boldA$ under which, for every $\boldb \in \RR^N$, there exists exactly one $\boldx \in \RR^N$ such that $\boldA \boldx = \boldb$.

The best way to understand this problem is as follows. Let $T$ be the linear map $T \boldx = \boldA \boldx$. The question we are asking here is: When does each point $\boldb \in \RR^N$ have one and only one preimage under $T$? In other words, when is $T$ a bijection?

To answer this question, we can refer back to the discussion in {ref}`ss-lfb`. We saw there that linear bijections are called nonsingular functions, and we discussed a number of equivalences for this property. The next fact replicates these equivalences in the language of matrices.

```{prf:property}
:label: fa-exinv

For $N \times N$ matrix $\boldA$, the following are equivalent:

(i) The columns of $\boldA$ are linearly independent.

(ii) The columns of $\boldA$ form a basis of $\RR^N$.

(iii) $\rank \boldA = N$.

(iv) $\colspace \boldA = \RR^N$.

(v) $\boldA \boldx = \boldA \boldy \; \implies \; \boldx = \boldy$.

(vi) $\boldA \boldx = \boldzero \; \implies \; \boldx = \boldzero$.

(vii) For each $\boldb \in \RR^N$, the equation $\boldA \boldx = \boldb$ has a solution.

(viii) For each $\boldb \in \RR^N$, the equation $\boldA \boldx = \boldb$ has a unique solution.
```


These results can all be verified using [theorem %s](#t-linfunc) by checking the corresponding implications for $T \boldx = \boldA \boldx$. For example, it's easy to see that if $\bolde_n$ is the $n$th canonical basis vector of $\RR^N$, then

```{math}
:enumerated: false

T\bolde_n = \boldA \bolde_n = \col_n \boldA
```

and hence the condition of linear independence of $\{T\bolde_1, \ldots, T\bolde_N\}$ in [theorem %s](#t-linfunc) translates to linear independence of the columns of $\boldA$.

Following common usage, if any of the equivalent conditions in [fact %s](#fa-exinv) are true we will call not just the map $T$ but also the underlying matrix $\boldA$ **nonsingular**. If any one—and hence all—of these conditions fail, then $\boldA$ is called **singular**.

Any bijection has an inverse (see {ref}`ss-functions`). In fact any nonsingular map $T$ has a nonsingular inverse $T^{-1}$ (see [fact %s](#fa-nsins)). In the present setting, where $T$ is generated by a matrix $\boldA$, the inverse $T^{-1}$ is also associated with a matrix, called the inverse of $\boldA$. The next theorem gives details:

```{prf:theorem}
:label: t-nonsigmat

For nonsingular $\boldA$ the following statements are true:

(i) There exists a square matrix $\boldB$ such that $\boldA \boldB = \boldB \boldA = \boldI$, where $\boldI$ is the identity matrix. The matrix $\boldB$ is called the **inverse** of $\boldA$, and written as $\boldA^{-1}$.

(ii) For each $\boldb \in \RR^N$, the unique solution to $\boldA \boldx = \boldb$ is given by

$$
\boldx = \boldA^{-1} \boldb
$$ (eq-sofox)
```


For this reason, nonsingular matrices are also referred to as **invertible**. The proof of [theorem %s](#t-nonsigmat) is a solved exercise ([ex. %s](#ex-nonsigmat)).

````{prf:example}

Consider the $N$ good linear demand system

```{math}
:enumerated: false

q_n = \sum_{k=1}^N a_{nk} p_k + b_n,
        \quad n = 1, \ldots, N
```

where $q_n$ and $p_n$ are quantity and price of the $n$th good respectively. We want to compute the inverse demand function, which gives prices in terms of quantities. To do so, we write our system in matrix form as $\boldq = \boldA \boldp + \boldb$. If the columns of $\boldA$ are linearly independent, then we can invert the system: a unique solution exists for each fixed $\boldq$ and $\boldb$. That solution is given by $\boldp = \boldA^{-1} (\boldq - \boldb)$.
````


As stated in the next fact, to show that $\boldA$ is invertible and $\boldB$ is the inverse of $\boldA$, it suffices to show that $\boldB$ is either a **left inverse**, in the sense that $\boldB \boldA = \boldI$, or a **right inverse**, in the sense that $\boldA \boldB = \boldI$. A proof is given in {ref}`ss-determ`.

```{prf:property}
:label: fa-lrinv

Let $\boldA$ and $\boldB$ be $N \times N$ square matrices. If $\boldB$ is either a left or a right inverse of $\boldA$, then $\boldA$ is nonsingular and $\boldB$ is its inverse.
```


The next fact collects more useful results about inverse matrices.

```{prf:property}
:label: fa-bpiv

If $\boldA$ and $\boldB$ are nonsingular and $\alpha \neq 0$, then

(i) $(\boldA^{-1})^{-1} = \boldA$,

(ii) $(\alpha \boldA)^{-1} = \alpha^{-1} \boldA^{-1}$, and

(iii) $(\boldA \boldB)^{-1} = \boldB^{-1} \boldA^{-1}$.
```


The relation $(\boldA \boldB)^{-1} = \boldB^{-1} \boldA^{-1}$ is just a special case of the analogous rule for inversion of bijections (see [fact %s](#fa-bij)). But you can also prove the equality directly by confirming that $\boldB^{-1} \boldA^{-1}$ is a right inverse (or left inverse) of $\boldA \boldB$.

(ss-determ)=
### Determinants
To each square matrix $\boldA$, we can associate a unique number $\det \boldA$ called the determinant of $\boldA$. To define it, let $S(N)$ be the set of all bijections from $\{1, \ldots, N\}$ to itself. For $\pi \in S(N)$ we define the **signature** of $\pi$ as

```{math}
:enumerated: false

\sgn(\pi) := \prod_{m < n} \frac{\pi(m) - \pi(n)}{m - n}
```

The **determinant** of $N \times N$ matrix $\boldA$ is then given as

```{math}
:enumerated: false

\det \boldA 
    := \sum_{\pi \in S(N)} \sgn(\pi) \prod_{n=1}^N a_{\pi(n) n}
```

We won't concern ourselves with the details of this definition. For now it's enough to know the following facts:

```{prf:property}
:label: fa-detfacts

If $\boldI$ is the $N \times N$ identity, $\boldA$ and $\boldB$ are $N \times N$ matrices and $\alpha \in \RR$, then

(i) $\det \boldI = 1$,

(ii) $\boldA$ is nonsingular if and only if $\det \boldA \neq 0$,

(iii) $\det(\boldA\boldB) = \det(\boldA) \det(\boldB)$,

(iv) $\det(\alpha \boldA) = \alpha^N \det(\boldA)$, and

(v) $\det(\boldA^{-1}) = (\det(\boldA))^{-1}$.
```


In the $2 \times 2$ case one can show that the determinant satisfies

$$
\det 
    \left(
    \begin{array}{cc}
        a & b  \\
        c & d  \\
    \end{array}
    \right)
    = ad - bc
$$ (eq-det22)

As an example of how these results can be useful, let's go back and prove [fact %s](#fa-lrinv). Fix square matrix $\boldA$ and suppose that a right inverse $\boldB$ exists, in the sense that $\boldA \boldB = \boldI$. This equality implies that both $\boldA$ and $\boldB$ are nonsingular, since applying $\det$ to both sides of $\boldA \boldB = \boldI$ and using the rules in [fact %s](#fa-detfacts) gives $\det(\boldA) \det(\boldB) = 1$. It follows that both $\det \boldA$ and $\det \boldB$ are nonzero. Hence both matrices are nonsingular.

To show that $\boldB$ is the inverse of $\boldA$, we just need to check that, in addition to $\boldA \boldB = \boldI$, we have $\boldB \boldA = \boldI$. To obtain the latter equality from the former, premultiply the former by $\boldB$ to get $\boldB \boldA \boldB = \boldB$ and then postmultiply by $\boldB^{-1}$ to get $\boldB \boldA = \boldI$. The proof for the left inverse case is similar.

## Properties of Matrices

Let's look at some special kinds of matrices and their role in linear algebra.

(ss-triang)=
### Diagonal and Triangular Matrices
A square matrix is called

- **lower triangular** if each element strictly above the principal diagonal is zero,

- **upper triangular** if every element strictly below the principal diagonal is zero, and

- **triangular** if it is either upper or lower triangular.

For example, if we define

```{math}
:label: eq-ltute
:enumerated: false

\boldL :=
    \left(
    \begin{array}{ccc}
        1 & 0 & 0  \\
        2 & 5 & 0 \\
        3 & 6 & 1
    \end{array}
    \right)
    \quad \text{and} \quad
    \boldU :=
    \left(
    \begin{array}{ccc}
        1 & 2 & 3  \\
        0 & 5 & 6 \\
        0 & 0 & 1
    \end{array}
    \right)
```

then $\boldL$ and $\boldU$ are lower and upper triangular respectively. The great advantage of triangular matrices is that the associated linear equations are simple to solve using either forward or backward substitution. For example, with the system

```{math}
:enumerated: false

\left(
    \begin{array}{ccc}
        1 & 0 & 0  \\
        2 & 5 & 0 \\
        3 & 6 & 1
    \end{array}
    \right)
    \left(
    \begin{array}{ccc}
        x_1 \\
        x_2 \\
        x_3 
    \end{array}
    \right)
    =
    \left(
    \begin{array}{c}
        x_1   \\
        2x_1 + 5x_2 \\
        3x_1 + 6x_2 + x_3
    \end{array}
    \right)
    = 
    \left(
    \begin{array}{c}
        b_1 \\
        b_2 \\
        b_3 \\
    \end{array}
    \right)
```

the top equation involves only $x_1$, so we can solve for its value directly. Plugging this value into the second equation, we can solve for $x_2$ and so on.

```{prf:property}
:label: fa-dtmat

If $\boldA = (a_{mn})$ is triangular, then $\det \boldA = \prod_{n=1}^N a_{nn}$.
```


An even better situation arises when our matrix is *both* upper and lower triangular. These are the **diagonal** matrices. In other words, a square matrix $\boldA$ is diagonal if all entries off the principal diagonal are zero. For example, the identity matrix is diagonal.

The following notation is often used to define diagonal matrices:

```{math}
:enumerated: false

\diag(d_1, \ldots, d_N) := 
    \left(
    \begin{array}{cccc}
        d_1 & 0 & \cdots & 0 \\
        0 & d_2 & \cdots & 0 \\
        \vdots & \vdots &  & \vdots \\
        0 & 0 & \cdots & d_N \\
    \end{array}
    \right)
```

With diagonal matrices it is trivial to compute powers, roots, inverses and products:

```{prf:property}
:label: fa-powdm

Let $\boldC = \diag(c_1, \ldots, c_N)$ and $\boldD = \diag(d_1, \ldots, d_N)$. The following statements are true:

(i) $\boldC + \boldD = \diag(c_1 + d_1, \ldots, c_N + d_N)$.

(ii) $\boldC \boldD = \diag(c_1 d_1, \ldots, c_N d_N)$.

(iii) $\boldD^k = \diag(d^k_1, \ldots, d^k_N)$ for any $k \in \NN$.

(iv) If $d_n \geq 0$ for all $n$, then $\sqrt{\boldD}$ exists and equals $\diag(\sqrt{d_1}, \ldots, \sqrt{d_N})$.

(v) If $d_n \neq 0$ for all $n$, then $\boldD$ is nonsingular and $\boldD^{-1} = \diag(d_1^{-1}, \ldots, d_N^{-1})$.
```


You can check parts (i) and (ii) directly. The other claims follow from (i) and (ii).

### Trace, Transpose, and Symmetry

The **trace** of an $N \times N$ matrix $\boldA$ is the sum of the elements on its principal diagonal:

```{math}
:enumerated: false

\trace \boldA = \sum_{n=1}^N a_{nn}
```

````{prf:property}
:label: fa-lmpt

If $\boldA$ and $\boldB$ are $N \times N$ matrices and $\alpha$ and $\beta$ are two scalars, then

```{math}
:enumerated: false

\trace(\alpha \boldA + \beta \boldB) 
        = \alpha \trace(\boldA) + \beta \trace(\boldB)
```

Moreover, if $\boldA$ is $N \times M$ and $\boldB$ is $M \times N$, then $\trace(\boldA \boldB) = \trace(\boldB \boldA)$.
````


The **transpose** of $N \times K$ matrix $\boldA$ is a $K \times N$ matrix $\boldA^\T$ such that $\col_n(\boldA^\T) = \row_n \boldA$. For example,

$$
\boldA := 
    \left(
    \begin{array}{cc}
        10 & 40  \\
        20 & 50  \\
        30 & 60
    \end{array}
    \right)
    \quad \implies \quad
    \boldA^\T = 
    \left(
    \begin{array}{ccc}
        10 & 20 & 30 \\
        40 & 50 & 60 
    \end{array}
    \right)
$$ (eq-aandb)

and

```{math}
:enumerated: false

\boldB := 
    \left(
    \begin{array}{ccc}
        1 & 3 & 5 \\
        2 & 4 & 6 \\
    \end{array}
    \right)
    \quad \implies \quad
    \boldB^\T = 
    \left(
    \begin{array}{cc}
        1 & 2  \\
        3 & 4  \\
        5 & 6 
    \end{array}
    \right)
```

A square matrix $\boldA$ is called **symmetric** if $\boldA^\T = \boldA$, or, equivalently, $a_{nk} = a_{kn}$ for every $k$ and $n$. Note that $\boldA^\T \boldA$ and $\boldA \boldA^\T$ are always well-defined and symmetric.

```{prf:property}
:label: fa-owtra

For conformable matrices $\boldA$ and $\boldB$, transposition satisfies

(i) $(\boldA^\T)^\T = \boldA$,

(ii) $(\boldA \boldB)^\T = \boldB^\T \boldA^\T$,

(iii) $(\boldA + \boldB)^\T = \boldA^\T +  \boldB^\T$, and

(iv) $(c \boldA)^\T = c \boldA^\T$ for any constant $c$.
```


```{prf:property}
:label: fa-owtra2

For each square matrix $\boldA$, we have

(i) $\trace(\boldA) = \trace(\boldA^\T)$ and

(ii) $\det(\boldA^\T) = \det(\boldA)$.

(iii) If $\boldA$ is nonsingular, then so is $\boldA^\T$, and its inverse is $(\boldA^\T)^{-1}= (\boldA^{-1})^\T$.
```


Note that if $\bolda$ and $\boldb$ are $N \times 1$ vectors, then the matrix product $\bolda^\T \boldb = \boldb^\T \bolda$ is equal to $\sum_{n=1}^N a_n b_n$, which is the same as the inner product $\inner{\bolda, \boldb}$. In what follows we'll often work with column vectors and use matrix product rather than inner product notation.

(ss-eiei)=
### Eigenvalues and Eigenvectors
Let $\boldA$ be $N \times N$. As in {ref}`ss-rank`, think of $\boldA$ as a linear map, so that $\boldA \boldx$ is the image of $\boldx$ under $\boldA$. In general, $\boldA$ will map $\boldx$ to some arbitrary new location but sometimes $\boldx$ will only be *scaled*. That is,

$$
\boldA \boldx = \lambda \boldx
$$ (eq-eiei)

for some scalar $\lambda$. If $\boldx$ and $\lambda$ satisfy {eq}`eq-eiei` and $\boldx$ is nonzero, then $\boldx$ is called an **eigenvector** of $\boldA$, $\lambda$ is called an **eigenvalue**, and $(\boldx, \lambda)$ is called an **eigenpair**.

````{prf:example}

If

```{math}
:enumerated: false

\boldA :=
        \begin{pmatrix}
            1 & -1 \\
            3 & 5
        \end{pmatrix},
        \quad
        \boldx :=
        \begin{pmatrix}
            1 \\
            -1
        \end{pmatrix}
        \quad \text{and} \quad
        \lambda := 2
```

then $(\boldx, \lambda)$ is an eigenpair of $\boldA$, since $\boldx \neq \boldzero$ and

```{math}
:enumerated: false

\boldA \boldx =
        \begin{pmatrix}
            1 & -1 \\
            3 & 5
        \end{pmatrix}
        \begin{pmatrix}
            1 \\
            -1
        \end{pmatrix}
        =
        \begin{pmatrix}
            2 \\
            -2
        \end{pmatrix}
        = 2
        \begin{pmatrix}
            1 \\
            -1
        \end{pmatrix}
        =
        \lambda \boldx
```
````


```{prf:example}

If $\boldI$ is the $N \times N$ identity then $(\boldx, 1)$ is an eigenpair of $\boldI$ for every nonzero $\boldx \in \RR^N$.
```


Consider now the matrix

```{math}
:enumerated: false

\boldR := 
    \left(
    \begin{array}{cc}
        0 & -1  \\
        1 & 0  \\
    \end{array}
    \right)
```

that rotates any point counterclockwise by $90^{\circ}$, as shown in [figure %s](#f-rotation_1). This rotation causes the scaling in {eq}`eq-eiei` to fail for any $\lambda \in \RR$ and nonzero $\boldx \in \RR^2$. If, however, we admit the possibility that $\lambda$ and the elements of $\boldx$ can be complex, then {eq}`eq-eiei` can hold. For example, direct calculation confirms that $\lambda = i$ and $\boldx = (1, -i)^\T$ is an eigenpair for $\boldR$. It turns out that contemplation of complex eigenpairs is useful. Here eigenpairs are taken to be complex-valued unless explicitly stated to be real.

```{figure} figures/rotation_1.svg
:name: f-rotation_1

$\boldR$ rotates points by $90^{\circ}$
```

````{prf:property}
:label: fa-chei

Given $N \times N$ matrix $\boldA$, the scalar $\lambda$ is an eigenvalue of $\boldA$ if and only if

```{math}
:enumerated: false

\det(\boldA - \lambda \boldI) = 0
```
````


Here $\boldI$ is the $N \times N$ identity. {prf:ref}`ex-chei` asks you to confirm this important fact. For the $2 \times 2$ matrix in {eq}`eq-det22`, the rule for the $2 \times 2$ determinant in {eq}`eq-det22`, [fact %s](#fa-chei) and a little bit of algebra imply that its eigenvalues are given by the two roots of the polynomial expression

```{math}
:enumerated: false

\lambda^2 - (a + d) \lambda + (ad - bc) = 0
```

More generally, given any $N \times N$ matrix $\boldA$, it can be shown via the fundamental theorem of algebra that there exist complex numbers $\lambda_1, \ldots, \lambda_N$, not necessarily distinct, such that

$$
\det(\boldA - \lambda \boldI) = \prod_{n=1}^N (\lambda_n - \lambda)
$$ (eq-loei)

It is clear that each $\lambda_n$ satisfies $\det(\boldA - \lambda_n \boldI) = 0$ and hence is an eigenvalue of $\boldA$. In particular, $\lambda_1, \ldots, \lambda_N$ is the set of eigenvalues of $\boldA$, although it's worth repeating that these numbers are not necessarily distinct.

```{prf:property}
:label: fa-eiprop

Let $\boldA$ be $N \times N$ and let $\lambda_1, \ldots, \lambda_N$ be the eigenvalues defined in {eq}`eq-loei`. The following statements are true:

(i) $\det \boldA = \prod_{n=1}^N \lambda_n$.

(ii) $\trace \boldA = \sum_{n=1}^N \lambda_n$.

(iii) If $\boldA$ is symmetric, then $\lambda_n \in \RR$ for all $n$.

(iv) If $\boldA$ is nonsingular, then the eigenvalues of $\boldA^{-1}$ are $1/\lambda_1, \ldots, 1/\lambda_N$.

(v) If $\boldA$ is triangular, then its eigenvalues coincide with the elements on the principle diagonal.
```


It is immediate from (i) that $\boldA$ is nonsingular $\iff$ all its eigenvalues are nonzero.

(ss-qaudf)=
### Quadratic Forms
In statistics and econometrics we often encounter quadratic expressions. In general, given symmetric $N \times N$ matrix $\boldA$, the **quadratic function** or **quadratic form** on $\RR^N$ associated with $\boldA$ is the map $Q$ defined by

```{math}
:enumerated: false

Q(\boldx) := \boldx^\T \boldA \boldx = \sum_{j=1}^N \sum_{i=1}^N a_{ij} x_i x_j
```

To give a simple illustration, let $N = 2$ and let $\boldA$ be the identity matrix $\boldI$. In this case,

```{math}
:enumerated: false

Q(\boldx) = \| \boldx \|^2 = x_1^2 + x_2^2
```

A 3D graph of this function is shown in [figure %s](#f-qform_pd).

```{figure} figures/qform_pd.svg
:name: f-qform_pd

Quadratic function $Q(\boldx) = x_1^2 + x_2^2$
```

One thing you'll notice about this function is that its graph lies everywhere above zero, or $Q(\boldx) \geq 0$. In fact we know that $\| \boldx \|^2$ is nonnegative and will be zero only when $\boldx = \boldzero$. Hence the graph touches zero only at the point $\boldx = \boldzero$. Many other choices of $\boldA$ yield a quadratic form with this property. Such $\boldA$ are said to be positive definite. More generally, an $N \times N$ symmetric matrix $\boldA$ is called

- **nonnegative definite** if $\boldx^\T \boldA \boldx \geq 0$ for all $\boldx \in \RR^N$,

- **positive definite** if $\boldx^\T \boldA \boldx > 0$ for all $\boldx \in \RR^N$ with $\boldx \neq \boldzero$,

- **nonpositive definite** if $\boldx^\T \boldA \boldx \leq 0$ for all $\boldx \in \RR^N$, and

- **negative definite** if $\boldx^\T \boldA \boldx < 0$ for all $\boldx \in \RR^N$ with $\boldx \neq \boldzero$.

If $\boldA$ fits none of these categories, then $\boldA$ is called **indefinite**. {numref}`f-qform_nd` shows the graph of a negative definite quadratic function. Now the function is hill-shaped, and $\boldzero$ is the unique global maximum. {numref}`f-qform_indef` shows an indefinite form.

```{figure} figures/qform_nd.svg
:name: f-qform_nd

Quadratic function $Q(\boldx) = -x_1^2 - x_2^2$
```

```{figure} figures/qform_indef.svg
:name: f-qform_indef

Quadratic function $Q(\boldx) = x_1^2/2 + 8 x_1 x_2 + x_2^2/2$
```

The easiest case for detecting definiteness is when the matrix $\boldA$ is diagonal, since

```{math}
:enumerated: false

\boldA = \diag(d_1, \ldots, d_N)
    \quad \text{implies} \quad
    Q(\boldx) = d_1 x_1^2 + \cdots + d_N x_N^2
```

From the right-hand expression we see that a diagonal matrix is positive definite if and only if all diagonal elements are positive. Analogous statements are true for nonnegative, nonpositive and negative definite matrices. The next fact generalizes this idea and is proved in {ref}`ss-simmat`.

```{prf:property}
:label: fa-eigdef

Let $\boldA$ be any symmetric matrix. $\boldA$ is

(i) positive definite if and only if its eigenvalues are all positive,

(ii) negative definite if and only if its eigenvalues are all negative,

and similarly for nonpositive and nonnegative definite.
```


It follows from [fact %s](#fa-eigdef) that

```{prf:property}
:label: fa-ipde

If $\boldA$ is positive definite, then $\boldA$ is nonsingular, with $\det \boldA > 0$.
```


Finally, here's a necessary (but not sufficient) condition for each kind of definiteness.

```{prf:property}
:label: fa-ipde0

If $\boldA$ is positive definite, then each element $a_{nn}$ on the principal diagonal is positive, and the same for nonnegative, nonpositive and negative.
```


(s-pd)=
## Projection and Decomposition
This section collects some essential results on projection and decomposition of matrices. The projection material takes our abstract projection theory from {ref}`s-orthog` and translates it into the more concrete language of matrices.

(s-pontob)=
### Projection Matrices
As stated in [theorem %s](#t-opt2), given any subspace $S$, the corresponding projection $\boldP = \proj S$ is a linear map from $\RR^N$ to $\RR^N$. In view of [theorem %s](#t-lmaeq), it follows that there exists an $N \times N$ matrix $\hat \boldP$ such that $\boldP \boldx = \hat \boldP \boldx$ for all $\boldx \in \RR^N$. In fact we've anticipated this in the notation $\boldP$, and from now on $\boldP$ will also represent the corresponding matrix. But what does this matrix look like?

```{prf:theorem}
:label: t-projon

Let $S$ be a subspace of $\RR^N$. If $\boldP = \proj S$, then

$$
\boldP = \boldB (\boldB^\T\boldB)^{-1} \boldB^\T
$$ (eq-projb)
for every matrix $\boldB$ such that the columns of $\boldB$ form a basis of $S$.
```


This result generalizes [fact %s](#fa-projon), which pertains to projection onto the span of an orthonormal basis. The expression in {eq}`eq-projb` is more complex, but at the same time it applies to any basis, orthonormal or otherwise. We further explore the connection between the two results in {ref}`ss-gso`. For the proof of [theorem %s](#t-projon) see [exercise %s](#ex-projon) and its solution.

The construction of $\boldP$ in {eq}`eq-projb` implicitly assumes that $\boldB^\T \boldB$ is nonsingular. This is justified because $\boldB$ has full column rank (see [ex. %s](#ex-fr)). As usual, we let $\boldM = \boldI - \boldP$ denote the residual projection.

```{prf:example}
:label: eg-pvonepm

Recall [example %s](#eg-pvones), where we found that the projection of $\boldy \in \RR^N$ onto $\Span\{\boldone\}$ is $\bar y \boldone$. We can get this from {eq}`eq-projb` as well. Since $\boldone$ is a basis for $\Span\{\boldone\}$, we have

$$
\boldP = \proj \; \Span\{\boldone\} 
        \; \implies \;
        \boldP 
        = \boldone (\boldone^\T\boldone)^{-1} \boldone^\T  
        = \frac{1}{N} \boldone \boldone^\T
$$ (eq-pcz)
This leads us back to $\boldP \boldy = \bar y \boldone$, as expected. The corresponding residual projection is

$$
\boldM_c
        = \boldI - \frac{1}{N} \boldone \boldone^\T
$$ (eq-pczm)
The reason for the subscripted $c$ was discussed in [example %s](#eg-pvonesr).
```


```{prf:property}
:label: fa-anox

In the setting of [theorem %s](#t-projon), we have

(i) $\boldM \boldB = \boldzero$

(ii) $\boldP \boldB = \boldB$
```


For example, it's easy to see that $\boldM_c$ in {eq}`eq-pczm` maps $\boldone$ to $\boldzero$. {prf:ref}`ex-anox` asks you to prove [fact %s](#fa-anox).

A square matrix $\boldA$ is called **idempotent** if $\boldA \boldA = \boldA$.

```{prf:property}
:label: fa-pmsi

Both $\boldP$ and $\boldM$ are symmetric and idempotent.
```


Idempotence of $\boldP$ and $\boldM$ can be checked by direct calculation. A better understanding is obtained by reflecting on the fact that projecting onto a subspace twice is the same as projecting once. After one projection the vector is already in the subspace. See [fact %s](#fa-opt3).

```{prf:property}
:label: fa-idtrra

If $\boldA$ is any idempotent matrix, then $\rank \boldA = \trace \boldA$.
```


For orthogonal projections, we can say more:

```{prf:property}
:label: fa-trop

Let $S$ be a linear subspace of $\RR^N$. If $\boldP = \proj S$ and $\boldM$ is the residual projection, then

(i) $\rank \boldP  = \trace \boldP  = \dim S$ and

(ii) $\rank \boldM  = \trace \boldM  = N - \dim S$.
```


To see why $\rank \boldP = \dim S$, recall that the rank of a linear map is the dimension of its range. When $\boldP = \proj S$, the range of $\boldP$ is exactly $S$.

To show that $\trace \boldP = \dim S$ also holds, we can appeal to [fact %s](#fa-idtrra). There's also a nice direct proof. See [exercise %s](#ex-dptpek) and its solution. Once we know that $\trace\boldP = \dim S$, it's clear that $\trace\boldM = N - \dim S$ because

```{math}
:enumerated: false

\trace\boldM 
    = \trace(\boldI - \boldP) 
    = \trace \boldI - \trace\boldP 
    = N - \dim S
```

(ss-odse)=
### Overdetermined Systems of Equations
In {ref}`ss-smi` we talked about solving equations of the form $\boldA \boldx = \boldb$ when $\boldA$ is square. In statistics and econometrics, we often work the case where $\boldA$ is $N \times K$ and $K \leq N$. When the inequality is strict, the system of equations is said to be **overdetermined**.

Consider the problem of whether or not there exists a vector $\boldx$ satisfying $\boldA \boldx = \boldb$ in this setting. Intuitively, when the number of equations is larger than the number of unknowns ($N > K$) we may not be able to find an $\boldx$ that satisfies all $N$ equations. There are several equivalent ways to formalize this intuition. The linear map $T \colon \RR^K \to \RR^N$ corresponding to $\boldA$ is $T\boldx = \boldA \boldx$ (see {ref}`ss-rank`). We know the following statements to be equivalent:

(i) there exists an $\boldx \in \RR^K$ with $\boldA \boldx = \boldb$.

(ii) $\boldb \in \colspace \boldA$.

(iii) $\boldb \in \range T$.

We also know from [theorem %s](#t-lfoc) that when $K < N$, the function $T$ cannot be onto, and hence it's possible that $\boldb$ lies outside the range of $T$.

In fact we can say more than this. When $K < N$, the scenario $\boldb \in \colspace \boldA$ is in some sense very rare. The reason is that the dimension of $\colspace \boldA$, which is precisely the rank of $\boldA$, is less than or equal to $K$ (see {ref}`ss-rank`) and hence strictly less than $N$. There is a sense in which $K$-dimensional subspaces of $\RR^N$ are negligible, and hence the "chance" of $\boldb$ happening to lie in this subspace is tiny.[^1]

As a result the standard approach is to admit that an exact solution may not exist, and to focus instead on finding a $\boldx \in \RR^K$ such that $\boldA \boldx$ is as close to $\boldb$ as possible. It's natural that "close to\" is defined in terms of ordinary Euclidean norm, which leads us to the minimization problem

$$
\min_{\boldx \in \RR^K} \| \boldb - \boldA \boldx \|
$$ (eq-mxls)

This is called a **least squares problem** because solving {eq}`eq-mxls` is the same as minimizing $\| \boldb - \boldA \boldx \|^2$ with respect to $\boldx$, and the squared norm is, by definition, a sum of squares. Assuming as before that $\boldA$ is $N \times K$ with $K \leq N$ and $\boldb$ is $N \times 1$, we can use the orthogonal projection theorem to solve {eq}`eq-mxls` as follows:

```{prf:theorem}
:label: t-lssol

If $\boldA$ has full column rank, then {eq}`eq-mxls` has the unique solution

$$
\hboldx := (\boldA^\T \boldA)^{-1} \boldA^\T \boldb
$$ (eq-dlse)
```


````{prf:proof}
:nonumber:

Let $\boldA$ and $\boldb$ be as in the statement of the theorem. Let $\hboldx$ be as in {eq}`eq-dlse` and let $S := \colspace \boldA$. By the full column rank assumption, the columns of $\boldA$ form a basis for $S$. Hence, applying [theorem %s](#t-projon), the orthogonal projection of $\boldb$ onto $S$ is

$$
\boldP\boldb 
        := \boldA (\boldA^\T \boldA)^{-1} \boldA^\T \boldb  
         = \boldA \hboldx
$$ (eq-ibprojon)
Moreover, since the orthogonal projection theorem gives a unique minimizer in terms of the closest point in $S$ to $\boldb$, we must have

$$
\| \boldb - \boldA\hboldx \| < \| \boldb - \boldy\|
        \quad \text{for all } \boldy \in S, \; \boldy \neq \boldA \hboldx
$$ (eq-sfy)
Pick any $\boldx \in \RR^K$ such that $\boldx \neq \hboldx$. By the definition of $S$ we have $\boldA \boldx \in S$. In addition, since $\boldx \neq \hboldx$, and since $\boldA$ has full column rank, it must be that $\boldA \boldx \neq \boldA \hboldx$ ([ex. %s](#ex-ui)). Hence

```{math}
:enumerated: false

\| \boldb - \boldA\hboldx \| < \| \boldb - \boldA \boldx\|
        \quad \text{for all } \boldx \in \RR^K, \; \boldx \neq \hboldx
```

In other words, $\hboldx$ is the unique solution to {eq}`eq-mxls`. ◻
````


In the expression for $\hboldx$ in {eq}`eq-dlse`, the matrix $(\boldA^\T \boldA)^{-1} \boldA^\T$ is called the **pseudoinverse** of $\boldA$. If $K = N$—that is, if $\boldA$ is in fact square—then, under our full rank assumption, the pseudoinverse reduces to the inverse $\boldA^{-1}$, and the least squares solution $\hboldx$ in {eq}`eq-dlse` reduces to the expression given in {eq}`eq-sofox`.

What happens if we drop the assumption that the columns of $\boldA$ are linearly independent? The set $\colspace \boldA$ is still a linear subspace, and the orthogonal projection theorem still gives us a closest point $\boldP \boldb$ to $\boldb$ in $\colspace \boldA$. Since $\boldP \boldb \in \colspace \boldA$, there still exists a vector $\hboldx$ such that $\boldP \boldb = \boldA \hboldx$. The problem is that now there exists an infinity of such vectors. {prf:ref}`ex-lsiosol` asks you to prove this.

(ss-gso)=
### QR Decomposition
The **QR decomposition** of a given matrix $\boldA$ is a product of the form $\boldQ \boldR$, where the first matrix has orthonormal columns and the second is upper triangular. This factorization has many applications, including least squares problems and the computation of eigenvalues. The decomposition is based on Gram–Schmidt orthogonalization.

```{prf:theorem}
:label: t-qr

If $\boldA$ is an $N \times K$ matrix with full column rank, then there exists a factorization $\boldA = \boldQ \boldR$ where

(i) $\boldR$ is $K \times K$, upper triangular and nonsingular, and

(ii) $\boldQ$ is $N \times K$, with orthonormal columns.
```


````{prf:proof}
:nonumber:

Let $\boldA$ be as stated, and let $\bolda_k := \col_k \boldA$. {prf:ref}`t-gso` gives us existence of an orthonormal set $\{\boldu_1, \ldots, \boldu_K\}$ such that the span of $\{\boldu_1, \ldots, \boldu_k\}$ equals that of $\{\bolda_1, \ldots, \bolda_k\}$ for $k = 1, \ldots, K$. In particular, $\bolda_k$ is in the span of $\{\boldu_1, \ldots, \boldu_k\}$, and hence, appealing to [fact %s](#fa-bdons), we can write

```{math}
:enumerated: false

\begin{aligned}
\bolda_1 & = (\bolda_1^\T \boldu_1) \boldu_1
        \\
        \bolda_2 & = (\bolda_2^\T \boldu_1) \boldu_1 + (\bolda_2^\T \boldu_2) \boldu_2
        \\
        \bolda_3 & = (\bolda_3^\T \boldu_1) \boldu_1 
            + (\bolda_3^\T \boldu_2) \boldu_2
            + (\bolda_3^\T \boldu_3) \boldu_3
\end{aligned}
```

and so on. Sticking to the $3 \times 3$ case to simplify expressions, we can stack these equations horizontally to get

```{math}
:enumerated: false

\left(
        \begin{array}{ccc}
            |        & |        & |        \\
            \bolda_1 & \bolda_2 & \bolda_3 \\
            |        & |        & |        
        \end{array}
        \right)
        =
        \left(
        \begin{array}{ccc}
            |        & |        & |        \\
            \boldu_1 & \boldu_2 & \boldu_3 \\
            |        & |        & |        
        \end{array}
        \right)
        \left(
        \begin{array}{ccc}
            (\bolda_1^\T \boldu_1)   & (\bolda_2^\T \boldu_1)  & (\bolda_3^\T \boldu_1)       \\
            0 & (\bolda_2^\T \boldu_2) & (\bolda_3^\T \boldu_2) \\
            0 & 0      & (\bolda_3^\T \boldu_3) 
        \end{array}
        \right)
```

or $\boldA = \boldQ \boldR$. This is our QR decomposition.

It remains to show that $\boldR$ is invertible. This will be so if each term $\bolda_k^\T \boldu_k$ is nonzero, since the determinant is the product of these elements (see [fact %s](#fa-eiprop)). Suppose, to the contrary, that $\bolda_k^\T \boldu_k = 0$ for some $k$. Then $\bolda_k$ lies in the span of $\{\boldu_1, \ldots, \boldu_{k-1}\}$, which, by construction, agrees with the span of $\{\bolda_1, \ldots, \bolda_{k-1}\}$. This contradicts linear independence of the columns of $\boldA$. ◻
````


Given the decomposition $\boldA = \boldQ \boldR$, the least squares solution $\hboldx$ defined in {eq}`eq-dlse` can also be written as $\hboldx = \boldR^{-1} \boldQ^\T \boldb$ ([ex. %s](#ex-qrlsq)). Premultiplying by $\boldR$ converts this expression to $\boldR \hboldx = \boldQ^\T \boldb$, which is easy to solve because $\boldR$ is triangular (see {ref}`ss-triang`).

(ss-simmat)=
### Diagonalization and Spectral Theory
One important concept in dynamic systems and related fields is topological conjugacy. If $f \colon A \to A$ and $g \colon B \to B$, then $g$ is said to be **topologically conjugate** to $f$ whenever there exists a continuous bijection $\tau \colon B \to A$ such that $f = \tau \circ g \circ \tau^{-1}$. The idea is that the action of $f$ can be replicated by transporting a point into the domain of $g$, applying $g$, and then transporting it back. This can be beneficial if $g$ is somehow simpler than $f$.

In the case of linear maps—that is, matrices—it is natural to study conjugacy in a setting where the bijection is also required to be linear. In linear algebra this is called *similarity*. In particular, a square matrix $\boldA$ is said to be **similar** to another matrix $\boldB$ if there exists an invertible matrix $\boldP$ such that $\boldA = \boldP \boldB \boldP^{-1}$. {numref}`f-diagonalize` shows the conjugate relationship of the two matrices when thought of as maps.

```{figure} figures/tikz_diagonalize.svg
:name: f-diagonalize

$\boldA$ is similar to $\boldB$
```

The next fact is a fun exercise.

```{prf:property}
:label: fa-powsm

If $\boldA$ is similar to $\boldB$, then $\boldA^t$ is similar to $\boldB^t$ for all $t \in \NN$.
```


As discussed above, similarity of $\boldA$ to a given matrix $\boldB$ is most useful when $\boldB$ is somehow simpler than $\boldA$, or more amenable to a given operation. The simplest kind of matrices we work with are diagonal matrices, so similarity to a diagonal matrix is particularly desirable. If $\boldA$ is similar to a diagonal matrix, then $\boldA$ is called **diagonalizable**.

```{prf:example}
:label: eg-eapo

Suppose that we want to calculate $\boldA^t$ for some given $t \in \NN$. If $\boldA = \boldP \boldLambda \boldP^{-1}$ for some $\boldLambda = \diag(\lambda_1, \ldots, \lambda_N)$, then by [fact %s](#fa-powsm) and [fact %s](#fa-powdm), we have $\boldA^t = \boldP \diag(\lambda_1^t, \ldots, \lambda_N^t) \boldP^{-1}$. Aside from mapping backward and forward with $\boldP$, the only effort required to take the $t$th power of $\boldA$ is to take the $t$th power of $N$ scalars.
```


If $\boldA = \boldP \boldLambda \boldP^{-1}$ where $\boldLambda$ is diagonal, it follows immediately that the elements on the principal diagonal of $\boldLambda$ are the eigenvalues of $\boldA$, and that the columns of $\boldP$ are eigenvectors. To see this, observe that $\boldA = \boldP \boldLambda \boldP^{-1}$ implies $\boldA \boldP = \boldP \boldLambda$. Equating the $n$th column on each side gives $\boldA \boldp_n = \lambda_n \boldp_n$, where $\boldp_n := \col_n \boldP$. Finally, $\boldp_n$ is not the zero vector because otherwise $\boldP$ would not be invertible. Let's summarize:

```{prf:property}
:label: fa-diagee

If $\boldA = \boldP \boldLambda \boldP^{-1}$ for some $\boldLambda = \diag(\lambda_1, \ldots, \lambda_N)$, then $(\col_n \boldP, \lambda_n)$ is an eigenpair of $\boldA$ for each $n$.
```


But when is $\boldA$ diagonalizable? The only thing we need to make the expression $\boldA = \boldP \boldLambda \boldP^{-1}$ work is for $\boldP$ to be invertible. This explains the next fact.

```{prf:property}

An $N \times N$ matrix $\boldA$ is diagonalizable if and only if it has $N$ linearly independent eigenvectors.[^2]
```


If $\boldA$ does have $N$ linearly independent eigenvectors, then we are in good shape—diagonalization is possible. But can we do any better? The matrix $\boldLambda$ cannot really be simplified, but in some cases $\boldP$ can be. In particular, things are even nicer if $\boldP$ has orthonormal columns—that is, its columns form an orthonormal set, and hence an orthonormal basis of $\RR^N$. These kinds of matrices are called **orthogonal matrices**. Here are some nice properties of such matrices.

```{prf:property}
:label: fa-orthmat

If $\boldQ$ and $\boldP$ are $N \times N$ orthogonal matrices, then

(i) $\boldQ^\T$ is orthogonal and $\boldQ^{-1} = \boldQ^\T$,

(ii) $\boldQ \boldP$ is orthogonal, and

(iii) $\det \boldQ \in \{-1, 1\}$.
```


The first result tells us is that if $\boldA = \boldQ \boldLambda \boldQ^{-1}$ and $\boldQ$ has orthonormal columns, then $\boldA = \boldQ \boldLambda \boldQ^\T$. It's easy to see from this expression that $\boldA$ must be symmetric, so if we hope to diagonalize with this extra orthonormal property, then we shouldn't look beyond symmetric matrices. The following fundamental theorem tells us that this form of diagonalization is possible precisely when $\boldA$ is symmetric.

```{prf:theorem}
:label: t-dismat

If $\boldA$ is symmetric, then $\boldA$ can be diagonalized as $\boldA = \boldQ \boldLambda \boldQ^\T$, where $\boldQ$ is an orthogonal matrix and $\boldLambda$ is the diagonal matrix formed from the eigenvalues of $\boldA$.
```


This is one version of the **spectral decomposition theorem**. See, for example, §10.3 of {cite}`jan`. $\boldQ \boldLambda \boldQ^\T$ is called the **symmetric eigenvalue decomposition** of $\boldA$. It's not hard to see that its action on a given $N \times 1$ vector $\boldx$ can be written as

```{math}
:enumerated: false

\boldA\boldx = \sum_{n=1}^N \lambda_n (\boldu_n^\T \boldx) \boldu_n
```

where $\lambda_n$ is the $n$th eigenvalue of $\boldA$ and $\boldu_n = \col_n \boldQ$. Compare with $\boldx = \sum_{n=1}^N (\boldu_n^\T \boldx) \boldu_n$, which is true by [fact %s](#fa-bdons).

One nice application of the spectral theorem is a proof of [fact %s](#fa-eigdef), which states among other things that symmetric matrix $\boldA$ is positive definite if and only if its eigenvalues are all positive. {prf:ref}`ex-eigdef` and its solution step you through the arguments.

```{prf:property}
:label: fa-sqrtsnd

If $\boldA$ is nonnegative definite, then $\sqrt{\boldA}$ exists and equals $\boldQ \sqrt{\boldLambda} \boldQ^\T$. The matrix $\sqrt{\boldLambda}$ is given by $\diag(\sqrt{\lambda_1}, \ldots, \sqrt{\lambda_N})$.
```


Here $\boldQ \boldLambda \boldQ^\T$ is the spectral decomposition of $\boldA$. To check that $\boldQ \sqrt{\boldLambda} \boldQ^\T$ is a square root just multiply. $\sqrt{\boldLambda}$ exists because $\boldA$ is nonnegative definite, and hence its eigenvalues are nonnegative (see [fact %s](#fa-eigdef) and [fact %s](#fa-powdm)).

By combining our results on spectral decomposition with the QR decomposition discussed in {ref}`ss-gso`, we can prove the following well-known fact:

```{prf:property}
:label: fa-cholesky

If $\boldA$ is positive definite, then there exists a nonsingular, upper triangular matrix $\boldR$ such that $\boldA = \boldR^\T \boldR$.
```


This decomposition is called the **Cholesky decomposition**. The proof is obtained by writing

```{math}
:enumerated: false

\boldA 
    = \boldQ \boldLambda \boldQ^\T
    = \boldQ \sqrt{\boldLambda} \sqrt{ \boldLambda} \boldQ^\T
    = (\sqrt{ \boldLambda} \boldQ^\T)^\T \sqrt{ \boldLambda} \boldQ^\T
```

and applying the QR decomposition to $\sqrt{ \boldLambda} \boldQ^\T$. This allows us to write $\sqrt{ \boldLambda} \boldQ^\T = \tilde \boldQ \boldR$, where $\boldR$ is nonsingular and upper triangular, and $\tilde \boldQ$ has orthonormal columns. Because the columns of $\tilde \boldQ$ are orthonormal,

```{math}
:enumerated: false

\boldA 
    = (\tilde \boldQ \boldR)^\T \tilde \boldQ \boldR
    =  \boldR^\T \tilde \boldQ^\T \tilde \boldQ \boldR
    =  \boldR^\T \boldR
```

Our decomposition has the properties stated in [fact %s](#fa-cholesky).

(ss-mn)=
### Norms and Continuity
Given vector sequence $\{\boldx_n\}$ in $\RR^K$ and any point $\boldx \in \RR^K$, we say that $\{\boldx_n\}$ **converges** to $\boldx$ and write **$\boldx_n \to \boldx$** if, for any $\epsilon > 0$, there exists an $N \in \NN$ such that $\|\boldx_n - \boldx\| < \epsilon$ whenever $n \geq N$. Another way to say this is that the real-valued sequence $z_n := \|\boldx_n - \boldx\|$ converges to zero in $\RR$ as $n \to \infty$.

```{prf:property}

The following results hold:

(i) If $\boldx_n \to \boldx$ and $\boldy_n \to \boldy$, then $\boldx_n + \boldy_n \to \boldx +  \boldy$.

(ii) If $\boldx_n \to \boldx$ and $\alpha \in \RR$, then $\alpha \boldx_n \to \alpha \boldx$.

(iii) $\boldx_n \to \boldx$ if and only if $\bolda^\T \boldx_n \to \bolda^\T \boldx$ for all $\bolda \in \RR^K$.
```


For econometrics it is helpful if we extend the notion of convergence to matrices. We can do this in a parallel fashion by defining norms over matrices. The **matrix norm** of $N \times K$ matrix $\boldA$ is defined as

$$
\| \boldA \| :=
    \max \left\{ 
        \| \boldA \boldx \| \,: \, 
        \boldx \in \RR^K, \; \| \boldx \| = 1
        \right\}
$$ (eq-mn)

In this expression there are two different norms in play: The left-hand side is a matrix norm, while on the right we have the ordinary vector norm.

Given a sequence of $N \times K$ matrices $\boldA_n$ and an $N \times K$ matrix $\boldA$, we say that $\boldA_n$ **converges** to $\boldA$ if the matrix norm deviation $\|\boldA_n - \boldA\|$ converges to zero in $\RR$.

While the value in {eq}`eq-mn` is not particularly easy to solve for in general, the definition is entirely standard, and one can show that the matrix norm behaves like the vector norm in many ways. For example,

```{prf:property}

For any conformable matrices $\boldA$ and $\boldB$, the matrix norm satisfies

(i) $\| \boldA \| \geq 0$ and $\| \boldA \| = 0$ if and only if all entries of $\boldA$ are zero,

(ii) $\| \alpha \boldA \| = |\alpha| \| \boldA \|$ for any scalar $\alpha$,

(iii) $\| \boldA + \boldB \| \leq \| \boldA \| + \| \boldB \|$, and

(iv) $\| \boldA \boldB \| \leq \| \boldA \| \| \boldB \|$.
```


Compare these results with those of [fact %s](#fa-profno).

````{prf:property}
:label: fa-mnbd

For any $J \times K$ matrix $\boldA$ with elements $a_{jk}$, we have

```{math}
:enumerated: false

\| \boldA \| \leq \sqrt{JK} \max_{jk} |a_{jk}|
```
````


This bound is handy. For example, it tells us that if every element of $\boldA$ is close to zero then $\| \boldA \|$ must also be close to zero.

(ss-mnns)=
#### Neumann Series
Let's look at an important result in analysis that uses matrix norms. Starting in [chapter 7](#c-tsm), we will investigate dynamic systems such as $\boldx_{t+1} = \boldA \boldx_t + \boldb$, where $\boldx_t$ represents the values of some variables of interest and $\boldA$ and $\boldb$ form the parameters in the law of motion for $\boldx_t$. One question we might ask in this setting is whether or not there exists a "stationary\" vector $\boldx \in \RR^N$, in the sense that $\boldx_t = \boldx$ implies $\boldx_{t+1} = \boldx$. In other words, we seek an $\boldx \in \RR^N$ that solves the system of equations

$$
\boldx = \boldA \boldx + \boldb 
    \qquad (\boldA \text{ is } N \times N \text{ and } \boldb \text{ is } N
    \times 1)
$$ (eq-lpfde)

We can get some insight by reflecting on the scalar case $x = a x + b$. If $|a| < 1$, then this equation has the unique solution

```{math}
:enumerated: false

\bar x 
    = \frac{b}{1-a} 
    = b \sum_{k=0}^{\infty} a^k
```

The second equality follows from elementary results on geometric series.

It turns out that a similar result is true in $\RR^N$ if we replace the condition $|a| < 1$ with an analogous result for matrices based around the matrix norm. We begin with the **Neumann series lemma**, which states the following:

```{prf:theorem}
:label: t-nms

If $\boldA$ is square and $\| \boldA^j \| < 1$ for some $j \in \NN$, then $\boldI - \boldA$ is invertible, and

$$
(\boldI - \boldA)^{-1} = \sum_{i=0}^{\infty} \boldA^i
$$ (eq-nmse)
```


Here $\boldI$ is the identity. The equality in {eq}`eq-nmse` means that the matrix sum $\sum_{i=0}^t \boldA^i$ converges to the left-hand side in matrix norm as $t \to \infty$. The sum is called the **Neumann series** associated with $\boldA$. The condition in the theorem ensures that it converges. When the condition holds, {eq}`eq-lpfde` has the unique solution

```{math}
:enumerated: false

\bar \boldx 
    = (\boldI - \boldA)^{-1} \boldb 
    = \sum_{i=0}^{\infty} \boldA^i  \boldb
```

How to test the condition in [theorem %s](#t-nms)? The most commonly used sufficient condition involves the **spectral radius** of $\boldA$, which is defined as

$$
\rho(\boldA) := \max \setntn{ |\lambda| }
    { \lambda \text{ is an eigenvalue of } \boldA}
$$ (eq-srmn)

In this definition, $|\lambda|$ is the **modulus** of the possibly complex number $\lambda$.[^3]

```{prf:property}

If $\rho(\boldA) < 1$, then $\| \boldA^j \| < 1$ for some $j \in \NN$.
```


To understand why $\rho(\boldA) < 1$ is sufficient for the result in the Neumann series lemma, consider [theorem %s](#t-nms) in this light: The claim is that $\sum_{i=0}^{\infty} \boldA^i$ is the inverse of $\boldI - \boldA$, so $\sum_{i=0}^t \boldA^i (\boldI - \boldA)$ should be close to $\boldI$ for large $t$. Evidently

```{math}
:enumerated: false

\sum_{i=0}^t \boldA^i (\boldI - \boldA)
    =
    \sum_{i=0}^t \boldA^i - \sum_{i=0}^t \boldA^{i+1}
    =
    \boldI - \boldA^{t+1}
```

Hence, for the result to go through, we require $\boldA^t \to \boldzero$ as $t \to \infty$. The case where $\boldA$ is diagonalizable gives the clearest insight here, so let's suppose, as in {ref}`ss-simmat`, that $\boldA = \boldP \boldLambda \boldP^{-1}$ where $\boldLambda$ is a diagonal matrix containing the eigenvalues $\lambda_1, \ldots, \lambda_N$ of $\boldA$ on its principal diagonal. As discussed in [fact %s](#fa-powsm), we then have

```{math}
:enumerated: false

\boldA^t
    = \boldP
    \left(
    \begin{matrix}
        \lambda_1^t & 0 & \cdots & 0 \\
        0 & \lambda_2^t & \cdots & 0 \\
        \vdots & \vdots &  & \vdots \\
        0 & 0 & \cdots & \lambda_N^t \\
    \end{matrix}
    \right)
    \boldP^{-1}
```

If $\rho(\boldA) < 1$, then $|\lambda_n| < 1$ for all $n$, and hence $\lambda_n^t \to 0$ as $t \to \infty$. It follows that $\boldA^t \to \boldzero$ as required.

## Further Reading

Good treatments of matrix and linear algebra include {cite}`jan` and {cite}`axler2015linear`.

## Exercises

```{exercise}
:label: ex-iai

Prove that the inverse of a nonsingular matrix is always nonsingular.
```


```{exercise}
:label: ex-ipde

Prove the claim in [fact %s](#fa-ipde) that if $\boldA$ is positive definite, then $\boldA$ is nonsingular. If you can, prove it without invoking positivity of its eigenvalues.
```


```{exercise}
:label: ex-ipde0

Prove [fact %s](#fa-ipde0).
```


```{exercise}
:label: ex-ui

Let $\boldA$ be $N \times K$ and full column rank. Show that $\boldx, \boldz \in \RR^K$ and $\boldx \neq \boldz$ implies $\boldA \boldx \neq \boldA \boldz$.
```


```{exercise}
:label: ex-nonsigmat

Prove [theorem %s](#t-nonsigmat). Before doing so, you might like to review [fact %s](#fa-nsins) and [theorem %s](#t-lmaeq).
```


```{exercise}
:label: ex-ch03_linear_algebra-auto-1

Let $\boldA$ be square. Assuming existence of the inverse $\boldA^{-1}$, show that $(\boldA^\T)^{-1} = (\boldA^{-1})^\T$.
```


```{exercise}
:label: ex-qfcv

Show that if $\bolde_i$ and $\bolde_j$ are the $i$th and $j$th canonical basis vectors of $\RR^N$ respectively, and $\boldA$ is an $N \times N$ matrix, then $\bolde_i^\T \boldA \bolde_j = a_{ij}$, the $i,j$th element of $\boldA$.
```


```{exercise}
:label: ex-iosol

Suppose that $\boldA$ is $N \times K$, the equation $\boldA \boldx = \boldb$ has a solution, and $K > N$. Show that the same equation has an infinity of solutions.[^4]
```


````{exercise}
:label: ex-ch03_linear_algebra-auto-2

Let

```{math}
:enumerated: false

\boldA := 
        \left(
        \begin{array}{cc}
            1 & -1  \\
            -1 & 1  
        \end{array}
        \right),
        \qquad
        \boldB := 
        \left(
        \begin{array}{cc}
            1 & 2 \\
            2 & 1 
        \end{array}
        \right)
```

Show that

(i) $\boldA$ is nonnegative definite.

(ii) $\boldB$ is *not* positive definite.
````


````{exercise}
:label: ex-ch03_linear_algebra-auto-3

Let $\boldA_1, \ldots, \boldA_J$ be invertible matrices. Use induction and [fact %s](#fa-bpiv) to show that the product of these matrices is invertible, and, in particular, that

```{math}
:enumerated: false

(\boldA_1 \boldA_2 \cdots \boldA_J)^{-1} 
        = \boldA_J^{-1} \cdots \boldA_2^{-1} \boldA_1^{-1}
```
````


```{exercise}
:label: ex-ch03_linear_algebra-auto-4

Show that for any matrix $\boldA$, the matrix $\boldA^\T \boldA$ is well-defined (i.e., multiplication is possible), square, and nonnegative definite.
```


```{exercise}
:label: ex-ch03_linear_algebra-auto-5

Show that if $\boldA$ and $\boldB$ are positive definite and $\boldA + \boldB$ is well-defined, then $\boldA + \boldB$ is also positive definite.
```


```{exercise}
:label: ex-ictz

Let $\boldA$ be $N \times K$. Show that if $\boldA \boldx = \boldzero$ for all $K \times 1$ vectors $\boldx$, then $\boldA = \boldzero$ (i.e., every element of $\boldA$ is zero). Show as a corollary that if $\boldA$ and $\boldB$ are $N \times K$ and $\boldA \boldx = \boldB \boldx$ for all $K \times 1$ vectors $\boldx$, then $\boldA = \boldB$.
```


```{exercise}
:label: ex-aaid

Let $\boldI$ be the $N \times N$ identity matrix.

(i) Explain why $\boldI$ is full column rank.

(ii) Show that $\boldI$ is the inverse of itself.

(iii) Let $\boldA := \alpha \boldI$. Give a condition on $\alpha$ such that $\boldA$ is positive definite.
```


```{exercise}
:label: ex-kwg

Let $\boldX := \boldI - 2 \boldu \boldu^\T$, where $\boldu$ is an $N \times 1$ vector with $\|\boldu\| =1$. Show that $\boldX$ is symmetric and $\boldX \boldX = \boldI$.
```


```{exercise}
:label: ex-ch03_linear_algebra-auto-6

Recall the definition of similarity of matrices, as given in {ref}`ss-simmat`. Let's write $\boldA \sim \boldB$ if $\boldA$ is similar to $\boldB$. Show that $\sim$ is an **equivalence relation** on the set of $N \times N$ matrices. In particular, show that, for any $N \times N$ matrices $\boldA, \boldB$, and $\boldC$, we have (i) $\boldA \sim \boldA$, (ii) $\boldA \sim \boldB$ implies $\boldB \sim \boldA$ and (iii) $\boldA \sim \boldB$ and $\boldB \sim \boldC$ implies $\boldA \sim \boldC$.
```


```{exercise}
:label: ex-chei

Confirm the claim in [fact %s](#fa-chei).
```


```{exercise}
:label: ex-ch03_linear_algebra-auto-7

Show that $\boldA$ is nonsingular if and only if $0$ is *not* an eigenvalue for $\boldA$.
```


```{exercise}
:label: ex-onsiei

Show that the only nonsingular idempotent matrix is the identity matrix.
```


```{exercise}
:label: ex-ch03_linear_algebra-auto-8

Let $\boldone$ be $N \times 1$ and let $\boldP := \frac{1}{N} \boldone \boldone^\T$. Verify that $\boldP$ is idempotent.
```


````{exercise}
:label: ex-wood

Show that, for conformable and suitably invertible matrices $\boldA$, $\boldU$ and $\boldW$, we have

```{math}
:enumerated: false

\left(\boldA+\boldU \boldW \right)^{-1} =
         \boldA^{-1} - \boldA^{-1}\boldU \left(\boldI+\boldW\boldA^{-1}\boldU
         \right)^{-1} \boldW\boldA^{-1}
```
````


```{exercise}
:label: ex-othmbp

Let $\boldQ$ be an orthogonal matrix. Show that $\boldQ^{-1} = \boldQ^\T$ and $\det \boldQ \in \{-1, 1\}$ both hold.
```


```{exercise}
:label: ex-eigdef

Use [theorem %s](#t-dismat) to prove the following part of [fact %s](#fa-eigdef): A symmetric matrix $\boldA$ is positive definite if and only if its eigenvalues are all positive.
```


```{exercise}
:label: ex-orthiis

Show that if $\boldQ$ is an $N \times N$ orthogonal matrix, then $\boldQ$ is an isometry on $\RR^N$. That is, for any $\boldx, \boldy \in \RR^N$, we have $\| \boldQ \boldx - \boldQ \boldy \| = \| \boldx - \boldy \|$.
```


```{exercise}
:label: ex-ch03_linear_algebra-auto-9

Let $\boldP$ be square, symmetric and idempotent. Let $S := \colspace \boldP$. Show that $\boldP = \proj S$.
```


```{exercise}
:label: ex-lssolnn

Consider [theorem %s](#t-lssol). If $N = K$, what does $\hboldx$ reduce to? Interpret.
```


```{exercise}
:label: ex-dptpek

Let $S$ be a linear subspace of $\RR^N$ and let $\boldP = \proj S$. Show that $\trace \boldP = \dim S$ without using the idempotence connection in [fact %s](#fa-idtrra).
```


```{exercise}
:label: ex-fr

Show that when $N \times K$ matrix $\boldB$ is full column rank, the matrix $\boldB^\T \boldB$ is nonsingular.[^5]
```


```{exercise}
:label: ex-ch03_linear_algebra-auto-10

Let $\boldA$ be an $N \times N$ matrix.

(i) Show that if $\boldI - \boldA$ is idempotent, then $\boldA$ is idempotent.

(ii) Show that if $\boldA$ is both symmetric and idempotent, then the matrix $\boldI - 2 \boldA$ is orthogonal.
```


```{exercise}
:label: ex-projon

Prove [theorem %s](#t-projon). (This takes a bit of work, of course.)
```


```{exercise}
:label: ex-anox

Let $\boldP = \boldB (\boldB^\T \boldB)^{-1} \boldB^\T$, as in [theorem %s](#t-projon), and let $\boldM$ be the residual projection. Show that $\boldM \boldB = \boldzero$ using matrix algebra.
```


```{exercise}
:label: ex-qrlsq

Let $\boldA$ be an $N \times K$ matrix with linearly independent columns and QR factorization $\boldA = \boldQ \boldR$ (see {ref}`ss-gso`). Fix $\boldb \in \RR^N$. Show that $\hboldx$ defined in {eq}`eq-dlse` can also be written as $\hboldx = \boldR^{-1} \boldQ^\T \boldb$.
```


````{exercise}
:label: ex-pcsin

Let's prove the Cauchy–Schwarz inequality $| \boldx^\T \boldy | \leq  \| \boldx \| \| \boldy \|$ from [fact %s](#fa-profno) via the orthogonal projection theorem. Let $\boldy$ and $\boldx$ be nonzero vectors in $\RR^N$ (since if either equals zero then the inequality is trivial), and let $\Span\{\boldx\}$ be all vectors of the form $\alpha \boldx$ for $\alpha \in \RR$.

(i) Letting $\boldP$ be the orthogonal projection onto $\Span\{\boldx\}$, show that

```{math}
:enumerated: false

\boldP \boldy = \frac{\boldx^\T\boldy}{\boldx^\T\boldx} \boldx
```

(ii) Using this expression and any relevant properties of orthogonal projections (see [theorem %s](#t-opt2)), confirm the Cauchy–Schwarz inequality.
````


```{exercise}
:label: ex-lsiosol

Prove the claim made after [theorem %s](#t-lssol) regarding failure of uniqueness without full column rank. In particular, let $\boldA$ be $N \times K$ with linearly *dependent* columns, and let $\boldP\boldb$ be the closest point to $\boldb$ in $\colspace \boldA$. Prove that there are infinitely many $\boldx \in \RR^K$ such that $\boldP\boldb = \boldA \boldx$.
```


```{exercise}
:label: ex-eigid

Let $\boldA$ be symmetric and idempotent. Show that every eigenvalue of $\boldA$ is either $0$ or $1$.
```


```{exercise}
:label: ex-nsqrtsnd

Show that if $\boldA$ is positive definite, then there exists a symmetric matrix $\boldC$ such that $\boldC \boldA \boldC = \boldI$.[^6]
```


### Solutions to Selected Exercises

```{solution} ex-iai

Let $\boldA$ be a nonsingular matrix. Being nonsingular, $\boldA$ is invertible, with inverse $\boldA^{-1}$. By the definition of the inverse, we have $\boldA \boldA^{-1} = \boldA^{-1} \boldA = \boldI$, where $\boldI$ is the identity. This tells us directly that $\boldA$ is an inverse for $\boldA^{-1}$. Hence $\boldA^{-1}$ is invertible, which is equivalent to nonsingularity.
```


```{solution} ex-ipde

Let $\boldA$ be positive definite and consider the following: If $\boldA$ is singular, then there exists nonzero $\boldx$ with $\boldA \boldx = \boldzero$ (see [fact %s](#fa-exinv)). But then $\boldx^\T \boldA \boldx = 0$ for nonzero $\boldx$. Contradiction.
```


```{solution} ex-ipde0

If $\boldx = \bolde_n$, then $\boldx^\T\boldA\boldx = a_{nn}$. The claim follows.
```


```{solution} ex-ui

Let $\boldA$, $\boldx$ and $\boldz$ be as stated in the question. Suppose, to the contrary, that $\boldA \boldx = \boldA \boldz$. Then $\boldA (\boldx -  \boldz) = \boldzero$, which, by [fact %s](#fa-frie), implies $\boldx - \boldz = \boldzero$, or $\boldx = \boldz$. Contradiction.
```


```{solution} ex-nonsigmat

Let $\boldA$ be a nonsingular matrix and let $T$ be the linear map associated with $\boldA$ via $T\boldx = \boldA \boldx$. Since $A$ is nonsingular, $T$ is also, by definition, nonsingular and hence, by [fact %s](#fa-nsins), has a nonsingular inverse $T^{-1}$. Being nonsingular, $T^{-1}$ is necessarily linear, and hence, by [theorem %s](#t-lmaeq), there exists a matrix $\boldB$ such that $T^{-1} \boldx = \boldB \boldx$ for all $\boldx$. By the definition of the inverse, we have $\boldA \boldB \boldx = T ( T^{-1} (\boldx )) = \boldx = \boldI \boldx$. Since this holds for any $\boldx$, we have $\boldA \boldB = \boldI$ (see [ex. %s](#ex-ictz)). A similar argument shows that $\boldB \boldA = \boldI$.

Regarding the second claim, $\boldA^{-1} \boldb$ is a solution to $\boldA \boldx = \boldb$, since $\boldA \boldA^{-1} \boldb = \boldI \boldb = \boldb$. Uniqueness follows from [fact %s](#fa-exinv) (and, in particular, the fact that nonsingularity of $\boldA$ includes the implication that the map $\boldx \mapsto \boldA \boldx$ is one-to-one).
```


```{solution} ex-iosol

Since the columns of $\boldA$ consist of $K$ vectors in $\RR^N$, the fact that $K > N$ implies that not all of the columns of $\boldA$ are linearly independent. (Recall [theorem %s](#t-exth).) It follows that $\boldA \boldz = \boldzero$ for some nonzero $\boldz$ in $\RR^K$, and hence $\boldA \lambda \boldz = \boldzero$ for any scalar $\lambda$. Now suppose that some $\boldx$ solves $\boldA \boldx = \boldb$. Then, for any $\lambda \in \RR$, we have $\boldA \boldx + \boldA \lambda \boldz = \boldA (\boldx + \lambda \boldz) = \boldb$. This proves the claim.
```


```{solution} ex-aaid

The solutions are as follows: (1) $\boldI$ is full column rank because its columns are the canonical basis vectors, which are independent. (2) By definition, $\boldB$ is the inverse of $\boldA$ if $\boldB \boldA = \boldA \boldB = \boldI$. It follows immediately that $\boldI$ is the inverse of itself. (3) A sufficient condition is $\alpha > 0$. If this holds, then given $\boldx \neq \boldzero$, we have $\boldx^\T \alpha \boldI \boldx = \alpha \| \boldx \|^2 > 0$.
```


````{solution} ex-kwg

First, $\boldX$ is symmetric because

```{math}
:enumerated: false

\boldX^\T 
        = (\boldI - 2 \boldu \boldu^\T)^\T
        = \boldI - 2 (\boldu \boldu^\T)^\T
        = \boldI - 2 (\boldu^\T)^\T \boldu^\T
        = \boldI - 2 \boldu \boldu^\T
        = \boldX
```

Second, $\boldX \boldX = \boldI$ because

```{math}
:enumerated: false

\begin{aligned}
\boldX \boldX 
        = (\boldI - 2 \boldu \boldu^\T) (\boldI^\T - 2 \boldu \boldu^\T)
        & = \boldI \boldI - 2 \boldI 2 \boldu \boldu^\T
         + (2 \boldu \boldu^\T)(2 \boldu \boldu^\T)
         \\
        & = \boldI  - 4 \boldu \boldu^\T + 4 \boldu \boldu^\T\boldu \boldu^\T
        = \boldI  - 4 \boldu \boldu^\T + 4 \boldu  \boldu^\T
        = \boldI
\end{aligned}
```

The second last equality is due to the assumption that $\boldu^\T \boldu = \| \boldu \|^2 = 1$.
````


````{solution} ex-chei

Let $\boldA$ be $N \times N$ and let $\boldI$ be the $N \times N$ identity. We have

```{math}
:enumerated: false

\begin{aligned}
\det(\boldA - \lambda \boldI) = 0
        & \iff \boldA - \lambda \boldI \text{ is singular}
        \\
        & \iff \exists \, \boldx \neq \boldzero \st
            (\boldA - \lambda \boldI) \boldx = \boldzero
        \iff \exists \, \boldx \neq \boldzero \st \boldA \boldx = \lambda \boldx
\end{aligned}
```

In other words, $\lambda$ is an eigenvalue of $\boldA$.
````


```{solution} ex-onsiei

Suppose that $\boldA$ is both idempotent and nonsingular. From idempotence we have $\boldA \boldA = \boldA$. Premultiplying by $\boldA^{-1}$ gives $\boldA = \boldI$.
```


````{solution} ex-wood

The claim is true because

```{math}
:enumerated: false

\begin{aligned}
&\left(\boldA+\boldU\boldW  \right) 
            \left[ 
                \boldA^{-1} - \boldA^{-1}\boldU 
                \left(\boldI+\boldW\boldA^{-1}\boldU \right)^{-1} 
                \boldW\boldA^{-1} 
            \right] \\
        & \quad = \boldI + \boldU\boldW \boldA^{-1} -
                (\boldU+\boldU \boldW \boldA^{-1}\boldU)(\boldI +
                \boldW\boldA^{-1}\boldU)^{-1}\boldW\boldA^{-1} \\
        & \quad = \boldI + \boldU\boldW \boldA^{-1} - \boldU(\boldI+
                \boldW\boldA^{-1}\boldU)(\boldI +
                \boldW\boldA^{-1}\boldU)^{-1}\boldW\boldA^{-1} \\
        & \quad = \boldI + \boldU\boldW \boldA^{-1} 
                - \boldU\boldW \boldA^{-1} = \boldI
        
\end{aligned}
```
````


```{solution} ex-othmbp

Let $\boldQ$ be an orthogonal matrix with columns $\boldu_1, \ldots, \boldu_N$. By the definition of matrix multiplication, the $m, n$th element of $\boldQ^\T \boldQ$ is $\boldu_m^\T \boldu_n$, which is $1$ if $m=n$ and zero otherwise. Hence $\boldQ^\T \boldQ = \boldI$. It follows from [fact %s](#fa-lrinv) that $\boldQ^\T$ is the inverse of $\boldQ$.

To see that $\det \boldQ \in \{-1, 1\}$, apply the results of [fact %s](#fa-detfacts) and [fact %s](#fa-owtra2) to the equality $\boldQ^\T \boldQ = \boldI$ to obtain $\det(\boldQ)^2 = 1$. The claim follows.
```


````{solution} ex-eigdef

Suppose that $\boldA$ is symmetric with eigenvalues $\lambda_1, \ldots, \lambda_N$. By [theorem %s](#t-dismat) we can decompose it as $\boldA = \boldQ \boldLambda \boldQ^\T$, where $\boldLambda$ is the diagonal matrix formed from eigenvalues and $\boldQ$ is an orthogonal matrix. Fixing $\boldx \in \RR^N$ and letting $\boldy := \boldQ^\T \boldx$, we have

$$
\boldx^\T \boldA \boldx 
        = (\boldQ^\T \boldx)^\T \boldLambda (\boldQ^\T \boldx)
        = \boldy^\T \boldLambda \boldy
        = \lambda_1 y_1^2 + \cdots + \lambda_N y_N^2
$$ (eq-eigdefs)
Suppose that all eigenvalues are positive. Take $\boldx$ to be nonzero. The vector $\boldy$ must be nonzero (why?), and it follows from {eq}`eq-eigdefs` that $\boldx^\T \boldA \boldx > 0$. Hence $\boldA$ is positive definite as claimed.

Conversely, suppose that $\boldA$ is positive definite. Fix $n \leq N$ and set $\boldx = \boldQ \bolde_n$. Evidently $\boldx$ is nonzero (why?). Hence $\boldx^\T \boldA \boldx > 0$. Since $\boldQ^\T$ is the inverse of $\boldQ$, it follows that

```{math}
:enumerated: false

\lambda_n 
        = \bolde_n^\T \boldLambda \bolde_n
        = (\boldQ^\T \boldx)^\T \boldLambda \boldQ^\T \boldx
        = \boldx^\T \boldQ \boldLambda \boldQ^\T \boldx
        = \boldx^\T \boldA \boldx
        > 0
```

Since $n$ was arbitrary, all eigenvalues are positive.
````


````{solution} ex-orthiis

Fixing $\boldx, \boldy \in \RR^N$ and letting $\boldz := \boldx - \boldy$ we have

```{math}
:enumerated: false

\| \boldQ \boldx - \boldQ \boldy \|^2 
        = \| \boldQ \boldz \|^2 
        = (\boldQ \boldz)^\T \boldQ \boldz
        = \boldz^\T \boldQ^\T \boldQ \boldz
        = \boldz^\T \boldz
        = \| \boldz \|^2
        = \| \boldx - \boldy \|^2
        
```
````


````{solution} ex-lssolnn

If $N = K$, then, in view of the full column rank assumption and [theorem %s](#t-nonsigmat), the matrix $\boldA$ is nonsingular. By [fact %s](#fa-owtra2), $\boldA^\T$ is likewise nonsingular. Applying the usual rule for inverse of products ([fact %s](#fa-bpiv)), we have

```{math}
:enumerated: false

\hboldx
        = (\boldA^\T \boldA)^{-1} \boldA^\T \boldb
        = \boldA^{-1} (\boldA^\T)^{-1} \boldA^\T \boldb
        = \boldA^{-1} \boldb
```

This is the solution to the system $\boldA \boldx = \boldb$ when $\boldA$ is square and invertible.
````


````{solution} ex-dptpek

Let $S$ and $\boldP$ be as in the statement of the exercise and let $K := \dim S$. We aim to show that $\trace \boldP = K$. Let $\boldB$ be a matrix such that its columns form a basis for $S$. By the definition of dimension, $\boldB$ has $K$ columns. Applying {eq}`eq-projb`, we have

```{math}
:enumerated: false

\trace \boldP =\trace(\boldB (\boldB^\T\boldB)^{-1} \boldB^\T )
```

Recalling [fact %s](#fa-lmpt), we rearrange to get

```{math}
:enumerated: false

\trace[\boldB (\boldB^\T\boldB)^{-1} \boldB^\T]
        = \trace[(\boldB^\T\boldB)^{-1} \boldB^\T\boldB ]
        = \trace \boldI
```

where $\boldI$ is the $K \times K$ identity. The claim follows.
````


````{solution} ex-fr

Let $\boldA = \boldB^\T \boldB$. It suffices to show that $\boldA$ is positive definite, since this implies that its determinant is strictly positive, and any matrix with nonzero determinant is nonsingular. To see that $\boldA$ is positive definite, pick any $\boldb \neq \boldzero$. We must show that $\boldb^\T \boldA \boldb > 0$. To see this, observe that

```{math}
:enumerated: false

\boldb^\T \boldA \boldb 
        = \boldb^\T \boldB^\T \boldB \boldb 
        = (\boldB \boldb)^\T \boldB \boldb = \| \boldB \boldb \|^2
```

By the properties of norms, the last term is zero only when $\boldB \boldb = \boldzero$. But this is not true because $\boldb \neq \boldzero$ and $\boldB$ is full column rank (see [theorem %s](#t-fdlinin)).
````


````{solution} ex-projon

Let $S$ and $\boldP$ be as stated in [theorem %s](#t-projon) and let $\boldB$ be a matrix such that the columns of $\boldB$ form a basis of $S$. Fix $\boldy \in \RR^N$. The claim is that $\hboldy := \boldB (\boldB^\T\boldB)^{-1} \boldB^\T \boldy$ is the orthogonal projection of $\boldy$ onto $S$. To verify this, we need to show that

(i) $\hboldy \in S$ and

(ii) $\boldy - \hboldy \perp S$.

Part (i) is true because $\hboldy$ can be written as $\hboldy = \boldB \boldx$ where $\boldx := (\boldB^\T\boldB)^{-1} \boldB^\T \boldy$. The vector $\boldB \boldx$ is a linear combination of the columns of $\boldB$. Since these columns form a basis of $S$, they must lie in $S$. Hence $\hboldy \in S$ as claimed.

Regarding (ii), from the assumption that $\boldB$ gives a basis for $S$, all points in $S$ have the form $\boldB \boldx$ for some $\boldx \in \RR^K$. Thus (ii) translates to the claim that

```{math}
:enumerated: false

\boldy - \boldB (\boldB^\T \boldB)^{-1} \boldB^\T \boldy \perp \boldB
        \boldx
        \quad \text{for all } \boldx \in \RR^K
```

This is true because if $\boldx \in \RR^K$, then

```{math}
:enumerated: false

(\boldB \boldx)^\T [\boldy - \boldB (\boldB^\T \boldB)^{-1} \boldB^\T
        \boldy] 
        =
        \boldx^\T [\boldB^\T \boldy - \boldB^\T \boldB (\boldB^\T \boldB)^{-1}
        \boldB^\T \boldy] 
        =
        \boldx^\T [\boldB^\T \boldy - \boldB^\T \boldy] 
        = 0
        
```
````


````{solution} ex-anox

These are straightforward. For example,

```{math}
:enumerated: false

\boldM \boldB 
        = \boldB - \boldB (\boldB^\T \boldB)^{-1} \boldB^\T \boldB 
        = \boldzero
        
```
````


````{solution} ex-qrlsq

Let $\boldA$, $\boldQ$, $\boldR$ and $\boldb \in \RR^N$ be as in the statement of the exercise. The claim is that $\hboldx$ defined in {eq}`eq-dlse` is equal to $\tilde \boldx := \boldR^{-1} \boldQ^\T \boldb$. To show this, in view of linear independence of the columns of $\boldA$, it suffices to show that $\boldA \tilde \boldx = \boldA \hboldx$, or

```{math}
:enumerated: false

\boldA (\boldA^\T \boldA)^{-1} \boldA^\T \boldb 
        = \boldQ \boldR \boldR^{-1} \boldQ^\T \boldb
```

After simplifying, we see it suffices to show that $\boldA (\boldA^\T \boldA)^{-1} \boldA^\T = \boldQ \boldQ^\T$. Since $\boldA$ and $\boldQ$ have the same column space, this follows from [theorem %s](#t-projon).
````


```{solution} ex-pcsin

Regarding (i), the expression for $\boldP \boldy$ given in [exercise %s](#ex-pcsin) can also be written as $\boldx (\boldx^\T\boldx)^{-1} \boldx^\T\boldy$. Since $\boldx$ is a basis for $\Span\{\boldx\}$, the validity of this expression as the projection onto $\Span\{\boldx\}$ follows immediately from [theorem %s](#t-projon). Regarding (ii), recall that orthogonal projections contract norms, so that, in particular, $\| \boldP \boldy \| \leq \| \boldy \|$ must hold. Using our expression for $\boldP \boldy$ from (i) and rearranging gives the desired bound $|\boldx^\T \boldy | \leq  \| \boldx \| \| \boldy \|$.
```


```{solution} ex-lsiosol

By the definition of orthogonal projection, $\boldP\boldb \in \colspace \boldA$, and hence there exists a vector $\boldx$ such that $\boldP\boldb = \boldA \boldx$. Since $\boldA$ has linearly dependent columns, there exists a nonzero vector $\bolda$ such that $\boldA \bolda = \boldzero$. Hence $\boldA \lambda \bolda = \boldzero$ for all $\lambda \in \RR$. For each such $\lambda$ we have $\boldP\boldb = \boldA \boldx = \boldA \boldx + \boldA \lambda \bolda = \boldA (\boldx + \lambda \bolda)$. This proves the claim.
```


```{solution} ex-eigid

By the spectral decomposition theorem, we know that $\boldA = \boldQ \boldLambda \boldQ^\T$ where $\boldQ$ is an orthogonal matrix and $\boldLambda$ is the diagonal matrix formed from the eigenvalues of $\boldA$. It follows that $\boldA^2 = \boldQ \boldLambda^2 \boldQ^\T$ and, since $\boldA$ is idempotent, that $\boldQ \boldLambda \boldQ^\T =  \boldQ \boldLambda^2 \boldQ^\T$. From this we obtain $\boldLambda = \boldLambda^2$. For diagonal matrices, powers are obtained by taking powers of the diagonal elements, from which we get $\lambda_n = \lambda^2_n$ for any eigenvalue $\lambda_n$. Hence $\lambda_n \in \{0, 1\}$ as claimed.
```


[^1]: Formally, $K$ dimensional subspaces have measure zero in $\RR^N$ whenever $K < N$. Hence every probability measure that is absolutely continuous with respect to Lebesgue measure puts zero mass on such a set. Less formally, consider the case where $N = 3$ and $K=2$. Then $\colspace \boldA$ forms at most a 2-dimensional plane in $\RR^3$. Intuitively, this set has no volume in $\RR^3$ because planes have no "thickness." Similarly, while we might visualize $\boldb$ as a dot in $\RR^3$, as a point it is in fact infinitesimally small. Hence the chance of a randomly chosen $\boldb$ lying in $\colspace \boldA$ is zero.

[^2]: If we permit the eigenvectors to be complex, then the requirement is that they form a basis of $\CC^N$, the set of complex $N$-vectors. In this case there is a sense in which almost all matrices are diagonalizable.

[^3]: The modulus of $a + ib \in \CC$ is $(a^2 + b^2)^{1/2}$. If the imaginary part is zero, this reduces to the usual notion of absolute value.

[^4]: This is the so-called "underdetermined\" case, where the number of equations is less than the number of unknowns. Intuitively, we do not have enough restrictions to pin down values uniquely.

[^5]: Hint: In view of [fact %s](#fa-ipde), it suffices to show that $\boldB^\T \boldB$ is positive definite.

[^6]: Hint: Look at [fact %s](#fa-sqrtsnd) and the argument that follows it.
