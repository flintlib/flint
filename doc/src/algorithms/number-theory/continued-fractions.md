---
title: Continued fractions of rational numbers
summary: Computes the continued fraction of a rational number with a half-GCD algorithm that finds the terms from leading bits, evaluates continued fractions by balanced matrix products, and uses the same machinery for the simplest fraction in an interval, Farey neighbours and Dedekind sums.
complexity: $O(M(b) \log b)$ for the expansion of a fraction with $b$-bit numerator and denominator, and for its evaluation from the terms.
order: 72
---

The continued fraction of a rational number is the sequence of quotients
of the Euclidean algorithm. FLINT computes it in quasi-linear time by the
half-GCD method [@Sch1971]. For the basic theory see [@HarWri2008,
Chapter X] and [@GKP1994, §6.7].

## Continued fractions and matrices {#matrices}

For $x = p/q$ with $q > 0$, the Euclidean algorithm $r_{-1} = p$, $r_0 =
q$, $r_{i-1} = c_i r_i + r_{i+1}$ with $c_i = \lfloor r_{i-1}/r_i \rfloor$
gives the expansion $x = [c_0; c_1, \ldots, c_{k-1}]$, with $c_0 =
\lfloor x \rfloor$ and $c_i \ge 1$ for $i \ge 1$, and $c_{k-1} \ge 2$ if
$k \ge 2$ (the shorter of the two expansions). Writing

$$
E_c = \begin{pmatrix} c & 1 \\ 1 & 0 \end{pmatrix}, \qquad
E_{c_0} E_{c_1} \cdots E_{c_{i-1}} =
\begin{pmatrix} p_{i-1} & p_{i-2} \\ q_{i-1} & q_{i-2} \end{pmatrix},
$$

the fractions $p_i/q_i$ are the *convergents* of $x$, and for every $y >
0$

$$
[c_0; c_1, \ldots, c_{i-1}, y] = \frac{p_{i-1} y + p_{i-2}}{q_{i-1} y +
q_{i-2}}.
$$

So the first $i$ terms of $x$ are a matrix $M$ of determinant $(-1)^i$
with nonnegative entries (apart from the first row when $c_0 < 0$), and
the rest of the expansion is that of $y = M^{-1}(x)$, where $M$ acts as a
Möbius transformation. The entries of $M$ have about as many bits in
total as the terms consumed, and $y$ has correspondingly fewer bits than
$x$.

## Terms from leading bits {#lehmer}

The key observation, due to Lehmer [@Leh1938], is that the first terms of
$x$ depend only on its leading bits.

::: lemma {#lem-interval}
Let $I = [l, r]$ be an interval with $l > 1$, and let $c = \lfloor l
\rfloor = \lfloor r \rfloor$. Then every $x \in I$ has first term $c$,
and $E_c^{-1}(I) = [1/(r - c), 1/(l - c)]$ is again an interval.
More generally, the terms common to the expansions of $l$ and $r$ (as
long as both remainders stay above 1) are terms of every $x \in I$.
:::

::: proof
For $x \in I$, $c \le l \le x \le r < c + 1$, so $\lfloor x \rfloor = c$,
and $x \mapsto 1/(x - c)$ is decreasing on $I$. The general statement
follows by induction on the number of terms.
:::

If $x = a/b$ with $a > b > 0$ of $n$ bits and $a', b'$ are the leading $k$
bits of $a$ and $b$ (so $a = 2^{n-k} a' + \alpha$, $0 \le \alpha <
2^{n-k}$, and similarly for $b$), then $x$ lies in the interval $(a'/(b' +
1), (a'+1)/b')$, whose endpoints have $O(k)$ bits. Expanding both
endpoints in parallel until they disagree gives about $k/2$ bits worth of
terms of $x$ at the cost of a computation on $k$-bit numbers. Applying
the resulting matrix $M^{-1}$ to the full $(a, b)$ then reduces it to a
pair of about $n - k/2$ bits.

## The half-GCD algorithm {#hgcd}

::: algorithm Half-GCD for continued fractions {#alg-hgcd}
**Input:** $a > b > 0$ of $n$ bits.

**Output:** a matrix $M$, a product of elementary matrices $E_{c_i}$, and
$(a', b') = M^{-1}(a, b)$ of about $n/2$ bits, such that the terms $c_i$
are terms of $a/b$ (and of every number in $(a/(b+1), (a+1)/b)$).

1. If $n$ is small, run Lehmer's algorithm with word-size matrices and
   return.
2. Recursively apply the algorithm to the top $n/2$ bits of $(a, b)$,
   giving $M_1$ with entries of about $n/4$ bits; set $(a, b) \gets
   M_1^{-1}(a, b)$, now of about $3n/4$ bits.
3. Perform one division step, so that the next quotient is not split.
4. Recursively apply the algorithm to the top $n/2$ bits of the new
   pair, giving $M_2$; set $(a, b) \gets M_2^{-1}(a, b)$.
5. Return $M = M_1 E_c M_2$ and $(a, b)$.
:::

Each level costs a constant number of multiplications of $n$-bit numbers,
so the cost is $O(M(n) \log n)$, against $O(n^2)$ for the Euclidean
algorithm. The half-GCD algorithm for integer GCDs, which this one
parallels, is described on [Integer GCD](algorithms/integers/gcd#hgcd). As with all half-GCD algorithms, the delicate part is the
correction when the leading bits produce a term that is not a term of the
full number; FLINT avoids it by working with intervals throughout
([](#lem-interval)), so that every term written is guaranteed.

## Evaluation {#evaluation}

Evaluating $[c_0; \ldots, c_{n-1}]$ amounts to computing the matrix
product $E_{c_0} \cdots E_{c_{n-1}}$. Multiplying the factors in order
costs $O(nb)$ for a $b$-bit result, since the partial products grow
linearly. Splitting the product in halves recursively balances the sizes
of the operands, and the cost becomes $O(M(b) \log n)$, as in
[](algorithms/integers/binary-splitting).

## The simplest fraction in an interval {#simplest}

The *Stern–Brocot tree* contains every positive rational once: its root is
$1/1$, and the children of a node are obtained as mediants
$(a+c)/(b+d)$ of the node and its nearest ancestors $a/b < c/d$
[@GKP1994, §4.5]. A path in the tree is a continued fraction expansion,
read as run lengths of left and right moves.

::: proposition {#prop-simplest}
Let $0 < l \le r$, and let $p/q$ be the first node on the path of the
Stern–Brocot tree from the root towards $l$ that lies in $[l, r]$. Then
every fraction $p'/q' \in [l, r]$ satisfies $p' \ge p$ and $q' \ge q$. In
particular $p/q$ is the simplest fraction in $[l, r]$.
:::

::: proof
The search starts with the bounds $0/1$ and $1/0$ and keeps two bounds
$a/b < c/d$ with $bc - ad = 1$, replacing one of them by the mediant
$(a+c)/(b+d)$ as long as the mediant is outside $[l, r]$; so $[l, r]$
stays inside the open interval $(a/b, c/d)$. Since $bc - ad = 1$, every
fraction strictly between $a/b$ and $c/d$ can be written $(\lambda a +
\mu c)/(\lambda b + \mu d)$ with integers $\lambda, \mu \ge 1$, so its
numerator is at least $a + c$ and its denominator at least $b + d$. When
the mediant $p/q = (a+c)/(b+d)$ lies in $[l, r]$, it attains both
bounds.
:::

In terms of continued fractions, the common terms of $l$ and $r$ are the
common part of the two paths. FLINT computes them with the interval
version of the half-GCD algorithm, which handles a ball $[l, r]$ directly
(`_fmpq_ball_get_cfrac`), and then completes the expansion with the
smallest integer in the transformed interval. Negative and mixed-sign
intervals are reduced to this case by subtracting an integer or noting
that $0$ or the smallest integer in the interval is the answer.

## Farey neighbours {#farey}

Two fractions $a/b < c/d$ with $bc - ad = 1$ are adjacent in every Farey
sequence of order $Q$ with $\max(b, d) \le Q < b + d$. Given $x = p/q$
and $Q \ge q$, the left neighbour $p_l/q_l$ is therefore determined by $p
q_l \equiv 1 \pmod q$ and $Q - q < q_l \le Q$, which `fmpq_farey_neighbors`
solves with one modular inverse; the right neighbour follows from the
mediant relation $p_l/q_l + p_r/q_r = v \cdot p/q$ with $v = \lfloor (Q +
q_l)/q \rfloor$ [@HarWri2008, Chapter III].

## Dedekind sums {#dedekind}

By the reciprocity law of Dedekind sums [@RadGro1972], $s(h, k)$ is a
simple function of the continued fraction of $k/h$: with the partial
quotients $c_1, \ldots, c_n$ and the convergent matrix $M = (m_{ij})$ of
$k/(h \bmod k)$,

$$
12\, s(h, k) = \sum_{i=1}^{n} (-1)^{i+1} c_i + \frac{m_{21} - (-1)^{n}
m_{12}}{m_{11}} - 3 \cdot [n \text{ odd}]
$$

[@Hickerson1977]. FLINT evaluates it with word-size quotients when $k$
fits in a word and otherwise with the half-GCD code, which accumulates the
alternating sum of the terms instead of storing them. The cost is that of
a GCD.

## Use in FLINT {#use}

- `fmpq_get_cfrac` computes the terms with `_fmpq_ball_get_cfrac`, the
  interval half-GCD (with `_fmpq_hgcd` as its recursive kernel and
  Lehmer's method on one and two words below), then finishes with
  division steps; `fmpq_get_cfrac_naive` is the Euclidean algorithm.
- `fmpq_set_cfrac` multiplies the elementary matrices in a balanced tree,
  with plain products below 32 terms.
- `fmpq_simplest_between`, `fmpq_dedekind_sum` and the large-modulus path
  of rational reconstruction (`_fmpq_reconstruct_fmpz_2`, see
  [](algorithms/number-theory/rational-reconstruction)) use the same
  half-GCD code.
- `fmpq_farey_neighbors` uses the modular inverse formula above.

::: functions fmpq_get_cfrac fmpq_set_cfrac fmpq_simplest_between fmpq_farey_neighbors fmpq_dedekind_sum :::
