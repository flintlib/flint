---
title: Squarefree factorisation
summary: Yun's algorithm in characteristic zero and the gcd iteration with $p$-th roots over finite fields, which split a polynomial into squarefree, pairwise coprime factors with distinct multiplicities using only gcds and exact divisions.
complexity: $O(M(n) \log n)$ field operations for a polynomial of degree $n$ with fast gcds (a constant number of gcds of total degree $O(n)$), plus $p$-th roots of at most $n$ coefficients in characteristic $p$.
order: 72
---

A polynomial $f$ over a field $K$ is *squarefree* if it is not divisible
by the square of a nonconstant polynomial. Every nonzero $f$ has a unique
*squarefree factorisation*
$$
f = c \prod_{i \ge 1} g_i^{\,i}
$$
with $c \in K^\times$ and monic, squarefree, pairwise coprime $g_i$, where
$g_i$ is the product of the irreducible factors of multiplicity exactly
$i$. It is the first step of every factorisation algorithm
([](algorithms/polynomials/factoring-finite-fields),
[](algorithms/polynomials/factoring-integers)) and of root finding, and
it costs only gcds and exact divisions. FLINT returns the nontrivial
$g_i$ together with their exponents $i$, made monic over a field and
primitive with positive leading coefficient over $\Z$.

The tool is the derivative. If $f = \prod_g g^{e_g}$ over the monic
irreducible $g$, then
$$
f' = \sum_g e_g\, g' \, \frac{f}{g}.
$$

::: lemma {#lem-gcd-derivative}
Suppose that every irreducible $g$ is separable ($\gcd(g, g') = 1$), which
holds in characteristic 0 and over finite fields. Then the multiplicity of
$g$ in $\gcd(f, f')$ is $e_g - 1$ if $e_g \ne 0$ in $K$, and $e_g$ if $e_g
= 0$ in $K$ (that is, if the characteristic $p$ divides $e_g$).
:::

::: proof
In the sum, every term except the one for $g$ is divisible by $g^{e_g}$.
The term for $g$ is divisible by exactly $g^{e_g - 1}$ if $e_g \ne 0$ in
$K$, since $g \nmid e_g g'$, and vanishes otherwise. Hence $g^{e_g - 1}
\| f'$ in the first case and $g^{e_g} \mid f'$ in the second.
:::

## Characteristic zero: Yun's algorithm {#yun}

In characteristic 0, [](#lem-gcd-derivative) gives $\gcd(f, f') =
\prod_i g_i^{\,i-1}$, so $f / \gcd(f, f') = \prod_i g_i$ is the
*squarefree part* of $f$. The multiplicities can be separated by
repeating this, but Yun's algorithm [@Yun1976] does it with gcds whose
degrees add up to $O(n)$.

::: algorithm Yun's algorithm {#alg-yun}
**Input:** $f$ of degree $n \ge 1$ over a field of characteristic 0,
monic (or primitive over a UFD).

**Output:** the pairs $(g_i, i)$ with $g_i \ne 1$.

1. $d \gets \gcd(f, f')$. If $d = 1$, return $(f, 1)$.
2. $v \gets f / d$, $w \gets f' / d$, $i \gets 1$.
3. Loop:
    1. $s \gets w - v'$. If $s = 0$, output $(v, i)$ if $v \ne 1$ and stop.
    2. $u \gets \gcd(v, s)$; output $(u, i)$ if $u \ne 1$.
    3. $v \gets v / u$, $w \gets s / u$, $i \gets i + 1$.
:::

::: theorem {#thm-yun}
[](#alg-yun) is correct. At the start of iteration $i$,
$$
v = \prod_{j \ge i} g_j, \qquad w = \sum_{j \ge i} (j - i + 1)\, g_j'\,
\frac{v}{g_j},
$$
and the gcd of iteration $i$ is $g_i$.
:::

::: proof
For $i = 1$, $d = \prod_j g_j^{\,j-1}$ gives $v = \prod_j g_j$ and $w =
f'/d = \sum_j j g_j' v / g_j$. Given the invariant,
$$
s = w - v' = \sum_{j > i} (j - i)\, g_j' \frac{v}{g_j}.
$$
Every term is divisible by $g_i$. For $j > i$ with $g_j \ne 1$, all terms
but the $j$-th are divisible by $g_j$, and the $j$-th is not, since $j - i
\ne 0$ and $\gcd(g_j, g_j') = 1$. Hence $\gcd(v, s) = g_i$, and $v / g_i$,
$s / g_i$ satisfy the invariant for $i + 1$. If $s = 0$, then $g_j = 1$
for all $j > i$, so $v = g_i$.
:::

The degrees of the polynomials in iteration $i$ are at most $\deg v \le
\sum_{j \ge i} \deg g_j$, and $\sum_i \sum_{j \ge i} \deg g_j = \sum_j j
\deg g_j = n$. With gcds in $O(M(m) \log m)$ for degree $m$
([](algorithms/polynomials/half-gcd)), the whole algorithm therefore
costs $O(M(n) \log n)$, the same as one gcd of degree $n$ up to a
constant factor [@Yun1976] [@GG2003, §14.6].

## Finite fields {#finite-fields}

In characteristic $p$ the derivative of $g^{e}$ vanishes when $p \mid e$,
and $f' = 0$ if and only if $f \in K[x^p]$. Over a perfect field, in
particular over $\mathbb{F}_q$ with $q = p^k$, such an $f$ is a $p$-th
power: if $f = \sum_i c_i x^{ip}$, then $f = (\sum_i c_i^{1/p} x^i)^p$,
where $c^{1/p} = c^{p^{k-1}}$ (and $c^{1/p} = c$ in a prime field). FLINT
uses the following algorithm, the classical gcd iteration [@Mus1971]
[@Knu1997, §4.6.2] combined with $p$-th roots [@GG2003, §14.6].

::: algorithm Squarefree factorisation over $\mathbb{F}_q$ {#alg-sqf-fq}
**Input:** $f$ of degree $n \ge 1$ over $\mathbb{F}_q$.

**Output:** the pairs $(g_i, i)$ with $g_i \ne 1$, the $g_i$ monic.

1. If $f' = 0$: compute $h$ with $h^p = f$ by taking $p$-th roots of the
   coefficients, factor $h$ recursively, multiply all exponents by $p$ and
   return.
2. $d \gets \gcd(f, f')$. If $d = 1$, return $(f, 1)$ (made monic).
3. $v \gets f / d$, $i \gets 1$. While $\deg v \ge 1$:
    1. $s \gets \gcd(v, d)$; $w \gets v / s$; output $(w, i)$ if $w \ne
       1$.
    2. $d \gets d / s$, $v \gets s$, $i \gets i + 1$.
4. If $d \ne 1$: $d$ is a $p$-th power; factor $d^{1/p}$ recursively and
   output its factors with exponents multiplied by $p$.
:::

::: theorem {#thm-sqf-fq}
[](#alg-sqf-fq) is correct. The exponents output in step 3 are not
divisible by $p$, those output in step 4 are, and all are distinct.
:::

::: proof
Write $f = \prod_g g^{e_g}$. By [](#lem-gcd-derivative), after step 2,
$$
d = \prod_{p \nmid e_g} g^{e_g - 1} \prod_{p \mid e_g} g^{e_g}, \qquad v =
\prod_{p \nmid e_g} g.
$$
By induction on $i$, at the start of iteration $i$ we have $v =
\prod_{p \nmid e_g,\, e_g \ge i} g$ and $d = \prod_{p \nmid e_g,\, e_g \ge
i} g^{e_g - i} \prod_{p \mid e_g} g^{e_g}$. Since $v$ is squarefree and
contains only factors with $p \nmid e_g$, $s = \gcd(v, d) = \prod_{p \nmid
e_g,\, e_g > i} g$, so $w = v/s$ is the product of the factors of
multiplicity exactly $i$ (which is 1 if $p \mid i$), and the updated $v$
and $d$ satisfy the invariant for $i + 1$. When $v = 1$, the remaining $d
= \prod_{p \mid e_g} g^{e_g}$ is a $p$-th power. In step 1 and step 4 the
recursion is applied to a polynomial of degree at most $n/p$, whose
squarefree factorisation, with exponents multiplied by $p$, gives the
factors of $f$ with multiplicities divisible by $p$.
:::

As in Yun's algorithm, the degrees of the gcds in step 3 add up to at
most $n$, so the cost is $O(M(n) \log n)$, plus at most $n/p$ $p$-th
roots of field elements per level of recursion.

The algorithm gives the factorisation with distinct exponents, not
necessarily with exponents $1, 2, 3, \ldots$ in order; FLINT returns the
factors in the order in which they are found.

## Squarefree part and squarefreeness test {#part}

The squarefree part $\prod_i g_i$ is $f / \gcd(f, f')$ (up to a unit) in
characteristic 0, while in characteristic $p$ it is computed from the
squarefree factorisation. To test whether $f$ is squarefree it suffices,
over a perfect field, to test $\gcd(f, f') = 1$: if $f' = 0$ and $\deg f
\ge 1$ then $f$ is a $p$-th power, hence not squarefree.

## Use in FLINT {#use}

`gr_poly_factor_squarefree` implements [](#alg-yun) for unique
factorisation domains of characteristic 0, after making $f$ monic (over a
field) or removing its content (otherwise), and
[](#alg-sqf-fq) for finite fields, where the $p$-th roots of the
coefficients are computed with `gr_fq_pth_root`.
`gr_poly_squarefree_part` and `gr_poly_is_squarefree` use the methods of
the previous section. The finite field modules wrap the generic code;
`fmpz_poly_factor_squarefree` has its own implementation of Yun's
algorithm with integer gcds, and the squarefreeness tests for `nmod_poly`
compute $\gcd(f, f')$ directly.

::: functions gr_poly_factor_squarefree gr_poly_squarefree_part gr_poly_is_squarefree nmod_poly_factor_squarefree nmod_poly_is_squarefree fmpz_mod_poly_factor_squarefree fmpz_mod_poly_is_squarefree fmpz_poly_factor_squarefree :::
