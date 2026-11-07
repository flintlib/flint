---
title: Barrett and Lemire reduction
summary: Reduce a single word modulo $n$ with two products and at most one correction, using $\lfloor 2^{\bits}/n \rfloor$ (Barrett) or $\lceil 2^{\bits}/n \rceil$ (Lemire, Kaser and Kurz).
complexity: Barrett: one high and one low product and one conditional subtraction. Lemire: one low and one high product, no correction.
order: 40
---

The methods of Barrett and of Lemire, Kaser and Kurz reduce a *single*
word $x$ modulo $n$, rather than a two-word product, using a precomputed
approximation of $B/n$. They are cheaper than a two-word
[Möller–Granlund](algorithms/arithmetic/moller-granlund) reduction when
the value to be reduced is known to fit in one word, as in hashing,
sieving and table lookups. Throughout, $B = 2^{\bits}$.

## Barrett reduction {#barrett}

Barrett's method [@Bar1987] estimates the quotient as follows.

::: proposition {#prop-barrett}
Let $2 \le n < B$, $m = \lfloor B/n \rfloor$ and $0 \le x < B$, and let
$\tilde q = \lfloor x m / B \rfloor$. Then
$\lfloor x/n \rfloor - 1 \le \tilde q \le \lfloor x/n \rfloor$, so
$0 \le x - \tilde q n < 2n$.
:::

::: proof
From $B/n - 1 < m \le B/n$ we get
$$
\frac{x}{n} - 1 < \frac{x}{n} - \frac{x}{B} < \frac{x m}{B} \le \frac{x}{n},
$$
using $x < B$. Taking floors gives the claim.
:::

The value $x - \tilde q n$ is at most $x$, so it is computed exactly with
word arithmetic even when $2n > B$. `n_mod_barrett_lazy` returns it, and
`n_mod_barrett` subtracts $n$ once if needed. `n_barrett_precomp` computes
$m$ with one division.

## Barrett reduction with an `nmod_t` {#nmod}

`NMOD_RED` reduces a word with the same kind of quotient estimate, but its
constant comes from the `nmod_t`, so it needs no division and no extra
field. Let $\nu$ be the number of bits of $n$, that is
$2^{\nu - 1} \le n < 2^{\nu}$ and $\mathit{norm} = \bits - \nu$. The
`ninv` field is the [Möller–Granlund reciprocal](algorithms/arithmetic/moller-granlund#reciprocal)
of $n 2^{\mathit{norm}}$.

::: lemma {#lem-nmod-m}
Let $m = 2^{\bits - 1} + \lfloor \mathit{ninv}/2 \rfloor$. Then
$m = \lfloor (2^{\nu + \bits - 1} - 1)/n \rfloor$, and
$2^{\bits - 1} \le m < B$.
:::

::: proof
By definition $B + \mathit{ninv} = \lfloor (B^2 - 1)/(n 2^{\mathit{norm}}) \rfloor
= \lfloor (2^{\nu + \bits} - 1)/n \rfloor$, using
$\lfloor \lfloor x \rfloor / k \rfloor = \lfloor x/k \rfloor$ for integers $k \ge 1$ and
$\lfloor (B^2 - 1)/2^{\mathit{norm}} \rfloor = 2^{\nu + \bits} - 1$. Hence
$m = \lfloor (B + \mathit{ninv})/2 \rfloor = \lfloor (2^{\nu + \bits} - 1)/(2n) \rfloor$.
With $y = 2^{\nu + \bits - 1}$, an integer $k$ satisfies $2nk \le 2y - 1$ if
and only if $2nk \le 2y - 2$, since $2nk$ is even; so this is
$\lfloor (y - 1)/n \rfloor$. The bounds follow from $\mathit{ninv} < B$.
:::

::: theorem {#thm-nmod-red}
Let $1 \le n < B$, $0 \le a < B$, $m$ as in [](#lem-nmod-m), and
$\tilde q = \lfloor a m / 2^{\nu + \bits - 1} \rfloor$. Then
$\lfloor a/n \rfloor - 1 \le \tilde q \le \lfloor a/n \rfloor$, so
$0 \le a - \tilde q n < 2n$.
:::

::: proof
Write $mn = 2^{\nu + \bits - 1} - 1 - \rho$ with $0 \le \rho < n$
([](#lem-nmod-m)). Then
$$
\frac{a m}{2^{\nu + \bits - 1}} = \frac{a}{n} - \frac{a (1 + \rho)}{n\, 2^{\nu + \bits - 1}},
$$
and the subtracted term is in $[0, 1)$ because $a < B$, $1 + \rho \le n$
and $\nu \ge 1$. Taking floors gives the claim.
:::

So `NMOD_RED` computes $\tilde q$ as the high word of $a m$ shifted right
by $\bits - 1 - \mathit{norm}$ bits, subtracts $\tilde q n$ (at most $a$, so
there is no overflow) and corrects once. This holds for every modulus
$1 \le n < B$, including powers of two and $n \ge 2^{\bits - 1}$, and every
word $a$. The cost is one high product, one low product and one comparison,
against a full product, a low product and two corrections for a
[Möller–Granlund](algorithms/arithmetic/moller-granlund) step on the
two-word number $(0, a)$.

When $n$ has `FLINT_BITS` bits (`mod.norm` $= 0$), $a < 2n$ for every word $a$,
so `NMOD_RED` skips the product and only corrects once. The test on `mod.norm`
is moved out of loops by the compiler, and the remaining loop of comparisons
can then be trivially vectorized.

Shifting $m$ gives the constants of the same family for other input sizes:
$\lfloor m / 2^k \rfloor = \lfloor (2^{\nu + \bits - 1 - k} - 1)/n \rfloor$
for $0 \le k < \nu + \bits - 1$. For instance $k = \nu - 1$ gives
$\lfloor (B - 1)/n \rfloor$, which is the constant of `n_barrett_precomp`
unless $n$ is a power of two.

## Two and three words {#two-words}

For $n < B/2$, the same estimate reduces a two-word $x < B^2$ with the
two-word constant $Q = \lfloor B^2/n \rfloor$: [](#prop-barrett) with
$B^2$ in place of $B$ shows that $\tilde q = \lfloor x Q / B^2 \rfloor$
satisfies $0 \le x - \tilde q n < 2n < B$. So the remainder is computed
from the low words of $x$ and $\tilde q n$ alone, and one subtraction of
$n$ finishes. `n_ll_mod_barrett` takes the low word of $\tilde q$, the
third word of $Q x$, from the low three words of this product
(`FLINT_MPN_MUL_3P2X2`, three word products).

A three-word $y = y_2 B^2 + y_1 B + y_0$ is first folded into two words,
$x = y_2 b_2 + y_1 b_1 + y_0 \equiv y \pmod n$ with $b_1 = B \bmod n$
and $b_2 = B^2 \bmod n$. Since $b_1, b_2 \le n - 1$,
$x \le (B - 1)(2n - 1) < B^2$ for $n < B/2$, whatever the words $y_i$
(`n_lll_mod_barrett`).

The three constants are the quotient and the remainders of the division of
$2^{\mathit{norm}} B^2$ by the normalized $d = 2^{\mathit{norm}} n$, in two
steps of [](algorithms/arithmetic/moller-granlund#alg-div) with the inverse
`ninv`: the first gives the high word $\lfloor B/n \rfloor$ of $Q$ and the
remainder $2^{\mathit{norm}} b_1$, the second the low word of $Q$ and
$2^{\mathit{norm}} b_2$ (`n_ll_barrett_precomp_preinv`).

## Lemire–Kaser–Kurz reduction {#lemire}

The method of Lemire, Kaser and Kurz [@LKK2019] computes the remainder
directly from the fractional part of $x/n$. With $c = \lceil B/n \rceil$,
the low word $xc \bmod B$ approximates $B \{x/n\}$, where $\{\cdot\}$
denotes the fractional part; multiplying it by $n$ and keeping the high
word gives $x \bmod n$. The approximation is exact when $n$ and $x$ are
below $2^{\bits/2}$, which is the condition required by `n_mod_lemire`.
The method needs no correction step, but its range restriction makes it
useful only for half-word values.

::: functions n_mod_barrett n_mod_barrett_lazy n_barrett_precomp n_ll_barrett_precomp_preinv n_ll_mod_barrett n_lll_mod_barrett n_mod_lemire n_lemire_precomp NMOD_RED :::
