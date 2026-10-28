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

Shifting $m$ gives the constants of the same family for other input sizes:
$\lfloor m / 2^k \rfloor = \lfloor (2^{\nu + \bits - 1 - k} - 1)/n \rfloor$
for $0 \le k < \nu + \bits - 1$. For instance $k = \nu - 1$ gives
$\lfloor (B - 1)/n \rfloor$, which is the constant of `n_barrett_precomp`
unless $n$ is a power of two.

## Lemire–Kaser–Kurz reduction {#lemire}

The method of Lemire, Kaser and Kurz [@LKK2019] computes the remainder
directly from the fractional part of $x/n$. With $c = \lceil B/n \rceil$,
the low word $xc \bmod B$ approximates $B \{x/n\}$, where $\{\cdot\}$
denotes the fractional part; multiplying it by $n$ and keeping the high
word gives $x \bmod n$. The approximation is exact when $n$ and $x$ are
below $2^{\bits/2}$, which is the condition required by `n_mod_lemire`.
The method needs no correction step, but its range restriction makes it
useful only for half-word values.

::: functions n_mod_barrett n_mod_barrett_lazy n_barrett_precomp n_mod_lemire n_lemire_precomp NMOD_RED :::
