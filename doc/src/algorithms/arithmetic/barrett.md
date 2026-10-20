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

## Lemire–Kaser–Kurz reduction {#lemire}

The method of Lemire, Kaser and Kurz [@LKK2019] computes the remainder
directly from the fractional part of $x/n$. With $c = \lceil B/n \rceil$,
the low word $xc \bmod B$ approximates $B \{x/n\}$, where $\{\cdot\}$
denotes the fractional part; multiplying it by $n$ and keeping the high
word gives $x \bmod n$. The approximation is exact when $n$ and $x$ are
below $2^{\bits/2}$, which is the condition required by `n_mod_lemire`.
The method needs no correction step, but its range restriction makes it
useful only for half-word values.

::: functions n_mod_barrett n_mod_barrett_lazy n_barrett_precomp n_mod_lemire n_lemire_precomp :::
