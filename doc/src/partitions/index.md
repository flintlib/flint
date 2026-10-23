---
title: partitions.h
header: partitions.h
description: The partition function p(n) for huge n by the Hardy–Ramanujan–Rademacher formula, with rigorous error bounds.
category: Number theory
order: 3
---

This module computes the integer partition function $p(n)$, the number of
ways of writing $n$ as a sum of positive integers without regard to
order: $p(4) = 5$ ($4, 3+1, 2+2, 2+1+1, 1+1+1+1$). Small values come from
a table and Euler's pentagonal recurrence; large values from the
Hardy–Ramanujan–Rademacher formula, a convergent series whose leading
term is about $e^{\pi \sqrt{2n/3}}/(4\sqrt 3\, n)$, so that $p(n)$ has
about $3.7 \sqrt n$ bits. The implementation [@Joh2012] evaluates a
truncation of the series in ball arithmetic with per-term precision,
factors the exponential sums symbolically, and adds Rademacher's bound
for the remainder, so the result is proved correct; $p(10^9)$, a number of
117 thousand bits, takes a few hundredths of a second.

The evaluation itself is done by `mp_real_partitions_hrr` in `mp_real`,
which uses hardware double-based balls (`dfloat`) for the many
low-precision terms and several threads for large $n$. For all values
$p(0), \ldots, p(n-1)$ at once, use `arith_number_of_partitions_vec`, and
for approximations `arb_partitions_fmpz`. See [](partitions/p) for usage
and [](algorithms/number-theory/partitions) for the mathematics.

::: summary algorithms/number-theory/partitions :::

## The partition function {#p}

::: function
void partitions_fmpz_fmpz(fmpz_t p, const fmpz_t n, int use_doubles)
void partitions_fmpz_ui(fmpz_t p, ulong n)

Set $p$ to $p(n)$, which is 0 for $n < 0$. The argument `use_doubles` is
ignored. The result is exact: the formula is evaluated in ball arithmetic
with a proved bound for the truncation error, and the function throws if
the final ball does not determine a unique integer (which is not
expected to happen). See [](partitions/p).

::: details
Values with $n < 128$ are read from a table, values with $n < 417$
(where $p(n) < 2^{64}$, on 64-bit machines) are computed with Euler's
pentagonal recurrence in single words, and larger values with
`mp_real_partitions_hrr`. The number of threads set with
`flint_set_num_threads` is used from about $n = 10^8$.
`partitions_fmpz_fmpz` accepts $n < 2^{128}$ and throws for larger $n$;
the result must also fit in an `fmpz`, which limits $n$ to about $1.4
\cdot 10^{21}$ on 64-bit machines. Beyond, use `mp_real_partitions_hrr`
directly.
:::
:::

::: macro
partitions_fmpz_ui_using_doubles(p, n)

Deprecated alias of `partitions_fmpz_ui`.
:::

## Bounds and approximations {#bounds}

::: function
void partitions_rademacher_bound(arf_t b, const fmpz_t n, ulong N)

Sets $b$ to an upper bound for

$$
M(n,N) = \frac{44 \pi^2}{225 \sqrt 3} N^{-1/2}
+ \frac{\pi \sqrt{2}}{75} \left( \frac{N}{n-1} \right)^{1/2}
\sinh\left(\frac{\pi}{N} \sqrt{\frac{2n}{3}}\right),
$$

Rademacher's bound for the error of the Hardy–Ramanujan–Rademacher
series truncated after the term $k = N$ [@Rademacher1937]. Requires $n
\ge 2$ and $N \ge 1$. The constants are rounded upwards and the bound is
evaluated with directed rounding at about 30 bits.
:::

::: function
void partitions_leading_fmpz(arb_t res, const fmpz_t n, slong prec)

Sets `res` to the dominant part of the leading term $k = 1$ of the
Hardy–Ramanujan–Rademacher series for $p(n)$,

$$
\frac{\sqrt{12}}{24n - 1} \left(1 - \frac{1}{t}\right) e^{t}, \qquad t =
\frac{\pi \sqrt{24n - 1}}{6},
$$

which omits the part $e^{-t}(1 + 1/t) \sqrt{12}/(24n-1)$ of the term
and approximates $p(n)$ with relative error about $e^{-t/2}$. The
precision of the exponential is increased by half the bit length of $n$
for the argument reduction. Requires $n \ge 1$.
:::
