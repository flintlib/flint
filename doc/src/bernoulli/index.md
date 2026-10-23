---
title: bernoulli.h
header: bernoulli.h
description: Exact Bernoulli numbers — isolated values by the zeta function and Harvey's multimodular algorithm, ranges by a reverse zeta iterator, residues modulo primes, and a global cache.
category: Number theory
order: 2
---

This module computes the Bernoulli numbers, defined by the exponential
generating function

$$
\frac{x}{e^x - 1} = \sum_{n=0}^{\infty} B_n \frac{x^n}{n!},
$$

as exact fractions: $B_0 = 1$, $B_1 = -1/2$, $B_2 = 1/6$, and $B_n = 0$
for odd $n \ge 3$. It provides isolated values (from the Riemann zeta
function, or for very large $n$ combined with Harvey's multimodular
algorithm), $B_n$ modulo a prime, ranges of consecutive values generated
in reverse order, and a cache used by the special functions of `arb` and
`acb` (gamma function, zeta function, Euler–Maclaurin summation) that
need many Bernoulli numbers repeatedly. The `arith` functions
`arith_bernoulli_number` and `arith_bernoulli_number_vec` are wrappers
and alternatives.

All algorithms are described in [](algorithms/number-theory/bernoulli),
and [](bernoulli/isolated) explains the choice of method for one value.

::: summary algorithms/number-theory/bernoulli :::

## Isolated Bernoulli numbers {#isolated}

::: function
void bernoulli_fmpq_ui(fmpq_t b, ulong n)
void _bernoulli_fmpq_ui(fmpz_t num, fmpz_t den, ulong n)

Set the result to $B_n$ in lowest terms with a positive denominator. The
value is copied from the cache if it is there; the cache is not
extended. Otherwise odd $n$ and $n < 18000$ use
`_bernoulli_fmpq_ui_zeta`, and larger even $n$ use
`_bernoulli_fmpq_ui_multi_mod` with the default tuning.
:::

::: function
void _bernoulli_fmpq_ui_zeta(fmpz_t num, fmpz_t den, ulong n)
void _bernoulli_fmpq_ui_multi_mod(fmpz_t num, fmpz_t den, ulong n, double alpha)

Set `num` and `den` to the numerator and denominator of $B_n$ in lowest
terms, by two methods.

::: details
The `zeta` version computes the denominator $d$ with
`arith_bernoulli_number_denom` (von Staudt–Clausen), approximates $B_n$
with `arb_bernoulli_ui_zeta` from $|B_n| = 2\, n!\, \zeta(n)/(2\pi)^n$ at
a precision slightly above the size of $d B_n$, and rounds $d B_n$ to
the unique integer in the resulting ball, increasing the precision if
needed.

The `multi_mod` version computes the numerator modulo the product $M$ of
the primes $5 \le p$ with $(p - 1) \nmid n$ (residues by
`bernoulli_mod_p_harvey`, in parallel with the threads of
`flint_set_num_threads`) up to a fraction `alpha` of the bits of the
numerator, and the remaining high bits from the zeta function: the
numerator is $r + tM$ where $r$ is the residue and $t$ the unique integer
in a ball for $(d B_n - r)/M$. The parameter $0 \le \alpha \le 1$
balances the two parts; a negative value selects the default, 0 for $n <
18000$, $0.005 + 3.6 \cdot 10^{-6} n$ for $n < 60000$ and $\min(0.18 +
0.5 \cdot 10^{-6} n, 0.28)$ above. For $n < 10$ or odd $n$ the zeta
version is called.
:::
:::

::: function
ulong bernoulli_mod_p_harvey(ulong n, ulong p)

Returns $B_n \bmod p$ for a prime $p < 2^{\mathtt{FLINT\_D\_BITS}}$
($2^{53}$ on 64-bit machines), computed with Harvey's algorithm
[@Har2010] in time $O(p \log^{O(1)} p)$, independent of $n$. If $p$
divides the denominator of $B_n$ (that is, $n \ge 2$ even with $(p - 1)
\mid n$, or $n = 1$ and $p = 2$), returns `UWORD_MAX`.

::: details
A return value of 0 means that $p$ divides the numerator of $B_n$; for
$n$ even, $p \ge 5$ and $(p-1) \nmid n$ this happens if and only if $p$
is irregular with an index $\equiv n \pmod{p-1}$, or $p \mid n$. The
function reduces $n$ modulo $p - 1$ by Kummer's congruence and computes
$B_m/m \bmod p$ for $m = n \bmod (p-1)$.
:::
:::

## Ranges of Bernoulli numbers {#ranges}

::: type
bernoulli_rev_t

An iterator that generates the even-indexed Bernoulli numbers $B_n,
B_{n-2}, B_{n-4}, \ldots, B_0$ in decreasing order, without caching. An
array of length one of `bernoulli_rev_struct`.

::: details
The numbers are computed from the zeta series $\zeta(k) = \sum_{j \ge
1} j^{-k}$ summed directly in fixed-point arithmetic, with a table of the
odd powers $j^{-n}$ that is updated by multiplication by $j^2$ when $n$
decreases by 2 [@Blo2009]. The reverse order is what makes the updates
multiplications rather than divisions, and lets the precision and the
number of terms decrease as $n$ decreases. See
[](algorithms/number-theory/bernoulli#zeta).
:::
:::

::: function
void bernoulli_rev_init(bernoulli_rev_t iter, ulong n)

Initialises `iter` so that the first call to `bernoulli_rev_next` returns
$B_n$. An odd $n$ is replaced by $n - 1$.
:::

::: function
void bernoulli_rev_next(fmpz_t numer, fmpz_t denom, bernoulli_rev_t iter)

Sets `numer` and `denom` to the numerator and denominator of $B_k$ in
lowest terms, where $B_k$ is the next number of the sequence, and
advances `iter` to $B_{k-2}$. Must not be called after $B_0$ has been
returned.
:::

::: function
void bernoulli_rev_clear(bernoulli_rev_t iter)

Frees the memory used by `iter`.
:::

::: function
void bernoulli_fmpq_vec_no_cache(fmpq * res, ulong a, slong num)

Sets `res` to the `num` consecutive Bernoulli numbers $B_a, B_{a+1},
\ldots, B_{a + \mathit{num} - 1}$, without reading or writing the cache.
To reuse a vector, use `bernoulli_cache_compute` and `bernoulli_fmpq_ui`
instead. Throws if $a > 2^{31}$ or $\mathit{num} > 10^9$.

::: details
For $a = 0$ and $\mathit{num} \le 128$ the function calls
`arith_bernoulli_number_vec`. Otherwise it runs reverse iterators over
blocks of indices; for $\mathit{num} \ge 200$ the blocks (of size about
$(a + \mathit{num})/32$, at least 128) are distributed over the threads
of `flint_set_num_threads`. The function is not optimised for small
`num` at large $a$, where each block costs about as much as $B_{a +
\mathit{num}}$ alone.
:::
:::

## Cache {#cache}

The cache holds $B_0, \ldots, B_{N-1}$ as an array of `fmpq`. It is
thread-local if FLINT is built with thread-local storage, so that each
thread has its own copy, and it is freed by `flint_cleanup`.

::: variable
slong bernoulli_cache_num
fmpq * bernoulli_cache

The number $N$ of cached values and the array of cached values $B_0,
\ldots, B_{N-1}$.
:::

::: function
void bernoulli_cache_compute(slong n)

Makes sure that $B_0, \ldots, B_{n-1}$ are cached. The cache grows by at
least 32 entries (below 128) or 128 entries at a time; new entries are
computed with `arith_bernoulli_number_vec` up to 128 and with
`bernoulli_fmpq_vec_no_cache` beyond.
:::

::: macro
BERNOULLI_ENSURE_CACHED(n)

Makes sure that $B_n$ is cached, by calling `bernoulli_cache_compute(n +
1)` if needed.
:::

## Bounds {#bounds}

::: function
slong bernoulli_bound_2exp_si(ulong n)

Returns an integer $b$ with $|B_n| \le 2^b$, using a table for $n <
512$ and otherwise the inequality $|B_n| < 4\, n!/(2\pi)^n < 4
(n+1)^{n+1} e^{-n}/(2\pi)^n$ evaluated in integer arithmetic with a
table of logarithms. Returns `WORD_MIN` if $B_n = 0$. The bound exceeds
the true exponent by at most one percent. Throws if $n > 2^{\bits - 6}$
or the computation would overflow; the function is meant for $n$ small
enough that $B_n$ could be computed exactly.
:::

## Internal functions {#internal}

These inline helpers estimate the precision and number of terms for the
zeta series used by the reverse iterator.

::: function
slong bernoulli_denom_size(slong n)

Returns a crude estimate of the number of bits of the denominator of
$B_n$: by von Staudt–Clausen the denominator is at most $n^{\sigma_0(n)}$,
and the estimate is $\tfrac12 \log_2(n)\, n^{1.066/\log\log n}$, using
$\sigma_0(n) < n^{1.066/\log\log n}$ for large $n$.
:::

::: function
slong bernoulli_global_prec(ulong nmax)

Returns the working precision of the iterator for $B_{\mathit{nmax}}$:
`arith_bernoulli_number_size(nmax) + bernoulli_denom_size(nmax)`. No
guard bits are added, as the denominator estimate is pessimistic.
:::

::: function
slong bernoulli_zeta_terms(ulong s, slong prec)

Returns the number $N$ of terms of the zeta series $\sum_{j \le N}
j^{-s}$ needed for an error below $2^{-\mathit{prec}}$: about
$2^{(\mathit{prec} + 1)/(s - 1)}$, made odd.
:::

::: function
slong bernoulli_power_prec(slong i, ulong s1, slong wp)

Returns the precision for the term $i^{-s_1}$ in a sum computed to
absolute precision $\mathit{wp}$: $\mathit{wp} - s_1 \log_2 i$, at least
10.
:::

::: macro
BERNOULLI_REV_MIN

The reverse iterator uses a table for $n <$ `BERNOULLI_REV_MIN` (32).
:::
