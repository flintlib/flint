---
title: Binary splitting of hypergeometric series in mp_real
summary: The general hypergeometric binary splitting of mp_real: exact leaves, removal of common factors, truncation near the root, rigorous tail bounds, transform reuse between the products of a merge and precomputed transforms of power tables, and parallelism.
complexity: $O(M(N)\log^2 N)$ for $N$ bits of a geometrically convergent series with polynomial coefficients.
order: 60
---

Binary splitting evaluates a sum of $N$ terms of a hypergeometric
series exactly as a ratio of two integers, organised as a balanced tree
so that the products at each level have similar sizes
([](algorithms/integers/binary-splitting)) [@Bre1976] [@HP1998]. It is
the method of choice for constants such as $\pi$, $e$, $\log 2$ and
$\zeta(3)$, and for the Taylor series of elementary functions at a
rational point (the bit-burst algorithm).

## The recursion {#recursion}

`mp_real_hypgeom_series` sums
$$
S = \frac{c_Q + c_P\sum_{k \ge 1} \frac{P(k)}{Q(k)}\prod_{j=1}^{k-1}\frac{R(j)}{Q(j)}}{c_D}
$$
(or its reciprocal) for integer polynomials $P$, $Q$, $R$, in the
format of y-cruncher's formula files [@Yee2025]. Over a range
$[a, b)$ of terms let
$$
T(a, b) = \sum_{k=a}^{b-1} P(k)\prod_{j=a}^{k-1} R(j)\prod_{j=k+1}^{b-1} Q(j),
\qquad Q(a, b) = \prod_{k=a}^{b-1} Q(k), \qquad R(a, b) = \prod_{k=a}^{b-1} R(k),
$$
so that the partial sum over $[a, b)$, relative to the product of the
ratios before $a$, is $T(a, b)/Q(a, b)$.

::: proposition {#prop-merge}
For $a < m < b$, with the left values $T_1, Q_1, R_1$ over $[a, m)$ and
the right values $T_2, Q_2, R_2$ over $[m, b)$,
$$
T(a, b) = T_1 Q_2 + R_1 T_2, \qquad Q(a, b) = Q_1 Q_2, \qquad R(a, b) = R_1 R_2.
$$
:::

::: proof
The terms with $k < m$ of $T(a, b)$ are those of $T_1$ multiplied by
$\prod_{j=m}^{b-1} Q(j) = Q_2$; the terms with $k \ge m$ are those of
$T_2$ multiplied by $\prod_{j=a}^{m-1} R(j) = R_1$.
:::

$R$ is needed only by left children, so the right spine of the tree
skips it. The result is $T(1, N+1)/Q(1, N+1)$ followed by one division
and the scalar factors.

## Leaves {#leaves}

The tree splits at leaf boundaries, and a leaf of $L$ terms (24 to 32,
or about 24 limbs of $T$) is evaluated exactly in one scratch buffer by
the backward recurrence
$$
T \gets P(k)\,Q_s + R(k)\,T, \qquad Q_s \gets Q(k)\,Q_s, \qquad R_s \gets R(k)\,R_s,
$$
for $k$ from the end of the leaf down. When $R$ divides $P$, $P = AR$
(as for $\pi$, $\log 2$ and most formulas), the step becomes
$T \gets R(k)(T + A(k)Q_s)$, where $A(k)$ is typically a single word even
when $P(k)$ is not. The polynomial values are computed by Horner's rule
in two's complement with a word count fixed by a bound on the values up
to the last term: whole leaves at a time for single-word values, then
`mpn_mul_1` and `mpn_addmul_1` passes with the signs folded in.

## Common factors {#content}

**Content powers.** When $Q$ (or $R$) has a large constant content $c$,
such as $q^2$ in the series of $\operatorname{atan}(p/q)$, the products
are carried without it, $Q(a, b) = c^{b-a}Q'(a, b)$, and the merges use
$T = T_1 Q_2' c^{b-m} + R_1 T_2$ with the pure powers $c^{L2^i}$ taken
from a table built by squarings over a perfectly balanced tree of $2^J$
leaves. This replaces a full-size product $Q_1Q_2$ by a short one but
adds a product by a power, so it pays only when the content dominates:
the table is used when the content has at least eight times the bits per
term of the rest of $Q$.

**Content removal.** A node only represents the ratios $T/Q$ and $R/Q$,
and $g = \gcd(R_1, Q_2)$ divides $T$, $Q$ and $R$ of the merge, so $R_1$
and $Q_2$ can be divided by $g$ before the products [@CHTZZ2007]. The
cancellation is large: the final $Q$ shrinks from 230 to 50 bits per
term for $\zeta(3)$, from 64 to 19 for Catalan's constant, from 34 to 14
for $\log 2$ and from 91 to 58 for $\pi$. The gcd is not computed by an
integer gcd but read off partial factorisations of $Q$ and $R$ over the
primes below $2^{32}$, kept with every node: $Q$ and $R$ are factored
over $\Z$, the values of their linear factors are factored by a
segmented sieve in $O(\sqrt N)$ memory, and a node's lists are merged
from its children's. A prime shared by $R(j)$ (left) and $Q(k)$ (right)
divides a fixed linear combination of $k - j$ for each pair of linear
factors, so only primes below a bound linear in the node's size need to
be compared, and larger ones wait until the level at which they can
cancel. The removal costs exact divisions, so it is done only where it
pays (not near the root, and not at small precisions).

## Truncation near the root {#truncation}

All values are `mp_real_t` balls, and exact ones stay exact as long as
they fit the working precision $w$; near the root the products exceed
$w$ limbs and are rounded, so the top levels of the tree are
floating-point computations with rigorous radii. Content removal stops
there, since it needs exact values.

::: lemma {#lem-vprec}
Let $e_1 = \operatorname{exp}(T_1) + \operatorname{exp}(Q_2)$ and
$e_2 = \operatorname{exp}(R_1) + \operatorname{exp}(T_2)$ be the sums of
the limb exponents of the factors, $g = e_1 - e_2$, and let $T_1Q_2$ be
rounded to $w$ limbs. Then at most $w - g + 1$ limbs of $R_1T_2$ lie
above the last limb of the rounded $T_1Q_2$.
:::

::: proof
A product of balls with limb exponents $e$ and $e'$ has its top limb at
position $e + e' - 1$ or $e + e' - 2$ (in units of $B$), so the rounded
$T_1Q_2$ ends at or above position $e_1 - 2 - w + 1 = e_1 - w - 1$, and
$R_1T_2$ starts at or below position $e_2 - 1$. The number of limb
positions from $e_1 - w - 1$ to $e_2 - 1$ is $e_2 - e_1 + w + 1 =
w - g + 1$.
:::

The sum takes nothing below the last limb of an inexact operand
([](algorithms/floating-point/mp-real-balls#sums)), so $R_1T_2$ is
formed at $w - g + 1$ limbs. For a convergent series, $R_1T_2/Q_1Q_2$ is the right partial sum times
the product of the left ratios, which is small: at the root $g$ is about
half the working precision. Computing $R_1T_2$ at $w - g + 1$ limbs cuts
that product by a third at $10^7$ digits of $\pi$.

## Tail bound {#tail}

The number of terms $N$ is the least for which a rigorous bound on the
tail is below the target. The product $\prod_{j<k}|R(j)/Q(j)|$ is
accumulated exactly (with bounded double rounding) for small $j$ and
bounded in closed form beyond the point where every polynomial is
dominated by its leading term; the terms beyond $N$ are bounded by a
geometric series from leading-term bounds. The bound is added to the
radius of the result, which therefore encloses the exact $S$. A series
whose ratio does not tend to a limit below 1 in absolute value, or whose
bound cannot be established, throws.

## Transform reuse in the merges {#transform-reuse}

At large sizes the products of a merge are FFT products, each costing
two forward transforms, a pointwise product and an inverse transform.
The four products of [](#prop-merge) come in two pairs with a common
operand:
$$
\{T_1 Q_2,\; Q_1 Q_2\} \text{ share } Q_2, \qquad
\{R_1 T_2,\; R_1 R_2\} \text{ share } R_1.
$$
Each pair is formed by `_fft_small_mpn_mul_shared`
([](fft_small/mpn_shared); see also
[multiplication with several primes](algorithms/fft/multimodular#shared)),
which computes the transform of the common operand once: three forward and two inverse transforms instead of four
and two. The pairs are those that the two threads of a parallel merge
take, so the sharing combines with parallelism.

A pair is shared when all three operands are exact and nonzero, both
products are kept whole (at most $w$ limbs), and
`_flint_mpn_mul_shared_use_fft` accepts the sizes (the longer product in
the FFT range of `flint_mpn_mul`, the shorter in half of it). Near the
root, where the products are truncated to middle products, separate
products are used. The second product of each pair is written over its
operand $Q_1$ (resp. $R_1$), and a shared pair needs no more memory than
one product, so the peak memory is unchanged. Measured on Emerald
Rapids for $\pi$, this makes the whole computation 6 to 10% faster at
$10^5$ and $10^6$ digits, 4 to 5% at $10^7$ and 3% at $10^8$ on one
thread, and 3 to 7% on two threads.

Summing $T_1Q_2 + R_1T_2$ in transform space, which would save one
more inverse transform, is not done: $Q_2$ lacks the low zero limbs
stripped from its mantissa (the 2-adic part of $Q$, for $\pi$ a third
of its length), so the two products are aligned differently, and
padding them to a common bottom lengthens every transform by more than
the inverse transform saves.

## Power tables with precomputed transforms {#powtab}

The bit-burst splittings of the exponential and of the sine and cosine
([](algorithms/floating-point/elementary-functions#bit-burst)) have
$P(k) = x$ (resp. $x^2$), a constant, so $R_1$ of a merge over
$[a, b)$ with midpoint $m$ is the power $x^{m-a}$. The tree splits at
fixed midpoints $m = a + \lfloor(b - a)/2\rfloor$, so only $O(\log N)$
distinct steps $s = m - a$ occur, and the powers $x^s$ are tabulated
once (a power table, as in `arb`). Every merge with step $s$ multiplies
$x^s$ by its right child, so the `fft_small` transform of $x^s$ can be
computed once per tree, and each such product saves the conversion and
forward transform of the table entry, about 30% of a two-prime product
(`_fft_small_mpn_mul_precomp`, [](fft_small/mpn_shared)).

::: lemma {#lem-step-uses}
In a binary splitting of $N$ terms with fixed midpoints, the ranges of
any level have at most two distinct lengths, differing by one.
Consequently the number of merges using each step $s$, and the longest
right child among them, can be computed in $O(\log N)$ operations from
$N$ alone.
:::

::: proof
By induction on the level: if all ranges of a level have length $\ell$
or $\ell + 1$, their halves have lengths $\lfloor\ell/2\rfloor$,
$\lceil\ell/2\rceil$, $\lfloor(\ell+1)/2\rfloor$ and
$\lceil(\ell+1)/2\rceil$, which take at most two consecutive values.
Tracking the (at most two) lengths and their multiplicities level by
level gives the counts.
:::

`_mp_real_bs_step_uses` computes these counts, and
`_mp_real_powtab_pre_init` transforms the entries used by at least two
merges, for partners up to the longest right child; the root's entry,
the largest, is used once and is not transformed, since that would only
add the storage of its transform. The transforms are made before the
recursion and only read during it, so the two halves of a merge can use
them on two threads. Their storage is one block, about four times the
transform of the largest entry kept; entries whose products exceed
$2^{20}$ limbs are not transformed, which bounds the extra memory to
about 50 MB per tree. Measured on Emerald Rapids, the exponential
becomes 8% faster at $10^6$ digits and 5% at $3 \cdot 10^5$, the sine and
cosine 4% and 3%; at $10^7$ digits the gain is about 3%, with 17% (exp)
and 8% (sine and cosine) more peak memory.

## Parallelism {#threads}

A splitting with a thread budget gives each half of a node half the
threads, as `flint_parallel_binary_splitting` does, but only while the
exact values of the node (estimated from the bits per term) stay within
four times the working precision (`MP_REAL_PAR_CAP`). Above that, two
concurrent halves would each hold full-size numbers, so they run one
after the other with the whole budget, and the merge runs its two pairs
of products on two threads, each pair's FFT using its half of the
threads. The result does not depend on the number of threads.

## Use in FLINT {#use}

- `mp_real_hypgeom_series` and `mp_real_hypgeom_series_int64`
  implement the general engine, with transform reuse in the merges;
  the constants of [](algorithms/floating-point/constants) and Zuniga's
  series for logarithms and arctangents
  ([](algorithms/floating-point/diophantine-reduction#zuniga)) are
  computed with it.
- `_mp_real_exp_sum_bs_powtab` and `_mp_real_sin_cos_sum_bs_powtab`
  implement the bit-burst splittings with transformed power tables
  (`_mp_real_bs_step_uses`, `_mp_real_powtab_pre_init`,
  `_mp_real_powtab_pre_mul`).
- `_mp_real_atan_frac_bsplit` and `_mp_real_atanh_frac_bsplit` are a
  dedicated splitting of the Taylor series of $\operatorname{atan}(p/q)$,
  using the content-power table for $q^2$.
- The shared and precomputed transforms are `_fft_small_mpn_mul_shared`
  and `_fft_small_mpn_mul_precomp` of `fft_small`; `flint_mpn_mul_shared`
  in `mpn_extras` is the user-level form of the first.

::: functions mp_real_hypgeom_series _mp_real_exp_sum_bs_powtab _mp_real_powtab_pre_init _fft_small_mpn_mul_shared _fft_small_mpn_mul_precomp :::
