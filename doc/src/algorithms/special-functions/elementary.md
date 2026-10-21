---
title: Elementary functions
summary: How arb evaluates exp, log, sin, cos and atan with rigorous error bounds - argument reduction modulo log 2 and π/4, precomputed tables, Taylor series by rectangular splitting on fixed-point numbers, reduction by logarithms of primes, Newton iteration and the bit-burst algorithm - and how the algorithm is chosen by precision.
complexity: $O(M(n) n^{1/2})$ or better at medium precision $n$ (rectangular splitting after table reduction); $O(M(n) \log^2 n)$ with the bit-burst algorithm at high precision.
order: 10
---

The real elementary functions of `arb` are built on four core functions:
the exponential, the logarithm, sine and cosine (computed together) and
the arctangent. All other real and complex elementary functions are
derived from these by algebraic identities. Each core function has
several algorithms, chosen by precision ([](#choice)). The design is
described in [@Joh2014c] (medium precision) and [@HJ2024] (high
precision); [@Bre1976b], [@Smi1989], [@BZ2011] and [@Mul2016] are
general references. The `mp_real` module implements its own versions of
the same functions with further refinements, described in
[](algorithms/floating-point/elementary-functions).

## The ball layer {#ball}

A function such as `arb_exp` receives a ball $[m \pm r]$ and proceeds as
follows.

1. Special values (zero, infinities, NaN, an infinite radius) are handled
   directly.
2. If $r > 0$, the precision is lowered to about the relative accuracy
   of the input plus 30 bits; more bits would be hidden by the propagated
   error. If the input has less than about 20 bits of accuracy, the
   function is evaluated at the endpoints of the ball instead, using its
   monotonicity (for $\sin$ and $\cos$, the endpoints and the extrema in
   between, in `double` arithmetic: `arb_sin_cos_wide`).
3. Otherwise the function of the exact midpoint is computed by a function
   such as `arb_exp_arf`, `arb_log_arf`, `arb_atan_arf` or the internal
   kernel of `arb_sin_cos`, which returns a ball with an exact error
   bound, and the propagated error
   ([](algorithms/ball-arithmetic/error-propagation)) is added.

Huge arguments are capped to keep the cost polynomial in the precision
(see [](guide/types/ball-arithmetic#precision)): `arb_exp` returns $[0 \pm
\infty]$ or a tiny ball when $|m| > 2^{\max(128, 2p)}$, and `arb_sin_cos`
returns $[\pm 1]$ when $|m| > 2^{\max(65536, 4p)}$. Arguments close to 0
are handled by Taylor polynomials of degree one or two with an explicit
remainder, which also preserves relative accuracy: for instance
$|e^x - (1 + x)| \le x^2$ for $|x| \le 1/2$.

## Argument reduction {#reduction}

At medium precision, the computations are done on fixed-point numbers:
arrays of $n$ limbs representing a number in $[0, 1)$ with an error
counted in units of $2^{-64n}$ (ulps). The tables below are stored to
4608 bits, which bounds the precision of this code path.

For the exponential, $x$ is reduced modulo $\log 2$:

$$
x = q \log 2 + w, \quad 0 \le w < \log 2, \quad e^x = 2^q e^w,
$$

where $q$ is an integer (`_arb_get_mpn_fixed_mod_log2`). Since $q$ is
computed exactly and $\log 2$ is read from a table with an error below
one ulp, the error of $w$ is at most $|q| + 1$ ulps; the precision is
raised by the number of bits of $q$, and if the table is too short for
that, the function falls back to the generic code. For sine and cosine,
$|x| = q \pi/4 + w$ with $0 \le w < \pi/4$
(`_arb_get_mpn_fixed_mod_pi4`), and the octant $q \bmod 8$ determines
which of $\pm \sin w$, $\pm \cos w$ and $\pm\cos(\pi/4 - w)$ etc. are
needed. For the logarithm, $x = 2^e (1 + t)$ with $0 \le t < 1$ and
$\log x = e \log 2 + \log(1 + t)$. For the arctangent, $|x| > 1$ is
replaced by $1/|x|$ using $\operatorname{atan} x = \pm\pi/2 -
\operatorname{atan}(1/x)$.

## Table lookups {#tables}

The reduced argument $w \in [0, 1)$ is reduced further with one or two
table lookups. Let $p$ be the integer formed by the top $b$ bits of $w$,
so that $w = p/2^b + v$ with $0 \le v < 2^{-b}$. The identities

$$
\begin{aligned}
e^w &= e^{p/2^b} e^{v}, \\
\sin w &= \sin(p/2^b) \cos v + \cos(p/2^b) \sin v, \quad
\cos w = \cos(p/2^b) \cos v - \sin(p/2^b) \sin v, \\
\log(1 + w) &= \log(1 + p/2^b) + \log\left(1 + \frac{2^b w - p}{2^b +
p}\right), \\
\operatorname{atan} w &= \operatorname{atan}(p/2^b) +
\operatorname{atan}\left(\frac{2^b w - p}{2^b + p w}\right)
\end{aligned}
$$

reduce the problem to evaluating the function at a small argument and
reading precomputed values at $p/2^b$. For $\log$ and $\operatorname{atan}$
the new argument is again below $2^{-b}$, since $0 \le 2^b w - p <
1$ and the denominators exceed $2^b$. Up to 512 bits a single lookup with
$b = 8$ is used (`arb_exp_tab1`, `arb_sin_cos_tab1`, `arb_atan_tab1`; for
the logarithm two 7-bit lookups, `arb_log_tab11` and `arb_log_tab12`); up
to 4608 bits, two successive lookups with $b = 5$ and then 5 more bits
(`arb_exp_tab21` and `arb_exp_tab22`, and so on). After the reduction the
argument is below $2^{-8}$ or $2^{-10}$, and the tables are small: the
exponential tables have 178 entries of 512 bits and $23 + 32$ entries of
4608 bits.

For the logarithm, the series that is summed is not that of $\log(1 +
v)$ but the faster one of
$$
\log(1 + v) = 2 \operatorname{atanh}\frac{v}{2 + v} = 2 \sum_{k=0}^{\infty}
\frac{1}{2k+1} \left(\frac{v}{2 + v}\right)^{2k+1},
$$
which has only odd powers and an argument half as large; the second
table lookup is fused with this transformation.

## Taylor series by rectangular splitting {#series}

After reduction, the Taylor series of $e^v$, $\sin v$ and $\cos v$, or
$\operatorname{atan} v$ and $\operatorname{atanh} v$, is summed on
fixed-point numbers by `_arb_exp_taylor_rs`, `_arb_sin_cos_taylor_rs`
and `_arb_atan_taylor_rs`. With $|v| \le 2^{-r}$ and a target of $p$
bits, about $p/r$ terms are needed for the exponential (with the number
of terms from `_arb_exp_taylor_bound`) and $p/(2r)$ for the odd series.

::: lemma Truncation of the exponential series {#lem-exp-tail}
If $0 \le v \le 1/2$ and $N \ge 1$, then $\sum_{k=N}^{\infty} v^k / k!
\le 2 v^N / N!$.
:::

::: proof
Consecutive terms have ratio $v/(k+1) \le 1/2$ for $k \ge N - 1 \ge 0$,
so the tail is bounded by the geometric series $\frac{v^N}{N!} (1 + 1/2 +
1/4 + \cdots)$.
:::

For the alternating series of $\sin$, $\cos$ and $\operatorname{atan}$,
the tail is bounded by its first term, and for $\operatorname{atanh}$ by
the first term times $1/(1 - v^2)$.

The sum $\sum_{k<N} c_k v^k$ is evaluated by rectangular splitting
[@PS1973] [@Smi1989]: the powers $v^2, \ldots, v^m$ with $m \approx
\sqrt{N}$ are computed once, and the polynomial is evaluated as a
polynomial in $v^m$ whose coefficients are linear combinations of $1, v,
\ldots, v^{m-1}$. This needs about $m + N/m \approx 2\sqrt{N}$ full
multiplications, and $N$ cheap scalar multiplications by the coefficients.
The scalar steps are made cheap in two ways. The coefficients $1/k!$ (or
$1/(2k+1)$) of a block are written over a common denominator that fits in
one word, with numerators that also fit in one word (the tables
`factorial_tab_numer` and `factorial_tab_denom` in
`arb/exp_taylor_rs.c`), so each block costs a few `mpn_addmul_1` calls
and one division by a single word; and because the arithmetic is
fixed-point, there are no normalisations or roundings in between. Each
operation adds at most one ulp of error, and the total is returned as an
explicit bound.

Two further tricks halve the work. When many terms are needed, the
exponential is computed as $e^v = \sinh v + \sqrt{1 + \sinh^2 v}$ from the
odd series of $\sinh v$ (`alternating` set to 0), and similarly $\cos v =
\sqrt{1 - \sin^2 v}$, each at the cost of one square root. Finally the
table values are multiplied in (one or two multiplications, or the
addition formulas for $\sin$ and $\cos$), $e \log 2$ or $\pi/2$ is added
for $\log$ and $\operatorname{atan}$, and the fixed-point result is
rounded into the midpoint with the accumulated error in the radius.

## High precision {#high}

Above the range of the fixed-point code, or when its tables are not
precise enough, the following methods are used.

*Reduction by logarithms of primes.* For $e^x$ between 2240 bits and 4
million bits (`arb_exp_arf_log_reduction`), the argument is reduced by an
integer combination of the logarithms of the first 13 primes,
$$
x = \sum_{i=1}^{13} c_i \log p_i + w, \qquad e^x = \Bigl(\prod_i
p_i^{c_i}\Bigr) e^w,
$$
where the $c_i$ are chosen with precomputed short integer relations
between the $\log p_i$ so that $|w|$ is about $2^{-256}$ to $2^{-768}$
[@HJ2024]. The logarithms are computed once per thread by binary
splitting of Machin-type formulas (`arb_log_primes_vec_bsplit`) and
cached, and the rational number $\prod p_i^{c_i}$ is computed exactly.
Sine and cosine between 2600 bits and 4 million bits use the analogous
reduction by the arguments $\theta_i$ of 13 Gaussian primes $a_i + b_i
i$ and multiply by the Gaussian rational $\prod (a_i + b_i i)^{c_i}$
normalised to modulus 1 (`arb_sin_cos_arf_atan_reduction`).
The relations are the same as those used by `mp_real`
([](algorithms/floating-point/diophantine-reduction)).

*Newton iteration.* Above 2800 bits the logarithm is computed from the
exponential: with $t$ an approximation of $\log x$ to $1/n$ of the
target precision,
$$
\log x = t + \log(1 + \delta), \qquad \delta = x e^{-t} - 1,
$$
and $\log(1 + \delta)$ is evaluated by its Taylor series to order $n$
(between 4 and 12), so that one step multiplies the precision by $n$
(`arb_log_newton`). Above 3400 bits the arctangent is computed in the same
way from $\sin$ and $\cos$, with $\operatorname{atan} x = t +
\operatorname{atan}((x \cos t - \sin t)/(\cos t + x \sin t))$
(`arb_atan_newton`). Logarithms of small smooth integers at high precision
are assembled from the cached logarithms of primes instead.

*Generic algorithms.* `arb_exp_arf_generic` and
`arb_sin_cos_arf_generic` work at any precision without tables. They
reduce by $\log 2$ or $\pi/4$ using `arb_const_log2` and `arb_const_pi`,
then either halve the argument $q$ times, sum the Taylor series by
rectangular splitting in ball arithmetic and square (or apply the double
angle formula) $q$ times (`arb_exp_arf_rs_generic`,
`arb_sin_cos_arf_rs_generic`), or use the bit-burst algorithm. The choice
depends on the precision, the size of the argument, the length of its
mantissa (a short argument favours binary splitting) and the number of
available threads.

::: algorithm Bit-burst evaluation of the exponential {#alg-bb}
**Input:** a fixed-point $x$ with $|x| < 2^{-d}$ and $p$ fractional bits.

**Output:** $e^x$ with an error bound.

1. Write $x = x_0 + x_1 + \cdots + x_K$, where $x_0$ consists of the
   first $d$ significant bits of $x$ and $x_k$ of the next $2^k d$ bits.
2. For each $k$, compute $e^{x_k}$ by summing its Taylor series, a
   series with a rational argument of $2^k d$ bits and denominator
   $2^{2^{k+1} d}$, by binary splitting (`_arb_exp_sum_bs_powtab`),
   with enough terms for $p$ bits.
3. Return $\prod_k e^{x_k}$.
:::

Since $|x_k| < 2^{-2^k d}$, the series for $e^{x_k}$ needs about
$p/(2^k d)$ terms, whose numerators grow by $2^k d$ bits each, so each
binary splitting costs $O(M(p) \log p)$ and the whole algorithm $O(M(p)
\log^2 p)$ [@Bre1976b]. The same scheme gives $\sin$ and $\cos$
(`arb_sin_cos_arf_bb`, summing the series of $e^{ix_k}$) and, with the
addition formula $\operatorname{atan} x = \operatorname{atan}(p/q) +
\operatorname{atan}((qx - p)/(q + px))$ in place of the product, the
arctangent (`arb_atan_arf_bb`). The terms of the exponentials can be
computed in parallel on several threads.

## Choice of algorithm {#choice}

The following table summarises the default choices for the midpoint
functions at precision $p$ bits (thresholds from `arb.h`).

| Function | $p < 512$ | $512 \le p \lesssim 2240$-$4608$ | to $4 \cdot 10^6$ | above |
|---|---|---|---|---|
| exp | tables ($b = 8$) + series | tables ($b = 5 + 5$) + series, to 2240 bits | log-prime reduction | generic (rs or bit-burst) |
| sin, cos | tables ($b = 8$) + series | tables ($b = 5 + 5$) + series, to 2600 bits | Gaussian-prime reduction | generic (rs or bit-burst) |
| log | tables ($7 + 7$) + atanh series | tables ($5 + 5$) + series, to 2800 bits | Newton on exp | Newton on exp |
| atan | tables ($b = 8$) + series | tables ($5 + 5$) + series, to 3400 bits | Newton on sin, cos | Newton on sin, cos |

## Derived functions {#derived}

- `arb_expm1` and `arb_log1p` use the same code with the subtraction of 1
  built in, adding the bits that would be lost to cancellation for small
  arguments.
- `arb_sinh_cosh`, `arb_tanh` and `arb_coth` are computed from `arb_exp`
  or `arb_expm1`; `arb_asinh`, `arb_acosh` and `arb_atanh` from
  `arb_log`, `arb_log1p` and `arb_sqrt`.
- `arb_asin` and `arb_acos` use $\operatorname{atan}(x / \sqrt{1 -
  x^2})$, and `arb_atan2` combines `arb_atan` with $\pi$.
- `arb_sin_pi` and `arb_cos_pi` reduce $x$ modulo $1/2$ exactly (by
  subtracting the integer nearest to $2x$) before multiplying by $\pi$, so
  they are exact at integers and half-integers and accurate near the
  zeros. `arb_sin_cos_pi_fmpq` evaluates small denominators as algebraic
  numbers, by refining a root of the minimal polynomial of $\cos(2\pi/q)$
  with Newton iteration.
- The complex functions of `acb` (`acb_exp`, `acb_log`, `acb_sin_cos`,
  `acb_atan`, ...) are built from the real ones: for example $e^{a + bi}
  = e^a (\cos b + i \sin b)$ and $\log z = \log|z| + i \arg z$ with
  `arb_log_hypot` and `arb_atan2`.

## Use in FLINT {#use}

::: functions arb_exp arb_log arb_sin_cos arb_atan arb_exp_arf arb_log_arf arb_atan_arf _arb_exp_taylor_rs _arb_sin_cos_taylor_rs _arb_atan_taylor_rs arb_exp_arf_log_reduction arb_sin_cos_arf_atan_reduction arb_log_newton arb_atan_newton arb_exp_arf_bb arb_sin_cos_arf_bb arb_atan_arf_bb :::
