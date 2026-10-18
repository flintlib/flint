---
title: double_extras.h
header: double_extras.h
description: Helper functions for double-precision floating-point numbers
category: Floating-point arithmetic
order: 5
---

This module collects small helpers for IEEE 754 binary64 numbers
(`double`) that are used throughout FLINT: random generation for tests,
exact scaling by powers of two, Horner evaluation, the Lambert W
function and a base-2 logarithm. The functions assume round-to-nearest
arithmetic. Vectors and matrices of doubles are in `d_vec` and `d_mat`,
interval arithmetic with doubles in `double_interval`, and double-word
to quad-double arithmetic in `dfloat`;
[](algorithms/floating-point/overview) explains how these relate to
FLINT's other floating-point types.

## Constants and bit-level access {#constants}

::: macro
D_BITS
D_EPS

The number of bits of a `double` mantissa, 53, and the distance
$2^{-52}$ from 1 to the next larger double (`DBL_EPSILON`). The unit
roundoff of round-to-nearest arithmetic is $u = 2^{-53}$, half of
`D_EPS`.
:::

::: macro
D_INF
D_NAN

Positive infinity (`HUGE_VAL`) and a quiet NaN (computed as
`HUGE_VAL - HUGE_VAL`).
:::

::: macro
D_MIN_NORMAL_EXPONENT
D_MAX_NORMAL_EXPONENT
D_EXPONENT_BIAS
D_EXPONENT_SHIFT

The range $-1022 \le e \le 1023$ of the exponents $e$ of normal doubles
$\pm 2^e \cdot 1.f$, the exponent bias 1023, and the position 52 of the
exponent field in the bit pattern.
:::

::: type
double_uint64_u

A union of a `double` (member `f`) and a `uint64_t` (member `i`), for
reading and writing the bit pattern of a double.
:::

## Random generation {#random}

These functions are intended for test code. They are not uniform: the
mantissa bits come from `n_randtest`, which favours long runs of zeros
and ones.

::: function
double d_randtest(flint_rand_t state)

Returns a random double in $[0.5, 1)$ with a random 53-bit mantissa.
:::

::: function
double d_randtest_signed(flint_rand_t state, slong minexp, slong maxexp)

Returns $m 2^e$, $-m 2^e$ or zero with probability 1/3 each, where $m$
is a value of `d_randtest` and $e$ is uniform in
$[\mathit{minexp}, \mathit{maxexp}]$, so that a nonzero result has
absolute value in $[2^{e-1}, 2^e)$. Requires
$\mathit{minexp} \le \mathit{maxexp}$; values outside the double range
overflow to infinity or underflow (through `ldexp`).
:::

::: function
double d_randtest_special(flint_rand_t state, slong minexp, slong maxexp)

Like `d_randtest_signed`, but returns a positive number, a negative
number, zero or a special value with probability 1/4 each. The special
value is NaN with probability 1/2 and $+\infty$ or $-\infty$ with
probability 1/4 each.
:::

## Arithmetic {#arithmetic}

::: function
double d_polyval(const double * poly, int len, double x)

Evaluates the polynomial $\sum_{i=0}^{\mathit{len}-1} \mathit{poly}[i]\, x^i$
at $x$ by Horner's rule, rounding after each operation. Requires
$\mathit{len} \ge 1$.
:::

::: function
double d_mul_2exp(double x, int i)
double d_mul_2exp_inrange(double x, int i)
double d_mul_2exp_inrange2(double x, int i)

Return $x \cdot 2^i$. The result is exact unless it overflows or is
subnormal. `d_mul_2exp` accepts any $i$ and calls `ldexp` when $2^i$
is not a normal double. The other two are faster special cases with
unchecked preconditions.

::: details
`d_mul_2exp_inrange` requires $-1022 \le i \le 1023$. It builds the
double $2^i$ from its bit pattern and performs one multiplication, so
overflow, underflow and special values of $x$ behave as for an ordinary
multiplication.

`d_mul_2exp_inrange2` adds $i$ to the exponent field of $x$ with an
integer addition. It additionally requires that $x$ and $x \cdot 2^i$
are normal doubles, in particular $x \ne 0$; otherwise the result is
meaningless. The preconditions of both functions are checked only in
builds configured with `--enable-assert`.
:::
:::

## Special functions {#special}

::: function
double d_lambertw(double x)

Returns the principal branch $W_0(x)$ of the Lambert W function, the
solution $w \ge -1$ of $w e^w = x$. Returns NaN for $x < -1/e$, where
the solution is complex, and for $x = -\infty$; returns $x$ itself for
$x = 0$, $x = +\infty$ and NaN.

::: details
Close to zero, the result is $x - x^2$ (for $|x| < 10^{-9}$) or a
rational approximation. Close to the branch point $-1/e$, it is a
polynomial in $\sqrt{x + 1/e}$, where $x + 1/e$ is computed with a
correction term for the rounding of $1/e$ (the most negative double in
the domain is $-6627126856707895 \cdot 2^{-54}$). Elsewhere, an initial
value from a piecewise rational approximation, or for $x > 40$ from the
asymptotic expansion in $L_1 = \log x$ and $L_2 = \log \log x$
[@CGHJK1996], is refined by one or two steps of Halley's method,
rescaled to avoid overflow of $e^w$ near the top of the exponent range.

A test of $10^7$ random inputs showed a maximum relative error below
$0.95 \cdot 2^{-52}$ for positive $x$; for negative $x$ the error is
larger, up to about $10 \cdot 2^{-52}$ close to $-1/e$. The accuracy
depends on the compiler flags and on the system `exp` and `log`. For a
rigorous enclosure at any precision, use `arb_lambertw`.
:::
:::

::: function
double d_log2(double x)

Returns an approximation of $\log_2 x$, computed as `log(x)` divided
by the double nearest to $\log 2$, so with an error of a few ulps. The
result is not exact for all powers of two: for $x = 2^{-1023}$ it is
$-1023.0000000000001$. Special values and domain errors are those of
`log`.
:::

::: function
int d_is_nan(double x)

Returns nonzero if $x$ is a NaN (with any payload) and 0 otherwise.
:::
