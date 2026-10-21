---
title: dfloat.h
header: dfloat.h
description: Double-word to quad-double floating-point numbers and balls built on IEEE doubles
category: Floating-point arithmetic
order: 2
---

This module implements real and complex numbers represented as
unevaluated sums of one to four IEEE doubles (double, double-double,
triple-double and quad-double, about $53N$ bits for $N$ doubles),
together with rigorous ball arithmetic on top of them. All arithmetic is
built from double operations and correctly rounded fused multiply-adds,
with SIMD code for vectors and dot products, and is intended to be much
faster than `arb_t` at precisions up to about 200 bits while remaining
rigorous. The formats are available directly, through functions such as
`d2_mul` and `d3b_exp`, and as generic rings created with
`gr_ctx_init_dfloat`, so that `gr_poly`, `gr_mat` and the other generic
modules work over them.

There are sixteen formats, with $N = 1, 2, 3, 4$:

| Format | Type | Represents |
|---|---|---|
| `dN` | `dN_t` | a real number $d_0 + d_1 + \cdots + d_{N-1}$ (approximate, like `double`) |
| `dNb` | `dNb_t` | a real ball: a `dN` midpoint and a `double` radius (rigorous) |
| `dNc` | `dNc_t` | a complex number with `dN` real and imaginary parts |
| `dNcb` | `dNcb_t` | a rectangular complex ball with `dNb` real and imaginary parts |

The plain formats are approximate arithmetic in the style of `double`:
results are rounded (to a few ulps at $53N$ bits), overflow gives
infinity, and infinities and NaN propagate. The ball formats are
rigorous: every operation returns a ball containing the exact image of
the input balls, including where the double arithmetic underflows or
overflows, and exact operations on exact balls stay exact. The guide
[](dfloat/using) shows both interfaces; the conventions of the balls for
infinities, NaN and partial functions are set out in [](dfloat/balls);
[](dfloat/performance) compares speeds with `nfloat` and `arb_t`. The
algorithms (error-free transformations, renormalisation, division and
square root, and the error analysis) are described in
[](algorithms/floating-point/double-word), and the elementary functions
in [](algorithms/floating-point/dfloat-elementary).

In the text, `dN_exp` stands for the four functions `d1_exp`, ...,
`d4_exp`, and similarly for the other formats. Aliasing of inputs and
outputs is allowed everywhere.

::: summary algorithms/floating-point/overview :::

## Types {#types}

::: type
d1_struct
d2_struct
d3_struct
d4_struct
d1_t
d2_t
d3_t
d4_t
d1_ptr
d2_ptr
d3_ptr
d4_ptr
d1_srcptr
d2_srcptr
d3_srcptr
d4_srcptr

A `dN_struct` is a struct with a single member `double d[N]`, the
components of the expansion $x = d_0 + \cdots + d_{N-1}$, and a `dN_t`
is a length-one array of it. The components are roughly nonoverlapping
($|d_{i+1}|$ is at most about $\operatorname{ulp}(d_i)$), the head $d_0$
is the double nearest to $x$ (up to ties), and $d_0 = 0$ implies
$x = 0$. Results of the arithmetic are *weakly* normalised; see
`dN_canonicalise` for the canonical form. `d1` is a plain double,
provided for uniformity.
:::

::: type
d1b_struct
d2b_struct
d3b_struct
d4b_struct
d1b_t
d2b_t
d3b_t
d4b_t
d1b_ptr
d2b_ptr
d3b_ptr
d4b_ptr
d1b_srcptr
d2b_srcptr
d3b_srcptr
d4b_srcptr

A `dNb_struct` has the members `double d[N]` (the midpoint, a `dN`
expansion) and `double rad` (the radius), and represents every real
number $x$ with $|x - \mathrm{mid}| \le \mathrm{rad}$. In a ball
returned by this module, the midpoint components are finite and the
radius is finite and nonnegative, or the ball is the whole real line in
the normal form $[0 \pm \infty]$; any ball with a nonfinite component or
radius is accepted as the whole real line (see [](dfloat/balls)).
:::

::: macro
DFLOAT_MAX_N

The largest number of components, 4.
:::

## Support and generic rings {#rings}

::: function
int dfloat_is_supported(void)

Returns whether the module works on this target. It requires that the
module was compiled without contraction of floating-point expressions
(`-ffp-contract=off` or the equivalent pragma), with doubles evaluated
in double precision (`FLT_EVAL_METHOD` 0, so not with x87 arithmetic),
and that `fma` is correctly rounded, which is tested once at run time.
Where it returns 0, the other functions of this module must not be
called; `gr_ctx_init_dfloat` returns `GR_UNABLE` and the power sums
return 0, so that their callers fall back to other code.
:::

::: function
int gr_ctx_init_dfloat(gr_ctx_t ctx, int n, int flags)

Initialises `ctx` to a generic ring of $n$-term expansions, $1 \le n \le 4$.
The flags select the format: without `DFLOAT_BALL` the elements are
plain `dN` numbers (an approximate floating-point ring), with it `dNb`
balls (a rigorous representation of $\mathbb{R}$); with `DFLOAT_COMPLEX` the
complex formats `dNc` or `dNcb`. Passing 0 or 1 selects the plain or
the ball real ring with the defaults. Returns `GR_UNABLE` if $n$ is out
of range or `dfloat_is_supported` returns 0, and `GR_SUCCESS`
otherwise.

::: details
The element size is $8(n + b)$ bytes, with $b = 1$ for balls and 0
otherwise, and twice that for the complex rings. The real ball ring
reports itself as a field, ordered, real, inexact and non-canonical,
like `gr_ctx_init_real_arb`, and the complex ball ring as an
algebraically closed field like `gr_ctx_init_complex_acb`. The
precision reported by `gr_ctx_get_real_prec` is $53n$ bits and cannot
be changed.

With `DFLOAT_STRONG`, the arithmetic operations (`add`, `sub`, `mul`,
`sqr`, `div`, `inv`, `sqrt`, `rsqrt`, `exp`, `log`, `sin`, `cos`,
`sin_cos`, `atan` and the vector operations) canonicalise their results;
this makes `equal` on plain numbers a comparison of values and costs up
to 60% on the cheapest operations for $n = 3, 4$. With `DFLOAT_FAST`
(balls only), products bound their rounding errors instead of tracking
them exactly when both inputs are exact; this is slightly cheaper but
loses exactness ($(2^{60} + 1) \cdot 3$ then gets a radius of about
$2^{-46}$).

The rings implement only operations with native code: the arithmetic,
the elementary functions documented below and their vector versions.
Special functions such as `gr_gamma` and `gr_zeta` return `GR_UNABLE`;
use an `arb_t` ring for them. The constants `gr_pi`, `gr_euler` and
`gr_catalan` are table lookups. Conversions with `gr_set_other` exist
between the `dfloat` rings and from and to `arb_t`, `acb_t`, `fmpz`
and `fmpq`; see [](dfloat/balls) for the rules on nonfinite values and
status codes.
:::
:::

::: macro
DFLOAT_BALL
DFLOAT_STRONG
DFLOAT_FAST
DFLOAT_COMPLEX

Flags for `gr_ctx_init_dfloat`, combined with bitwise or: ball
arithmetic instead of plain numbers (1), canonical results (2), bounded
instead of tracked rounding errors in ball products (4), and complex
instead of real numbers (8).
:::

## Assignment and constants {#assignment}

::: function
void d1_set_d(d1_t res, double x)
void d2_set_d(d2_t res, double x)
void d3_set_d(d3_t res, double x)
void d4_set_d(d4_t res, double x)
void d1b_set_d(d1b_t res, double x)
void d2b_set_d(d2b_t res, double x)
void d3b_set_d(d3b_t res, double x)
void d4b_set_d(d4b_t res, double x)

Sets `res` to the double $x$, exactly. In the ball version a nonfinite
$x$ gives the whole real line.
:::

::: function
void d1_set_dd(d1_t res, double x0, double x1)
void d2_set_dd(d2_t res, double x0, double x1)
void d3_set_dd(d3_t res, double x0, double x1)
void d4_set_dd(d4_t res, double x0, double x1)
void d1b_set_dd(d1b_t res, double x0, double x1)
void d2b_set_dd(d2b_t res, double x0, double x1)
void d3b_set_dd(d3b_t res, double x0, double x1)
void d4b_set_dd(d4b_t res, double x0, double x1)

Sets `res` to $x_0 + x_1$, renormalised. In the ball version, anything
that does not fit in $N$ components is added to the radius.
:::

::: function
void d1b_set_d_rad(d1b_t res, double x, double rad)
void d2b_set_d_rad(d2b_t res, double x, double rad)
void d3b_set_d_rad(d3b_t res, double x, double rad)
void d4b_set_d_rad(d4b_t res, double x, double rad)

Sets `res` to the ball $[x \pm \mathit{rad}]$; the whole real line if
$x$ or $\mathit{rad}$ is not finite or $\mathit{rad}$ is negative.
:::

::: function
void d1_zero(d1_t res)
void d2_zero(d2_t res)
void d3_zero(d3_t res)
void d4_zero(d4_t res)
void d1_one(d1_t res)
void d2_one(d2_t res)
void d3_one(d3_t res)
void d4_one(d4_t res)
void d1b_zero(d1b_t res)
void d2b_zero(d2b_t res)
void d3b_zero(d3b_t res)
void d4b_zero(d4b_t res)
void d1b_one(d1b_t res)
void d2b_one(d2b_t res)
void d3b_one(d3b_t res)
void d4b_one(d4b_t res)
void d1b_indeterminate(d1b_t res)
void d2b_indeterminate(d2b_t res)
void d3b_indeterminate(d3b_t res)
void d4b_indeterminate(d4b_t res)

Set `res` to 0, to 1, or (for `indeterminate`) to the whole real line
in its normal form $[0 \pm \infty]$.
:::

::: function
void d1_set(d1_t res, d1_srcptr x)
void d2_set(d2_t res, d2_srcptr x)
void d3_set(d3_t res, d3_srcptr x)
void d4_set(d4_t res, d4_srcptr x)
void d1_neg(d1_t res, d1_srcptr x)
void d2_neg(d2_t res, d2_srcptr x)
void d3_neg(d3_t res, d3_srcptr x)
void d4_neg(d4_t res, d4_srcptr x)
void d1_abs(d1_t res, d1_srcptr x)
void d2_abs(d2_t res, d2_srcptr x)
void d3_abs(d3_t res, d3_srcptr x)
void d4_abs(d4_t res, d4_srcptr x)
void d1b_set(d1b_t res, d1b_srcptr x)
void d2b_set(d2b_t res, d2b_srcptr x)
void d3b_set(d3b_t res, d3b_srcptr x)
void d4b_set(d4b_t res, d4b_srcptr x)
void d1b_neg(d1b_t res, d1b_srcptr x)
void d2b_neg(d2b_t res, d2b_srcptr x)
void d3b_neg(d3b_t res, d3b_srcptr x)
void d4b_neg(d4b_t res, d4b_srcptr x)
void d1b_abs(d1b_t res, d1b_srcptr x)
void d2b_abs(d2b_t res, d2b_srcptr x)
void d3b_abs(d3b_t res, d3b_srcptr x)
void d4b_abs(d4b_t res, d4b_srcptr x)

Copy, negation and absolute value, all exact. `dNb_set` copies any form
of the whole real line as it is; `dNb_neg` and `dNb_abs` return its
normal form.
:::

## Arithmetic {#arithmetic}

::: function
void d1_add(d1_t res, d1_srcptr x, d1_srcptr y)
void d2_add(d2_t res, d2_srcptr x, d2_srcptr y)
void d3_add(d3_t res, d3_srcptr x, d3_srcptr y)
void d4_add(d4_t res, d4_srcptr x, d4_srcptr y)
void d1_sub(d1_t res, d1_srcptr x, d1_srcptr y)
void d2_sub(d2_t res, d2_srcptr x, d2_srcptr y)
void d3_sub(d3_t res, d3_srcptr x, d3_srcptr y)
void d4_sub(d4_t res, d4_srcptr x, d4_srcptr y)
void d1_mul(d1_t res, d1_srcptr x, d1_srcptr y)
void d2_mul(d2_t res, d2_srcptr x, d2_srcptr y)
void d3_mul(d3_t res, d3_srcptr x, d3_srcptr y)
void d4_mul(d4_t res, d4_srcptr x, d4_srcptr y)
void d1_sqr(d1_t res, d1_srcptr x)
void d2_sqr(d2_t res, d2_srcptr x)
void d3_sqr(d3_t res, d3_srcptr x)
void d4_sqr(d4_t res, d4_srcptr x)
void d1_mul_d(d1_t res, d1_srcptr x, double y)
void d2_mul_d(d2_t res, d2_srcptr x, double y)
void d3_mul_d(d3_t res, d3_srcptr x, double y)
void d4_mul_d(d4_t res, d4_srcptr x, double y)
void d1_div(d1_t res, d1_srcptr x, d1_srcptr y)
void d2_div(d2_t res, d2_srcptr x, d2_srcptr y)
void d3_div(d3_t res, d3_srcptr x, d3_srcptr y)
void d4_div(d4_t res, d4_srcptr x, d4_srcptr y)
void d1_inv(d1_t res, d1_srcptr x)
void d2_inv(d2_t res, d2_srcptr x)
void d3_inv(d3_t res, d3_srcptr x)
void d4_inv(d4_t res, d4_srcptr x)
void d1_sqrt(d1_t res, d1_srcptr x)
void d2_sqrt(d2_t res, d2_srcptr x)
void d3_sqrt(d3_t res, d3_srcptr x)
void d4_sqrt(d4_t res, d4_srcptr x)
void d1_rsqrt(d1_t res, d1_srcptr x)
void d2_rsqrt(d2_t res, d2_srcptr x)
void d3_rsqrt(d3_t res, d3_srcptr x)
void d4_rsqrt(d4_t res, d4_srcptr x)
void d1_mul_2exp_si(d1_t res, d1_srcptr x, slong e)
void d2_mul_2exp_si(d2_t res, d2_srcptr x, slong e)
void d3_mul_2exp_si(d3_t res, d3_srcptr x, slong e)
void d4_mul_2exp_si(d4_t res, d4_srcptr x, slong e)

Arithmetic on the plain formats, accurate to a few ulps at $53N$ bits:
sum, difference, product, square, product with a double, quotient,
inverse, square root, reciprocal square root $x^{-1/2}$ and
multiplication by $2^e$. Overflow gives infinity, and nonfinite values
propagate as in double arithmetic.

::: details
Addition, subtraction and multiplication compute the exact terms of
each order with error-free transformations and renormalise; terms of
order $N$ and beyond are dropped in multiplication. Division is a long
division whose digits are single double divisions, and the square roots
are Newton iterations with tapered precision; `dN_rsqrt` doubles the
number of correct components per step and is accurate to about
$8 u^N$, $u = 2^{-53}$. See [](algorithms/floating-point/double-word).

The results are *compact* weakly normalised expansions: a component
that comes out exactly zero after a cancellation is squeezed out, so
that the magnitudes of the components always decrease geometrically.
The tapered kernels (division, square roots and the elementary
functions) rely on this and lose accuracy on an expansion with an
internal zero such as $(2^{-86}, 0, 2^{-140}, 0)$, which the module
never produces; an expansion built by hand in this form should be
passed through `dN_canonicalise` first (the elementary functions do
this themselves, but `dN_div`, `dN_sqrt` and `dN_rsqrt` do not).
:::
:::

::: function
void d1b_add(d1b_t res, d1b_srcptr x, d1b_srcptr y)
void d2b_add(d2b_t res, d2b_srcptr x, d2b_srcptr y)
void d3b_add(d3b_t res, d3b_srcptr x, d3b_srcptr y)
void d4b_add(d4b_t res, d4b_srcptr x, d4b_srcptr y)
void d1b_sub(d1b_t res, d1b_srcptr x, d1b_srcptr y)
void d2b_sub(d2b_t res, d2b_srcptr x, d2b_srcptr y)
void d3b_sub(d3b_t res, d3b_srcptr x, d3b_srcptr y)
void d4b_sub(d4b_t res, d4b_srcptr x, d4b_srcptr y)
void d1b_mul(d1b_t res, d1b_srcptr x, d1b_srcptr y)
void d2b_mul(d2b_t res, d2b_srcptr x, d2b_srcptr y)
void d3b_mul(d3b_t res, d3b_srcptr x, d3b_srcptr y)
void d4b_mul(d4b_t res, d4b_srcptr x, d4b_srcptr y)
void d1b_sqr(d1b_t res, d1b_srcptr x)
void d2b_sqr(d2b_t res, d2b_srcptr x)
void d3b_sqr(d3b_t res, d3b_srcptr x)
void d4b_sqr(d4b_t res, d4b_srcptr x)
void d1b_mul_d(d1b_t res, d1b_srcptr x, double y)
void d2b_mul_d(d2b_t res, d2b_srcptr x, double y)
void d3b_mul_d(d3b_t res, d3b_srcptr x, double y)
void d4b_mul_d(d4b_t res, d4b_srcptr x, double y)
void d1b_mul_2exp_si(d1b_t res, d1b_srcptr x, slong e)
void d2b_mul_2exp_si(d2b_t res, d2b_srcptr x, slong e)
void d3b_mul_2exp_si(d3b_t res, d3b_srcptr x, slong e)
void d4b_mul_2exp_si(d4b_t res, d4b_srcptr x, slong e)

Rigorous ball arithmetic. The rounding errors of the midpoint operation
are added to the radius, so that exact operations on exact balls give
exact results: $[3 \pm 0] \cdot [3 \pm 0] = [9 \pm 0]$, and
$(2^{60} + 1) \cdot 3$ is exact in `d2b`.

::: details
The midpoint terms are computed with error-free transformations and the
magnitudes of the dropped terms are added to the radius; the radius
itself is computed in round-to-nearest arithmetic and multiplied by
$1 + 2^{-45}$ at the end, which covers up to $2^7$ rounding errors.
When both inputs are exact the products track all rounding errors
exactly; otherwise (or always, in a context with `DFLOAT_FAST`) the
last-order terms are bounded by $u$ times the accumulated values, which
is cheaper and loses only a constant factor on a radius that is nonzero
anyway. Products whose terms could underflow take a slower path that
scales the operands to unit exponent (`_dfloat_ball_mul_scaled`). A
midpoint that overflows makes the result the whole real line.
:::
:::

::: function
int d1b_div(d1b_t res, d1b_srcptr x, d1b_srcptr y)
int d2b_div(d2b_t res, d2b_srcptr x, d2b_srcptr y)
int d3b_div(d3b_t res, d3b_srcptr x, d3b_srcptr y)
int d4b_div(d4b_t res, d4b_srcptr x, d4b_srcptr y)
int d1b_inv(d1b_t res, d1b_srcptr x)
int d2b_inv(d2b_t res, d2b_srcptr x)
int d3b_inv(d3b_t res, d3b_srcptr x)
int d4b_inv(d4b_t res, d4b_srcptr x)

Division and inverse. Return `GR_DOMAIN` if the divisor is exactly zero
(whatever the dividend), `GR_UNABLE` if it contains zero (in particular
if it is the whole real line), and `GR_SUCCESS` otherwise, in which
case `res` contains the quotient.

::: details
For the computed quotient $q$, the radius is
$(|x - qy| + r_x + |q| r_y)/(|y| - r_y)$, with the residual $x - qy$
enclosed by the tracked kernels, so that an exact quotient of exact
balls is exact. Heads outside $[2^{-150}, 2^{150}]$ are handled by
scaling (`_dfloat_ball_div`).
:::
:::

::: function
int d1b_sqrt(d1b_t res, d1b_srcptr x)
int d2b_sqrt(d2b_t res, d2b_srcptr x)
int d3b_sqrt(d3b_t res, d3b_srcptr x)
int d4b_sqrt(d4b_t res, d4b_srcptr x)

Square root. Returns `GR_DOMAIN` if $x$ is entirely negative,
`GR_UNABLE` if it contains negative numbers (and positive ones or zero),
and `GR_SUCCESS` otherwise; $\sqrt{[0 \pm 0]} = [0 \pm 0]$.
:::

::: function
int d1b_rsqrt(d1b_t res, d1b_srcptr x)
int d2b_rsqrt(d2b_t res, d2b_srcptr x)
int d3b_rsqrt(d3b_t res, d3b_srcptr x)
int d4b_rsqrt(d4b_t res, d4b_srcptr x)

Reciprocal square root. Returns `GR_DOMAIN` if $x$ is entirely
nonpositive, `GR_UNABLE` if it contains zero, and `GR_SUCCESS`
otherwise.

::: details
For the computed $r$ and $e' = 1 - x' r^2$, $x'^{-1/2} = r (1 - e')^{-1/2}$,
so that $|x'^{-1/2} - r| \le |r| (|e'|/2)(1 + |e'|)$ for
$|e'| \le 1/4$, where $|e'| \le |1 - x r^2| + r_x r^2$ and $1 - x r^2$
is enclosed by the tracked kernels (two products and a subtraction,
which makes the ball version relatively more expensive than the ball
square root, whose residual needs one product).
:::
:::

::: function
void d1_canonicalise(d1_t res, d1_srcptr x)
void d2_canonicalise(d2_t res, d2_srcptr x)
void d3_canonicalise(d3_t res, d3_srcptr x)
void d4_canonicalise(d4_t res, d4_srcptr x)
void d1b_canonicalise(d1b_t res, d1b_srcptr x)
void d2b_canonicalise(d2b_t res, d2b_srcptr x)
void d3b_canonicalise(d3b_t res, d3b_srcptr x)
void d4b_canonicalise(d4b_t res, d4b_srcptr x)

Sets `res` to $x$ (exactly, with the radius unchanged) with the
components in canonical form: $|x_{i+1}| \le \operatorname{ulp}(x_i)/2$
and zeros only at the end, which is the form obtained by rounding the
value to nearest term by term. Works for any $x$ whose components
overlap by at most a few bits, in particular any result of this module.
A nonfinite head is copied as it is by `dN_canonicalise`;
`dNb_canonicalise` gives the whole real line in its normal form.

::: details
The arithmetic renormalises with a single bottom-up pass (*weak*
renormalisation), which is all the error analysis needs, but after a
cancellation a component can be smaller than the one following it, so
the representation is not unique and `dN_equal` compares values rather
than components. For $N = 2$ the weak form is already canonical; for
$N = 3, 4$ about 4 to 7% of the results of random operations are not.
Canonicalisation is a top-down TwoSum chain with a check, and a branchy
repair pass in about $10^{-5}$ to $10^{-4}$ of the cases.
:::
:::

::: function
void d1_renorm(d1_t res, double * t, slong len)
void d2_renorm(d2_t res, double * t, slong len)
void d3_renorm(d3_t res, double * t, slong len)
void d4_renorm(d4_t res, double * t, slong len)

Sets `res` to the sum of the `len` doubles $t_0, \ldots, t_{len-1}$,
rounded to $N$ nonoverlapping components; the dropped part is
discarded. The terms may be in any order, but for an accurate result
they should be roughly sorted by decreasing magnitude. The array `t` is
used as scratch space and overwritten; pass a copy
if the terms are needed afterwards. `res` must not overlap `t`.
:::

::: function
void d1b_add_error_d(d1b_t res, double err)
void d2b_add_error_d(d2b_t res, double err)
void d3b_add_error_d(d3b_t res, double err)
void d4b_add_error_d(d4b_t res, double err)

Adds $|\mathit{err}|$ to the radius of `res`, rounding up; a nonfinite
`err` gives the whole real line.
:::

## Elementary functions {#elementary}

The elementary functions of the plain formats are floating-point kernels
with static error bounds, proved once for all arguments in a range; the
observed errors on random arguments are 4 to 100 times smaller. The ball versions apply the kernel to the
midpoint, add the static bound, and propagate the input radius with a
bound on the derivative. The algorithms (table-based argument
reduction, polynomials evaluated in several chains with tapered
precision) are described in
[](algorithms/floating-point/dfloat-elementary). The functions are not
correctly rounded; for $N = 1$ they are about as fast and as accurate
as the C library functions.

::: summary algorithms/floating-point/dfloat-elementary :::

::: macro
DFLOAT_EXP_EPS_1
DFLOAT_EXP_EPS_2
DFLOAT_EXP_EPS_3
DFLOAT_EXP_EPS_4
DFLOAT_EXP_EPS(n)

Relative error bounds of the exponential kernels, about $2^{-51.0}$,
$2^{-100.2}$, $2^{-151.0}$ and $2^{-201.6}$: for
$-746 \le x_0 \le 709.9$, `dN_exp` returns $\exp(x)(1 + \delta)$ with
$|\delta| \le$ `DFLOAT_EXP_EPS_N`, up to an additional absolute error of
$2^{-1070}$ when the result is subnormal.
:::

::: function
void d1_exp(d1_t res, d1_srcptr x)
void d2_exp(d2_t res, d2_srcptr x)
void d3_exp(d3_t res, d3_srcptr x)
void d4_exp(d4_t res, d4_srcptr x)

The exponential function, with the error bound `DFLOAT_EXP_EPS_N`.
Nonfinite arguments and arguments beyond the range give the double
result (infinity, 0 or NaN); `dN_exp` of 0 is exactly 1.
:::

::: function
void d1b_exp(d1b_t res, d1b_srcptr x)
void d2b_exp(d2b_t res, d2b_srcptr x)
void d3b_exp(d3b_t res, d3b_srcptr x)
void d4b_exp(d4b_t res, d4b_srcptr x)

The exponential of a ball. The radius is
$|y| (1 + 2\varepsilon)((e^r - 1) + 2\varepsilon)$ for the kernel value
$y$ at the midpoint, the input radius $r$ and the static bound
$\varepsilon$. Where this overflows for a wide ball whose upper end $h$
does not, the result is $[0 \pm e^h]$. Arguments beyond the double
range give the whole real line (overflow) or a ball of radius at most
$2^{-1074}$ around zero (underflow); $\exp([0 \pm 0]) = [1 \pm 0]$.
:::

::: macro
DFLOAT_TRIG_EPS_1
DFLOAT_TRIG_EPS_2
DFLOAT_TRIG_EPS_3
DFLOAT_TRIG_EPS_4
DFLOAT_TRIG_EPS(n)

Error bounds of the sine and cosine kernels, about $2^{-50.4}$,
$2^{-98.5}$, $2^{-148.5}$ and $2^{-200.1}$: for
$|x| \le 2^{20} \pi/2 - 1$, `dN_sin` has an absolute error of at most
$\varepsilon (|x| + \min(|x|, 1))$ and `dN_cos` of at most
$\varepsilon (|x| + 1)$, with $\varepsilon$ = `DFLOAT_TRIG_EPS_N`.
The bounds are absolute and scaled by $|x|$: the relative accuracy near
the zeros of the functions is not preserved, since the argument is
reduced with $\pi/2$ to about $53N + 48$ bits and the reduced argument
keeps the absolute accuracy $u^N |x|$ of the input. For $|x| \le 1$ the
sine bound $2\varepsilon|x|$ is relative.
:::

::: function
void d1_sin(d1_t res, d1_srcptr x)
void d2_sin(d2_t res, d2_srcptr x)
void d3_sin(d3_t res, d3_srcptr x)
void d4_sin(d4_t res, d4_srcptr x)
void d1_cos(d1_t res, d1_srcptr x)
void d2_cos(d2_t res, d2_srcptr x)
void d3_cos(d3_t res, d3_srcptr x)
void d4_cos(d4_t res, d4_srcptr x)
void d1_sin_cos(d1_t sn, d1_t cs, d1_srcptr x)
void d2_sin_cos(d2_t sn, d2_t cs, d2_srcptr x)
void d3_sin_cos(d3_t sn, d3_t cs, d3_srcptr x)
void d4_sin_cos(d4_t sn, d4_t cs, d4_srcptr x)

The sine and cosine, with the bounds `DFLOAT_TRIG_EPS_N`. `dN_sin_cos`
computes both from one argument reduction, at about 10% more than one
of them, with results identical to those of `dN_sin` and `dN_cos`.
Arguments beyond $2^{20}\pi/2 - 1$ in absolute value are evaluated
through `arb_t` with the same guarantee; nonfinite arguments give NaN.
:::

::: function
void d1b_sin(d1b_t res, d1b_srcptr x)
void d2b_sin(d2b_t res, d2b_srcptr x)
void d3b_sin(d3b_t res, d3b_srcptr x)
void d4b_sin(d4b_t res, d4b_srcptr x)
void d1b_cos(d1b_t res, d1b_srcptr x)
void d2b_cos(d2b_t res, d2b_srcptr x)
void d3b_cos(d3b_t res, d3b_srcptr x)
void d4b_cos(d4b_t res, d4b_srcptr x)
void d1b_sin_cos(d1b_t sn, d1b_t cs, d1b_srcptr x)
void d2b_sin_cos(d2b_t sn, d2b_t cs, d2b_srcptr x)
void d3b_sin_cos(d3b_t sn, d3b_t cs, d3b_srcptr x)
void d4b_sin_cos(d4b_t sn, d4b_t cs, d4b_srcptr x)

The sine and cosine of a ball. The radius is $\varepsilon(|m| + a) + r$
with $a = \min(|m|, 1)$ for the sine and $a = 1$ for the cosine, for
the midpoint $m$ and the input radius $r$ (both functions are
1-Lipschitz); a radius of 1 or more is replaced by the enclosure
$[0 \pm 1]$. Infinite radii, nonfinite midpoints and midpoints beyond
the kernel range go through `arb_t`. $\sin([0 \pm 0]) = [0 \pm 0]$
and $\cos([0 \pm 0]) = [1 \pm 0]$.
:::

::: macro
DFLOAT_LOG_EPS_1
DFLOAT_LOG_EPS_2
DFLOAT_LOG_EPS_3
DFLOAT_LOG_EPS_4
DFLOAT_LOG_EPS(n)
DFLOAT_ATAN_EPS_1
DFLOAT_ATAN_EPS_2
DFLOAT_ATAN_EPS_3
DFLOAT_ATAN_EPS_4
DFLOAT_ATAN_EPS(n)

Relative error bounds of the logarithm and arctangent kernels:
$|\mathtt{dN\_log}(x) - \log x| \le \varepsilon |\log x|$ for
$2^{-1022} \le x < 2^{1022}$, with $\varepsilon$ about $2^{-49.3}$,
$2^{-99.3}$, $2^{-150.2}$ and $2^{-200.8}$, and
$|\mathtt{dN\_atan}(x) - \operatorname{atan} x| \le \varepsilon |\operatorname{atan} x|$
for all finite $x$, with $\varepsilon$ about $2^{-49.8}$, $2^{-99.6}$,
$2^{-150.5}$ and $2^{-201.3}$. Both bounds are relative everywhere,
including near $x = 1$ for the logarithm and near 0 for the arctangent.
:::

::: function
void d1_log(d1_t res, d1_srcptr x)
void d2_log(d2_t res, d2_srcptr x)
void d3_log(d3_t res, d3_srcptr x)
void d4_log(d4_t res, d4_srcptr x)
void d1_atan(d1_t res, d1_srcptr x)
void d2_atan(d2_t res, d2_srcptr x)
void d3_atan(d3_t res, d3_srcptr x)
void d4_atan(d4_t res, d4_srcptr x)

The natural logarithm and the arctangent, with the bounds above. The
logarithm of a nonpositive number is $-\infty$ (for zero) or NaN, of
$+\infty$ it is $+\infty$, and $\log 1 = 0$ exactly; arguments outside
$[2^{-1022}, 2^{1022})$ go through `arb_t`. The arctangent of
$\pm\infty$ is $\pm \pi/2$, and $\operatorname{atan} 0 = 0$ exactly.
:::

::: function
int d1b_log(d1b_t res, d1b_srcptr x)
int d2b_log(d2b_t res, d2b_srcptr x)
int d3b_log(d3b_t res, d3b_srcptr x)
int d4b_log(d4b_t res, d4b_srcptr x)
void d1b_atan(d1b_t res, d1b_srcptr x)
void d2b_atan(d2b_t res, d2b_srcptr x)
void d3b_atan(d3b_t res, d3b_srcptr x)
void d4b_atan(d4b_t res, d4b_srcptr x)

The logarithm and arctangent of a ball $[m \pm r]$, with the radii
$\varepsilon |y| + r/(m - r)$ and $\varepsilon |y| + r$. `dNb_log`
returns `GR_DOMAIN` if the ball is entirely nonpositive, `GR_UNABLE` if
it is otherwise not contained in $(0, \infty)$, and `GR_SUCCESS`
otherwise. Infinite radii, nonfinite midpoints and midpoints beyond the
kernel ranges go through `arb_t`.
:::

::: macro
DFLOAT_EXPM1_EPS_1
DFLOAT_EXPM1_EPS_2
DFLOAT_EXPM1_EPS_3
DFLOAT_EXPM1_EPS_4
DFLOAT_EXPM1_EPS(n)
DFLOAT_LOG1P_EPS_1
DFLOAT_LOG1P_EPS_2
DFLOAT_LOG1P_EPS_3
DFLOAT_LOG1P_EPS_4
DFLOAT_LOG1P_EPS(n)

Relative error bounds of `dN_expm1`, valid for $-746 \le x \le 709.9$
(about $2^{-48.9}$, $2^{-98.3}$, $2^{-149.2}$, $2^{-199.8}$), and of
`dN_log1p`, valid for $1 + x \in [2^{-1022}, 2^{1022})$ (about
$2^{-49.3}$, $2^{-99.1}$, $2^{-149.9}$, $2^{-200.5}$). They are relative
everywhere, in particular near 0, where both functions are
$x(1 + O(x))$.
:::

::: function
void d1_expm1(d1_t res, d1_srcptr x)
void d2_expm1(d2_t res, d2_srcptr x)
void d3_expm1(d3_t res, d3_srcptr x)
void d4_expm1(d4_t res, d4_srcptr x)
void d1_log1p(d1_t res, d1_srcptr x)
void d2_log1p(d2_t res, d2_srcptr x)
void d3_log1p(d3_t res, d3_srcptr x)
void d4_log1p(d4_t res, d4_srcptr x)

$e^x - 1$ and $\log(1 + x)$, with the bounds above, accurate also for
small $|x|$ where `dN_exp` and `dN_log` would cancel.
:::

::: function
void d1b_expm1(d1b_t res, d1b_srcptr x)
void d2b_expm1(d2b_t res, d2b_srcptr x)
void d3b_expm1(d3b_t res, d3b_srcptr x)
void d4b_expm1(d4b_t res, d4b_srcptr x)
int d1b_log1p(d1b_t res, d1b_srcptr x)
int d2b_log1p(d2b_t res, d2b_srcptr x)
int d3b_log1p(d3b_t res, d3b_srcptr x)
int d4b_log1p(d4b_t res, d4b_srcptr x)

$e^x - 1$ and $\log(1 + x)$ for balls, propagating the input radius as
`dNb_exp` and `dNb_log` do. `dNb_log1p` returns `GR_DOMAIN` if the ball
is entirely in $(-\infty, -1]$, `GR_UNABLE` if it is otherwise not
contained in $(-1, \infty)$, and `GR_SUCCESS` otherwise.
:::

::: function
void d1_tan(d1_t res, d1_srcptr x)
void d2_tan(d2_t res, d2_srcptr x)
void d3_tan(d3_t res, d3_srcptr x)
void d4_tan(d4_t res, d4_srcptr x)
void d1_sinh(d1_t res, d1_srcptr x)
void d2_sinh(d2_t res, d2_srcptr x)
void d3_sinh(d3_t res, d3_srcptr x)
void d4_sinh(d4_t res, d4_srcptr x)
void d1_cosh(d1_t res, d1_srcptr x)
void d2_cosh(d2_t res, d2_srcptr x)
void d3_cosh(d3_t res, d3_srcptr x)
void d4_cosh(d4_t res, d4_srcptr x)
void d1_tanh(d1_t res, d1_srcptr x)
void d2_tanh(d2_t res, d2_srcptr x)
void d3_tanh(d3_t res, d3_srcptr x)
void d4_tanh(d4_t res, d4_srcptr x)
void d1_asin(d1_t res, d1_srcptr x)
void d2_asin(d2_t res, d2_srcptr x)
void d3_asin(d3_t res, d3_srcptr x)
void d4_asin(d4_t res, d4_srcptr x)
void d1_acos(d1_t res, d1_srcptr x)
void d2_acos(d2_t res, d2_srcptr x)
void d3_acos(d3_t res, d3_srcptr x)
void d4_acos(d4_t res, d4_srcptr x)
void d1_atan2(d1_t res, d1_srcptr y, d1_srcptr x)
void d2_atan2(d2_t res, d2_srcptr y, d2_srcptr x)
void d3_atan2(d3_t res, d3_srcptr y, d3_srcptr x)
void d4_atan2(d4_t res, d4_srcptr y, d4_srcptr x)
void d1_pow(d1_t res, d1_srcptr x, d1_srcptr y)
void d2_pow(d2_t res, d2_srcptr x, d2_srcptr y)
void d3_pow(d3_t res, d3_srcptr x, d3_srcptr y)
void d4_pow(d4_t res, d4_srcptr x, d4_srcptr y)

Further elementary functions, as compositions of the kernels with the
arithmetic (the arguments are canonicalised first). Their relative
errors are about the sums of those of their parts (about $2^{-49}$,
$2^{-98}$, $2^{-149}$, $2^{-199}$), except that `dN_pow` amplifies the
error of the logarithm by $|y \log x|$ (up to $2^{9.5}$ at the end of
the range of the exponential) and that `tan`, `sinh` and `tanh` inherit
absolute errors near their zeros. `dN_atan2(res, y, x)` follows the
double `atan2` for nonfinite arguments, and `dN_pow` gives NaN for
$x < 0$ and noninteger $y$.

::: details
The formulas are $\tan x = \sin x / \cos x$;
$\sinh |x| = e(e + 2)/(2(e + 1))$ with $e = \operatorname{expm1}|x|$
for $|x| < 1$ and $(E - 1/E)/2$ with $E = e^{|x|}$ beyond;
$\cosh x = (E + 1/E)/2$; $\tanh |x| = e/(e + 2)$ with
$e = \operatorname{expm1}(2|x|)$ for $|x| < 1$ and $1 - 2/(E + 1)$ with
$E = e^{2|x|}$ beyond; $\operatorname{asin} x = \operatorname{atan}(x/\sqrt{(1 - x)(1 + x)})$
and $\operatorname{acos} x = 2 \operatorname{atan}\sqrt{(1 - x)/(1 + x)}$
(where $1 \pm x$ cancel they are exact, so both stay relative up to
$\pm 1$); $\operatorname{atan2}(y, x) = \operatorname{atan}(y/x)$
adjusted by $\pm\pi$ for $x < 0$ and $\pm\pi/2$ for $x = 0$; and
$x^y = \exp(y \log x)$ for $x > 0$, binary powering of $x$ or $1/x$ for
an integer $y$ with $|y| < 2^{40}$ (any base).
:::
:::

::: function
int d1b_tan(d1b_t res, d1b_srcptr x)
int d2b_tan(d2b_t res, d2b_srcptr x)
int d3b_tan(d3b_t res, d3b_srcptr x)
int d4b_tan(d4b_t res, d4b_srcptr x)
void d1b_sinh(d1b_t res, d1b_srcptr x)
void d2b_sinh(d2b_t res, d2b_srcptr x)
void d3b_sinh(d3b_t res, d3b_srcptr x)
void d4b_sinh(d4b_t res, d4b_srcptr x)
void d1b_cosh(d1b_t res, d1b_srcptr x)
void d2b_cosh(d2b_t res, d2b_srcptr x)
void d3b_cosh(d3b_t res, d3b_srcptr x)
void d4b_cosh(d4b_t res, d4b_srcptr x)
void d1b_tanh(d1b_t res, d1b_srcptr x)
void d2b_tanh(d2b_t res, d2b_srcptr x)
void d3b_tanh(d3b_t res, d3b_srcptr x)
void d4b_tanh(d4b_t res, d4b_srcptr x)
int d1b_asin(d1b_t res, d1b_srcptr x)
int d2b_asin(d2b_t res, d2b_srcptr x)
int d3b_asin(d3b_t res, d3b_srcptr x)
int d4b_asin(d4b_t res, d4b_srcptr x)
int d1b_acos(d1b_t res, d1b_srcptr x)
int d2b_acos(d2b_t res, d2b_srcptr x)
int d3b_acos(d3b_t res, d3b_srcptr x)
int d4b_acos(d4b_t res, d4b_srcptr x)
void d1b_atan2(d1b_t res, d1b_srcptr y, d1b_srcptr x)
void d2b_atan2(d2b_t res, d2b_srcptr y, d2b_srcptr x)
void d3b_atan2(d3b_t res, d3b_srcptr y, d3b_srcptr x)
void d4b_atan2(d4b_t res, d4b_srcptr y, d4b_srcptr x)
int d1b_pow(d1b_t res, d1b_srcptr x, d1b_srcptr y)
int d2b_pow(d2b_t res, d2b_srcptr x, d2b_srcptr y)
int d3b_pow(d3b_t res, d3b_srcptr x, d3b_srcptr y)
int d4b_pow(d4b_t res, d4b_srcptr x, d4b_srcptr y)

The same functions on balls, as the same compositions of the ball
functions. `dNb_tan`, `dNb_asin`, `dNb_acos` and `dNb_pow` return a
status: `GR_SUCCESS` if the input is inside the domain, `GR_DOMAIN` if
it is outside, and `GR_UNABLE` otherwise. The domains are
$x \notin \pi/2 + \pi\Z$ for the tangent (never `GR_DOMAIN`, since the
poles are irrational), $|x| \le 1$ for `asin` and `acos`, and for
`pow`, $x > 0$, or $x = 0$ and $y \ge 0$, or $x < 0$ and $y \in \Z$.

::: details
Where the composition fails for a ball inside the domain, `dNb_asin`
and `dNb_acos` (at and near $\pm 1$, where $1 - x^2$ contains zero)
evaluate the monotone functions at the exact endpoints through `arb_t`,
and `dNb_tan` (near a pole, where the absolute error of the cosine
covers zero) calls `arb_tan`; `dNb_pow` with an exact integer exponent
of $2^{40}$ or more uses `arb_pow_fmpz`. `dNb_atan2` always returns a
finite enclosure: $\operatorname{atan}(y/x)$ adjusted by $\pm\pi$ for
$x$ away from zero, $\pi$ when $x < 0$ and $y = 0$ exactly,
$[0 \pm \pi]$ when $x < 0$ and $y$ otherwise contains zero (the branch
cut), $\pm\pi/2 - \operatorname{atan}(x/y)$ for $x$ containing zero and
$y$ away from it, and $[0 \pm \pi]$ when both contain zero.
:::
:::

## Predicates and comparison {#predicates}

::: function
int d1_is_zero(d1_srcptr x)
int d2_is_zero(d2_srcptr x)
int d3_is_zero(d3_srcptr x)
int d4_is_zero(d4_srcptr x)
int d1_is_one(d1_srcptr x)
int d2_is_one(d2_srcptr x)
int d3_is_one(d3_srcptr x)
int d4_is_one(d4_srcptr x)
int d1_is_finite(d1_srcptr x)
int d2_is_finite(d2_srcptr x)
int d3_is_finite(d3_srcptr x)
int d4_is_finite(d4_srcptr x)
int d1_is_nan(d1_srcptr x)
int d2_is_nan(d2_srcptr x)
int d3_is_nan(d3_srcptr x)
int d4_is_nan(d4_srcptr x)
int d1_equal(d1_srcptr x, d1_srcptr y)
int d2_equal(d2_srcptr x, d2_srcptr y)
int d3_equal(d3_srcptr x, d3_srcptr y)
int d4_equal(d4_srcptr x, d4_srcptr y)
int d1_cmp(d1_srcptr x, d1_srcptr y)
int d2_cmp(d2_srcptr x, d2_srcptr y)
int d3_cmp(d3_srcptr x, d3_srcptr y)
int d4_cmp(d4_srcptr x, d4_srcptr y)
int d1_sgn(d1_srcptr x)
int d2_sgn(d2_srcptr x)
int d3_sgn(d3_srcptr x)
int d4_sgn(d4_srcptr x)

Predicates and comparisons of the exact values of plain expansions:
`dN_equal` compares values, not components, and `dN_cmp` and `dN_sgn`
return $-1$, 0 or 1. `dN_is_finite` is true if every component is
finite.
:::

::: function
int d1b_is_exact(d1b_srcptr x)
int d2b_is_exact(d2b_srcptr x)
int d3b_is_exact(d3b_srcptr x)
int d4b_is_exact(d4b_srcptr x)
int d1b_is_finite(d1b_srcptr x)
int d2b_is_finite(d2b_srcptr x)
int d3b_is_finite(d3b_srcptr x)
int d4b_is_finite(d4b_srcptr x)
int d1b_is_zero(d1b_srcptr x)
int d2b_is_zero(d2b_srcptr x)
int d3b_is_zero(d3b_srcptr x)
int d4b_is_zero(d4b_srcptr x)
int d1b_contains_zero(d1b_srcptr x)
int d2b_contains_zero(d2b_srcptr x)
int d3b_contains_zero(d3b_srcptr x)
int d4b_contains_zero(d4b_srcptr x)
int d1b_contains_d(d1b_srcptr x, double c)
int d2b_contains_d(d2b_srcptr x, double c)
int d3b_contains_d(d3b_srcptr x, double c)
int d4b_contains_d(d4b_srcptr x, double c)
int d1b_contains(d1b_srcptr x, d1b_srcptr y)
int d2b_contains(d2b_srcptr x, d2b_srcptr y)
int d3b_contains(d3b_srcptr x, d3b_srcptr y)
int d4b_contains(d4b_srcptr x, d4b_srcptr y)
int d1b_overlaps(d1b_srcptr x, d1b_srcptr y)
int d2b_overlaps(d2b_srcptr x, d2b_srcptr y)
int d3b_overlaps(d3b_srcptr x, d3b_srcptr y)
int d4b_overlaps(d4b_srcptr x, d4b_srcptr y)
int d1b_equal(d1b_srcptr x, d1b_srcptr y)
int d2b_equal(d2b_srcptr x, d2b_srcptr y)
int d3b_equal(d3b_srcptr x, d3b_srcptr y)
int d4b_equal(d4b_srcptr x, d4b_srcptr y)
int d1b_cmp(int * res, d1b_srcptr x, d1b_srcptr y)
int d2b_cmp(int * res, d2b_srcptr x, d2b_srcptr y)
int d3b_cmp(int * res, d3b_srcptr x, d3b_srcptr y)
int d4b_cmp(int * res, d4b_srcptr x, d4b_srcptr y)

Ball predicates. `is_exact` is true for a finite ball of radius zero
and `is_finite` for a ball that is not the whole real line. `is_zero`,
`equal`, `contains` and `contains_d` answer 1 only when certain
(`is_zero` and `equal` only for exact balls, `contains_d` only for a
finite $c$); `contains_zero` and `overlaps` answer 0 only when certain.
The whole real line contains and overlaps everything, and is contained
only in the whole real line. `dNb_cmp` sets `*res` to $-1$, 0 or 1 and
returns `GR_SUCCESS` when every pair of points of the balls compares
the same way (or both are the same exact number), and sets `*res` to 0
and returns `GR_UNABLE` otherwise.
:::

## Conversions, input and output {#conversions}

::: function
double d1_get_d(d1_srcptr x)
double d2_get_d(d2_srcptr x)
double d3_get_d(d3_srcptr x)
double d4_get_d(d4_srcptr x)
double d1b_get_d(d1b_srcptr x)
double d2b_get_d(d2b_srcptr x)
double d3b_get_d(d3b_srcptr x)
double d4b_get_d(d4b_srcptr x)

Return the head $d_0$ of the midpoint, the double nearest to it (up to
ties); 0 for the whole real line.
:::

::: function
void d1b_get_mid_d(double * mid, d1b_srcptr x)
void d2b_get_mid_d(double * mid, d2b_srcptr x)
void d3b_get_mid_d(double * mid, d3b_srcptr x)
void d4b_get_mid_d(double * mid, d4b_srcptr x)

Sets `mid[0]`, ..., `mid[N-1]` to the components of the midpoint (zeros
for the whole real line).
:::

::: function
void d1_get_arf(arf_t res, d1_srcptr x)
void d2_get_arf(arf_t res, d2_srcptr x)
void d3_get_arf(arf_t res, d3_srcptr x)
void d4_get_arf(arf_t res, d4_srcptr x)
void d1_set_arf(d1_t res, const arf_t x)
void d2_set_arf(d2_t res, const arf_t x)
void d3_set_arf(d3_t res, const arf_t x)
void d4_set_arf(d4_t res, const arf_t x)
void d1_get_arb(arb_t res, d1_srcptr x)
void d2_get_arb(arb_t res, d2_srcptr x)
void d3_get_arb(arb_t res, d3_srcptr x)
void d4_get_arb(arb_t res, d4_srcptr x)
void d1_set_arb(d1_t res, const arb_t x)
void d2_set_arb(d2_t res, const arb_t x)
void d3_set_arb(d3_t res, const arb_t x)
void d4_set_arb(d4_t res, const arb_t x)
void d1b_get_arb(arb_t res, d1b_srcptr x)
void d2b_get_arb(arb_t res, d2b_srcptr x)
void d3b_get_arb(arb_t res, d3b_srcptr x)
void d4b_get_arb(arb_t res, d4b_srcptr x)
void d1b_set_arb(d1b_t res, const arb_t x)
void d2b_set_arb(d2b_t res, const arb_t x)
void d3b_set_arb(d3b_t res, const arb_t x)
void d4b_set_arb(d4b_t res, const arb_t x)
void d1b_set_arf(d1b_t res, const arf_t x)
void d2b_set_arf(d2b_t res, const arf_t x)
void d3b_set_arf(d3b_t res, const arf_t x)
void d4b_set_arf(d4b_t res, const arf_t x)

Conversions from and to `arf_t` and `arb_t`. The `get` functions are
exact (the radius is rounded up to a `mag_t`). The `set` functions
round the midpoint to $N$ doubles by greedy rounding and, for balls,
add the remainder to the radius. A nonfinite `arb_t` or `arf_t` gives
the whole real line, and `dNb_get_arb` gives the whole real line as
$[0 \pm \infty]$.
:::

::: function
void d1_set_fmpz(d1_t res, const fmpz_t x)
void d2_set_fmpz(d2_t res, const fmpz_t x)
void d3_set_fmpz(d3_t res, const fmpz_t x)
void d4_set_fmpz(d4_t res, const fmpz_t x)
void d1_set_fmpq(d1_t res, const fmpq_t x)
void d2_set_fmpq(d2_t res, const fmpq_t x)
void d3_set_fmpq(d3_t res, const fmpq_t x)
void d4_set_fmpq(d4_t res, const fmpq_t x)
void d1b_set_fmpz(d1b_t res, const fmpz_t x)
void d2b_set_fmpz(d2b_t res, const fmpz_t x)
void d3b_set_fmpz(d3b_t res, const fmpz_t x)
void d4b_set_fmpz(d4b_t res, const fmpz_t x)
void d1b_set_fmpq(d1b_t res, const fmpq_t x)
void d2b_set_fmpq(d2b_t res, const fmpq_t x)
void d3b_set_fmpq(d3b_t res, const fmpq_t x)
void d4b_set_fmpq(d4b_t res, const fmpq_t x)

Set `res` to an integer or a rational number, rounded to $N$ doubles
through an `arb_t` (rigorously, for the balls).
:::

::: function
int d1_set_str(d1_t res, const char * s)
int d2_set_str(d2_t res, const char * s)
int d3_set_str(d3_t res, const char * s)
int d4_set_str(d4_t res, const char * s)
int d1b_set_str(d1b_t res, const char * s)
int d2b_set_str(d2b_t res, const char * s)
int d3b_set_str(d3b_t res, const char * s)
int d4b_set_str(d4b_t res, const char * s)

Parse a decimal string with `arb_set_str` (at $53N + 64$ bits) and
convert, rigorously for the balls. Return `GR_SUCCESS`, or `GR_UNABLE`
if the string could not be parsed.
:::

::: function
char * d1_get_str(d1_srcptr x, slong digits)
char * d2_get_str(d2_srcptr x, slong digits)
char * d3_get_str(d3_srcptr x, slong digits)
char * d4_get_str(d4_srcptr x, slong digits)
void d1_print(d1_srcptr x)
void d2_print(d2_srcptr x)
void d3_print(d3_srcptr x)
void d4_print(d4_srcptr x)
char * d1b_get_str(d1b_srcptr x, slong digits)
char * d2b_get_str(d2b_srcptr x, slong digits)
char * d3b_get_str(d3b_srcptr x, slong digits)
char * d4b_get_str(d4b_srcptr x, slong digits)
void d1b_print(d1b_srcptr x)
void d2b_print(d2b_srcptr x)
void d3b_print(d3b_srcptr x)
void d4b_print(d4b_srcptr x)

Convert to a decimal string with `digits` significant digits through
`arb_get_str` (`digits` $\le 0$ selects $\lfloor 53N \log_{10} 2 \rfloor + 1$
digits), or print to standard output with the default number of digits.
Balls are printed with their radius, and every form of the whole real
line as `[+/- inf]`. The string must be freed with `flint_free`.
:::

::: function
void d1_randtest(d1_t res, flint_rand_t state)
void d2_randtest(d2_t res, flint_rand_t state)
void d3_randtest(d3_t res, flint_rand_t state)
void d4_randtest(d4_t res, flint_rand_t state)
void d1_randtest_special(d1_t res, flint_rand_t state)
void d2_randtest_special(d2_t res, flint_rand_t state)
void d3_randtest_special(d3_t res, flint_rand_t state)
void d4_randtest_special(d4_t res, flint_rand_t state)
void d1b_randtest(d1b_t res, flint_rand_t state)
void d2b_randtest(d2b_t res, flint_rand_t state)
void d3b_randtest(d3b_t res, flint_rand_t state)
void d4b_randtest(d4b_t res, flint_rand_t state)
void d1b_randtest_special(d1b_t res, flint_rand_t state)
void d2b_randtest_special(d2b_t res, flint_rand_t state)
void d3b_randtest_special(d3b_t res, flint_rand_t state)
void d4b_randtest_special(d4b_t res, flint_rand_t state)

Random values for testing; the `special` versions also produce
subnormal numbers, values at the ends of the double range, nonfinite
values and radii from zero to infinity.
:::

## Vector operations {#vectors}

These functions take vectors of `len` elements and are used by the
generic rings for `_gr_vec_add`, `_gr_vec_dot` and so on. With AVX2 and
FMA on x86-64, or NEON on AArch64, they process four elements at a time
with the same results as the scalar functions; elsewhere they are loops
over the scalar code.

::: function
void _d1_vec_add(d1_ptr res, d1_srcptr x, d1_srcptr y, slong len)
void _d2_vec_add(d2_ptr res, d2_srcptr x, d2_srcptr y, slong len)
void _d3_vec_add(d3_ptr res, d3_srcptr x, d3_srcptr y, slong len)
void _d4_vec_add(d4_ptr res, d4_srcptr x, d4_srcptr y, slong len)
void _d1_vec_sub(d1_ptr res, d1_srcptr x, d1_srcptr y, slong len)
void _d2_vec_sub(d2_ptr res, d2_srcptr x, d2_srcptr y, slong len)
void _d3_vec_sub(d3_ptr res, d3_srcptr x, d3_srcptr y, slong len)
void _d4_vec_sub(d4_ptr res, d4_srcptr x, d4_srcptr y, slong len)
void _d1_vec_mul(d1_ptr res, d1_srcptr x, d1_srcptr y, slong len)
void _d2_vec_mul(d2_ptr res, d2_srcptr x, d2_srcptr y, slong len)
void _d3_vec_mul(d3_ptr res, d3_srcptr x, d3_srcptr y, slong len)
void _d4_vec_mul(d4_ptr res, d4_srcptr x, d4_srcptr y, slong len)
void _d1_vec_mul_scalar(d1_ptr res, d1_srcptr x, slong len, d1_srcptr c)
void _d2_vec_mul_scalar(d2_ptr res, d2_srcptr x, slong len, d2_srcptr c)
void _d3_vec_mul_scalar(d3_ptr res, d3_srcptr x, slong len, d3_srcptr c)
void _d4_vec_mul_scalar(d4_ptr res, d4_srcptr x, slong len, d4_srcptr c)
void _d1_vec_addmul_scalar(d1_ptr res, d1_srcptr x, slong len, d1_srcptr c)
void _d2_vec_addmul_scalar(d2_ptr res, d2_srcptr x, slong len, d2_srcptr c)
void _d3_vec_addmul_scalar(d3_ptr res, d3_srcptr x, slong len, d3_srcptr c)
void _d4_vec_addmul_scalar(d4_ptr res, d4_srcptr x, slong len, d4_srcptr c)
void _d1_vec_submul_scalar(d1_ptr res, d1_srcptr x, slong len, d1_srcptr c)
void _d2_vec_submul_scalar(d2_ptr res, d2_srcptr x, slong len, d2_srcptr c)
void _d3_vec_submul_scalar(d3_ptr res, d3_srcptr x, slong len, d3_srcptr c)
void _d4_vec_submul_scalar(d4_ptr res, d4_srcptr x, slong len, d4_srcptr c)
void _d1_vec_div(d1_ptr res, d1_srcptr x, d1_srcptr y, slong len)
void _d2_vec_div(d2_ptr res, d2_srcptr x, d2_srcptr y, slong len)
void _d3_vec_div(d3_ptr res, d3_srcptr x, d3_srcptr y, slong len)
void _d4_vec_div(d4_ptr res, d4_srcptr x, d4_srcptr y, slong len)
void _d1_vec_sqrt(d1_ptr res, d1_srcptr x, slong len)
void _d2_vec_sqrt(d2_ptr res, d2_srcptr x, slong len)
void _d3_vec_sqrt(d3_ptr res, d3_srcptr x, slong len)
void _d4_vec_sqrt(d4_ptr res, d4_srcptr x, slong len)
void _d1_vec_rsqrt(d1_ptr res, d1_srcptr x, slong len)
void _d2_vec_rsqrt(d2_ptr res, d2_srcptr x, slong len)
void _d3_vec_rsqrt(d3_ptr res, d3_srcptr x, slong len)
void _d4_vec_rsqrt(d4_ptr res, d4_srcptr x, slong len)

Elementwise sum, difference, product, product with the scalar $c$,
accumulation of $\pm c\, x_i$ into `res`, quotient, square root and
reciprocal square root.
:::

::: function
void _d1_vec_dot(d1_t res, d1_srcptr initial, int subtract, d1_srcptr x, d1_srcptr y, slong len)
void _d2_vec_dot(d2_t res, d2_srcptr initial, int subtract, d2_srcptr x, d2_srcptr y, slong len)
void _d3_vec_dot(d3_t res, d3_srcptr initial, int subtract, d3_srcptr x, d3_srcptr y, slong len)
void _d4_vec_dot(d4_t res, d4_srcptr initial, int subtract, d4_srcptr x, d4_srcptr y, slong len)
void _d1_vec_dot_rev(d1_t res, d1_srcptr initial, int subtract, d1_srcptr x, d1_srcptr y, slong len)
void _d2_vec_dot_rev(d2_t res, d2_srcptr initial, int subtract, d2_srcptr x, d2_srcptr y, slong len)
void _d3_vec_dot_rev(d3_t res, d3_srcptr initial, int subtract, d3_srcptr x, d3_srcptr y, slong len)
void _d4_vec_dot_rev(d4_t res, d4_srcptr initial, int subtract, d4_srcptr x, d4_srcptr y, slong len)

Set `res` to $s \pm \sum_i x_i y_i$ (`vec_dot`) or
$s \pm \sum_i x_i y_{len-1-i}$ (`vec_dot_rev`), where $s$ is `initial`
(or zero if `initial` is `NULL`) and the sign is negative if `subtract`
is nonzero, as `_gr_vec_dot` and `_gr_vec_dot_rev`.

::: details
The SIMD versions keep four lane accumulators and combine them at the
end; `vec_dot_rev` (the kernel of classical polynomial multiplication
and of the power series base cases) loads the blocks of $y$ in reverse
lane order, so that it gives the same result as `vec_dot` on a reversed
copy of $y$. Dot products shorter than 32 terms for `d1`, 10 for `d2`
and `d3`, and 8 for the balls use a lean scalar loop.
:::
:::

::: function
void _d1_vec_exp(d1_ptr res, d1_srcptr x, slong len)
void _d2_vec_exp(d2_ptr res, d2_srcptr x, slong len)
void _d3_vec_exp(d3_ptr res, d3_srcptr x, slong len)
void _d4_vec_exp(d4_ptr res, d4_srcptr x, slong len)
void _d1_vec_sin(d1_ptr res, d1_srcptr x, slong len)
void _d2_vec_sin(d2_ptr res, d2_srcptr x, slong len)
void _d3_vec_sin(d3_ptr res, d3_srcptr x, slong len)
void _d4_vec_sin(d4_ptr res, d4_srcptr x, slong len)
void _d1_vec_cos(d1_ptr res, d1_srcptr x, slong len)
void _d2_vec_cos(d2_ptr res, d2_srcptr x, slong len)
void _d3_vec_cos(d3_ptr res, d3_srcptr x, slong len)
void _d4_vec_cos(d4_ptr res, d4_srcptr x, slong len)
void _d1_vec_sin_cos(d1_ptr sres, d1_ptr cres, d1_srcptr x, slong len)
void _d2_vec_sin_cos(d2_ptr sres, d2_ptr cres, d2_srcptr x, slong len)
void _d3_vec_sin_cos(d3_ptr sres, d3_ptr cres, d3_srcptr x, slong len)
void _d4_vec_sin_cos(d4_ptr sres, d4_ptr cres, d4_srcptr x, slong len)
void _d1_vec_log(d1_ptr res, d1_srcptr x, slong len)
void _d2_vec_log(d2_ptr res, d2_srcptr x, slong len)
void _d3_vec_log(d3_ptr res, d3_srcptr x, slong len)
void _d4_vec_log(d4_ptr res, d4_srcptr x, slong len)
void _d1_vec_atan(d1_ptr res, d1_srcptr x, slong len)
void _d2_vec_atan(d2_ptr res, d2_srcptr x, slong len)
void _d3_vec_atan(d3_ptr res, d3_srcptr x, slong len)
void _d4_vec_atan(d4_ptr res, d4_srcptr x, slong len)

Elementwise elementary functions with the same static error bounds as
the scalar kernels, which they run on four arguments at a time (3 to 5
times faster per element than scalar calls). Elements that hit a
special case are redone with the scalar code.
:::

::: function
void _d1b_vec_add(d1b_ptr res, d1b_srcptr x, d1b_srcptr y, slong len)
void _d2b_vec_add(d2b_ptr res, d2b_srcptr x, d2b_srcptr y, slong len)
void _d3b_vec_add(d3b_ptr res, d3b_srcptr x, d3b_srcptr y, slong len)
void _d4b_vec_add(d4b_ptr res, d4b_srcptr x, d4b_srcptr y, slong len)
void _d1b_vec_sub(d1b_ptr res, d1b_srcptr x, d1b_srcptr y, slong len)
void _d2b_vec_sub(d2b_ptr res, d2b_srcptr x, d2b_srcptr y, slong len)
void _d3b_vec_sub(d3b_ptr res, d3b_srcptr x, d3b_srcptr y, slong len)
void _d4b_vec_sub(d4b_ptr res, d4b_srcptr x, d4b_srcptr y, slong len)
void _d1b_vec_mul(d1b_ptr res, d1b_srcptr x, d1b_srcptr y, slong len)
void _d2b_vec_mul(d2b_ptr res, d2b_srcptr x, d2b_srcptr y, slong len)
void _d3b_vec_mul(d3b_ptr res, d3b_srcptr x, d3b_srcptr y, slong len)
void _d4b_vec_mul(d4b_ptr res, d4b_srcptr x, d4b_srcptr y, slong len)
void _d1b_vec_mul_scalar(d1b_ptr res, d1b_srcptr x, slong len, d1b_srcptr c)
void _d2b_vec_mul_scalar(d2b_ptr res, d2b_srcptr x, slong len, d2b_srcptr c)
void _d3b_vec_mul_scalar(d3b_ptr res, d3b_srcptr x, slong len, d3b_srcptr c)
void _d4b_vec_mul_scalar(d4b_ptr res, d4b_srcptr x, slong len, d4b_srcptr c)
void _d1b_vec_addmul_scalar(d1b_ptr res, d1b_srcptr x, slong len, d1b_srcptr c)
void _d2b_vec_addmul_scalar(d2b_ptr res, d2b_srcptr x, slong len, d2b_srcptr c)
void _d3b_vec_addmul_scalar(d3b_ptr res, d3b_srcptr x, slong len, d3b_srcptr c)
void _d4b_vec_addmul_scalar(d4b_ptr res, d4b_srcptr x, slong len, d4b_srcptr c)
void _d1b_vec_submul_scalar(d1b_ptr res, d1b_srcptr x, slong len, d1b_srcptr c)
void _d2b_vec_submul_scalar(d2b_ptr res, d2b_srcptr x, slong len, d2b_srcptr c)
void _d3b_vec_submul_scalar(d3b_ptr res, d3b_srcptr x, slong len, d3b_srcptr c)
void _d4b_vec_submul_scalar(d4b_ptr res, d4b_srcptr x, slong len, d4b_srcptr c)
void _d1b_vec_div(d1b_ptr res, d1b_srcptr x, d1b_srcptr y, slong len)
void _d2b_vec_div(d2b_ptr res, d2b_srcptr x, d2b_srcptr y, slong len)
void _d3b_vec_div(d3b_ptr res, d3b_srcptr x, d3b_srcptr y, slong len)
void _d4b_vec_div(d4b_ptr res, d4b_srcptr x, d4b_srcptr y, slong len)
void _d1b_vec_dot(d1b_t res, d1b_srcptr initial, int subtract, d1b_srcptr x, d1b_srcptr y, slong len)
void _d2b_vec_dot(d2b_t res, d2b_srcptr initial, int subtract, d2b_srcptr x, d2b_srcptr y, slong len)
void _d3b_vec_dot(d3b_t res, d3b_srcptr initial, int subtract, d3b_srcptr x, d3b_srcptr y, slong len)
void _d4b_vec_dot(d4b_t res, d4b_srcptr initial, int subtract, d4b_srcptr x, d4b_srcptr y, slong len)
void _d1b_vec_dot_rev(d1b_t res, d1b_srcptr initial, int subtract, d1b_srcptr x, d1b_srcptr y, slong len)
void _d2b_vec_dot_rev(d2b_t res, d2b_srcptr initial, int subtract, d2b_srcptr x, d2b_srcptr y, slong len)
void _d3b_vec_dot_rev(d3b_t res, d3b_srcptr initial, int subtract, d3b_srcptr x, d3b_srcptr y, slong len)
void _d4b_vec_dot_rev(d4b_t res, d4b_srcptr initial, int subtract, d4b_srcptr x, d4b_srcptr y, slong len)

The ball versions of the vector arithmetic, with the same results as
the scalar ball operations.
:::

::: function
int _d1b_vec_sqrt(d1b_ptr res, d1b_srcptr x, slong len)
int _d2b_vec_sqrt(d2b_ptr res, d2b_srcptr x, slong len)
int _d3b_vec_sqrt(d3b_ptr res, d3b_srcptr x, slong len)
int _d4b_vec_sqrt(d4b_ptr res, d4b_srcptr x, slong len)
int _d1b_vec_rsqrt(d1b_ptr res, d1b_srcptr x, slong len)
int _d2b_vec_rsqrt(d2b_ptr res, d2b_srcptr x, slong len)
int _d3b_vec_rsqrt(d3b_ptr res, d3b_srcptr x, slong len)
int _d4b_vec_rsqrt(d4b_ptr res, d4b_srcptr x, slong len)
void _d1b_vec_exp(d1b_ptr res, d1b_srcptr x, slong len)
void _d2b_vec_exp(d2b_ptr res, d2b_srcptr x, slong len)
void _d3b_vec_exp(d3b_ptr res, d3b_srcptr x, slong len)
void _d4b_vec_exp(d4b_ptr res, d4b_srcptr x, slong len)
void _d1b_vec_sin(d1b_ptr res, d1b_srcptr x, slong len)
void _d2b_vec_sin(d2b_ptr res, d2b_srcptr x, slong len)
void _d3b_vec_sin(d3b_ptr res, d3b_srcptr x, slong len)
void _d4b_vec_sin(d4b_ptr res, d4b_srcptr x, slong len)
void _d1b_vec_cos(d1b_ptr res, d1b_srcptr x, slong len)
void _d2b_vec_cos(d2b_ptr res, d2b_srcptr x, slong len)
void _d3b_vec_cos(d3b_ptr res, d3b_srcptr x, slong len)
void _d4b_vec_cos(d4b_ptr res, d4b_srcptr x, slong len)
void _d1b_vec_sin_cos(d1b_ptr sres, d1b_ptr cres, d1b_srcptr x, slong len)
void _d2b_vec_sin_cos(d2b_ptr sres, d2b_ptr cres, d2b_srcptr x, slong len)
void _d3b_vec_sin_cos(d3b_ptr sres, d3b_ptr cres, d3b_srcptr x, slong len)
void _d4b_vec_sin_cos(d4b_ptr sres, d4b_ptr cres, d4b_srcptr x, slong len)
int _d1b_vec_log(d1b_ptr res, d1b_srcptr x, slong len)
int _d2b_vec_log(d2b_ptr res, d2b_srcptr x, slong len)
int _d3b_vec_log(d3b_ptr res, d3b_srcptr x, slong len)
int _d4b_vec_log(d4b_ptr res, d4b_srcptr x, slong len)
void _d1b_vec_atan(d1b_ptr res, d1b_srcptr x, slong len)
void _d2b_vec_atan(d2b_ptr res, d2b_srcptr x, slong len)
void _d3b_vec_atan(d3b_ptr res, d3b_srcptr x, slong len)
void _d4b_vec_atan(d4b_ptr res, d4b_srcptr x, slong len)

Elementwise square roots and elementary functions of balls.
`_dNb_vec_sqrt`, `_dNb_vec_rsqrt` and `_dNb_vec_log` return the bitwise
or of the statuses of the elements; where an element fails, its output
is unspecified.
:::

::: function
void dfloat_set_dfloat(double * res, int nres, const double * x, int nx)
void dfloat_ball_set_ball(double * res, int nres, const double * x, int nx)

Convert the `nx`-term expansion `x` to `nres` terms. Converting to a
longer format is exact; converting to a shorter one rounds, and for
balls the dropped part is added to the radius. In the ball version the
radius is stored after the midpoint, at `x[nx]` and `res[nres]`, and a
result whose magnitude overflows becomes the whole real line. Requires
$1 \le$ `nx`, `nres` $\le 4$.
:::

## Complex numbers {#complex}

The complex formats `dNc` are pairs of `dN` numbers and the complex ball
formats `dNcb` pairs of `dNb` balls, representing every $a + bi$ with
$a$ in the real ball and $b$ in the imaginary ball (rectangular balls,
as `acb_t`). The parts are stored consecutively, so that an array of
$n$ complex numbers is also an array of $2n$ real numbers of the same
format. Each part follows the conventions of the real balls; see
[](dfloat/balls#complex) for the domains of the complex functions and
their status codes. The plain complex formats are approximate, with
results normwise accurate to a few ulps where the functions are well
conditioned.

::: type
d1c_struct
d2c_struct
d3c_struct
d4c_struct
d1c_t
d2c_t
d3c_t
d4c_t
d1c_ptr
d2c_ptr
d3c_ptr
d4c_ptr
d1c_srcptr
d2c_srcptr
d3c_srcptr
d4c_srcptr

A `dNc_struct` has two `dN_struct` members `re` and `im`.
:::

::: type
d1cb_struct
d2cb_struct
d3cb_struct
d4cb_struct
d1cb_t
d2cb_t
d3cb_t
d4cb_t
d1cb_ptr
d2cb_ptr
d3cb_ptr
d4cb_ptr
d1cb_srcptr
d2cb_srcptr
d3cb_srcptr
d4cb_srcptr

A `dNcb_struct` has two `dNb_struct` members `re` and `im`.
:::

::: macro
d1c_realref(x)
d2c_realref(x)
d3c_realref(x)
d4c_realref(x)
d1c_imagref(x)
d2c_imagref(x)
d3c_imagref(x)
d4c_imagref(x)
d1cb_realref(x)
d2cb_realref(x)
d3cb_realref(x)
d4cb_realref(x)
d1cb_imagref(x)
d2cb_imagref(x)
d3cb_imagref(x)
d4cb_imagref(x)

Pointers to the real and imaginary parts of `x`.
:::

::: function
void d1c_zero(d1c_t res)
void d2c_zero(d2c_t res)
void d3c_zero(d3c_t res)
void d4c_zero(d4c_t res)
void d1c_one(d1c_t res)
void d2c_one(d2c_t res)
void d3c_one(d3c_t res)
void d4c_one(d4c_t res)
void d1c_onei(d1c_t res)
void d2c_onei(d2c_t res)
void d3c_onei(d3c_t res)
void d4c_onei(d4c_t res)
void d1c_set(d1c_t res, d1c_srcptr x)
void d2c_set(d2c_t res, d2c_srcptr x)
void d3c_set(d3c_t res, d3c_srcptr x)
void d4c_set(d4c_t res, d4c_srcptr x)
void d1c_set_d(d1c_t res, double x)
void d2c_set_d(d2c_t res, double x)
void d3c_set_d(d3c_t res, double x)
void d4c_set_d(d4c_t res, double x)
void d1c_set_d_d(d1c_t res, double x, double y)
void d2c_set_d_d(d2c_t res, double x, double y)
void d3c_set_d_d(d3c_t res, double x, double y)
void d4c_set_d_d(d4c_t res, double x, double y)
void d1c_set_re_im(d1c_t res, d1_srcptr x, d1_srcptr y)
void d2c_set_re_im(d2c_t res, d2_srcptr x, d2_srcptr y)
void d3c_set_re_im(d3c_t res, d3_srcptr x, d3_srcptr y)
void d4c_set_re_im(d4c_t res, d4_srcptr x, d4_srcptr y)
void d1cb_zero(d1cb_t res)
void d2cb_zero(d2cb_t res)
void d3cb_zero(d3cb_t res)
void d4cb_zero(d4cb_t res)
void d1cb_one(d1cb_t res)
void d2cb_one(d2cb_t res)
void d3cb_one(d3cb_t res)
void d4cb_one(d4cb_t res)
void d1cb_onei(d1cb_t res)
void d2cb_onei(d2cb_t res)
void d3cb_onei(d3cb_t res)
void d4cb_onei(d4cb_t res)
void d1cb_indeterminate(d1cb_t res)
void d2cb_indeterminate(d2cb_t res)
void d3cb_indeterminate(d3cb_t res)
void d4cb_indeterminate(d4cb_t res)
void d1cb_set(d1cb_t res, d1cb_srcptr x)
void d2cb_set(d2cb_t res, d2cb_srcptr x)
void d3cb_set(d3cb_t res, d3cb_srcptr x)
void d4cb_set(d4cb_t res, d4cb_srcptr x)
void d1cb_set_d(d1cb_t res, double x)
void d2cb_set_d(d2cb_t res, double x)
void d3cb_set_d(d3cb_t res, double x)
void d4cb_set_d(d4cb_t res, double x)
void d1cb_set_d_d(d1cb_t res, double x, double y)
void d2cb_set_d_d(d2cb_t res, double x, double y)
void d3cb_set_d_d(d3cb_t res, double x, double y)
void d4cb_set_d_d(d4cb_t res, double x, double y)
void d1cb_set_re_im(d1cb_t res, d1b_srcptr x, d1b_srcptr y)
void d2cb_set_re_im(d2cb_t res, d2b_srcptr x, d2b_srcptr y)
void d3cb_set_re_im(d3cb_t res, d3b_srcptr x, d3b_srcptr y)
void d4cb_set_re_im(d4cb_t res, d4b_srcptr x, d4b_srcptr y)

Assignments: 0, 1, $i$, a copy, a double, $x + yi$ from two doubles,
and $x + yi$ from two real elements; `indeterminate` is the whole
complex plane. `set` and `set_re_im` copy the parts as they are.
:::

::: function
void d1c_neg(d1c_t res, d1c_srcptr x)
void d2c_neg(d2c_t res, d2c_srcptr x)
void d3c_neg(d3c_t res, d3c_srcptr x)
void d4c_neg(d4c_t res, d4c_srcptr x)
void d1c_conj(d1c_t res, d1c_srcptr x)
void d2c_conj(d2c_t res, d2c_srcptr x)
void d3c_conj(d3c_t res, d3c_srcptr x)
void d4c_conj(d4c_t res, d4c_srcptr x)
void d1c_add(d1c_t res, d1c_srcptr x, d1c_srcptr y)
void d2c_add(d2c_t res, d2c_srcptr x, d2c_srcptr y)
void d3c_add(d3c_t res, d3c_srcptr x, d3c_srcptr y)
void d4c_add(d4c_t res, d4c_srcptr x, d4c_srcptr y)
void d1c_sub(d1c_t res, d1c_srcptr x, d1c_srcptr y)
void d2c_sub(d2c_t res, d2c_srcptr x, d2c_srcptr y)
void d3c_sub(d3c_t res, d3c_srcptr x, d3c_srcptr y)
void d4c_sub(d4c_t res, d4c_srcptr x, d4c_srcptr y)
void d1c_mul(d1c_t res, d1c_srcptr x, d1c_srcptr y)
void d2c_mul(d2c_t res, d2c_srcptr x, d2c_srcptr y)
void d3c_mul(d3c_t res, d3c_srcptr x, d3c_srcptr y)
void d4c_mul(d4c_t res, d4c_srcptr x, d4c_srcptr y)
void d1c_sqr(d1c_t res, d1c_srcptr x)
void d2c_sqr(d2c_t res, d2c_srcptr x)
void d3c_sqr(d3c_t res, d3c_srcptr x)
void d4c_sqr(d4c_t res, d4c_srcptr x)
void d1c_mul_re(d1c_t res, d1c_srcptr x, d1_srcptr y)
void d2c_mul_re(d2c_t res, d2c_srcptr x, d2_srcptr y)
void d3c_mul_re(d3c_t res, d3c_srcptr x, d3_srcptr y)
void d4c_mul_re(d4c_t res, d4c_srcptr x, d4_srcptr y)
void d1c_mul_onei(d1c_t res, d1c_srcptr x)
void d2c_mul_onei(d2c_t res, d2c_srcptr x)
void d3c_mul_onei(d3c_t res, d3c_srcptr x)
void d4c_mul_onei(d4c_t res, d4c_srcptr x)
void d1c_mul_2exp_si(d1c_t res, d1c_srcptr x, slong e)
void d2c_mul_2exp_si(d2c_t res, d2c_srcptr x, slong e)
void d3c_mul_2exp_si(d3c_t res, d3c_srcptr x, slong e)
void d4c_mul_2exp_si(d4c_t res, d4c_srcptr x, slong e)
void d1c_div(d1c_t res, d1c_srcptr x, d1c_srcptr y)
void d2c_div(d2c_t res, d2c_srcptr x, d2c_srcptr y)
void d3c_div(d3c_t res, d3c_srcptr x, d3c_srcptr y)
void d4c_div(d4c_t res, d4c_srcptr x, d4c_srcptr y)
void d1c_inv(d1c_t res, d1c_srcptr x)
void d2c_inv(d2c_t res, d2c_srcptr x)
void d3c_inv(d3c_t res, d3c_srcptr x)
void d4c_inv(d4c_t res, d4c_srcptr x)
void d1c_sqrt(d1c_t res, d1c_srcptr x)
void d2c_sqrt(d2c_t res, d2c_srcptr x)
void d3c_sqrt(d3c_t res, d3c_srcptr x)
void d4c_sqrt(d4c_t res, d4c_srcptr x)
void d1c_rsqrt(d1c_t res, d1c_srcptr x)
void d2c_rsqrt(d2c_t res, d2c_srcptr x)
void d3c_rsqrt(d3c_t res, d3c_srcptr x)
void d4c_rsqrt(d4c_t res, d4c_srcptr x)
void d1c_abs(d1_t res, d1c_srcptr x)
void d2c_abs(d2_t res, d2c_srcptr x)
void d3c_abs(d3_t res, d3c_srcptr x)
void d4c_abs(d4_t res, d4c_srcptr x)
void d1c_arg(d1_t res, d1c_srcptr x)
void d2c_arg(d2_t res, d2c_srcptr x)
void d3c_arg(d3_t res, d3c_srcptr x)
void d4c_arg(d4_t res, d4c_srcptr x)
void d1cb_neg(d1cb_t res, d1cb_srcptr x)
void d2cb_neg(d2cb_t res, d2cb_srcptr x)
void d3cb_neg(d3cb_t res, d3cb_srcptr x)
void d4cb_neg(d4cb_t res, d4cb_srcptr x)
void d1cb_conj(d1cb_t res, d1cb_srcptr x)
void d2cb_conj(d2cb_t res, d2cb_srcptr x)
void d3cb_conj(d3cb_t res, d3cb_srcptr x)
void d4cb_conj(d4cb_t res, d4cb_srcptr x)
void d1cb_add(d1cb_t res, d1cb_srcptr x, d1cb_srcptr y)
void d2cb_add(d2cb_t res, d2cb_srcptr x, d2cb_srcptr y)
void d3cb_add(d3cb_t res, d3cb_srcptr x, d3cb_srcptr y)
void d4cb_add(d4cb_t res, d4cb_srcptr x, d4cb_srcptr y)
void d1cb_sub(d1cb_t res, d1cb_srcptr x, d1cb_srcptr y)
void d2cb_sub(d2cb_t res, d2cb_srcptr x, d2cb_srcptr y)
void d3cb_sub(d3cb_t res, d3cb_srcptr x, d3cb_srcptr y)
void d4cb_sub(d4cb_t res, d4cb_srcptr x, d4cb_srcptr y)
void d1cb_mul(d1cb_t res, d1cb_srcptr x, d1cb_srcptr y)
void d2cb_mul(d2cb_t res, d2cb_srcptr x, d2cb_srcptr y)
void d3cb_mul(d3cb_t res, d3cb_srcptr x, d3cb_srcptr y)
void d4cb_mul(d4cb_t res, d4cb_srcptr x, d4cb_srcptr y)
void d1cb_sqr(d1cb_t res, d1cb_srcptr x)
void d2cb_sqr(d2cb_t res, d2cb_srcptr x)
void d3cb_sqr(d3cb_t res, d3cb_srcptr x)
void d4cb_sqr(d4cb_t res, d4cb_srcptr x)
void d1cb_mul_re(d1cb_t res, d1cb_srcptr x, d1b_srcptr y)
void d2cb_mul_re(d2cb_t res, d2cb_srcptr x, d2b_srcptr y)
void d3cb_mul_re(d3cb_t res, d3cb_srcptr x, d3b_srcptr y)
void d4cb_mul_re(d4cb_t res, d4cb_srcptr x, d4b_srcptr y)
void d1cb_mul_onei(d1cb_t res, d1cb_srcptr x)
void d2cb_mul_onei(d2cb_t res, d2cb_srcptr x)
void d3cb_mul_onei(d3cb_t res, d3cb_srcptr x)
void d4cb_mul_onei(d4cb_t res, d4cb_srcptr x)
void d1cb_mul_2exp_si(d1cb_t res, d1cb_srcptr x, slong e)
void d2cb_mul_2exp_si(d2cb_t res, d2cb_srcptr x, slong e)
void d3cb_mul_2exp_si(d3cb_t res, d3cb_srcptr x, slong e)
void d4cb_mul_2exp_si(d4cb_t res, d4cb_srcptr x, slong e)
int d1cb_div(d1cb_t res, d1cb_srcptr x, d1cb_srcptr y)
int d2cb_div(d2cb_t res, d2cb_srcptr x, d2cb_srcptr y)
int d3cb_div(d3cb_t res, d3cb_srcptr x, d3cb_srcptr y)
int d4cb_div(d4cb_t res, d4cb_srcptr x, d4cb_srcptr y)
int d1cb_inv(d1cb_t res, d1cb_srcptr x)
int d2cb_inv(d2cb_t res, d2cb_srcptr x)
int d3cb_inv(d3cb_t res, d3cb_srcptr x)
int d4cb_inv(d4cb_t res, d4cb_srcptr x)
void d1cb_sqrt(d1cb_t res, d1cb_srcptr x)
void d2cb_sqrt(d2cb_t res, d2cb_srcptr x)
void d3cb_sqrt(d3cb_t res, d3cb_srcptr x)
void d4cb_sqrt(d4cb_t res, d4cb_srcptr x)
int d1cb_rsqrt(d1cb_t res, d1cb_srcptr x)
int d2cb_rsqrt(d2cb_t res, d2cb_srcptr x)
int d3cb_rsqrt(d3cb_t res, d3cb_srcptr x)
int d4cb_rsqrt(d4cb_t res, d4cb_srcptr x)
void d1cb_abs(d1b_t res, d1cb_srcptr x)
void d2cb_abs(d2b_t res, d2cb_srcptr x)
void d3cb_abs(d3b_t res, d3cb_srcptr x)
void d4cb_abs(d4b_t res, d4cb_srcptr x)
void d1cb_arg(d1b_t res, d1cb_srcptr x)
void d2cb_arg(d2b_t res, d2cb_srcptr x)
void d3cb_arg(d3b_t res, d3cb_srcptr x)
void d4cb_arg(d4b_t res, d4cb_srcptr x)

Arithmetic: negation, conjugate, sum, difference, product, square,
product with a real number (`mul_re`) or with $i$ (`mul_onei`),
multiplication by $2^e$, quotient, inverse, principal square root and
reciprocal square root, absolute value and argument
$\operatorname{atan2}(b, a) \in (-\pi, \pi]$ (0 for 0) as real
elements. The ball versions of `div`, `inv` and `rsqrt` return a status
(`GR_DOMAIN` only for an exactly zero divisor or argument, `GR_UNABLE`
if it contains zero).

::: details
The product is $(ac - bd) + (ad + bc)i$, or two real products when
either factor is exactly real or imaginary (which is also tighter for
balls); the square is $(a + b)(a - b) + 2abi$ for the plain formats (a
relative real part) and $(a^2 - b^2) + 2abi$ for the balls (a smaller
radius). The quotient is $x\bar{y}/|y|^2$ after scaling $y$ by a power
of two; for balls, where this fails for a divisor away from zero (a
wide rectangle), the result is the enclosure $|x/y| \le \max|x|/\min|y|$
of both parts. The square root is $t + ui$ with
$t = \sqrt{(|z| + a)/2}$ and $u = b/(2t)$ for $a \ge 0$, and the
symmetric form for $a < 0$, so that neither part cancels. The plain
division, inverse, square roots and absolute value have inlined fast
paths for operands whose parts have heads in $[2^{-250}, 2^{250}]$.
:::
:::

::: function
int d1c_is_zero(d1c_srcptr x)
int d2c_is_zero(d2c_srcptr x)
int d3c_is_zero(d3c_srcptr x)
int d4c_is_zero(d4c_srcptr x)
int d1c_is_one(d1c_srcptr x)
int d2c_is_one(d2c_srcptr x)
int d3c_is_one(d3c_srcptr x)
int d4c_is_one(d4c_srcptr x)
int d1c_is_real(d1c_srcptr x)
int d2c_is_real(d2c_srcptr x)
int d3c_is_real(d3c_srcptr x)
int d4c_is_real(d4c_srcptr x)
int d1c_is_finite(d1c_srcptr x)
int d2c_is_finite(d2c_srcptr x)
int d3c_is_finite(d3c_srcptr x)
int d4c_is_finite(d4c_srcptr x)
int d1c_equal(d1c_srcptr x, d1c_srcptr y)
int d2c_equal(d2c_srcptr x, d2c_srcptr y)
int d3c_equal(d3c_srcptr x, d3c_srcptr y)
int d4c_equal(d4c_srcptr x, d4c_srcptr y)
int d1cb_is_zero(d1cb_srcptr x)
int d2cb_is_zero(d2cb_srcptr x)
int d3cb_is_zero(d3cb_srcptr x)
int d4cb_is_zero(d4cb_srcptr x)
int d1cb_is_one(d1cb_srcptr x)
int d2cb_is_one(d2cb_srcptr x)
int d3cb_is_one(d3cb_srcptr x)
int d4cb_is_one(d4cb_srcptr x)
int d1cb_is_real(d1cb_srcptr x)
int d2cb_is_real(d2cb_srcptr x)
int d3cb_is_real(d3cb_srcptr x)
int d4cb_is_real(d4cb_srcptr x)
int d1cb_is_exact(d1cb_srcptr x)
int d2cb_is_exact(d2cb_srcptr x)
int d3cb_is_exact(d3cb_srcptr x)
int d4cb_is_exact(d4cb_srcptr x)
int d1cb_is_finite(d1cb_srcptr x)
int d2cb_is_finite(d2cb_srcptr x)
int d3cb_is_finite(d3cb_srcptr x)
int d4cb_is_finite(d4cb_srcptr x)
int d1cb_contains_zero(d1cb_srcptr x)
int d2cb_contains_zero(d2cb_srcptr x)
int d3cb_contains_zero(d3cb_srcptr x)
int d4cb_contains_zero(d4cb_srcptr x)
int d1cb_contains(d1cb_srcptr x, d1cb_srcptr y)
int d2cb_contains(d2cb_srcptr x, d2cb_srcptr y)
int d3cb_contains(d3cb_srcptr x, d3cb_srcptr y)
int d4cb_contains(d4cb_srcptr x, d4cb_srcptr y)
int d1cb_overlaps(d1cb_srcptr x, d1cb_srcptr y)
int d2cb_overlaps(d2cb_srcptr x, d2cb_srcptr y)
int d3cb_overlaps(d3cb_srcptr x, d3cb_srcptr y)
int d4cb_overlaps(d4cb_srcptr x, d4cb_srcptr y)
int d1cb_equal(d1cb_srcptr x, d1cb_srcptr y)
int d2cb_equal(d2cb_srcptr x, d2cb_srcptr y)
int d3cb_equal(d3cb_srcptr x, d3cb_srcptr y)
int d4cb_equal(d4cb_srcptr x, d4cb_srcptr y)

Predicates, combining those of the parts; `is_real` means that the
imaginary part is exactly zero. For balls, answers of 1 are certain
except for `contains_zero` and `overlaps`, whose answers of 0 are.
:::

::: function
void d1c_exp(d1c_t res, d1c_srcptr x)
void d2c_exp(d2c_t res, d2c_srcptr x)
void d3c_exp(d3c_t res, d3c_srcptr x)
void d4c_exp(d4c_t res, d4c_srcptr x)
void d1c_log(d1c_t res, d1c_srcptr x)
void d2c_log(d2c_t res, d2c_srcptr x)
void d3c_log(d3c_t res, d3c_srcptr x)
void d4c_log(d4c_t res, d4c_srcptr x)
void d1c_sin(d1c_t res, d1c_srcptr x)
void d2c_sin(d2c_t res, d2c_srcptr x)
void d3c_sin(d3c_t res, d3c_srcptr x)
void d4c_sin(d4c_t res, d4c_srcptr x)
void d1c_cos(d1c_t res, d1c_srcptr x)
void d2c_cos(d2c_t res, d2c_srcptr x)
void d3c_cos(d3c_t res, d3c_srcptr x)
void d4c_cos(d4c_t res, d4c_srcptr x)
void d1c_sin_cos(d1c_t sn, d1c_t cs, d1c_srcptr x)
void d2c_sin_cos(d2c_t sn, d2c_t cs, d2c_srcptr x)
void d3c_sin_cos(d3c_t sn, d3c_t cs, d3c_srcptr x)
void d4c_sin_cos(d4c_t sn, d4c_t cs, d4c_srcptr x)
void d1c_tan(d1c_t res, d1c_srcptr x)
void d2c_tan(d2c_t res, d2c_srcptr x)
void d3c_tan(d3c_t res, d3c_srcptr x)
void d4c_tan(d4c_t res, d4c_srcptr x)
void d1c_sinh(d1c_t res, d1c_srcptr x)
void d2c_sinh(d2c_t res, d2c_srcptr x)
void d3c_sinh(d3c_t res, d3c_srcptr x)
void d4c_sinh(d4c_t res, d4c_srcptr x)
void d1c_cosh(d1c_t res, d1c_srcptr x)
void d2c_cosh(d2c_t res, d2c_srcptr x)
void d3c_cosh(d3c_t res, d3c_srcptr x)
void d4c_cosh(d4c_t res, d4c_srcptr x)
void d1c_tanh(d1c_t res, d1c_srcptr x)
void d2c_tanh(d2c_t res, d2c_srcptr x)
void d3c_tanh(d3c_t res, d3c_srcptr x)
void d4c_tanh(d4c_t res, d4c_srcptr x)
void d1c_pow(d1c_t res, d1c_srcptr x, d1c_srcptr y)
void d2c_pow(d2c_t res, d2c_srcptr x, d2c_srcptr y)
void d3c_pow(d3c_t res, d3c_srcptr x, d3c_srcptr y)
void d4c_pow(d4c_t res, d4c_srcptr x, d4c_srcptr y)
void d1cb_exp(d1cb_t res, d1cb_srcptr x)
void d2cb_exp(d2cb_t res, d2cb_srcptr x)
void d3cb_exp(d3cb_t res, d3cb_srcptr x)
void d4cb_exp(d4cb_t res, d4cb_srcptr x)
int d1cb_log(d1cb_t res, d1cb_srcptr x)
int d2cb_log(d2cb_t res, d2cb_srcptr x)
int d3cb_log(d3cb_t res, d3cb_srcptr x)
int d4cb_log(d4cb_t res, d4cb_srcptr x)
void d1cb_sin(d1cb_t res, d1cb_srcptr x)
void d2cb_sin(d2cb_t res, d2cb_srcptr x)
void d3cb_sin(d3cb_t res, d3cb_srcptr x)
void d4cb_sin(d4cb_t res, d4cb_srcptr x)
void d1cb_cos(d1cb_t res, d1cb_srcptr x)
void d2cb_cos(d2cb_t res, d2cb_srcptr x)
void d3cb_cos(d3cb_t res, d3cb_srcptr x)
void d4cb_cos(d4cb_t res, d4cb_srcptr x)
void d1cb_sin_cos(d1cb_t sn, d1cb_t cs, d1cb_srcptr x)
void d2cb_sin_cos(d2cb_t sn, d2cb_t cs, d2cb_srcptr x)
void d3cb_sin_cos(d3cb_t sn, d3cb_t cs, d3cb_srcptr x)
void d4cb_sin_cos(d4cb_t sn, d4cb_t cs, d4cb_srcptr x)
int d1cb_tan(d1cb_t res, d1cb_srcptr x)
int d2cb_tan(d2cb_t res, d2cb_srcptr x)
int d3cb_tan(d3cb_t res, d3cb_srcptr x)
int d4cb_tan(d4cb_t res, d4cb_srcptr x)
void d1cb_sinh(d1cb_t res, d1cb_srcptr x)
void d2cb_sinh(d2cb_t res, d2cb_srcptr x)
void d3cb_sinh(d3cb_t res, d3cb_srcptr x)
void d4cb_sinh(d4cb_t res, d4cb_srcptr x)
void d1cb_cosh(d1cb_t res, d1cb_srcptr x)
void d2cb_cosh(d2cb_t res, d2cb_srcptr x)
void d3cb_cosh(d3cb_t res, d3cb_srcptr x)
void d4cb_cosh(d4cb_t res, d4cb_srcptr x)
int d1cb_tanh(d1cb_t res, d1cb_srcptr x)
int d2cb_tanh(d2cb_t res, d2cb_srcptr x)
int d3cb_tanh(d3cb_t res, d3cb_srcptr x)
int d4cb_tanh(d4cb_t res, d4cb_srcptr x)
int d1cb_pow(d1cb_t res, d1cb_srcptr x, d1cb_srcptr y)
int d2cb_pow(d2cb_t res, d2cb_srcptr x, d2cb_srcptr y)
int d3cb_pow(d3cb_t res, d3cb_srcptr x, d3cb_srcptr y)
int d4cb_pow(d4cb_t res, d4cb_srcptr x, d4cb_srcptr y)

Elementary functions (principal branches, with the branch cut of
`log`, `sqrt` and `pow` on $(-\infty, 0]$, continuous from above), as
compositions of the real functions. The ball versions of `log`, `tan`,
`tanh` and `pow` return a status following the domains in
[](dfloat/balls#complex).

::: details
$\exp(a + bi) = e^a(\cos b + i\sin b)$;
$\log z = \log|z| + i\operatorname{atan2}(b, a)$, where $\log|z|$ is
$\operatorname{log1p}((p - 1)(p + 1) + q^2)/2$ ($p$ the larger part)
when $|z|^2$ is near 1, so that it stays relative; $\sin$, $\cos$,
$\sinh$, $\cosh$ by the addition formulas; $\tan z = (t + ui)/(1 - tui)$
with $t = \tan a$ and $u = \tanh b$, or, where $\tan a$ fails and $b$ is
away from zero, $i(1 - w)/(1 + w)$ with $w = e^{2iz}$ taken with
$|w| < 1$; $\tanh z = -i\tan(iz)$; and $x^y$ by binary powering for an
integer $y$ with $|y| < 2^{40}$ and $\exp(y \log x)$ otherwise.
:::
:::

::: function
void d1c_get_acb(acb_t res, d1c_srcptr x)
void d2c_get_acb(acb_t res, d2c_srcptr x)
void d3c_get_acb(acb_t res, d3c_srcptr x)
void d4c_get_acb(acb_t res, d4c_srcptr x)
void d1c_set_acb(d1c_t res, const acb_t x)
void d2c_set_acb(d2c_t res, const acb_t x)
void d3c_set_acb(d3c_t res, const acb_t x)
void d4c_set_acb(d4c_t res, const acb_t x)
char * d1c_get_str(d1c_srcptr x, slong digits)
char * d2c_get_str(d2c_srcptr x, slong digits)
char * d3c_get_str(d3c_srcptr x, slong digits)
char * d4c_get_str(d4c_srcptr x, slong digits)
void d1c_print(d1c_srcptr x)
void d2c_print(d2c_srcptr x)
void d3c_print(d3c_srcptr x)
void d4c_print(d4c_srcptr x)
void d1c_randtest(d1c_t res, flint_rand_t state)
void d2c_randtest(d2c_t res, flint_rand_t state)
void d3c_randtest(d3c_t res, flint_rand_t state)
void d4c_randtest(d4c_t res, flint_rand_t state)
void d1c_randtest_special(d1c_t res, flint_rand_t state)
void d2c_randtest_special(d2c_t res, flint_rand_t state)
void d3c_randtest_special(d3c_t res, flint_rand_t state)
void d4c_randtest_special(d4c_t res, flint_rand_t state)
void d1cb_get_acb(acb_t res, d1cb_srcptr x)
void d2cb_get_acb(acb_t res, d2cb_srcptr x)
void d3cb_get_acb(acb_t res, d3cb_srcptr x)
void d4cb_get_acb(acb_t res, d4cb_srcptr x)
void d1cb_set_acb(d1cb_t res, const acb_t x)
void d2cb_set_acb(d2cb_t res, const acb_t x)
void d3cb_set_acb(d3cb_t res, const acb_t x)
void d4cb_set_acb(d4cb_t res, const acb_t x)
char * d1cb_get_str(d1cb_srcptr x, slong digits)
char * d2cb_get_str(d2cb_srcptr x, slong digits)
char * d3cb_get_str(d3cb_srcptr x, slong digits)
char * d4cb_get_str(d4cb_srcptr x, slong digits)
void d1cb_print(d1cb_srcptr x)
void d2cb_print(d2cb_srcptr x)
void d3cb_print(d3cb_srcptr x)
void d4cb_print(d4cb_srcptr x)
void d1cb_randtest(d1cb_t res, flint_rand_t state)
void d2cb_randtest(d2cb_t res, flint_rand_t state)
void d3cb_randtest(d3cb_t res, flint_rand_t state)
void d4cb_randtest(d4cb_t res, flint_rand_t state)
void d1cb_randtest_special(d1cb_t res, flint_rand_t state)
void d2cb_randtest_special(d2cb_t res, flint_rand_t state)
void d3cb_randtest_special(d3cb_t res, flint_rand_t state)
void d4cb_randtest_special(d4cb_t res, flint_rand_t state)

Conversions from and to `acb_t`, part by part as for the real formats;
strings in the format of the `acb_t` ring (`re`, `im*I` or
`(re + im*I)`); and random elements (including real, imaginary and zero
ones).
:::

::: function
void _d1c_vec_add(d1c_ptr res, d1c_srcptr x, d1c_srcptr y, slong len)
void _d2c_vec_add(d2c_ptr res, d2c_srcptr x, d2c_srcptr y, slong len)
void _d3c_vec_add(d3c_ptr res, d3c_srcptr x, d3c_srcptr y, slong len)
void _d4c_vec_add(d4c_ptr res, d4c_srcptr x, d4c_srcptr y, slong len)
void _d1c_vec_sub(d1c_ptr res, d1c_srcptr x, d1c_srcptr y, slong len)
void _d2c_vec_sub(d2c_ptr res, d2c_srcptr x, d2c_srcptr y, slong len)
void _d3c_vec_sub(d3c_ptr res, d3c_srcptr x, d3c_srcptr y, slong len)
void _d4c_vec_sub(d4c_ptr res, d4c_srcptr x, d4c_srcptr y, slong len)
void _d1c_vec_mul(d1c_ptr res, d1c_srcptr x, d1c_srcptr y, slong len)
void _d2c_vec_mul(d2c_ptr res, d2c_srcptr x, d2c_srcptr y, slong len)
void _d3c_vec_mul(d3c_ptr res, d3c_srcptr x, d3c_srcptr y, slong len)
void _d4c_vec_mul(d4c_ptr res, d4c_srcptr x, d4c_srcptr y, slong len)
void _d1c_vec_mul_scalar(d1c_ptr res, d1c_srcptr x, slong len, d1c_srcptr c)
void _d2c_vec_mul_scalar(d2c_ptr res, d2c_srcptr x, slong len, d2c_srcptr c)
void _d3c_vec_mul_scalar(d3c_ptr res, d3c_srcptr x, slong len, d3c_srcptr c)
void _d4c_vec_mul_scalar(d4c_ptr res, d4c_srcptr x, slong len, d4c_srcptr c)
void _d1c_vec_addmul_scalar(d1c_ptr res, d1c_srcptr x, slong len, d1c_srcptr c)
void _d2c_vec_addmul_scalar(d2c_ptr res, d2c_srcptr x, slong len, d2c_srcptr c)
void _d3c_vec_addmul_scalar(d3c_ptr res, d3c_srcptr x, slong len, d3c_srcptr c)
void _d4c_vec_addmul_scalar(d4c_ptr res, d4c_srcptr x, slong len, d4c_srcptr c)
void _d1c_vec_submul_scalar(d1c_ptr res, d1c_srcptr x, slong len, d1c_srcptr c)
void _d2c_vec_submul_scalar(d2c_ptr res, d2c_srcptr x, slong len, d2c_srcptr c)
void _d3c_vec_submul_scalar(d3c_ptr res, d3c_srcptr x, slong len, d3c_srcptr c)
void _d4c_vec_submul_scalar(d4c_ptr res, d4c_srcptr x, slong len, d4c_srcptr c)
void _d1c_vec_dot(d1c_t res, d1c_srcptr initial, int subtract, d1c_srcptr x, d1c_srcptr y, slong len)
void _d2c_vec_dot(d2c_t res, d2c_srcptr initial, int subtract, d2c_srcptr x, d2c_srcptr y, slong len)
void _d3c_vec_dot(d3c_t res, d3c_srcptr initial, int subtract, d3c_srcptr x, d3c_srcptr y, slong len)
void _d4c_vec_dot(d4c_t res, d4c_srcptr initial, int subtract, d4c_srcptr x, d4c_srcptr y, slong len)
void _d1c_vec_dot_rev(d1c_t res, d1c_srcptr initial, int subtract, d1c_srcptr x, d1c_srcptr y, slong len)
void _d2c_vec_dot_rev(d2c_t res, d2c_srcptr initial, int subtract, d2c_srcptr x, d2c_srcptr y, slong len)
void _d3c_vec_dot_rev(d3c_t res, d3c_srcptr initial, int subtract, d3c_srcptr x, d3c_srcptr y, slong len)
void _d4c_vec_dot_rev(d4c_t res, d4c_srcptr initial, int subtract, d4c_srcptr x, d4c_srcptr y, slong len)
void _d1c_vec_exp(d1c_ptr res, d1c_srcptr x, slong len)
void _d2c_vec_exp(d2c_ptr res, d2c_srcptr x, slong len)
void _d3c_vec_exp(d3c_ptr res, d3c_srcptr x, slong len)
void _d4c_vec_exp(d4c_ptr res, d4c_srcptr x, slong len)
void _d1cb_vec_add(d1cb_ptr res, d1cb_srcptr x, d1cb_srcptr y, slong len)
void _d2cb_vec_add(d2cb_ptr res, d2cb_srcptr x, d2cb_srcptr y, slong len)
void _d3cb_vec_add(d3cb_ptr res, d3cb_srcptr x, d3cb_srcptr y, slong len)
void _d4cb_vec_add(d4cb_ptr res, d4cb_srcptr x, d4cb_srcptr y, slong len)
void _d1cb_vec_sub(d1cb_ptr res, d1cb_srcptr x, d1cb_srcptr y, slong len)
void _d2cb_vec_sub(d2cb_ptr res, d2cb_srcptr x, d2cb_srcptr y, slong len)
void _d3cb_vec_sub(d3cb_ptr res, d3cb_srcptr x, d3cb_srcptr y, slong len)
void _d4cb_vec_sub(d4cb_ptr res, d4cb_srcptr x, d4cb_srcptr y, slong len)
void _d1cb_vec_mul(d1cb_ptr res, d1cb_srcptr x, d1cb_srcptr y, slong len)
void _d2cb_vec_mul(d2cb_ptr res, d2cb_srcptr x, d2cb_srcptr y, slong len)
void _d3cb_vec_mul(d3cb_ptr res, d3cb_srcptr x, d3cb_srcptr y, slong len)
void _d4cb_vec_mul(d4cb_ptr res, d4cb_srcptr x, d4cb_srcptr y, slong len)
void _d1cb_vec_mul_scalar(d1cb_ptr res, d1cb_srcptr x, slong len, d1cb_srcptr c)
void _d2cb_vec_mul_scalar(d2cb_ptr res, d2cb_srcptr x, slong len, d2cb_srcptr c)
void _d3cb_vec_mul_scalar(d3cb_ptr res, d3cb_srcptr x, slong len, d3cb_srcptr c)
void _d4cb_vec_mul_scalar(d4cb_ptr res, d4cb_srcptr x, slong len, d4cb_srcptr c)
void _d1cb_vec_addmul_scalar(d1cb_ptr res, d1cb_srcptr x, slong len, d1cb_srcptr c)
void _d2cb_vec_addmul_scalar(d2cb_ptr res, d2cb_srcptr x, slong len, d2cb_srcptr c)
void _d3cb_vec_addmul_scalar(d3cb_ptr res, d3cb_srcptr x, slong len, d3cb_srcptr c)
void _d4cb_vec_addmul_scalar(d4cb_ptr res, d4cb_srcptr x, slong len, d4cb_srcptr c)
void _d1cb_vec_submul_scalar(d1cb_ptr res, d1cb_srcptr x, slong len, d1cb_srcptr c)
void _d2cb_vec_submul_scalar(d2cb_ptr res, d2cb_srcptr x, slong len, d2cb_srcptr c)
void _d3cb_vec_submul_scalar(d3cb_ptr res, d3cb_srcptr x, slong len, d3cb_srcptr c)
void _d4cb_vec_submul_scalar(d4cb_ptr res, d4cb_srcptr x, slong len, d4cb_srcptr c)
void _d1cb_vec_dot(d1cb_t res, d1cb_srcptr initial, int subtract, d1cb_srcptr x, d1cb_srcptr y, slong len)
void _d2cb_vec_dot(d2cb_t res, d2cb_srcptr initial, int subtract, d2cb_srcptr x, d2cb_srcptr y, slong len)
void _d3cb_vec_dot(d3cb_t res, d3cb_srcptr initial, int subtract, d3cb_srcptr x, d3cb_srcptr y, slong len)
void _d4cb_vec_dot(d4cb_t res, d4cb_srcptr initial, int subtract, d4cb_srcptr x, d4cb_srcptr y, slong len)
void _d1cb_vec_dot_rev(d1cb_t res, d1cb_srcptr initial, int subtract, d1cb_srcptr x, d1cb_srcptr y, slong len)
void _d2cb_vec_dot_rev(d2cb_t res, d2cb_srcptr initial, int subtract, d2cb_srcptr x, d2cb_srcptr y, slong len)
void _d3cb_vec_dot_rev(d3cb_t res, d3cb_srcptr initial, int subtract, d3cb_srcptr x, d3cb_srcptr y, slong len)
void _d4cb_vec_dot_rev(d4cb_t res, d4cb_srcptr initial, int subtract, d4cb_srcptr x, d4cb_srcptr y, slong len)
void _d1cb_vec_exp(d1cb_ptr res, d1cb_srcptr x, slong len)
void _d2cb_vec_exp(d2cb_ptr res, d2cb_srcptr x, slong len)
void _d3cb_vec_exp(d3cb_ptr res, d3cb_srcptr x, slong len)
void _d4cb_vec_exp(d4cb_ptr res, d4cb_srcptr x, slong len)

Vector operations, implemented with the real vector operations:
`vec_add` and `vec_sub` are the real ones on $2\,\mathit{len}$ elements;
the dot products are two real dot products per block of 32 terms (3 to
4 times faster than a loop of complex products); `vec_exp` uses the
real `vec_exp` and `vec_sin_cos` on blocks of 32 elements.
:::

## Sieved power sums and Platt's multi-evaluation {#powsum}

These functions use the `dfloat` balls to accelerate the Riemann–Siegel
and Platt methods for the Riemann zeta function on the critical line
(`acb_dirichlet_zeta_rs_r`, `acb_dirichlet_platt_multieval`). The
algorithms are described in [](dfloat/powsum).

::: function
int gr_powsum_sieved(acb_t res, const acb_t s, ulong N, gr_ctx_t T, gr_ctx_t P, gr_ctx_t F, double dm, double da, ulong M)

Sets `res` to an enclosure of $\sum_{n=1}^N n^{-s}$, as
`acb_dirichlet_powsum_sieved` with `len = 1`, computed with the ball
ring $T$ for the terms, the ball ring $P$ for the phases and optionally
the plain ring $F$ for the composite terms and block sums (`NULL`: all
balls), whose products and sums have relative errors at most `dm` and
`da`. `M` is the table bound (0 for the default). Returns a `gr` status.
:::

::: function
int dfloat_powsum_sieved(acb_t res, const acb_t s, ulong N, slong prec)
int _dfloat_powsum_sieved(acb_t res, const acb_t s, ulong N, int tn, int pn, ulong M)
int _dfloat_powsum_sieved_ball(acb_t res, const acb_t s, ulong N, int tn, int pn, ulong M)

The `dfloat` instantiations of `gr_powsum_sieved`, with terms in
`dTb` and phases in `dPb` for $(T, P)$ = (`tn`, `pn`) one of $(1, 2)$,
$(1, 3)$, $(2, 3)$, $(2, 4)$, $(3, 4)$. `dfloat_powsum_sieved` chooses
$T$ and $P$ for about `prec` bits after the binary point. They return 1
on success and 0, leaving `res` undefined, for unsupported arguments.
:::

::: macro
DFLOAT_POWSUM_MAX_BYTES

The memory limit for the tables of `dfloat_powsum_sieved`: $4 \cdot 10^9$
bytes on 64-bit systems and $5 \cdot 10^8$ on 32-bit systems.
:::

::: function
int _dfloat_platt_smk(acb_ptr S_table, const fmpz * smk_points, const arb_t t0, slong A, slong B, ulong J, slong K, int pn, ulong M, slong prec)
int _dfloat_platt_smk_dd(double * S5, const fmpz * smk_points, const arb_t t0, slong A, slong B, ulong J, slong K, int pn, ulong M)

The sums over $j$ in Platt's multi-evaluation of the scaled Lambda
function (`acb_dirichlet_platt_multieval`): for the buckets $m < AB$ of
the $j \le J$ with `smk_points[m]` $\le j <$ `smk_points[m+1]`, sets
entry $kAB + m$ to $\sum_j z_j b_j^k$ for $0 \le k < K$, where
$z_j = j^{-1/2} e^{-i t_0 \log(j\sqrt{\pi})}$ and
$b_j = \log(j\sqrt{\pi})/(2\pi) - m/B$. Return 1 on success and 0 if
`pn` is not 3 or 4 or the module is not supported.

::: details
`_dfloat_platt_smk_dd` omits the common factor
$c = e^{-it_0\log\sqrt{\pi}}$ and writes compact balls of five doubles
per entry (the real and imaginary midpoints as double-words and a
common radius): `S5[5 (k N + m) + 0..4]` holds Re hi, Re lo, Im hi,
Im lo and the radius, with $N = AB$. `M` bounds the table of terms (0
for the default $\min(J/2, 5 \cdot 10^6)$).

The rest of the multi-evaluation amplifies an error in the moment $k$
by about $2^{24 + 8.6k}$, while $|S_k|$ decreases like $2^{-13k}$. So
the terms $j^{-s}$ and $b_j$ are `d3b` balls (with phases in
`d3b` or `d4b`), $S_0$ is summed in `d3b` balls, and the moments
$1 \le k < 14$ are computed in plain double-double arithmetic and
$k \ge 14$ in double arithmetic, vectorised over the $j$ of a block,
with an a priori error bound added per bucket. The terms are obtained
by sieving: $f(j) = f(d) f(j/d)$ and $\log j = \log d + \log(j/d)$ for
the smallest prime factor $d$, from a table for $j \le M$ and by
splitting $j = xy$ with $x, y \le M$ beyond. The work is split into
chunks computed in parallel and summed in a fixed order, so the result
does not depend on the number of threads. At the height of the
$10^{15}$-th zero this takes about 0.42 microseconds per $j$ on one
thread, against 7.7 for the `arb_t` version at 228 bits.
:::
:::

## Internal functions {#internal}

### Error-free transformations and radius arithmetic {#internal-eft}

::: macro
DFLOAT_TWO_SUM(s, e, a, b)
DFLOAT_FAST_TWO_SUM(s, e, a, b)
DFLOAT_TWO_PROD(p, e, a, b)

Error-free transformations: set $s$ and $e$ such that
$s + e = a + b$ exactly with $s = \operatorname{fl}(a + b)$ (Knuth's
TwoSum, no precondition; Dekker's FastTwoSum, requiring
$\operatorname{exponent}(a) \ge \operatorname{exponent}(b)$ or $a = 0$),
or $p + e = ab$ exactly with $p = \operatorname{fl}(ab)$ (with `fma`;
requires that $ab$ neither overflows nor underflows, that is
$|ab| \ge 2^{-969}$ or $ab = 0$). See
[](algorithms/floating-point/double-word#eft).
:::

::: macro
DFLOAT_RAD_SLACK
DFLOAT_RAD_UNDERFLOW_THRESHOLD
DFLOAT_RAD_UNDERFLOW_BUMP

The factor $1 + 2^{-45}$ by which an accumulated radius is multiplied
to make it an upper bound (it covers up to $2^7$ nearest roundings of
normal numbers), the threshold $2^{-1000}$ below which radius products
and quotients are checked for underflow, and the subnormal $2^{-1070}$
added to such a result.
:::

::: function
double _dfloat_rad_mul_slow(double a, double b, double p)
double _dfloat_rad_mul(double a, double b)
double _dfloat_rad_div(double a, double b)
double _dfloat_rad_finish(double r)
double _dfloat_lower(double r)

Radius arithmetic on nonnegative doubles: an upper bound (up to the
final slack factor) for $ab$ and for $a/b$ ($b > 0$), with products and
quotients that may have underflowed bumped up by $2^{-1070}$ and with
$0 \cdot \infty = 0$ (`_dfloat_rad_mul_slow` is the underflow path of
`_dfloat_rad_mul`, given $p = \operatorname{fl}(ab)$); the final upper
bound $r(1 + 2^{-45})$; and a lower bound $r(1 - 2^{-45})$ for a
positive quantity computed with a few roundings (0 for $r < 2^{-1000}$).
:::

::: function
int _dfloat_ball_is_whole(const double * d, int n, double rad)
double _dfloat_abs_sum(const double * d, int n)

Return whether a ball with midpoint components `d[0..n-1]` and radius
`rad` is the whole real line (some component or the radius is
infinite or NaN), and the sum $\sum |d_i|$, an upper bound (up to the
slack factor) for the absolute value of an expansion.
:::

### Generic kernels {#internal-generic}

::: function
void _dfloat_renorm(double * res, int n, double * err, double * t, int len)

Renormalises the sum of the `len` doubles `t` (overwritten) to `n`
nonoverlapping components `res`, adding the magnitude of the dropped
part to `*err` if `err` is not `NULL`. This is a bottom-up VecSum
followed by a top-down pass, with a branchy fallback (VecSumErrBranch
[@JMPT2016]) when a gap would otherwise be left. The terms should be
roughly sorted by decreasing magnitude.
:::

::: function
void _dfloat_ball_mul_scaled(double * res, double * rrad, int n, const double * x, double rx, const double * y, double ry)
void _dfloat_ball_div(double * res, double * rrad, int n, const double * x, double rx, const double * y, double ry)
int _dfloat_ball_rsqrt(double * res, double * rrad, int n, const double * x, double rx)
int _dfloat_ball_sqrt(double * res, double * rrad, int n, const double * x, double rx)

Rigorous product, quotient, reciprocal square root and square root of
balls with midpoints of `n` components and radii `rx`, `ry`, for inputs
of any magnitude, by scaling the heads to $[1, 2)$, applying the
fixed-$N$ kernels and scaling back. Components that would become
subnormal are flushed into the radius. The product is the slow path of
`dNb_mul`; the others are the only implementations.
:::

::: function
void _dfloat_canonicalise(double * res, int n, const double * x)
void _dfloat_add_n(double * res, double * err, int n, const double * x, const double * y)
void _dfloat_sub_n(double * res, double * err, int n, const double * x, const double * y)
void _dfloat_mul_n(double * res, double * err, int n, const double * x, const double * y)
void _dfloat_sqr_n(double * res, double * err, int n, const double * x)
void _dfloat_div_approx_n(double * q, int n, const double * x, const double * y)
void _dfloat_sqrt_approx_n(double * s, int n, const double * x)
void _dfloat_rsqrt_approx_n(double * r, int n, const double * x)

Dispatch to the fixed-$N$ kernels below for a run-time `n`.
:::

::: function
void _d1_add_tracked(double * res, double * err, const double * x, const double * y)
void _d2_add_tracked(double * res, double * err, const double * x, const double * y)
void _d3_add_tracked(double * res, double * err, const double * x, const double * y)
void _d4_add_tracked(double * res, double * err, const double * x, const double * y)
void _d1_sub_tracked(double * res, double * err, const double * x, const double * y)
void _d2_sub_tracked(double * res, double * err, const double * x, const double * y)
void _d3_sub_tracked(double * res, double * err, const double * x, const double * y)
void _d4_sub_tracked(double * res, double * err, const double * x, const double * y)
void _d1_mul_tracked(double * res, double * err, const double * x, const double * y)
void _d2_mul_tracked(double * res, double * err, const double * x, const double * y)
void _d3_mul_tracked(double * res, double * err, const double * x, const double * y)
void _d4_mul_tracked(double * res, double * err, const double * x, const double * y)
void _d1_sqr_tracked(double * res, double * err, const double * x)
void _d2_sqr_tracked(double * res, double * err, const double * x)
void _d3_sqr_tracked(double * res, double * err, const double * x)
void _d4_sqr_tracked(double * res, double * err, const double * x)

The tracked kernels: set `res` to the $N$-term rounded sum, difference,
product or square and add a bound for the magnitude of the dropped part
to `*err`.
:::

::: function
void _d1_div_approx(double * q, const double * x, const double * y)
void _d2_div_approx(double * q, const double * x, const double * y)
void _d3_div_approx(double * q, const double * x, const double * y)
void _d4_div_approx(double * q, const double * x, const double * y)
void _d1_sqrt_approx(double * s, const double * x)
void _d2_sqrt_approx(double * s, const double * x)
void _d3_sqrt_approx(double * s, const double * x)
void _d4_sqrt_approx(double * s, const double * x)
void _d1_rsqrt_approx(double * r, const double * x)
void _d2_rsqrt_approx(double * r, const double * x)
void _d3_rsqrt_approx(double * r, const double * x)
void _d4_rsqrt_approx(double * r, const double * x)

The approximate quotient, square root and reciprocal square root
kernels of the plain formats.
:::

### Conversions and random generation {#internal-conv}

::: function
void _dfloat_get_arf(arf_t res, const double * x, int n)
void _dfloat_set_arf(double * res, int n, double * err, const arf_t x)
void _dfloat_get_arb(arb_t res, const double * x, int n, double rad)
void _dfloat_set_arb(double * res, double * rad, int n, const arb_t x)
int _dfloat_part_get_arb(arb_t res, const double * x, int n, int ball)
int _dfloat_set_str(double * res, double * rad, int n, const char * s)
char * _dfloat_get_str(const double * x, int n, double rad, slong digits)
char * _dfloat_complex_get_str(const double * re, const double * im, int n, double rrad, double irad, slong digits)
void _dfloat_randtest(double * res, int n, flint_rand_t state)
void _dfloat_randtest_special(double * res, int n, flint_rand_t state)
double _dfloat_randtest_rad(flint_rand_t state)

The implementations of the conversions, strings and random generation
for a run-time number of components `n`. A negative radius `rad` in the
string functions means a plain (non-ball) value.
:::

### Elementary function internals {#internal-elementary}

The kernels use tables in `src/dfloat/exp_tables.c`, generated by
`dev/gen_dfloat_exp_tables.py`, of 5-term expansions: `_dfloat_exp_T1`,
`_dfloat_exp_T2` ($e^{i/256}$, $e^{j/65536}$), `_dfloat_exp_E1`,
`_dfloat_exp_E2` (the same minus 1), `_dfloat_exp_ln2`,
`_dfloat_exp_ln2_42`, `_dfloat_exp_invfac`, `_dfloat_trig_A`,
`_dfloat_trig_B`, `_dfloat_trig_A1`, `_dfloat_trig_coef`,
`_dfloat_trig_pi2_33`, `_dfloat_log_R1`, `_dfloat_log_L1`,
`_dfloat_log_R2`, `_dfloat_log_L2`, `_dfloat_log_coef`,
`_dfloat_atan_C`, `_dfloat_pi_2`, `_dfloat_atan_coef` and the
constants `_dfloat_const`, with the bound `_dfloat_exp_tab_err` on the
neglected tails.

::: macro
DFLOAT_CONST_PI
DFLOAT_CONST_EULER
DFLOAT_CONST_CATALAN
DFLOAT_CONST_LOG2

Indices of $\pi$, Euler's constant, Catalan's constant and $\log 2$ in
the table `_dfloat_const`.
:::

::: function
void _dfloat_exp(double * res, int n, const double * x)
void _dfloat_expm1(double * res, int n, const double * x)
void _dfloat_log(double * res, int n, const double * x)
void _dfloat_log1p(double * res, int n, const double * x)
void _dfloat_atan(double * res, int n, const double * x)
void _dfloat_sin_cos(double * sn, double * cs, int n, const double * x, int want)

Dispatch to `dN_exp`, `dN_expm1`, `dN_log`, `dN_log1p`, `dN_atan` and
the sine and cosine kernels for a run-time `n`; `want` is 1 for the
sine, 2 for the cosine and 3 for both.
:::

::: function
void _dfloat_ball_exp(double * res, double * rrad, int n, const double * x, double rx)
void _dfloat_ball_expm1(double * res, double * rrad, int n, const double * x, double rx)
int _dfloat_ball_log(double * res, double * rrad, int n, const double * x, double rx)
int _dfloat_ball_log1p(double * res, double * rrad, int n, const double * x, double rx)
void _dfloat_ball_atan(double * res, double * rrad, int n, const double * x, double rx)
void _dfloat_ball_sin_cos(double * sn, double * srad, double * cs, double * crad, int n, const double * x, double rx, int want)

The ball versions for a run-time `n`; the logarithms return a `gr`
status.
:::

::: function
void _dfloat_exp_slow(double * res, int n, const double * x)
void _dfloat_expm1_slow(double * res, int n, const double * x)
void _dfloat_log_slow(double * res, int n, const double * x)
void _dfloat_log1p_slow(double * res, int n, const double * x)
void _dfloat_atan_slow(double * res, int n, const double * x)
void _dfloat_trig_slow(double * sn, double * cs, int n, const double * x, int want)
int _dfloat_ball_log_slow(double * res, double * rrad, int n, const double * x, double rx)
int _dfloat_ball_log1p_slow(double * res, double * rrad, int n, const double * x, double rx)
void _dfloat_ball_atan_slow(double * res, double * rrad, int n, const double * x, double rx)
void _dfloat_ball_trig_slow(double * sn, double * srad, double * cs, double * crad, int n, const double * x, double rx, int want)

Fallbacks through `arb_t` for any argument, with the same guarantees
as the kernels; used for arguments outside the kernel ranges.
:::

### Complex internals {#internal-complex}

::: macro
DFLOAT_MODE_AUTO
DFLOAT_MODE_FAST

Error modes of ball products and quotients: track rounding errors
exactly for exact inputs (0), or always bound them (2).
:::

::: type
_dfloat_ops_struct
_dfloat_ops_t

A table of the real operations of one format (`n` components, `ball`
flag), as function pointers on arrays of doubles, through which the
complex functions that are not performance-critical are implemented
once for all formats.
:::

::: variable
const _dfloat_ops_struct _d1_ops
const _dfloat_ops_struct _d2_ops
const _dfloat_ops_struct _d3_ops
const _dfloat_ops_struct _d4_ops
const _dfloat_ops_struct _d1b_ops
const _dfloat_ops_struct _d2b_ops
const _dfloat_ops_struct _d3b_ops
const _dfloat_ops_struct _d4b_ops

The operation tables of the eight real formats.
:::

::: function
_dfloat_ops_t _dfloat_ops(int n, int ball)

Returns the operation table of `dN` (`ball` = 0) or `dNb` (`ball` = 1).
:::

::: function
int _dfloat_complex_div(_dfloat_ops_t ops, double * res, const double * x, const double * y, int mode)
int _dfloat_complex_inv(_dfloat_ops_t ops, double * res, const double * x, int mode)
void _dfloat_complex_sqrt(_dfloat_ops_t ops, double * res, const double * x, int mode)
int _dfloat_complex_rsqrt(_dfloat_ops_t ops, double * res, const double * x, int mode)
void _dfloat_complex_abs(_dfloat_ops_t ops, double * res, const double * x)
void _dfloat_complex_arg(_dfloat_ops_t ops, double * res, const double * x)
void _dfloat_complex_exp(_dfloat_ops_t ops, double * res, const double * x)
int _dfloat_complex_log(_dfloat_ops_t ops, double * res, const double * x)
void _dfloat_complex_sin_cos(_dfloat_ops_t ops, double * sn, double * cs, const double * x)
void _dfloat_complex_sinh_cosh(_dfloat_ops_t ops, double * sh, double * ch, const double * x)
int _dfloat_complex_tan(_dfloat_ops_t ops, double * res, const double * x)
int _dfloat_complex_tanh(_dfloat_ops_t ops, double * res, const double * x)
int _dfloat_complex_pow(_dfloat_ops_t ops, double * res, const double * x, const double * y)
void _dfloat_complex_get_acb(_dfloat_ops_t ops, acb_t res, const double * x)
void _dfloat_complex_set_acb(_dfloat_ops_t ops, double * res, const acb_t x)
char * _dfloat_complex_get_str_ops(_dfloat_ops_t ops, const double * x, slong digits)
void _dfloat_complex_randtest(_dfloat_ops_t ops, double * res, flint_rand_t state, int special)

The complex functions for any format, given its operation table; `res`,
`x` and `y` are complex elements (two real elements each). They return
the status of the ball function (`GR_SUCCESS`, `GR_DOMAIN` or
`GR_UNABLE`); the plain versions always compute a result and return
`GR_DOMAIN` only for an exact singularity.
:::

### Generic ring contexts {#internal-ctx}

::: type
_dfloat_ctx_struct

The context data of a `dfloat` ring: the number of components `n` and
the flags `ball`, `strong`, `fast` and `cplx`.
:::

::: macro
DFLOAT_CTX(ctx)
DFLOAT_CTX_N(ctx)
DFLOAT_CTX_BALL(ctx)
DFLOAT_CTX_STRONG(ctx)
DFLOAT_CTX_FAST(ctx)
DFLOAT_CTX_COMPLEX(ctx)
DFLOAT_CTX_PREC(ctx)

Access to the context data; `DFLOAT_CTX_PREC` is $53n$.
:::
