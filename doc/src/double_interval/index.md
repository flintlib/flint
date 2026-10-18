---
title: double_interval.h
header: double_interval.h
description: Fast interval arithmetic with double endpoints
category: Floating-point arithmetic
order: 6
---

This module implements interval arithmetic with `double` endpoints,
for computing cheap but rigorous enclosures where the precision of a
double is enough, such as bounds on the magnitude of an integrand or of
the terms of a series. It is used by the numerical integration code of
`arb_hypgeom` (for example `arb_hypgeom_1f1_integration`). Unlike the
IEEE 1788 approach, which switches the rounding mode of the processor,
the operations here compute the endpoints in the default
round-to-nearest mode and then widen them by a few ulps, which is
valid whether or not the operation was exact. The method and its
correctness are explained in [](algorithms/floating-point/double-intervals).

For rigorous arithmetic at higher precision, use the `arb_t` balls, or
the double-word to quad-double balls `d2b_t`, ..., `d4b_t` of `dfloat`.

## Types {#types}

::: type
di_t

A struct with two `double` members `a` and `b`, representing the
extended real interval $[a, b]$ (the endpoints may be infinite). Unless
stated otherwise, functions require $a \le b$ and that neither endpoint
is NaN. The type is passed and returned by value.
:::

::: macro
DI_CHECK(x)

Throws an exception unless `x.a <= x.b` (in particular, if an endpoint
is NaN).
:::

## Construction and conversion {#construction}

::: function
di_t di_interval(double a, double b)

Returns the interval $[a, b]$. Throws an exception unless $a \le b$
(in particular if $a$ or $b$ is NaN).
:::

::: function
di_t arb_get_di(const arb_t x)

Returns an interval containing the ball $x$: its lower and upper
endpoints, rounded down and up to doubles. A ball with a NaN midpoint
gives $[-\infty, +\infty]$, and an infinite radius gives infinite
endpoints.
:::

::: function
void arb_set_di(arb_t res, di_t x, slong prec)

Sets $\mathit{res}$ to a ball containing the interval $x$, with the
midpoint rounded to $\mathit{prec}$ bits (by `arb_set_interval_arf`).
:::

::: function
void di_print(di_t x)

Prints $x$ to standard output as `[a, b]` with 17 significant digits
per endpoint. The decimal endpoints are not rounded outward, so the
printed interval need not contain $x$.
:::

::: function
double d_randtest2(flint_rand_t state)

Returns a random non-NaN double $\pm m 2^e$ with $m$ from `d_randtest`
and $e$ uniform in $[-1200, 1199]$, so that the result can also be
zero (after underflow), subnormal or infinite.
:::

::: function
di_t di_randtest(flint_rand_t state)

Returns an interval with two random endpoints from `d_randtest2`,
ordered.
:::

## Arithmetic {#arithmetic}

The `di_fast_` functions compute the endpoints with the default
rounding (assumed to be to nearest) and then move the lower endpoint
down and the upper endpoint up by a relative amount $2^{-51}$ plus an
absolute amount of about $2^{-1047}$, regardless of whether the
operations were exact. The results are therefore valid enclosures, but
slightly wider than optimal even for exact operands: the sum of
$[1, 1]$ and $[1, 1]$ is an interval of width about $2^{-49}$
around 2.

::: function
di_t di_neg(di_t x)

Returns $[-b, -a]$, the exact negation of $x = [a, b]$.
:::

::: function
di_t di_fast_add(di_t x, di_t y)
di_t di_fast_sub(di_t x, di_t y)
di_t di_fast_mul(di_t x, di_t y)
di_t di_fast_div(di_t x, di_t y)

Return an enclosure of the sum, difference, product or quotient of $x$
and $y$. The product and quotient select the endpoint products by the
signs of the operands; a product of the form $0 \cdot \infty$ gives
$[-\infty, +\infty]$, and so does division by an interval containing
zero.
:::

::: function
di_t di_fast_sqr(di_t x)

Returns an enclosure of $\{t^2 : t \in x\}$. The lower endpoint is
nonnegative, and exactly 0 if $x$ contains zero.
:::

::: function
di_t di_fast_add_d(di_t x, double y)
di_t di_fast_sub_d(di_t x, double y)
di_t di_fast_mul_d(di_t x, double y)
di_t di_fast_div_d(di_t x, double y)

The same operations with the exact point interval $[y, y]$ as the
second operand. Throws an exception if $y$ is NaN.
:::

::: function
di_t di_fast_log_nonnegative(di_t x)

Returns an enclosure of $\log x$ for $x \subseteq [0, +\infty]$; a
nonpositive lower endpoint gives the lower endpoint $-\infty$. The
endpoints are computed with `mag_d_log_lower_bound` and
`mag_d_log_upper_bound`.
:::

::: function
di_t di_fast_mid(di_t x)

Returns an enclosure of the midpoint $(a + b)/2$ of $x$, or
$[-\infty, +\infty]$ if an endpoint is infinite.
:::

::: function
double di_fast_ubound_radius(di_t x)

Returns an upper bound for the radius $(b - a)/2$ of $x$.
:::

## Internal functions {#internal}

::: function
double _di_below(double x)
double _di_above(double x)

Return a lower respectively upper bound for any real number whose
round-to-nearest value is $x$. `_di_below` returns
$x - (|x| + 10^{-300}) \cdot 2^{-51}$ for $x \le 10^{300}$ (rounded to
nearest), $10^{300}$ for larger $x$ including $+\infty$, and $-\infty$
for NaN (which arises from $\infty - \infty$ in the endpoint
arithmetic); `_di_above` is symmetric.
:::
