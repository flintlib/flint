---
title: arb_calc.h
header: arb_calc.h
description: Calculus with real functions given by ball evaluators: rigorous isolation of the simple real roots of an analytic function by bisection, and refinement by bisection or interval Newton iteration.
category: Real and complex numbers
order: 11
---

This module contains algorithms of real calculus for functions that the
user provides as ball evaluators: a C function that, given a ball $x$,
computes enclosures of $f(x)$ and of its first derivatives. The
algorithms isolate the real roots of $f$ on an interval by bisection,
using the signs of $f$ at points and the range of $f'$ on subintervals,
and refine isolated roots by bisection or by Newton's method in ball
arithmetic, with quadratic convergence and a rigorous error bound at
each step ([](algorithms/ball-arithmetic/real-root-isolation)). All
results are certified: an interval reported to contain exactly one root
does so. For integration and complex functions, see [](acb_calc); for
polynomials with integer coefficients, `arb_fmpz_poly_real_roots` is
much faster. The tutorial [](tutorials/numerics/real-roots) isolates and
refines the roots of a function with this module, and
[](tutorials/numerics/interval-newton) develops the Newton iteration from
scratch.

::: summary algorithms/ball-arithmetic/real-root-isolation :::

## Types and constants {#types}

::: type
arb_calc_func_t

The type of a function evaluator, a pointer to a function

```c
int func(arb_ptr out, const arb_t inp, void * param, slong order, slong prec)
```

which writes to `out[0], ..., out[order - 1]` enclosures of the first
`order` Taylor coefficients $f(x), f'(x), f''(x)/2, \ldots$ of $f$ at
every point of the ball `inp`, computed at a precision of about `prec`
bits. `param` passes user data through. The return value is reserved
for an error code and should be 0.

::: details
`out` and `inp` are never aliased, and `order` is at least 1. The
algorithms call the function with `order` 1 (value only), 2 (value and
first derivative) and 3 (for Newton convergence factors). If $f$ is
built from arithmetic operations and functions of [](arb), the series
versions of [](arb_poly) give the Taylor coefficients, for example
`_arb_poly_sin_series` applied to the series $x + t$ (see
[](arb_calc/isolate_roots)). The enclosures must be valid for the whole
ball `inp`, which in root isolation is an entire subinterval.
:::
:::

::: macro
ARB_CALC_SUCCESS
ARB_CALC_IMPRECISE_INPUT
ARB_CALC_NO_CONVERGENCE

Return codes (0, 1 and 2): success; the input probably needs to be
computed more accurately; the algorithm did not converge, because the
problem has no solution, the algorithm does not apply or the precision
is insufficient.
:::

::: variable
int arb_calc_verbose

If nonzero, the functions of this module print information about the
computation to standard output. The variable is thread-local.
:::

## Intervals {#intervals}

Root isolation works with intervals $[a, b]$ given by exact endpoints,
since a ball cannot represent an interval exactly after bisection.

::: type
arf_interval_struct
arf_interval_t
arf_interval_ptr
arf_interval_srcptr

An `arf_interval_struct` holds two `arf_struct` endpoints `a` and `b`
and represents the interval $[a, b]$. `arf_interval_t` is an array of
one `arf_interval_struct`, and `arf_interval_ptr` and
`arf_interval_srcptr` are pointers used for vectors of intervals.
:::

::: function
void arf_interval_init(arf_interval_t v)
void arf_interval_clear(arf_interval_t v)

Initialise `v` (both endpoints zero), and free its memory.
:::

::: function
arf_interval_ptr _arf_interval_vec_init(slong n)
void _arf_interval_vec_clear(arf_interval_ptr v, slong n)

Allocate and initialise a vector of $n$ intervals, and clear and free
it.
:::

::: function
void arf_interval_set(arf_interval_t v, const arf_interval_t u)
void arf_interval_swap(arf_interval_t v, arf_interval_t u)

Set `v` to a copy of `u`, and swap `v` and `u` efficiently.
:::

::: function
void arf_interval_get_arb(arb_t x, const arf_interval_t v, slong prec)

Sets `x` to a ball that contains the interval $[a, b]$ (using
`arb_set_interval_arf`); the ball is exact if the interval is a point
and is otherwise rounded outward to `prec` bits.
:::

::: function
void arf_interval_printd(const arf_interval_t v, slong n)
void arf_interval_fprintd(FILE * file, const arf_interval_t v, slong n)

Print `v` as `[a, b]` with the endpoints printed to $n$ decimal digits,
to standard output or to `file`.
:::

## Root isolation by bisection {#isolation}

::: function
slong arb_calc_isolate_roots(arf_interval_ptr * found, int ** flags, arb_calc_func_t func, void * param, const arf_interval_t interval, slong maxdepth, slong maxeval, slong maxfound, slong prec)

Isolates the simple roots of a real analytic function $f$ in the
interior of `interval`. Allocates and writes to `*found` an array of $n$
subintervals of `interval` and to `*flags` an array of $n$ flags, and
returns $n$. The caller must free them with `_arf_interval_vec_clear`
and `flint_free`. The output satisfies:

- $f$ has no root in `interval` outside the output subintervals;
- the subintervals are sorted in increasing order and do not overlap,
  except that consecutive ones may share an endpoint;
- a subinterval with flag 1 contains exactly one root, which is simple;
- a subinterval with another flag may or may not contain roots.

If all flags are 1, all roots in the interior of `interval` have been
isolated. Roots of multiplicity greater than one and roots at the
endpoints of the subdivision cannot be isolated; in their place, a
subinterval of minimal size with flag other than 1 is returned. The
other subintervals can be searched again with a higher precision or
larger limits.

::: details
The search stops subdividing after `maxdepth` levels of bisection, so
the smallest subintervals have about $2^{-\mathit{maxdepth}}$ times the
width of `interval` (20 to 50 is typical; `maxdepth` should not exceed
`prec`); after `maxeval` subintervals have been examined, the
remaining ones are returned unexamined, so that the number of calls of
`func` is a small multiple of `maxeval` (100 to 100000 is typical); and
after `maxfound` roots have been isolated (pass `WORD_MAX` to find all
roots, 1 to find one). `prec` is the precision used to evaluate $f$.
The bisection points must be exactly representable, which holds unless
the endpoints have enormous exponents.

A subinterval $I$ is discarded if $f(I)$ excludes zero. It is reported
as containing exactly one root if $f$ has opposite nonzero signs at its
endpoints and $f'(I)$ excludes zero ($f$ is then strictly monotonic on
$I$). Otherwise it is bisected at its midpoint (`arb_calc_partition`)
([](algorithms/ball-arithmetic/real-root-isolation)).
:::
:::

::: function
int arb_calc_partition(arf_interval_t L, arf_interval_t R, arb_calc_func_t func, void * param, const arf_interval_t block, slong prec)

Splits `block` $= [a, b]$ at its exact midpoint $m$ into $L = [a, m]$
and $R = [m, b]$ and returns the sign of $f(m)$: $\pm 1$ if it is
determined, and 0 if $f(m)$ contains zero.
:::

::: function
int arb_calc_refine_root_bisect(arf_interval_t r, arb_calc_func_t func, void * param, const arf_interval_t start, slong iter, slong prec)

Given an interval `start` with $f$ of opposite nonzero signs at its
endpoints, performs `iter` bisection steps, keeping the half with a sign
change, and sets `r` to the result, which contains a root (the only one
if `start` contains only one). Returns `ARB_CALC_SUCCESS`;
`ARB_CALC_IMPRECISE_INPUT` (with `r` unchanged) if the signs at the
endpoints of `start` are not determined or equal; or
`ARB_CALC_NO_CONVERGENCE` if the sign at a bisection point cannot be
determined, in which case `r` is the last valid interval.
:::

## Newton iteration {#newton}

::: function
void arb_calc_newton_conv_factor(arf_t conv_factor, arb_calc_func_t func, void * param, const arb_t conv_region, slong prec)

Sets `conv_factor` to an upper bound for $C = \sup_{t, u \in I}
\tfrac12 |f''(t)| / |f'(u)|$ on the ball $I$ = `conv_region`, obtained
by evaluating the Taylor coefficients $f'(I)$ and $f''(I)/2$ with one
call to `func` and bounding the absolute value of their quotient. If
$f$ is ill-conditioned, $I$ may need to be extremely narrow for the
bound to be finite.

::: note
In FLINT 3.6 and earlier
the function returned only $C/2$, and `arb_calc_newton_step` and
`arb_calc_refine_root_newton` given this factor could return a ball
that excludes the root. Doubling the result is a safe workaround for
older versions.
:::
:::

::: function
int arb_calc_newton_step(arb_t xnew, arb_calc_func_t func, void * param, const arb_t x, const arb_t conv_region, const arf_t conv_factor, slong prec)

Performs one step of interval Newton iteration for a ball $x = [m \pm
r]$ that contains exactly one root of $f$, given a ball $I \supseteq x$
(`conv_region`) and the bound $C$ for $\tfrac12 |f''| / |f'|$ on $I$
(`conv_factor`). Computes $m' = m - f(m)/f'(m)$ in ball arithmetic at
`prec` bits and adds $C r^2$ to its radius. If the new ball $x'$ lies in
$I$ and has a smaller radius than $x$, sets `xnew` to $x'$ and returns
`ARB_CALC_SUCCESS`; otherwise sets `xnew` to $x$ and returns
`ARB_CALC_NO_CONVERGENCE`.

::: details
For the root $\xi \in x$, Taylor's theorem gives $m' - \xi =
\frac{f''(\eta)}{2 f'(m)} (m - \xi)^2$ for some $\eta$ between $m$ and
$\xi$, so $|m' - \xi| \le C r^2$, and the rounding errors in $m'$ are
included in the ball arithmetic. This is only valid if $C$ really
bounds $\tfrac12 |f''|/|f'|$ on $I$.
:::
:::

::: function
int arb_calc_refine_root_newton(arb_t r, arb_calc_func_t func, void * param, const arb_t start, const arb_t conv_region, const arf_t conv_factor, slong eval_extra_prec, slong prec)

Refines the enclosure `start` of a simple root to a relative accuracy
of about `prec` bits by Newton steps at working precisions that roughly
double from the accuracy of `start` up to `prec` (plus a padding of
about $\log_2 C$ bits). `eval_extra_prec` is the number of guard bits
needed to evaluate $f$ accurately near the root (for a polynomial with
large coefficients of alternating signs evaluated by Horner's rule,
about the bit size of the coefficients). Returns `ARB_CALC_SUCCESS` if
every step succeeded; this does not guarantee `prec` accurate bits,
only an improved enclosure. On failure returns
`ARB_CALC_IMPRECISE_INPUT` or `ARB_CALC_NO_CONVERGENCE`, and `r` is a
valid but possibly less accurate enclosure (possibly `start`).
:::
