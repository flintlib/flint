---
title: acb_calc.h
header: acb_calc.h
description: Rigorous numerical integration of complex functions along straight-line paths, combining direct interval enclosures, Gauss–Legendre quadrature with error bounds from analyticity, and adaptive subdivision.
category: Real and complex numbers
order: 12
---

This module computes rigorous enclosures of integrals $\int_a^b f(t)\,
dt$ of complex functions along the segment from $a$ to $b$ in the
complex plane. The function is given by an evaluator that encloses $f$
on a complex ball and reports when $f$ is not holomorphic there. The
main function `acb_calc_integrate` bisects the path adaptively; on each
segment it uses either the trivial enclosure $(b - a) f([a, b])$ or a
Gauss–Legendre quadrature rule whose error is bounded through a bound
for $|f|$ on an ellipse around the segment. This gives exponential
convergence for piecewise holomorphic integrands and correct results
for integrands with discontinuities, branch points and singularities
near the path [@Joh2018a]. The guide [](acb_calc/integrate) explains
how to write integrands and choose the parameters,
[](algorithms/ball-arithmetic/integration) gives the algorithm and its
error bounds, and the tutorial [](tutorials/numerics/integrals) computes
several integrals in complete programs.

::: summary algorithms/ball-arithmetic/integration :::

## Integrands {#integrands}

::: type
acb_calc_func_t

The type of an integrand, a pointer to a function

```c
int func(acb_ptr out, const acb_t inp, void * param, slong order, slong prec)
```

which computes $f$ on the complex ball `inp` at a precision of about
`prec` bits; `param` passes user data through. The return value is
reserved for an error code and should be 0. `out` and `inp` are not
aliased.

::: details
With `order` = 0, `func` writes an enclosure of $f(z)$ for all $z$ in
`inp` to `out[0]`. Here $f$ may be any function, also discontinuous or
non-holomorphic. With `order` = $n \ge 1$, `func` writes the first $n$
Taylor coefficients of $f$ at `inp`, and *must verify that $f$ is
holomorphic on `inp`*, writing non-finite values otherwise. The
integration algorithms only use `order` 0 and 1, so derivatives are not
needed: with `order` 1, `func` computes $f$ and checks holomorphy.

Integrands built from field operations and meromorphic functions need
no special care, since division by zero and poles give non-finite
enclosures automatically. Functions with branch cuts must detect when
`inp` meets the cut; the functions `acb_sqrt_analytic`,
`acb_rsqrt_analytic`, `acb_log_analytic` and `acb_pow_analytic` do this
when their `analytic` flag is set, and `acb_real_abs`, `acb_real_sgn`,
`acb_real_heaviside`, `acb_real_floor`, `acb_real_ceil`, `acb_real_max`,
`acb_real_min` and `acb_real_sqrtpos` give piecewise holomorphic
extensions of real functions. For example, a piecewise holomorphic
extension of $\sqrt{\lfloor z \rfloor}$ is

```c
int func(acb_ptr out, const acb_t inp, void * param, slong order, slong prec)
{
    if (order > 1)
        flint_abort();  /* derivatives not implemented */

    acb_real_floor(out, inp, order != 0, prec);
    acb_sqrt_analytic(out, out, order != 0, prec);
    return 0;
}
```

and a custom function with branch cut detection follows the pattern

```c
/* Square root with detection of the branch cut (-inf, 0]. */
void sqrt_analytic(acb_t res, const acb_t z, int analytic, slong prec)
{
    if (analytic && arb_contains_zero(acb_imagref(z))
                 && arb_contains_nonpositive(acb_realref(z)))
        acb_indeterminate(res);
    else
        acb_sqrt(res, z, prec);
}
```

The example program `examples/integrals.c` contains many more
integrands.
:::
:::

## Integration {#integration}

::: function
int acb_calc_integrate(acb_t res, acb_calc_func_t func, void * param, const acb_t a, const acb_t b, slong rel_goal, const mag_t abs_tol, const acb_calc_integrate_opt_t options, slong prec)

Sets `res` to an enclosure of $\int_a^b f(t)\, dt$ along the straight
line from $a$ to $b$, where $f$ is given by `func` and `param`. Returns
`ARB_CALC_SUCCESS` if the integration converged to the requested
accuracy on all segments and `ARB_CALC_NO_CONVERGENCE` otherwise; in
both cases `res` is a valid enclosure, which is finite only if $a$ and
$b$ are finite and $f$ is bounded on the path. Improper integrals must
be truncated or transformed by the caller.

::: details
`rel_goal` $r \ge 0$ and `abs_tol` $\varepsilon_{\mathrm{abs}}$ set the
target error $\max(\varepsilon_{\mathrm{abs}}, M 2^{-r})$ on each
segment, where $M$ is the magnitude of the integral, estimated during
the computation; `prec` is the working precision for evaluating $f$ and
for the endpoints, which the algorithm does not change by itself. The
usual choice is $r = \mathit{prec}$ and $\varepsilon_{\mathrm{abs}} =
2^{-\mathit{prec}}$. The total error may exceed both goals because of
the number of segments, cancellation between them and rounding errors.
To compute an integral much smaller than 1 to high relative accuracy,
set $\varepsilon_{\mathrm{abs}}$ to about $M 2^{-r}$ for an estimate
$M$, or to 0 (then only the relative goal applies, which may need extra
subdivisions to estimate $M$). The absolute tolerance is raised
automatically when the partial integrals turn out to be large.
`options` may be `NULL` for the defaults. See [](acb_calc/integrate)
for the algorithm and the parameters.
:::
:::

::: type
acb_calc_integrate_opt_struct
acb_calc_integrate_opt_t

Options for `acb_calc_integrate`. `acb_calc_integrate_opt_t` is an
array of one `acb_calc_integrate_opt_struct`. The fields are:

- `slong deg_limit`: the largest Gauss–Legendre degree on a segment;
  default (if 0 or negative) $0.5 \min(\mathit{prec},
  \mathit{rel\_goal}) + 60$;
- `slong eval_limit`: the approximate maximal number of evaluations of
  $f$; default $1000\,\mathit{prec} + \mathit{prec}^2$;
- `slong depth_limit`: the maximal number of segments waiting to be
  processed (memory use is proportional to it); default $2\,\mathit{prec}$;
- `int use_heap`: if 1, process the segment with the largest error
  first (a priority queue) instead of the most recently created one (a
  stack);
- `int verbose`: 1 prints a summary, 2 also information on each
  segment.

::: details
A higher `deg_limit` can help for functions that are holomorphic in a
large neighbourhood of the path but oscillate, such as entire functions
with many zeros; the precomputation of the nodes grows with the degree.
A lower `eval_limit` makes the function give up sooner when it
converges slowly; integrands with many discontinuities or singularities
close to the path may need a higher one. The heap order sometimes gives
better results when convergence fails, but may need many more queued
segments (a larger `depth_limit`).
:::
:::

::: function
void acb_calc_integrate_opt_init(acb_calc_integrate_opt_t options)

Initialises `options` with all fields 0, which selects the defaults.
:::

## Local integration algorithms {#local}

::: function
int acb_calc_integrate_gl_auto_deg(acb_t res, slong * num_eval, acb_calc_func_t func, void * param, const acb_t a, const acb_t b, const mag_t tol, slong deg_limit, int verbose, slong prec)

Attempts to compute $\int_a^b f(t)\, dt$ with one Gauss–Legendre rule
of degree at most `deg_limit`, chosen automatically so that the error
bound is below `tol`. Returns `ARB_CALC_SUCCESS` and the enclosure in
`res` if this is possible, and `ARB_CALC_NO_CONVERGENCE` with `res`
indeterminate otherwise. The number of evaluations of $f$ is written to
`num_eval`. If `verbose` is nonzero, the chosen degree is printed.

::: details
With $\Delta = (b - a)/2$ and $m = (a + b)/2$, the integral is
$\int_{-1}^{1} g$ with $g(t) = \Delta f(\Delta t + m)$. If $g$ is
holomorphic and $|g| \le M$ inside the ellipse with foci $\pm 1$ and
semi-axes $X$ and $Y = \sqrt{X^2 - 1}$, $\rho = X + Y$, the error of the
$n$-point Gauss–Legendre rule is at most
$$\frac{64 M}{15 (\rho - 1) \rho^{2n - 1}},$$
a slightly weakened form of a bound of Trefethen [@Tre2008]. $M$ is
obtained by evaluating $f$ with `order` 1 on the rectangle $[\pm X] +
[\pm Y] i$ mapped to the segment. The function tries $X = 2, 4, 8, 32,
512, \ldots$ (doubling the exponent) as long as $M$ is finite, and for
each $X$ the smallest degree $n$ from the sequence $1, 2, 4, 6, 8, 12,
16, 22, 32, \ldots$ (steps of about $\sqrt{2}$) that meets the
tolerance; the smallest $n$ found is used. The nodes and weights are
computed by `arb_hypgeom_legendre_p_ui_root` [@JM2018]
([](arb_hypgeom/legendre-roots)) and cached
(per thread), and the evaluations are distributed over threads when
`flint_get_num_threads()` is greater than one.
:::
:::

## Older integration functions {#old}

::: function
void acb_calc_cauchy_bound(arb_t bound, acb_calc_func_t func, void * param, const acb_t x, const arb_t radius, slong maxdepth, slong prec)

Sets `bound` to a ball containing
$$C(x, r) = \frac{1}{2\pi r} \oint_{|z - x| = r} |f(z)|\, |dz| = \int_0^1 |f(x + r e^{2\pi i t})|\, dt,$$
where $r$ is `radius`, by a step sum over $n = 16, 32, 64, \ldots$
arcs of the circle (at most `maxdepth` doublings), stopping when the
enclosure excludes zero. The whole circle is subdivided uniformly,
without adaptivity.
:::

::: function
int acb_calc_integrate_taylor(acb_t res, acb_calc_func_t func, void * param, const acb_t a, const acb_t b, const arf_t inner_radius, const arf_t outer_radius, slong accuracy_goal, slong prec)

Computes $\int_a^b f(t)\, dt$ along the segment from $a$ to $b$ (both
finite) by Taylor polynomials of $f$ centred at uniformly spaced points,
with truncation errors bounded by the Cauchy integral formula. The
singularities of $f$ must lie at distance greater than `outer_radius`
from the path; the steps do not exceed `inner_radius`, which must be
smaller (half of `outer_radius` gives about one bit per term). The
truncation error targets $2^{-\mathit{accuracy\_goal}}$ in absolute
terms; arithmetic is done at `prec` bits. This function uses
derivatives of $f$ (`order` greater than 1) and is superseded by
`acb_calc_integrate`.

::: details
If $f(m + x) = \sum_n a_n x^n$, the integral of the tail satisfies
$$\left|\sum_{n=N}^{\infty} a_n \frac{x^{n+1}}{n+1}\right| \le
\frac{C(m, R)\, R\, |x|}{(R - |x|)(N + 1)} \left(\frac{|x|}{R}\right)^N$$
for $|x| < R$, with $C(m, R)$ from `acb_calc_cauchy_bound`.
:::
:::
