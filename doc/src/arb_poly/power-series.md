---
title: Power series and Taylor expansions
summary: How truncated power series with ball coefficients are used to compute rigorous Taylor coefficients and derivatives of composite functions, with the conventions for lengths, precision and aliasing of the series functions.
order: 1
---

An `arb_poly_t` of length at most $n$ also represents a power series
modulo $x^n$. The functions whose name ends in `_series` take a target
length $n$ and return the first $n$ coefficients of the power series of
their result: `arb_poly_mullow(c, a, b, n, prec)` gives $ab \bmod x^n$,
`arb_poly_exp_series(f, h, n, prec)` gives $\exp(h) \bmod x^n$, and so
on. Since the first $n$ coefficients of $F(h)$ depend only on the first
$n$ coefficients of $h$, these operations compose: a chain of series
operations computes the Taylor coefficients of the composite function
with rigorous error bounds. This is forward-mode automatic
differentiation of any order, and it is how FLINT computes derivatives
of special functions. The algorithms are described in
[](algorithms/ball-arithmetic/power-series).

## Derivatives of a composite function {#derivatives}

To expand $f(x)$ at a point $x_0$, start from the series $x_0 + t$ of
length 2 and apply the series versions of the operations that define
$f$. The coefficient of $t^k$ of the result is $f^{(k)}(x_0)/k!$. The
program computes the first eight derivatives of
$$f(x) = \frac{\exp(\sin x)}{1 + x^2}$$
at $x_0 = 1/2$:

```c example name=arb_poly_series_derivatives
#include <flint/arb_poly.h>

int main(void)
{
    arb_poly_t x, s, e, d, f;
    arb_t t;
    slong n = 8, prec = 64, k;

    arb_poly_init(x);
    arb_poly_init(s);
    arb_poly_init(e);
    arb_poly_init(d);
    arb_poly_init(f);
    arb_init(t);

    /* the series x0 + t with x0 = 1/2 */
    arb_set_d(t, 0.5);
    arb_poly_set_coeff_arb(x, 0, t);
    arb_poly_set_coeff_si(x, 1, 1);

    arb_poly_sin_series(s, x, n, prec);         /* sin(x0 + t) */
    arb_poly_exp_series(e, s, n, prec);         /* exp(sin(x0 + t)) */
    arb_poly_mullow(d, x, x, n, prec);          /* (x0 + t)^2 */
    arb_poly_add_si(d, d, 1, prec);             /* 1 + (x0 + t)^2 */
    arb_poly_div_series(f, e, d, n, prec);

    for (k = 0; k < n; k++)
    {
        /* f^(k)(x0) = k! * coefficient k */
        arb_fac_ui(t, k, prec);
        arb_mul(t, t, f->coeffs + k, prec);
        flint_printf("f^(%wd)(1/2) = ", k);
        arb_printn(t, 15, 0);
        flint_printf("\n");
    }

    arb_poly_clear(x);
    arb_poly_clear(s);
    arb_poly_clear(e);
    arb_poly_clear(d);
    arb_poly_clear(f);
    arb_clear(t);
    return 0;
}
```
```output
f^(0)(1/2) = [1.29211703715367 +/- 3.01e-15]
f^(1)(1/2) = [0.100245750004579 +/- 3.93e-16]
f^(2)(1/2) = [-1.85212894004061 +/- 3.03e-15]
f^(3)(1/2) = [2.07237670339131 +/- 4.52e-15]
f^(4)(1/2) = [6.58264514459575 +/- 3.32e-15]
f^(5)(1/2) = [-58.5378725284893 +/- 5.01e-14]
f^(6)(1/2) = [157.535450408399 +/- 1.89e-13]
f^(7)(1/2) = [1154.76846929024 +/- 3.86e-12]
```

The cost of each operation is that of a few multiplications of length
$n$ (the elementary functions use Newton iteration, $O(M(n))$), so
derivatives of high order are cheap. The radii grow
roughly in proportion to the coefficients, so that the relative accuracy
decreases only slowly with the order: the derivative of order 999 of
this function, about $-8.63 \cdot 10^{2515}$, still has 13 correct
digits.

If $x_0$ is a ball rather than a point, the coefficients enclose the
Taylor coefficients at every point of the ball. Combined with the
`arb_poly_majorant` of the coefficients, this gives bounds for $f$ and
its derivatives on a whole interval, which is the basis of the
root-finding and integration algorithms in [](arb_calc) and
[](acb_calc).

## Lengths, precision and aliasing {#conventions}

- The length argument of a series function is the number of
  coefficients to compute. Inputs may be longer (they are truncated) or
  shorter (they are padded with zeros). For the underscore versions the
  output array must have room for exactly that many coefficients, and
  both lengths must be positive.
- A series function is only as good as the function at the constant
  term. `arb_poly_log_series` needs a constant term that excludes zero,
  `arb_poly_inv_series` and `arb_poly_div_series` a divisor whose
  constant term excludes zero, `arb_poly_sqrt_series` a positive
  constant term, and so on; otherwise the coefficients are not finite.
  `arb_poly_compose_series` and `arb_poly_revert_series` need an exactly
  zero constant term in the inner series.
- The output can be the same object as an input for all functions
  without underscore. Many underscore functions also allow aliasing;
  the overview states for each function whether it does.
- Functions that take two series outputs, such as
  `arb_poly_sin_cos_series`, compute both at the cost of little more
  than one. Prefer them when both values are needed.
- Squaring is detected when the same object is passed twice to
  `arb_poly_mullow`; this is faster and gives tighter enclosures than
  multiplying two copies.

## Special functions and inverse functions {#special}

Many special functions accept series arguments directly: in this module
the gamma, digamma, zeta, Lambert W and Riemann–Siegel functions, and
in [](arb_hypgeom) the error functions, exponential and trigonometric
integrals, Bessel and Airy functions and others. The compositional
inverse of a series is computed by `arb_poly_revert_series`. The
following program expands $\Gamma(1 + x)$, whose coefficients are
$1$, $-\gamma$, $(\gamma^2 + \zeta(2))/2$, ..., and computes the series
of the Lambert W function as the inverse of $x e^x$; its coefficients
are $(-k)^{k-1}/k!$:

```c example name=arb_poly_series_special
#include <flint/arb_poly.h>

int main(void)
{
    arb_poly_t h, g, w;
    slong k, n = 6, prec = 64;

    arb_poly_init(h);
    arb_poly_init(g);
    arb_poly_init(w);

    /* Gamma(1 + x) */
    arb_poly_set_coeff_si(h, 0, 1);
    arb_poly_set_coeff_si(h, 1, 1);
    arb_poly_gamma_series(g, h, n, prec);
    flint_printf("Gamma(1+x) = ");
    arb_poly_printd(g, 10);
    flint_printf("\n");

    /* W(x) is the inverse of x exp(x) */
    arb_poly_zero(h);
    arb_poly_set_coeff_si(h, 1, 1);           /* x */
    arb_poly_exp_series(g, h, n, prec);       /* exp(x) */
    arb_poly_shift_left(g, g, 1);             /* x exp(x) */
    arb_poly_revert_series(w, g, n, prec);
    for (k = 1; k < n; k++)
    {
        flint_printf("W: x^%wd: ", k);
        arb_printn(w->coeffs + k, 10, 0);
        flint_printf("\n");
    }

    arb_poly_clear(h);
    arb_poly_clear(g);
    arb_poly_clear(w);
    return 0;
}
```
```output
Gamma(1+x) = [(1.000000000 +/- 0), (-0.5772156649 +/- 5.4422e-20), (0.9890559953 +/- 5.5302e-20), (-0.9074790761 +/- 5.5981e-20), (0.9817280868 +/- 5.6223e-20), (-0.9819950689 +/- 2.8364e-21)]
W: x^1: 1.000000000
W: x^2: -1.000000000
W: x^3: 1.500000000
W: x^4: [-2.666666667 +/- 3.34e-10]
W: x^5: [5.208333333 +/- 3.34e-10]
```

The radii $3.34 \cdot 10^{-10}$ printed for $W$ come from the
conversion of $-8/3$ and $125/24$ to ten digits; the computed balls are
accurate to about 60 bits. The constant term of $x e^x$ is exactly zero
here, as `arb_poly_revert_series` requires: `arb_poly_shift_left`
inserts an exact zero.

For complex arguments, [](acb_poly) provides the same functions, and
the guide [](acb_poly/taylor-series) continues with the zeta function
and polylogarithms.
