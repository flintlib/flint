---
title: Ball arithmetic
kind: topic
summary: How FLINT represents real and complex numbers as midpoint-radius balls, how each operation computes a midpoint and a bound for the propagated and rounding errors, why the radius has its own fast 30-bit type, and which pages describe the algorithms.
order: 0
---

Ball arithmetic represents a real number by a ball $[m \pm r]$ and
replaces each operation $f$ by a function that returns a ball containing
$f(t)$ for every $t$ in the input ball. FLINT implements it in the
modules [mag](mag), [arf](arf), [arb](arb) and [acb](acb), following the
design of Arb [@Joh2017]; the general theory is in [@Hoe2009].

## Midpoint and radius {#representation}

An `arb_t` has an `arf_t` midpoint, a binary floating-point number with
a mantissa of arbitrary length and an `fmpz` exponent, and a `mag_t`
radius, an unsigned floating-point number with a 30-bit mantissa and an
`fmpz` exponent. An `acb_t` is a pair of `arb_t` for the real and
imaginary parts, so it represents a rectangle in the complex plane; disks
would be a little tighter for some operations, but rectangles make the
real and imaginary parts available directly and reuse all of the real
code.

The asymmetry between midpoint and radius is the main point of the
design. The midpoint needs arbitrary precision and correct rounding, the
radius only an upper bound with a few significant bits. With a 30-bit
mantissa, all radius operations use single-word integer arithmetic (with
an upward correction instead of a test for exactness), and exponents are
almost always inline `fmpz` values, so an operation on radii costs a few
nanoseconds. At high precision the cost of a ball operation is therefore
that of the midpoint operation, and even at 64 bits it is only a small
factor more. Inf-sup interval arithmetic, by contrast, needs two
full-precision operations with opposite rounding for every operation.
[](algorithms/ball-arithmetic/representation) describes the formats and
proves that the basic radius operations round upward.

## How an operation is evaluated {#evaluation}

For an operation $f$ applied to $[m \pm r]$, FLINT computes

1. an approximation $\tilde y$ of $f(m)$, rounded to `prec` bits, with a
   bound $\varepsilon$ for its error;
2. a bound $\rho$ for the propagated error $\sup_{|a| \le r} |f(m + a) -
   f(m)|$;

and returns $[\tilde y \pm (\varepsilon + \rho)]$. For the arithmetic
operations $\tilde y$ is a single `arf` operation, and $\varepsilon$ is
one unit in the last place, added only if the operation was inexact, so
exact input with an exact result gives an exact ball. The propagated
error bounds for arithmetic and elementary functions, with proofs, are described in
[](algorithms/ball-arithmetic/error-propagation).

Two refinements make this efficient. First, the midpoint is not computed
to more than about 30 bits beyond the relative accuracy of the input:
`arb_exp` of a ball whose radius is $2^{-100}$ times its midpoint
computes $e^m$ to about 130 bits even if `prec` is $10^6$, since more bits
would be hidden by $\rho$. Second, when the input ball is wide (relative
accuracy below about 20 bits), the first-order bound for $\rho$
overestimates badly, so the elementary functions and division evaluate the
range of $f$ from the endpoints, using the monotonicity or the extrema of
the function, and convert the resulting interval back to a ball.

Complex operations combine the real ones: a product of rectangles applies
the product bound to the real and imaginary parts (computing the
midpoints with `arf_complex_mul`, which uses three real products at high
precision), and functions such as $\exp(a + bi) = e^a (\cos b + i \sin b)$
reduce to real functions with the errors propagated through the formulas.

## Sums, dot products and series {#sums}

Long sums are where ball arithmetic is most easily slower or less
accurate than necessary. `arb_dot` and `acb_dot` sum many products in one
fixed-point accumulator with a single rounding and compute all radius
contributions with word arithmetic; see
[](algorithms/ball-arithmetic/dot-products). Polynomials and power series
of balls (`arb_poly`, `acb_poly`) are multiplied by scaling the
coefficients to a common exponent, multiplying the midpoints exactly as
integer polynomials and bounding the radius contribution with a separate
product of upper bounds. Infinite series of special functions are
truncated where a tail bound, computed with `mag_t` arithmetic, falls
below the target accuracy; common tail bounds are described in
[](algorithms/ball-arithmetic/error-propagation).

## Precision and the inclusion principle {#precision}

The inclusion principle is guaranteed; the size of the radius is a matter
of quality. For the basic operations the result is within a few units in
the last place of the best possible ball, so the accuracy of a computation
can be predicted from the precision and the conditioning of the problem.
For higher functions the radius can be larger, and convergence as the
precision grows is the practical goal. The user guide
[](guide/types/ball-arithmetic#precision) explains how to choose the precision,
and [](guide/types/ball-arithmetic#conventions) lists technical conventions such as
the treatment of aliased arguments.

## Overview of the algorithm pages {#pages}

| Topic | Page | Main functions |
|---|---|---|
| Formats of radii and midpoints, rounding | [](algorithms/ball-arithmetic/representation) | `mag_mul`, `arf_add`, `arf_mag_add_ulp` |
| Propagated error of operations, tail bounds | [](algorithms/ball-arithmetic/error-propagation) | `arb_mul`, `arb_div`, `arb_exp`, `mag_geom_series` |
| Dot products | [](algorithms/ball-arithmetic/dot-products) | `arb_dot`, `acb_dot` |
| Elementary functions | [](algorithms/special-functions/elementary) | `arb_exp`, `arb_log`, `arb_sin_cos`, `arb_atan` |
| Mathematical constants | [](algorithms/special-functions/constants) | `arb_const_pi`, `arb_const_euler` |
| Other special functions | [](algorithms/special-functions/overview) | `arb_hypgeom_gamma`, `acb_hypgeom_2f1` |
| Polynomial multiplication | [](algorithms/ball-arithmetic/polynomial-multiplication) | `arb_poly_mullow`, `acb_poly_mullow` |
| Power series | [](algorithms/ball-arithmetic/power-series) | `arb_poly_inv_series`, `arb_poly_exp_series` |
| Complex roots of polynomials | [](algorithms/ball-arithmetic/polynomial-roots) | `acb_poly_find_roots`, `arb_fmpz_poly_complex_roots` |
| Real root isolation | [](algorithms/ball-arithmetic/real-root-isolation) | `arb_calc_isolate_roots`, `arb_calc_refine_root_newton` |
| Matrix multiplication | [](algorithms/ball-arithmetic/matrix-multiplication) | `arb_mat_mul`, `acb_mat_mul` |
| Linear systems, inverses, determinants | [](algorithms/ball-arithmetic/linear-systems) | `arb_mat_solve`, `arb_mat_det` |
| Eigenvalues and eigenvectors | [](algorithms/ball-arithmetic/eigenvalues) | `acb_mat_eig_simple`, `acb_mat_eig_multiple` |
| Integration | [](algorithms/ball-arithmetic/integration) | `acb_calc_integrate` |

The other floating-point, ball and interval types of FLINT (`nfloat`,
`mp_real`, `dfloat` and intervals of doubles) are compared with `arb` in
[](algorithms/floating-point/overview).
