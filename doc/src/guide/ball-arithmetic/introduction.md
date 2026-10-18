---
title: What ball arithmetic is
summary: Midpoint-radius arithmetic as automatic error analysis, how it compares with floating-point and interval arithmetic, and the modules that implement it in FLINT.
order: 1
---

Ball arithmetic, also called midpoint-radius interval arithmetic, extends
floating-point arithmetic by attaching an error bound to every number. A
real number is represented by a ball $[m \pm r]$, the set of real numbers
$t$ with $|t - m| \le r$, where the midpoint $m$ is a binary
floating-point number of arbitrary precision and the radius $r$ is a
nonnegative floating-point number of low precision. Every operation
returns a ball that contains all possible exact results for inputs in the
input balls. A computation in ball arithmetic therefore produces a
mathematically rigorous enclosure of the exact result, with no error
analysis by the user.

## Why balls {#why}

With plain floating-point arithmetic, the result of a computation is an
approximation with an unknown error. To know how many digits are correct,
one has to carry out an error analysis by hand, which is laborious and
error-prone, or rely on heuristics such as repeating the computation at
a higher precision and comparing the results, which can fail. Ball
arithmetic makes the error analysis automatic: the radius of the output
is a proof of the accuracy of the midpoint.

The classical alternative is inf-sup interval arithmetic, which stores
both endpoints $[a, b]$ of an interval as full-precision numbers and so
costs at least twice as much as floating-point arithmetic. In ball
arithmetic only the midpoint has full precision; the radius needs only a
few bits, and FLINT stores it as a 30-bit `mag_t` whose arithmetic is
much cheaper than the arithmetic on the midpoint. At high precision a
ball operation therefore costs about the same as the corresponding
floating-point operation. A ball is also more natural than an interval for
numbers that are known to many digits, and the mathematics of error
propagation is often simpler in the midpoint-radius form. Van der
Hoeven's report [@Hoe2009] is a good introduction to ball arithmetic,
the paper [@Joh2017] describes the design of the implementation in
FLINT (then the separate library Arb), and Rump's survey [@Rum2010] covers rigorous
numerical methods built on interval arithmetic in general.

Other implementations of ball arithmetic include iRRAM and Mathemagix.
FLINT's implementation differs in technical aspects that make certain
computations faster, in giving full access to the internals through a
low-level interface, and in covering a large part of the special functions
of standard references such as the DLMF [@NIST2012], with rigorous error
bounds.

## What FLINT provides {#modules}

The ball arithmetic part of FLINT consists of the following modules.

- [mag](mag): unsigned floating-point numbers with a 30-bit mantissa and
  an arbitrary exponent, used for radii and other error bounds; every
  operation rounds in a fixed direction.
- [arf](arf): arbitrary-precision binary floating-point numbers with
  correct rounding and unbounded exponents, used for midpoints.
- [arb](arb): real balls, with arithmetic, elementary functions,
  constants and some special functions.
- [acb](acb): complex numbers as pairs of real balls (a rectangle in the
  complex plane), with the corresponding functions.
- [acf](acf): complex floating-point numbers, used by approximate
  algorithms.
- [arb_poly](arb_poly) and [acb_poly](acb_poly): polynomials and power
  series with ball coefficients, with fast multiplication, composition,
  root finding and power series of special functions.
- [arb_mat](arb_mat) and [acb_mat](acb_mat): matrices of balls, with
  rigorous linear solving, determinants, eigenvalues and matrix
  exponentials.
- [arb_hypgeom](arb_hypgeom) and [acb_hypgeom](acb_hypgeom):
  hypergeometric and related special functions (gamma, error functions,
  Bessel functions, orthogonal polynomials, and many more).
- [acb_elliptic](acb_elliptic), [acb_modular](acb_modular),
  [acb_theta](acb_theta), [acb_dirichlet](acb_dirichlet): elliptic and
  modular functions, theta functions, Dirichlet L-functions and the
  Riemann zeta function.
- [arb_calc](arb_calc), [acb_calc](acb_calc): root isolation of real
  functions and rigorous numerical integration.
- [arb_fmpz_poly](arb_fmpz_poly): root isolation and evaluation of
  integer polynomials using balls, and [arb_fpwrap](arb_fpwrap):
  correctly rounded `double` special functions built on ball arithmetic.

The balls are also available as generic rings, see
[gr_domains](gr_domains), which is how matrices and polynomials over
$\mathbb{R}$ and $\mathbb{C}$ in the generic modules are implemented.

## Where to go next {#next}

[](guide/ball-arithmetic/using) introduces the programming interface
through examples. [](guide/ball-arithmetic/precision) explains how the
precision argument relates to the accuracy of results.
[](guide/ball-arithmetic/pitfalls) describes the conventions specific
to ball arithmetic: precision limits, aliasing and algebraic
simplification, cached constants and the use of hardware floating
point.
