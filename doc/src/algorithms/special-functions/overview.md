---
title: Evaluating special functions
kind: topic
summary: How FLINT evaluates special functions with rigorous error bounds — the common principles (enclosures, error propagation, cancellation-aware working precision, argument reduction, convergent and asymptotic expansions, analytic continuation, caching) and a map of the methods to the functions and modules that use them.
order: 0
---

FLINT's special functions take balls as input (`arb_t`, `acb_t`, or
power series with ball coefficients) and return balls that contain the
exact value of the function at every point of the input. Every
algorithm therefore comes with an error bound: for a truncated series,
for the propagation of the input radius, and for the rounding errors,
which ball arithmetic tracks automatically. Most of the code goes back
to Arb [@Joh2017; @Joh2014a].

## Principles {#principles}

**Enclosures, not approximations.** A result $[m \pm r]$ is correct if
it contains the true value; it is useful if $r$ is small. The functions
are written so that $r$ is close to the best possible for exact
inputs at the target precision, but they never return an incorrect
enclosure to obtain a smaller one. When no finite enclosure can be
given — a pole in the input, an input that straddles a branch cut, a
missing algorithm, or too low a precision for the size of the
arguments — the result is indeterminate (infinite radius). The
semantics and the cost model of ball arithmetic are explained in
[](guide/types/ball-arithmetic) and [](algorithms/ball-arithmetic/overview).

**Truncation bounds.** Every series is truncated with an explicit bound
for the remainder: tail bounds by geometric series for convergent
hypergeometric series ([](algorithms/special-functions/hypergeometric#tail-bound)),
Olver-type bounds for asymptotic expansions
([](algorithms/special-functions/hypergeometric#asymptotic),
[](algorithms/special-functions/gamma#stirling)), Cauchy majorants for
Taylor series of solutions of differential equations
([](algorithms/special-functions/gauss-2f1#corner)), and
Euler–Maclaurin remainder bounds for zeta functions
([](algorithms/special-functions/zeta)). The bounds are computed with
the `mag_t` type, which rounds upward, so that bounding costs little
compared to the evaluation.

**Working precision.** Ball arithmetic propagates rounding errors
correctly but does not prevent cancellation: summing $\sum (-x)^k/k!$
for $x = 100$ at 64 bits gives a correct but useless ball. The
functions therefore estimate the cancellation of the method they use,
usually in double precision, and add that many guard bits; conversely
they lower the working precision to the accuracy of the input, since
more precision cannot improve the result.

**Midpoint evaluation and propagated error.** For an input ball $x = [m
\pm r]$, evaluating an algorithm directly in ball arithmetic often
overestimates the result (each operation widens the ball, and
correlations are lost). Many functions instead evaluate at the exact
midpoint $m$ and add $r \cdot \sup_{|t - m| \le r} |f'(t)|$, using a
cheap bound for the derivative: $e^{-z^2}$ for the error function,
$\psi(x)$ for the gamma function, Bessel and Airy function bounds,
Legendre polynomial derivative bounds. This also allows the evaluation
at the midpoint to use a lower precision.

**Choosing among expansions.** Most functions have several
representations — a convergent series at the origin, an asymptotic
expansion at infinity, transformations of the argument, integral
representations — with very different costs depending on the argument,
the parameters and the precision. FLINT chooses among them with cheap
floating-point estimates of the number of terms and of the
cancellation, and some functions retry with another method if the
result is inaccurate.

**Argument reduction and functional equations.** Recurrences and
symmetries move the argument to a region where an expansion converges
quickly: $\Gamma(z) = \Gamma(z+r)/(z)_r$ and the reflection formula,
the linear fractional transformations of ${}_2F_1$ and the dilogarithm,
the modular group for theta functions and modular forms, $\exp(x) =
\exp(x/2^r)^{2^r}$ and table lookups for the elementary functions.

**Summation algorithms.** The cost of a series is dominated by the
multiplications involving full-precision numbers. Rectangular
splitting reduces their number to $O(\sqrt N)$ for $N$ terms when the
term ratios have small coefficients; binary splitting computes series
with short arguments and rational coefficients in quasi-linear time;
the bit-burst algorithm extends this to arbitrary arguments for
holonomic functions ([](algorithms/special-functions/hypergeometric#summation)).

**Caching.** Constants ($\pi$, $\log 2$, Euler's $\gamma$,
$\Gamma(1/3)$, ...) and Bernoulli numbers are computed once at the
highest precision requested so far and cached per thread or globally
([](algorithms/special-functions/constants)).

**Power series.** Functions of power series arguments (`_series`) and
Taylor expansions at a point (`_jet`) give derivatives with rigorous
bounds. They are computed either from the differential equation of the
function (Bessel, Airy, Coulomb functions), from the derivative and
power series integration (error functions, exponential integrals), or
by running the scalar algorithm with power series arithmetic (gamma
function, hypergeometric series with power series parameters, which is
also how limits at removable singularities are computed).

## Map of functions and methods {#map}

| Functions | Main methods | Modules | Algorithm page |
|---|---|---|---|
| $\exp$, $\log$, trigonometric and hyperbolic functions and their inverses | argument reduction, tables, rectangular splitting on fixed-point numbers, bit-burst | [`arb`](arb), [`acb`](acb) | [](algorithms/special-functions/elementary) |
| $\pi$, $e$, $\log 2$, $\gamma$, $\zeta(3)$, Catalan, ... | hypergeometric series by binary splitting, caching | [`arb`](arb), [`mp_real`](mp_real) | [](algorithms/special-functions/constants) |
| $\Gamma$, $1/\Gamma$, $\log\Gamma$, rising factorials, digamma | Stirling series, Taylor series of $1/\Gamma$, rational arguments, rectangular and binary splitting | [`arb_hypgeom`](arb_hypgeom), [`acb_hypgeom`](acb_hypgeom), [`arb`](arb), [`acb`](acb) | [](algorithms/special-functions/gamma) |
| ${}_pF_q$, ${}_0F_1$, ${}_1F_1$, $U$ | convergent series with tail bounds, asymptotic expansion of $U$, connection formulas, numerical integration | [`acb_hypgeom`](acb_hypgeom) | [](algorithms/special-functions/hypergeometric) |
| sums with rational parameters | forward recurrence, rectangular splitting, binary splitting | [`arb_hypgeom`](arb_hypgeom) | [](algorithms/special-functions/hypergeometric-sums) |
| $\operatorname{erf}$, $\operatorname{erfc}$, Fresnel, $\Gamma(s,z)$, $\gamma(s,z)$, $B(a,b;z)$, $E_s$, $\operatorname{Ei}$, $\operatorname{Si}$, $\operatorname{Ci}$, $\operatorname{li}$ | ${}_1F_1$ and $U$ representations, real series with rectangular splitting, bit-burst for $\Gamma(a, z)$ | [`arb_hypgeom`](arb_hypgeom), [`acb_hypgeom`](acb_hypgeom) | [](algorithms/special-functions/error-functions) |
| Bessel $J, Y, I, K$, Airy, Coulomb wave functions | ${}_0F_1$ series, asymptotic expansions of $U$, dedicated Airy bounds, interval Newton for zeros | [`arb_hypgeom`](arb_hypgeom), [`acb_hypgeom`](acb_hypgeom) | [](algorithms/special-functions/bessel-airy) |
| ${}_2F_1$, incomplete beta | transformations, limits by power series, analytic continuation, integration | [`acb_hypgeom`](acb_hypgeom) | [](algorithms/special-functions/gauss-2f1) |
| Legendre polynomials and functions, Gauss–Legendre nodes, orthogonal polynomials, spherical harmonics | series at 0 and 1, asymptotic expansion, recurrences, ${}_2F_1$ | [`arb_hypgeom`](arb_hypgeom), [`acb_hypgeom`](acb_hypgeom) | [](algorithms/special-functions/legendre) |
| dilogarithm | transformations, Bernoulli series, Taylor continuation, bit-burst | [`acb_hypgeom`](acb_hypgeom) | [](algorithms/special-functions/dilogarithm) |
| Riemann and Hurwitz zeta functions, polylogarithms, Dirichlet $L$-functions, zeta zeros | Euler–Maclaurin summation, Riemann–Siegel, Euler products, Borwein's algorithm | [`arb`](arb), [`acb`](acb), [`acb_dirichlet`](acb_dirichlet) | [](algorithms/special-functions/zeta), [](algorithms/special-functions/zeta-zeros) |
| elliptic integrals and functions, modular forms, theta functions | AGM, Carlson's algorithm, reduction to the fundamental domain, $q$-series | [`acb_elliptic`](acb_elliptic), [`acb_modular`](acb_modular), [`acb_theta`](acb_theta) | [](algorithms/special-functions/elliptic-modular), [](algorithms/special-functions/theta-functions), [](algorithms/special-functions/theta-quasi-linear) |
| Lambert $W$ | Halley iteration from asymptotic and series approximations, a posteriori error bounds and a branch check | [`arb`](arb), [`acb`](acb) | [](algorithms/special-functions/lambert-w) |
| Bernoulli numbers, partitions | zeta function or multimodular algorithms; Hardy–Ramanujan–Rademacher formula | [`bernoulli`](bernoulli), [`partitions`](partitions) | [](algorithms/number-theory/bernoulli), [](algorithms/number-theory/partitions) |

The power series versions in [`arb_poly`](arb_poly) and
[`acb_poly`](acb_poly) (`arb_poly_gamma_series`,
`acb_poly_zeta_series`, ...) use the same algorithms with power series
arithmetic.

## Interfaces {#interfaces}

The functions are available at three levels.

- **Ball arithmetic.** The modules in the table, with `arb_t` and
  `acb_t` arguments, are the primary interface. The real functions of
  [`arb_hypgeom`](arb_hypgeom) mostly wrap the complex ones and return
  an indeterminate result when the value is not real.
- **Generic rings.** [`gr_special`](gr_special) defines `gr_gamma`,
  `gr_erf`, `gr_hypgeom_2f1`, `gr_bessel_j`, `gr_zeta` and many more for
  every ring; the contexts of real and complex balls call the functions
  above, the floating-point contexts of [`nfloat`](nfloat) and
  [`mp_real`](mp_real) have their own implementations of the
  elementary functions, and exact contexts such as
  [`ca`](ca) return exact values where they can (for example
  $\Gamma(1/2) = \sqrt\pi$). Some power series algorithms exist only in
  generic form, such as `_gr_poly_bessel_j_series`.
- **Floating-point wrappers.** [`arb_fpwrap`](arb_fpwrap) evaluates
  the functions for `double` and `complex double` arguments with
  correct rounding or a guaranteed relative accuracy, increasing the
  precision until the ball determines the result.

## Use in FLINT {#use}

The module pages give the conventions (precision, branch cuts,
regularization, optional outputs); the guide [](acb_hypgeom/using)
explains them for the hypergeometric functions with examples.
