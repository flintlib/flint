---
title: arb_fmpz_poly.h
header: arb_fmpz_poly.h
description: Functions for integer polynomials that use real and complex balls: evaluation at balls, certified real and complex root isolation, deflation, and minimal polynomials of Gaussian periods.
category: Real and complex numbers
order: 8
---

This module contains functions for integer polynomials (`fmpz_poly_t`)
that use ball arithmetic, either to produce real or complex numbers or
internally to produce an exact result. The central functions are
`arb_fmpz_poly_real_roots` and `arb_fmpz_poly_complex_roots`, which
compute disjoint enclosures of all real or complex roots of a squarefree
integer polynomial to a requested accuracy, combining exact real root
isolation in [](fmpz_poly) ([](algorithms/polynomials/real-roots-integers)), refinement by interval Newton iteration and
bisection, and the certified numerical root finder of [](acb_poly); see
[](arb_fmpz_poly/roots) and the tutorial
[](tutorials/numerics/polynomial-roots). Conversions between integer polynomials and
ball polynomials are in [](arb_poly) and [](acb_poly)
(`arb_poly_set_fmpz_poly`, `arb_poly_get_unique_fmpz_poly`, ...).

The header declares its functions with the types from `acb_types.h` but
does not include `arb.h` or `acb.h`; programs that create balls must
include those headers too.

::: summary algorithms/ball-arithmetic/polynomial-roots :::

## Evaluation {#evaluation}

::: function
void _arb_fmpz_poly_evaluate_arb_horner(arb_t res, const fmpz * poly, slong len, const arb_t x, slong prec)
void arb_fmpz_poly_evaluate_arb_horner(arb_t res, const fmpz_poly_t poly, const arb_t x, slong prec)
void _arb_fmpz_poly_evaluate_arb_rectangular(arb_t res, const fmpz * poly, slong len, const arb_t x, slong prec)
void arb_fmpz_poly_evaluate_arb_rectangular(arb_t res, const fmpz_poly_t poly, const arb_t x, slong prec)
void _arb_fmpz_poly_evaluate_arb(arb_t res, const fmpz * poly, slong len, const arb_t x, slong prec)
void arb_fmpz_poly_evaluate_arb(arb_t res, const fmpz_poly_t poly, const arb_t x, slong prec)

Set `res` to the value of the integer polynomial (given as a
polynomial or as `len` coefficients) at the real ball $x$, using
Horner's rule, rectangular splitting [@PS1973] (which multiplies the
integer coefficients by precomputed powers of $x$ in dot products), or
an automatic choice.

::: details
The default uses rectangular splitting for $\mathit{len} \ge
\max(6, 5 + 2500 / (\max(\mathit{prec}, 64) + 64))$, unless
$\mathit{prec} > 1024$ and the coefficients have more than
$\mathit{prec}/2$ bits, and Horner's rule otherwise.
:::
:::

::: function
void _arb_fmpz_poly_evaluate_acb_horner(acb_t res, const fmpz * poly, slong len, const acb_t x, slong prec)
void arb_fmpz_poly_evaluate_acb_horner(acb_t res, const fmpz_poly_t poly, const acb_t x, slong prec)
void _arb_fmpz_poly_evaluate_acb_rectangular(acb_t res, const fmpz * poly, slong len, const acb_t x, slong prec)
void arb_fmpz_poly_evaluate_acb_rectangular(acb_t res, const fmpz_poly_t poly, const acb_t x, slong prec)
void _arb_fmpz_poly_evaluate_acb(acb_t res, const fmpz * poly, slong len, const acb_t x, slong prec)
void arb_fmpz_poly_evaluate_acb(acb_t res, const fmpz_poly_t poly, const acb_t x, slong prec)

Set `res` to the value of the integer polynomial at the complex ball
$x$, using Horner's rule, rectangular splitting, or an automatic choice
(the real evaluation if $x$ is real, rectangular splitting otherwise).
:::

## Deflation {#deflation}

::: function
ulong arb_fmpz_poly_deflation(const fmpz_poly_t poly)

Returns the largest $d$ such that $\mathit{poly}(x) = g(x^d)$ for an
integer polynomial $g$, that is, the GCD of the exponents of the
nonconstant terms. For a constant polynomial the result is 1, and for
the zero polynomial 0.
:::

::: function
void arb_fmpz_poly_deflate(fmpz_poly_t res, const fmpz_poly_t poly, ulong deflation)

Sets `res` to $g$ where $\mathit{poly}(x) = g(x^d)$ with $d$ =
`deflation`, which must divide the exponent of every nonzero term of
`poly`.
:::

## Polynomial roots {#roots}

The functions in this section assume that the polynomial is
squarefree; this is not checked. For a general polynomial, isolate the
roots of the factors of a squarefree factorisation, which also gives
the multiplicities:

```c
fmpz_poly_factor_t fac;
fmpz_poly_factor_init(fac);
fmpz_poly_factor_squarefree(fac, poly);

for (i = 0; i < fac->num; i++)
{
    deg = fmpz_poly_degree(fac->p + i);
    flint_printf("%wd roots of multiplicity %wd\n", deg, fac->exp[i]);
    roots = _acb_vec_init(deg);
    arb_fmpz_poly_complex_roots(roots, fac->p + i, 0, prec);
    _acb_vec_clear(roots, deg);
}

fmpz_poly_factor_clear(fac);
```

::: macro
ARB_FMPZ_POLY_ROOTS_VERBOSE
ARB_FMPZ_POLY_ROOTS_RATIONAL

Flags for the root isolation functions (values 1 and 2): print
information about the progress of the computation, and (only for
`arb_fmpz_poly_real_roots`) compute the rational roots exactly first.
:::

::: function
void arb_fmpz_poly_refine_root_arb(arb_t res, const fmpz_poly_t poly, const arb_t initial, slong prec)

Given a ball `initial` that contains exactly one root of the
squarefree polynomial `poly`, in its interior unless `initial` is
exact, sets `res` to an enclosure of that root with a relative accuracy
of at least `prec` bits that is contained in `initial`. If the root is
an endpoint of `initial`, the function may fail to terminate.

::: details
The function maintains an exact interval $[\ell, h]$ with the root in
its interior and tries, at each step, an interval Newton step $N(z) = m -
f(m)/f'(z)$ (which proves uniqueness and gives quadratic convergence
when $N(z) \subseteq z$), then a floating-point Newton step whose
result is confirmed by sign changes at the ends of a small candidate
interval, and finally a bisection. The working precision follows the
accuracy of the current enclosure.
:::
:::

::: function
slong arb_fmpz_poly_real_roots(arb_ptr roots, const fmpz_poly_t poly, int flags, slong prec)
void arb_fmpz_poly_complex_roots(acb_ptr roots, const fmpz_poly_t poly, int flags, slong prec)

Write to `roots` disjoint enclosures of all real roots, or of all
complex roots, of the squarefree polynomial `poly`, each with a relative
accuracy of at least `prec` bits. `arb_fmpz_poly_real_roots` returns the
number of real roots and writes them in increasing order;
`arb_fmpz_poly_complex_roots` writes $\deg(\mathit{poly})$ roots: first
the real roots in increasing order (with imaginary parts exactly zero),
then the nonreal roots in conjugate pairs, the root in the upper
half-plane first. The output array must have room for the degree of
`poly` entries. A polynomial of degree less than 1 has no roots.

::: details
`arb_fmpz_poly_real_roots` deflates the polynomial ($f(x) = g(x^d)$, a
root at 0 removed), isolates the real roots of $g$ exactly with
`fmpz_poly_isolate_real_roots` (or of $g$ on $(0, \infty)$ for even $d$)
([](algorithms/polynomials/real-roots-integers)), divides out the exact
dyadic roots, refines the isolating intervals with
`arb_fmpz_poly_refine_root_arb` and takes $d$-th roots. With
`ARB_FMPZ_POLY_ROOTS_RATIONAL`, the rational roots are first found
exactly with `fmpz_poly_roots_fmpq` and divided out; this is much faster
for polynomials with many rational roots and adds some overhead
otherwise. Rational roots that are dyadic are returned as
exact balls.

`arb_fmpz_poly_complex_roots` deflates the polynomial, isolates the
real roots, and then runs `acb_poly_find_roots` with doubling precision,
starting at 53 bits, until all roots are isolated and accurate to `prec`
bits. If the real roots are numerous (a quarter of the roots or more)
or clustered, they are divided out in ball arithmetic first and only
the nonreal roots are computed numerically. Otherwise the reality of
the real roots is proved with `acb_poly_validate_real_roots`. See
[](arb_fmpz_poly/roots).
:::
:::

## Special polynomials {#special}

::: function
void arb_fmpz_poly_gauss_period_minpoly(fmpz_poly_t res, ulong q, ulong n)

Sets `res` to the minimal polynomial of the Gaussian periods
$\sum_{a \in H} \zeta^a$, where $\zeta = e^{2 \pi i / q}$ and $H$ runs
over the cosets of the subgroup of order $d = (q-1)/n$ of $(\Z /
q\Z)^{\times}$. The result has degree $n$; for $d = 1$ it is the
cyclotomic polynomial $\Phi_q$. Requires $q$ prime, $n \mid q - 1$ and
$\gcd(n, d) = 1$; otherwise `res` is set to zero.

::: details
The periods are computed numerically with `acb_dirichlet_roots`, the
product $\prod (x - \eta_k)$ is formed in ball arithmetic, and the
integer coefficients are recovered with `arb_poly_get_unique_fmpz_poly`,
doubling the precision until this succeeds. If additionally $(q-1)/f$
is coprime to $n$, where $f$ is the multiplicative order of a prime $p$
modulo $q$, the polynomial is irreducible over $\mathbb{F}_p$
[@CraPom2005], which gives a fast way to construct finite field extensions
of a given degree (see [](algorithms/finite-fields/defining-polynomials)).
:::
:::
