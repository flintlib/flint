---
title: arb_fpwrap.h
header: arb_fpwrap.h
description: Floating-point wrappers of the real and complex special functions: double and complex double input and output with guaranteed accuracy, obtained by ball arithmetic with automatically increasing precision.
category: Real and complex numbers
order: 13
---

This module provides the special functions of FLINT's ball arithmetic
modules as ordinary floating-point functions of `double` and
`complex_double` arguments. Each wrapper converts its arguments to
balls, evaluates the function with the corresponding `arb` or `acb`
function at 64, 128, 256, ... bits until the result is accurate enough
to be rounded to a `double`, and returns that `double` with a status
code; if the target accuracy cannot be reached, it returns NaN and
`FPWRAP_UNABLE`. Since the error bounds of ball arithmetic are
rigorous, a successful result is accurate to within one unit in the
last place, or correctly rounded if requested. The guide
[](arb_fpwrap/usage) shows how to call the wrappers from C and from
other languages, and the tutorial section
[](tutorials/numerics/special-functions#fpwrap) uses them in a complete
program.

The module has the following limitations:

- Only finite input and points where the function value is finite are
  handled; for example, the wrappers do not know that $\log 0 =
  -\infty$ or that $\exp(-\infty) = 0$. Infinite or NaN input, poles
  and other singularities give `FPWRAP_UNABLE` and NaN.
- `-0.0` is treated as `+0.0`; signs of zero on branch cuts must be
  handled by the caller.
- With `FPWRAP_CORRECT_ROUNDING`, the wrappers can fail in exact or
  asymptotic cases, where deciding the rounding needs special methods.
- A value that is computed accurately but is too small or too large for
  a `double` is returned as $\pm 0$ or $\pm\infty$ with `FPWRAP_SUCCESS`,
  this being the best floating-point approximation.
- Only rounding to nearest is implemented.

The module was introduced in Arb 2.21 and its interface may still
change.

## Option and return flags {#flags}

The functions return an `int` status, and take a `flags` argument
that selects options; 0 gives the defaults.

::: macro
FPWRAP_SUCCESS
FPWRAP_UNABLE

Return values: `FPWRAP_SUCCESS` (0) indicates an accurate result (up to
unavoidable underflow or overflow in the final conversion);
`FPWRAP_UNABLE` (1) indicates that the target accuracy could not be
reached, possibly because of a singularity, and the output is NaN.
:::

::: macro
FPWRAP_ACCURATE_PARTS

Flag (value 1): for complex output, compute both the real part and the
imaginary part to full relative accuracy. By default, a complex result
is accurate to 53 bits relative to its absolute value, so that a part
much smaller than the other may have a large relative error. Setting
this flag can make the evaluation slower or make it fail.
:::

::: macro
FPWRAP_CORRECT_ROUNDING

Flag (value 2): guarantee correct rounding (to nearest) of the result,
or of both parts of a complex result. By default, the result is
accurate to within one unit in the last place but not necessarily
correctly rounded. Setting this flag can make the evaluation slower or
make it fail.
:::

::: macro
FPWRAP_WORK_LIMIT

The unit (value 65536) for the maximal working precision: with $n
\cdot$`FPWRAP_WORK_LIMIT` added to `flags`, the wrappers try the
precisions $64, 128, \ldots, 64 \cdot 2^{n}$ bits (at most $64 \cdot
2^{24}$). The default $n = 0$ tries $64, 128, \ldots, 8192$ bits.

::: details
The implementation stops after the precision $64 \cdot 2^n$, that
is, after $n + 1$ precision levels for $n \ge 1$, while the default uses
8 levels, the same as $n = 7$. With `2 * FPWRAP_WORK_LIMIT` the wrappers
thus try 64, 128 and 256 bits (the documentation of FLINT 3.6 and
earlier stated $n$ levels).
:::
:::

## Types {#types}

The outputs are passed by reference, so that the status can be
returned and functions with several outputs have a uniform interface.

::: type
complex_double

A structure with two `double` fields, `real` and `imag`, representing a
complex number. It is used instead of the C99 complex types; its memory
layout is the same as that of `double _Complex`, so conversion is
straightforward.
:::

## Accuracy {#accuracy}

A real result is accepted when the ball has a relative accuracy of at
least 54 bits (or, with `FPWRAP_CORRECT_ROUNDING`, when `arb_can_round_arf`
shows that its midpoint rounds to the same `double` as every point of the
ball), or when it is certainly smaller than $2^{-1077}$ or larger than
$2^{1024}$ in absolute value, in which case it rounds to zero or
infinity. A complex result is accepted under the same conditions
applied to the complex number as a whole, or, with
`FPWRAP_ACCURATE_PARTS` or `FPWRAP_CORRECT_ROUNDING`, to each part. The
inputs are exact (`double` values are converted to balls without
error), so the result is the correctly computed value of the function
at the floating-point input; ill-conditioning of the function near the
input only increases the precision needed.

## Elementary functions {#elementary}

::: function
int arb_fpwrap_double_exp(double * res, double x, int flags)
int arb_fpwrap_cdouble_exp(complex_double * res, complex_double x, int flags)

Compute $e^x$, using `arb_exp` and `acb_exp`.
:::

::: function
int arb_fpwrap_double_expm1(double * res, double x, int flags)
int arb_fpwrap_cdouble_expm1(complex_double * res, complex_double x, int flags)

Compute $e^x - 1$, accurately also for small $x$ (`arb_expm1`, `acb_expm1`).
:::

::: function
int arb_fpwrap_double_log(double * res, double x, int flags)
int arb_fpwrap_cdouble_log(complex_double * res, complex_double x, int flags)

Compute the natural logarithm (the principal branch for complex input), using `arb_log` and `acb_log`. The real version fails for $x \le 0$.
:::

::: function
int arb_fpwrap_double_log1p(double * res, double x, int flags)
int arb_fpwrap_cdouble_log1p(complex_double * res, complex_double x, int flags)

Compute $\log(1 + x)$, accurately also for small $x$ (`arb_log1p`, `acb_log1p`).
:::

::: function
int arb_fpwrap_double_pow(double * res, double x, double y, int flags)
int arb_fpwrap_cdouble_pow(complex_double * res, complex_double x, complex_double y, int flags)

Compute $x^y$ (`arb_pow`, `acb_pow`); the complex version uses the principal branch $e^{y \log x}$.
:::

::: function
int arb_fpwrap_double_sqrt(double * res, double x, int flags)
int arb_fpwrap_cdouble_sqrt(complex_double * res, complex_double x, int flags)

Compute the square root (principal branch), using `arb_sqrt` and `acb_sqrt`.
:::

::: function
int arb_fpwrap_double_rsqrt(double * res, double x, int flags)
int arb_fpwrap_cdouble_rsqrt(complex_double * res, complex_double x, int flags)

Compute the reciprocal square root $x^{-1/2}$ (`arb_rsqrt`, `acb_rsqrt`).
:::

::: function
int arb_fpwrap_double_cbrt(double * res, double x, int flags)
int arb_fpwrap_cdouble_cbrt(complex_double * res, complex_double x, int flags)

Compute the cube root: the real cube root $x^{1/3}$ (`arb_root_ui`; this fails for negative $x$) and the principal complex cube root.
:::

::: function
int arb_fpwrap_double_sin(double * res, double x, int flags)
int arb_fpwrap_cdouble_sin(complex_double * res, complex_double x, int flags)

Compute $\sin x$ (`arb_sin`, `acb_sin`).
:::

::: function
int arb_fpwrap_double_cos(double * res, double x, int flags)
int arb_fpwrap_cdouble_cos(complex_double * res, complex_double x, int flags)

Compute $\cos x$ (`arb_cos`, `acb_cos`).
:::

::: function
int arb_fpwrap_double_tan(double * res, double x, int flags)
int arb_fpwrap_cdouble_tan(complex_double * res, complex_double x, int flags)

Compute $\tan x$ (`arb_tan`, `acb_tan`).
:::

::: function
int arb_fpwrap_double_cot(double * res, double x, int flags)
int arb_fpwrap_cdouble_cot(complex_double * res, complex_double x, int flags)

Compute $\cot x$ (`arb_cot`, `acb_cot`).
:::

::: function
int arb_fpwrap_double_sec(double * res, double x, int flags)
int arb_fpwrap_cdouble_sec(complex_double * res, complex_double x, int flags)

Compute $\sec x = 1/\cos x$ (`arb_sec`, `acb_sec`).
:::

::: function
int arb_fpwrap_double_csc(double * res, double x, int flags)
int arb_fpwrap_cdouble_csc(complex_double * res, complex_double x, int flags)

Compute $\csc x = 1/\sin x$ (`arb_csc`, `acb_csc`).
:::

::: function
int arb_fpwrap_double_sinc(double * res, double x, int flags)
int arb_fpwrap_cdouble_sinc(complex_double * res, complex_double x, int flags)

Compute $\operatorname{sinc} x = \sin(x)/x$ (`arb_sinc`, `acb_sinc`).
:::

::: function
int arb_fpwrap_double_sin_pi(double * res, double x, int flags)
int arb_fpwrap_cdouble_sin_pi(complex_double * res, complex_double x, int flags)

Compute $\sin(\pi x)$ (`arb_sin_pi`, `acb_sin_pi`), which is exact at integers.
:::

::: function
int arb_fpwrap_double_cos_pi(double * res, double x, int flags)
int arb_fpwrap_cdouble_cos_pi(complex_double * res, complex_double x, int flags)

Compute $\cos(\pi x)$ (`arb_cos_pi`, `acb_cos_pi`).
:::

::: function
int arb_fpwrap_double_tan_pi(double * res, double x, int flags)
int arb_fpwrap_cdouble_tan_pi(complex_double * res, complex_double x, int flags)

Compute $\tan(\pi x)$ (`arb_tan_pi`, `acb_tan_pi`).
:::

::: function
int arb_fpwrap_double_cot_pi(double * res, double x, int flags)
int arb_fpwrap_cdouble_cot_pi(complex_double * res, complex_double x, int flags)

Compute $\cot(\pi x)$ (`arb_cot_pi`, `acb_cot_pi`).
:::

::: function
int arb_fpwrap_double_sinc_pi(double * res, double x, int flags)
int arb_fpwrap_cdouble_sinc_pi(complex_double * res, complex_double x, int flags)

Compute $\operatorname{sinc}(\pi x) = \sin(\pi x)/(\pi x)$ (`arb_sinc_pi`, `acb_sinc_pi`).
:::

::: function
int arb_fpwrap_double_asin(double * res, double x, int flags)
int arb_fpwrap_cdouble_asin(complex_double * res, complex_double x, int flags)

Compute $\operatorname{asin} x$ (`arb_asin`, `acb_asin`).
:::

::: function
int arb_fpwrap_double_acos(double * res, double x, int flags)
int arb_fpwrap_cdouble_acos(complex_double * res, complex_double x, int flags)

Compute $\operatorname{acos} x$ (`arb_acos`, `acb_acos`).
:::

::: function
int arb_fpwrap_double_atan(double * res, double x, int flags)
int arb_fpwrap_cdouble_atan(complex_double * res, complex_double x, int flags)

Compute $\operatorname{atan} x$ (`arb_atan`, `acb_atan`).
:::

::: function
int arb_fpwrap_double_atan2(double * res, double x1, double x2, int flags)

Computes the two-argument arctangent $\operatorname{atan2}(x_1, x_2)$, the argument of the point $(x_2, x_1)$, as `arb_atan2` with $y = x_1$ and $x = x_2$. It fails at the origin.
:::

::: function
int arb_fpwrap_double_asinh(double * res, double x, int flags)
int arb_fpwrap_cdouble_asinh(complex_double * res, complex_double x, int flags)

Compute $\operatorname{asinh} x$ (`arb_asinh`, `acb_asinh`).
:::

::: function
int arb_fpwrap_double_acosh(double * res, double x, int flags)
int arb_fpwrap_cdouble_acosh(complex_double * res, complex_double x, int flags)

Compute $\operatorname{acosh} x$ (`arb_acosh`, `acb_acosh`).
:::

::: function
int arb_fpwrap_double_atanh(double * res, double x, int flags)
int arb_fpwrap_cdouble_atanh(complex_double * res, complex_double x, int flags)

Compute $\operatorname{atanh} x$ (`arb_atanh`, `acb_atanh`).
:::

::: function
int arb_fpwrap_double_lambertw(double * res, double x, slong branch, int flags)
int arb_fpwrap_cdouble_lambertw(complex_double * res, complex_double x, slong branch, int flags)

Compute the Lambert W function on the given branch: for real input, `branch` must be 0 (the principal branch $W_0$) or $-1$ ($W_{-1}$); for complex input, any branch $k$ (`arb_lambertw`, `acb_lambertw`).
:::

## Gamma, zeta and related functions {#gamma}

::: function
int arb_fpwrap_double_rising(double * res, double x, double n, int flags)
int arb_fpwrap_cdouble_rising(complex_double * res, complex_double x, complex_double n, int flags)

Compute the rising factorial $(x)_n = \Gamma(x + n)/\Gamma(x)$ (`arb_rising`, `acb_rising`).
:::

::: function
int arb_fpwrap_double_gamma(double * res, double x, int flags)
int arb_fpwrap_cdouble_gamma(complex_double * res, complex_double x, int flags)

Compute $\Gamma(x)$ (`arb_gamma`, `acb_gamma`).
:::

::: function
int arb_fpwrap_double_rgamma(double * res, double x, int flags)
int arb_fpwrap_cdouble_rgamma(complex_double * res, complex_double x, int flags)

Compute $1/\Gamma(x)$ (`arb_rgamma`, `acb_rgamma`), an entire function.
:::

::: function
int arb_fpwrap_double_lgamma(double * res, double x, int flags)
int arb_fpwrap_cdouble_lgamma(complex_double * res, complex_double x, int flags)

Compute $\log \Gamma(x)$ (`arb_lgamma`, `acb_lgamma`); the real version requires $x > 0$, the complex version gives the principal branch of the log-gamma function.
:::

::: function
int arb_fpwrap_double_digamma(double * res, double x, int flags)
int arb_fpwrap_cdouble_digamma(complex_double * res, complex_double x, int flags)

Compute the digamma function $\psi(x)$ (`arb_digamma`, `acb_digamma`).
:::

::: function
int arb_fpwrap_double_zeta(double * res, double x, int flags)
int arb_fpwrap_cdouble_zeta(complex_double * res, complex_double x, int flags)

Compute the Riemann zeta function $\zeta(x)$ (`arb_zeta`, `acb_zeta`).
:::

::: function
int arb_fpwrap_double_hurwitz_zeta(double * res, double s, double z, int flags)
int arb_fpwrap_cdouble_hurwitz_zeta(complex_double * res, complex_double s, complex_double z, int flags)

Compute the Hurwitz zeta function $\zeta(s, z)$ (`arb_hurwitz_zeta`, `acb_hurwitz_zeta`).
:::

::: function
int arb_fpwrap_double_lerch_phi(double * res, double z, double s, double a, int flags)
int arb_fpwrap_cdouble_lerch_phi(complex_double * res, complex_double z, complex_double s, complex_double a, int flags)

Compute the Lerch transcendent $\Phi(z, s, a)$ (`acb_dirichlet_lerch_phi`); the real version fails if the value is not real.
:::

::: function
int arb_fpwrap_double_barnes_g(double * res, double x, int flags)
int arb_fpwrap_cdouble_barnes_g(complex_double * res, complex_double x, int flags)

Compute the Barnes G-function (`acb_barnes_g`; the real version fails if the result is not real).
:::

::: function
int arb_fpwrap_double_log_barnes_g(double * res, double x, int flags)
int arb_fpwrap_cdouble_log_barnes_g(complex_double * res, complex_double x, int flags)

Compute the logarithm of the Barnes G-function (`acb_log_barnes_g`).
:::

::: function
int arb_fpwrap_double_polygamma(double * res, double s, double z, int flags)
int arb_fpwrap_cdouble_polygamma(complex_double * res, complex_double s, complex_double z, int flags)

Compute the polygamma function $\psi^{(s)}(z)$ (`acb_polygamma`).
:::

::: function
int arb_fpwrap_double_polylog(double * res, double s, double z, int flags)
int arb_fpwrap_cdouble_polylog(complex_double * res, complex_double s, complex_double z, int flags)

Compute the polylogarithm $\operatorname{Li}_s(z)$ (`arb_polylog`, `acb_polylog`).
:::

::: function
int arb_fpwrap_cdouble_dirichlet_eta(complex_double * res, complex_double s, int flags)

Computes the Dirichlet eta function $\eta(s) = (1 - 2^{1-s}) \zeta(s)$ (`acb_dirichlet_eta`).
:::

::: function
int arb_fpwrap_cdouble_riemann_xi(complex_double * res, complex_double s, int flags)

Computes the Riemann xi function $\xi(s)$ (`acb_dirichlet_xi`).
:::

::: function
int arb_fpwrap_cdouble_hardy_theta(complex_double * res, complex_double z, int flags)

Computes the Hardy (Riemann–Siegel) theta function $\theta(z)$ (`acb_dirichlet_hardy_theta`).
:::

::: function
int arb_fpwrap_cdouble_hardy_z(complex_double * res, complex_double z, int flags)

Computes the Hardy Z-function $Z(z)$ (`acb_dirichlet_hardy_z`).
:::

::: function
int arb_fpwrap_cdouble_zeta_zero(complex_double * res, ulong n, int flags)

Computes the $n$-th nontrivial zero of the Riemann zeta function with positive imaginary part, $n \ge 1$ (`acb_dirichlet_zeta_zero`).
:::

## Error functions and exponential integrals {#error-functions}

::: function
int arb_fpwrap_double_erf(double * res, double x, int flags)
int arb_fpwrap_cdouble_erf(complex_double * res, complex_double x, int flags)

Compute the error function $\operatorname{erf}(x)$ (`arb_hypgeom_erf`, `acb_hypgeom_erf`).
:::

::: function
int arb_fpwrap_double_erfc(double * res, double x, int flags)
int arb_fpwrap_cdouble_erfc(complex_double * res, complex_double x, int flags)

Compute the complementary error function $\operatorname{erfc}(x) = 1 - \operatorname{erf}(x)$, accurately also where it is small (`arb_hypgeom_erfc`, `acb_hypgeom_erfc`).
:::

::: function
int arb_fpwrap_double_erfi(double * res, double x, int flags)
int arb_fpwrap_cdouble_erfi(complex_double * res, complex_double x, int flags)

Compute the imaginary error function $\operatorname{erfi}(x) = -i \operatorname{erf}(ix)$ (`arb_hypgeom_erfi`, `acb_hypgeom_erfi`).
:::

::: function
int arb_fpwrap_double_erfinv(double * res, double x, int flags)

Computes the inverse error function on $(-1, 1)$ (`arb_hypgeom_erfinv`).
:::

::: function
int arb_fpwrap_double_erfcinv(double * res, double x, int flags)

Computes the inverse complementary error function on $(0, 2)$ (`arb_hypgeom_erfcinv`).
:::

::: function
int arb_fpwrap_double_fresnel_s(double * res, double x, int normalized, int flags)
int arb_fpwrap_cdouble_fresnel_s(complex_double * res, complex_double x, int normalized, int flags)

Compute the Fresnel integral $S(x)$, normalised ($\int_0^x \sin(\pi t^2/2)\,dt$) if `normalized` is nonzero and $\int_0^x \sin(t^2)\,dt$ otherwise (`arb_hypgeom_fresnel`, `acb_hypgeom_fresnel`).
:::

::: function
int arb_fpwrap_double_fresnel_c(double * res, double x, int normalized, int flags)
int arb_fpwrap_cdouble_fresnel_c(complex_double * res, complex_double x, int normalized, int flags)

Compute the Fresnel integral $C(x)$, with the same normalisations as `arb_fpwrap_double_fresnel_s`.
:::

::: function
int arb_fpwrap_double_gamma_upper(double * res, double s, double z, int regularized, int flags)
int arb_fpwrap_cdouble_gamma_upper(complex_double * res, complex_double s, complex_double z, int regularized, int flags)

Compute the upper incomplete gamma function $\Gamma(s, z)$ (`arb_hypgeom_gamma_upper`, `acb_hypgeom_gamma_upper`); `regularized` selects the variants as in those functions (0: plain, 1: divided by $\Gamma(s)$, 2: the generalised exponential integral normalisation).
:::

::: function
int arb_fpwrap_double_gamma_lower(double * res, double s, double z, int regularized, int flags)
int arb_fpwrap_cdouble_gamma_lower(complex_double * res, complex_double s, complex_double z, int regularized, int flags)

Compute the lower incomplete gamma function $\gamma(s, z)$ (`arb_hypgeom_gamma_lower`, `acb_hypgeom_gamma_lower`), with `regularized` as for those functions.
:::

::: function
int arb_fpwrap_double_beta_lower(double * res, double a, double b, double z, int regularized, int flags)
int arb_fpwrap_cdouble_beta_lower(complex_double * res, complex_double a, complex_double b, complex_double z, int regularized, int flags)

Compute the incomplete beta function $B(a, b; z)$ (`arb_hypgeom_beta_lower`, `acb_hypgeom_beta_lower`), regularized if `regularized` is nonzero.
:::

::: function
int arb_fpwrap_double_exp_integral_e(double * res, double s, double z, int flags)
int arb_fpwrap_cdouble_exp_integral_e(complex_double * res, complex_double s, complex_double z, int flags)

Compute the generalised exponential integral $E_s(z)$ (`arb_hypgeom_expint`, `acb_hypgeom_expint`).
:::

::: function
int arb_fpwrap_double_exp_integral_ei(double * res, double x, int flags)
int arb_fpwrap_cdouble_exp_integral_ei(complex_double * res, complex_double x, int flags)

Compute the exponential integral $\operatorname{Ei}(x)$ (`arb_hypgeom_ei`, `acb_hypgeom_ei`).
:::

::: function
int arb_fpwrap_double_sin_integral(double * res, double x, int flags)
int arb_fpwrap_cdouble_sin_integral(complex_double * res, complex_double x, int flags)

Compute the sine integral $\operatorname{Si}(x)$ (`arb_hypgeom_si`, `acb_hypgeom_si`).
:::

::: function
int arb_fpwrap_double_cos_integral(double * res, double x, int flags)
int arb_fpwrap_cdouble_cos_integral(complex_double * res, complex_double x, int flags)

Compute the cosine integral $\operatorname{Ci}(x)$ (`arb_hypgeom_ci`, `acb_hypgeom_ci`).
:::

::: function
int arb_fpwrap_double_sinh_integral(double * res, double x, int flags)
int arb_fpwrap_cdouble_sinh_integral(complex_double * res, complex_double x, int flags)

Compute the hyperbolic sine integral $\operatorname{Shi}(x)$ (`arb_hypgeom_shi`, `acb_hypgeom_shi`).
:::

::: function
int arb_fpwrap_double_cosh_integral(double * res, double x, int flags)
int arb_fpwrap_cdouble_cosh_integral(complex_double * res, complex_double x, int flags)

Compute the hyperbolic cosine integral $\operatorname{Chi}(x)$ (`arb_hypgeom_chi`, `acb_hypgeom_chi`).
:::

::: function
int arb_fpwrap_double_log_integral(double * res, double x, int offset, int flags)
int arb_fpwrap_cdouble_log_integral(complex_double * res, complex_double x, int offset, int flags)

Compute the logarithmic integral $\operatorname{li}(x)$, or the offset integral $\operatorname{Li}(x) = \operatorname{li}(x) - \operatorname{li}(2)$ if `offset` is nonzero (`arb_hypgeom_li`, `acb_hypgeom_li`).
:::

::: function
int arb_fpwrap_double_dilog(double * res, double x, int flags)
int arb_fpwrap_cdouble_dilog(complex_double * res, complex_double x, int flags)

Compute the dilogarithm $\operatorname{Li}_2(x)$ (`arb_hypgeom_dilog`, `acb_hypgeom_dilog`).
:::

## Bessel, Airy and Coulomb functions {#bessel}

::: function
int arb_fpwrap_double_bessel_j(double * res, double nu, double x, int flags)
int arb_fpwrap_cdouble_bessel_j(complex_double * res, complex_double nu, complex_double x, int flags)

Compute the Bessel function of the first kind $J_\nu(x)$ (`arb_hypgeom_bessel_j`, `acb_hypgeom_bessel_j`).
:::

::: function
int arb_fpwrap_double_bessel_y(double * res, double nu, double x, int flags)
int arb_fpwrap_cdouble_bessel_y(complex_double * res, complex_double nu, complex_double x, int flags)

Compute the Bessel function of the second kind $Y_\nu(x)$ (`arb_hypgeom_bessel_y`, `acb_hypgeom_bessel_y`).
:::

::: function
int arb_fpwrap_double_bessel_i(double * res, double nu, double x, int flags)
int arb_fpwrap_cdouble_bessel_i(complex_double * res, complex_double nu, complex_double x, int flags)

Compute the modified Bessel function of the first kind $I_\nu(x)$ (`arb_hypgeom_bessel_i`, `acb_hypgeom_bessel_i`).
:::

::: function
int arb_fpwrap_double_bessel_k(double * res, double nu, double x, int flags)
int arb_fpwrap_cdouble_bessel_k(complex_double * res, complex_double nu, complex_double x, int flags)

Compute the modified Bessel function of the second kind $K_\nu(x)$ (`arb_hypgeom_bessel_k`, `acb_hypgeom_bessel_k`).
:::

::: function
int arb_fpwrap_double_bessel_k_scaled(double * res, double nu, double x, int flags)
int arb_fpwrap_cdouble_bessel_k_scaled(complex_double * res, complex_double nu, complex_double x, int flags)

Compute $e^x K_\nu(x)$, which avoids underflow for large $x$ (`arb_hypgeom_bessel_k_scaled`, `acb_hypgeom_bessel_k_scaled`).
:::

::: function
int arb_fpwrap_double_airy_ai(double * res, double x, int flags)
int arb_fpwrap_cdouble_airy_ai(complex_double * res, complex_double x, int flags)

Compute the Airy function $\operatorname{Ai}(x)$ (`arb_hypgeom_airy`, `acb_hypgeom_airy`).
:::

::: function
int arb_fpwrap_double_airy_ai_prime(double * res, double x, int flags)
int arb_fpwrap_cdouble_airy_ai_prime(complex_double * res, complex_double x, int flags)

Compute $\operatorname{Ai}'(x)$.
:::

::: function
int arb_fpwrap_double_airy_bi(double * res, double x, int flags)
int arb_fpwrap_cdouble_airy_bi(complex_double * res, complex_double x, int flags)

Compute the Airy function $\operatorname{Bi}(x)$.
:::

::: function
int arb_fpwrap_double_airy_bi_prime(double * res, double x, int flags)
int arb_fpwrap_cdouble_airy_bi_prime(complex_double * res, complex_double x, int flags)

Compute $\operatorname{Bi}'(x)$.
:::

::: function
int arb_fpwrap_double_airy_ai_zero(double * res, ulong n, int flags)

Computes the $n$-th zero $a_n$ of $\operatorname{Ai}$, $n \ge 1$ (`arb_hypgeom_airy_zero`); the zeros are negative.
:::

::: function
int arb_fpwrap_double_airy_ai_prime_zero(double * res, ulong n, int flags)

Computes the $n$-th zero $a'_n$ of $\operatorname{Ai}'$, $n \ge 1$.
:::

::: function
int arb_fpwrap_double_airy_bi_zero(double * res, ulong n, int flags)

Computes the $n$-th zero $b_n$ of $\operatorname{Bi}$, $n \ge 1$.
:::

::: function
int arb_fpwrap_double_airy_bi_prime_zero(double * res, ulong n, int flags)

Computes the $n$-th zero $b'_n$ of $\operatorname{Bi}'$, $n \ge 1$.
:::

::: function
int arb_fpwrap_double_coulomb_f(double * res, double l, double eta, double x, int flags)
int arb_fpwrap_cdouble_coulomb_f(complex_double * res, complex_double l, complex_double eta, complex_double x, int flags)

Compute the regular Coulomb wave function $F_\ell(\eta, x)$ (`arb_hypgeom_coulomb`, `acb_hypgeom_coulomb`).
:::

::: function
int arb_fpwrap_double_coulomb_g(double * res, double l, double eta, double x, int flags)
int arb_fpwrap_cdouble_coulomb_g(complex_double * res, complex_double l, complex_double eta, complex_double x, int flags)

Compute the irregular Coulomb wave function $G_\ell(\eta, x)$.
:::

::: function
int arb_fpwrap_cdouble_coulomb_hpos(complex_double * res, complex_double l, complex_double eta, complex_double x, int flags)
int arb_fpwrap_cdouble_coulomb_hneg(complex_double * res, complex_double l, complex_double eta, complex_double x, int flags)

Compute the complex Coulomb wave functions $H^{\pm}_\ell(\eta, z) = G_\ell \pm i F_\ell$ (`acb_hypgeom_coulomb`).
:::

## Orthogonal polynomials {#orthogonal}

::: function
int arb_fpwrap_double_chebyshev_t(double * res, double n, double x, int flags)
int arb_fpwrap_cdouble_chebyshev_t(complex_double * res, complex_double n, complex_double x, int flags)

Compute the Chebyshev function of the first kind $T_n(x)$ for real or complex $n$ (`arb_hypgeom_chebyshev_t`, `acb_hypgeom_chebyshev_t`).
:::

::: function
int arb_fpwrap_double_chebyshev_u(double * res, double n, double x, int flags)
int arb_fpwrap_cdouble_chebyshev_u(complex_double * res, complex_double n, complex_double x, int flags)

Compute the Chebyshev function of the second kind $U_n(x)$.
:::

::: function
int arb_fpwrap_double_jacobi_p(double * res, double n, double a, double b, double x, int flags)
int arb_fpwrap_cdouble_jacobi_p(complex_double * res, complex_double n, complex_double a, complex_double b, complex_double x, int flags)

Compute the Jacobi function $P_n^{(a,b)}(x)$ (`arb_hypgeom_jacobi_p`, `acb_hypgeom_jacobi_p`).
In FLINT 3.6 and earlier the real
version computed $P_n^{(a,x)}(0)$ instead; the complex version was
correct.
:::

::: function
int arb_fpwrap_double_gegenbauer_c(double * res, double n, double m, double x, int flags)
int arb_fpwrap_cdouble_gegenbauer_c(complex_double * res, complex_double n, complex_double m, complex_double x, int flags)

Compute the Gegenbauer function $C_n^{(m)}(x)$.
:::

::: function
int arb_fpwrap_double_laguerre_l(double * res, double n, double m, double x, int flags)
int arb_fpwrap_cdouble_laguerre_l(complex_double * res, complex_double n, complex_double m, complex_double x, int flags)

Compute the associated Laguerre function $L_n^{(m)}(x)$.
:::

::: function
int arb_fpwrap_double_hermite_h(double * res, double n, double x, int flags)
int arb_fpwrap_cdouble_hermite_h(complex_double * res, complex_double n, complex_double x, int flags)

Compute the Hermite function $H_n(x)$.
:::

::: function
int arb_fpwrap_double_legendre_p(double * res, double n, double m, double x, int type, int flags)
int arb_fpwrap_cdouble_legendre_p(complex_double * res, complex_double n, complex_double m, complex_double x, int type, int flags)

Compute the Legendre function of the first kind $P_n^m(x)$; `type` selects the branch cut convention as in `arb_hypgeom_legendre_p` (0 or 1).
:::

::: function
int arb_fpwrap_double_legendre_q(double * res, double n, double m, double x, int type, int flags)
int arb_fpwrap_cdouble_legendre_q(complex_double * res, complex_double n, complex_double m, complex_double x, int type, int flags)

Compute the Legendre function of the second kind $Q_n^m(x)$, with `type` as in `arb_hypgeom_legendre_q`.
:::

::: function
int arb_fpwrap_double_legendre_root(double * res1, double * res2, ulong n, ulong k, int flags)

Sets `res1` to the root of index $k$ (with $0 \le k < n$, in decreasing order) of the Legendre polynomial $P_n$ and `res2` to the corresponding Gauss–Legendre quadrature weight (`arb_hypgeom_legendre_p_ui_root`). Fails if $k \ge n$.
:::

::: function
int arb_fpwrap_cdouble_spherical_y(complex_double * res, slong n, slong m, complex_double x1, complex_double x2, int flags)

Computes the spherical harmonic $Y_n^m(\theta, \phi)$ with $\theta$ = `x1` and $\phi$ = `x2` (`acb_hypgeom_spherical_y`).
:::

## Hypergeometric functions {#hypergeometric}

::: function
int arb_fpwrap_double_hypgeom_0f1(double * res, double a, double x, int regularized, int flags)
int arb_fpwrap_cdouble_hypgeom_0f1(complex_double * res, complex_double a, complex_double x, int regularized, int flags)

Compute ${}_0F_1(;a;x)$, or its regularized version ${}_0F_1(;a;x)/\Gamma(a)$ if `regularized` is nonzero (`arb_hypgeom_0f1`, `acb_hypgeom_0f1`).
:::

::: function
int arb_fpwrap_double_hypgeom_1f1(double * res, double a, double b, double x, int regularized, int flags)
int arb_fpwrap_cdouble_hypgeom_1f1(complex_double * res, complex_double a, complex_double b, complex_double x, int regularized, int flags)

Compute the confluent hypergeometric function ${}_1F_1(a;b;x)$, regularized if `regularized` is nonzero (`arb_hypgeom_1f1`, `acb_hypgeom_1f1`).
:::

::: function
int arb_fpwrap_double_hypgeom_u(double * res, double a, double b, double x, int flags)
int arb_fpwrap_cdouble_hypgeom_u(complex_double * res, complex_double a, complex_double b, complex_double x, int flags)

Compute the confluent hypergeometric function $U(a, b, x)$ (`arb_hypgeom_u`, `acb_hypgeom_u`).
:::

::: function
int arb_fpwrap_double_hypgeom_2f1(double * res, double a, double b, double c, double x, int regularized, int flags)
int arb_fpwrap_cdouble_hypgeom_2f1(complex_double * res, complex_double a, complex_double b, complex_double c, complex_double x, int regularized, int flags)

Compute the Gauss hypergeometric function ${}_2F_1(a, b; c; x)$, regularized if `regularized` is nonzero (`arb_hypgeom_2f1`, `acb_hypgeom_2f1`).
:::

::: function
int arb_fpwrap_double_hypgeom_pfq(double * res, const double * a, slong p, const double * b, slong q, double z, int regularized, int flags)
int arb_fpwrap_cdouble_hypgeom_pfq(complex_double * res, const complex_double * a, slong p, const complex_double * b, slong q, complex_double z, int regularized, int flags)

Compute the generalised hypergeometric function ${}_pF_q(a_1, \ldots, a_p; b_1, \ldots, b_q; z)$ with the parameters given as arrays of length $p$ and $q$, regularized if `regularized` is nonzero (`arb_hypgeom_pfq`, `acb_hypgeom_pfq`).
:::

## Elliptic integrals, elliptic functions and modular forms {#elliptic}

::: function
int arb_fpwrap_double_agm(double * res, double x, double y, int flags)
int arb_fpwrap_cdouble_agm(complex_double * res, complex_double x, complex_double y, int flags)

Compute the arithmetic-geometric mean $M(x, y)$ (`arb_agm`, `acb_agm`).
:::

::: function
int arb_fpwrap_cdouble_elliptic_k(complex_double * res, complex_double m, int flags)

Computes the complete elliptic integral of the first kind $K(m)$ in terms of the parameter $m = k^2$ (`acb_elliptic_k`).
:::

::: function
int arb_fpwrap_cdouble_elliptic_e(complex_double * res, complex_double m, int flags)

Computes the complete elliptic integral of the second kind $E(m)$ (`acb_elliptic_e`).
:::

::: function
int arb_fpwrap_cdouble_elliptic_pi(complex_double * res, complex_double n, complex_double m, int flags)

Computes the complete elliptic integral of the third kind $\Pi(n, m)$ (`acb_elliptic_pi`).
:::

::: function
int arb_fpwrap_cdouble_elliptic_f(complex_double * res, complex_double phi, complex_double m, int pi, int flags)

Computes the incomplete elliptic integral of the first kind $F(\phi, m)$; if `pi` is nonzero, the argument is $\pi\phi$ instead of $\phi$ (`acb_elliptic_f`).
:::

::: function
int arb_fpwrap_cdouble_elliptic_e_inc(complex_double * res, complex_double phi, complex_double m, int pi, int flags)

Computes the incomplete elliptic integral of the second kind $E(\phi, m)$, with `pi` as for `arb_fpwrap_cdouble_elliptic_f` (`acb_elliptic_e_inc`).
:::

::: function
int arb_fpwrap_cdouble_elliptic_pi_inc(complex_double * res, complex_double n, complex_double phi, complex_double m, int pi, int flags)

Computes the incomplete elliptic integral of the third kind $\Pi(n, \phi, m)$, with `pi` as above (`acb_elliptic_pi_inc`).
:::

::: function
int arb_fpwrap_cdouble_elliptic_rf(complex_double * res, complex_double x, complex_double y, complex_double z, int option, int flags)

Computes the Carlson symmetric elliptic integral $R_F(x, y, z)$ (`acb_elliptic_rf`); `option` is passed as its `flags` and should be 0.
:::

::: function
int arb_fpwrap_cdouble_elliptic_rg(complex_double * res, complex_double x, complex_double y, complex_double z, int option, int flags)

Computes the Carlson symmetric elliptic integral $R_G(x, y, z)$ (`acb_elliptic_rg`), with `option` as above.
:::

::: function
int arb_fpwrap_cdouble_elliptic_rj(complex_double * res, complex_double x, complex_double y, complex_double z, complex_double w, int option, int flags)

Computes the Carlson symmetric elliptic integral $R_J(x, y, z, w)$ (`acb_elliptic_rj`), with `option` as above.
:::

::: function
int arb_fpwrap_cdouble_elliptic_p(complex_double * res, complex_double z, complex_double tau, int flags)

Computes the Weierstrass elliptic function $\wp(z, \tau)$ for the lattice with periods $1$ and $\tau$ (`acb_elliptic_p`).
:::

::: function
int arb_fpwrap_cdouble_elliptic_p_prime(complex_double * res, complex_double z, complex_double tau, int flags)

Computes $\wp'(z, \tau)$ (`acb_elliptic_p_prime`).
:::

::: function
int arb_fpwrap_cdouble_elliptic_inv_p(complex_double * res, complex_double z, complex_double tau, int flags)

Computes the inverse Weierstrass elliptic function $\wp^{-1}(z, \tau)$ (`acb_elliptic_inv_p`).
:::

::: function
int arb_fpwrap_cdouble_elliptic_zeta(complex_double * res, complex_double z, complex_double tau, int flags)

Computes the Weierstrass zeta function $\zeta(z, \tau)$ (`acb_elliptic_zeta`).
:::

::: function
int arb_fpwrap_cdouble_elliptic_sigma(complex_double * res, complex_double z, complex_double tau, int flags)

Computes the Weierstrass sigma function $\sigma(z, \tau)$ (`acb_elliptic_sigma`).
:::

::: function
int arb_fpwrap_cdouble_jacobi_theta_1(complex_double * res, complex_double z, complex_double tau, int flags)

Computes the Jacobi theta function $\theta_1(z, \tau)$ in the convention of `acb_modular_theta`.
:::

::: function
int arb_fpwrap_cdouble_jacobi_theta_2(complex_double * res, complex_double z, complex_double tau, int flags)

Computes $\theta_2(z, \tau)$.
:::

::: function
int arb_fpwrap_cdouble_jacobi_theta_3(complex_double * res, complex_double z, complex_double tau, int flags)

Computes $\theta_3(z, \tau)$.
:::

::: function
int arb_fpwrap_cdouble_jacobi_theta_4(complex_double * res, complex_double z, complex_double tau, int flags)

Computes $\theta_4(z, \tau)$.
:::

::: function
int arb_fpwrap_cdouble_dedekind_eta(complex_double * res, complex_double tau, int flags)

Computes the Dedekind eta function $\eta(\tau)$ (`acb_modular_eta`).
:::

::: function
int arb_fpwrap_cdouble_modular_j(complex_double * res, complex_double tau, int flags)

Computes the modular $j$-invariant $j(\tau)$ (`acb_modular_j`).
:::

::: function
int arb_fpwrap_cdouble_modular_lambda(complex_double * res, complex_double tau, int flags)

Computes the modular lambda function $\lambda(\tau)$ (`acb_modular_lambda`).
:::

::: function
int arb_fpwrap_cdouble_modular_delta(complex_double * res, complex_double tau, int flags)

Computes the modular discriminant $\Delta(\tau)$ (`acb_modular_delta`).
:::
