---
title: gr_special.h
header: gr_special.h
description: Mathematical constants, elementary functions, combinatorial numbers and special functions over generic rings.
category: Generic rings
order: 4
---

This module gives generic access to mathematical constants, elementary
and special functions and combinatorial numbers. Each function is a
method of the ring: it computes in the context `ctx`, with arguments and
result in that ring, and returns a status code ([`gr`](gr#status)).
The same call therefore evaluates $\Gamma(x)$ in real or complex ball
arithmetic, in floating-point arithmetic, as a power series, or returns
an exact value: `gr_fac_ui` gives $10!$ exactly in $\Z$, modulo $n$ in
$\Z/n\Z$ and as a ball in $\mathbb{R}$.

Which rings implement which functions:

* Real and complex balls (`gr_ctx_init_real_arb`,
  `gr_ctx_init_complex_acb`) implement nearly all functions, through the
  `arb`, `acb`, `arb_hypgeom`, `acb_hypgeom`, `acb_modular`,
  `acb_elliptic` and `acb_dirichlet` modules, with rigorous error bounds
  ([](algorithms/special-functions/overview)). Over the reals, a
  function whose value is not real returns `GR_DOMAIN` (for example
  $\log(-1)$) or `GR_UNABLE` when this cannot be decided.
* Floating-point contexts (`arf`, `acf`, `nfloat`, `dfloat`, `decimal`)
  implement the elementary functions and many special functions, with
  approximate results.
* Exact real and complex numbers (`ca`, `qqbar`) implement the
  elementary functions where the result is exact (`qqbar` for roots of
  unity and algebraic values such as `gr_sin_pi` of a rational number).
* Power series rings (`gr_series`) implement the functions of power
  series, via `gr_poly` (for example `_gr_poly_exp_series`).
* $\Z$, $\mathbb{Q}$, $\Z/n\Z$ and finite fields implement the combinatorial
  numbers, many through the generic implementations documented at the
  end, which work in any ring.

A function that is not implemented for a ring returns `GR_UNABLE`. The
defaults of the elementary functions handle only trivial arguments
($\exp(0) = 1$, $\log(1) = 0$, $\sin(0) = 0$, ...) and reduce some
functions to others ($\exp_2(x) = 2^x$ with `gr_pow`,
$\operatorname{acot}(x) = \operatorname{atan}(1/x)$).

Arguments that are parameters of a function precede the main argument,
in the same order as in the corresponding `arb` and `acb` functions:
`gr_bessel_j(res, nu, z, ctx)` computes $J_\nu(z)$. In the prototypes,
`x` and `y` are the first and second argument in this order.

## Mathematical constants {#constants}


::: function
int gr_pi(gr_ptr res, gr_ctx_t ctx)
int gr_euler(gr_ptr res, gr_ctx_t ctx)
int gr_catalan(gr_ptr res, gr_ctx_t ctx)
int gr_khinchin(gr_ptr res, gr_ctx_t ctx)
int gr_glaisher(gr_ptr res, gr_ctx_t ctx)

Set `res` to $\pi$, Euler's constant $\gamma$, Catalan's constant $G$, Khinchin's constant $K$ or Glaisher's constant $A$. Rings that do not contain the constant return `GR_DOMAIN` (`gr_pi` in $\mathbb{Q}$) or `GR_UNABLE`.
:::

## Exponentials and logarithms {#exp-log}


::: function
int gr_exp(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)
int gr_expm1(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)
int gr_exp2(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)
int gr_exp10(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)

Set `res` to $e^x$, $e^x - 1$ (accurate for small $x$), $2^x$ or $10^x$.
:::

::: function
int gr_exp_pi_i(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)

Sets `res` to $e^{\pi i x}$. For rational $x$ this is a root of unity, computed exactly by `qqbar`.
:::

::: function
int gr_log(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)
int gr_log1p(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)
int gr_log2(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)
int gr_log10(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)

Set `res` to the principal natural logarithm $\log x$, to $\log(1 + x)$ (accurate for small $x$), or to the logarithms in base 2 and 10. `GR_DOMAIN` for $x = 0$.
:::

::: function
int gr_log_pi_i(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)

Sets `res` to $\log(x) / (\pi i)$, the inverse of `gr_exp_pi_i` on the principal branch.
:::

::: function
int gr_lambertw(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)
int gr_lambertw_fmpz(gr_ptr res, gr_srcptr x, const fmpz_t k, gr_ctx_t ctx)

Set `res` to the principal branch $W_0(x)$ of the Lambert W function, or to the branch $W_k(x)$ with index $k$ given by the `fmpz_t` argument [@Joh2017b].
:::

## Trigonometric and hyperbolic functions {#trig}


::: function
int gr_sin(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)
int gr_cos(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)
int gr_sin_cos(gr_ptr res1, gr_ptr res2, gr_srcptr x, gr_ctx_t ctx)
int gr_tan(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)
int gr_cot(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)
int gr_sec(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)
int gr_csc(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)

Set `res` to $\sin x$, $\cos x$, both at once (`gr_sin_cos` sets `res1` to the sine and `res2` to the cosine), $\tan x$, $\cot x$, $\sec x = 1/\cos x$ or $\csc x = 1/\sin x$. Poles give `GR_DOMAIN`.
:::

::: function
int gr_sin_pi(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)
int gr_cos_pi(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)
int gr_sin_cos_pi(gr_ptr res1, gr_ptr res2, gr_srcptr x, gr_ctx_t ctx)
int gr_tan_pi(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)
int gr_cot_pi(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)
int gr_sec_pi(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)
int gr_csc_pi(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)

The same functions evaluated at $\pi x$: $\sin(\pi x)$ and so on. They are exact at integers and half-integers in ball arithmetic, and algebraic for rational $x$, which `qqbar` and `ca` compute exactly.
:::

::: function
int gr_sinc(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)
int gr_sinc_pi(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)

Set `res` to $\operatorname{sinc}(x) = \sin(x)/x$, respectively $\sin(\pi x)/(\pi x)$, with the value 1 at $x = 0$.
:::

::: function
int gr_sinh(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)
int gr_cosh(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)
int gr_sinh_cosh(gr_ptr res1, gr_ptr res2, gr_srcptr x, gr_ctx_t ctx)
int gr_tanh(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)
int gr_coth(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)
int gr_sech(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)
int gr_csch(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)

Set `res` to the hyperbolic functions $\sinh x$, $\cosh x$ (both at once for `gr_sinh_cosh`), $\tanh x$, $\coth x$, $\operatorname{sech} x$ and $\operatorname{csch} x$.
:::

## Inverse trigonometric and hyperbolic functions {#inverse-trig}

All functions use the principal branches of DLMF chapter 4.


::: function
int gr_asin(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)
int gr_acos(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)
int gr_atan(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)
int gr_acot(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)
int gr_asec(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)
int gr_acsc(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)

Set `res` to $\operatorname{asin} x$, $\operatorname{acos} x$, $\operatorname{atan} x$, $\operatorname{acot} x = \operatorname{atan}(1/x)$, $\operatorname{asec} x = \operatorname{acos}(1/x)$ and $\operatorname{acsc} x = \operatorname{asin}(1/x)$.
:::

::: function
int gr_atan2(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_ctx_t ctx)

Sets `res` to $\operatorname{atan2}(y, x)$, the argument of the point $(x, y)$, where the first argument is $y$ and the second $x$, as in C. Defined for real rings.
:::

::: function
int gr_asin_pi(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)
int gr_acos_pi(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)
int gr_atan_pi(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)
int gr_acot_pi(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)
int gr_asec_pi(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)
int gr_acsc_pi(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)

The inverse functions divided by $\pi$: $\operatorname{asin}(x)/\pi$ and so on, which are rational for many algebraic $x$.
:::

::: function
int gr_asinh(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)
int gr_acosh(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)
int gr_atanh(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)
int gr_acoth(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)
int gr_asech(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)
int gr_acsch(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)

Set `res` to the inverse hyperbolic functions.
:::

## Factorials and the gamma function {#gamma}


::: function
int gr_fac(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)
int gr_fac_ui(gr_ptr res, ulong x, gr_ctx_t ctx)
int gr_fac_fmpz(gr_ptr res, const fmpz_t x, gr_ctx_t ctx)
int gr_fac_vec(gr_ptr res, slong len, gr_ctx_t ctx)

Set `res` to the factorial $x! = \Gamma(x + 1)$. The *vec* version writes the first `len` factorials $0!, 1!, \ldots, (\mathit{len} - 1)!$ to the array `res`.
:::

::: function
int gr_rfac(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)
int gr_rfac_ui(gr_ptr res, ulong x, gr_ctx_t ctx)
int gr_rfac_fmpz(gr_ptr res, const fmpz_t x, gr_ctx_t ctx)
int gr_rfac_vec(gr_ptr res, slong len, gr_ctx_t ctx)

Set `res` to the reciprocal factorial $1/x! = 1/\Gamma(x+1)$, and the *vec* version to $1/0!, \ldots, 1/(\mathit{len} - 1)!$; `GR_DOMAIN` if a factorial is not invertible in the ring.
:::

::: function
int gr_bin(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_ctx_t ctx)
int gr_bin_ui(gr_ptr res, gr_srcptr x, ulong y, gr_ctx_t ctx)
int gr_bin_uiui(gr_ptr res, ulong x, ulong y, gr_ctx_t ctx)
int gr_bin_vec(gr_ptr res, gr_srcptr x, slong len, gr_ctx_t ctx)
int gr_bin_ui_vec(gr_ptr res, ulong x, slong len, gr_ctx_t ctx)

Set `res` to the binomial coefficient $\binom{x}{y}$. The *vec* versions write $\binom{x}{0}, \binom{x}{1}, \ldots, \binom{x}{\mathit{len}-1}$ to `res`. For a two-dimensional table (Pascal's triangle), `gr_mat_pascal` is faster than repeated calls.
:::

::: function
int gr_rising(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_ctx_t ctx)
int gr_rising_ui(gr_ptr res, gr_srcptr x, ulong y, gr_ctx_t ctx)
int gr_falling(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_ctx_t ctx)
int gr_falling_ui(gr_ptr res, gr_srcptr x, ulong y, gr_ctx_t ctx)

Set `res` to the rising factorial $x (x+1) \cdots (x+y-1)$ and the falling factorial $x (x-1) \cdots (x-y+1)$, or their extensions $\Gamma(x+y)/\Gamma(x)$ and $\Gamma(x+1)/\Gamma(x-y+1)$ to non-integer $y$.
:::

::: function
int gr_gamma(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)
int gr_gamma_fmpz(gr_ptr res, const fmpz_t x, gr_ctx_t ctx)
int gr_gamma_fmpq(gr_ptr res, const fmpq_t x, gr_ctx_t ctx)
int gr_rgamma(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)
int gr_lgamma(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)
int gr_digamma(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)

Set `res` to the gamma function $\Gamma(x)$, its reciprocal $1/\Gamma(x)$ (entire, zero at the poles of $\Gamma$), the log-gamma function $\log \Gamma(x)$ (the principal branch continuous away from the negative real axis) or the digamma function $\psi(x) = \Gamma'(x)/\Gamma(x)$. `gr_gamma` returns `GR_DOMAIN` at the poles $0, -1, -2, \ldots$. See [](algorithms/special-functions/gamma).
:::

::: function
int gr_beta(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_ctx_t ctx)

Sets `res` to the beta function $B(x, y) = \Gamma(x) \Gamma(y) / \Gamma(x + y)$.
:::

::: function
int gr_barnes_g(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)
int gr_log_barnes_g(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)

Set `res` to the Barnes G-function $G(x)$, respectively its logarithm $\log G(x)$.
:::

::: function
int gr_doublefac(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)
int gr_doublefac_ui(gr_ptr res, ulong x, gr_ctx_t ctx)

Set `res` to the double factorial $x!! = x (x - 2) (x - 4) \cdots$, ending at 1 or 2.
:::

::: function
int gr_harmonic(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)
int gr_harmonic_ui(gr_ptr res, ulong x, gr_ctx_t ctx)

Set `res` to the harmonic number $H_x = \sum_{k=1}^{x} 1/k$, or $\psi(x+1) + \gamma$ for non-integer $x$.
:::

## Combinatorial numbers {#combinatorial}

The *vec* functions for a sequence $c_n$ write the `len` values
$c_0, c_1, \ldots, c_{\mathit{len}-1}$ to the array `res`. Exact rings
receive exact values (or `GR_DOMAIN` if a denominator is not invertible),
rings with a real precision receive numerical values.


::: function
int gr_bernoulli_ui(gr_ptr res, ulong x, gr_ctx_t ctx)
int gr_bernoulli_fmpz(gr_ptr res, const fmpz_t x, gr_ctx_t ctx)
int gr_bernoulli_vec(gr_ptr res, slong len, gr_ctx_t ctx)

Set `res` to the Bernoulli number $B_n$ (with $B_1 = -1/2$). See [`bernoulli`](bernoulli) and [](algorithms/number-theory/bernoulli).
:::

::: function
int gr_eulernum_ui(gr_ptr res, ulong x, gr_ctx_t ctx)
int gr_eulernum_fmpz(gr_ptr res, const fmpz_t x, gr_ctx_t ctx)
int gr_eulernum_vec(gr_ptr res, slong len, gr_ctx_t ctx)

Set `res` to the Euler number $E_n$ (the coefficients of $\operatorname{sech}(t) = \sum E_n t^n / n!$).
:::

::: function
int gr_bernpoly_ui(gr_ptr res, ulong n, gr_srcptr x, gr_ctx_t ctx)
int gr_eulerpoly_ui(gr_ptr res, ulong n, gr_srcptr x, gr_ctx_t ctx)

Set `res` to the Bernoulli polynomial $B_n(x)$ or the Euler polynomial $E_n(x)$. These methods are implemented for real and complex balls (`arb_bernoulli_poly_ui`, `acb_bernoulli_poly_ui`); other rings return `GR_UNABLE`.

::: details
The ball versions compute $E_n(x) = 2 (B_{n+1}(x) - 2^{n+1} B_{n+1}(x/2)) / (n+1)$, with extra working precision to absorb the cancellation between the two terms for large $|x|$.
:::
:::

::: function
int gr_fib_ui(gr_ptr res, ulong x, gr_ctx_t ctx)
int gr_fib_fmpz(gr_ptr res, const fmpz_t x, gr_ctx_t ctx)
int gr_fib_vec(gr_ptr res, slong len, gr_ctx_t ctx)

Set `res` to the Fibonacci number $F_n$, with $F_0 = 0$ and $F_1 = 1$; `gr_fib_fmpz` accepts negative $n$, with $F_{-n} = (-1)^{n+1} F_n$.
:::

::: function
int gr_stirling_s1u_uiui(gr_ptr res, ulong x, ulong y, gr_ctx_t ctx)
int gr_stirling_s1_uiui(gr_ptr res, ulong x, ulong y, gr_ctx_t ctx)
int gr_stirling_s2_uiui(gr_ptr res, ulong x, ulong y, gr_ctx_t ctx)

Set `res` to the unsigned Stirling number of the first kind $\left[{x \atop y}\right]$, the signed Stirling number of the first kind $s(x, y)$, or the Stirling number of the second kind $\left\{{x \atop y}\right\}$ ([](algorithms/number-theory/stirling-bell)).
:::

::: function
int gr_stirling_s1u_ui_vec(gr_ptr res, ulong x, slong len, gr_ctx_t ctx)
int gr_stirling_s1_ui_vec(gr_ptr res, ulong x, slong len, gr_ctx_t ctx)
int gr_stirling_s2_ui_vec(gr_ptr res, ulong x, slong len, gr_ctx_t ctx)

Write the row $S(x, 0), S(x, 1), \ldots, S(x, \mathit{len} - 1)$ of the respective Stirling numbers to `res`. For a two-dimensional table, `gr_mat_stirling` is faster than repeated calls.
:::

::: function
int gr_bellnum_ui(gr_ptr res, ulong x, gr_ctx_t ctx)
int gr_bellnum_fmpz(gr_ptr res, const fmpz_t x, gr_ctx_t ctx)
int gr_bellnum_vec(gr_ptr res, slong len, gr_ctx_t ctx)

Set `res` to the Bell number $B_n$, the number of partitions of a set of $n$ elements ([](algorithms/number-theory/stirling-bell)).
:::

::: function
int gr_partitions_ui(gr_ptr res, ulong x, gr_ctx_t ctx)
int gr_partitions_fmpz(gr_ptr res, const fmpz_t x, gr_ctx_t ctx)
int gr_partitions_vec(gr_ptr res, slong len, gr_ctx_t ctx)

Set `res` to the number of partitions $p(n)$ of the integer $n$ ([](algorithms/number-theory/partitions)).
:::

## Error functions and exponential integrals {#erf}

For balls, the algorithms are described in
[](algorithms/special-functions/error-functions).

::: function
int gr_erf(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)
int gr_erfc(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)
int gr_erfcx(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)
int gr_erfi(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)

Set `res` to the error function $\operatorname{erf}(x)$, the complementary error function $\operatorname{erfc}(x) = 1 - \operatorname{erf}(x)$, the scaled complementary error function $\operatorname{erfcx}(x) = e^{x^2} \operatorname{erfc}(x)$, or the imaginary error function $\operatorname{erfi}(x) = -i \operatorname{erf}(ix)$.
:::

::: function
int gr_erfinv(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)
int gr_erfcinv(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)

Set `res` to the inverse error function and the inverse complementary error function, for real arguments in $(-1, 1)$, respectively $(0, 2)$.
:::

::: function
int gr_fresnel_s(gr_ptr res, gr_srcptr x, int normalized, gr_ctx_t ctx)
int gr_fresnel_c(gr_ptr res, gr_srcptr x, int normalized, gr_ctx_t ctx)
int gr_fresnel(gr_ptr res1, gr_ptr res2, gr_srcptr x, int normalized, gr_ctx_t ctx)

Set `res` to the Fresnel integral $S(x) = \int_0^x \sin(\pi t^2 / 2) \, dt$ or $C(x) = \int_0^x \cos(\pi t^2 / 2) \, dt$ if `normalized` is nonzero, and to the unnormalised integrals $\int_0^x \sin(t^2) \, dt$ and $\int_0^x \cos(t^2) \, dt$ otherwise; `gr_fresnel` sets `res1` to $S$ and `res2` to $C$.
:::

::: function
int gr_gamma_upper(gr_ptr res, gr_srcptr x, gr_srcptr y, int regularized, gr_ctx_t ctx)
int gr_gamma_lower(gr_ptr res, gr_srcptr x, gr_srcptr y, int regularized, gr_ctx_t ctx)

Set `res` to the incomplete gamma functions with parameter $s = x$ and argument $z = y$: if `regularized` is 0, the upper function $\Gamma(s, z) = \int_z^\infty t^{s-1} e^{-t} \, dt$ and the lower function $\gamma(s, z) = \int_0^z t^{s-1} e^{-t} \, dt$; if 1, the regularized versions $Q(s, z) = \Gamma(s,z)/\Gamma(s)$ and $P(s,z) = \gamma(s,z)/\Gamma(s)$; if 2, $z^{-s} \Gamma(s, z)$ and $\gamma^*(s, z) = z^{-s} P(s, z)$, as for `acb_hypgeom_gamma_upper` and `acb_hypgeom_gamma_lower`.
:::

::: function
int gr_beta_lower(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_srcptr z, int regularized, gr_ctx_t ctx)

Sets `res` to the incomplete beta function $B(a, b; z) = \int_0^z t^{a-1} (1-t)^{b-1} \, dt$ with $a = x$, $b = y$, or to the regularized version $I(a, b; z) = B(a, b; z)/B(a, b)$ if `regularized` is nonzero.
:::

::: function
int gr_exp_integral(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_ctx_t ctx)
int gr_exp_integral_ei(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)

Set `res` to the generalised exponential integral $E_s(z) = \int_1^\infty e^{-zt} t^{-s} \, dt$ with $s = x$ and $z = y$, respectively the exponential integral $\operatorname{Ei}(x)$.
:::

::: function
int gr_sin_integral(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)
int gr_cos_integral(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)
int gr_sinh_integral(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)
int gr_cosh_integral(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)

Set `res` to the sine integral $\operatorname{Si}(x)$, the cosine integral $\operatorname{Ci}(x)$, the hyperbolic sine integral $\operatorname{Shi}(x)$ or the hyperbolic cosine integral $\operatorname{Chi}(x)$.
:::

::: function
int gr_log_integral(gr_ptr res, gr_srcptr x, int offset, gr_ctx_t ctx)

Sets `res` to the logarithmic integral $\operatorname{li}(x) = \int_0^x dt / \log t$, or to the offset logarithmic integral $\operatorname{Li}(x) = \operatorname{li}(x) - \operatorname{li}(2)$ if `offset` is nonzero.
:::

::: function
int gr_dilog(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)

Sets `res` to the dilogarithm $\operatorname{Li}_2(x)$.
:::

## Bessel, Airy and Coulomb functions {#bessel}

For balls, the algorithms are described in
[](algorithms/special-functions/bessel-airy).

::: function
int gr_bessel_j(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_ctx_t ctx)
int gr_bessel_y(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_ctx_t ctx)
int gr_bessel_i(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_ctx_t ctx)
int gr_bessel_k(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_ctx_t ctx)
int gr_bessel_j_y(gr_ptr res1, gr_ptr res2, gr_srcptr x, gr_srcptr y, gr_ctx_t ctx)

Set `res` to the Bessel functions $J_\nu(z)$, $Y_\nu(z)$, $I_\nu(z)$ or $K_\nu(z)$ of order $\nu = x$ and argument $z = y$; `gr_bessel_j_y` sets `res1` to $J_\nu(z)$ and `res2` to $Y_\nu(z)$.
:::

::: function
int gr_bessel_i_scaled(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_ctx_t ctx)
int gr_bessel_k_scaled(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_ctx_t ctx)

Set `res` to the scaled Bessel functions $e^{-z} I_\nu(z)$ and $e^{z} K_\nu(z)$, which do not overflow for large $z$.
:::

::: function
int gr_bessel_j_jet(gr_ptr res, gr_srcptr nu, gr_srcptr z, slong len, gr_ctx_t ctx)

Sets `res` to the first `len` Taylor coefficients of $J_\nu(z + t)$ at $t = 0$, that is, $J_\nu^{(k)}(z) / k!$ for $0 \le k < \mathit{len}$ (the *jet* of $J_\nu$ at $z$). This is a generic function, not a method: it computes $J_\nu(z)$ and $J_{\nu+1}(z)$ and then the higher coefficients by a three-term recurrence (https://fungrim.org/entry/9b2f38/), and uses the series expansion at $z = 0$ for integer $\nu$.
:::

::: function
int gr_airy(gr_ptr res1, gr_ptr res2, gr_ptr res3, gr_ptr res4, gr_srcptr x, gr_ctx_t ctx)
int gr_airy_ai(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)
int gr_airy_bi(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)
int gr_airy_ai_prime(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)
int gr_airy_bi_prime(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)

Set `res` to the Airy functions $\operatorname{Ai}(x)$, $\operatorname{Bi}(x)$ and their derivatives; `gr_airy` sets `res1`, `res2`, `res3`, `res4` to $\operatorname{Ai}(x)$, $\operatorname{Ai}'(x)$, $\operatorname{Bi}(x)$, $\operatorname{Bi}'(x)$ (as `acb_hypgeom_airy`).
:::

::: function
int gr_airy_ai_zero(gr_ptr res, const fmpz_t n, gr_ctx_t ctx)
int gr_airy_bi_zero(gr_ptr res, const fmpz_t n, gr_ctx_t ctx)
int gr_airy_ai_prime_zero(gr_ptr res, const fmpz_t n, gr_ctx_t ctx)
int gr_airy_bi_prime_zero(gr_ptr res, const fmpz_t n, gr_ctx_t ctx)

Set `res` to the $n$-th real zero $a_n$, $b_n$, $a_n'$ or $b_n'$ of $\operatorname{Ai}$, $\operatorname{Bi}$, $\operatorname{Ai}'$ or $\operatorname{Bi}'$, counted from $n = 1$ and ordered by decreasing value; `GR_DOMAIN` for $n < 1$.
:::

::: function
int gr_coulomb(gr_ptr res1, gr_ptr res2, gr_ptr res3, gr_ptr res4, gr_srcptr x, gr_srcptr y, gr_srcptr z, gr_ctx_t ctx)
int gr_coulomb_f(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_srcptr z, gr_ctx_t ctx)
int gr_coulomb_g(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_srcptr z, gr_ctx_t ctx)
int gr_coulomb_hpos(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_srcptr z, gr_ctx_t ctx)
int gr_coulomb_hneg(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_srcptr z, gr_ctx_t ctx)

Set `res` to the Coulomb wave functions $F_\ell(\eta, z)$, $G_\ell(\eta, z)$, $H^{+}_\ell(\eta, z)$ and $H^{-}_\ell(\eta, z)$ with $\ell = x$, $\eta = y$ and argument $z$; `gr_coulomb` computes all four (as `acb_hypgeom_coulomb`).
:::

## Orthogonal polynomials {#orthogonal}


::: function
int gr_chebyshev_t(gr_ptr res, gr_srcptr n, gr_srcptr x, gr_ctx_t ctx)
int gr_chebyshev_t_fmpz(gr_ptr res, const fmpz_t n, gr_srcptr x, gr_ctx_t ctx)
int gr_chebyshev_u(gr_ptr res, gr_srcptr n, gr_srcptr x, gr_ctx_t ctx)
int gr_chebyshev_u_fmpz(gr_ptr res, const fmpz_t n, gr_srcptr x, gr_ctx_t ctx)

Set `res` to the Chebyshev polynomials $T_n(x)$ and $U_n(x)$, for an integer $n$ (the `_fmpz` versions, which work in any ring) or $n$ in the ring (as hypergeometric functions in ball arithmetic).
:::

::: function
int gr_jacobi_p(gr_ptr res, gr_srcptr n, gr_srcptr a, gr_srcptr b, gr_srcptr z, gr_ctx_t ctx)
int gr_gegenbauer_c(gr_ptr res, gr_srcptr n, gr_srcptr m, gr_srcptr z, gr_ctx_t ctx)
int gr_laguerre_l(gr_ptr res, gr_srcptr n, gr_srcptr m, gr_srcptr z, gr_ctx_t ctx)
int gr_hermite_h(gr_ptr res, gr_srcptr n, gr_srcptr z, gr_ctx_t ctx)

Set `res` to the Jacobi polynomial $P_n^{(a,b)}(z)$, the Gegenbauer polynomial $C_n^{(m)}(z)$, the associated Laguerre polynomial $L_n^{(m)}(z)$ or the Hermite polynomial $H_n(z)$, extended to non-integer $n$ as hypergeometric functions (as the corresponding `acb_hypgeom` functions).
:::

::: function
int gr_legendre_p(gr_ptr res, gr_srcptr n, gr_srcptr m, gr_srcptr z, int type, gr_ctx_t ctx)
int gr_legendre_q(gr_ptr res, gr_srcptr n, gr_srcptr m, gr_srcptr z, int type, gr_ctx_t ctx)

Set `res` to the associated Legendre functions of the first and second kind $P_n^m(z)$ and $Q_n^m(z)$; `type` selects the branch cut convention: 0 for Ferrers' functions on $(-1, 1)$ and 1 for the functions on $\mathbb{C} \setminus (-\infty, 1]$ (types 2 and 3 of *Mathematica*), as in `acb_hypgeom_legendre_p`.
:::

::: function
int gr_spherical_y_si(gr_ptr res, slong n, slong m, gr_srcptr theta, gr_srcptr phi, gr_ctx_t ctx)

Sets `res` to the spherical harmonic $Y_n^m(\theta, \phi)$ for integers $n$ and $m$.
:::

::: function
int gr_legendre_p_root_ui(gr_ptr root, gr_ptr weight, ulong n, ulong k, gr_ctx_t ctx)

Sets `root` to the root of index $k$ of the Legendre polynomial $P_n(x)$, with $0 \le k < n$ and the roots ordered by decreasing value, and, if `weight` is not `NULL`, `weight` to the corresponding weight of Gauss–Legendre quadrature (as `arb_hypgeom_legendre_p_ui_root`). `GR_DOMAIN` if $k \ge n$.
:::

## Hypergeometric functions {#hypergeometric}

The `flags` argument is passed to the `acb_hypgeom` functions: bit 0
(`flags = 1`) selects the regularized function, which divides by
$\Gamma(b)$ (by $\prod \Gamma(b_i)$ for ${}_pF_q$) and is entire in the
lower parameters. For ${}_2F_1$ the other bits are those of
`acb_hypgeom_2f1`. See [](algorithms/special-functions/hypergeometric)
and [](algorithms/special-functions/gauss-2f1).


::: function
int gr_hypgeom_0f1(gr_ptr res, gr_srcptr a, gr_srcptr z, int flags, gr_ctx_t ctx)
int gr_hypgeom_1f1(gr_ptr res, gr_srcptr a, gr_srcptr b, gr_srcptr z, int flags, gr_ctx_t ctx)
int gr_hypgeom_2f1(gr_ptr res, gr_srcptr a, gr_srcptr b, gr_srcptr c, gr_srcptr z, int flags, gr_ctx_t ctx)

Set `res` to the hypergeometric functions ${}_0F_1(a; z)$, ${}_1F_1(a; b; z)$ (Kummer's function $M(a, b, z)$) and ${}_2F_1(a, b; c; z)$, or to their regularized versions.
:::

::: function
int gr_hypgeom_u(gr_ptr res, gr_srcptr a, gr_srcptr b, gr_srcptr z, int flags, gr_ctx_t ctx)

Sets `res` to the confluent hypergeometric function of the second kind $U(a, b, z)$ (Tricomi's function). `flags` is ignored.
:::

::: function
int gr_hypgeom_pfq(gr_ptr res, const gr_vec_t a, const gr_vec_t b, gr_srcptr z, int flags, gr_ctx_t ctx)

Sets `res` to the generalised hypergeometric function ${}_pF_q(a_1, \ldots, a_p; b_1, \ldots, b_q; z)$ with the parameters given by the vectors `a` and `b`, or its regularized version.
:::

::: macro
GR_PFQ_OP(ctx, NAME)

The method of `gr_hypgeom_pfq`, cast to its signature
(`gr_method_pfq_op_op`).
:::

## Zeta and L-functions {#zeta}


::: function
int gr_zeta(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)
int gr_zeta_ui(gr_ptr res, ulong x, gr_ctx_t ctx)

Set `res` to the Riemann zeta function $\zeta(s)$; `GR_DOMAIN` at the pole $s = 1$. The method of `gr_zeta_ui` is implemented for real and complex balls (`arb_zeta_ui`); other rings return `GR_UNABLE`, and `gr_zeta` with an integer argument can be used instead.
:::

::: function
int gr_hurwitz_zeta(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_ctx_t ctx)

Sets `res` to the Hurwitz zeta function $\zeta(s, a) = \sum_{k \ge 0} (k + a)^{-s}$ with $s = x$, $a = y$ ([](algorithms/special-functions/zeta)).
:::

::: function
int gr_polylog(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_ctx_t ctx)

Sets `res` to the polylogarithm $\operatorname{Li}_s(z)$ with $s = x$, $z = y$.
:::

::: function
int gr_polygamma(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_ctx_t ctx)

Sets `res` to the polygamma function $\psi^{(s)}(z)$ with $s = x$, $z = y$, extended to complex $s$ as in `acb_polygamma`.
:::

::: function
int gr_lerch_phi(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_srcptr z, gr_ctx_t ctx)

Sets `res` to the Lerch transcendent $\Phi(z, s, a) = \sum_{k \ge 0} z^k (k + a)^{-s}$ with $z = x$, $s = y$, $a = z$.
:::

::: function
int gr_stieltjes(gr_ptr res, const fmpz_t x, gr_srcptr y, gr_ctx_t ctx)

Sets `res` to the generalised Stieltjes constant $\gamma_n(a)$, the coefficients of the Laurent expansion of $\zeta(s, a)$ at $s = 1$.
:::

::: function
int gr_dirichlet_eta(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)
int gr_dirichlet_beta(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)

Set `res` to the Dirichlet eta function $\eta(s) = (1 - 2^{1-s}) \zeta(s)$ and the Dirichlet beta function $\beta(s) = \sum_{k \ge 0} (-1)^k (2k+1)^{-s}$. `gr_dirichlet_beta` is implemented for real and complex balls, as the Dirichlet L-function of the nontrivial character modulo 4 (`acb_dirichlet_l`).
:::

::: function
int gr_riemann_xi(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)

Sets `res` to the Riemann xi function $\xi(s) = \tfrac{1}{2} s (s-1) \pi^{-s/2} \Gamma(s/2) \zeta(s)$.
:::

::: function
int gr_zeta_zero(gr_ptr res, const fmpz_t n, gr_ctx_t ctx)
int gr_zeta_zero_vec(gr_ptr res, const fmpz_t n, slong len, gr_ctx_t ctx)
int gr_zeta_nzeros(gr_ptr res, gr_srcptr t, gr_ctx_t ctx)

Set `res` to the $n$-th nontrivial zero $\rho_n = 1/2 + i t_n$ of $\zeta(s)$ in the upper half-plane, counted from $n = 1$ by increasing height; to the `len` consecutive zeros $\rho_n, \rho_{n+1}, \ldots$; or to the number $N(t)$ of zeros with imaginary part in $(0, t]$. The ball implementations return `GR_UNABLE` for $n < 1$ and for non-real or non-finite $t$.
:::

::: function
int gr_dirichlet_chi_fmpz(gr_ptr res, const dirichlet_group_t G, const dirichlet_char_t chi, const fmpz_t n, gr_ctx_t ctx)
int gr_dirichlet_chi_vec(gr_ptr res, const dirichlet_group_t G, const dirichlet_char_t chi, slong len, gr_ctx_t ctx)

Set `res` to the value $\chi(n)$ of the Dirichlet character `chi` of the group `G`, respectively the values $\chi(0), \ldots, \chi(\mathit{len}-1)$, as complex numbers. These are not methods: they are implemented for complex balls (`acb`) only and return `GR_UNABLE` otherwise. Declared only if `dirichlet.h` is included first.
:::

::: function
int gr_dirichlet_l(gr_ptr res, const dirichlet_group_t G, const dirichlet_char_t chi, gr_srcptr s, gr_ctx_t ctx)

Sets `res` to the Dirichlet L-function $L(s, \chi)$ for the character `chi` of the group `G`, in complex balls, in real balls (`GR_UNABLE` unless the value is real), or as a power series in `gr_series` rings (`gr_series_dirichlet_l`); `GR_UNABLE` in other rings.
:::

::: function
int gr_dirichlet_hardy_theta(gr_ptr res, const dirichlet_group_t G, const dirichlet_char_t chi, gr_srcptr t, gr_ctx_t ctx)
int gr_dirichlet_hardy_z(gr_ptr res, const dirichlet_group_t G, const dirichlet_char_t chi, gr_srcptr t, gr_ctx_t ctx)

Set `res` to the Hardy theta function $\theta(t)$ and Hardy Z-function $Z(t)$ of the character `chi`, defined so that $Z(t) = e^{i \theta(t)} L(1/2 + it, \chi)$ is real for real $t$, with the same rings as `gr_dirichlet_l`.
:::

## Modular forms and theta functions {#modular}

Here $\tau$ is in the upper half-plane and $q = e^{\pi i \tau}$, with the
conventions of [`acb_modular`](acb_modular).


::: function
int gr_jacobi_theta(gr_ptr res1, gr_ptr res2, gr_ptr res3, gr_ptr res4, gr_srcptr z, gr_srcptr tau, gr_ctx_t ctx)
int gr_jacobi_theta_1(gr_ptr res, gr_srcptr z, gr_srcptr tau, gr_ctx_t ctx)
int gr_jacobi_theta_2(gr_ptr res, gr_srcptr z, gr_srcptr tau, gr_ctx_t ctx)
int gr_jacobi_theta_3(gr_ptr res, gr_srcptr z, gr_srcptr tau, gr_ctx_t ctx)
int gr_jacobi_theta_4(gr_ptr res, gr_srcptr z, gr_srcptr tau, gr_ctx_t ctx)

Set `res` to the Jacobi theta functions $\theta_j(z, \tau)$, with the argument scaled by $\pi$ as in `acb_modular_theta` (for example $\theta_3(z, \tau) = \sum_{n \in \Z} q^{n^2} e^{2\pi i n z}$); `gr_jacobi_theta` computes all four.
:::

::: function
int gr_modular_j(gr_ptr res, gr_srcptr tau, gr_ctx_t ctx)
int gr_modular_lambda(gr_ptr res, gr_srcptr tau, gr_ctx_t ctx)
int gr_modular_delta(gr_ptr res, gr_srcptr tau, gr_ctx_t ctx)

Set `res` to the modular $j$-invariant $j(\tau)$, the modular lambda function $\lambda(\tau)$ or the modular discriminant $\Delta(\tau)$.
:::

::: function
int gr_dedekind_eta(gr_ptr res, gr_srcptr tau, gr_ctx_t ctx)
int gr_dedekind_eta_q(gr_ptr res, gr_srcptr q, gr_ctx_t ctx)

Set `res` to the Dedekind eta function $\eta(\tau)$. `gr_dedekind_eta_q` takes the nome $q = e^{2 \pi i \tau}$ instead and computes $\eta$ without the leading factor $q^{1/24}$, that is, $q^{-1/24} \eta(\tau) = \prod_{k \ge 1} (1 - q^k)$ (as `acb_modular_eta_sum` does), which is a single-valued function of $q$ for $|q| < 1$. It returns `GR_DOMAIN` if $|q| \ge 1$, and is implemented for real and complex balls.

::: details
For $|q| \le 1/2$ (or a ball $q$ containing 0) the ball versions sum the series with `acb_modular_eta_sum`; otherwise they compute $\tau = \log(q) / (2 \pi i)$ and use `acb_modular_eta`, which moves $\tau$ to the fundamental domain, and multiply by $q^{-1/24} = e^{-\log(q)/24}$. A ball $q$ that contains points of modulus both below and at least 1 gives `GR_UNABLE`.
:::
:::

::: function
int gr_hilbert_class_poly(gr_ptr res, slong D, gr_srcptr x, gr_ctx_t ctx)

Sets `res` to the value $H_D(x)$ of the Hilbert class polynomial for the discriminant $D$ at `x`. If `ctx` is a `gr_poly` ring and `x` its generator, the polynomial itself is returned. `GR_DOMAIN` if $D$ is not a negative discriminant. The default computes $H_D$ with `acb_modular_hilbert_class_poly` and evaluates it, so it works in any ring.
:::

::: function
int gr_eisenstein_e(gr_ptr res, ulong n, gr_srcptr tau, gr_ctx_t ctx)
int gr_eisenstein_g(gr_ptr res, ulong n, gr_srcptr tau, gr_ctx_t ctx)
int gr_eisenstein_g_vec(gr_ptr res, gr_srcptr tau, slong len, gr_ctx_t ctx)

Set `res` to the Eisenstein series $G_k(\tau) = \sum_{(m,n) \ne (0,0)} (m + n\tau)^{-k}$ or its normalisation $E_k(\tau) = G_k(\tau) / (2 \zeta(k))$, for even $k \ge 2$ (`GR_DOMAIN` for odd $k$ or $k = 0$); `gr_eisenstein_g_vec` writes $G_4, G_6, \ldots, G_{2\mathit{len}+2}$.
:::

## Elliptic integrals and functions {#elliptic}

The conventions are those of [`acb_elliptic`](acb_elliptic): the
parameter $m = k^2$, and the lattice $(1, \tau)$ for the elliptic
functions.


::: function
int gr_agm(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_ctx_t ctx)
int gr_agm1(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)

Set `res` to the arithmetic-geometric mean $M(x, y)$, respectively $M(1, x)$.
:::

::: function
int gr_elliptic_k(gr_ptr res, gr_srcptr m, gr_ctx_t ctx)
int gr_elliptic_e(gr_ptr res, gr_srcptr m, gr_ctx_t ctx)
int gr_elliptic_pi(gr_ptr res, gr_srcptr n, gr_srcptr m, gr_ctx_t ctx)

Set `res` to the complete elliptic integrals of the first, second and third kind $K(m)$, $E(m)$ and $\Pi(n, m)$.
:::

::: function
int gr_elliptic_f(gr_ptr res, gr_srcptr phi, gr_srcptr m, int pi, gr_ctx_t ctx)
int gr_elliptic_e_inc(gr_ptr res, gr_srcptr phi, gr_srcptr m, int pi, gr_ctx_t ctx)
int gr_elliptic_pi_inc(gr_ptr res, gr_srcptr n, gr_srcptr phi, gr_srcptr m, int pi, gr_ctx_t ctx)

Set `res` to the incomplete elliptic integrals $F(\phi, m)$, $E(\phi, m)$ and $\Pi(n, \phi, m)$; if `pi` is nonzero, the amplitude is $\pi \phi$ instead of $\phi$.
:::

::: function
int gr_carlson_rf(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_srcptr z, int flags, gr_ctx_t ctx)
int gr_carlson_rc(gr_ptr res, gr_srcptr x, gr_srcptr y, int flags, gr_ctx_t ctx)
int gr_carlson_rd(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_srcptr z, int flags, gr_ctx_t ctx)
int gr_carlson_rg(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_srcptr z, int flags, gr_ctx_t ctx)
int gr_carlson_rj(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_srcptr z, gr_srcptr w, int flags, gr_ctx_t ctx)

Set `res` to the Carlson symmetric elliptic integrals $R_F(x, y, z)$, $R_C(x, y) = R_F(x, y, y)$, $R_D(x, y, z) = R_J(x, y, z, z)$, $R_G(x, y, z)$ and $R_J(x, y, z, w)$. `flags` is passed to `acb_elliptic_rf` and `acb_elliptic_rj` and should be 0.
:::

::: function
int gr_elliptic_invariants(gr_ptr res1, gr_ptr res2, gr_srcptr tau, gr_ctx_t ctx)
int gr_elliptic_roots(gr_ptr res1, gr_ptr res2, gr_ptr res3, gr_srcptr tau, gr_ctx_t ctx)

Set `res1`, `res2` to the lattice invariants $g_2$, $g_3$, respectively `res1`, `res2`, `res3` to the roots $e_1$, $e_2$, $e_3$ of $4t^3 - g_2 t - g_3$, of the lattice $(1, \tau)$.
:::

::: function
int gr_weierstrass_p(gr_ptr res, gr_srcptr z, gr_srcptr tau, gr_ctx_t ctx)
int gr_weierstrass_p_prime(gr_ptr res, gr_srcptr z, gr_srcptr tau, gr_ctx_t ctx)
int gr_weierstrass_p_inv(gr_ptr res, gr_srcptr z, gr_srcptr tau, gr_ctx_t ctx)
int gr_weierstrass_zeta(gr_ptr res, gr_srcptr z, gr_srcptr tau, gr_ctx_t ctx)
int gr_weierstrass_sigma(gr_ptr res, gr_srcptr z, gr_srcptr tau, gr_ctx_t ctx)

Set `res` to the Weierstrass elliptic function $\wp(z, \tau)$, its derivative $\wp'(z, \tau)$, its inverse $\wp^{-1}(z, \tau)$ (an elliptic integral), the Weierstrass zeta function $\zeta(z, \tau)$ or the Weierstrass sigma function $\sigma(z, \tau)$ for the lattice $(1, \tau)$.
:::

## Generic implementations {#generic}

These are the defaults installed in every method table
([`gr_generic`](gr_generic)). They work in any ring, using the
arithmetic of the ring and, where useful, exact computation in $\Z$ or
$\mathbb{Q}$ followed by conversion, or computation with balls followed by
conversion in rings with a real precision.


::: function
int gr_generic_exp(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)
int gr_generic_log(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)
int gr_generic_sin(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)
int gr_generic_cos(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)
int gr_generic_sin_cos(gr_ptr res1, gr_ptr res2, gr_srcptr x, gr_ctx_t ctx)
int gr_generic_tan(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)
int gr_generic_tanh(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)
int gr_generic_asin(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)
int gr_generic_asinh(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)
int gr_generic_atan(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)
int gr_generic_atanh(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)

Handle the trivial argument ($0$, or $1$ for the logarithm) exactly and return `GR_UNABLE` otherwise.
:::

::: function
int gr_generic_expm1(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)
int gr_generic_exp2(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)
int gr_generic_exp10(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)
int gr_generic_log1p(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)
int gr_generic_log2(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)
int gr_generic_log10(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)

Reduce to other functions: $e^x - 1$, $2^x$ and $10^x$ with `gr_pow`, $\log(1+x)$, and $\log(x)/\log(2)$ or $\log(x)/\log(10)$ (any failure of the latter two is reported as `GR_UNABLE`).
:::

::: function
int gr_generic_sec_pi(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)
int gr_generic_csc_pi(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)

Compute $1/\cos(\pi x)$ and $1/\sin(\pi x)$.
:::

::: function
int gr_generic_asin_pi(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)
int gr_generic_acos_pi(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)
int gr_generic_atan_pi(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)
int gr_generic_acot_pi(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)
int gr_generic_asec_pi(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)
int gr_generic_acsc_pi(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)

Compute the inverse function and divide by $\pi$.
:::

::: function
int gr_generic_acot(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)
int gr_generic_asec(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)
int gr_generic_acsc(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)
int gr_generic_acoth(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)
int gr_generic_asech(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)
int gr_generic_acsch(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)

Compute $\operatorname{atan}(1/x)$, $\operatorname{acos}(1/x)$, $\operatorname{asin}(1/x)$ and the hyperbolic analogues; the inversion fails for $x = 0$ even where the function has a finite value.
:::

::: function
int gr_generic_fac(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)
int gr_generic_fac_ui(gr_ptr res, ulong n, gr_ctx_t ctx)
int gr_generic_fac_fmpz(gr_ptr res, const fmpz_t n, gr_ctx_t ctx)
int gr_generic_fac_vec(gr_ptr res, slong len, gr_ctx_t ctx)

Factorials: `fmpz_fac_ui` in $\Z$; a table for small $n$; in rings of finite characteristic, products of several consecutive integers per multiplication, and `n_factorial_fast_mod2_preinv` for `nmod`; numerically with `arb_fac_ui` in rings with a real precision; otherwise exactly in $\Z$ and converted. `gr_generic_fac` uses `gr_fac_ui` for integer arguments and $\Gamma(x + 1)$ otherwise.
:::

::: function
int gr_generic_rfac(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)
int gr_generic_rfac_ui(gr_ptr res, ulong n, gr_ctx_t ctx)
int gr_generic_rfac_fmpz(gr_ptr res, const fmpz_t n, gr_ctx_t ctx)
int gr_generic_rfac_vec(gr_ptr res, slong len, gr_ctx_t ctx)

Reciprocal factorials, computed as inverses of factorials.
:::

::: function
int gr_generic_rising(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_ctx_t ctx)
int gr_generic_rising_ui(gr_ptr res, gr_srcptr x, ulong y, gr_ctx_t ctx)
int gr_generic_falling(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_ctx_t ctx)
int gr_generic_falling_ui(gr_ptr res, gr_srcptr x, ulong y, gr_ctx_t ctx)

Rising factorials by binary splitting of the product $x (x+1) \cdots (x+n-1)$, with a basecase of 16 factors; falling factorials as rising factorials of $x - n + 1$. Non-integer $y$ gives `GR_UNABLE`.
:::

::: function
int gr_generic_bin(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_ctx_t ctx)
int gr_generic_bin_ui(gr_ptr res, gr_srcptr x, ulong y, gr_ctx_t ctx)
int gr_generic_bin_uiui(gr_ptr res, ulong x, ulong y, gr_ctx_t ctx)
int gr_generic_bin_vec(gr_ptr res, gr_srcptr x, slong len, gr_ctx_t ctx)
int gr_generic_bin_ui_vec(gr_ptr res, ulong x, slong len, gr_ctx_t ctx)

Binomial coefficients: `gr_generic_bin_uiui` uses `fmpz_bin_uiui` in $\Z$, and in other rings computes the integer exactly and converts it when $n \le 100$ or the ring is exact of characteristic zero, and otherwise divides a falling factorial by a factorial in the ring (falling back to the exact integer if the division fails). `gr_generic_bin_ui` divides the falling factorial $x (x-1) \cdots (x-k+1)$ by $k!$; `gr_generic_bin` requires an integer $y$. The *vec* versions use recurrences between consecutive coefficients, taking care in rings of finite characteristic where divisions may not be possible.
:::

::: function
int gr_generic_doublefac(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)
int gr_generic_doublefac_ui(gr_ptr res, ulong n, gr_ctx_t ctx)

Double factorials: $(2k)!! = 2^k k!$ and $(2k+1)!! = (2k+1)!/(2k)!!$; for non-integer $x$ the extension $2^{x/2} (\pi/2)^{(\cos(\pi x) - 1)/4} \Gamma(1 + x/2)$.
:::

::: function
int gr_generic_harmonic(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)
int gr_generic_harmonic_ui(gr_ptr res, ulong n, gr_ctx_t ctx)

Harmonic numbers: exactly with `fmpq_harmonic_ui` and converted, or numerically in rings with a real precision.
:::

::: function
int gr_generic_beta(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_ctx_t ctx)

Computes $B(x, y) = \Gamma(x) \Gamma(y) / \Gamma(x + y)$.
:::

::: function
int gr_generic_fib2_fmpz(gr_ptr v, gr_ptr u, const fmpz_t n, gr_ctx_t ctx)
int gr_generic_fib_fmpz(gr_ptr res, const fmpz_t n, gr_ctx_t ctx)
int gr_generic_fib_ui(gr_ptr res, ulong n, gr_ctx_t ctx)
int gr_generic_fib_vec(gr_ptr res, slong len, gr_ctx_t ctx)

Fibonacci numbers by a table for small $n$ and otherwise the doubling formulas $F_{2k-1} = F_k^2 + F_{k-1}^2$, $F_{2k+1} = 4F_k^2 - F_{k-1}^2 + 2(-1)^k$ along the binary expansion of $n$. `gr_generic_fib2_fmpz` sets `v` to $F_n$ and `u` to $F_{n-1}$. The *vec* version uses the recurrence.
:::

::: function
int gr_generic_bellnum_fmpz(gr_ptr res, const fmpz_t n, gr_ctx_t ctx)
int gr_generic_bellnum_ui(gr_ptr res, ulong n, gr_ctx_t ctx)
int gr_generic_bellnum_vec(gr_ptr res, slong len, gr_ctx_t ctx)

Bell numbers: exactly in $\Z$ with the functions of `arith`, numerically in rings with a real precision, directly modulo $n$ for `nmod` (`arith_bell_number_nmod` and `arith_bell_number_nmod_vec`), and otherwise exactly in $\Z$ and converted.
:::

::: function
int gr_generic_partitions_fmpz(gr_ptr res, const fmpz_t n, gr_ctx_t ctx)
int gr_generic_partitions_ui(gr_ptr res, ulong n, gr_ctx_t ctx)
int gr_generic_partitions_vec(gr_ptr res, slong len, gr_ctx_t ctx)

Partition numbers: `partitions_fmpz_ui` (the Hardy–Ramanujan–Rademacher formula [@Joh2012]) in $\Z$, numerically in rings with a real precision, and otherwise exactly and converted. The *vec* version uses `arith_number_of_partitions_vec` in $\Z$ and `arith_number_of_partitions_nmod_vec` for `nmod` (power series inversion of $\prod (1 - x^k)$), numerical values in rings with a real precision, and otherwise exact values converted.
:::

::: function
int gr_generic_chebyshev_t2_fmpz(gr_ptr a, gr_ptr b, const fmpz_t n, gr_srcptr x, gr_ctx_t ctx)
int gr_generic_chebyshev_t_fmpz(gr_ptr y, const fmpz_t n, gr_srcptr x, gr_ctx_t ctx)
int gr_generic_chebyshev_u2_fmpz(gr_ptr a, gr_ptr b, const fmpz_t n, gr_srcptr x, gr_ctx_t ctx)
int gr_generic_chebyshev_u_fmpz(gr_ptr y, const fmpz_t n, gr_srcptr x, gr_ctx_t ctx)

Chebyshev polynomials $T_n(x)$ and $U_n(x)$ for $n \ge 0$ by the doubling formulas $T_{2k} = 2T_k^2 - 1$ and $T_{2k+1} = 2 T_k T_{k+1} - x$ (and analogues for $U$) along the binary expansion of $n$, using $O(\log n)$ multiplications. The *2* versions set `a` to $T_n(x)$ and `b` to $T_{n-1}(x)$ (respectively $U$). Negative $n$ gives `GR_UNABLE`.
:::

::: function
int gr_generic_erfcx(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)

Computes $e^{x^2} \operatorname{erfc}(x)$ directly from `gr_exp` and `gr_erfc`, which is not numerically stable for large $x$; failures are reported as `GR_UNABLE`.
:::

::: function
int gr_generic_hilbert_class_poly(gr_ptr res, slong D, gr_srcptr x, gr_ctx_t ctx)

Computes $H_D$ with `acb_modular_hilbert_class_poly` and evaluates it at `x` (or returns it in a polynomial ring); `GR_DOMAIN` if $D$ is not a negative discriminant.
:::
