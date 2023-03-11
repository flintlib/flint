.. _arb_fpwrap:

**arb_fpwrap.h** -- floating-point wrappers of Arb mathematical functions
=========================================================================================

This module provides wrappers of Arb functions intended users who
want accurate floating-point mathematical functions
without necessarily caring about ball arithmetic.
The wrappers take floating-point input, give floating-point output,
and automatically increase the internal working precision
to ensure that the output is accurate
(in the rare case of failure, they output NaN along with an error code).

**Warning:** This module is experimental (as of Arb 2.21). It has not
been extensively tested, and interfaces may change in the future.

Supported types:

* ``double`` and ``complex_double`` (53-bit precision)

Limitations:

* The wrappers currently only handle finite input and points where function
  value is finite. For example,
  they do not know that `\log(0) = -\infty` or that `\exp(-\infty) = 0`.
  Singular input or output result in ``FPWRAP_UNABLE`` and a NaN output value.
  Evaluation of limit values may be implemented in the future for some functions.
* The wrappers currently treat ``-0.0`` as ``+0.0``. Users who need to
  distinguish signs of zero, e.g. on branch cuts, currently need to do so
  manually.
* When requesting *correct rounding*, the wrappers can fail to converge
  in asymptotic or exact cases (where special algorithms are required).
* If the value is computed accurately internally but is too small to represent
  as a floating-point number, the result will be ``-0.0`` or ``+0.0`` (on underflow)
  or ``-Inf`` or ``+Inf`` (on overflow). Since the underflowed or overflowed
  result is the best possible floating-point approximation of the true value,
  this outcome is considered correct and the flag ``FPWRAP_SUCCESS`` is returned.
  In the future, return status flags may be added to indicate that underflow
  or overflow has occurred.
* Different rounding modes are not yet implemented.

Option and return flags
-------------------------------------------------------------------------------

Functions return an ``int`` flag indicating the status.

.. macro:: FPWRAP_SUCCESS

    Indicates an accurate result. (Up to inevitable underflow or
    overflow in the final conversion to a floating-point result; see above.)

    This flag has the numerical value 0.

.. macro:: FPWRAP_UNABLE

    Indicates failure (unable to achieve to target accuracy,
    possibly because of a singularity). The output is set to NaN.

    This flag has the numerical value 1.

Functions take a *flags* parameter specifying optional rounding and termination
behavior. This can be set to 0 to use defaults.

.. macro:: FPWRAP_ACCURATE_PARTS

    For complex output, compute both real and imaginary parts to full relative accuracy.
    By default (if this flag is not set), complex results are computed to
    at least 53-bit accuracy as a whole, but if either the real or imaginary
    part is much smaller than the other, that part can have a large relative error.
    Setting this flag can result in slower evaluation or failure to converge
    in some cases.

    This flag has the numerical value 1.

.. macro:: FPWRAP_CORRECT_ROUNDING

    Guarantees *correct rounding*.
    By default (if this flag is not set), real results are accurate up to the
    rounding of the last bit, but the last bit is not guaranteed to
    be rounded optimally.
    Setting this flag can result in slower
    evaluation or failure to converge in some cases.
    Correct rounding automatically applies to both real and imaginary parts
    of complex numbers, so it is unnecessary to set both this flag and
    *FPWRAP_ACCURATE_PARTS*.

    This flag has the numerical value 2.

.. macro:: FPWRAP_WORK_LIMIT

    Multiplied by an integer, specifies the maximum working precision to use
    before giving up. With ``n * FPWRAP_WORK_LIMIT`` added to *flags*,
    `n` levels of precision will be used.
    The default `n = 0` is equivalent to `n = 8`, which for ``double``
    means trying with a working precision of 64, 128, 256, 512, 1024, 2048,
    4096, 8192 bits.
    With ``flags = 2 * FPWRAP_WORK_LIMIT``, we only try 64 and 128
    bits, and with ``flags = 16 * FPWRAP_WORK_LIMIT`` we
    go up to 2097152 bits.

    This flag has the numerical value 65536.

Types
-------------------------------------------------------------------------------

Outputs are passed by reference so that we can return status
flags and so that the interface is uniform for functions with
multiple outputs.

.. type:: complex_double

    A struct of two ``double`` components (``real`` and ``imag``), used to
    represent a machine-precision complex number. We use this custom type
    instead of the complex types defined in ``<complex.h>`` since Arb
    does not depend on C99. Users should easily be able to convert
    to the C99 complex type since the layout in memory is identical.

Functions
-------------------------------------------------------------------------------

Elementary functions
...............................................................................

.. function:: int arb_fpwrap_double_exp(double * res, double x, int flags)
              int arb_fpwrap_cdouble_exp(complex_double * res, complex_double x, int flags)

.. function:: int arb_fpwrap_double_expm1(double * res, double x, int flags)
              int arb_fpwrap_cdouble_expm1(complex_double * res, complex_double x, int flags)

.. function:: int arb_fpwrap_double_log(double * res, double x, int flags)
              int arb_fpwrap_cdouble_log(complex_double * res, complex_double x, int flags)

.. function:: int arb_fpwrap_double_log1p(double * res, double x, int flags)
              int arb_fpwrap_cdouble_log1p(complex_double * res, complex_double x, int flags)

.. function:: int arb_fpwrap_double_pow(double * res, double x, double y, int flags)
              int arb_fpwrap_cdouble_pow(complex_double * res, complex_double x, complex_double y, int flags)

.. function:: int arb_fpwrap_double_sqrt(double * res, double x, int flags)
              int arb_fpwrap_cdouble_sqrt(complex_double * res, complex_double x, int flags)

.. function:: int arb_fpwrap_double_rsqrt(double * res, double x, int flags)
              int arb_fpwrap_cdouble_rsqrt(complex_double * res, complex_double x, int flags)

.. function:: int arb_fpwrap_double_cbrt(double * res, double x, int flags)
              int arb_fpwrap_cdouble_cbrt(complex_double * res, complex_double x, int flags)

.. function:: int arb_fpwrap_double_sin(double * res, double x, int flags)
              int arb_fpwrap_cdouble_sin(complex_double * res, complex_double x, int flags)

.. function:: int arb_fpwrap_double_cos(double * res, double x, int flags)
              int arb_fpwrap_cdouble_cos(complex_double * res, complex_double x, int flags)

.. function:: int arb_fpwrap_double_tan(double * res, double x, int flags)
              int arb_fpwrap_cdouble_tan(complex_double * res, complex_double x, int flags)

.. function:: int arb_fpwrap_double_cot(double * res, double x, int flags)
              int arb_fpwrap_cdouble_cot(complex_double * res, complex_double x, int flags)

.. function:: int arb_fpwrap_double_sec(double * res, double x, int flags)
              int arb_fpwrap_cdouble_sec(complex_double * res, complex_double x, int flags)

.. function:: int arb_fpwrap_double_csc(double * res, double x, int flags)
              int arb_fpwrap_cdouble_csc(complex_double * res, complex_double x, int flags)

.. function:: int arb_fpwrap_double_sinc(double * res, double x, int flags)
              int arb_fpwrap_cdouble_sinc(complex_double * res, complex_double x, int flags)

.. function:: int arb_fpwrap_double_sin_pi(double * res, double x, int flags)
              int arb_fpwrap_cdouble_sin_pi(complex_double * res, complex_double x, int flags)

.. function:: int arb_fpwrap_double_cos_pi(double * res, double x, int flags)
              int arb_fpwrap_cdouble_cos_pi(complex_double * res, complex_double x, int flags)

.. function:: int arb_fpwrap_double_tan_pi(double * res, double x, int flags)
              int arb_fpwrap_cdouble_tan_pi(complex_double * res, complex_double x, int flags)

.. function:: int arb_fpwrap_double_cot_pi(double * res, double x, int flags)
              int arb_fpwrap_cdouble_cot_pi(complex_double * res, complex_double x, int flags)

.. function:: int arb_fpwrap_double_sinc_pi(double * res, double x, int flags)
              int arb_fpwrap_cdouble_sinc_pi(complex_double * res, complex_double x, int flags)

.. function:: int arb_fpwrap_double_asin(double * res, double x, int flags)
              int arb_fpwrap_cdouble_asin(complex_double * res, complex_double x, int flags)

.. function:: int arb_fpwrap_double_acos(double * res, double x, int flags)
              int arb_fpwrap_cdouble_acos(complex_double * res, complex_double x, int flags)

.. function:: int arb_fpwrap_double_atan(double * res, double x, int flags)
              int arb_fpwrap_cdouble_atan(complex_double * res, complex_double x, int flags)

.. function:: int arb_fpwrap_double_atan2(double * res, double x1, double x2, int flags)

.. function:: int arb_fpwrap_double_asinh(double * res, double x, int flags)
              int arb_fpwrap_cdouble_asinh(complex_double * res, complex_double x, int flags)

.. function:: int arb_fpwrap_double_acosh(double * res, double x, int flags)
              int arb_fpwrap_cdouble_acosh(complex_double * res, complex_double x, int flags)

.. function:: int arb_fpwrap_double_atanh(double * res, double x, int flags)
              int arb_fpwrap_cdouble_atanh(complex_double * res, complex_double x, int flags)

.. function:: int arb_fpwrap_double_lambertw(double * res, double x, slong branch, int flags)
              int arb_fpwrap_cdouble_lambertw(complex_double * res, complex_double x, slong branch, int flags)

Gamma, zeta and related functions
...............................................................................

.. function:: int arb_fpwrap_double_rising(double * res, double x, double n, int flags)
              int arb_fpwrap_cdouble_rising(complex_double * res, complex_double x, complex_double n, int flags)

    Rising factorial.

.. function:: int arb_fpwrap_double_gamma(double * res, double x, int flags)
              int arb_fpwrap_cdouble_gamma(complex_double * res, complex_double x, int flags)

    Gamma function.

.. function:: int arb_fpwrap_double_rgamma(double * res, double x, int flags)
              int arb_fpwrap_cdouble_rgamma(complex_double * res, complex_double x, int flags)

    Reciprocal gamma function.

.. function:: int arb_fpwrap_double_lgamma(double * res, double x, int flags)
              int arb_fpwrap_cdouble_lgamma(complex_double * res, complex_double x, int flags)

    Log-gamma function.

.. function:: int arb_fpwrap_double_digamma(double * res, double x, int flags)
              int arb_fpwrap_cdouble_digamma(complex_double * res, complex_double x, int flags)

    Digamma function.

.. function:: int arb_fpwrap_double_zeta(double * res, double x, int flags)
              int arb_fpwrap_cdouble_zeta(complex_double * res, complex_double x, int flags)

    Riemann zeta function.

.. function:: int arb_fpwrap_double_hurwitz_zeta(double * res, double s, double z, int flags)
              int arb_fpwrap_cdouble_hurwitz_zeta(complex_double * res, complex_double s, complex_double z, int flags)

    Hurwitz zeta function.

.. function:: int arb_fpwrap_double_lerch_phi(double * res, double z, double s, double a, int flags)
              int arb_fpwrap_cdouble_lerch_phi(complex_double * res, complex_double z, complex_double s, complex_double a, int flags)

    Lerch transcendent.

.. function:: int arb_fpwrap_double_barnes_g(double * res, double x, int flags)
              int arb_fpwrap_cdouble_barnes_g(complex_double * res, complex_double x, int flags)

    Barnes G-function.

.. function:: int arb_fpwrap_double_log_barnes_g(double * res, double x, int flags)
              int arb_fpwrap_cdouble_log_barnes_g(complex_double * res, complex_double x, int flags)

    Logarithmic Barnes G-function.

.. function:: int arb_fpwrap_double_polygamma(double * res, double s, double z, int flags)
              int arb_fpwrap_cdouble_polygamma(complex_double * res, complex_double s, complex_double z, int flags)

    Polygamma function.

.. function:: int arb_fpwrap_double_polylog(double * res, double s, double z, int flags)
              int arb_fpwrap_cdouble_polylog(complex_double * res, complex_double s, complex_double z, int flags)

    Polylogarithm.

.. function:: int arb_fpwrap_cdouble_dirichlet_eta(complex_double * res, complex_double s, int flags)

.. function:: int arb_fpwrap_cdouble_riemann_xi(complex_double * res, complex_double s, int flags)

.. function:: int arb_fpwrap_cdouble_hardy_theta(complex_double * res, complex_double z, int flags)

.. function:: int arb_fpwrap_cdouble_hardy_z(complex_double * res, complex_double z, int flags)

.. function:: int arb_fpwrap_cdouble_zeta_zero(complex_double * res, ulong n, int flags)

Error functions and exponential integrals
...............................................................................

.. function:: int arb_fpwrap_double_erf(double * res, double x, int flags)
              int arb_fpwrap_cdouble_erf(complex_double * res, complex_double x, int flags)

.. function:: int arb_fpwrap_double_erfc(double * res, double x, int flags)
              int arb_fpwrap_cdouble_erfc(complex_double * res, complex_double x, int flags)

.. function:: int arb_fpwrap_double_erfi(double * res, double x, int flags)
              int arb_fpwrap_cdouble_erfi(complex_double * res, complex_double x, int flags)

.. function:: int arb_fpwrap_double_erfinv(double * res, double x, int flags)
.. function:: int arb_fpwrap_double_erfcinv(double * res, double x, int flags)

.. function:: int arb_fpwrap_double_fresnel_s(double * res, double x, int normalized, int flags)
              int arb_fpwrap_cdouble_fresnel_s(complex_double * res, complex_double x, int normalized, int flags)

.. function:: int arb_fpwrap_double_fresnel_c(double * res, double x, int normalized, int flags)
              int arb_fpwrap_cdouble_fresnel_c(complex_double * res, complex_double x, int normalized, int flags)

.. function:: int arb_fpwrap_double_gamma_upper(double * res, double s, double z, int regularized, int flags)
              int arb_fpwrap_cdouble_gamma_upper(complex_double * res, complex_double s, complex_double z, int regularized, int flags)

.. function:: int arb_fpwrap_double_gamma_lower(double * res, double s, double z, int regularized, int flags)
              int arb_fpwrap_cdouble_gamma_lower(complex_double * res, complex_double s, complex_double z, int regularized, int flags)

.. function:: int arb_fpwrap_double_beta_lower(double * res, double a, double b, double z, int regularized, int flags)
              int arb_fpwrap_cdouble_beta_lower(complex_double * res, complex_double a, complex_double b, complex_double z, int regularized, int flags)

.. function:: int arb_fpwrap_double_exp_integral_e(double * res, double s, double z, int flags)
              int arb_fpwrap_cdouble_exp_integral_e(complex_double * res, complex_double s, complex_double z, int flags)

.. function:: int arb_fpwrap_double_exp_integral_ei(double * res, double x, int flags)
              int arb_fpwrap_cdouble_exp_integral_ei(complex_double * res, complex_double x, int flags)

.. function:: int arb_fpwrap_double_sin_integral(double * res, double x, int flags)
              int arb_fpwrap_cdouble_sin_integral(complex_double * res, complex_double x, int flags)

.. function:: int arb_fpwrap_double_cos_integral(double * res, double x, int flags)
              int arb_fpwrap_cdouble_cos_integral(complex_double * res, complex_double x, int flags)

.. function:: int arb_fpwrap_double_sinh_integral(double * res, double x, int flags)
              int arb_fpwrap_cdouble_sinh_integral(complex_double * res, complex_double x, int flags)

.. function:: int arb_fpwrap_double_cosh_integral(double * res, double x, int flags)
              int arb_fpwrap_cdouble_cosh_integral(complex_double * res, complex_double x, int flags)

.. function:: int arb_fpwrap_double_log_integral(double * res, double x, int offset, int flags)
              int arb_fpwrap_cdouble_log_integral(complex_double * res, complex_double x, int offset, int flags)

.. function:: int arb_fpwrap_double_dilog(double * res, double x, int flags)
              int arb_fpwrap_cdouble_dilog(complex_double * res, complex_double x, int flags)

Bessel, Airy and Coulomb functions
...............................................................................

.. function:: int arb_fpwrap_double_bessel_j(double * res, double nu, double x, int flags)
              int arb_fpwrap_cdouble_bessel_j(complex_double * res, complex_double nu, complex_double x, int flags)

.. function:: int arb_fpwrap_double_bessel_y(double * res, double nu, double x, int flags)
              int arb_fpwrap_cdouble_bessel_y(complex_double * res, complex_double nu, complex_double x, int flags)

.. function:: int arb_fpwrap_double_bessel_i(double * res, double nu, double x, int flags)
              int arb_fpwrap_cdouble_bessel_i(complex_double * res, complex_double nu, complex_double x, int flags)

.. function:: int arb_fpwrap_double_bessel_k(double * res, double nu, double x, int flags)
              int arb_fpwrap_cdouble_bessel_k(complex_double * res, complex_double nu, complex_double x, int flags)

.. function:: int arb_fpwrap_double_bessel_k_scaled(double * res, double nu, double x, int flags)
              int arb_fpwrap_cdouble_bessel_k_scaled(complex_double * res, complex_double nu, complex_double x, int flags)

.. function:: int arb_fpwrap_double_airy_ai(double * res, double x, int flags)
              int arb_fpwrap_cdouble_airy_ai(complex_double * res, complex_double x, int flags)

.. function:: int arb_fpwrap_double_airy_ai_prime(double * res, double x, int flags)
              int arb_fpwrap_cdouble_airy_ai_prime(complex_double * res, complex_double x, int flags)

.. function:: int arb_fpwrap_double_airy_bi(double * res, double x, int flags)
              int arb_fpwrap_cdouble_airy_bi(complex_double * res, complex_double x, int flags)

.. function:: int arb_fpwrap_double_airy_bi_prime(double * res, double x, int flags)
              int arb_fpwrap_cdouble_airy_bi_prime(complex_double * res, complex_double x, int flags)

.. function:: int arb_fpwrap_double_airy_ai_zero(double * res, ulong n, int flags)

.. function:: int arb_fpwrap_double_airy_ai_prime_zero(double * res, ulong n, int flags)

.. function:: int arb_fpwrap_double_airy_bi_zero(double * res, ulong n, int flags)

.. function:: int arb_fpwrap_double_airy_bi_prime_zero(double * res, ulong n, int flags)

.. function:: int arb_fpwrap_double_coulomb_f(double * res, double l, double eta, double x, int flags)
              int arb_fpwrap_cdouble_coulomb_f(complex_double * res, complex_double l, complex_double eta, complex_double x, int flags)

.. function:: int arb_fpwrap_double_coulomb_g(double * res, double l, double eta, double x, int flags)
              int arb_fpwrap_cdouble_coulomb_g(complex_double * res, complex_double l, complex_double eta, complex_double x, int flags)

.. function:: int arb_fpwrap_cdouble_coulomb_hpos(complex_double * res, complex_double l, complex_double eta, complex_double x, int flags)
              int arb_fpwrap_cdouble_coulomb_hneg(complex_double * res, complex_double l, complex_double eta, complex_double x, int flags)

Orthogonal polynomials
...............................................................................

.. function:: int arb_fpwrap_double_chebyshev_t(double * res, double n, double x, int flags)
              int arb_fpwrap_cdouble_chebyshev_t(complex_double * res, complex_double n, complex_double x, int flags)

.. function:: int arb_fpwrap_double_chebyshev_u(double * res, double n, double x, int flags)
              int arb_fpwrap_cdouble_chebyshev_u(complex_double * res, complex_double n, complex_double x, int flags)

.. function:: int arb_fpwrap_double_jacobi_p(double * res, double n, double a, double b, double x, int flags)
              int arb_fpwrap_cdouble_jacobi_p(complex_double * res, complex_double n, complex_double a, complex_double b, complex_double x, int flags)

.. function:: int arb_fpwrap_double_gegenbauer_c(double * res, double n, double m, double x, int flags)
              int arb_fpwrap_cdouble_gegenbauer_c(complex_double * res, complex_double n, complex_double m, complex_double x, int flags)

.. function:: int arb_fpwrap_double_laguerre_l(double * res, double n, double m, double x, int flags)
              int arb_fpwrap_cdouble_laguerre_l(complex_double * res, complex_double n, complex_double m, complex_double x, int flags)

.. function:: int arb_fpwrap_double_hermite_h(double * res, double n, double x, int flags)
              int arb_fpwrap_cdouble_hermite_h(complex_double * res, complex_double n, complex_double x, int flags)

.. function:: int arb_fpwrap_double_legendre_p(double * res, double n, double m, double x, int type, int flags)
              int arb_fpwrap_cdouble_legendre_p(complex_double * res, complex_double n, complex_double m, complex_double x, int type, int flags)

.. function:: int arb_fpwrap_double_legendre_q(double * res, double n, double m, double x, int type, int flags)
              int arb_fpwrap_cdouble_legendre_q(complex_double * res, complex_double n, complex_double m, complex_double x, int type, int flags)

.. function:: int arb_fpwrap_double_legendre_root(double * res1, double * res2, ulong n, ulong k, int flags)

    Sets *res1* to the index *k* root of the Legendre polynomial `P_n(x)`,
    and simultaneously sets *res2* to the corresponding weight for
    Gauss-Legendre quadrature.

.. function:: int arb_fpwrap_cdouble_spherical_y(complex_double * res, slong n, slong m, complex_double x1, complex_double x2, int flags)

Hypergeometric functions
...............................................................................

.. function:: int arb_fpwrap_double_hypgeom_0f1(double * res, double a, double x, int regularized, int flags)
              int arb_fpwrap_cdouble_hypgeom_0f1(complex_double * res, complex_double a, complex_double x, int regularized, int flags)

.. function:: int arb_fpwrap_double_hypgeom_1f1(double * res, double a, double b, double x, int regularized, int flags)
              int arb_fpwrap_cdouble_hypgeom_1f1(complex_double * res, complex_double a, complex_double b, complex_double x, int regularized, int flags)

.. function:: int arb_fpwrap_double_hypgeom_u(double * res, double a, double b, double x, int flags)
              int arb_fpwrap_cdouble_hypgeom_u(complex_double * res, complex_double a, complex_double b, complex_double x, int flags)

.. function:: int arb_fpwrap_double_hypgeom_2f1(double * res, double a, double b, double c, double x, int regularized, int flags)
              int arb_fpwrap_cdouble_hypgeom_2f1(complex_double * res, complex_double a, complex_double b, complex_double c, complex_double x, int regularized, int flags)

.. function:: int arb_fpwrap_double_hypgeom_pfq(double * res, const double * a, slong p, const double * b, slong q, double z, int regularized, int flags)
              int arb_fpwrap_cdouble_hypgeom_pfq(complex_double * res, const complex_double * a, slong p, const complex_double * b, slong q, complex_double z, int regularized, int flags)


Elliptic integrals, elliptic functions and modular forms
...............................................................................

.. function:: int arb_fpwrap_double_agm(double * res, double x, double y, int flags)
              int arb_fpwrap_cdouble_agm(complex_double * res, complex_double x, complex_double y, int flags)

    Arithmetic-geometric mean.

.. function:: int arb_fpwrap_cdouble_elliptic_k(complex_double * res, complex_double m, int flags)

.. function:: int arb_fpwrap_cdouble_elliptic_e(complex_double * res, complex_double m, int flags)

.. function:: int arb_fpwrap_cdouble_elliptic_pi(complex_double * res, complex_double n, complex_double m, int flags)

.. function:: int arb_fpwrap_cdouble_elliptic_f(complex_double * res, complex_double phi, complex_double m, int pi, int flags)

.. function:: int arb_fpwrap_cdouble_elliptic_e_inc(complex_double * res, complex_double phi, complex_double m, int pi, int flags)

.. function:: int arb_fpwrap_cdouble_elliptic_pi_inc(complex_double * res, complex_double n, complex_double phi, complex_double m, int pi, int flags)

    Complete and incomplete elliptic integrals.

.. function:: int arb_fpwrap_cdouble_elliptic_rf(complex_double * res, complex_double x, complex_double y, complex_double z, int option, int flags)

.. function:: int arb_fpwrap_cdouble_elliptic_rg(complex_double * res, complex_double x, complex_double y, complex_double z, int option, int flags)

.. function:: int arb_fpwrap_cdouble_elliptic_rj(complex_double * res, complex_double x, complex_double y, complex_double z, complex_double w, int option, int flags)

    Carlson symmetric elliptic integrals.

.. function:: int arb_fpwrap_cdouble_elliptic_p(complex_double * res, complex_double z, complex_double tau, int flags)

.. function:: int arb_fpwrap_cdouble_elliptic_p_prime(complex_double * res, complex_double z, complex_double tau, int flags)

.. function:: int arb_fpwrap_cdouble_elliptic_inv_p(complex_double * res, complex_double z, complex_double tau, int flags)

.. function:: int arb_fpwrap_cdouble_elliptic_zeta(complex_double * res, complex_double z, complex_double tau, int flags)

.. function:: int arb_fpwrap_cdouble_elliptic_sigma(complex_double * res, complex_double z, complex_double tau, int flags)

    Weierstrass elliptic functions.

.. function:: int arb_fpwrap_cdouble_jacobi_theta_1(complex_double * res, complex_double z, complex_double tau, int flags)

.. function:: int arb_fpwrap_cdouble_jacobi_theta_2(complex_double * res, complex_double z, complex_double tau, int flags)

.. function:: int arb_fpwrap_cdouble_jacobi_theta_3(complex_double * res, complex_double z, complex_double tau, int flags)

.. function:: int arb_fpwrap_cdouble_jacobi_theta_4(complex_double * res, complex_double z, complex_double tau, int flags)

    Jacobi theta functions.

.. function:: int arb_fpwrap_cdouble_dedekind_eta(complex_double * res, complex_double tau, int flags)

.. function:: int arb_fpwrap_cdouble_modular_j(complex_double * res, complex_double tau, int flags)

.. function:: int arb_fpwrap_cdouble_modular_lambda(complex_double * res, complex_double tau, int flags)

.. function:: int arb_fpwrap_cdouble_modular_delta(complex_double * res, complex_double tau, int flags)

Calling from C
-------------------------------------------------------------------------------

The program ``examples/fpwrap.c`` provides a usage example::

    #include "arb_fpwrap.h"

    int main()
    {
        double x, y;
        complex_double cx, cy;
        int flags = 0;    /* default options */

        x = 2.0;
        cx.real = 0.5;
        cx.imag = 123.0;

        arb_fpwrap_double_zeta(&y, x, flags);
        arb_fpwrap_cdouble_zeta(&cy, cx, flags);

        printf("zeta(%g) = %.16g\n", x, y);
        printf("zeta(%g + %gi) = %.16g + %.16gi\n", cx.real, cx.imag, cy.real, cy.imag);

        flint_cleanup();
        return 0;
    }

.. highlight:: text

This should print::

    > build/examples/fpwrap 
    zeta(2) = 1.644934066848226
    zeta(0.5 + 123i) = 0.006252861175594465 + 0.08206030514520983i

Note that this program does not check the return flag
to perform error handling.


Interfacing from Python
-------------------------------------------------------------------------------

.. highlight:: python

This illustrates how to call functions from Python using ``ctypes``::

    import ctypes
    import ctypes.util

    libarb_path = ctypes.util.find_library('arb')
    libarb = ctypes.CDLL(libarb_path)

    class _complex_double(ctypes.Structure):
        _fields_ = [('real', ctypes.c_double),
                    ('imag', ctypes.c_double)]

    def wrap_double_fun(fun):
        def f(x):
            y = ctypes.c_double()
            if fun(ctypes.byref(y), ctypes.c_double(x), 0):
                raise ValueError(f"unable to evaluate function accurately at {x}")
            return y.value
        return f

    def wrap_cdouble_fun(fun):
        def f(x):
            x = complex(x)
            cx = _complex_double()
            cy = _complex_double()
            cx.real = x.real
            cx.imag = x.imag
            if fun(ctypes.byref(cy), cx, 0):
                raise ValueError(f"unable to evaluate function accurately at {x}")
            return complex(cy.real, cy.imag)
        return f

    zeta = wrap_double_fun(libarb.arb_fpwrap_double_zeta)
    czeta = wrap_cdouble_fun(libarb.arb_fpwrap_cdouble_zeta)

    print(zeta(2.0))
    print(czeta(0.5+1e9j))
    print(zeta(1.0))       # pole, where wrapper throws exception

.. highlight:: text

This should print::

    1.6449340668482264
    (-2.761748029838061-1.6775122409894598j)
    Traceback (most recent call last):
      ...
    ValueError: unable to evaluate function accurately at 1.0

Interfacing from Julia
-------------------------------------------------------------------------------

.. highlight:: julia

This illustrates how to call functions from Julia using ``ccall``::

    using Libdl

    dlopen("/home/fredrik/src/arb/libarb.so")

    function zeta(x::Float64)
        cy = Ref{Float64}()
        if Bool(ccall((:arb_fpwrap_double_zeta, :libarb), Cint, (Ptr{Float64}, Float64, Cint), cy, x, 0))
            error("unable to evaluate accurately at ", x)
        end
        return cy[]
    end

    function zeta(x::Complex{Float64})
        cy = Ref{Complex{Float64}}()
        if Bool(ccall((:arb_fpwrap_cdouble_zeta, :libarb), Cint, (Ptr{Complex{Float64}}, Complex{Float64}, Cint), cy, x, 0))
            error("unable to evaluate accurately at ", x)
        end
        return cy[]
    end

    println(zeta(2.0))
    println(zeta(0.5 + 1e9im))
    println(zeta(1.0))       # pole, where wrapper throws exception

.. highlight:: text

This should print::

    1.6449340668482264
    -2.761748029838061 - 1.6775122409894598im
    ERROR: unable to evaluate accurately at 1.0
    Stacktrace:
     ...


.. highlight:: c

