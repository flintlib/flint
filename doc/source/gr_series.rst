.. _gr-series:

**gr_series.h** -- formal power series over generic rings
===============================================================================

We provide two implementations of formal power series:

* ``gr_series_mod`` - truncated power series `R[[x]] / x^n`, with
  exact reresentation of elements.

* ``gr_series`` - infinite power series `R[[x]]`, with elements
  represented as finite truncations with an error term `O(x^n)`
  where *n* is tracked per-element. Polynomials of small degree can be
  represented exactly without error term.

It should be noted that these are different algebraic structures.
For example, `\mathbb{Q}[[x]]` is an integral domain, but 
`\mathbb{Q}[[x]] / x^n` is not.
Some operations are admissible for ``gr_series`` but not for
``gr_series_mod``, and vice versa.
Note that for :func:`gr_equal`, `0 = 0` is true in both `R[[x]]` and `R[[x]] / x^n`, but
`0 + O(x^n) = 0 + O(x^n)` (or even `0 + O(x^n) = 0`) is unknown in `R[[x]]`
since we do not know whether the big-O term hides some nonzero term.
For a computation where either representation may be used,
it is recommended to choose ``gr_series_mod`` as this will have less
overhead.

Context constructors
---------------------------------------------------------------------------------

.. function:: void gr_series_mod_ctx_init(gr_ctx_t ctx, gr_ctx_t base_ring, slong prec)

    Initializes *ctx* to a ring of truncated power series `R[[x]] / \langle x^n \rangle`
    over the given *base_ring*.
    Elements have type :type:`gr_poly_t`.
    It is assumed that all inputs are already truncated to length *n*,
    and this invariant is enforced for all outputs.

    The truncation order *n* of a ``gr_series_mod`` context is fixed.
    If `m \le n`, then elements belonging to a context object for `R[[x]] / x^m` can be given as
    input to a context for `R[[x]] / x^n`, but not vice versa.

.. function:: void gr_series_ctx_init(gr_ctx_t ctx, gr_ctx_t base_ring, slong prec)

    Initializes *ctx* to a ring of power series `R[[x]]` over the given *base_ring*.
    Elements have type :type:`gr_series_t`.

    Elements are generally inexact, having an error term `O(x^n)` where *n*
    may be different for each element.
    The context parameter *prec* defines the default precision: tails `c_n x^n + c_{n+1} x^{n+1} + \ldots`
    with `n \ge prec` will be replaced by the error term `O(x^n)`.

    The default precision is mutable. Elements can be shared freely between
    two context objects with the same base ring, even if those context
    objects have different precision.


Truncated power series
---------------------------------------------------------------------------------

Types and macros
.................................................................................

.. macro:: GR_SERIES_MOD_ELEM_CTX(ctx)

    Given a ``gr_series_mod`` context object, accesses the context object
    of the base ring.

.. macro:: GR_SERIES_MOD_N(ctx)

    Given a ``gr_series_mod`` context object, accesses the truncation
    order *n*.


Generic methods
.................................................................................

Elements of a ``gr_series_mod`` ring can be manipulated using the standard
``gr`` interface. The following is a list of implemented functions such
that ``gr_series_mod_foo`` overloads ``gr_foo`` in the method table
for ``gr_series_mod``; the user can optionally call ``gr_series_mod_foo``
directly.

.. function:: void gr_series_mod_ctx_clear(gr_ctx_t ctx)
              int gr_series_mod_ctx_write(gr_stream_t out, gr_ctx_t ctx)
              truth_t gr_series_mod_ctx_is_ring(gr_ctx_t ctx)
              truth_t gr_series_mod_ctx_is_commutative_ring(gr_ctx_t ctx)
              truth_t gr_series_mod_ctx_is_integral_domain(gr_ctx_t ctx)
              truth_t gr_series_mod_ctx_is_rational_vector_space(gr_ctx_t ctx)
              truth_t gr_series_mod_ctx_is_real_vector_space(gr_ctx_t ctx)
              truth_t gr_series_mod_ctx_is_complex_vector_space(gr_ctx_t ctx)
              truth_t gr_series_mod_ctx_is_field(gr_ctx_t ctx)
              int gr_series_mod_ctx_set_gen_name(gr_ctx_t ctx, const char * s)
              int gr_series_mod_ctx_set_gen_names(gr_ctx_t ctx, const char ** s)
              int gr_series_mod_gens_recursive(gr_vec_t vec, gr_ctx_t ctx)
              void gr_series_mod_init(gr_poly_t res, gr_ctx_t ctx)
              void gr_series_mod_clear(gr_poly_t res, gr_ctx_t ctx)
              void gr_series_mod_swap(gr_poly_t x, gr_poly_t y, gr_ctx_t ctx)
              int gr_series_mod_randtest(gr_poly_t res, flint_rand_t state, gr_ctx_t ctx)
              int gr_series_mod_write(gr_stream_t out, const gr_poly_t x, gr_ctx_t ctx)
              int gr_series_mod_zero(gr_poly_t res, gr_ctx_t ctx)
              int gr_series_mod_one(gr_poly_t res, gr_ctx_t ctx)
              int gr_series_mod_gen(gr_poly_t res, gr_ctx_t ctx)
              int gr_series_mod_set(gr_poly_t res, const gr_poly_t x, gr_ctx_t ctx)
              int gr_series_mod_set_other(gr_poly_t res, gr_srcptr x, gr_ctx_t x_ctx, gr_ctx_t ctx)
              void gr_series_mod_set_shallow(gr_poly_t res, const gr_poly_t x, gr_ctx_t ctx)
              int gr_series_mod_set_si(gr_poly_t res, slong c, gr_ctx_t ctx)
              int gr_series_mod_set_ui(gr_poly_t res, ulong c, gr_ctx_t ctx)
              int gr_series_mod_set_fmpz(gr_poly_t res, const fmpz_t c, gr_ctx_t ctx)
              int gr_series_mod_set_fmpq(gr_poly_t res, const fmpq_t c, gr_ctx_t ctx)
              truth_t gr_series_mod_is_zero(const gr_poly_t x, gr_ctx_t ctx)
              truth_t gr_series_mod_is_one(const gr_poly_t x, gr_ctx_t ctx)
              truth_t gr_series_mod_equal(const gr_poly_t x, const gr_poly_t y, gr_ctx_t ctx)
              int gr_series_mod_neg(gr_poly_t res, const gr_poly_t x, gr_ctx_t ctx)
              int gr_series_mod_add(gr_poly_t res, const gr_poly_t x, const gr_poly_t y, gr_ctx_t ctx)
              int gr_series_mod_sub(gr_poly_t res, const gr_poly_t x, const gr_poly_t y, gr_ctx_t ctx)
              int gr_series_mod_mul(gr_poly_t res, const gr_poly_t x, const gr_poly_t y, gr_ctx_t ctx)
              int gr_series_mod_inv(gr_poly_t res, const gr_poly_t x, gr_ctx_t ctx)
              int gr_series_mod_div(gr_poly_t res, const gr_poly_t x, const gr_poly_t y, gr_ctx_t ctx)
              int gr_series_mod_exp(gr_poly_t res, const gr_poly_t x, gr_ctx_t ctx)
              int gr_series_mod_log(gr_poly_t res, const gr_poly_t x, gr_ctx_t ctx)
              int gr_series_mod_sqrt(gr_poly_t res, const gr_poly_t x, gr_ctx_t ctx)
              int gr_series_mod_rsqrt(gr_poly_t res, const gr_poly_t x, gr_ctx_t ctx)
              int gr_series_mod_tan(gr_poly_t res, const gr_poly_t x, gr_ctx_t ctx)
              int gr_series_mod_asin(gr_poly_t res, const gr_poly_t x, gr_ctx_t ctx)
              int gr_series_mod_acos(gr_poly_t res, const gr_poly_t x, gr_ctx_t ctx)
              int gr_series_mod_atan(gr_poly_t res, const gr_poly_t x, gr_ctx_t ctx)
              int gr_series_mod_asinh(gr_poly_t res, const gr_poly_t x, gr_ctx_t ctx)
              int gr_series_mod_acosh(gr_poly_t res, const gr_poly_t x, gr_ctx_t ctx)
              int gr_series_mod_atanh(gr_poly_t res, const gr_poly_t x, gr_ctx_t ctx)

Power series
---------------------------------------------------------------------------------

Types and macros
.................................................................................

.. type:: gr_series_struct
          gr_series_t

    A structure containing a :type:`gr_poly_t` followed by a :type:`slong`
    representing the exponent in the error term.

.. type:: gr_series_vec_struct
          gr_series_vec_t

    Alias for :type:`gr_vec_t` with :type:`gr_series_t` elements,
    provided for convenience.

.. macro:: GR_SERIES_POLY(x)

    Macro accessing the polynomial part of a :type:`gr_series_t`.

.. macro:: GR_SERIES_ERROR(x)

    Macro accessing the error of a :type:`gr_series_t`.

.. macro:: GR_SERIES_ERR_MAX

    The maximum allowed *n* in an error term `O(x^n)`.

.. macro:: GR_SERIES_ERR_EXACT

    A special value of *n* used to indicate that a series is exact.

.. macro:: GR_SERIES_ELEM_CTX(ctx)

    Given a ``gr_series`` context object, accesses the context object
    of the base ring.

.. macro:: GR_SERIES_PREC(ctx)

    Given a ``gr_series`` context object, accesses the default precision.

Error term manipulation
.................................................................................

.. function:: slong _gr_series_get_error(const gr_series_t f, gr_ctx_t ctx)

    Return the exponent `n` of the error term `O(x^n)` of *f*.
    If *f* is exact, returns ``GR_SERIES_ERR_EXACT``.

.. function:: truth_t _gr_series_is_exact(const gr_series_t f, gr_ctx_t ctx)

    Returns whether *f* is exact as a power series, i.e. lacks error term
    `O(x^n)`. This does not recursively check exactness of the underlying elements.

.. function:: void _gr_series_set_error(gr_series_t f, slong n, gr_ctx_t ctx)

    Add an error term `O(x^n)` in-place to *f*. The exponent *n* is
    clamped between 0 and ``GR_SERIES_ERR_MAX``. Terms of order higher
    than or equal to *n* will be truncated.

.. function:: void _gr_series_make_exact(gr_series_t f, gr_ctx_t ctx)

    Remove the `O(x^n)` error term (if any) from *f*.

Generic methods
.................................................................................

Elements of a ``gr_series`` ring can be manipulated using the standard
``gr`` interface. The following is a list of implemented functions such
that ``gr_series_foo`` overloads ``gr_foo`` in the method table
for ``gr_series``; the user can optionally call ``gr_series_foo``
directly.

.. function:: void gr_series_ctx_clear(gr_ctx_t ctx)
              int gr_series_ctx_write(gr_stream_t out, gr_ctx_t ctx)
              truth_t gr_series_ctx_is_ring(gr_ctx_t ctx)
              truth_t gr_series_ctx_is_commutative_ring(gr_ctx_t ctx)
              truth_t gr_series_ctx_is_integral_domain(gr_ctx_t ctx)
              truth_t gr_series_ctx_is_rational_vector_space(gr_ctx_t ctx)
              truth_t gr_series_ctx_is_real_vector_space(gr_ctx_t ctx)
              truth_t gr_series_ctx_is_complex_vector_space(gr_ctx_t ctx)
              int gr_series_ctx_set_gen_name(gr_ctx_t ctx, const char * s)
              int gr_series_ctx_set_gen_names(gr_ctx_t ctx, const char ** s)
              void gr_series_init(gr_series_t res, gr_ctx_t ctx)
              void gr_series_clear(gr_series_t res, gr_ctx_t ctx)
              void gr_series_swap(gr_series_t x, gr_series_t y, gr_ctx_t ctx)
              void gr_series_set_shallow(gr_series_t res, const gr_series_t x, gr_ctx_t ctx)
              int gr_series_randtest(gr_series_t res, flint_rand_t state, gr_ctx_t ctx)
              int gr_series_write(gr_stream_t out, const gr_series_t x, gr_ctx_t ctx)
              int gr_series_zero(gr_series_t res, gr_ctx_t ctx)
              int gr_series_one(gr_series_t res, gr_ctx_t ctx)
              int gr_series_set(gr_series_t res, const gr_series_t x, gr_ctx_t ctx)
              int gr_series_gen(gr_series_t res, gr_ctx_t ctx)
              int gr_series_gens_recursive(gr_vec_t vec, gr_ctx_t ctx)
              int gr_series_neg(gr_series_t res, const gr_series_t x, gr_ctx_t ctx)
              int gr_series_set_gr_poly(gr_series_t res, const gr_poly_t x, gr_ctx_t ctx)
              int gr_series_set_scalar(gr_series_t res, gr_srcptr x, gr_ctx_t ctx)
              int gr_series_set_si(gr_series_t res, slong c, gr_ctx_t ctx)
              int gr_series_set_ui(gr_series_t res, ulong c, gr_ctx_t ctx)
              int gr_series_set_fmpz(gr_series_t res, const fmpz_t c, gr_ctx_t ctx)
              int gr_series_set_fmpq(gr_series_t res, const fmpq_t c, gr_ctx_t ctx)
              int gr_series_set_other(gr_series_t res, gr_srcptr x, gr_ctx_t x_ctx, gr_ctx_t ctx)
              truth_t gr_series_is_zero(const gr_series_t x, gr_ctx_t ctx)
              truth_t gr_series_is_one(const gr_series_t x, gr_ctx_t ctx)
              truth_t gr_series_coeff_is_zero(const gr_series_t x, slong i, gr_ctx_t ctx)
              truth_t gr_series_equal(const gr_series_t x, const gr_series_t y, gr_ctx_t ctx)
              int gr_series_add(gr_series_t res, const gr_series_t x, const gr_series_t y, gr_ctx_t ctx)
              int gr_series_sub(gr_series_t res, const gr_series_t x, const gr_series_t y, gr_ctx_t ctx)
              int gr_series_mul(gr_series_t res, const gr_series_t x, const gr_series_t y, gr_ctx_t ctx)
              int gr_series_inv(gr_series_t res, const gr_series_t x, gr_ctx_t ctx)
              int gr_series_div(gr_series_t res, const gr_series_t x, const gr_series_t y, gr_ctx_t ctx)
              int gr_series_divexact(gr_series_t res, const gr_series_t x, const gr_series_t y, gr_ctx_t ctx)
              int gr_series_sqrt(gr_series_t res, const gr_series_t x, gr_ctx_t ctx)
              int gr_series_rsqrt(gr_series_t res, const gr_series_t x, gr_ctx_t ctx)
              int gr_series_exp(gr_series_t res, const gr_series_t x, gr_ctx_t ctx)
              int gr_series_log(gr_series_t res, const gr_series_t x, gr_ctx_t ctx)
              int gr_series_tan(gr_series_t res, const gr_series_t x, gr_ctx_t ctx)
              int gr_series_asin(gr_series_t res, const gr_series_t x, gr_ctx_t ctx)
              int gr_series_acos(gr_series_t res, const gr_series_t x, gr_ctx_t ctx)
              int gr_series_atan(gr_series_t res, const gr_series_t x, gr_ctx_t ctx)
              int gr_series_asinh(gr_series_t res, const gr_series_t x, gr_ctx_t ctx)
              int gr_series_acosh(gr_series_t res, const gr_series_t x, gr_ctx_t ctx)
              int gr_series_atanh(gr_series_t res, const gr_series_t x, gr_ctx_t ctx)
              int gr_series_gamma(gr_series_t res, const gr_series_t x, gr_ctx_t ctx)
              int gr_series_rgamma(gr_series_t res, const gr_series_t x, gr_ctx_t ctx)
              int gr_series_lgamma(gr_series_t res, const gr_series_t x, gr_ctx_t ctx)
              int gr_series_digamma(gr_series_t res, const gr_series_t x, gr_ctx_t ctx)
              int gr_series_erf(gr_series_t res, const gr_series_t x, gr_ctx_t ctx)
              int gr_series_erfc(gr_series_t res, const gr_series_t x, gr_ctx_t ctx)
              int gr_series_erfi(gr_series_t res, const gr_series_t x, gr_ctx_t ctx)
              int gr_series_exp_integral_ei(gr_series_t res, const gr_series_t x, gr_ctx_t ctx)
              int gr_series_cos_integral(gr_series_t res, const gr_series_t x, gr_ctx_t ctx)
              int gr_series_cosh_integral(gr_series_t res, const gr_series_t x, gr_ctx_t ctx)
              int gr_series_sin_integral(gr_series_t res, const gr_series_t x, gr_ctx_t ctx)
              int gr_series_sinh_integral(gr_series_t res, const gr_series_t x, gr_ctx_t ctx)
              int gr_series_fresnel(gr_series_t res1, gr_series_t res2, const gr_series_t x, int normalized, gr_ctx_t ctx)
              int gr_series_fresnel_s(gr_series_t res, const gr_series_t x, int normalized, gr_ctx_t ctx)
              int gr_series_fresnel_c(gr_series_t res, const gr_series_t x, int normalized, gr_ctx_t ctx)
              int gr_series_airy(gr_series_t res1, gr_series_t res2, gr_series_t res3, gr_series_t res4, const gr_series_t x, gr_ctx_t ctx)
              int gr_series_airy_ai(gr_series_t res, const gr_series_t x, gr_ctx_t ctx)
              int gr_series_airy_ai_prime(gr_series_t res, const gr_series_t x, gr_ctx_t ctx)
              int gr_series_airy_bi(gr_series_t res, const gr_series_t x, gr_ctx_t ctx)
              int gr_series_airy_bi_prime(gr_series_t res, const gr_series_t x, gr_ctx_t ctx)
              int gr_series_log_integral(gr_series_t res, const gr_series_t x, int offset, gr_ctx_t ctx)
              int gr_series_gamma_upper(gr_series_t res, const gr_series_t s, const gr_series_t x, int regularized, gr_ctx_t ctx)
              int gr_series_gamma_lower(gr_series_t res, const gr_series_t s, const gr_series_t x, int regularized, gr_ctx_t ctx)
              int gr_series_beta_lower(gr_series_t res, const gr_series_t a, const gr_series_t b, const gr_series_t x, int regularized, gr_ctx_t ctx)
              int gr_series_polylog(gr_series_t res, const gr_series_t s, const gr_series_t z, gr_ctx_t ctx)
              int gr_series_hurwitz_zeta(gr_series_t res, const gr_series_t s, const gr_series_t z, gr_ctx_t ctx)
              int gr_series_dirichlet_l(gr_series_t res, const dirichlet_group_t G, const dirichlet_char_t chi, const gr_series_t x, gr_ctx_t ctx)
              int gr_series_dirichlet_hardy_theta(gr_series_t res, const dirichlet_group_t G, const dirichlet_char_t chi, const gr_series_t x, gr_ctx_t ctx)
              int gr_series_dirichlet_hardy_z(gr_series_t res, const dirichlet_group_t G, const dirichlet_char_t chi, const gr_series_t x, gr_ctx_t ctx)
              int gr_series_jacobi_theta(gr_series_t res1, gr_series_t res2, gr_series_t res3, gr_series_t res4, const gr_series_t x, const gr_series_t tau, gr_ctx_t ctx)
              int gr_series_jacobi_theta_1(gr_series_t res, const gr_series_t x, const gr_series_t tau, gr_ctx_t ctx)
              int gr_series_jacobi_theta_2(gr_series_t res, const gr_series_t x, const gr_series_t tau, gr_ctx_t ctx)
              int gr_series_jacobi_theta_3(gr_series_t res, const gr_series_t x, const gr_series_t tau, gr_ctx_t ctx)
              int gr_series_jacobi_theta_4(gr_series_t res, const gr_series_t x, const gr_series_t tau, gr_ctx_t ctx)
              int gr_series_agm1(gr_series_t res, const gr_series_t x, gr_ctx_t ctx)
              int gr_series_elliptic_k(gr_series_t res, const gr_series_t x, gr_ctx_t ctx)
              int gr_series_weierstrass_p(gr_series_t res, const gr_series_t x, const gr_series_t tau, gr_ctx_t ctx)
              int gr_series_hypgeom_pfq(gr_series_t res, const gr_series_vec_t a, const gr_series_vec_t b, const gr_series_t x, int regularized, gr_ctx_t ctx)




.. raw:: latex

    \newpage
