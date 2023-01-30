.. _gr-poly:

**gr_poly.h** -- dense univariate polynomials over generic rings
===============================================================================

A :type:`gr_poly_t` represents a univariate polynomial `f \in R[X]`
implemented as a dense array of coefficients in a generic ring *R*.

In this module, the context object ``ctx`` always represents the
coefficient ring *R* unless otherwise stated.
Creating a context object representing the polynomial ring `R[X]`
only becomes necessary when one
wants to manipulate polynomials using generic ring methods
like ``gr_add`` instead of the designated polynomial
methods like ``gr_poly_add``.

Most functions are provided in two versions: an underscore method which
operates directly on pre-allocated arrays of coefficients and generally
has some restrictions (often requiring the lengths to be nonzero
and not supporting aliasing of the input and output arrays),
and a non-underscore method which performs automatic memory
management and handles degenerate cases.

Type compatibility
-------------------------------------------------------------------------------

The ``gr_poly`` type has the same data layout as the following
polynomial types: ``fmpz_poly``, ``fq_poly``, ``fq_nmod_poly``,
``fq_zech_poly``, ``arb_poly``, ``acb_poly``, ``ca_poly``.
Methods in this module can therefore be mixed freely with
methods in the corresponding Flint, Arb and Calcium modules
when the underlying coefficient type is the same.

It is not directly compatible with the following types:
``fmpq_poly`` (coefficients are stored with a common denominator),
``nmod_poly`` (modulus data is stored as part of the polynomial object).

Weak normalization
-------------------------------------------------------------------------------

A :type:`gr_poly_t` is always normalised by removing leading zeros.
For rings without decidable equality (e.g. rings with inexact
representation), only coefficients that are provably zero will be
removed, and there can thus be spurious leading zeros in the
internal representation.
Methods that depend on knowing the exact degree of a polynomial
will act appropriately, typically by returning ``GR_UNABLE``
when it is unknown whether the leading stored coefficient is nonzero.

Types, macros and constants
-------------------------------------------------------------------------------

.. type:: gr_poly_struct

.. type:: gr_poly_t

    Contains a pointer to an array of coefficients (``coeffs``), the used
    length (``length``), and the allocated size of the array (``alloc``).

    A ``gr_poly_t`` is defined as an array of length one of type
    ``gr_poly_struct``, permitting a ``gr_poly_t`` to
    be passed by reference.

Memory management
-------------------------------------------------------------------------------

.. function:: void gr_poly_init(gr_poly_t poly, gr_ctx_t ctx)

.. function:: void gr_poly_init2(gr_poly_t poly, slong len, gr_ctx_t ctx)

.. function:: void gr_poly_clear(gr_poly_t poly, gr_ctx_t ctx)

.. function:: gr_ptr gr_poly_entry_ptr(gr_poly_t poly, slong i, gr_ctx_t ctx)

.. function:: slong gr_poly_length(const gr_poly_t poly, gr_ctx_t ctx)

.. function:: void gr_poly_swap(gr_poly_t poly1, gr_poly_t poly2, gr_ctx_t ctx)

.. function:: void gr_poly_fit_length(gr_poly_t poly, slong len, gr_ctx_t ctx)

.. function:: void _gr_poly_set_length(gr_poly_t poly, slong len, gr_ctx_t ctx)

Basic manipulation
-------------------------------------------------------------------------------

.. function:: void _gr_poly_normalise(gr_poly_t poly, gr_ctx_t ctx)

.. function:: int gr_poly_set(gr_poly_t res, const gr_poly_t src, gr_ctx_t ctx)
              int gr_poly_get_fmpz_poly(gr_poly_t res, const fmpz_poly_t src, gr_ctx_t ctx)
              int gr_poly_set_fmpq_poly(gr_poly_t res, const fmpq_poly_t src, gr_ctx_t ctx)
              int gr_poly_set_gr_poly_other(gr_poly_t res, const gr_poly_t x, gr_ctx_t x_ctx, gr_ctx_t ctx)

.. function:: int _gr_poly_reverse(gr_ptr res, gr_srcptr poly, slong len, slong n, gr_ctx_t ctx)
              int gr_poly_reverse(gr_poly_t res, const gr_poly_t poly, slong n, gr_ctx_t ctx)

.. function:: int gr_poly_truncate(gr_poly_t poly, slong newlen, gr_ctx_t ctx)

.. function:: int gr_poly_zero(gr_poly_t poly, gr_ctx_t ctx)
              int gr_poly_one(gr_poly_t poly, gr_ctx_t ctx)
              int gr_poly_neg_one(gr_poly_t poly, gr_ctx_t ctx)

.. function:: int gr_poly_write(gr_stream_t out, const gr_poly_t poly, gr_ctx_t ctx)
              int gr_poly_print(const gr_poly_t poly, gr_ctx_t ctx)

.. function:: int gr_poly_randtest(gr_poly_t poly, flint_rand_t state, slong len, gr_ctx_t ctx)

.. function:: truth_t _gr_poly_equal(gr_srcptr poly1, slong len1, gr_srcptr poly2, slong len2, gr_ctx_t ctx)
              truth_t gr_poly_equal(const gr_poly_t poly1, const gr_poly_t poly2, gr_ctx_t ctx)

.. function:: truth_t gr_poly_is_zero(const gr_poly_t poly, gr_ctx_t ctx)
              truth_t gr_poly_is_one(const gr_poly_t poly, gr_ctx_t ctx)

.. function:: int gr_poly_set_scalar(gr_poly_t poly, gr_srcptr c, gr_ctx_t ctx)
              int gr_poly_set_si(gr_poly_t poly, slong c, gr_ctx_t ctx)
              int gr_poly_set_ui(gr_poly_t poly, slong c, gr_ctx_t ctx)
              int gr_poly_set_fmpz(gr_poly_t poly, const fmpz_t c, gr_ctx_t ctx)
              int gr_poly_set_fmpq(gr_poly_t poly, const fmpq_t c, gr_ctx_t ctx)

.. function:: int gr_poly_set_coeff_scalar(gr_poly_t poly, slong n, gr_srcptr c, gr_ctx_t ctx)
              int gr_poly_set_coeff_si(gr_poly_t poly, slong n, slong c, gr_ctx_t ctx)
              int gr_poly_set_coeff_ui(gr_poly_t poly, slong n, ulong c, gr_ctx_t ctx)
              int gr_poly_set_coeff_fmpz(gr_poly_t poly, slong n, const fmpz_t c, gr_ctx_t ctx)
              int gr_poly_set_coeff_fmpq(gr_poly_t poly, slong n, const fmpq_t c, gr_ctx_t ctx)

.. function:: int gr_poly_get_coeff_scalar(gr_ptr res, const gr_poly_t poly, slong n, gr_ctx_t ctx)

Arithmetic
-------------------------------------------------------------------------------

.. function:: int gr_poly_neg(gr_poly_t res, const gr_poly_t src, gr_ctx_t ctx)

.. function:: int _gr_poly_add(gr_ptr res, gr_srcptr poly1, slong len1, gr_srcptr poly2, slong len2, gr_ctx_t ctx)
              int gr_poly_add(gr_poly_t res, const gr_poly_t poly1, const gr_poly_t poly2, gr_ctx_t ctx)

.. function:: int _gr_poly_sub(gr_ptr res, gr_srcptr poly1, slong len1, gr_srcptr poly2, slong len2, gr_ctx_t ctx)
              int gr_poly_sub(gr_poly_t res, const gr_poly_t poly1, const gr_poly_t poly2, gr_ctx_t ctx)

.. function:: int _gr_poly_mul(gr_ptr res, gr_srcptr poly1, slong len1, gr_srcptr poly2, slong len2, gr_ctx_t ctx)
              int gr_poly_mul(gr_poly_t res, const gr_poly_t poly1, const gr_poly_t poly2, gr_ctx_t ctx)

.. function:: int _gr_poly_mullow_generic(gr_ptr res, gr_srcptr poly1, slong len1, gr_srcptr poly2, slong len2, slong len, gr_ctx_t ctx)
              int gr_poly_mullow(gr_poly_t res, const gr_poly_t poly1, const gr_poly_t poly2, slong len, gr_ctx_t ctx)

.. function:: int gr_poly_mul_scalar(gr_poly_t res, const gr_poly_t poly, gr_srcptr c, gr_ctx_t ctx)

Powering
--------------------------------------------------------------------------------

.. function:: int _gr_poly_pow_series_ui_binexp(gr_ptr res, gr_srcptr f, slong flen, ulong exp, slong len, gr_ctx_t ctx)
              int gr_poly_pow_series_ui_binexp(gr_poly_t res, const gr_poly_t poly, ulong exp, slong len, gr_ctx_t ctx)

.. function:: int _gr_poly_pow_series_ui(gr_ptr res, gr_srcptr f, slong flen, ulong exp, slong len, gr_ctx_t ctx)
              int gr_poly_pow_series_ui(gr_poly_t res, const gr_poly_t poly, ulong exp, slong len, gr_ctx_t ctx)

.. function:: int _gr_poly_pow_ui_binexp(gr_ptr res, gr_srcptr f, slong flen, ulong exp, gr_ctx_t ctx);
              int gr_poly_pow_ui_binexp(gr_poly_t res, const gr_poly_t poly, ulong exp, gr_ctx_t ctx);

.. function:: int _gr_poly_pow_ui(gr_ptr res, gr_srcptr f, slong flen, ulong exp, gr_ctx_t ctx);
              int gr_poly_pow_ui(gr_poly_t res, const gr_poly_t poly, ulong exp, gr_ctx_t ctx);

.. function:: int gr_poly_pow_fmpz(gr_poly_t res, const gr_poly_t poly, const fmpz_t exp, gr_ctx_t ctx);

.. function:: int _gr_poly_pow_series_fmpq_recurrence(gr_ptr h, gr_srcptr f, slong flen, const fmpq_t exp, slong len, gr_ctx_t ctx)
              int gr_poly_pow_series_fmpq_recurrence(gr_poly_t res, const gr_poly_t poly, const fmpq_t exp, slong len, gr_ctx_t ctx)

Division
--------------------------------------------------------------------------------

TODO: algorithms handle allocation for R differently

.. function:: int _gr_poly_divrem_divconquer_preinv(gr_ptr Q, gr_ptr R, gr_srcptr A, slong lenA, gr_srcptr B, slong lenB, gr_srcptr invB, slong cutoff, gr_ctx_t ctx)
              int _gr_poly_divrem_divconquer(gr_ptr Q, gr_ptr R, gr_srcptr A, slong lenA, gr_srcptr B, slong lenB, slong cutoff, gr_ctx_t ctx)
              int gr_poly_divrem_divconquer(gr_poly_t Q, gr_poly_t R, const gr_poly_t A, const gr_poly_t B, slong cutoff, gr_ctx_t ctx)
              int _gr_poly_divrem_basecase_preinv(gr_ptr Q, gr_ptr R, gr_srcptr A, slong lenA, gr_srcptr B, slong lenB, gr_srcptr invB, gr_ctx_t ctx)
              int _gr_poly_divrem_basecase_noinv(gr_ptr Q, gr_ptr R, gr_srcptr A, slong lenA, gr_srcptr B, slong lenB, gr_ctx_t ctx)
              int _gr_poly_divrem_basecase(gr_ptr Q, gr_ptr R, gr_srcptr A, slong lenA, gr_srcptr B, slong lenB, gr_ctx_t ctx)
              int gr_poly_divrem_basecase(gr_poly_t Q, gr_poly_t R, const gr_poly_t A, const gr_poly_t B, gr_ctx_t ctx)
              int _gr_poly_divrem(gr_ptr Q, gr_ptr R, gr_srcptr A, slong lenA, gr_srcptr B, slong lenB, gr_ctx_t ctx)
              int gr_poly_divrem(gr_poly_t Q, gr_poly_t R, const gr_poly_t A, const gr_poly_t B, gr_ctx_t ctx)

    Polynomial division with remainder. ``GR_DOMAIN`` is returned when
    *B* is provably zero or when encountering an impossible division
    in the polynomial division algorithm.

    The underscore methods make the following assumptions:

    * *Q* has room for ``lenA - lenB + 1`` coefficients.
    * *R* has room for ``lenA`` coefficients.
    * ``lenA >= lenB >= 1``.
    * *Q* is not aliased with either *A* or *B*.
    * *R* is not aliased with *B*.
    * The divisor *B* is normalized to have nonzero leading coefficient.
      (The non-underscore methods check for leading coefficients that
      are not provably nonzero and return ``GR_UNABLE``)

    The ``preinv`` functions take a precomputed inverse of the
    leading coefficient as input.

.. function:: int _gr_poly_div_newton(gr_ptr Q, gr_srcptr A, slong lenA, gr_srcptr B, slong lenB, gr_ctx_t ctx);
              int _gr_poly_divrem_newton(gr_ptr Q, gr_ptr R, gr_srcptr A, slong lenA, gr_srcptr B, slong lenB, gr_ctx_t ctx);
              int gr_poly_divrem_newton(gr_poly_t Q, gr_poly_t R, const gr_poly_t A, const gr_poly_t B, gr_ctx_t ctx);


Power series division
--------------------------------------------------------------------------------

For divide-and-conquer (including Newton-like) algorithms, *cutoff* has the
following meaning: we use the basecase algorithm for lengths `n < \operatorname{cutoff}`
and the divide-and-conquer algorithm for `n \ge \operatorname{cutoff}`.
Using `\operatorname{cutoff} = n` thus results in exactly one divide-and-conquer
step with a basecase length of `\lceil n / 2 \rceil`.
One should **avoid** calling the Newton methods with `n < \operatorname{cutoff}`
as this may result in much worse performance if those methods
do not have a specific escape check for that case.

The *newton* versions uses Newton iteration, switching to a basecase
algorithm when the length is smaller than the specified *cutoff*.
Division uses the Karp-Markstein algorithm.

.. function:: int _gr_poly_inv_series_newton(gr_ptr res, gr_srcptr A, slong Alen, slong len, slong cutoff, gr_ctx_t ctx)
              int gr_poly_inv_series_newton(gr_poly_t res, const gr_poly_t A, slong len, slong cutoff, gr_ctx_t ctx)
              int _gr_poly_inv_series_basecase(gr_ptr res, gr_srcptr A, slong Alen, slong len, gr_ctx_t ctx)
              int gr_poly_inv_series_basecase(gr_poly_t res, const gr_poly_t A, slong len, gr_ctx_t ctx)
              int _gr_poly_inv_series(gr_ptr res, gr_srcptr A, slong Alen, slong len, gr_ctx_t ctx)
              int gr_poly_inv_series(gr_poly_t res, const gr_poly_t A, slong len, gr_ctx_t ctx)

.. function:: int _gr_poly_div_series_newton(gr_ptr res, gr_srcptr A, slong Alen, gr_srcptr B, slong Blen, slong len, gr_ctx_t ctx)
              int gr_poly_div_series_newton(gr_poly_t res, const gr_poly_t A, const gr_poly_t B, slong len, gr_ctx_t ctx)
              int _gr_poly_div_series_basecase(gr_ptr res, gr_srcptr A, slong Alen, gr_srcptr B, slong Blen, slong len, gr_ctx_t ctx)
              int gr_poly_div_series_basecase(gr_poly_t res, const gr_poly_t A, const gr_poly_t B, slong len, gr_ctx_t ctx)
              int _gr_poly_div_series(gr_ptr res, gr_srcptr A, slong Alen, gr_srcptr B, slong Blen, slong len, gr_ctx_t ctx)
              int gr_poly_div_series(gr_poly_t res, const gr_poly_t A, const gr_poly_t B, slong len, gr_ctx_t ctx)

Square roots
--------------------------------------------------------------------------------

.. function:: int _gr_poly_sqrt_series_newton(gr_ptr res, gr_srcptr f, slong flen, slong len, slong cutoff, gr_ctx_t ctx)
              int gr_poly_sqrt_series_newton(gr_poly_t res, const gr_poly_t f, slong len, slong cutoff, gr_ctx_t ctx)
              int _gr_poly_sqrt_series_basecase(gr_ptr res, gr_srcptr f, slong flen, slong len, gr_ctx_t ctx)
              int gr_poly_sqrt_series_basecase(gr_poly_t res, const gr_poly_t f, slong len, gr_ctx_t ctx)
              int _gr_poly_sqrt_series_miller(gr_ptr res, gr_srcptr f, slong flen, slong len, gr_ctx_t ctx)
              int gr_poly_sqrt_series_miller(gr_poly_t res, const gr_poly_t f, slong len, gr_ctx_t ctx)
              int _gr_poly_sqrt_series(gr_ptr res, gr_srcptr f, slong flen, slong len, gr_ctx_t ctx)
              int gr_poly_sqrt_series(gr_poly_t res, const gr_poly_t f, slong len, gr_ctx_t ctx)

.. function:: int _gr_poly_rsqrt_series_newton(gr_ptr res, gr_srcptr f, slong flen, slong len, slong cutoff, gr_ctx_t ctx)
              int gr_poly_rsqrt_series_newton(gr_poly_t res, const gr_poly_t f, slong len, slong cutoff, gr_ctx_t ctx)
              int _gr_poly_rsqrt_series_basecase(gr_ptr res, gr_srcptr f, slong flen, slong len, gr_ctx_t ctx)
              int gr_poly_rsqrt_series_basecase(gr_poly_t res, const gr_poly_t f, slong len, gr_ctx_t ctx)
              int _gr_poly_rsqrt_series_miller(gr_ptr res, gr_srcptr f, slong flen, slong len, gr_ctx_t ctx)
              int gr_poly_rsqrt_series_miller(gr_poly_t res, const gr_poly_t f, slong len, gr_ctx_t ctx)
              int _gr_poly_rsqrt_series(gr_ptr res, gr_srcptr f, slong flen, slong len, gr_ctx_t ctx)
              int gr_poly_rsqrt_series(gr_poly_t res, const gr_poly_t f, slong len, gr_ctx_t ctx)

Evaluation
-------------------------------------------------------------------------------

.. function:: int _gr_poly_evaluate_rectangular(gr_ptr res, gr_srcptr poly, slong len, gr_srcptr x, gr_ctx_t ctx)
              int gr_poly_evaluate_rectangular(gr_ptr res, const gr_poly_t poly, gr_srcptr x, gr_ctx_t ctx)

.. function:: int _gr_poly_evaluate_horner(gr_ptr res, gr_srcptr poly, slong len, gr_srcptr x, gr_ctx_t ctx)
              int gr_poly_evaluate_horner(gr_ptr res, const gr_poly_t poly, gr_srcptr x, gr_ctx_t ctx)

.. function:: int _gr_poly_evaluate(gr_ptr res, gr_srcptr poly, slong len, gr_srcptr x, gr_ctx_t ctx)
              int gr_poly_evaluate(gr_ptr res, const gr_poly_t poly, gr_srcptr x, gr_ctx_t ctx)

    Set *res* to *poly* evaluated at *x*.

.. function:: int _gr_poly_evaluate_other_horner(gr_ptr res, gr_srcptr f, slong len, const gr_srcptr x, gr_ctx_t x_ctx, gr_ctx_t ctx)
              int gr_poly_evaluate_other_horner(gr_ptr res, const gr_poly_t f, gr_srcptr x, gr_ctx_t x_ctx, gr_ctx_t ctx)
              int _gr_poly_evaluate_other_rectangular(gr_ptr res, gr_srcptr f, slong len, const gr_srcptr x, gr_ctx_t x_ctx, gr_ctx_t ctx)
              int gr_poly_evaluate_other_rectangular(gr_ptr res, const gr_poly_t f, gr_srcptr x, gr_ctx_t x_ctx, gr_ctx_t ctx)
              int _gr_poly_evaluate_other(gr_ptr res, gr_srcptr f, slong len, const gr_srcptr x, gr_ctx_t x_ctx, gr_ctx_t ctx)
              int gr_poly_evaluate_other(gr_ptr res, const gr_poly_t f, gr_srcptr x, gr_ctx_t x_ctx, gr_ctx_t ctx)

    Set *res* to *poly* evaluated at *x*, where the coefficients of *f*
    belong to *ctx* while both *x* and *res* belong to *x_ctx*.

Multipoint evaluation and interpolation
-------------------------------------------------------------------------------

.. function:: gr_ptr * _gr_poly_tree_alloc(slong len, gr_ctx_t ctx)

.. function:: void _gr_poly_tree_free(gr_ptr * tree, slong len, gr_ctx_t ctx)

.. function:: int _gr_poly_tree_build(gr_ptr * tree, gr_srcptr roots, slong len, gr_ctx_t ctx)

.. function:: int _gr_poly_evaluate_vec_fast_precomp(gr_ptr vs, gr_srcptr poly, slong plen, gr_ptr * tree, slong len, gr_ctx_t ctx)

.. function:: int _gr_poly_evaluate_vec_fast(gr_ptr ys, gr_srcptr poly, slong plen, gr_srcptr xs, slong n, gr_ctx_t ctx)

.. function:: int gr_poly_evaluate_vec_fast(gr_vec_t ys, const gr_poly_t poly, const gr_vec_t xs, gr_ctx_t ctx)

.. function:: int _gr_poly_evaluate_vec_iter(gr_ptr ys, gr_srcptr poly, slong plen, gr_srcptr xs, slong n, gr_ctx_t ctx)

.. function:: int gr_poly_evaluate_vec_iter(gr_vec_t ys, const gr_poly_t poly, const gr_vec_t xs, gr_ctx_t ctx)


Composition
-------------------------------------------------------------------------------

.. function:: int _gr_poly_taylor_shift_horner(gr_ptr res, gr_srcptr poly, slong len, gr_srcptr c, gr_ctx_t ctx)
              int gr_poly_taylor_shift_horner(gr_poly_t res, const gr_poly_t poly, gr_srcptr c, gr_ctx_t ctx)
              int _gr_poly_taylor_shift_divconquer(gr_ptr res, gr_srcptr poly, slong len, gr_srcptr c, gr_ctx_t ctx)
              int gr_poly_taylor_shift_divconquer(gr_poly_t res, const gr_poly_t poly, gr_srcptr c, gr_ctx_t ctx)
              int _gr_poly_taylor_shift(gr_ptr res, gr_srcptr poly, slong len, gr_srcptr c, gr_ctx_t ctx)
              int gr_poly_taylor_shift(gr_poly_t res, const gr_poly_t poly, gr_srcptr c, gr_ctx_t ctx)

.. function:: int _gr_poly_compose_horner(gr_ptr res, gr_srcptr poly1, slong len1, gr_srcptr poly2, slong len2, gr_ctx_t ctx)
              int gr_poly_compose_horner(gr_poly_t res, const gr_poly_t poly1, const gr_poly_t poly2, gr_ctx_t ctx)
              int _gr_poly_compose_divconquer(gr_ptr res, gr_srcptr poly1, slong len1, gr_srcptr poly2, slong len2, gr_ctx_t ctx)
              int gr_poly_compose_divconquer(gr_poly_t res, const gr_poly_t poly1, const gr_poly_t poly2, gr_ctx_t ctx)
              int _gr_poly_compose(gr_ptr res, gr_srcptr poly1, slong len1, gr_srcptr poly2, slong len2, gr_ctx_t ctx)
              int gr_poly_compose(gr_poly_t res, const gr_poly_t poly1, const gr_poly_t poly2, gr_ctx_t ctx)

Derivative and integral
-------------------------------------------------------------------------------

.. function:: int _gr_poly_derivative(gr_ptr res, gr_srcptr poly, slong len, gr_ctx_t ctx)
              int gr_poly_derivative(gr_poly_t res, const gr_poly_t poly, gr_ctx_t ctx)

.. function:: int _gr_poly_integral(gr_ptr res, gr_srcptr poly, slong len, gr_ctx_t ctx)
              int gr_poly_integral(gr_poly_t res, const gr_poly_t poly, gr_ctx_t ctx)

Monic polynomials
-------------------------------------------------------------------------------

.. function:: int _gr_poly_make_monic(gr_ptr res, gr_srcptr poly, slong len, gr_ctx_t ctx)
              int gr_poly_make_monic(gr_poly_t res, const gr_poly_t src, gr_ctx_t ctx)

.. function:: truth_t _gr_poly_is_monic(gr_srcptr poly, slong len, gr_ctx_t ctx)
              truth_t gr_poly_is_monic(const gr_poly_t res, gr_ctx_t ctx)

GCD
-------------------------------------------------------------------------------

.. function:: int _gr_poly_gcd_euclidean(gr_ptr G, slong * lenG, gr_srcptr A, slong lenA, gr_srcptr B, slong lenB, gr_ctx_t ctx)
              int gr_poly_gcd_euclidean(gr_poly_t G, const gr_poly_t A, const gr_poly_t B, gr_ctx_t ctx)
              int _gr_poly_gcd(gr_ptr G, slong * lenG, gr_srcptr A, slong lenA, gr_srcptr B, slong lenB, gr_ctx_t ctx)
              int gr_poly_gcd(gr_poly_t G, const gr_poly_t A, const gr_poly_t B, gr_ctx_t ctx)

    Polynomial GCD. Currently only useful over fields.

    The underscore methods assume ``lenA >= lenB >= 1`` and that both
    *A* and *B* have nonzero leading coefficient.

Roots
-------------------------------------------------------------------------------

.. function:: int gr_poly_roots(gr_vec_t roots, gr_vec_t mult, const gr_poly_t poly, int flags, gr_ctx_t ctx)
              int gr_poly_roots_other(gr_vec_t roots, gr_vec_t mult, const gr_poly_t poly, gr_ctx_t poly_ctx, int flags, gr_ctx_t ctx)

    Finds all roots of the given polynomial in the ring defined by *ctx*,
    storing the roots without duplication in *roots* (a vector with
    elements of type ``ctx``) and the corresponding multiplicities in
    *mult* (a vector with elements of type ``fmpz``).

    If the target ring is not an algebraically closed field, then
    the sum of multiplicities can be smaller than the degree of the
    polynomial. For example, with ``fmpz`` coefficients, we only
    find integer roots.
    The *other* version of this function takes as input a polynomial
    with entries in a different ring ``poly_ctx``. For example,
    we can compute ``qqbar`` or ``arb`` roots for a polynomial
    with ``fmpz`` coefficients.

    Whether the roots are sorted in any particular order is
    ring-dependent (and possibly undefined for a given ring).

    We consider roots of the zero polynomial to be ill-defined and return
    ``GR_DOMAIN`` in that case.

Power series special functions
--------------------------------------------------------------------------------

.. function:: int _gr_poly_atan_series(gr_ptr res, gr_srcptr A, slong Alen, slong len, gr_ctx_t ctx)
              int gr_poly_atan_series(gr_poly_t res, const gr_poly_t A, slong len, gr_ctx_t ctx)
              int _gr_poly_atanh_series(gr_ptr res, gr_srcptr A, slong Alen, slong len, gr_ctx_t ctx)
              int gr_poly_atanh_series(gr_poly_t res, const gr_poly_t A, slong len, gr_ctx_t ctx)
              int _gr_poly_log_series(gr_ptr res, gr_srcptr f, slong flen, slong len, gr_ctx_t ctx)
              int gr_poly_log_series(gr_poly_t res, const gr_poly_t f, slong len, gr_ctx_t ctx)

.. raw:: latex

    \newpage

