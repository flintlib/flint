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

.. function:: int _gr_poly_reverse(gr_ptr res, gr_srcptr poly, slong len, slong n, gr_ctx_t ctx)
              int gr_poly_reverse(gr_poly_t res, const gr_poly_t poly, slong n, gr_ctx_t ctx)

.. function:: int gr_poly_zero(gr_poly_t poly, gr_ctx_t ctx)
              int gr_poly_one(gr_poly_t poly, gr_ctx_t ctx)
              int gr_poly_neg_one(gr_poly_t poly, gr_ctx_t ctx)

.. function:: int gr_poly_write(gr_stream_t out, const gr_poly_t poly, gr_ctx_t ctx)
              int gr_poly_print(const gr_poly_t poly, gr_ctx_t ctx)

.. function:: int gr_poly_randtest(gr_poly_t poly, flint_rand_t state, slong len, gr_ctx_t ctx)

.. function:: truth_t _gr_poly_equal(gr_srcptr poly1, slong len1, gr_srcptr poly2, slong len2, gr_ctx_t ctx)
              truth_t gr_poly_equal(const gr_poly_t poly1, const gr_poly_t poly2, gr_ctx_t ctx)

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

.. function:: int _gr_poly_inv_series(gr_ptr res, gr_srcptr f, slong flen, slong len, gr_ctx_t ctx)
              int gr_poly_inv_series(gr_poly_t res, const gr_poly_t f, slong len, gr_ctx_t ctx)

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

.. raw:: latex

    \newpage

