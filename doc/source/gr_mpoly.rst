.. _gr-mpoly:

**gr_mpoly.h** -- sparse multivariate polynomials over generic rings
===============================================================================

A :type:`gr_mpoly_t` represents a multivariate polynomial
`f \in R[X_1,\ldots,X_n]` implemented as an array of coefficients
in a generic ring *R* together with an array of packed exponents.

Weak normalization
-------------------------------------------------------------------------------

A :type:`gr_mpoly_t` is always normalised by removing zero
coefficients.
For rings without decidable equality (e.g. rings with inexact
representation), only coefficients that are provably zero will be
removed, and there can thus be spurious zeros in the
internal representation.
Methods that depend on knowing the exact structure of a polynomial
will act appropriately, typically by returning ``GR_UNABLE``
when it is unknown whether any stored coefficients are nonzero.

Types, macros and constants
-------------------------------------------------------------------------------

.. type:: gr_mpoly_struct

.. type:: gr_mpoly_t

    A ``gr_mpoly_t`` is defined as an array of length one of type
    ``gr_mpoly_struct``, permitting a ``gr_mpoly_t`` to
    be passed by reference.

Memory management
-------------------------------------------------------------------------------

.. function:: void gr_mpoly_init(gr_mpoly_t A, const mpoly_ctx_t mctx, gr_ctx_t cctx)

    Initializes and sets *A* to the zero polynomial.

.. function:: void gr_mpoly_init3(gr_mpoly_t A, slong alloc, flint_bitcnt_t bits, const mpoly_ctx_t mctx, gr_ctx_t cctx)
              void gr_mpoly_init2(gr_mpoly_t A, slong alloc, const mpoly_ctx_t mctx, gr_ctx_t cctx)

    Initializes *A* with space allocated for the given number
    of coefficients and exponents with the given number of bits.

.. function:: void gr_mpoly_clear(gr_mpoly_t A, const mpoly_ctx_t mctx, gr_ctx_t cctx)

    Clears *A*, freeing all allocated data.

Basic manipulation
-------------------------------------------------------------------------------

.. function:: void gr_mpoly_swap(gr_mpoly_t A, gr_mpoly_t B, const mpoly_ctx_t mctx, gr_ctx_t cctx)

    Swaps *A* and *B* efficiently.

.. function:: int gr_mpoly_set(gr_mpoly_t A, const gr_mpoly_t B, const mpoly_ctx_t mctx, gr_ctx_t cctx);

    Sets *A* to *B*.

.. function:: int gr_mpoly_zero(gr_mpoly_t A, const mpoly_ctx_t mctx, gr_ctx_t cctx)

    Sets *A* to the zero polynomial.

.. function:: truth_t gr_mpoly_is_zero(const gr_mpoly_t A, const mpoly_ctx_t mctx, gr_ctx_t cctx)

    Returns whether *A* is the zero polynomial.

.. function:: int gr_mpoly_gen(gr_mpoly_t A, slong var, const mpoly_ctx_t mctx, gr_ctx_t cctx)

    Sets *A* to the generator with index *var* (indexed from zero).

.. function:: truth_t gr_mpoly_is_gen(const gr_mpoly_t A, slong var, const mpoly_ctx_t mctx, gr_ctx_t cctx)

    Returns whether *A* is the generator with index *var* (indexed from zero).

Comparisons
-------------------------------------------------------------------------------

.. function:: truth_t gr_mpoly_equal(const gr_mpoly_t A, const gr_mpoly_t B, const mpoly_ctx_t mctx, gr_ctx_t cctx)

    Returns whether *A* and *B* are equal.

Input and output
-------------------------------------------------------------------------------

.. function:: int gr_mpoly_randtest_bits(gr_mpoly_t A, flint_rand_t state, slong length, flint_bitcnt_t exp_bits, const mpoly_ctx_t mctx, gr_ctx_t cctx)

    Sets *A* to a random polynomial with up to *length* terms
    and up to *exp_bits* bits in the exponents.

Input and output
-------------------------------------------------------------------------------

.. function:: int gr_mpoly_write_pretty(gr_stream_t out, const gr_mpoly_t A, const char ** x, const mpoly_ctx_t mctx, gr_ctx_t cctx)
              int gr_mpoly_print_pretty(const gr_mpoly_t A, const char ** x, const mpoly_ctx_t mctx, gr_ctx_t cctx)

    Prints *A* using the strings in *x* for the variables.
    If *x* is *NULL*, defaults are used.

Coefficient and exponent access
-------------------------------------------------------------------------------

.. function:: int gr_mpoly_get_coeff_scalar_fmpz(gr_ptr c, const gr_mpoly_t A, const fmpz * exp, const mpoly_ctx_t mctx, gr_ctx_t cctx)
              int gr_mpoly_get_coeff_scalar_ui(gr_ptr c, const gr_mpoly_t A, const ulong * exp, const mpoly_ctx_t mctx, gr_ctx_t cctx)

    Sets *c* to the coefficient in *A* with exponents *exp*.

.. function:: int gr_mpoly_set_coeff_scalar_fmpz(gr_mpoly_t A, gr_srcptr c, const fmpz * exp, const mpoly_ctx_t mctx, gr_ctx_t cctx)
              int gr_mpoly_set_coeff_ui_fmpz(gr_mpoly_t A, ulong c, const fmpz * exp, const mpoly_ctx_t mctx, gr_ctx_t cctx)
              int gr_mpoly_set_coeff_si_fmpz(gr_mpoly_t A, slong c, const fmpz * exp, const mpoly_ctx_t mctx, gr_ctx_t cctx)
              int gr_mpoly_set_coeff_fmpz_fmpz(gr_mpoly_t A, const fmpz_t c, const fmpz * exp, const mpoly_ctx_t mctx, gr_ctx_t cctx)
              int gr_mpoly_set_coeff_fmpq_fmpz(gr_mpoly_t A, const fmpq_t c, const fmpz * exp, const mpoly_ctx_t mctx, gr_ctx_t cctx)

.. function:: int gr_mpoly_set_coeff_scalar_ui(gr_mpoly_t poly, gr_srcptr c, const ulong * exp, const mpoly_ctx_t mctx, gr_ctx_t cctx)
              int gr_mpoly_set_coeff_ui_ui(gr_mpoly_t A, ulong c, const ulong * exp, const mpoly_ctx_t mctx, gr_ctx_t cctx)
              int gr_mpoly_set_coeff_si_ui(gr_mpoly_t A, slong c, const ulong * exp, const mpoly_ctx_t mctx, gr_ctx_t cctx)
              int gr_mpoly_set_coeff_fmpz_ui(gr_mpoly_t A, const fmpz_t c, const ulong * exp, const mpoly_ctx_t mctx, gr_ctx_t cctx)
              int gr_mpoly_set_coeff_fmpq_ui(gr_mpoly_t A, const fmpq_t c, const ulong * exp, const mpoly_ctx_t mctx, gr_ctx_t cctx)

    Sets the coefficient with exponents *exp* in *A* to the scalar *c*
    which must be an element of or coercible to the coefficient ring.

Arithmetic
-------------------------------------------------------------------------------

.. function:: int gr_mpoly_neg(gr_mpoly_t A, const gr_mpoly_t B, const mpoly_ctx_t mctx, gr_ctx_t cctx)

    Sets *A* to the negation of *B*.

.. function:: int gr_mpoly_add(gr_mpoly_t A, const gr_mpoly_t B, const gr_mpoly_t C, const mpoly_ctx_t mctx, gr_ctx_t cctx)

    Sets *A* to the difference of *B* and *C*.

.. function:: int gr_mpoly_sub(gr_mpoly_t A, const gr_mpoly_t B, const gr_mpoly_t C, const mpoly_ctx_t mctx, gr_ctx_t cctx)

    Sets *A* to the difference of *B* and *C*.

.. function:: int gr_mpoly_mul(gr_mpoly_t A, const gr_mpoly_t B, const gr_mpoly_t C, const mpoly_ctx_t mctx, gr_ctx_t cctx)
              int gr_mpoly_mul_johnson(gr_mpoly_t A, const gr_mpoly_t B, const gr_mpoly_t C, const mpoly_ctx_t mctx, gr_ctx_t cctx)
              int gr_mpoly_mul_monomial(gr_mpoly_t A, const gr_mpoly_t B, const gr_mpoly_t C, const mpoly_ctx_t mctx, gr_ctx_t cctx)

    Sets *A* to the product of *B* and *C*.
    The *monomial* version assumes that *C* is a monomial.

.. function:: int gr_mpoly_mul_scalar(gr_mpoly_t A, const gr_mpoly_t B, gr_srcptr c, const mpoly_ctx_t mctx, gr_ctx_t cctx)
              int gr_mpoly_mul_si(gr_mpoly_t A, const gr_mpoly_t B, slong c, const mpoly_ctx_t mctx, gr_ctx_t cctx)
              int gr_mpoly_mul_ui(gr_mpoly_t A, const gr_mpoly_t B, ulong c, const mpoly_ctx_t mctx, gr_ctx_t cctx)
              int gr_mpoly_mul_fmpz(gr_mpoly_t A, const gr_mpoly_t B, const fmpz_t c, const mpoly_ctx_t mctx, gr_ctx_t cctx)
              int gr_mpoly_mul_fmpq(gr_mpoly_t A, const gr_mpoly_t B, const fmpq_t c, const mpoly_ctx_t mctx, gr_ctx_t cctx)

    Sets *A* to *B* multiplied by the scalar *c* which must be
    an element of or coercible to the coefficient ring.

Container operations
-------------------------------------------------------------------------------

Mostly intended for internal use.

.. function:: void _gr_mpoly_fit_length(gr_ptr * coeffs, slong * coeffs_alloc, ulong ** exps, slong * exps_alloc, slong N, slong length, gr_ctx_t cctx)

.. function:: void gr_mpoly_fit_length(gr_mpoly_t A, slong len, const mpoly_ctx_t mctx, gr_ctx_t cctx);

    Ensures that *A* has space for *len* coefficients and exponents.

.. function:: void gr_mpoly_fit_bits(gr_mpoly_t A, flint_bitcnt_t bits, const mpoly_ctx_t mctx, gr_ctx_t cctx)

.. function:: void gr_mpoly_fit_length_fit_bits(gr_mpoly_t A, slong len, flint_bitcnt_t bits, const mpoly_ctx_t mctx, gr_ctx_t cctx)

.. function:: void gr_mpoly_fit_length_reset_bits(gr_mpoly_t A, slong len, flint_bitcnt_t bits, const mpoly_ctx_t mctx, gr_ctx_t cctx)

.. function:: void _gr_mpoly_set_length(gr_mpoly_t A, slong newlen, const mpoly_ctx_t mctx, gr_ctx_t cctx)

.. function:: void _gr_mpoly_push_exp_ui(gr_mpoly_t A, const ulong * exp, const mpoly_ctx_t mctx, gr_ctx_t cctx)

.. function:: int gr_mpoly_push_term_scalar_ui(gr_mpoly_t A, gr_srcptr c, const ulong * exp, const mpoly_ctx_t mctx, gr_ctx_t cctx)

.. function:: void _gr_mpoly_push_exp_fmpz(gr_mpoly_t A, const fmpz * exp, const mpoly_ctx_t mctx, gr_ctx_t cctx)

.. function:: int gr_mpoly_push_term_scalar_fmpz(gr_mpoly_t A, gr_srcptr c, const fmpz * exp, const mpoly_ctx_t mctx, gr_ctx_t cctx)

.. function:: void gr_mpoly_sort_terms(gr_mpoly_t A, const mpoly_ctx_t mctx, gr_ctx_t cctx)

.. function:: int gr_mpoly_combine_like_terms(gr_mpoly_t A, const mpoly_ctx_t mctx, gr_ctx_t cctx)

.. function:: truth_t gr_mpoly_is_canonical(const gr_mpoly_t A, const mpoly_ctx_t mctx, gr_ctx_t cctx)

.. function:: void gr_mpoly_assert_canonical(const gr_mpoly_t A, const mpoly_ctx_t mctx, gr_ctx_t cctx)



.. raw:: latex

    \newpage
