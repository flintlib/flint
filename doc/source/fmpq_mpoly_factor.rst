.. _fmpq-mpoly-factor:

**fmpq_mpoly_factor.h** -- factorisation of multivariate polynomials over the rational numbers
==============================================================================================

Types, macros and constants
-------------------------------------------------------------------------------

.. type:: fmpq_mpoly_factor_struct

    A struct for holding a factored rational polynomial. There is a
    single constant and a product of bases to corresponding exponents.

.. type:: fmpq_mpoly_factor_t

    An array of length `1` of ``fmpq_mpoly_factor_struct``.


Memory management
--------------------------------------------------------------------------------


.. function:: void fmpq_mpoly_factor_init(fmpq_mpoly_factor_t f, const fmpq_mpoly_ctx_t ctx)

    Initialise *f*.

.. function:: void fmpq_mpoly_factor_clear(fmpq_mpoly_factor_t f, const fmpq_mpoly_ctx_t ctx)

    Clear *f*.


Basic manipulation
--------------------------------------------------------------------------------


.. function:: void fmpq_mpoly_factor_swap(fmpq_mpoly_factor_t f, fmpq_mpoly_factor_t g, const fmpq_mpoly_ctx_t ctx)

    Efficiently swap *f* and *g*.

.. function:: slong fmpq_mpoly_factor_length(const fmpq_mpoly_factor_t f, const fmpq_mpoly_ctx_t ctx)

    Return the length of the product in *f*.

.. function:: void fmpq_mpoly_factor_get_constant_fmpq(fmpq_t c, const fmpq_mpoly_factor_t f, const fmpq_mpoly_ctx_t ctx)

    Set *c* to the constant of *f*.

.. function:: void fmpq_mpoly_factor_get_base(fmpq_mpoly_t B, const fmpq_mpoly_factor_t f, slong i, const fmpq_mpoly_ctx_t ctx)
              void fmpq_mpoly_factor_swap_base(fmpq_mpoly_t B, fmpq_mpoly_factor_t f, slong i, const fmpq_mpoly_ctx_t ctx)

    Set (resp. swap) *B* to (resp. with) the base of the term of index *i* in  *A*.

.. function:: slong fmpq_mpoly_factor_get_exp_si(fmpq_mpoly_factor_t f, slong i, const fmpq_mpoly_ctx_t ctx)

    Return the exponent of the term of index *i* in *A*. It is assumed to fit an ``slong``.

.. function:: void fmpq_mpoly_factor_sort(fmpq_mpoly_factor_t f, const fmpq_mpoly_ctx_t ctx)

    Sort the product of *f* first by exponent and then by base.

.. function:: int fmpq_mpoly_factor_make_monic(fmpq_mpoly_factor_t f, const fmpq_mpoly_ctx_t ctx)
              int fmpq_mpoly_factor_make_integral(fmpq_mpoly_factor_t f, const fmpq_mpoly_ctx_t ctx)

    Make the bases in *f* monic (resp. integral and primitive with positive leading coefficient).
    Return `1` for success, `0` for failure.


Factorisation
--------------------------------------------------------------------------------

    A return of `1` indicates that the function was successful. Otherwise,
    the return is `0` and *f* is undefined. None of these functions
    multiply *f* by *A*: *f* is simply set to a factorisation of *A*, and thus
    these functions should not depend on the initial value of the output *f*.
    The normalization of the factors is not yet specified: use :func:`fmpq_mpoly_factor_make_monic`
    or :func:`fmpq_mpoly_factor_make_integral` for common normalizations.

.. function:: int fmpq_mpoly_factor_squarefree(fmpq_mpoly_factor_t f, const fmpq_mpoly_t A, const fmpq_mpoly_ctx_t ctx)

    Set *f* to a factorization of *A* where the bases are primitive and
    pairwise relatively prime. If the product of all irreducible factors with
    a given exponent is desired, it is recommend to call :func:`fmpq_mpoly_factor_sort`
    and then multiply the bases with the desired exponent.

.. function:: int fmpq_mpoly_factor(fmpq_mpoly_factor_t f, const fmpq_mpoly_t A, const fmpq_mpoly_ctx_t ctx)

    Set *f* to a factorization of *A* where the bases are irreducible.

