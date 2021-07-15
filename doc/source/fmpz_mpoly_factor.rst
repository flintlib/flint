.. _fmpz-mpoly-factor:

**fmpz_mpoly_factor.h** -- factorisation of multivariate polynomials over the integers
======================================================================================

Types, macros and constants
-------------------------------------------------------------------------------

.. type:: fmpz_mpoly_factor_struct

    A struct for holding a factored integer polynomial. There is a
    single constant and a product of bases to corresponding exponents.

.. type:: fmpz_mpoly_factor_t

    An array of length `1` of ``fmpz_mpoly_factor_struct``.


Memory management
--------------------------------------------------------------------------------


.. function:: void fmpz_mpoly_factor_init(fmpz_mpoly_factor_t f, const fmpz_mpoly_ctx_t ctx)

    Initialise *f*.

.. function:: void fmpz_mpoly_factor_clear(fmpz_mpoly_factor_t f, const fmpz_mpoly_ctx_t ctx)

    Clear *f*.


Basic manipulation
--------------------------------------------------------------------------------


.. function:: void fmpz_mpoly_factor_swap(fmpz_mpoly_factor_t f, fmpz_mpoly_factor_t g, const fmpz_mpoly_ctx_t ctx)

    Efficiently swap *f* and *g*.

.. function:: slong fmpz_mpoly_factor_length(const fmpz_mpoly_factor_t f, const fmpz_mpoly_ctx_t ctx)

    Return the length of the product in *f*.

.. function:: void fmpz_mpoly_factor_get_constant_fmpz(fmpz_t c, const fmpz_mpoly_factor_t f, const fmpz_mpoly_ctx_t ctx)
              void fmpz_mpoly_factor_get_constant_fmpq(fmpq_t c, const fmpz_mpoly_factor_t f, const fmpz_mpoly_ctx_t ctx)

    Set `c` to the constant of *f*.

.. function:: void fmpz_mpoly_factor_get_base(fmpz_mpoly_t B, const fmpz_mpoly_factor_t f, slong i, const fmpz_mpoly_ctx_t ctx)
              void fmpz_mpoly_factor_swap_base(fmpz_mpoly_t B, fmpz_mpoly_factor_t f, slong i, const fmpz_mpoly_ctx_t ctx)

    Set (resp. swap) *B* to (resp. with) the base of the term of index `i` in  *A*.

.. function:: slong fmpz_mpoly_factor_get_exp_si(fmpz_mpoly_factor_t f, slong i, const fmpz_mpoly_ctx_t ctx)

    Return the exponent of the term of index `i` in *A*. It is assumed to fit an ``slong``.

.. function:: void fmpz_mpoly_factor_sort(fmpz_mpoly_factor_t f, const fmpz_mpoly_ctx_t ctx)

    Sort the product of *f* first by exponent and then by base.


Factorisation
--------------------------------------------------------------------------------

    A return of `1` indicates that the function was successful. Otherwise,
    the return is `0` and *f* is undefined. None of these functions
    multiply *f* by *A*: *f* is simply set to a factorisation of *A*, and thus
    these functions should not depend on the initial value of the output *f*.

.. function:: int fmpz_mpoly_factor_squarefree(fmpz_mpoly_factor_t f, const fmpz_mpoly_t A, const fmpz_mpoly_ctx_t ctx)

    Set *f* to a factorization of *A* where the bases are primitive and
    pairwise relatively prime. If the product of all irreducible factors with
    a given exponent is desired, it is recommend to call :func:`fmpz_mpoly_factor_sort`
    and then multiply the bases with the desired exponent.

.. function:: int fmpz_mpoly_factor(fmpz_mpoly_factor_t f, const fmpz_mpoly_t A, const fmpz_mpoly_ctx_t ctx)

    Set *f* to a factorization of *A* where the bases are irreducible.

