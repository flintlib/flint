.. _fmpz-mod-mpoly-factor:

**fmpz_mod_mpoly_factor.h** -- factorisation of multivariate polynomials over the integers mod n
================================================================================================

Types, macros and constants
-------------------------------------------------------------------------------

.. type:: fmpz_mod_mpoly_factor_struct

    A struct for holding a factored polynomial over the integers mod n. There is a
    single constant and a product of bases to corresponding exponents.

.. type:: fmpz_mod_mpoly_factor_t

    An array of length `1` of ``fmpz_mod_mpoly_factor_struct``.


Memory management
--------------------------------------------------------------------------------


.. function:: void fmpz_mod_mpoly_factor_init(fmpz_mod_mpoly_factor_t f, const fmpz_mod_mpoly_ctx_t ctx)

    Initialise *f*.

.. function:: void fmpz_mod_mpoly_factor_clear(fmpz_mod_mpoly_factor_t f, const fmpz_mod_mpoly_ctx_t ctx)

    Clear *f*.


Basic manipulation
--------------------------------------------------------------------------------


.. function:: void fmpz_mod_mpoly_factor_swap(fmpz_mod_mpoly_factor_t f, fmpz_mod_mpoly_factor_t g, const fmpz_mod_mpoly_ctx_t ctx)

    Efficiently swap *f* and *g*.

.. function:: slong fmpz_mod_mpoly_factor_length(const fmpz_mod_mpoly_factor_t f, const fmpz_mod_mpoly_ctx_t ctx)

    Return the length of the product in *f*.

.. function:: void fmpz_mod_mpoly_factor_get_constant_fmpz(fmpz_t c, const fmpz_mod_mpoly_factor_t f, const fmpz_mod_mpoly_ctx_t ctx)

    Set *c* to the constant of *f*.

.. function:: void fmpz_mod_mpoly_factor_get_base(fmpz_mod_mpoly_t B, const fmpz_mod_mpoly_factor_t f, slong i, const fmpz_mod_mpoly_ctx_t ctx)
              void fmpz_mod_mpoly_factor_swap_base(fmpz_mod_mpoly_t B, fmpz_mod_mpoly_factor_t f, slong i, const fmpz_mod_mpoly_ctx_t ctx)

    Set (resp. swap) *B* to (resp. with) the base of the term of index *i* in  *f*.

.. function:: slong fmpz_mod_mpoly_factor_get_exp_si(fmpz_mod_mpoly_factor_t f, slong i, const fmpz_mod_mpoly_ctx_t ctx)

    Return the exponent of the term of index *i* in *f*. It is assumed to fit an ``slong``.

.. function:: void fmpz_mod_mpoly_factor_sort(fmpz_mod_mpoly_factor_t f, const fmpz_mod_mpoly_ctx_t ctx)

    Sort the product of *f* first by exponent and then by base.


Factorisation
--------------------------------------------------------------------------------

    A return of `1` indicates that the function was successful. Otherwise,
    the return is `0` and *f* is undefined. None of these functions
    multiply *f* by *A*: *f* is simply set to a factorisation of *A*, and thus
    these functions should not depend on the initial value of the output *f*.

.. function:: int fmpz_mod_mpoly_factor_squarefree(fmpz_mod_mpoly_factor_t f, const fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_ctx_t ctx)

    Set *f* to a factorization of *A* where the bases are primitive and
    pairwise relatively prime. If the product of all irreducible factors with
    a given exponent is desired, it is recommend to call :func:`fmpz_mod_mpoly_factor_sort`
    and then multiply the bases with the desired exponent.

.. function:: int fmpz_mod_mpoly_factor(fmpz_mod_mpoly_factor_t f, const fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_ctx_t ctx)

    Set *f* to a factorization of *A* where the bases are irreducible.

