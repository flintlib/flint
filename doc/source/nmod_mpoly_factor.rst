.. _nmod-mpoly-factor:

**nmod_mpoly_factor.h** -- factorisation of multivariate polynomials over integers mod n (word-size n)
======================================================================================================

Types, macros and constants
-------------------------------------------------------------------------------

.. type:: nmod_mpoly_factor_struct

    A struct for holding a factored polynomial. There is a
    single constant and a product of bases to corresponding exponents.

.. type:: nmod_mpoly_factor_t

    An array of length `1` of ``nmod_mpoly_factor_struct``.


Memory management
--------------------------------------------------------------------------------


.. function:: void nmod_mpoly_factor_init(nmod_mpoly_factor_t f, const nmod_mpoly_ctx_t ctx)

    Initialise *f*.

.. function:: void nmod_mpoly_factor_clear(nmod_mpoly_factor_t f, const nmod_mpoly_ctx_t ctx)

    Clear *f*.


Basic manipulation
--------------------------------------------------------------------------------


.. function:: void nmod_mpoly_factor_swap(nmod_mpoly_factor_t f, nmod_mpoly_factor_t g, const nmod_mpoly_ctx_t ctx)

    Efficiently swap *f* and `*g*`.

.. function:: slong nmod_mpoly_factor_length(const nmod_mpoly_factor_t f, const nmod_mpoly_ctx_t ctx)

    Return the length of the product in *f*.

.. function:: void nmod_mpoly_factor_get_constant_ui(const nmod_mpoly_factor_t f, const nmod_mpoly_ctx_t ctx)

    Return the constant of *f*.

.. function:: void nmod_mpoly_factor_get_base(nmod_mpoly_t p, const nmod_mpoly_factor_t f, slong i, const nmod_mpoly_ctx_t ctx)
              void nmod_mpoly_factor_swap_base(nmod_mpoly_t p, nmod_mpoly_factor_t f, slong i, const nmod_mpoly_ctx_t ctx)

    Set (resp. swap) *B* to (resp. with) the base of the term of index `i` in  *A*.

.. function:: slong nmod_mpoly_factor_get_exp_si(nmod_mpoly_factor_t f, slong i, const nmod_mpoly_ctx_t ctx)

    Return the exponent of the term of index `i` in *A*. It is assumed to fit an ``slong``.

.. function:: void nmod_mpoly_factor_sort(nmod_mpoly_factor_t f, const nmod_mpoly_ctx_t ctx)

    Sort the product of *f* first by exponent and then by base.


Factorisation
--------------------------------------------------------------------------------

    A return of `1` indicates that the function was successful. Otherwise,
    the return is `0` and *f* is undefined. None of these functions
    multiply *f* by *A*: *f* is simply set to a factorisation of *A*, and thus
    these functions should not depend on the initial value of the output *f*.

.. function:: int nmod_mpoly_factor_squarefree(nmod_mpoly_factor_t f, const nmod_mpoly_t A, const nmod_mpoly_ctx_t ctx)

    Set *f* to a factorization of *A* where the bases are primitive and
    pairwise relatively prime. If the product of all irreducible factors with
    a given exponent is desired, it is recommend to call :func:`nmod_mpoly_factor_sort`
    and then multiply the bases with the desired exponent.

.. function:: int nmod_mpoly_factor(nmod_mpoly_factor_t f, const nmod_mpoly_t A, const nmod_mpoly_ctx_t ctx)

    Set *f* to a factorization of *A* where the bases are irreducible.

