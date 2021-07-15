.. _fq_nmod-mpoly-factor:

**fq_nmod_mpoly_factor.h** -- factorisation of multivariate polynomials over finite fields of word-sized characteristic
========================================================================================================================

Types, macros and constants
-------------------------------------------------------------------------------

.. type:: fq_nmod_mpoly_factor_struct

    A struct for holding a factored polynomial. There is a
    single constant and a product of bases to corresponding exponents.

.. type:: fq_nmod_mpoly_factor_t

    An array of length `1` of ``fq_nmod_mpoly_factor_struct``.


Memory management
--------------------------------------------------------------------------------


.. function:: void fq_nmod_mpoly_factor_init(fq_nmod_mpoly_factor_t f, const fq_nmod_mpoly_ctx_t ctx)

    Initialise *f*.

.. function:: void fq_nmod_mpoly_factor_clear(fq_nmod_mpoly_factor_t f, const fq_nmod_mpoly_ctx_t ctx)

    Clear *f*.


Basic manipulation
--------------------------------------------------------------------------------


.. function:: void fq_nmod_mpoly_factor_swap(fq_nmod_mpoly_factor_t f, fq_nmod_mpoly_factor_t g, const fq_nmod_mpoly_ctx_t ctx)

    Efficiently swap *f* and *g*.

.. function:: slong fq_nmod_mpoly_factor_length(const fq_nmod_mpoly_factor_t f, const fq_nmod_mpoly_ctx_t ctx)

    Return the length of the product in *f*.

.. function:: void fq_nmod_mpoly_factor_get_constant_fq_nmod(fq_nmod_t c, const fq_nmod_mpoly_factor_t f, const fq_nmod_mpoly_ctx_t ctx)

    Set `c` to the the constant of *f*.

.. function:: void fq_nmod_mpoly_factor_get_base(fq_nmod_mpoly_t p, const fq_nmod_mpoly_factor_t f, slong i, const fq_nmod_mpoly_ctx_t ctx)
              void fq_nmod_mpoly_factor_swap_base(fq_nmod_mpoly_t p, fq_nmod_mpoly_factor_t f, slong i, const fq_nmod_mpoly_ctx_t ctx)

    Set (resp. swap) *B* to (resp. with) the base of the term of index *i* in *A*.

.. function:: slong fq_nmod_mpoly_factor_get_exp_si(fq_nmod_mpoly_factor_t f, slong i, const fq_nmod_mpoly_ctx_t ctx)

    Return the exponent of the term of index *i* in *A*. It is assumed to fit an ``slong``.

.. function:: void fq_nmod_mpoly_factor_sort(fq_nmod_mpoly_factor_t f, const fq_nmod_mpoly_ctx_t ctx)

    Sort the product of *f* first by exponent and then by base.


Factorisation
--------------------------------------------------------------------------------

    A return of `1` indicates that the function was successful. Otherwise,
    the return is `0` and *f* is undefined. None of these functions
    multiply *f* by *A*: *f* is simply set to a factorisation of *A*, and thus
    these functions should not depend on the initial value of the output *f*.

.. function:: int fq_nmod_mpoly_factor_squarefree(fq_nmod_mpoly_factor_t f, const fq_nmod_mpoly_t A, const fq_nmod_mpoly_ctx_t ctx)

    Set *f* to a factorization of *A* where the bases are primitive and
    pairwise relatively prime. If the product of all irreducible factors with
    a given exponent is desired, it is recommend to call :func:`fq_nmod_mpoly_factor_sort`
    and then multiply the bases with the desired exponent.

.. function:: int fq_nmod_mpoly_factor(fq_nmod_mpoly_factor_t f, const fq_nmod_mpoly_t A, const fq_nmod_mpoly_ctx_t ctx)

    Set *f* to a factorization of *A* where the bases are irreducible.

