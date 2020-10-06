.. _utils-flint:

**utils_flint.h** -- extra methods for Flint types
===============================================================================

General methods for multivariate polynomials
-------------------------------------------------------------------------------

.. function:: void fmpz_mpoly_primitive_part(fmpz_mpoly_t res, const fmpz_mpoly_t f, const fmpz_mpoly_ctx_t ctx)

    Sets *res* to the primitive part of *f*, obtained by dividing
    out the content of all coefficients and normalizing the leading
    coefficient to be positive. The zero polynomial is unchanged.

.. function:: void fmpz_mpoly_symmetric_gens(fmpz_mpoly_t res, ulong k, slong * vars, slong n, const fmpz_mpoly_ctx_t ctx)

.. function:: void fmpz_mpoly_symmetric(fmpz_mpoly_t res, ulong k, const fmpz_mpoly_ctx_t ctx)

    Sets *res* to the elementary symmetric polynomial
    `e_k(X_1,\ldots,X_n)`.

    The *gens* version takes `X_1,\ldots,X_n` to be the subset of
    generators given by *vars* and *n*.
    The indices in *vars* start from zero.
    Currently, the indices in *vars* must be distinct.

.. type:: fmpz_mpoly_vec_struct

.. type:: fmpz_mpoly_vec_t

    A type holding a vector of :type:`fmpz_mpoly_t`.

.. macro::  fmpz_mpoly_vec_entry(vec, i)

    Macro for accessing the entry at position *i* in *vec*.

.. function:: void fmpz_mpoly_vec_init(fmpz_mpoly_vec_t vec, const fmpz_mpoly_ctx_t ctx)

    Initializes *vec* for use, setting it to the empty vector.

.. function::void fmpz_mpoly_vec_clear(fmpz_mpoly_vec_t vec, const fmpz_mpoly_ctx_t ctx)

    Clears *vec*, freeing its allocated memory.

.. function:: void fmpz_mpoly_vec_print(const fmpz_mpoly_vec_t vec, const fmpz_mpoly_ctx_t ctx)

    Prints *vec* to standard output.

.. function:: void fmpz_mpoly_vec_swap(fmpz_mpoly_vec_t x, fmpz_mpoly_vec_t y, const fmpz_mpoly_ctx_t ctx)

    Swaps *x* and *y* efficiently.

.. function:: void fmpz_mpoly_vec_fit_length(fmpz_mpoly_vec_t vec, slong len, const fmpz_mpoly_ctx_t ctx)

    Allocates room for *len* entries in *vec*.

.. function:: void fmpz_mpoly_vec_set(fmpz_mpoly_vec_t dest, const fmpz_mpoly_vec_t src, const fmpz_mpoly_ctx_t ctx)

    Sets *dest* to a copy of *src*.

.. function:: void fmpz_mpoly_vec_append(fmpz_mpoly_vec_t vec, const fmpz_mpoly_t f, const fmpz_mpoly_ctx_t ctx)

    Appends *f* to the end of *vec*.

.. function:: slong fmpz_mpoly_vec_insert_unique(fmpz_mpoly_vec_t vec, const fmpz_mpoly_t f, const fmpz_mpoly_ctx_t ctx)

    Inserts *f* without duplication into *vec* and returns its index.
    If this polynomial already exists, *vec* is unchanged. If this
    polynomial does not exist in *vec*, it is appended.

.. function:: void fmpz_mpoly_vec_set_length(fmpz_mpoly_vec_t vec, slong len, const fmpz_mpoly_ctx_t ctx)

    Sets the length of *vec* to *len*, truncating or zero-extending
    as needed.

.. function:: void fmpz_mpoly_vec_randtest_not_zero(fmpz_mpoly_vec_t vec, flint_rand_t state, slong len, slong poly_len, slong bits, ulong exp_bound, fmpz_mpoly_ctx_t ctx)

    Sets *vec* to a random vector with exactly *len* entries, all nonzero,
    with random parameters defined by *poly_len*, *bits* and *exp_bound*.

.. function:: void fmpz_mpoly_vec_set_primitive_unique(fmpz_mpoly_vec_t res, const fmpz_mpoly_vec_t src, const fmpz_mpoly_ctx_t ctx)

    Sets *res* to a vector containing all polynomials in *src* reduced
    to their primitive parts, without duplication. The zero polynomial
    is skipped if present. The output order is arbitrary.

Ideals and Gröbner bases
-------------------------------------------------------------------------------

The following methods deal with ideals in `\mathbb{Q}[X_1,\ldots,X_n]`.
We use primitive integer polynomials as normalised generators
in place of monic rational polynomials.

.. function:: void fmpz_mpoly_spoly(fmpz_mpoly_t res, const fmpz_mpoly_t f, const fmpz_mpoly_t g, const fmpz_mpoly_ctx_t ctx)

    Sets *res* to the *S*-polynomial of *f* and *g*, scaled to
    an integer polynomial by computing the LCM of the leading coefficients.

.. function:: void fmpz_mpoly_reduction_primitive_part(fmpz_mpoly_t res, const fmpz_mpoly_t f, const fmpz_mpoly_vec_t vec, const fmpz_mpoly_ctx_t ctx)

    Sets *res* to the primitive part of the reduction (remainder of multivariate
    quasidivision with remainder) with respect to the polynomials *vec*.

.. function:: int fmpz_mpoly_vec_is_groebner(const fmpz_mpoly_vec_t G, const fmpz_mpoly_vec_t F, const fmpz_mpoly_ctx_t ctx)

    If *F* is *NULL*, checks if *G* is a Gröbner basis. If *F* is not *NULL*,
    checks if *G* is a Gröbner basis for *F*.

.. function:: int fmpz_mpoly_vec_is_autoreduced(const fmpz_mpoly_vec_t F, const fmpz_mpoly_ctx_t ctx);

    Checks whether the vector *F* is autoreduced (or inter-reduced).

.. function:: void fmpz_mpoly_vec_autoreduction(fmpz_mpoly_vec_t H, const fmpz_mpoly_vec_t F, const fmpz_mpoly_ctx_t ctx);

    Sets *H* to the autoreduction (inter-reduction) of *F*.

.. function:: void fmpz_mpoly_vec_autoreduction_groebner(fmpz_mpoly_vec_t H, const fmpz_mpoly_vec_t G, const fmpz_mpoly_ctx_t ctx)

    Sets *H* to the autoreduction (inter-reduction) of *G*.
    Assumes that *G* is a Gröbner basis.
    This produces a reduced Gröbner basis, which is unique
    (up to the sort order of the entries in the vector).

.. function:: pair_t fmpz_mpoly_select_pop_pair(pairs_t pairs, const fmpz_mpoly_vec_t G, const fmpz_mpoly_ctx_t ctx)

    Given a vector *pairs* of indices `(i, j)` into *G*, selects one pair
    for elimination in Buchberger's algorithm. The pair is removed
    from *pairs* and returned.

.. function:: void fmpz_mpoly_buchberger_naive(fmpz_mpoly_vec_t G, const fmpz_mpoly_vec_t F, const fmpz_mpoly_ctx_t ctx)

    Sets *G* to a Gröbner basis for *F*, computed using
    a naive implementation of Buchberger's algorithm.

.. function:: int fmpz_mpoly_buchberger_naive_with_limits(fmpz_mpoly_vec_t G, const fmpz_mpoly_vec_t F, slong ideal_len_limit, slong poly_len_limit, slong poly_bits_limit, const fmpz_mpoly_ctx_t ctx)

    As :func:`fmpz_mpoly_buchberger_naive`, but halts if during the
    execution of Buchberger's algorithm the length of the
    ideal basis set exceeds *ideal_len_limit*, the length of any
    polynomial exceeds *poly_len_limit*, or the size of the
    coefficients of any polynomial exceeds *poly_bits_limit*.
    Returns 1 for success and 0 for failure. On failure, *G* is
    a valid basis for *F* but it might not be a Gröbner basis.

Index pairs
-------------------------------------------------------------------------------

.. type:: pair_t

    A pair of *slong* indices *a* and *b*.

.. type:: pairs_struct

.. type:: pairs_t

    A type holding a vector of :type:`pair_t`.

.. function:: void pairs_init(pairs_t vec)

    Initializes *vec* for use, setting it to the empty vector of pairs.

.. function:: void pairs_fit_length(pairs_t vec, slong len)

    Allocates space for *len* elements in *vec*.

.. function:: void pairs_clear(pairs_t vec)

    Frees *vec*.

.. function:: void pairs_append(pairs_t vec, slong i, slong j)

    Appends the pair `(i, j)` to the end of *vec*.

.. function:: void pairs_insert_unique(pairs_t vec, slong i, slong j)

    Inserts `(i, j)` without duplication into *vec*. If this pair
    already exists, *vec* is unchanged. If this pair does not exist
    in *vec*, it is appended.

.. raw:: latex

    \newpage
