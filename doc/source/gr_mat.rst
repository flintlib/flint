.. _gr-mat:

**gr_mat.h** -- dense matrices over generic rings
===============================================================================

A :type:`gr_mat_t` represents a matrix implemented as a dense
array of entries in a generic ring *R*.

In this module, the context object ``ctx`` always represents the
coefficient ring *R* unless otherwise stated.
Creating a context object representing a matrix
space only becomes necessary when one
wants to manipulate matrices using generic ring methods
like ``gr_add`` instead of the designated matrix
methods like ``gr_mat_add``.

Warning: matrix functions generally assume that input as well
as output operands have compatible shapes.
Shape errors are not usually handled (this may change).

Type compatibility
-------------------------------------------------------------------------------

The ``gr_mat`` type has the same data layout as most
Flint, Arb and Calcium matrix types.
Methods in this module can therefore be mixed freely with
methods in the corresponding Flint, Arb and Calcium modules
when the underlying coefficient type is the same.

It is not directly compatible with the ``nmod_mat`` type
(modulus data is stored as part of the matrix object).

Types, macros and constants
-------------------------------------------------------------------------------

.. type:: gr_mat_struct

.. type:: gr_mat_t

    Contains a pointer to an array of coefficients (``entries``), the
    number of rows (``r``), the number of columns (``c``),
    and an array to pointers marking the start of each row (``rows``).

    A ``gr_mat_t`` is defined as an array of length one of type
    ``gr_mat_struct``, permitting a ``gr_mat_t`` to
    be passed by reference.

Basic operations
-------------------------------------------------------------------------------

.. macro:: GR_MAT_ENTRY(mat, i, j, sz)

    Macro to access the entry at row *i* and column *j* of the
    matrix *mat* whose entries have size *sz* bytes.

.. macro:: gr_mat_nrows(mat, ctx)

    Macro accessing the number of rows of *mat*.

.. macro:: gr_mat_ncols(mat, ctx)

    Macro accessing the number of columns of *mat*.

.. function:: void gr_mat_init(gr_mat_t mat, slong rows, slong cols, gr_ctx_t ctx)

    Initializes *mat* to a matrix with the given number of rows and
    columns.

.. function:: void gr_mat_clear(gr_mat_t mat, gr_ctx_t ctx)

    Clears the matrix.

.. function:: void gr_mat_window_init(gr_mat_t window, const gr_mat_t mat, slong r1, slong c1, slong r2, slong c2, gr_ctx_t ctx)

    Initializes *window* to a window matrix into the submatrix of *mat*
    starting at the corner at row *r1* and column *c1* (inclusive) and ending
    at row *r2* and column *c2* (exclusive).
    The indices must be within bounds.

.. function:: void gr_mat_window_clear(gr_mat_t window, gr_ctx_t ctx)

    Frees the window matrix.

.. function:: void gr_mat_swap(gr_mat_t mat1, gr_mat_t mat2, gr_ctx_t ctx)

    Swaps *mat1* and *mat12* efficiently.

.. function:: int gr_mat_swap_entrywise(gr_mat_t mat1, const gr_mat_t mat2, gr_ctx_t ctx)

    Performs a deep swap of *mat1* and *mat2*, swapping the individual
    entries rather than the top-level structures.

.. function:: int gr_mat_write(gr_stream_t out, const gr_mat_t mat, gr_ctx_t ctx)

    Write *mat* to the stream *out*.

.. function:: int gr_mat_print(const gr_mat_t mat, gr_ctx_t ctx)

    Prints *mat* to standard output.

.. function:: int gr_mat_randtest(gr_mat_t mat, flint_rand_t state, void * options, gr_ctx_t ctx)

    Sets *mat* to a random matrix.

.. function:: truth_t gr_mat_is_empty(const gr_mat_t mat, gr_ctx_t ctx)

    Returns whether *mat* is an empty matrix, having either zero
    rows or zero column. This predicate is always decidable (even if
    the underlying ring is not computable), returning
    ``T_TRUE`` or ``T_FALSE``.

.. function:: truth_t gr_mat_is_square(const gr_mat_t mat, gr_ctx_t ctx)

    Returns whether *mat* is a square matrix, having the same number
    of rows as columns (not the same thing as being a perfect square!).
    This predicate is always decidable (even if the underlying ring
    is not computable), returning ``T_TRUE`` or ``T_FALSE``.

.. function:: truth_t gr_mat_equal(const gr_mat_t mat1, const gr_mat_t mat2, gr_ctx_t ctx)

    Returns whether *mat1* and *mat2* are equal.

.. function:: truth_t gr_mat_is_zero(const gr_mat_t mat, gr_ctx_t ctx)
              truth_t gr_mat_is_one(const gr_mat_t mat, gr_ctx_t ctx)
              truth_t gr_mat_is_neg_one(const gr_mat_t mat, gr_ctx_t ctx)

    Returns whether *mat* respectively is the zero matrix or
    the scalar matrix with 1 or -1 on the main diagonal.

.. function:: int gr_mat_zero(gr_mat_t res, gr_ctx_t ctx)

    Sets *res* to the zero matrix.

.. function:: int gr_mat_one(gr_mat_t res, gr_ctx_t ctx)

    Sets *res* to the scalar matrix with 1 on the main diagonal
    and zero elsewhere.

.. function:: int gr_mat_set(gr_mat_t res, const gr_mat_t mat, gr_ctx_t ctx)

    Sets *res* to the value of *mat*.

.. function:: int gr_mat_set_scalar(gr_mat_t res, gr_srcptr c, gr_ctx_t ctx)
              int gr_mat_set_ui(gr_mat_t res, ulong c, gr_ctx_t ctx)
              int gr_mat_set_si(gr_mat_t res, slong c, gr_ctx_t ctx)
              int gr_mat_set_fmpz(gr_mat_t res, const fmpz_t c, gr_ctx_t ctx)
              int gr_mat_set_fmpq(gr_mat_t res, const fmpq_t c, gr_ctx_t ctx)

    Set *res* to the scalar matrix with *c* on the main diagonal
    and zero elsewhere.

.. function:: int gr_mat_transpose(gr_mat_t B, const gr_mat_t A, gr_ctx_t ctx)

Arithmetic
-------------------------------------------------------------------------------

.. function:: int gr_mat_neg(gr_mat_t res, const gr_mat_t mat, gr_ctx_t ctx)

.. function:: int gr_mat_add(gr_mat_t res, const gr_mat_t mat1, const gr_mat_t mat2, gr_ctx_t ctx)

.. function:: int gr_mat_sub(gr_mat_t res, const gr_mat_t mat1, const gr_mat_t mat2, gr_ctx_t ctx)

.. function:: int gr_mat_mul_classical(gr_mat_t res, const gr_mat_t mat1, const gr_mat_t mat2, gr_ctx_t ctx)
              int gr_mat_mul(gr_mat_t res, const gr_mat_t mat1, const gr_mat_t mat2, gr_ctx_t ctx)

.. function:: int gr_mat_sqr(gr_mat_t res, const gr_mat_t mat, gr_ctx_t ctx)

.. function:: int gr_mat_add_scalar(gr_mat_t res, const gr_mat_t mat, gr_srcptr c, gr_ctx_t ctx)
              int gr_mat_sub_scalar(gr_mat_t res, const gr_mat_t mat, gr_srcptr c, gr_ctx_t ctx)
              int gr_mat_mul_scalar(gr_mat_t res, const gr_mat_t mat, gr_srcptr c, gr_ctx_t ctx)
              int gr_mat_addmul_scalar(gr_mat_t res, const gr_mat_t mat, gr_srcptr c, gr_ctx_t ctx)
              int gr_mat_submul_scalar(gr_mat_t res, const gr_mat_t mat, gr_srcptr c, gr_ctx_t ctx)

LU decomposition
-------------------------------------------------------------------------------

.. function:: int gr_mat_lu_classical(slong * rank, slong * P, gr_mat_t LU, const gr_mat_t A, int full_rank_check, gr_ctx_t ctx)

Determinant and trace
-------------------------------------------------------------------------------

.. function:: int gr_mat_trace(gr_ptr res, const gr_mat_t mat, gr_ctx_t ctx)

Solving
-------------------------------------------------------------------------------

.. function:: int gr_mat_solve_tril_classical(gr_mat_t X, const gr_mat_t L, const gr_mat_t B, int unit, gr_ctx_t ctx)
              int gr_mat_solve_tril_recursive(gr_mat_t X, const gr_mat_t L, const gr_mat_t B, int unit, gr_ctx_t ctx)
              int gr_mat_solve_tril(gr_mat_t X, const gr_mat_t L, const gr_mat_t B, int unit, gr_ctx_t ctx)
              int gr_mat_solve_triu_classical(gr_mat_t X, const gr_mat_t U, const gr_mat_t B, int unit, gr_ctx_t ctx)
              int gr_mat_solve_triu_recursive(gr_mat_t X, const gr_mat_t U, const gr_mat_t B, int unit, gr_ctx_t ctx)
              int gr_mat_solve_triu(gr_mat_t X, const gr_mat_t U, const gr_mat_t B, int unit, gr_ctx_t ctx)

    Solves the lower triangular system `LX = B` or the upper triangular system
    `UX = B`, respectively. Division by the the diagonal entries must
    be possible; if not a division fails, ``GR_DOMAIN`` is returned
    even if the system is solvable.
    If *unit* is set, the main diagonal of *L* or *U*
    is taken to consist of all ones, and in that case the actual entries on
    the diagonal are not read at all and can contain other data.

    The *classical* versions perform the computations iteratively while the
    *recursive* versions perform the computations in a block recursive
    way to benefit from fast matrix multiplication. The default versions
    choose an algorithm automatically.


Characteristic polynomial
-------------------------------------------------------------------------------

.. function:: int _gr_mat_charpoly_berkowitz(gr_ptr res, const gr_mat_t mat, gr_ctx_t ctx)
              int gr_mat_charpoly_berkowitz(gr_poly_t res, const gr_mat_t mat, gr_ctx_t ctx)

    Sets *res* to the characteristic polynomial of the square matrix
    *mat*, computed using the division-free Berkowitz algorithm.
    The number of operations is `O(n^4)` where *n* is the
    size of the matrix. The
    underscore method assumes that *res* is a preallocated
    array of `n + 1` coefficients.

.. function:: int _gr_mat_charpoly_danilevsky_inplace(gr_ptr res, gr_mat_t mat, gr_ctx_t ctx)
              int _gr_mat_charpoly_danilevsky(gr_ptr res, const gr_mat_t mat, gr_ctx_t ctx)
              int gr_mat_charpoly_danilevsky(gr_poly_t res, const gr_mat_t mat, gr_ctx_t ctx)

    Sets *res* to the characteristic polynomial of the square matrix
    *mat*, computed using the Danilevsky algorithm.
    The number of operations is `O(n^3)` where *n* is the
    size of the matrix. The
    underscore method assumes that *res* is a preallocated
    array of `n + 1` coefficients.
    The *inplace* version overwrites the input matrix.

    This method requires divisions and can therefore fail when the
    ring is not a field, but will sometimes succeed anywyay. It
    also requires testing for zero. It returns
    the ``GR_UNABLE`` or ``GR_DOMAIN`` flag when an impossible division
    is encountered or when a comparison cannot be performed.

.. function:: int _gr_mat_charpoly_hessenberg(gr_ptr res, const gr_mat_t mat, gr_ctx_t ctx)
              int gr_mat_charpoly_hessenberg(gr_poly_t res, const gr_mat_t mat, gr_ctx_t ctx)

    Sets *res* to the characteristic polynomial of the square matrix
    *mat*, which is assumed to be in Hessenberg form (this is
    currently not checked).

.. function:: int _gr_mat_charpoly_faddeev(gr_ptr res, gr_mat_t adj, const gr_mat_t mat, gr_ctx_t ctx)
              int gr_mat_charpoly_faddeev(gr_poly_t res, gr_mat_t adj, const gr_mat_t mat, gr_ctx_t ctx)
              int _gr_mat_charpoly_faddeev_bsgs(gr_ptr res, gr_mat_t adj, const gr_mat_t mat, gr_ctx_t ctx)
              int gr_mat_charpoly_faddeev_bsgs(gr_poly_t res, gr_mat_t adj, const gr_mat_t mat, gr_ctx_t ctx)

    Sets *res* to the characteristic polynomial of the square matrix
    *mat*, computed using the Faddeev-LeVerrier algorithm.
    If the optional output argument *adj* is not *NULL*, it is
    set to the adjugate matrix, which is computed free of charge.

    The *bsgs* version uses a baby-step giant-step strategy,
    also known as the Preparata-Sarwate algorithm.
    This reduces the complexity from `O(n^4)` to `O(n^{3.5})` operations
    at the cost of requiring `n^{0.5}` temporary matrices to be
    stored.

    This method requires divisions by small integers and can
    therefore fail (returning the ``GR_UNABLE`` or ``GR_DOMAIN`` flags)
    in finite characteristic or when the underlying ring does
    not implement a division algorithm.

Hessenberg form
-------------------------------------------------------------------------------

.. function:: truth_t gr_mat_is_hessenberg(const gr_mat_t mat, gr_ctx_t ctx)

    Returns whether *mat* is in upper Hessenberg form.

.. function:: int gr_mat_hessenberg_gauss(gr_mat_t res, const gr_mat_t mat, gr_ctx_t ctx)
              int gr_mat_hessenberg_householder(gr_mat_t res, const gr_mat_t mat, gr_ctx_t ctx)
              int gr_mat_hessenberg(gr_mat_t res, const gr_mat_t mat, gr_ctx_t ctx)

    Sets *res* to an upper Hessenberg form of *mat*.
    The *gauss* version uses Gaussian elimination.
    The *householder* version uses Householder reflections.

    These methods require divisions and zero testing
    and can therefore fail (returning ``GR_UNABLE`` or ``GR_DOMAIN``)
    when the ring is not a field.
    The *householder* version additionally requires complex
    conjugation and the ability to compute square roots.

Special matrices
-------------------------------------------------------------------------------

.. function:: int gr_mat_ones(gr_mat_t res, gr_ctx_t ctx)

    Sets all entries in *res* to one.

.. function:: int gr_mat_pascal(gr_mat_t res, int triangular, gr_ctx_t ctx)

    Sets *res* to a Pascal matrix, whose entries are binomial coefficients.
    If *triangular* is 0, constructs a full symmetric matrix
    with the rows of Pascal's triangle as successive antidiagonals.
    If *triangular* is 1, constructs the upper triangular matrix with
    the rows of Pascal's triangle as columns, and if *triangular* is -1,
    constructs the lower triangular matrix with the rows of Pascal's
    triangle as rows.

.. function:: int gr_mat_stirling(gr_mat_t res, int kind, gr_ctx_t ctx)

    Sets *res* to a Stirling matrix, whose entries are Stirling numbers.
    If *kind* is 0, the entries are set to the unsigned Stirling numbers
    of the first kind. If *kind* is 1, the entries are set to the signed
    Stirling numbers of the first kind. If *kind* is 2, the entries are
    set to the Stirling numbers of the second kind.

.. function:: int gr_mat_hilbert(gr_mat_t mat, gr_ctx_t ctx)

    Sets *res* to the Hilbert matrix, which has entries `1/(i+j+1)`
    for `i, j \ge 0`.
