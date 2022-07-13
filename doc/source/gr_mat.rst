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

Warnings
-------------------------------------------------------------------------------

* Matrix functions generally assume that input as well
  as output operands have compatible shapes.
  Shape errors are not usually handled (this may change).
* Some operations (like rank, LU factorization) generally only make
  sense when the base ring is an integral domain.
  Typically the algorithms designed for integral domains also work
  over non-integral domains as long as all inversions of nonzero
  elements succeed. If an inversion fails, the algorithm will return
  the ``GR_DOMAIN`` or ``GR_UNABLE`` flag.
  This might not yet be entirely consistent.

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
              int gr_mat_set_fmpz_mat(gr_mat_t res, const fmpz_mat_t mat, gr_ctx_t ctx)
              int gr_mat_set_fmpq_mat(gr_mat_t res, const fmpq_mat_t mat, gr_ctx_t ctx)

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
              int gr_mat_mul_generic(gr_mat_t C, const gr_mat_t A, const gr_mat_t B, gr_ctx_t ctx)
              int gr_mat_mul(gr_mat_t res, const gr_mat_t mat1, const gr_mat_t mat2, gr_ctx_t ctx)

    Matrix multiplication. The default function can be overloaded by specific rings;
    otherwise, it falls back to :func:`gr_mat_mul_generic` which currently
    only performs classical multiplication.

.. function:: int gr_mat_sqr(gr_mat_t res, const gr_mat_t mat, gr_ctx_t ctx)

.. function:: int gr_mat_add_scalar(gr_mat_t res, const gr_mat_t mat, gr_srcptr c, gr_ctx_t ctx)
              int gr_mat_sub_scalar(gr_mat_t res, const gr_mat_t mat, gr_srcptr c, gr_ctx_t ctx)
              int gr_mat_mul_scalar(gr_mat_t res, const gr_mat_t mat, gr_srcptr c, gr_ctx_t ctx)
              int gr_mat_addmul_scalar(gr_mat_t res, const gr_mat_t mat, gr_srcptr c, gr_ctx_t ctx)
              int gr_mat_submul_scalar(gr_mat_t res, const gr_mat_t mat, gr_srcptr c, gr_ctx_t ctx)
              int gr_mat_div_scalar(gr_mat_t res, const gr_mat_t mat, gr_srcptr c, gr_ctx_t ctx)

.. function:: int _gr_mat_gr_poly_evaluate(gr_mat_t res, gr_srcptr poly, slong len, const gr_mat_t mat, gr_ctx_t ctx)
              int gr_mat_gr_poly_evaluate(gr_mat_t res, const gr_poly_t poly, const gr_mat_t mat, gr_ctx_t ctx)

    Sets *res* to the matrix obtained by evaluating the
    scalar polynomial *poly* with matrix argument *mat*.

Gaussian elimination
-------------------------------------------------------------------------------

.. function:: int gr_mat_find_nonzero_pivot(slong * pivot_row, gr_mat_t mat, slong start_row, slong end_row, slong column, gr_ctx_t ctx)

    Attempts to find a nonzero element in column number *column*
    of the matrix *mat* in a row between *start_row* (inclusive)
    and *end_row* (exclusive).
    On success, sets ``pivot_row`` to the row index and returns
    ``GR_SUCCESS``. If no nonzero pivot element exists, returns ``GR_DOMAIN``.
    If no nonzero pivot element exists and zero-testing fails for some
    element, returns the flag ``GR_UNABLE``.

    This function may be destructive: any elements that are nontrivially
    zero but can be certified zero may be overwritten by exact zeros.

.. function:: int gr_mat_lu_classical(slong * rank, slong * P, gr_mat_t LU, const gr_mat_t A, int rank_check, gr_ctx_t ctx)
              int gr_mat_lu_recursive(slong * rank, slong * P, gr_mat_t LU, const gr_mat_t A, int rank_check, gr_ctx_t ctx)
              int gr_mat_lu(slong * rank, slong * P, gr_mat_t LU, const gr_mat_t A, int rank_check, gr_ctx_t ctx)

    Computes a generalized LU decomposition `A = PLU` of a given
    matrix *A*, writing the rank of *A* to *rank*.

    If *A* is a nonsingular square matrix, *LU* will be set to
    a unit diagonal lower triangular matrix *L* and an upper
    triangular matrix *U* (the diagonal of *L* will not be stored
    explicitly).

    If *A* is an arbitrary matrix of rank *r*, *U* will be in row
    echelon form having *r* nonzero rows, and *L* will be lower
    triangular but truncated to *r* columns, having implicit ones on
    the *r* first entries of the main diagonal. All other entries will
    be zero.

    If a nonzero value for ``rank_check`` is passed, the function
    will abandon the output matrix in an undefined state and set
    the rank to 0 if *A* is detected to be rank-deficient.
    This currently only works as expected for square matrices.

    The algorithm can fail if it fails to certify that a pivot
    element is zero or nonzero, in which case the correct rank
    cannot be determined. It can also fail if a pivot element
    is not invertible. In these cases the ``GR_UNABLE`` and/or
    ``GR_DOMAIN`` flags will be returned. On failure,
    the data in the output variables
    ``rank``, ``P`` and ``LU`` will be meaningless.

    The *classical* version uses iterative Gaussian elimination.
    The *recursive* version uses a block recursive algorithm
    to take advantage of fast matrix multiplication.

.. function:: int gr_mat_fflu(slong * rank, slong * P, gr_mat_t LU, gr_ptr den, const gr_mat_t A, int rank_check, gr_ctx_t ctx)

    Similar to :func:`gr_mat_lu`, but computes a fraction-free
    LU decomposition using the Bareiss algorithm.
    The denominator is written to *den*.

Solving
-------------------------------------------------------------------------------

.. function:: int gr_mat_nonsingular_solve_tril_classical(gr_mat_t X, const gr_mat_t L, const gr_mat_t B, int unit, gr_ctx_t ctx)
              int gr_mat_nonsingular_solve_tril_recursive(gr_mat_t X, const gr_mat_t L, const gr_mat_t B, int unit, gr_ctx_t ctx)
              int gr_mat_nonsingular_solve_tril(gr_mat_t X, const gr_mat_t L, const gr_mat_t B, int unit, gr_ctx_t ctx)
              int gr_mat_nonsingular_solve_triu_classical(gr_mat_t X, const gr_mat_t U, const gr_mat_t B, int unit, gr_ctx_t ctx)
              int gr_mat_nonsingular_solve_triu_recursive(gr_mat_t X, const gr_mat_t U, const gr_mat_t B, int unit, gr_ctx_t ctx)
              int gr_mat_nonsingular_solve_triu(gr_mat_t X, const gr_mat_t U, const gr_mat_t B, int unit, gr_ctx_t ctx)

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

.. function:: int gr_mat_nonsingular_solve_fflu(gr_mat_t X, const gr_mat_t A, const gr_mat_t B, gr_ctx_t ctx)
              int gr_mat_nonsingular_solve_lu(gr_mat_t X, const gr_mat_t A, const gr_mat_t B, gr_ctx_t ctx)
              int gr_mat_nonsingular_solve(gr_mat_t X, const gr_mat_t A, const gr_mat_t B, gr_ctx_t ctx)

    Solves `AX = B`. If *A* is not invertible,
    returns ``GR_DOMAIN`` even if the system has a solution.

.. function:: int gr_mat_nonsingular_solve_fflu_precomp(gr_mat_t X, const slong * perm, const gr_mat_t LU, const gr_mat_t B, gr_ctx_t ctx)
              int gr_mat_nonsingular_solve_lu_precomp(gr_mat_t X, const slong * perm, const gr_mat_t LU, const gr_mat_t B, gr_ctx_t ctx)

    Solves `AX = B` given a precomputed FFLU or LU factorization of *A*.

.. function:: int gr_mat_nonsingular_solve_den_fflu(gr_mat_t X, gr_ptr den, const gr_mat_t A, const gr_mat_t B, gr_ctx_t ctx)
              int gr_mat_nonsingular_solve_den(gr_mat_t X, gr_ptr den, const gr_mat_t A, const gr_mat_t B, gr_ctx_t ctx)

    Solves `AX = B` over the fraction field of the present ring
    (assumed to be an integral domain), returning `X` with
    an implied denominator *den*.
    If *A* is not invertible over the fraction field, returns
    ``GR_DOMAIN`` even if the system has a solution.

Determinant and trace
-------------------------------------------------------------------------------

.. function:: int gr_mat_det_bareiss(gr_ptr res, const gr_mat_t mat, gr_ctx_t ctx)
              int gr_mat_det_berkowitz(gr_ptr res, const gr_mat_t mat, gr_ctx_t ctx)
              int gr_mat_det_lu(gr_ptr res, const gr_mat_t mat, gr_ctx_t ctx)
              int gr_mat_det_cofactor(gr_ptr res, const gr_mat_t mat, gr_ctx_t ctx)
              int gr_mat_det(gr_ptr res, const gr_mat_t mat, gr_ctx_t ctx)

    Sets *res* to the determinant of the square matrix *mat*.
    Various algorithms are available:

    * The *berkowitz* version uses the division-free Berkowitz algorithm
      performing `O(n^4)` operations. Since no zero tests are required, it
      is guaranteed to succeed if the ring arithmetic succeeds.

    * The *cofactor* version performs cofactor expansion. This is currently
      only supported for matrices up to size 4, and for larger
      matrices returns the ``GR_UNABLE`` flag.

    * The *lu* and *bareiss* versions use rational LU decomposition
      and fraction-free LU decomposition (Bareiss algorithm) respectively,
      requiring `O(n^3)` operations. These algorithms can fail if zero
      certification or inversion fails, in which case the ``GR_UNABLE``
      flag is returned.

    If the matrix is not square, ``GR_DOMAIN`` is returned.

.. function:: int gr_mat_trace(gr_ptr res, const gr_mat_t mat, gr_ctx_t ctx)

    Sets *res* to the trace (sum of entries on the main diagonal) of
    the square matrix *mat*.
    If the matrix is not square, ``GR_DOMAIN`` is returned.


Rank
-------------------------------------------------------------------------------

.. function:: int gr_mat_rank_fflu(slong * rank, const gr_mat_t mat, gr_ctx_t ctx)
              int gr_mat_rank_lu(slong * rank, const gr_mat_t mat, gr_ctx_t ctx)
              int gr_mat_rank(slong * rank, const gr_mat_t mat, gr_ctx_t ctx)

    Sets *res* to the rank of *mat*.
    The default method returns ``GR_DOMAIN`` if the element ring
    is not an integral domain, in which case the usual rank is
    not well-defined. The *fflu* and *lu* variants currently do
    not check the element domain, and simply return this flag if they
    encounter an impossible inverse in the execution of the
    respective algorithms.

Row echelon form
-------------------------------------------------------------------------------

.. function:: int gr_mat_rref_lu(slong * rank, gr_mat_t R, const gr_mat_t A, gr_ctx_t ctx)
              int gr_mat_rref_fflu(slong * rank, gr_mat_t R, const gr_mat_t A, gr_ctx_t ctx)
              int gr_mat_rref(slong * rank, gr_mat_t R, const gr_mat_t A, gr_ctx_t ctx)

    Sets *R* to the reduced row echelon form of *A*, also setting
    *rank* to its rank.

.. function:: int gr_mat_rref_den_fflu(slong * rank, gr_mat_t R, gr_ptr den, const gr_mat_t A, gr_ctx_t ctx)
              int gr_mat_rref_den(slong * rank, gr_mat_t R, gr_ptr den, const gr_mat_t A, gr_ctx_t ctx)

    Like *rref*, but computes the reduced row echelon multiplied
    by a common (not necessarily minimal) denominator which is written
    to *den*. This can be used to compute the rref over an integral
    domain which is not a field.

Nullspace
-------------------------------------------------------------------------------

.. function:: int gr_mat_nullspace(gr_mat_t X, const gr_mat_t A, gr_ctx_t ctx)

    Sets *X* to a basis for the (right) nullspace of *A*.
    On success, the output matrix will be resized to the correct
    number of columns.

    The basis is not guaranteed to be presented in a
    canonical or minimal form.

    If the ring is not a field, this is implied to compute a nullspace
    basis over the fraction field. The result may be meaningless
    if the ring is not an integral domain.

Inverse and adjugate
-------------------------------------------------------------------------------

.. function:: int gr_mat_inv(gr_mat_t res, const gr_mat_t mat, gr_ctx_t ctx)

    Sets *res* to the inverse of *mat*, computed by solving
    `A A^{-1} = I`.

    Returns ``GR_DOMAIN`` if it can be determined that *mat* is not
    invertible over the present ring (warning: this may not work
    over non-integral domains). If invertibility cannot be proved,
    returns ``GR_UNABLE``.

    To compute the inverse over the fraction field, one may use
    :func:`gr_mat_nonsingular_solve_den` or :func:`gr_mat_adjugate`.

.. function:: int gr_mat_adjugate_charpoly(gr_mat_t adj, gr_ptr det, const gr_mat_t mat, gr_ctx_t ctx)
              int gr_mat_adjugate_cofactor(gr_mat_t adj, gr_ptr det, const gr_mat_t mat, gr_ctx_t ctx)
              int gr_mat_adjugate(gr_mat_t adj, gr_ptr det, const gr_mat_t mat, gr_ctx_t ctx)

    Sets *adj* to the adjugate matrix of *mat*, simultaneously
    setting *det* to the determinant of *mat*. We have
    `\operatorname{adj}(A) A = A \operatorname{adj}(A) = \det(A) I`,
    and `A^{-1} = \operatorname{adj}(A) / \det(A)` when *A*
    is invertible.

    The *cofactor* version uses cofactor expansion, requiring the
    evaluation of `n^2` determinants.
    The *charpoly* version computes and then evaluates the
    characteristic polynomial, requiring `O(n^{1/2})`
    matrix multiplications plus `O(n^3)` or `O(n^4)` operations
    for the characteristic polynomial itself depending on the
    algorithm used.

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

Random matrices
-------------------------------------------------------------------------------

.. function:: int gr_mat_randtest(gr_mat_t res, flint_rand_t state, gr_ctx_t ctx)

    Sets *res* to a random matrix. The distribution is nonuniform.

.. function:: int gr_mat_randops(gr_mat_t mat, flint_rand_t state, slong count, gr_ctx_t ctx)

    Randomises *mat* in-place by performing elementary row or column
    operations. More precisely, at most *count* random additions or
    subtractions of distinct rows and columns will be performed.

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

.. function:: int gr_mat_hilbert(gr_mat_t res, gr_ctx_t ctx)

    Sets *res* to the Hilbert matrix, which has entries `1/(i+j+1)`
    for `i, j \ge 0`.

.. function:: int gr_mat_hadamard(gr_mat_t res, gr_ctx_t ctx)

    If possible, sets *res* to a Hadamard matrix of the provided size
    and returns ``GR_SUCCESS``. Returns ``GR_DOMAIN``
    if no Hadamard matrix of the given size exists,
    and ``GR_UNABLE`` if the implementation does
    not know how to construct a Hadamard matrix of the given
    size.

    A Hadamard matrix of size *n* can only exist if *n* is 0, 1, 2,
    or a multiple of 4. It is not known whether a
    Hadamard matrix exists for every size that is a multiple of 4.
    This function uses the Paley construction, which
    succeeds for all *n* of the form `n = 2^e` or `n = 2^e (q + 1)` where
    *q* is an odd prime power. Orders *n* for which Hadamard matrices are
    known to exist but for which this construction fails are
    92, 116, 156, ... (OEIS A046116).

