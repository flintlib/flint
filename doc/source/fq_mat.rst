.. _fq-mat:

**fq_mat.h** -- matrices over finite fields
===============================================================================

Types, macros and constants
-------------------------------------------------------------------------------

.. type:: fq_mat_struct

.. type:: fq_mat_t

Memory management
--------------------------------------------------------------------------------


.. function:: void fq_mat_init(fq_mat_t mat, slong rows, slong cols, const fq_ctx_t ctx)

    Initialises ``mat`` to a ``rows``-by-``cols`` matrix with
    coefficients in `\mathbf{F}_{q}` given by ``ctx``. All elements
    are set to zero.

.. function:: void fq_mat_init_set(fq_mat_t mat, const fq_mat_t src, const fq_ctx_t ctx)

    Initialises ``mat`` and sets its dimensions and elements to
    those of ``src``.

.. function:: void fq_mat_clear(fq_mat_t mat, const fq_ctx_t ctx)

    Clears the matrix and releases any memory it used. The matrix
    cannot be used again until it is initialised. This function must be
    called exactly once when finished using an ``fq_mat_t`` object.

.. function:: void fq_mat_set(fq_mat_t mat, const fq_mat_t src, const fq_ctx_t ctx)

    Sets ``mat`` to a copy of ``src``. It is assumed
    that ``mat`` and ``src`` have identical dimensions.


Basic properties and manipulation
--------------------------------------------------------------------------------


.. function:: fq_struct * fq_mat_entry(const fq_mat_t mat, slong i, slong j)

    Directly accesses the entry in ``mat`` in row `i` and column `j`,
    indexed from zero. No bounds checking is performed.

.. function:: void fq_mat_entry_set(fq_mat_t mat, slong i, slong j, const fq_t x, const fq_ctx_t ctx)

    Sets the entry in ``mat`` in row `i` and column `j` to ``x``.

.. function:: slong fq_mat_nrows(const fq_mat_t mat, const fq_ctx_t ctx)

    Returns the number of rows in ``mat``.

.. function:: slong fq_mat_ncols(const fq_mat_t mat, const fq_ctx_t ctx)

    Returns the number of columns in ``mat``.

.. function:: void fq_mat_swap(fq_mat_t mat1, fq_mat_t mat2, const fq_ctx_t ctx)

    Swaps two matrices. The dimensions of ``mat1`` and ``mat2``
    are allowed to be different.

.. function:: void fq_mat_swap_entrywise(fq_mat_t mat1, fq_mat_t mat2, const fq_ctx_t ctx)

    Swaps two matrices by swapping the individual entries rather than swapping
    the contents of the structs.

.. function:: void fq_mat_zero(fq_mat_t mat, const fq_ctx_t ctx)

    Sets all entries of ``mat`` to 0.

.. function:: void fq_mat_one(fq_mat_t mat, const fq_ctx_t ctx)

    Sets all the diagonal entries of ``mat`` to 1 and all other entries to 0.

.. function:: void fq_mat_swap_rows(fq_mat_t mat, slong * perm, slong r, slong s, const fq_ctx_t ctx)

    Swaps rows ``r`` and ``s`` of ``mat``.  If ``perm`` is non-``NULL``, the
    permutation of the rows will also be applied to ``perm``.

.. function:: void fq_mat_swap_cols(fq_mat_t mat, slong * perm, slong r, slong s, const fq_ctx_t ctx)

    Swaps columns ``r`` and ``s`` of ``mat``.  If ``perm`` is non-``NULL``, the
    permutation of the columns will also be applied to ``perm``.

.. function:: void fq_mat_invert_rows(fq_mat_t mat, slong * perm, const fq_ctx_t ctx)

    Swaps rows ``i`` and ``r - i`` of ``mat`` for ``0 <= i < r/2``, where
    ``r`` is the number of rows of ``mat``. If ``perm`` is non-``NULL``, the
    permutation of the rows will also be applied to ``perm``.

.. function:: void fq_mat_invert_cols(fq_mat_t mat, slong * perm, const fq_ctx_t ctx)

    Swaps columns ``i`` and ``c - i`` of ``mat`` for ``0 <= i < c/2``, where
    ``c`` is the number of columns of ``mat``. If ``perm`` is non-``NULL``, the
    permutation of the columns will also be applied to ``perm``.


Conversions
--------------------------------------------------------------------------------

.. function:: void fq_mat_set_nmod_mat(fq_mat_t mat1, const nmod_mat_t mat2, const fq_ctx_t ctx)

    Sets the matrix ``mat1`` to the matrix ``mat2``.

.. function:: void fq_mat_set_fmpz_mod_mat(fq_mat_t mat1, const fmpz_mod_mat_t mat2, const fq_ctx_t ctx)

    Sets the matrix ``mat1`` to the matrix ``mat2``.


Concatenate
--------------------------------------------------------------------------------


.. function:: void fq_mat_concat_vertical(fq_mat_t res, const fq_mat_t mat1, const fq_mat_t mat2, const fq_ctx_t ctx)

    Sets ``res`` to vertical concatenation of (``mat1``, ``mat2``) in that order. Matrix dimensions : ``mat1`` : `m \times n`, ``mat2`` : `k \times n`, ``res`` : `(m + k) \times n`.

.. function:: void fq_mat_concat_horizontal(fq_mat_t res, const fq_mat_t mat1, const fq_mat_t mat2, const fq_ctx_t ctx)

    Sets ``res`` to horizontal concatenation of (``mat1``, ``mat2``) in that order. Matrix dimensions : ``mat1`` : `m \times n`, ``mat2`` : `m \times k`, ``res``  : `m \times (n + k)`.


Printing
--------------------------------------------------------------------------------


.. function:: int fq_mat_print_pretty(const fq_mat_t mat, const fq_ctx_t ctx)

    Pretty-prints ``mat`` to ``stdout``. A header is printed
    followed by the rows enclosed in brackets.

.. function:: int fq_mat_fprint_pretty(FILE * file, const fq_mat_t mat, const fq_ctx_t ctx)

    Pretty-prints ``mat`` to ``file``. A header is printed
    followed by the rows enclosed in brackets.

    In case of success, returns a positive value.  In case of failure,
    returns a non-positive value.

.. function:: int fq_mat_print(const fq_mat_t mat, const fq_ctx_t ctx)

    Prints ``mat`` to ``stdout``. A header is printed followed
    by the rows enclosed in brackets.

.. function:: int fq_mat_fprint(FILE * file, const fq_mat_t mat, const fq_ctx_t ctx)

    Prints ``mat`` to ``file``. A header is printed followed by
    the rows enclosed in brackets.

    In case of success, returns a positive value.  In case of failure,
    returns a non-positive value.


Window
--------------------------------------------------------------------------------


.. function:: void fq_mat_window_init(fq_mat_t window, const fq_mat_t mat, slong r1, slong c1, slong r2, slong c2, const fq_ctx_t ctx)

     Initializes the matrix ``window`` to be an ``r2 - r1`` by
     ``c2 - c1`` submatrix of ``mat`` whose ``(0,0)`` entry
     is the ``(r1, c1)`` entry of ``mat``.  The memory for the
     elements of ``window`` is shared with ``mat``.


.. function:: void fq_mat_window_clear(fq_mat_t window, const fq_ctx_t ctx)

     Clears the matrix ``window`` and releases any memory that it
     uses.  Note that the memory to the underlying matrix that
     ``window`` points to is not freed.



Random matrix generation
--------------------------------------------------------------------------------


.. function:: void fq_mat_randtest(fq_mat_t mat, flint_rand_t state, const fq_ctx_t ctx)

    Sets the elements of ``mat`` to random elements of
    `\mathbf{F}_{q}`, given by ``ctx``.

.. function:: int fq_mat_randpermdiag(fq_mat_t mat, flint_rand_t state, fq_struct * diag, slong n, const fq_ctx_t ctx)

    Sets ``mat`` to a random permutation of the diagonal matrix
    with `n` leading entries given by the vector ``diag``. It is
    assumed that the main diagonal of ``mat`` has room for at
    least `n` entries.

    Returns `0` or `1`, depending on whether the permutation is even
    or odd respectively.

.. function:: void fq_mat_randrank(fq_mat_t mat, flint_rand_t state, slong rank, const fq_ctx_t ctx)

    Sets ``mat`` to a random sparse matrix with the given rank,
    having exactly as many non-zero elements as the rank, with the
    non-zero elements being uniformly random elements of
    `\mathbf{F}_{q}`.

    The matrix can be transformed into a dense matrix with unchanged
    rank by subsequently calling :func:`fq_mat_randops`.

.. function:: void fq_mat_randops(fq_mat_t mat, flint_rand_t state, slong count, const fq_ctx_t ctx)

    Randomises ``mat`` by performing elementary row or column
    operations. More precisely, at most ``count`` random additions
    or subtractions of distinct rows and columns will be performed.
    This leaves the rank (and for square matrices, determinant)
    unchanged.

.. function:: void fq_mat_randtril(fq_mat_t mat, flint_rand_t state, int unit, const fq_ctx_t ctx)

    Sets ``mat`` to a random lower triangular matrix. If
    ``unit`` is 1, it will have ones on the main diagonal,
    otherwise it will have random nonzero entries on the main
    diagonal.

.. function:: void fq_mat_randtriu(fq_mat_t mat, flint_rand_t state, int unit, const fq_ctx_t ctx)

    Sets ``mat`` to a random upper triangular matrix. If
    ``unit`` is 1, it will have ones on the main diagonal,
    otherwise it will have random nonzero entries on the main
    diagonal.


Comparison
--------------------------------------------------------------------------------


.. function:: int fq_mat_equal(const fq_mat_t mat1, const fq_mat_t mat2, const fq_ctx_t ctx)

    Returns nonzero if mat1 and mat2 have the same dimensions and elements,
    and zero otherwise.

.. function:: int fq_mat_is_zero(const fq_mat_t mat, const fq_ctx_t ctx)

    Returns a non-zero value if all entries of ``mat`` are zero, and
    otherwise returns zero.

.. function:: int fq_mat_is_one(const fq_mat_t mat, const fq_ctx_t ctx)

    Returns a non-zero value if all entries ``mat`` are zero except the
    diagonal entries which must be one, otherwise returns zero..

.. function:: int fq_mat_is_empty(const fq_mat_t mat, const fq_ctx_t ctx)

    Returns a non-zero value if the number of rows or the number of
    columns in ``mat`` is zero, and otherwise returns zero.

.. function:: int fq_mat_is_square(const fq_mat_t mat, const fq_ctx_t ctx)

    Returns a non-zero value if the number of rows is equal to the
    number of columns in ``mat``, and otherwise returns zero.




Addition and subtraction
--------------------------------------------------------------------------------


.. function:: void fq_mat_add(fq_mat_t C, const fq_mat_t A, const fq_mat_t B, const fq_ctx_t ctx)

    Computes `C = A + B`. Dimensions must be identical.

.. function:: void fq_mat_sub(fq_mat_t C, const fq_mat_t A, const fq_mat_t B, const fq_ctx_t ctx)

    Computes `C = A - B`. Dimensions must be identical.

.. function:: void fq_mat_neg(fq_mat_t A, const fq_mat_t B, const fq_ctx_t ctx)

    Sets `B = -A`. Dimensions must be identical.


Matrix multiplication
--------------------------------------------------------------------------------


.. function:: void fq_mat_mul(fq_mat_t C, const fq_mat_t A, const fq_mat_t B, const fq_ctx_t ctx)

    Sets `C = AB`. Dimensions must be compatible for matrix
    multiplication.  Aliasing is allowed. This function automatically chooses
    between classical and KS multiplication.

.. function:: void fq_mat_mul_classical(fq_mat_t C, const fq_mat_t A, const fq_mat_t B, const fq_ctx_t ctx)

    Sets `C = AB`. Dimensions must be compatible for matrix multiplication.
    `C` is not allowed to be aliased with `A` or `B`. Uses classical
    matrix multiplication.

.. function:: void fq_mat_mul_KS(fq_mat_t C, const fq_mat_t A, const fq_mat_t B, const fq_ctx_t ctx)

    Sets `C = AB`. Dimensions must be compatible for matrix
    multiplication.  `C` is not allowed to be aliased with `A` or
    `B`. Uses Kronecker substitution to perform the multiplication
    over the integers.

.. function:: void fq_mat_submul(fq_mat_t D, const fq_mat_t C, const fq_mat_t A, const fq_mat_t B, const fq_ctx_t ctx)

    Sets `D = C + AB`. `C` and `D` may be aliased with each other but
    not with `A` or `B`.

.. function:: void fq_mat_mul_vec(fq_struct * c, const fq_mat_t A, const fq_struct * b, slong blen, const fq_ctx_t ctx)
              void fq_mat_mul_vec_ptr(fq_struct * const * c, const fq_mat_t A, const fq_struct * const * b, slong blen, const fq_ctx_t ctx)

    Compute a matrix-vector product of ``A`` and ``(b, blen)`` and store the result in ``c``.
    The vector ``(b, blen)`` is either truncated or zero-extended to the number of columns of ``A``.
    The number entries written to ``c`` is always equal to the number of rows of ``A``.

.. function:: void fq_mat_vec_mul(fq_struct * c, const fq_struct * a, slong alen, const fq_mat_t B, const fq_ctx_t ctx)
              void fq_mat_vec_mul_ptr(fq_struct * const * c, const fq_struct * const * a, slong alen, const fq_mat_t B, const fq_ctx_t ctx)

    Compute a vector-matrix product of ``(a, alen)`` and ``B`` and and store the result in ``c``.
    The vector ``(a, alen)`` is either truncated or zero-extended to the number of rows of ``B``.
    The number entries written to ``c`` is always equal to the number of columns of ``B``.


Inverse
--------------------------------------------------------------------------------


.. function:: int fq_mat_inv(fq_mat_t B, fq_mat_t A, const fq_ctx_t ctx)

    Sets `B = A^{-1}` and returns `1` if `A` is invertible. If `A` is singular,
    returns `0` and sets the elements of `B` to undefined values.

    `A` and `B` must be square matrices with the same dimensions.


LU decomposition
--------------------------------------------------------------------------------


.. function:: slong fq_mat_lu(slong * P, fq_mat_t A, int rank_check, const fq_ctx_t ctx)

    Computes a generalised LU decomposition `LU = PA` of a given
    matrix `A`, returning the rank of `A`.

    If `A` is a nonsingular square matrix, it will be overwritten with
    a unit diagonal lower triangular matrix `L` and an upper
    triangular matrix `U` (the diagonal of `L` will not be stored
    explicitly).

    If `A` is an arbitrary matrix of rank `r`, `U` will be in row
    echelon form having `r` nonzero rows, and `L` will be lower
    triangular but truncated to `r` columns, having implicit ones on
    the `r` first entries of the main diagonal. All other entries will
    be zero.

    If a nonzero value for ``rank_check`` is passed, the function
    will abandon the output matrix in an undefined state and return 0
    if `A` is detected to be rank-deficient.

    This function calls ``fq_mat_lu_recursive``.

.. function:: slong fq_mat_lu_classical(slong * P, fq_mat_t A, int rank_check, const fq_ctx_t ctx)

    Computes a generalised LU decomposition `LU = PA` of a given
    matrix `A`, returning the rank of `A`. The behavior of this
    function is identical to that of ``fq_mat_lu``. Uses Gaussian
    elimination.

.. function:: slong fq_mat_lu_recursive(slong * P, fq_mat_t A, int rank_check, const fq_ctx_t ctx)

    Computes a generalised LU decomposition `LU = PA` of a given
    matrix `A`, returning the rank of `A`. The behavior of this
    function is identical to that of ``fq_mat_lu``. Uses recursive
    block decomposition, switching to classical Gaussian elimination
    for sufficiently small blocks.


Reduced row echelon form
--------------------------------------------------------------------------------


.. function:: slong fq_mat_rref(fq_mat_t B, const fq_mat_t A, const fq_ctx_t ctx)

    Puts `B` in reduced row echelon form and returns the rank of `A`.

    The rref is computed by first obtaining an unreduced row echelon
    form via LU decomposition and then solving an additional
    triangular system.

.. function:: slong fq_mat_reduce_row(fq_mat_t A, slong * P, slong * L, slong n, const fq_ctx_t ctx)

    Reduce row n of the matrix `A`, assuming the prior rows are in Gauss
    form. However those rows may not be in order. The entry `i` of the array
    `P` is the row of `A` which has a pivot in the `i`-th column. If no such
    row exists, the entry of `P` will be `-1`. The function returns the column
    in which the `n`-th row has a pivot after reduction. This will always be
    chosen to be the first available column for a pivot from the left. This
    information is also updated in `P`. Entry `i` of the array `L` contains the
    number of possibly nonzero columns of `A` row `i`. This speeds up reduction
    in the case that `A` is chambered on the right. Otherwise the entries of
    `L` can all be set to the number of columns of `A`. We require the entries
    of `L` to be monotonic increasing.


Triangular solving
--------------------------------------------------------------------------------


.. function:: void fq_mat_solve_tril(fq_mat_t X, const fq_mat_t L, const fq_mat_t B, int unit, const fq_ctx_t ctx)

    Sets `X = L^{-1} B` where `L` is a full rank lower triangular
    square matrix. If ``unit`` = 1, `L` is assumed to have ones on
    its main diagonal, and the main diagonal will not be read.  `X`
    and `B` are allowed to be the same matrix, but no other aliasing
    is allowed. Automatically chooses between the classical and
    recursive algorithms.

.. function:: void fq_mat_solve_tril_classical(fq_mat_t X, const fq_mat_t L, const fq_mat_t B, int unit, const fq_ctx_t ctx)

    Sets `X = L^{-1} B` where `L` is a full rank lower triangular
    square matrix. If ``unit`` = 1, `L` is assumed to have ones on
    its main diagonal, and the main diagonal will not be read.  `X`
    and `B` are allowed to be the same matrix, but no other aliasing
    is allowed. Uses forward substitution.

.. function:: void fq_mat_solve_tril_recursive(fq_mat_t X, const fq_mat_t L, const fq_mat_t B, int unit, const fq_ctx_t ctx)

    Sets `X = L^{-1} B` where `L` is a full rank lower triangular
    square matrix. If ``unit`` = 1, `L` is assumed to have ones on
    its main diagonal, and the main diagonal will not be read.  `X`
    and `B` are allowed to be the same matrix, but no other aliasing
    is allowed.

    Uses the block inversion formula

    .. math::
        \begin{pmatrix} A & 0 \\ C & D \end{pmatrix}^{-1}
        \begin{pmatrix} X \\ Y \end{pmatrix} =
        \begin{pmatrix} A^{-1} X \\ D^{-1} ( Y - C A^{-1} X ) \end{pmatrix}


    to reduce the problem to matrix multiplication and triangular
    solving of smaller systems.

.. function:: void fq_mat_solve_triu(fq_mat_t X, const fq_mat_t U, const fq_mat_t B, int unit, const fq_ctx_t ctx)

    Sets `X = U^{-1} B` where `U` is a full rank upper triangular
    square matrix. If ``unit`` = 1, `U` is assumed to have ones on
    its main diagonal, and the main diagonal will not be read.  `X`
    and `B` are allowed to be the same matrix, but no other aliasing
    is allowed. Automatically chooses between the classical and
    recursive algorithms.

.. function:: void fq_mat_solve_triu_classical(fq_mat_t X, const fq_mat_t U, const fq_mat_t B, int unit, const fq_ctx_t ctx)

    Sets `X = U^{-1} B` where `U` is a full rank upper triangular
    square matrix. If ``unit`` = 1, `U` is assumed to have ones on
    its main diagonal, and the main diagonal will not be read.  `X`
    and `B` are allowed to be the same matrix, but no other aliasing
    is allowed. Uses forward substitution.

.. function:: void fq_mat_solve_triu_recursive(fq_mat_t X, const fq_mat_t U, const fq_mat_t B, int unit, const fq_ctx_t ctx)

    Sets `X = U^{-1} B` where `U` is a full rank upper triangular
    square matrix. If ``unit`` = 1, `U` is assumed to have ones on
    its main diagonal, and the main diagonal will not be read.  `X`
    and `B` are allowed to be the same matrix, but no other aliasing
    is allowed.

    Uses the block inversion formula

    .. math::
        \begin{pmatrix} A & B \\ 0 & D \end{pmatrix}^{-1}
        \begin{pmatrix} X \\ Y \end{pmatrix} =
        \begin{pmatrix} A^{-1} (X - B D^{-1} Y) \\ D^{-1} Y \end{pmatrix}


    to reduce the problem to matrix multiplication and triangular
    solving of smaller systems.


Solving
--------------------------------------------------------------------------------


.. function:: int fq_mat_solve(fq_mat_t X, const fq_mat_t A, const fq_mat_t B, const fq_ctx_t ctx)

    Solves the matrix-matrix equation `AX = B`.

    Returns `1` if `A` has full rank; otherwise returns `0` and sets the
    elements of `X` to undefined values.

    The matrix `A` must be square.

.. function:: int fq_mat_can_solve(fq_mat_t X, const fq_mat_t A, const fq_mat_t B, const fq_ctx_t ctx)

    Solves the matrix-matrix equation `AX = B` over `Fq`.

    Returns `1` if a solution exists; otherwise returns `0` and sets the
    elements of `X` to zero. If more than one solution exists, one of the
    valid solutions is given.

    There are no restrictions on the shape of `A` and it may be singular.


Transforms
--------------------------------------------------------------------------------


.. function:: void fq_mat_similarity(fq_mat_t M, slong r, fq_t d, const fq_ctx_t ctx)

    Applies a similarity transform to the `n\times n` matrix `M` in-place.

    If `P` is the `n\times n` identity matrix the zero entries of whose row
    `r` (`0`-indexed) have been replaced by `d`, this transform is equivalent
    to `M = P^{-1}MP`.

    Similarity transforms preserve the determinant, characteristic polynomial
    and minimal polynomial.

    The value `d` is required to be reduced modulo the modulus of the entries
    in the matrix.


Characteristic polynomial
--------------------------------------------------------------------------------


.. function:: void fq_mat_charpoly_danilevsky(fq_poly_t p, const fq_mat_t M, const fq_ctx_t ctx)

    Compute the characteristic polynomial `p` of the matrix `M`. The matrix
    is assumed to be square.

.. function:: void fq_mat_charpoly(fq_poly_t p, const fq_mat_t M, const fq_ctx_t ctx)

    Compute the characteristic polynomial `p` of the matrix `M`. The matrix
    is required to be square, otherwise an exception is raised.


Minimal polynomial
--------------------------------------------------------------------------------


.. function:: void fq_mat_minpoly(fq_poly_t p, const fq_mat_t M, const fq_ctx_t ctx)

    Compute the minimal polynomial `p` of the matrix `M`. The matrix
    is required to be square, otherwise an exception is raised.
