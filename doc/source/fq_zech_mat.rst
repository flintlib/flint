.. _fq-zech-mat:

**fq_zech_mat.h** -- matrices over finite fields (Zech logarithm representation)
================================================================================

Types, macros and constants
-------------------------------------------------------------------------------

.. type:: fq_zech_mat_struct

.. type:: fq_zech_mat_t

Memory management
--------------------------------------------------------------------------------


.. function:: void fq_zech_mat_init(fq_zech_mat_t mat, slong rows, slong cols, const fq_zech_ctx_t ctx)

    Initialises ``mat`` to a ``rows``-by-``cols`` matrix with
    coefficients in `\mathbf{F}_{q}` given by ``ctx``. All elements
    are set to zero.

.. function:: void fq_zech_mat_init_set(fq_zech_mat_t mat, const fq_zech_mat_t src, const fq_zech_ctx_t ctx)

    Initialises ``mat`` and sets its dimensions and elements to
    those of ``src``.

.. function:: void fq_zech_mat_clear(fq_zech_mat_t mat, const fq_zech_ctx_t ctx)

    Clears the matrix and releases any memory it used. The matrix
    cannot be used again until it is initialised. This function must be
    called exactly once when finished using an ``fq_zech_mat_t`` object.

.. function:: void fq_zech_mat_set(fq_zech_mat_t mat, const fq_zech_mat_t src, const fq_zech_ctx_t ctx)

    Sets ``mat`` to a copy of ``src``. It is assumed
    that ``mat`` and ``src`` have identical dimensions.


Basic properties and manipulation
--------------------------------------------------------------------------------


.. function:: fq_zech_struct * fq_zech_mat_entry(const fq_zech_mat_t mat, slong i, slong j)

    Directly accesses the entry in ``mat`` in row `i` and column `j`,
    indexed from zero. No bounds checking is performed.

.. function:: void fq_zech_mat_entry_set(fq_zech_mat_t mat, slong i, slong j, const fq_zech_t x, const fq_zech_ctx_t ctx)

    Sets the entry in ``mat`` in row `i` and column `j` to ``x``.

.. function:: slong fq_zech_mat_nrows(const fq_zech_mat_t mat, const fq_zech_ctx_t ctx)

    Returns the number of rows in ``mat``.

.. function:: slong fq_zech_mat_ncols(const fq_zech_mat_t mat, const fq_zech_ctx_t ctx)

    Returns the number of columns in ``mat``.

.. function:: void fq_zech_mat_swap(fq_zech_mat_t mat1, fq_zech_mat_t mat2, const fq_zech_ctx_t ctx)

    Swaps two matrices. The dimensions of ``mat1`` and ``mat2``
    are allowed to be different.

.. function:: void fq_zech_mat_swap_entrywise(fq_zech_mat_t mat1, fq_zech_mat_t mat2, const fq_zech_ctx_t ctx)

    Swaps two matrices by swapping the individual entries rather than swapping
    the contents of the structs.

.. function:: void fq_zech_mat_zero(fq_zech_mat_t mat, const fq_zech_ctx_t ctx)

    Sets all entries of ``mat`` to 0.

.. function:: void fq_zech_mat_one(fq_zech_mat_t mat, const fq_zech_ctx_t ctx)

    Sets all diagonal entries of ``mat`` to 1 and all other entries to 0.


Conversions
--------------------------------------------------------------------------------

.. function:: void fq_zech_mat_set_nmod_mat(fq_zech_mat_t mat1, const nmod_mat_t mat2, const fq_zech_ctx_t ctx)

    Sets the matrix ``mat1`` to the matrix ``mat2``.

.. function:: void fq_zech_mat_set_fmpz_mod_mat(fq_zech_mat_t mat1, const fmpz_mod_mat_t mat2, const fq_zech_ctx_t ctx)

    Sets the matrix ``mat1`` to the matrix ``mat2``.


Concatenate
--------------------------------------------------------------------------------


.. function:: void fq_zech_mat_concat_vertical(fq_zech_mat_t res, const fq_zech_mat_t mat1, const fq_zech_mat_t mat2, const fq_zech_ctx_t ctx)

    Sets ``res`` to vertical concatenation of (``mat1``, ``mat2``) in that order. Matrix dimensions : ``mat1`` : `m \times n`, ``mat2`` : `k \times n`, ``res`` : `(m + k) \times n`.

.. function:: void fq_zech_mat_concat_horizontal(fq_zech_mat_t res, const fq_zech_mat_t mat1, const fq_zech_mat_t mat2, const fq_zech_ctx_t ctx)

    Sets ``res`` to horizontal concatenation of (``mat1``, ``mat2``) in that order. Matrix dimensions : ``mat1`` : `m \times n`, ``mat2`` : `m \times k`, ``res``  : `m \times (n + k)`.


Printing
--------------------------------------------------------------------------------


.. function:: int fq_zech_mat_print_pretty(const fq_zech_mat_t mat, const fq_zech_ctx_t ctx)

    Pretty-prints ``mat`` to ``stdout``. A header is printed
    followed by the rows enclosed in brackets.

.. function:: int fq_zech_mat_fprint_pretty(FILE * file, const fq_zech_mat_t mat, const fq_zech_ctx_t ctx)

    Pretty-prints ``mat`` to ``file``. A header is printed
    followed by the rows enclosed in brackets.

    In case of success, returns a positive value.  In case of failure,
    returns a non-positive value.

.. function:: int fq_zech_mat_print(const fq_zech_mat_t mat, const fq_zech_ctx_t ctx)

    Prints ``mat`` to ``stdout``. A header is printed followed
    by the rows enclosed in brackets.

.. function:: int fq_zech_mat_fprint(FILE * file, const fq_zech_mat_t mat, const fq_zech_ctx_t ctx)

    Prints ``mat`` to ``file``. A header is printed followed by
    the rows enclosed in brackets.

    In case of success, returns a positive value.  In case of failure,
    returns a non-positive value.


Window
--------------------------------------------------------------------------------


.. function:: void fq_zech_mat_window_init(fq_zech_mat_t window, const fq_zech_mat_t mat, slong r1, slong c1, slong r2, slong c2, const fq_zech_ctx_t ctx)

     Initializes the matrix ``window`` to be an ``r2 - r1`` by
     ``c2 - c1`` submatrix of ``mat`` whose ``(0,0)`` entry
     is the ``(r1, c1)`` entry of ``mat``.  The memory for the
     elements of ``window`` is shared with ``mat``.


.. function:: void fq_zech_mat_window_clear(fq_zech_mat_t window, const fq_zech_ctx_t ctx)

     Clears the matrix ``window`` and releases any memory that it
     uses.  Note that the memory to the underlying matrix that
     ``window`` points to is not freed.


Random matrix generation
--------------------------------------------------------------------------------


.. function:: void fq_zech_mat_randtest(fq_zech_mat_t mat, flint_rand_t state, const fq_zech_ctx_t ctx)

    Sets the elements of ``mat`` to random elements of
    `\mathbf{F}_{q}`, given by ``ctx``.

.. function:: int fq_zech_mat_randpermdiag(fq_zech_mat_t mat, flint_rand_t state, fq_zech_struct * diag, slong n, const fq_zech_ctx_t ctx)

    Sets ``mat`` to a random permutation of the diagonal matrix
    with `n` leading entries given by the vector ``diag``. It is
    assumed that the main diagonal of ``mat`` has room for at
    least `n` entries.

    Returns `0` or `1`, depending on whether the permutation is even
    or odd respectively.

.. function:: void fq_zech_mat_randrank(fq_zech_mat_t mat, flint_rand_t state, slong rank, const fq_zech_ctx_t ctx)

    Sets ``mat`` to a random sparse matrix with the given rank,
    having exactly as many non-zero elements as the rank, with the
    non-zero elements being uniformly random elements of
    `\mathbf{F}_{q}`.

    The matrix can be transformed into a dense matrix with unchanged
    rank by subsequently calling :func:`fq_zech_mat_randops`.

.. function:: void fq_zech_mat_randops(fq_zech_mat_t mat, flint_rand_t state, slong count, const fq_zech_ctx_t ctx)

    Randomises ``mat`` by performing elementary row or column
    operations. More precisely, at most ``count`` random additions
    or subtractions of distinct rows and columns will be performed.
    This leaves the rank (and for square matrices, determinant)
    unchanged.

.. function:: void fq_zech_mat_randtril(fq_zech_mat_t mat, flint_rand_t state, int unit, const fq_zech_ctx_t ctx)

    Sets ``mat`` to a random lower triangular matrix. If
    ``unit`` is 1, it will have ones on the main diagonal,
    otherwise it will have random nonzero entries on the main
    diagonal.

.. function:: void fq_zech_mat_randtriu(fq_zech_mat_t mat, flint_rand_t state, int unit, const fq_zech_ctx_t ctx)

    Sets ``mat`` to a random upper triangular matrix. If
    ``unit`` is 1, it will have ones on the main diagonal,
    otherwise it will have random nonzero entries on the main
    diagonal.


Comparison
--------------------------------------------------------------------------------


.. function:: int fq_zech_mat_equal(const fq_zech_mat_t mat1, const fq_zech_mat_t mat2, const fq_zech_ctx_t ctx)

    Returns nonzero if mat1 and mat2 have the same dimensions and elements,
    and zero otherwise.

.. function:: int fq_zech_mat_is_zero(const fq_zech_mat_t mat, const fq_zech_ctx_t ctx)

    Returns a non-zero value if all entries ``mat`` are zero, and
    otherwise returns zero.

.. function:: int fq_zech_mat_is_one(const fq_zech_mat_t mat, const fq_zech_ctx_t ctx)

    Returns a non-zero value if all entries ``mat`` are zero except the
    diagonal entries which must be one, otherwise returns zero.

.. function:: int fq_zech_mat_is_empty(const fq_zech_mat_t mat, const fq_zech_ctx_t ctx)

    Returns a non-zero value if the number of rows or the number of
    columns in ``mat`` is zero, and otherwise returns zero.

.. function:: int fq_zech_mat_is_square(const fq_zech_mat_t mat, const fq_zech_ctx_t ctx)

    Returns a non-zero value if the number of rows is equal to the
    number of columns in ``mat``, and otherwise returns zero.



Transpose
--------------------------------------------------------------------------------


.. function:: void fq_zech_mat_transpose(fq_zech_mat_t B, const fq_zech_mat_t A, const fq_zech_ctx_t ctx)

    Sets `B` to `A^T`, the transpose of `A`. Dimensions must be compatible.
    Aliasing is allowed for square matrices.



Addition and subtraction
--------------------------------------------------------------------------------


.. function:: void fq_zech_mat_add(fq_zech_mat_t C, const fq_zech_mat_t A, const fq_zech_mat_t B,  const fq_zech_ctx_t ctx)

    Computes `C = A + B`. Dimensions must be identical.

.. function:: void fq_zech_mat_sub(fq_zech_mat_t C, const fq_zech_mat_t A, const fq_zech_mat_t B, const fq_zech_ctx_t ctx)

    Computes `C = A - B`. Dimensions must be identical.

.. function:: void fq_zech_mat_neg(fq_zech_mat_t A, const fq_zech_mat_t B, const fq_zech_ctx_t ctx)

    Sets `B = -A`. Dimensions must be identical.


Matrix multiplication
--------------------------------------------------------------------------------


.. function:: void fq_zech_mat_mul(fq_zech_mat_t C, const fq_zech_mat_t A, const fq_zech_mat_t B,  const fq_zech_ctx_t ctx)

    Sets `C = AB`. Dimensions must be compatible for matrix
    multiplication.  `C` is not allowed to be aliased with `A` or
    `B`. This function automatically chooses between classical and
    KS multiplication.

.. function:: void fq_zech_mat_mul_classical(fq_zech_mat_t C, const fq_zech_mat_t A, const fq_zech_mat_t B, const fq_zech_ctx_t ctx)

    Sets `C = AB`. Dimensions must be compatible for matrix multiplication.
    `C` is not allowed to be aliased with `A` or `B`. Uses classical
    matrix multiplication.

.. function:: void fq_zech_mat_mul_KS(fq_zech_mat_t C, const fq_zech_mat_t A, const fq_zech_mat_t B, const fq_zech_ctx_t ctx)

    Sets `C = AB`. Dimensions must be compatible for matrix
    multiplication.  `C` is not allowed to be aliased with `A` or
    `B`. Uses Kronecker substitution to perform the multiplication
    over the integers.

.. function:: void fq_zech_mat_submul(fq_zech_mat_t D, const fq_zech_mat_t C, const fq_zech_mat_t A, const fq_zech_mat_t B, const fq_zech_ctx_t ctx)

    Sets `D = C + AB`. `C` and `D` may be aliased with each other but
    not with `A` or `B`.

.. function:: void fq_zech_mat_mul_vec(fq_zech_struct * c, const fq_zech_mat_t A, const fq_zech_struct * b, slong blen, const fq_zech_ctx_t ctx)
              void fq_zech_mat_mul_vec_ptr(fq_zech_struct * const * c, const fq_zech_mat_t A, const fq_zech_struct * const * b, slong blen, const fq_zech_ctx_t ctx)

    Compute a matrix-vector product of ``A`` and ``(b, blen)`` and store the result in ``c``.
    The vector ``(b, blen)`` is either truncated or zero-extended to the number of columns of ``A``.
    The number entries written to ``c`` is always equal to the number of rows of ``A``.

.. function:: void fq_zech_mat_vec_mul(fq_zech_struct * c, const fq_zech_struct * a, slong alen, const fq_zech_mat_t B, const fq_zech_ctx_t ctx)
              void fq_zech_mat_vec_mul_ptr(fq_zech_struct * const * c, const fq_zech_struct * const * a, slong alen, const fq_zech_mat_t B, const fq_zech_ctx_t ctx)

    Compute a vector-matrix product of ``(a, alen)`` and ``B`` and and store the result in ``c``.
    The vector ``(a, alen)`` is either truncated or zero-extended to the number of rows of ``B``.
    The number entries written to ``c`` is always equal to the number of columns of ``B``.


LU decomposition
--------------------------------------------------------------------------------


.. function:: slong fq_zech_mat_lu(slong * P, fq_zech_mat_t A, int rank_check, const fq_zech_ctx_t ctx)

    Computes a generalised LU decomposition `PLU = A` of a given
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


Reduced row echelon form
--------------------------------------------------------------------------------


.. function:: slong fq_zech_mat_rref(fq_zech_mat_t B, const fq_zech_mat_t A, const fq_zech_ctx_t ctx)

    Puts `B` in reduced row echelon form and returns the rank of `A`.

    The rref is computed by first obtaining an unreduced row echelon
    form via LU decomposition and then solving an additional
    triangular system.


Triangular solving
--------------------------------------------------------------------------------


.. function:: void fq_zech_mat_solve_tril(fq_zech_mat_t X, const fq_zech_mat_t L, const fq_zech_mat_t B, int unit, const fq_zech_ctx_t ctx)

    Sets `X = L^{-1} B` where `L` is a full rank lower triangular
    square matrix. If ``unit`` = 1, `L` is assumed to have ones on
    its main diagonal, and the main diagonal will not be read.  `X`
    and `B` are allowed to be the same matrix, but no other aliasing
    is allowed. Automatically chooses between the classical and
    recursive algorithms.

.. function:: void fq_zech_mat_solve_triu(fq_zech_mat_t X, const fq_zech_mat_t U, const fq_zech_mat_t B, int unit, const fq_zech_ctx_t ctx)

    Sets `X = U^{-1} B` where `U` is a full rank upper triangular
    square matrix. If ``unit`` = 1, `U` is assumed to have ones on
    its main diagonal, and the main diagonal will not be read.  `X`
    and `B` are allowed to be the same matrix, but no other aliasing
    is allowed. Automatically chooses between the classical and
    recursive algorithms.


Solving
--------------------------------------------------------------------------------

.. function:: int fq_zech_mat_solve(fq_zech_mat_t X, const fq_zech_mat_t A, const fq_zech_mat_t B, const fq_zech_ctx_t ctx)

    Solves the matrix-matrix equation `AX = B`.

    Returns `1` if `A` has full rank; otherwise returns `0` and sets the
    elements of `X` to undefined values.

    The matrix `A` must be square.

.. function:: int fq_zech_mat_can_solve(fq_zech_mat_t X, const fq_zech_mat_t A, const fq_zech_mat_t B, const fq_zech_ctx_t ctx)

    Solves the matrix-matrix equation `AX = B` over `Fq`.

    Returns `1` if a solution exists; otherwise returns `0` and sets the
    elements of `X` to zero. If more than one solution exists, one of the
    valid solutions is given.

    There are no restrictions on the shape of `A` and it may be singular.

Transforms
--------------------------------------------------------------------------------


.. function:: void fq_zech_mat_similarity(fq_zech_mat_t M, slong r, fq_zech_t d, const fq_zech_ctx_t ctx)

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


.. function:: void fq_zech_mat_charpoly(fq_zech_poly_t p, const fq_zech_mat_t M, const fq_zech_ctx_t ctx)

    Compute the characteristic polynomial `p` of the matrix `M`. The matrix
    is required to be square, otherwise an exception is raised.


Minimal polynomial
--------------------------------------------------------------------------------


.. function:: void fq_zech_mat_minpoly(fq_zech_poly_t p, const fq_zech_mat_t M, const fq_zech_ctx_t ctx)

    Compute the minimal polynomial `p` of the matrix `M`. The matrix
    is required to be square, otherwise an exception is raised.
