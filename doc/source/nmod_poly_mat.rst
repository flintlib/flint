.. _nmod-poly-mat:

**nmod_poly_mat.h** -- matrices of univariate polynomials over integers mod n (word-size n)
===========================================================================================

Description.

Types, macros and constants
-------------------------------------------------------------------------------

.. type:: nmod_poly_mat_struct

.. type:: nmod_poly_mat_t

    Description.


Memory management
--------------------------------------------------------------------------------


.. function:: void nmod_poly_mat_init(nmod_poly_mat_t mat, slong rows, slong cols, mp_limb_t n)

    Initialises a matrix with the given number of rows and columns for use.
    The modulus is set to `n`.

.. function:: void nmod_poly_mat_init_set(nmod_poly_mat_t mat, const nmod_poly_mat_t src)

    Initialises a matrix ``mat`` of the same dimensions and modulus
    as ``src``, and sets it to a copy of ``src``.

.. function:: void nmod_poly_mat_clear(nmod_poly_mat_t mat)

    Frees all memory associated with the matrix. The matrix must be
    reinitialised if it is to be used again.


Basic properties
--------------------------------------------------------------------------------


.. function:: slong nmod_poly_mat_nrows(const nmod_poly_mat_t mat)

    Returns the number of rows in ``mat``.

.. function:: slong nmod_poly_mat_ncols(const nmod_poly_mat_t mat)

    Returns the number of columns in ``mat``.

.. function:: mp_limb_t nmod_poly_mat_modulus(const nmod_poly_mat_t mat)

    Returns the modulus of ``mat``.


Basic assignment and manipulation
--------------------------------------------------------------------------------


.. function:: nmod_poly_struct * nmod_poly_mat_entry(const nmod_poly_mat_t mat, slong i, slong j)

    Gives a reference to the entry at row ``i`` and column ``j``.
    The reference can be passed as an input or output variable to any
    ``nmod_poly`` function for direct manipulation of the matrix element.
    No bounds checking is performed.

.. function:: void nmod_poly_mat_set(nmod_poly_mat_t mat1, const nmod_poly_mat_t mat2)

    Sets ``mat1`` to a copy of ``mat2``.

.. function:: void nmod_poly_mat_swap(nmod_poly_mat_t mat1, nmod_poly_mat_t mat2)

    Swaps ``mat1`` and ``mat2`` efficiently.



Input and output
--------------------------------------------------------------------------------


.. function:: void nmod_poly_mat_print(const nmod_poly_mat_t mat, const char * x)

    Prints the matrix ``mat`` to standard output, using the
    variable ``x``.


Random matrix generation
--------------------------------------------------------------------------------


.. function:: void nmod_poly_mat_randtest(nmod_poly_mat_t mat, flint_rand_t state, slong len)

    This is equivalent to applying ``nmod_poly_randtest`` to all entries
    in the matrix.

.. function:: void nmod_poly_mat_randtest_sparse(nmod_poly_mat_t A, flint_rand_t state, slong len, float density)

    Creates a random matrix with the amount of nonzero entries given
    approximately by the ``density`` variable, which should be a fraction
    between 0 (most sparse) and 1 (most dense).

    The nonzero entries will have random lengths between 1 and ``len``.


Special matrices
--------------------------------------------------------------------------------


.. function:: void nmod_poly_mat_zero(nmod_poly_mat_t mat)

    Sets ``mat`` to the zero matrix.

.. function:: void nmod_poly_mat_one(nmod_poly_mat_t mat)

    Sets ``mat`` to the unit or identity matrix of given shape,
    having the element 1 on the main diagonal and zeros elsewhere.
    If ``mat`` is nonsquare, it is set to the truncation of a unit matrix.


Basic comparison and properties
--------------------------------------------------------------------------------


.. function:: int nmod_poly_mat_equal(const nmod_poly_mat_t mat1, const nmod_poly_mat_t mat2)

    Returns nonzero if ``mat1`` and ``mat2`` have the same shape and
    all their entries agree, and returns zero otherwise.

.. function:: int nmod_poly_mat_is_zero(const nmod_poly_mat_t mat)

    Returns nonzero if all entries in ``mat`` are zero, and returns
    zero otherwise.

.. function:: int nmod_poly_mat_is_one(const nmod_poly_mat_t mat)

    Returns nonzero if all entry of ``mat`` on the main diagonal
    are the constant polynomial 1 and all remaining entries are zero,
    and returns zero otherwise. The matrix need not be square.

.. function:: int nmod_poly_mat_is_empty(const nmod_poly_mat_t mat)

    Returns a non-zero value if the number of rows or the number of
    columns in ``mat`` is zero, and otherwise returns
    zero.

.. function:: int nmod_poly_mat_is_square(const nmod_poly_mat_t mat)

    Returns a non-zero value if the number of rows is equal to the
    number of columns in ``mat``, and otherwise returns zero.



Norms
--------------------------------------------------------------------------------


.. function:: slong nmod_poly_mat_max_length(const nmod_poly_mat_t A)

    Returns the maximum polynomial length among all the entries in ``A``.



Evaluation
--------------------------------------------------------------------------------


.. function:: void nmod_poly_mat_evaluate_nmod(nmod_mat_t B, const nmod_poly_mat_t A, mp_limb_t x)

    Sets the ``nmod_mat_t`` ``B`` to ``A`` evaluated entrywise
    at the point ``x``.



Arithmetic
--------------------------------------------------------------------------------


.. function:: void nmod_poly_mat_scalar_mul_nmod_poly(nmod_poly_mat_t B, const nmod_poly_mat_t A, const nmod_poly_t c)

    Sets ``B`` to ``A`` multiplied entrywise by the polynomial ``c``.

.. function:: void nmod_poly_mat_scalar_mul_nmod(nmod_poly_mat_t B, const nmod_poly_mat_t A, mp_limb_t c)

    Sets ``B`` to ``A`` multiplied entrywise by the coefficient
    ``c``, which is assumed to be reduced modulo the modulus.

.. function:: void nmod_poly_mat_add(nmod_poly_mat_t C, const nmod_poly_mat_t A, const nmod_poly_mat_t B)

    Sets ``C`` to the sum of ``A`` and ``B``.
    All matrices must have the same shape. Aliasing is allowed.

.. function:: void nmod_poly_mat_sub(nmod_poly_mat_t C, const nmod_poly_mat_t A, const nmod_poly_mat_t B)

    Sets ``C`` to the sum of ``A`` and ``B``.
    All matrices must have the same shape. Aliasing is allowed.

.. function:: void nmod_poly_mat_neg(nmod_poly_mat_t B, const nmod_poly_mat_t A)

    Sets ``B`` to the negation of ``A``.
    The matrices must have the same shape. Aliasing is allowed.

.. function:: void nmod_poly_mat_mul(nmod_poly_mat_t C, const nmod_poly_mat_t A, const nmod_poly_mat_t B)

    Sets ``C`` to the matrix product of ``A`` and ``B``.
    The matrices must have compatible dimensions for matrix multiplication.
    Aliasing is allowed. This function automatically chooses between
    classical, KS and evaluation-interpolation multiplication.

.. function:: void nmod_poly_mat_mul_classical(nmod_poly_mat_t C, const nmod_poly_mat_t A, const nmod_poly_mat_t B)

    Sets ``C`` to the matrix product of ``A`` and ``B``, 
    computed using the classical algorithm. The matrices must have 
    compatible dimensions for matrix multiplication. Aliasing is allowed.

.. function:: void nmod_poly_mat_mul_KS(nmod_poly_mat_t C, const nmod_poly_mat_t A, const nmod_poly_mat_t B)

    Sets ``C`` to the matrix product of ``A`` and ``B``, 
    computed using Kronecker segmentation. The matrices must have 
    compatible dimensions for matrix multiplication. Aliasing is allowed.

.. function:: void nmod_poly_mat_mul_interpolate(nmod_poly_mat_t C, const nmod_poly_mat_t A, const nmod_poly_mat_t B)

    Sets ``C`` to the matrix product of ``A`` and ``B``,
    computed through evaluation and interpolation. The matrices must have
    compatible dimensions for matrix multiplication. For interpolation
    to be well-defined, we require that the modulus is a prime at least as
    large as `m + n - 1` where `m` and `n` are the maximum lengths of
    polynomials in the input matrices. Aliasing is allowed.

.. function:: void nmod_poly_mat_sqr(nmod_poly_mat_t B, const nmod_poly_mat_t A)

    Sets ``B`` to the square of ``A``, which must be a square matrix.
    Aliasing is allowed. This function automatically chooses between
    classical and KS squaring.

.. function:: void nmod_poly_mat_sqr_classical(nmod_poly_mat_t B, const nmod_poly_mat_t A)

    Sets ``B`` to the square of ``A``, which must be a square matrix.
    Aliasing is allowed. This function uses direct formulas for very small
    matrices, and otherwise classical matrix multiplication.

.. function:: void nmod_poly_mat_sqr_KS(nmod_poly_mat_t B, const nmod_poly_mat_t A)

    Sets ``B`` to the square of ``A``, which must be a square matrix.
    Aliasing is allowed. This function uses Kronecker segmentation.

.. function:: void nmod_poly_mat_sqr_interpolate(nmod_poly_mat_t B, const nmod_poly_mat_t A)

    Sets ``B`` to the square of ``A``, which must be a square matrix,
    computed through evaluation and interpolation. For interpolation
    to be well-defined, we require that the modulus is a prime at least as
    large as `2n - 1` where `n` is the maximum length of
    polynomials in the input matrix. Aliasing is allowed.

.. function:: void nmod_poly_mat_pow(nmod_poly_mat_t B, const nmod_poly_mat_t A, ulong exp)

    Sets ``B`` to ``A`` raised to the power ``exp``, where ``A``
    is a square matrix. Uses exponentiation by squaring. Aliasing is allowed.


Row reduction
--------------------------------------------------------------------------------


.. function:: slong nmod_poly_mat_find_pivot_any(const nmod_poly_mat_t mat, slong start_row, slong end_row, slong c)

    Attempts to find a pivot entry for row reduction.
    Returns a row index `r` between ``start_row`` (inclusive) and
    ``stop_row`` (exclusive) such that column `c` in ``mat`` has
    a nonzero entry on row `r`, or returns -1 if no such entry exists.

    This implementation simply chooses the first nonzero entry from
    it encounters. This is likely to be a nearly optimal choice if all
    entries in the matrix have roughly the same size, but can lead to
    unnecessary coefficient growth if the entries vary in size.

.. function:: slong nmod_poly_mat_find_pivot_partial(const nmod_poly_mat_t mat, slong start_row, slong end_row, slong c)

    Attempts to find a pivot entry for row reduction.
    Returns a row index `r` between ``start_row`` (inclusive) and
    ``stop_row`` (exclusive) such that column `c` in ``mat`` has
    a nonzero entry on row `r`, or returns -1 if no such entry exists.

    This implementation searches all the rows in the column and
    chooses the nonzero entry of smallest degree. This heuristic
    typically reduces coefficient growth when the matrix entries
    vary in size.

.. function:: slong nmod_poly_mat_fflu(nmod_poly_mat_t B, nmod_poly_t den, slong * perm, const nmod_poly_mat_t A, int rank_check)

    Uses fraction-free Gaussian elimination to set (``B``, ``den``) to a
    fraction-free LU decomposition of ``A`` and returns the
    rank of ``A``. Aliasing of ``A`` and ``B`` is allowed.

    Pivot elements are chosen with ``nmod_poly_mat_find_pivot_partial``.
    If ``perm`` is non-``NULL``, the permutation of
    rows in the matrix will also be applied to ``perm``.

    If ``rank_check`` is set, the function aborts and returns 0 if the
    matrix is detected not to have full rank without completing the
    elimination.

    The denominator ``den`` is set to `\pm \operatorname{det}(A)`, where
    the sign is decided by the parity of the permutation. Note that the
    determinant is not generally the minimal denominator.

.. function:: slong nmod_poly_mat_rref(nmod_poly_mat_t B, nmod_poly_t den, const nmod_poly_mat_t A)

    Sets (``B``, ``den``) to the reduced row echelon form of ``A``
    and returns the rank of ``A``.
    Aliasing of ``A`` and ``B`` is allowed.

    The denominator ``den`` is set to `\pm \operatorname{det}(A)`.
    Note that the determinant is not generally the minimal denominator.


Trace
--------------------------------------------------------------------------------


.. function:: void nmod_poly_mat_trace(nmod_poly_t trace, const nmod_poly_mat_t mat)

    Computes the trace of the matrix, i.e. the sum of the entries on
    the main diagonal. The matrix is required to be square.


Determinant and rank
--------------------------------------------------------------------------------


.. function:: void nmod_poly_mat_det(nmod_poly_t det, const nmod_poly_mat_t A)

    Sets ``det`` to the determinant of the square matrix ``A``. Uses
    a direct formula, fraction-free LU decomposition, or interpolation,
    depending on the size of the matrix.

.. function:: void nmod_poly_mat_det_fflu(nmod_poly_t det, const nmod_poly_mat_t A)

    Sets ``det`` to the determinant of the square matrix ``A``.
    The determinant is computed by performing a fraction-free LU
    decomposition on a copy of ``A``.

.. function:: void nmod_poly_mat_det_interpolate(nmod_poly_t det, const nmod_poly_mat_t A)

    Sets ``det`` to the determinant of the square matrix ``A``.
    The determinant is computed by determining a bound `n` for its length,
    evaluating the matrix at `n` distinct points, computing the determinant
    of each coefficient matrix, and forming the interpolating polynomial.

    If the coefficient ring does not contain `n` distinct points (that is,
    if working over `\mathbf{Z}/p\mathbf{Z}` where `p < n`),
    this function automatically falls back to ``nmod_poly_mat_det_fflu``.

.. function:: slong nmod_poly_mat_rank(const nmod_poly_mat_t A)

    Returns the rank of ``A``. Performs fraction-free LU decomposition
    on a copy of ``A``.



Inverse
--------------------------------------------------------------------------------


.. function:: int nmod_poly_mat_inv(nmod_poly_mat_t Ainv, nmod_poly_t den, const nmod_poly_mat_t A)

    Sets (``Ainv``, ``den``) to the inverse matrix of ``A``.
    Returns 1 if ``A`` is nonsingular and 0 if ``A`` is singular.
    Aliasing of ``Ainv`` and ``A`` is allowed.

    More precisely, ``det`` will be set to the determinant of ``A``
    and ``Ainv`` will be set to the adjugate matrix of ``A``.
    Note that the determinant is not necessarily the minimal denominator.

    Uses fraction-free LU decomposition, followed by solving for
    the identity matrix.



Nullspace
--------------------------------------------------------------------------------


.. function:: slong nmod_poly_mat_nullspace(nmod_poly_mat_t res, const nmod_poly_mat_t mat)

    Computes the right rational nullspace of the matrix ``mat`` and
    returns the nullity.

    More precisely, assume that ``mat`` has rank `r` and nullity `n`.
    Then this function sets the first `n` columns of ``res``
    to linearly independent vectors spanning the nullspace of ``mat``.
    As a result, we always have rank(``res``) `= n`, and
    ``mat`` `\times` ``res`` is the zero matrix.

    The computed basis vectors will not generally be in a reduced form.
    In general, the polynomials in each column vector in the result
    will have a nontrivial common GCD.


Solving
--------------------------------------------------------------------------------


.. function:: int nmod_poly_mat_solve(nmod_poly_mat_t X, nmod_poly_t den, const nmod_poly_mat_t A, const nmod_poly_mat_t B)

    Solves the equation `AX = B` for nonsingular `A`. More precisely, computes
    (``X``, ``den``) such that `AX = B \times \operatorname{den}`.
    Returns 1 if `A` is nonsingular and 0 if `A` is singular.
    The computed denominator will not generally be minimal.

    Uses fraction-free LU decomposition followed by fraction-free
    forward and back substitution.

.. function:: int nmod_poly_mat_solve_fflu(nmod_poly_mat_t X, nmod_poly_t den, const nmod_poly_mat_t A, const nmod_poly_mat_t B)

    Solves the equation `AX = B` for nonsingular `A`. More precisely, computes
    (``X``, ``den``) such that `AX = B \times \operatorname{den}`.
    Returns 1 if `A` is nonsingular and 0 if `A` is singular.
    The computed denominator will not generally be minimal.

    Uses fraction-free LU decomposition followed by fraction-free
    forward and back substitution.

.. function:: void nmod_poly_mat_solve_fflu_precomp(nmod_poly_mat_t X, const slong * perm, const nmod_poly_mat_t FFLU, const nmod_poly_mat_t B)

    Performs fraction-free forward and back substitution given a precomputed
    fraction-free LU decomposition and corresponding permutation.
