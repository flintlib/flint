.. _fmpq-mat:

**fmpq_mat.h** -- matrices over the rational numbers
===============================================================================

Description.

Types, macros and constants
-------------------------------------------------------------------------------

.. type:: fmpq_mat_struct

.. type:: fmpq_mat_t

    Description.


Memory management
--------------------------------------------------------------------------------


.. function:: void fmpq_mat_init(fmpq_mat_t mat, slong rows, slong cols)

    Initialises a matrix with the given number of rows and columns for use.

.. function:: void fmpq_mat_init_set(fmpq_mat_t mat1, const fmpq_mat_t mat2)

    Initialises ``mat1`` and sets it equal to ``mat2``.

.. function:: void fmpq_mat_clear(fmpq_mat_t mat)

    Frees all memory associated with the matrix. The matrix must be
    reinitialised if it is to be used again.

.. function:: void fmpq_mat_swap(fmpq_mat_t mat1, fmpq_mat_t mat2)

    Swaps two matrices. The dimensions of ``mat1`` and ``mat2``
    are allowed to be different.

.. function:: void fmpq_mat_swap_entrywise(fmpq_mat_t mat1, fmpq_mat_t mat2)

    Swaps two matrices by swapping the individual entries rather than swapping
    the contents of the structs.


Entry access
--------------------------------------------------------------------------------


.. function:: fmpq * fmpq_mat_entry(const fmpq_mat_t mat, slong i, slong j)

    Gives a reference to the entry at row ``i`` and column ``j``.
    The reference can be passed as an input or output variable to any
    ``fmpq`` function for direct manipulation of the matrix element.
    No bounds checking is performed.

.. function:: fmpz * fmpq_mat_entry_num(const fmpq_mat_t mat, slong i, slong j)

    Gives a reference to the numerator of the entry at row ``i`` and
    column ``j``. The reference can be passed as an input or output
    variable to any ``fmpz`` function for direct manipulation of the
    matrix element. No bounds checking is performed.

.. function:: fmpz * fmpq_mat_entry_den(const fmpq_mat_t mat, slong i, slong j)

    Gives a reference to the denominator of the entry at row ``i`` and
    column ``j``. The reference can be passed as an input or output
    variable to any ``fmpz`` function for direct manipulation of the
    matrix element. No bounds checking is performed.

.. function:: slong fmpq_mat_nrows(const fmpq_mat_t mat)

    Return the number of rows of the matrix ``mat``.

.. function:: slong fmpq_mat_ncols(const fmpq_mat_t mat)

    Return the number of columns of the matrix ``mat``.


Basic assignment
--------------------------------------------------------------------------------


.. function:: void fmpq_mat_set(fmpq_mat_t dest, const fmpq_mat_t src)

    Sets the entries in ``dest`` to the same values as in ``src``, 
    assuming the two matrices have the same dimensions.

.. function:: void fmpq_mat_zero(fmpq_mat_t mat)

    Sets ``mat`` to the zero matrix.

.. function:: void fmpq_mat_one(fmpq_mat_t mat)

    Let `m` be the minimum of the number of rows and columns 
    in the matrix ``mat``.  This function sets the first 
    `m \times m` block to the identity matrix, and the remaining 
    block to zero.

.. function:: void fmpq_mat_transpose(fmpq_mat_t rop, const fmpq_mat_t op)

    Sets the matrix ``rop`` to the transpose of the matrix ``op``, 
    assuming that their dimensions are compatible.

.. function:: void fmpq_mat_swap_rows(fmpq_mat_t mat, slong * perm, slong r, slong s)
    
    Swaps rows ``r`` and ``s`` of ``mat``.  If ``perm`` is non-``NULL``, the
    permutation of the rows will also be applied to ``perm``.

.. function:: void fmpq_mat_swap_cols(fmpq_mat_t mat, slong * perm, slong r, slong s)
    
    Swaps columns ``r`` and ``s`` of ``mat``.  If ``perm`` is non-``NULL``, the
    permutation of the columns will also be applied to ``perm``.

.. function:: void fmpq_mat_invert_rows(fmpq_mat_t mat, slong * perm)
    
    Swaps rows ``i`` and ``r - i`` of ``mat`` for ``0 <= i < r/2``, where
    ``r`` is the number of rows of ``mat``. If ``perm`` is non-``NULL``, the
    permutation of the rows will also be applied to ``perm``.

.. function:: void fmpq_mat_invert_cols(fmpq_mat_t mat, slong * perm)
    
    Swaps columns ``i`` and ``c - i`` of ``mat`` for ``0 <= i < c/2``, where
    ``c`` is the number of columns of ``mat``. If ``perm`` is non-``NULL``, the
    permutation of the columns will also be applied to ``perm``.

Addition, scalar multiplication
--------------------------------------------------------------------------------


.. function:: void fmpq_mat_add(fmpq_mat_t mat, const fmpq_mat_t mat1, const fmpq_mat_t mat2)

    Sets ``mat`` to the sum of ``mat1`` and ``mat2``, 
    assuming that all three matrices have the same dimensions.

.. function:: void fmpq_mat_sub(fmpq_mat_t mat, const fmpq_mat_t mat1, const fmpq_mat_t mat2)

    Sets ``mat`` to the difference of ``mat1`` and ``mat2``, 
    assuming that all three matrices have the same dimensions.

.. function:: void fmpq_mat_neg(fmpq_mat_t rop, const fmpq_mat_t op)

    Sets ``rop`` to the negative of ``op``, assuming that 
    the two matrices have the same dimensions.

.. function:: void fmpq_mat_scalar_mul_fmpq(fmpq_mat_t rop, const fmpq_mat_t op, const fmpq_t x)

    Sets ``rop`` to ``op`` multiplied by the rational `x`, 
    assuming that the two matrices have the same dimensions.

    Note that the rational ``x`` may not be aliased with any part of the
    entries of ``rop``.

.. function:: void fmpq_mat_scalar_mul_fmpz(fmpq_mat_t rop, const fmpq_mat_t op, const fmpz_t x)

    Sets ``rop`` to ``op`` multiplied by the integer `x`, 
    assuming that the two matrices have the same dimensions.

    Note that the integer `x` may not be aliased with any part of 
    the entries of ``rop``.

.. function:: void fmpq_mat_scalar_div_fmpz(fmpq_mat_t rop, const fmpq_mat_t op, const fmpz_t x)

    Sets ``rop`` to ``op`` divided by the integer `x`, 
    assuming that the two matrices have the same dimensions 
    and that `x` is non-zero.

    Note that the integer `x` may not be aliased with any part of 
    the entries of ``rop``.


Input and output
--------------------------------------------------------------------------------


.. function:: void fmpq_mat_print(const fmpq_mat_t mat)

    Prints the matrix ``mat`` to standard output.


Random matrix generation
--------------------------------------------------------------------------------


.. function:: void fmpq_mat_randbits(fmpq_mat_t mat, flint_rand_t state, flint_bitcnt_t bits)

    This is equivalent to applying ``fmpq_randbits`` to all entries
    in the matrix.

.. function:: void fmpq_mat_randtest(fmpq_mat_t mat, flint_rand_t state, flint_bitcnt_t bits)

    This is equivalent to applying ``fmpq_randtest`` to all entries
    in the matrix.


Window
--------------------------------------------------------------------------------


.. function:: void fmpq_mat_window_init(fmpq_mat_t window, const fmpq_mat_t mat, slong r1, slong c1, slong r2, slong c2)

    Initializes the matrix ``window`` to be an ``r2 - r1`` by
    ``c2 - c1`` submatrix of ``mat`` whose ``(0,0)`` entry
    is the ``(r1, c1)`` entry of ``mat``. The memory for the
    elements of ``window`` is shared with ``mat``.

.. function:: void fmpq_mat_window_clear(fmpq_mat_t window)

    Clears the matrix ``window`` and releases any memory that it
    uses. Note that the memory to the underlying matrix that
    ``window`` points to is not freed.


Concatenate
--------------------------------------------------------------------------------


.. function:: void fmpq_mat_concat_vertical(fmpq_mat_t res, const fmpq_mat_t mat1, const fmpq_mat_t mat2)

    Sets ``res`` to vertical concatenation of (``mat1``, ``mat2``) in that order. Matrix dimensions: ``mat1``: `m \times n`, ``mat2``: `k \times n`, ``res``: `(m + k) \times n`.

.. function:: void fmpq_mat_concat_horizontal(fmpq_mat_t res, const fmpq_mat_t mat1, const fmpq_mat_t mat2)

    Sets ``res`` to horizontal concatenation of (``mat1``, ``mat2``) in that order. Matrix dimensions: ``mat1``: `m \times n`, ``mat2``: `m \times k`, ``res``: `m \times (n + k)`.


Special matrices
--------------------------------------------------------------------------------


.. function:: void fmpq_mat_hilbert_matrix(fmpq_mat_t mat)

    Sets ``mat`` to a Hilbert matrix of the given size. That is,
    the entry at row `i` and column `j` is set to `1/(i+j+1)`.


Basic comparison and properties
--------------------------------------------------------------------------------


.. function:: int fmpq_mat_equal(const fmpq_mat_t mat1, const fmpq_mat_t mat2)

    Returns nonzero if ``mat1`` and ``mat2`` have the same shape and
    all their entries agree, and returns zero otherwise. Assumes the
    entries in both ``mat1`` and ``mat2`` are in canonical form.

.. function:: int fmpq_mat_is_integral(const fmpq_mat_t mat)

    Returns nonzero if all entries in ``mat`` are integer-valued, and
    returns zero otherwise. Assumes that the entries in ``mat``
    are in canonical form.

.. function:: int fmpq_mat_is_zero(const fmpq_mat_t mat)

    Returns nonzero if all entries in ``mat`` are zero, and returns
    zero otherwise.

.. function:: int fmpq_mat_is_one(const fmpq_mat_t mat)

    Returns nonzero if ``mat`` ones along the diagonal and zeros elsewhere,
    and returns zero otherwise.

.. function:: int fmpq_mat_is_empty(const fmpq_mat_t mat)

    Returns a non-zero value if the number of rows or the number of
    columns in ``mat`` is zero, and otherwise returns
    zero.

.. function:: int fmpq_mat_is_square(const fmpq_mat_t mat)

    Returns a non-zero value if the number of rows is equal to the
    number of columns in ``mat``, and otherwise returns zero.



Integer matrix conversion
--------------------------------------------------------------------------------


.. function:: int fmpq_mat_get_fmpz_mat(fmpz_mat_t dest, const fmpq_mat_t mat)

    Sets ``dest`` to ``mat`` and returns nonzero if all entries
    in ``mat`` are integer-valued. If not all entries in ``mat``
    are integer-valued, sets ``dest`` to an undefined matrix
    and returns zero. Assumes that the entries in ``mat`` are
    in canonical form.

.. function:: void fmpq_mat_get_fmpz_mat_entrywise(fmpz_mat_t num, fmpz_mat_t den, const fmpq_mat_t mat)

    Sets the integer matrices ``num`` and ``den`` respectively
    to the numerators and denominators of the entries in ``mat``.

.. function:: void fmpq_mat_get_fmpz_mat_matwise(fmpz_mat_t num, fmpz_t den, const fmpq_mat_t mat)

    Converts all entries in ``mat`` to a common denominator,
    storing the rescaled numerators in ``num`` and the
    denominator in ``den``. The denominator will be minimal
    if the entries in ``mat`` are in canonical form.

.. function:: void fmpq_mat_get_fmpz_mat_rowwise(fmpz_mat_t num, fmpz * den, const fmpq_mat_t mat)

    Clears denominators in ``mat`` row by row. The rescaled
    numerators are written to ``num``, and the denominator
    of row ``i`` is written to position ``i`` in ``den``
    which can be a preinitialised ``fmpz`` vector. Alternatively,
    ``NULL`` can be passed as the ``den`` variable, in which
    case the denominators will not be stored.

.. function:: void fmpq_mat_get_fmpz_mat_rowwise_2(fmpz_mat_t num, fmpz_mat_t num2, fmpz * den, const fmpq_mat_t mat, const fmpq_mat_t mat2)

    Clears denominators row by row of both ``mat`` and ``mat2``,
    writing the respective numerators to ``num`` and ``num2``.
    This is equivalent to concatenating ``mat`` and ``mat2``
    horizontally, calling ``fmpq_mat_get_fmpz_mat_rowwise``,
    and extracting the two submatrices in the result.

.. function:: void fmpq_mat_get_fmpz_mat_colwise(fmpz_mat_t num, fmpz * den, const fmpq_mat_t mat)

    Clears denominators in ``mat`` column by column. The rescaled
    numerators are written to ``num``, and the denominator
    of column ``i`` is written to position ``i`` in ``den``
    which can be a preinitialised ``fmpz`` vector. Alternatively,
    ``NULL`` can be passed as the ``den`` variable, in which
    case the denominators will not be stored.

.. function:: void fmpq_mat_set_fmpz_mat(fmpq_mat_t dest, const fmpz_mat_t src)

    Sets ``dest`` to ``src``.

.. function:: void fmpq_mat_set_fmpz_mat_div_fmpz(fmpq_mat_t mat, const fmpz_mat_t num, const fmpz_t den)

    Sets ``mat`` to the integer matrix ``num`` divided by the
    common denominator ``den``.


Modular reduction and rational reconstruction
--------------------------------------------------------------------------------


.. function:: void fmpq_mat_get_fmpz_mat_mod_fmpz(fmpz_mat_t dest, const fmpq_mat_t mat, const fmpz_t mod)

    Sets each entry in ``dest`` to the corresponding entry in ``mat``,
    reduced modulo ``mod``.

.. function:: int fmpq_mat_set_fmpz_mat_mod_fmpz(fmpq_mat_t X, const fmpz_mat_t Xmod, const fmpz_t mod)

    Sets ``X`` to the entrywise rational reconstruction integer matrix
    ``Xmod`` modulo ``mod``, and returns nonzero if the reconstruction
    is successful. If rational reconstruction fails for any element,
    returns zero and sets the entries in ``X`` to undefined values.


Matrix multiplication
--------------------------------------------------------------------------------


.. function:: void fmpq_mat_mul_direct(fmpq_mat_t C, const fmpq_mat_t A, const fmpq_mat_t B)

    Sets ``C`` to the matrix product ``AB``, computed
    naively using rational arithmetic. This is typically very slow and
    should only be used in circumstances where clearing denominators
    would consume too much memory.

.. function:: void fmpq_mat_mul_cleared(fmpq_mat_t C, const fmpq_mat_t A, const fmpq_mat_t B)

    Sets ``C`` to the matrix product ``AB``, computed
    by clearing denominators and multiplying over the integers.

.. function:: void fmpq_mat_mul(fmpq_mat_t C, const fmpq_mat_t A, const fmpq_mat_t B)

    Sets ``C`` to the matrix product ``AB``. This
    simply calls ``fmpq_mat_mul_cleared``.

.. function:: void fmpq_mat_mul_fmpz_mat(fmpq_mat_t C, const fmpq_mat_t A, const fmpz_mat_t B)

    Sets ``C`` to the matrix product ``AB``, with ``B``
    an integer matrix. This function works efficiently by clearing
    denominators of ``A``.

.. function:: void fmpq_mat_mul_r_fmpz_mat(fmpq_mat_t C, const fmpz_mat_t A, const fmpq_mat_t B)

    Sets ``C`` to the matrix product ``AB``, with ``A``
    an integer matrix. This function works efficiently by clearing
    denominators of ``B``.

.. function:: void fmpq_mat_mul_fmpq_vec(fmpq * c, const fmpq_mat_t A, const fmpq * b, slong blen)
              void fmpq_mat_mul_fmpz_vec(fmpq * c, const fmpq_mat_t A, const fmpz * b, slong blen)
              void fmpq_mat_mul_fmpq_vec_ptr(fmpq * const * c, const fmpq_mat_t A, const fmpq * const * b, slong blen)
              void fmpq_mat_mul_fmpz_vec_ptr(fmpq * const * c, const fmpq_mat_t A, const fmpz * const * b, slong blen)

    Compute a matrix-vector product of ``A`` and ``(b, blen)`` and store the result in ``c``.
    The vector ``(b, blen)`` is either truncated or zero-extended to the number of columns of ``A``.
    The number entries written to ``c`` is always equal to the number of rows of ``A``.

.. function:: void fmpq_mat_fmpq_vec_mul(fmpq * c, const fmpq * a, slong alen, const fmpq_mat_t B)
              void fmpq_mat_fmpz_vec_mul(fmpq * c, const fmpz * a, slong alen, const fmpq_mat_t B)
              void fmpq_mat_fmpq_vec_mul_ptr(fmpq * const * c, const fmpq * const * a, slong alen, const fmpq_mat_t B)
              void fmpq_mat_fmpz_vec_mul_ptr(fmpq * const * c, const fmpz * const * a, slong alen, const fmpq_mat_t B)

    Compute a vector-matrix product of ``(a, alen)`` and ``B`` and and store the result in ``c``.
    The vector ``(a, alen)`` is either truncated or zero-extended to the number of rows of ``B``.
    The number entries written to ``c`` is always equal to the number of columns of ``B``.


Kronecker product
--------------------------------------------------------------------------------


.. function:: void fmpq_mat_kronecker_product(fmpq_mat_t C, const fmpq_mat_t A, const fmpq_mat_t B)

    Sets ``C`` to the Kronecker product of ``A`` and ``B``.


Trace
--------------------------------------------------------------------------------


.. function:: void fmpq_mat_trace(fmpq_t trace, const fmpq_mat_t mat)

    Computes the trace of the matrix, i.e. the sum of the entries on
    the main diagonal. The matrix is required to be square.


Determinant
--------------------------------------------------------------------------------


.. function:: void fmpq_mat_det(fmpq_t det, const fmpq_mat_t mat)

    Sets ``det`` to the determinant of ``mat``. In the general case,
    the determinant is computed by clearing denominators and computing a
    determinant over the integers. Matrices of size 0, 1 or 2 are handled
    directly.


Nonsingular solving
--------------------------------------------------------------------------------


.. function:: int fmpq_mat_solve_fraction_free(fmpq_mat_t X, const fmpq_mat_t A, const fmpq_mat_t B)
              int fmpq_mat_solve_dixon(fmpq_mat_t X, const fmpq_mat_t A, const fmpq_mat_t B)
              int fmpq_mat_solve_multi_mod(fmpq_mat_t X, const fmpq_mat_t A, const fmpq_mat_t B)
              int fmpq_mat_solve(fmpq_mat_t X, const fmpq_mat_t A, const fmpq_mat_t B)

    Solves ``AX = B`` for nonsingular ``A``.
    Returns nonzero if ``A`` is nonsingular or if the right hand side
    is empty, and zero otherwise.

    All algorithms clear denominators to obtain a rescaled system over the integers.
    The *fraction_free* algorithm uses FFLU solving over the integers.
    The *dixon* and *multi_mod* algorithms use Dixon p-adic lifting
    or multimodular solving, followed by rational reconstruction
    with an adaptive stopping test. The *dixon* and *multi_mod* algorithms
    are generally the best choice for large systems.

    The default method chooses an algorithm automatically.

.. function:: int fmpq_mat_solve_fmpz_mat_fraction_free(fmpq_mat_t X, const fmpz_mat_t A, const fmpz_mat_t B)
              int fmpq_mat_solve_fmpz_mat_dixon(fmpq_mat_t X, const fmpz_mat_t A, const fmpz_mat_t B)
              int fmpq_mat_solve_fmpz_mat_multi_mod(fmpq_mat_t X, const fmpz_mat_t A, const fmpz_mat_t B)
              int fmpq_mat_solve_fmpz_mat(fmpq_mat_t X, const fmpz_mat_t A, const fmpz_mat_t B)

    Solves ``AX = B`` for nonsingular ``A``, where *A* and *B* are integer
    matrices. Returns nonzero if ``A`` is nonsingular or if the right hand side
    is empty, and zero otherwise.


.. function:: int fmpq_mat_can_solve_multi_mod(fmpq_mat_t X, const fmpq_mat_t A, const fmpq_mat_t B)

    Returns `1` if ``AX = B`` has a solution and if so, sets ``X`` to one such
    solution. The matrices can have any shape but must have the same number of
    rows.

.. function:: int fmpq_mat_can_solve_fraction_free(fmpq_mat_t X, const fmpq_mat_t A, const fmpq_mat_t B)

    Returns `1` if ``AX = B`` has a solution and if so, sets ``X`` to one such
    solution. The matrices can have any shape but must have the same number of
    rows.

int fmpq_mat_can_solve_fmpz_mat_dixon(fmpq_mat_t X, const fmpz_mat_t A, const fmpz_mat_t B)

    Returns `1` if ``AX = B`` has a solution and if so, sets ``X`` to one such
    solution. The matrices can have any shape but must have the same number of
    rows. The input matrices must have integer entries and `A` cannot be an
    empty matrix.

int fmpq_mat_can_solve_dixon(fmpq_mat_t X, const fmpq_mat_t A, const fmpq_mat_t B)

    Returns `1` if ``AX = B`` has a solution and if so, sets ``X`` to one such
    solution. The matrices can have any shape but must have the same number of
    rows.

.. function:: int fmpq_mat_can_solve(fmpq_mat_t X, const fmpq_mat_t A, const fmpq_mat_t B)

    Returns `1` if ``AX = B`` has a solution and if so, sets ``X`` to one such
    solution. The matrices can have any shape but must have the same number of
    rows.


Inverse
--------------------------------------------------------------------------------


.. function:: int fmpq_mat_inv(fmpq_mat_t B, const fmpq_mat_t A)

    Sets ``B`` to the inverse matrix of ``A`` and returns nonzero.
    Returns zero if ``A`` is singular. ``A`` must be a square matrix.



Echelon form
--------------------------------------------------------------------------------


.. function:: int fmpq_mat_pivot(slong * perm, fmpq_mat_t mat, slong r, slong c)

    Helper function for row reduction. Returns 1 if the entry of ``mat``
    at row `r` and column `c` is nonzero. Otherwise searches for a nonzero
    entry in the same column among rows `r+1, r+2, \ldots`. If a nonzero
    entry is found at row `s`, swaps rows `r` and `s` and the corresponding
    entries in ``perm`` (unless ``NULL``) and returns -1. If no
    nonzero pivot entry is found, leaves the inputs unchanged and returns 0.

.. function:: slong fmpq_mat_rref_classical(fmpq_mat_t B, const fmpq_mat_t A)

    Sets ``B`` to the reduced row echelon form of ``A`` and returns
    the rank. Performs Gauss-Jordan elimination directly over the rational
    numbers. This algorithm is usually inefficient and is mainly intended
    to be used for testing purposes.

.. function:: slong fmpq_mat_rref_fraction_free(fmpq_mat_t B, const fmpq_mat_t A)

    Sets ``B`` to the reduced row echelon form of ``A`` and returns
    the rank. Clears denominators and performs fraction-free Gauss-Jordan
    elimination using ``fmpz_mat`` functions.

.. function:: slong fmpq_mat_rref(fmpq_mat_t B, const fmpq_mat_t A)

    Sets ``B`` to the reduced row echelon form of ``A`` and returns
    the rank. This function automatically chooses between the classical and
    fraction-free algorithms depending on the size of the matrix.


Gram-Schmidt Orthogonalisation
--------------------------------------------------------------------------------


.. function:: void fmpq_mat_gso(fmpq_mat_t B, const fmpq_mat_t A)

    Takes a subset of `\mathbb{Q}^m` `S = \{a_1, a_2, \ldots ,a_n\}` (as the
    columns of a `m \times n` matrix ``A``) and generates an orthogonal set
    `S' = \{b_1, b_2, \ldots ,b_n\}` (as the columns of the `m \times n` matrix 
    ``B``) that spans the same subspace of `\mathbb{Q}^m` as `S`.


Transforms
--------------------------------------------------------------------------------


.. function:: void fmpq_mat_similarity(fmpq_mat_t A, slong r, fmpq_t d)

    Applies a similarity transform to the `n\times n` matrix `M` in-place.

    If `P` is the `n\times n` identity matrix the zero entries of whose row
    `r` (`0`-indexed) have been replaced by `d`, this transform is equivalent
    to `M = P^{-1}MP`.

    Similarity transforms preserve the determinant, characteristic polynomial
    and minimal polynomial.


Characteristic polynomial
--------------------------------------------------------------------------------


.. function:: void _fmpq_mat_charpoly(fmpz * coeffs, fmpz_t den, const fmpq_mat_t mat)

    Set ``(coeffs, den)`` to the characteristic polynomial of the given
    `n\times n` matrix.

.. function:: void fmpq_mat_charpoly(fmpq_poly_t pol, const fmpq_mat_t mat)

    Set ``pol`` to the characteristic polynomial of the given `n\times n`
    matrix. If ``mat`` is not square, an exception is raised.


Minimal polynomial
--------------------------------------------------------------------------------


.. function:: slong _fmpq_mat_minpoly(fmpz * coeffs, fmpz_t den, const fmpq_mat_t mat)

    Set ``(coeffs, den)`` to the minimal polynomial of the given
    `n\times n` matrix and return the length of the polynomial.

.. function:: void fmpq_mat_minpoly(fmpq_poly_t pol, const fmpq_mat_t mat)

    Set ``pol`` to the minimal polynomial of the given `n\times n`
    matrix. If ``mat`` is not square, an exception is raised.
