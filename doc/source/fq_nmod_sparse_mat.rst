.. _fq-nmod-sparse-mat:

**fq_nmod_sparse_mat.h** -- sparse matrixs over finite fields (word-size characteristic)
===============================================================================

Description.

Types, macros and constants
-------------------------------------------------------------------------------

.. type:: fq_nmod_sparse_mat_t

    Holds an array of (possibly empty) sparse vectors corresponding to rows in 
    the matrix

Memory management
--------------------------------------------------------------------------------

.. function:: void fq_nmod_sparse_mat_init(fq_nmod_sparse_mat_t M, slong rows, slong cols, const fq_nmod_ctx_t ctx)

    Initializes an empty sparse matrix ``M`` with given number of rows and columns

.. function:: void fq_nmod_sparse_mat_clear(fq_nmod_sparse_mat_t M, const fq_nmod_ctx_t ctx)

    Clears the entries of the matrix ``M`` and frees the space allocated for it

.. function:: void fq_nmod_sparse_mat_swap(fq_nmod_sparse_mat_t M1, fq_nmod_sparse_mat_t M2, const fq_nmod_ctx_t ctx)

    Swaps two matrices ``M1`` and ``M2`` (no reallocaton)


Instantiation
--------------------------------------------------------------------------------

.. function:: void fq_nmod_sparse_mat_zero(fq_nmod_sparse_mat_t M, const fq_nmod_ctx_t ctx)

    Sets matrix ``M`` to zero (the empty sparse matrix)

.. function:: void fq_nmod_sparse_mat_one(fq_nmod_sparse_mat_t M, const fq_nmod_ctx_t ctx)

    Sets matrix ``M`` to identity matrix (based on its number of rows)

.. function:: void fq_nmod_sparse_mat_set(fq_nmod_sparse_mat_t N, const fq_nmod_sparse_mat_t M, const fq_nmod_ctx_t ctx)

    Makes ``N`` a (deep) copy of ``M``

.. function:: void fq_nmod_sparse_mat_from_entries(fq_nmod_sparse_mat_t M, slong *rows, slong *inds, fq_nmod_struct *vals, slong nnz, const fq_nmod_ctx_t ctx)

    Constructs matrix ``M`` from a given sequence of ``rows``, ``cols``, and 
    corresponding ``vals``, all of length ``nnz``, assumes sorted by rows 
    with no duplicate (row, col) indices

.. function:: void fq_nmod_sparse_mat_append_col(fq_nmod_sparse_mat_t M, const fq_nmod_struct *v, const fq_nmod_ctx_t ctx)

    Add a dense column to the right of the matrix

.. function:: void fq_nmod_sparse_mat_append_row(fq_nmod_sparse_mat_t M, const fq_nmod_sparse_vec_t v, const fq_nmod_ctx_t ctx)

    Add a sparse row to the bottom of the matrix


Conversion to/from dense matrix
--------------------------------------------------------------------------------

.. function:: void fq_nmod_sparse_mat_from_dense(fq_nmod_sparse_mat_t M, const fq_nmod_mat_t dM, const fq_nmod_ctx_t ctx)

    Converts the dense matrix ``dM`` to a sparse matrix ``M``

.. function:: void fq_nmod_sparse_mat_to_dense(fq_nmod_mat_t dM, const fq_nmod_sparse_mat_t M, const fq_nmod_ctx_t ctx)

    Converts the sparse matrix ``M`` to a dense matrix ``dM``

Windows, concatenation, and splitting
--------------------------------------------------------------------------------

.. function:: void fq_nmod_sparse_mat_window_init(fq_nmod_sparse_mat_t window, const fq_nmod_sparse_mat_t M, slong r1, slong c1, slong r2, slong c2, const fq_nmod_ctx_t ctx)

    Constructs a window on the sparse matrix ``M`` between rows ``r1`` and ``r2`` 
    and cols ``c1`` and ``c2`` (valid as long as original matrix remains unmodified)

.. function:: void fq_nmod_sparse_mat_window_clear(fq_nmod_sparse_mat_t window, const fq_nmod_ctx_t ctx)

    Clears a window

.. function:: void fq_nmod_sparse_mat_concat_horizontal(fq_nmod_sparse_mat_t B, const fq_nmod_sparse_mat_t M1, const fq_nmod_sparse_mat_t M2, const fq_nmod_ctx_t ctx)

    Horizontally concatenates two matrices ``M1`` and ``M2`` into block matrix ``B``

.. function:: void fq_nmod_sparse_mat_concat_vertical(fq_nmod_sparse_mat_t B, const fq_nmod_sparse_mat_t M1, const fq_nmod_sparse_mat_t M2, const fq_nmod_ctx_t ctx)

    Vertically concatenates two matrices ``M1`` and ``M2`` into block matrix ``B``

.. function:: void fq_nmod_sparse_mat_split_horizontal(fq_nmod_sparse_mat_t M1, fq_nmod_sparse_mat_t M2, const fq_nmod_sparse_mat_t B, slong c, const fq_nmod_ctx_t ctx)

    Splits ``B`` horizontally into two submatrices ``M1`` and ``M2``, dividing at column ``c``

.. function:: void fq_nmod_sparse_mat_split_vertical(fq_nmod_sparse_mat_t M1, fq_nmod_sparse_mat_t M2, const fq_nmod_sparse_mat_t B, slong r, const fq_nmod_ctx_t ctx)

    Splits ``B`` vertically into two submatrices ``M1`` and ``M2``, dividing at row ``r``


Permutation
--------------------------------------------------------------------------------

.. function:: void fq_nmod_sparse_mat_permute_cols(fq_nmod_sparse_mat_t M, slong *Q, const fq_nmod_ctx_t ctx)

    Permutes the columns indices of ``M`` according to ``Q``, and re-sorts each row

.. function:: void fq_nmod_sparse_mat_permute_rows(fq_nmod_sparse_mat_t M, slong *P, const fq_nmod_ctx_t ctx)

    Permutes the rows of ``M`` according to ``P``


Randomization
--------------------------------------------------------------------------------


.. function:: void fq_nmod_sparse_mat_randtest(fq_nmod_sparse_mat_t M, flint_rand_t state, slong min_nnz, slong max_nnz, const fq_nmod_ctx_t ctx)

    Makes ``M`` a sparse matrix with between ``min_nnz`` and ``max_nnz`` nonzero 
    entries per row, with individual entries generated by fq_nmod_randtest


Output
--------------------------------------------------------------------------------

.. function:: void fq_nmod_sparse_mat_print_pretty(const fq_nmod_sparse_mat_t M, const fq_nmod_ctx_t ctx)

    Prints the matrix ``M`` to ``stdout`` in a human-readable format


Comparison
--------------------------------------------------------------------------------

.. function:: void fq_nmod_sparse_is_zero(fq_nmod_sparse_mat_t M, const fq_nmod_ctx_t ctx)

    Checks if the given matrix ``M`` is trivial (empty), returning `1` if so and `0` 
    otherwise

.. function:: void fq_nmod_sparse_mat_equal(const fq_nmod_sparse_mat_t M1, const fq_nmod_sparse_mat_t M2, slong ioff, const fq_nmod_ctx_t ctx)

    Checks if ``M1`` equals ``M2``, returning `1` if so and `0` otherwise


Transpose
--------------------------------------------------------------------------------

.. function:: void fq_sparse_mat_transpose(fq_nmod_sparse_mat_t N, const fq_nmod_sparse_mat_t M, const fq_nmod_ctx_t ctx)

    Transposes ``M`` into the matrix ``N`` (must have swapped rows and columns)


Arithmetic
--------------------------------------------------------------------------------

.. function:: void fq_nmod_sparse_mat_neg(fq_nmod_sparse_mat_t N, const fq_nmod_sparse_mat_t M, const fq_nmod_ctx_t ctx)

    Sets ``N`` to the negation of ``M``

.. function:: void fq_nmod_sparse_mat_scalar_mul_fq_nmod(fq_nmod_sparse_mat_t N, const fq_nmod_sparse_mat_t M, const fq_nmod_t c, const fq_nmod_ctx_t ctx)

    Sets ``N`` to the scalar multiple of ``M`` by ``c``

.. function:: void fq_nmod_sparse_mat_add(fq_nmod_sparse_mat_t O, const fq_nmod_sparse_mat_t M, const fq_nmod_sparse_mat_t N, const fq_nmod_ctx_t ctx)

    Sets ``O`` to the sum of ``M`` and ``N``

.. function:: void fq_nmod_sparse_mat_sub(fq_nmod_sparse_mat_t O, const fq_nmod_sparse_mat_t M, const fq_nmod_sparse_mat_t N, const fq_nmod_ctx_t ctx)

    Sets ``O`` to the difference of ``M`` and ``N``

.. function:: void fq_nmod_sparse_mat_scalar_addmul_fq_nmod(fq_nmod_sparse_mat_t O, const fq_nmod_sparse_mat_t M, const fq_nmod_sparse_mat_t N, const fq_nmod_t c, const fq_nmod_ctx_t ctx)

    Sets ``O`` to the sum of ``M`` and ``c` times ``N``

.. function:: void fq_nmod_sparse_mat_scalar_submul_fq_nmod(fq_nmod_sparse_mat_t O, const fq_nmod_sparse_mat_t M, const fq_nmod_sparse_mat_t N, const fq_nmod_t c, const fq_nmod_ctx_t ctx)

    Sets ``O`` to the difference of ``M`` and ``N` times ``v``

.. function:: void fq_nmod_sparse_mat_mul_vec(fq_nmod_struct *y, const fq_nmod_sparse_mat_t M, const fq_nmod_struct *x, const fq_nmod_ctx_t ctx)

    Sets ``y`` to the product of ``M`` and ``x``

.. function:: void fq_nmod_sparse_mat_mul_mat(fq_nmod_mat_t Y, const fq_nmod_sparse_mat_t M, const fq_nmod_mat_t X, const fq_nmod_ctx_t ctx)

    Sets ``Y`` to the product of ``M`` and ``X``

.. function:: slong fq_nmod_sparse_mat_inv(fq_nmod_sparse_mat_t N, const fq_nmod_sparse_mat_t M, const fq_nmod_ctx_t ctx)

    Sets ``N`` to the "inverse" of ``M``, i.e., the matrix such that NM is
    in reduced row-echelon form


Decomposition/reduction
--------------------------------------------------------------------------------

.. function:: slong fq_nmod_sparse_mat_lu(slong *P, slong *Q, fq_nmod_sparse_mat_t L, fq_nmod_sparse_mat_t U, const fq_nmod_sparse_mat_t M, const fq_nmod_ctx_t ctx)

    Computes the decomposition PMQ = LU for a given sparse matrix ``M``, where 
    ``P`` is a row permutation, ``Q`` is a column permutation, ``L``is a lower
    triangular matrix, and ``U`` is an upper triangular matrix

.. function:: void fq_nmod_sparse_mat_rref(fq_nmod_sparse_mat_t M, const fq_nmod_ctx_t ctx)

    Applies row reduction to put ``M`` in reduced row echelon form (in place)

Solving
--------------------------------------------------------------------------------

.. function:: int fq_nmod_sparse_mat_solve_lu(fq_nmod_struct *x, const fq_nmod_sparse_mat_t M, const fq_nmod_struct *b, const fq_nmod_ctx_t ctx)

    Given a matrix ``M`` and target vector ``b``, use LU decomposition to find
    a vector ``x`` such that Mx = b, returns `1` if successful and `0` if not

.. function:: int fq_nmod_sparse_mat_solve_rref(fq_nmod_struct *x, const fq_nmod_sparse_mat_t M, const fq_nmod_struct *b, const fq_nmod_ctx_t ctx)

    Given a matrix ``M`` and target vector ``b``, use the reduced row-echelon
    form to find a vector ``x`` such that Mx = b, returns `1` if successful and 
    `0` if not

.. function:: int fq_nmod_sparse_mat_solve_lanczos(fq_nmod_struct *x, const fq_nmod_sparse_mat_t M, const fq_nmod_struct *b, flint_rand_t state, const fq_nmod_ctx_t ctx)

    Given a matrix ``M`` and target vector ``b``, use the Lanczos algorithm to
    find a vector ``x`` such that Mx = b, returns `1` if successful and `0` if not

.. function:: int fq_nmod_sparse_mat_solve_wiedemann(fq_nmod_struct *x, const fq_nmod_sparse_mat_t M, const fq_nmod_struct *b, const fq_nmod_ctx_t ctx)

    Given a matrix ``M`` and target vector ``b``, use the Wiedemann algorithm to
    find a vector ``x`` such that Mx = b, returns `1` if successful and `0` if not

.. function:: int fq_nmod_sparse_mat_solve_block_lanczos(fq_nmod_struct *x, const fq_nmod_sparse_mat_t M, const fq_nmod_struct *b, slong block_size, flint_rand_t state, const fq_nmod_ctx_t ctx)

    Given a matrix ``M`` and target vector ``b``, use Coppersmith's block Lanczos 
    algorithm (with specified block size) to find a vector ``x`` such that Mx = b, 
    returns `1` if successful and `0` if not

.. function:: int fq_nmod_sparse_mat_solve_block_wiedemann(fq_nmod_struct *x, const fq_nmod_sparse_mat_t M, const fq_nmod_struct *b, slong block_size, flint_rand_t state, const fq_nmod_ctx_t ctx)

    Given a matrix ``M`` and target vector ``b``, use Coppersmith's block Wiedemann
    algorithm (with specified block size) to find a vector ``x`` such that Mx = b, 
    returns `1` if successful and `0` if not

Nullvector and nullspace computation
--------------------------------------------------------------------------------

.. function:: int fq_nmod_sparse_mat_nullvector_lanczos(fq_nmod_struct *x, const fq_nmod_sparse_mat_t M, flint_rand_t state, const fq_nmod_ctx_t ctx)

    Given a matrix ``M``, use the Lanczos algorithm to find a nullvector ``x`` 
    s.t. Mx = 0, returns `1` if successful and `0` if not

.. function:: int fq_nmod_sparse_mat_nullvector_wiedemann(fq_nmod_struct *x, const fq_nmod_sparse_mat_t M, flint_rand_t state, const fq_nmod_ctx_t ctx)

    Given a matrix ``M``, use the Wiedemann algorithm to find a nullvector ``x`` 
    s.t. Mx = 0, returns `1` if successful and `0` if not

.. function:: int fq_nmod_sparse_mat_nullvector_block_lanczos(fq_nmod_struct *x, const fq_nmod_sparse_mat_t M, slong block_size, flint_rand_t state, const fq_nmod_ctx_t ctx)

    Given a matrix ``M``, use Coppersmith's block Lanczos algorithm to find a 
    nullvector ``x`` s.t. Mx = 0, returns `1` if successful and `0` if not

.. function:: int fq_nmod_sparse_mat_nullvector_block_wiedemann(fq_nmod_struct *x, const fq_nmod_sparse_mat_t M, slong block_size, flint_rand_t state, const fq_nmod_ctx_t ctx)

    Given a matrix ``M``, use Coppersmith's block Wiedemann algorithm to find a 
    nullvector ``x`` s.t. Mx = 0, returns `1` if successful and `0` if not

.. function:: int fq_nmod_sparse_mat_nullspace_rref(fq_nmod_mat_t X, const fq_nmod_sparse_mat_t M, const fq_nmod_ctx_t ctx)

    Given a matrix ``M``, use the reduced row echelon form to construct the 
    nullspace ``X`` of M (initialized by this function), returns the nullity

.. function:: int fq_nmod_sparse_mat_nullspace_lu(fq_nmod_mat_t X, const fq_nmod_sparse_mat_t M, const fq_nmod_ctx_t ctx)

    Given a matrix ``M``, use the LU decomposition to construct the nullspace ``X``
    of M (initialized by this function), returns the nullity

.. function:: int fq_nmod_sparse_mat_nullspace_lanczos(fq_nmod_mat_t X, const fq_nmod_sparse_mat_t M, flint_rand_t state, const fq_nmod_ctx_t ctx)

    Given a matrix ``M``, use the Lanczos algorithm to find a nullspace ``X`` 
    of M (initialized by this function), returns the found nullity

.. function:: int fq_nmod_sparse_mat_nullspace_wiedemann(fq_nmod_mat_t X, const fq_nmod_sparse_mat_t M, flint_rand_t state, const fq_nmod_ctx_t ctx)

    Given a matrix ``M``, use the Wiedemann algorithm to find a nullspace ``X`` 
    of M (initialized by this function), returns the found nullity

.. function:: int fq_nmod_sparse_mat_nullspace_block_lanczos(fq_nmod_mat_t X, const fq_nmod_sparse_mat_t M, slong block_size, flint_rand_t state, const fq_nmod_ctx_t ctx)

    Given a matrix ``M``, use Coppersmith's block Lanczos algorithm to find a 
    nullspace ``X`` of M (initialized by this function), returns the found nullity

.. function:: int fq_nmod_sparse_mat_nullspace_block_wiedemann(fq_nmod_mat_t X, const fq_nmod_sparse_mat_t M, slong block_size, flint_rand_t state, const fq_nmod_ctx_t ctx)

    Given a matrix ``M``, use Coppersmith's block Wiedemann algorithm to find a 
    nullspace ``X`` of M (initialized by this function), returns the found nullity

