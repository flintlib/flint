.. _gr-sparse-vec:

**gr_sparse_mat.h** -- Sparse matrices over generic rings
===============================================================================

The types :type:`gr_csr_mat_t` and :type:`gr_lil_mat_t` are representions
of matrices over a generic
ring which are optimized for the situation where there are many entries which
are zero (the matrix is *sparse*, or has *low density*).  Internally, the
matrix records only the positions of the nonzeros and their values.
Although technically these types and all its functions will operate correctly
with any density, a regular :type:`gr_mat_t` may be more efficient in the number
of nonzeros is not small.

Types and basic access
--------------------------------------------------------------------------------

.. type:: gr_csr_mat_struct

.. type:: gr_csr_mat_t
        
    This struct presents a sparse matrix in "compressed sparse row" form, and contains:

    * An `slong` value ``r``, the number of rows in the matrix.
    * An `slong` value ``c``, the number of columns in the matrix.
    * An `slong` value ``nnz``, the number of nonzeroes in the matrix.
    * An `slong` value ``alloc``, the maximum number of nonzeroes currently allocated.
    * A `ulong` array ``rows`` of length ``r+1`` providing row offsets.
    * A `ulong` array ``cols`` of length ``nnz`` providing column indices.
    * A :type:`gr_ptr` array ``nzs`` of length ``nnz`` providing nonzero values.

    For a given row of index ``row``, its nonzero columns and associated values are given
    by the components of ``cols`` and ``entries`` associated to
    indices between ``rows[row]`` and ``rows[row+1]``. 
    We always require:

    * The ``cols`` for each row are sorted into strictly increasing order.
    * The ``entries`` are nonzero (meaning, ``gr_is_zero(GR_ENTRY(vec->nzs, i, ctx->sizeof_elem), ctx)`` returns ``T_FALSE`` for ``0 <= i < vec->nnz``).
    * We have ``nnz <= alloc <= length``.
    
    A ``gr_csr_mat_t`` is defined as an array of length one of type
    ``gr_csr_mat_struct``, permitting a ``gr_csr_mat_t`` to
    be passed by reference. Note that ``gr_csr_mat_t`` is the more efficient
    way to represent a sparse matrix, but is less convenient for modification.

.. type:: gr_lil_mat_struct

.. type:: gr_lil_mat_t
        
    This struct presents a sparse matrix in "list of lists" form, and contains:

    * An `slong` value ``r``, the number of rows in the matrix.
    * An `slong` value ``c``, the number of columns in the matrix.
    * An `slong` value ``nnz``, the number of nonzeroes in the matrix.
    * A :type:`gr_sparse_vec_t` array ``rows`` of length ``r``, the (sparse) rows of the matrix.

    A ``gr_lil_mat_t`` is defined as an array of length one of type
    ``gr_lil_mat_struct``, permitting a ``gr_lil_mat_t`` to
    be passed by reference. Note that ``gr_lil_mat_t`` is a less efficient
    way to represent a sparse matrix, but is more convenient for modification.

.. type:: gr_coo_mat_struct

.. type:: gr_coo_mat_t
        
    This struct presents a sparse matrix in (coo)rdinate form, and contains:

    * An `slong` value ``r``, the number of rows in the matrix.
    * An `slong` value ``c``, the number of columns in the matrix.
    * An `slong` value ``nnz``, the number of nonzeroes in the matrix.
    * An `slong` value ``alloc``, the maximum number of nonzeroes currently allocated.
    * A :type:`gr_sparse_vec_t` array ``rows`` of length ``nnz``, the list of rows for each entry.
    * A :type:`gr_sparse_vec_t` array ``cols`` of length ``nnz``, the list of columns for each entry.
    * A :type:`gr_sparse_vec_t` array ``nzs`` of length ``nnz``, the list of columns for each entry.
    * A `truth_t` value `is_canonical`, with is `T_TRUE` if the matrix is in "canonical" form (sorted by rows and columns, with unique non-zero entries).

    A ``gr_coo_mat_t`` is defined as an array of length one of type
    ``gr_coo_mat_struct``, permitting a ``gr_coo_mat_t`` to
    be passed by reference. Note that ``gr_coo_mat_t`` is the least efficient
    way to represent a sparse matrix (and has minimal functionality), but
    is convenient for constructing from arbitrary lists of entries.

.. function:: void gr_sparse_mat_nrows(gr_csr_mat_t mat, gr_ctx_t ctx)
              void gr_sparse_mat_ncols(gr_csr_mat_t mat, gr_ctx_t ctx)
              void gr_sparse_mat_nnz(gr_csr_mat_t mat, gr_ctx_t ctx)
              void gr_sparse_mat_nrows(gr_csr_mat_t mat, gr_ctx_t ctx)
              void gr_sparse_mat_ncols(gr_csr_mat_t mat, gr_ctx_t ctx)
              void gr_sparse_mat_nnz(gr_csr_mat_t mat, gr_ctx_t ctx)

    Get the number of rows/columns/nonzeroes of the given sparse matrix,
    in any representation.

.. macro:: GR_CSR_MAT_COL(mat, i, j)
           GR_LIL_MAT_COL(mat, i, j)

    Get a pointer to the column of the *j*-th nonzero in the *i*-th row of *mat*.
    There is no bounds checking.

.. macro:: GR_CSR_MAT_ENTRY(mat, i, j, sz)
           GR_LIL_MAT_ENTRY(mat, i, j, sz)

    Get a pointer to the value of the *j*-th nonzero entry in the *i*-th row of *mat*,
    given the size of each entry (typically obtained as ``ctx->sizeof_elem``).
    There is no bounds checking.

.. macro:: GR_COO_MAT_ROW(mat, i)
           GR_LIL_MAT_COL(mat, i)
           GR_LIL_MAT_ENTRY(mat, i)

    Get a pointer to the row, column, and value of the *i*-th nonzero in the matrix, repectively.
    There is no bounds checking.

.. function:: ulong * gr_csr_mat_col_ptr(gr_csr_mat_t mat, slong i, slong j)
              const ulong * gr_csr_mat_col_srcptr(const gr_csr_mat_t mat, slong i, slong j)
              ulong * gr_lil_mat_col_ptr(gr_lil_mat_t mat, slong i, slong j)
              const ulong * gr_lil_mat_col_srcptr(const gr_lil_mat_t mat, slong i, slong j)

    Get a (const) pointer to the column of the *j*-th nonzero in the *i*-th row of *mat*.
    If the location is out of bounds, the function returns NULL.

.. function:: gr_ptr gr_csr_mat_entry_ptr(gr_csr_mat_t mat, slong i, slong j, gr_ctx_t ctx)
              gr_srcptr gr_csr_mat_entry_srcptr(gr_csr_mat_t mat, slong i, slong j, gr_ctx_t ctx)
              gr_ptr gr_lil_mat_entry_ptr(gr_lil_mat_t mat, slong i, slong j, gr_ctx_t ctx)
              gr_srcptr gr_lil_mat_entry_srcptr(gr_lil_mat_t mat, slong i, slong j, gr_ctx_t ctx)

    Get a (const) pointer to the *j*-th nonzero entry in the *i*-th row of *mat*.
    If the location is out of bounds, the function returns NULL.

.. function:: ulong * gr_coo_mat_row_ptr(gr_coo_mat_t mat, slong i)
              const ulong * gr_coo_mat_row_srcptr(const gr_coo_mat_t mat, slong i)
              ulong * gr_coo_mat_col_ptr(gr_coo_mat_t mat, slong i)
              const ulong * gr_coo_mat_col_srcptr(const gr_coo_mat_t mat, slong i)
              gr_ptr gr_coo_mat_entry_ptr(gr_coo_mat_t mat, slong i)
              gr_srcptr gr_coo_mat_entry_ptr(const gr_coo_mat_t mat, slong i)

.. function:: void gr_csr_mat_init(gr_csr_mat_t mat, slong rows, slong cols, gr_ctx_t ctx)
              void gr_lil_mat_init(gr_lil_mat_t mat, slong rows, slong cols, gr_ctx_t ctx)
              void gr_coo_mat_init(gr_coo_mat_t mat, slong rows, slong cols, gr_ctx_t ctx)

    Initializes *mat* to a *rows* x *cols* matrix with no nonzeros.

.. function:: void gr_csr_mat_clear(gr_csr_mat_t mat, gr_ctx_t ctx)
              void gr_lil_mat_clear(gr_lil_mat_t mat, gr_ctx_t ctx)
              void gr_coo_mat_clear(gr_coo_mat_t mat, gr_ctx_t ctx)

    Clears the matrix *mat*.

.. function:: gr_csr_mat_swap(gr_csr_mat_t mat1, gr_csr_mat_t mat2, gr_ctx_t ctx)
              gr_lil_mat_swap(gr_lil_mat_t mat1, gr_lil_mat_t mat2, gr_ctx_t ctx)
              gr_coo_mat_swap(gr_coo_mat_t mat1, gr_coo_mat_t mat2, gr_ctx_t ctx)

    Swap data underlying two sparse matrices (no allocation or copying).

.. function:: void gr_csr_mat_fit_nnz(gr_csr_mat_t mat, slong nnz, gr_ctx_t ctx)
              void gr_lil_mat_fit_nnz(gr_lil_mat_t mat, slong *nnz, gr_ctx_t ctx)
              void gr_coo_mat_fit_nnz(gr_coo_mat_t mat, slong nnz, gr_ctx_t ctx)

    Ensure that *mat* has enough storage to hold at least *nnz* nonzeros.  This does
    not change the dimensions of the matrix or the number of nonzeros stored.
    Note that, for matrices in *lil* form, *nnz* is a vector of sizes, one for each row.

.. function:: void gr_csr_mat_shrink_to_nnz(gr_csr_mat_t mat, gr_ctx_t ctx)
              void gr_lil_mat_shrink_to_nnz(gr_lil_mat_t mat, gr_ctx_t ctx)
              void gr_coo_mat_shrink_to_nnz(gr_coo_mat_t mat, gr_ctx_t ctx)

    Reallocate the storage in *mat* down the current number of nonzeros.

    Note that, for matrices in *lil* form, the operation is performed
    on each row.

.. function:: void gr_csr_mat_set_cols(gr_csr_mat_t mat, slong cols, gr_ctx_t ctx)
              void gr_lil_mat_set_cols(gr_lil_mat_t mat, slong cols, gr_ctx_t ctx)
              void gr_coo_mat_set_cols(gr_coo_mat_t mat, slong cols, gr_ctx_t ctx)

    Set the nominal number of columns of the matrix *mat* to *cols*.  If *cols* is smaller than
    the current number of columns of *mat*, any entries whose columns are at least *cols*
    are truncated. That is, the number of nonzeros can change (but the allocation does not).

.. function:: int gr_coo_mat_from_entries(gr_coo_mat_t mat, ulong *rows, ulong *cols, gr_srcptr entries, slong nnz, truth_t is_canonical, gr_ctx_t ctx)

    Construct a sparse matrix in coordinate form from three list of corresponding
    rows, columns, and (presumably nonzero) values. If *is_canonical* is set to *T_TRUE*,
    the rows and columns are assumed to sorted in row-column order and unique, and the
    values are assumed to be non-zero. Otherwise, they may be in any order with zeroes,
    and duplicate elements are assumed to add.

.. function:: truth_t gr_coo_mat_is_canonical(gr_coo_mat_t mat, gr_ctx_t ctx)

    Check if a matrix in coordinate form is in canonical form, with sorted and unique
    (row, column) pairs and no (known) zero entries.

.. function:: int gr_coo_mat_canonicalize(gr_coo_mat_t mat, gr_ctx_t ctx)

    Put a coordinate-form sparse matrix in canonical form.

.. function:: int gr_coo_mat_randtest(gr_coo_mat_t mat, slong nnz, int replacement, truth_t is_canonical, flint_rand_t state, gr_ctx_t ctx)

    Construct a random coordinate-form sparse matrix *mat* with (approximately) *nnz* nonzeroes,
    in coordinate form iff *is_canonical* is set to *T_TRUE*. If replacement is set, the 
    (row, column) indices are chosen randomly with replacement, so the actual number of nonzeroes
    may be slightly less. Otherwise, the matrix is guaranteed to have *nnz* identically distributed
    nonzeroes, using reservoir sampling. The latter is better for sampling high density matrices,
    the former for low density ones. In both cases, the nonzero values are chosen using ``gr_randtest_nonzero``.

.. function:: int gr_coo_mat_randtest_prob(gr_coo_mat_t mat, double prob, flint_rand_t state, gr_ctx_t ctx)

    Construct a random coordinate-form sparse matrix *mat*, with each entry set to a random nonzero
    value with probability *prob*.

Getting, setting and conversion
--------------------------------------------------------------------------------

.. function:: gr_ptr gr_csr_mat_find_entry(gr_csr_mat_t mat, slong row, slong col, gr_ctx_t ctx)
              gr_ptr gr_lil_mat_find_entry(gr_lil_mat_t mat, slong row, slong col, gr_ctx_t ctx)
              gr_ptr gr_coo_mat_find_entry(gr_coo_mat_t mat, slong row, slong col, gr_ctx_t ctx)
              int gr_csr_mat_get_entry(gr_ptr res, gr_csr_mat_t mat, slong row, slong col, gr_ctx_t ctx)
              int gr_lil_mat_get_entry(gr_ptr res, gr_lil_mat_t mat, slong row, slong col, gr_ctx_t ctx)
              int gr_coo_mat_get_entry(gr_ptr res, gr_coo_mat_t mat, slong row, slong col, gr_ctx_t ctx)

    The *find* functions look for an entry at position (*row*, *col*), returning it if found and NULL
    otherwise; the corresponding *get* functions instead copy any value found into *res*, and set
    *res* to zero if that position is not found. In the case of non-canonical form COO
    sparse matrices, the value is taken to be the first nonzero found at that positiion.
    
    Because of the way sparse matrices are represented, this is 
    * logarithmic time in the number of nonzeros in the specified row, for CSR and LIL matrices,
    * logarithmic time in the number of nonzeros in the matrix, for canonical form COO matrices, and
    * linear time in the number of nonzeroes in the matrix, for non-canonical form COO matrices.

.. function:: int gr_csr_mat_set_entry(gr_csr_mat_t mat, slong row, slong col, gr_srcptr entry, gr_ctx_t ctx)
              int gr_lil_mat_set_entry(gr_lil_mat_t mat, slong row, slong col, gr_srcptr entry, gr_ctx_t ctx)
              int gr_coo_mat_set_entry(gr_coo_mat_t mat, slong row, slong col, gr_srcptr entry, gr_ctx_t ctx)

    Set the the entry at location (*row*, *col*) to be *entry*.
    
    Because of the way sparse vectors are represented, it is not efficient to call this function
    repeatedly: for the list of lists representation, it is linear time in the number of nonzeros
    in the associated row; for the sparse compressed row representation, it is linear time in the
    overall number of nonzeroes. If possible, the entries to update should be batched up and
    given using `gr_coo_mat_from_entries`.

.. function:: int gr_csr_mat_zero(gr_csr_mat_t mat, gr_ctx_t ctx)
              int gr_lil_mat_zero(gr_lil_mat_t mat, gr_ctx_t ctx)
              int gr_coo_mat_zero(gr_coo_mat_t mat, gr_ctx_t ctx)

    Set *mat* to the zero matrix (by setting *nnz* to 0; no actual element values are changed).

.. function:: int gr_csr_mat_set(gr_csr_mat_t res, const gr_csr_mat_t mat, gr_ctx_t ctx)
              int gr_lil_mat_set(gr_lil_mat_t res, const gr_lil_mat_t mat, gr_ctx_t ctx)
              int gr_coo_mat_set(gr_coo_mat_t res, const gr_lil_mat_t mat, gr_ctx_t ctx)
              int gr_csr_mat_set_lil_mat(gr_csr_mat_t res, const gr_lil_mat_t mat, gr_ctx_t ctx)
              int gr_lil_mat_set_csr_mat(gr_lil_mat_t res, const gr_csr_mat_t mat, gr_ctx_t ctx)
              int gr_coo_mat_set_lil_mat(gr_coo_mat_t res, const gr_lil_mat_t mat, gr_ctx_t ctx)
              int gr_coo_mat_set_csr_mat(gr_coo_mat_t res, const gr_csr_mat_t mat, gr_ctx_t ctx)
              int gr_csr_mat_set_coo_mat(gr_csr_mat_t res, const gr_coo_mat_t mat, gr_ctx_t ctx)
              int gr_lil_mat_set_coo_mat(gr_lil_mat_t res, const gr_coo_mat_t mat, gr_ctx_t ctx)

    Set *res* to a copy of *mat*, possibly changing the type of representation.

.. function:: int gr_csr_mat_set_mat(gr_csr_mat_t res, gr_mat_t mat, gr_ctx_t ctx)
              int gr_lil_mat_set_mat(gr_lil_mat_t res, gr_mat_t mat, gr_ctx_t ctx)
              int gr_coo_mat_set_mat(gr_coo_mat_t res, gr_mat_t mat, gr_ctx_t ctx)

    Set *res* from the (nominally) dense matrix *mat*.
    
.. function:: int gr_mat_set_csr_mat(gr_mat_t res, const gr_csr_mat_t mat, gr_ctx_t ctx)
              int gr_mat_set_lil_mat(gr_mat_t res, const gr_lil_mat_t mat, gr_ctx_t ctx)
              int gr_mat_set_coo_mat(gr_mat_t res, const gr_coo_mat_t mat, gr_ctx_t ctx)

    Set a dense matrix *res* from the sparse matrix *mat*, for some representation.

.. function:: int gr_csr_mat_init_set(gr_csr_mat_t res, const gr_csr_mat_t mat, gr_ctx_t ctx)
              int gr_lil_mat_init_set(gr_lil_mat_t res, const gr_lil_mat_t mat, gr_ctx_t ctx)
              int gr_csr_mat_init_set_lil_mat(gr_csr_mat_t res, const gr_lil_mat_t mat, gr_ctx_t ctx)
              int gr_lil_mat_init_set_csr_mat(gr_lil_mat_t res, const gr_csr_mat_t mat, gr_ctx_t ctx)
              int gr_coo_mat_init_set_lil_mat(gr_coo_mat_t res, const gr_lil_mat_t mat, gr_ctx_t ctx)
              int gr_coo_mat_init_set_csr_mat(gr_coo_mat_t res, const gr_csr_mat_t mat, gr_ctx_t ctx)
              int gr_csr_mat_init_set_coo_mat(gr_csr_mat_t res, const gr_coo_mat_t mat, gr_ctx_t ctx)
              int gr_lil_mat_init_set_coo_mat(gr_lil_mat_t res, const gr_coo_mat_t mat, gr_ctx_t ctx)
              int gr_csr_mat_init_set_mat(gr_csr_mat_t res, const gr_mat_t mat, gr_ctx_t ctx)
              int gr_lil_mat_init_set_mat(gr_lil_mat_t res, const gr_mat_t mat, gr_ctx_t ctx)
              int gr_coo_mat_init_set_mat(gr_coo_mat_t res, const gr_mat_t mat, gr_ctx_t ctx)
              int gr_mat_init_set_csr_mat(gr_mat_t res, const gr_csr_mat_t mat, gr_ctx_t ctx)
              int gr_mat_init_set_lil_mat(gr_mat_t res, const gr_lil_mat_t mat, gr_ctx_t ctx)
              int gr_mat_init_set_coo_mat(gr_mat_t res, const gr_coo_mat_t mat, gr_ctx_t ctx)

    Simultaneous initialize and setting.

.. function:: int gr_lil_mat_update(gr_lil_mat_t res, const gr_lil_mat_t_t src, gr_ctx_t ctx)

    Update *res* with the nonzeros in *src*.  That is, any (row, col) indices in *res* which also appear
    in *src* are overwritten with their values in *src*.  Any indices in *res* which do
    not appear in *src* are left unchanged.

.. function:: void gr_lil_mat_window_init(gr_lil_mat_t window, const gr_lil_mat_t mat, slong r1, slong c1, slong r2, slong c2, gr_ctx_t ctx)
              void gr_lil_mat_window_clear(gr_lil_mat_t window, gr_ctx_t ctx)

    A window is a view on a submatrix with a given interval of rows and columns, and is provided
    by using pointer offsets into the given matrix (so no copying is performed). The window
    produced is read-only. TO BE IMPLEMENTD

.. function:: int gr_scr_mat_permute_cols(gr_scr_mat_t mat, slong * perm, gr_ctx_t ctx)
              int gr_lil_mat_permute_cols(gr_lil_mat_t mat, slong * perm, gr_ctx_t ctx)
              int gr_coo_mat_permute_cols(gr_coo_mat_t mat, slong * perm, gr_ctx_t ctx)

    Permute the columns in *mat* according to the given permutation, i.e., ``mat[r][perm[i]] = mat[i]``.

.. function:: int gr_lil_mat_swap_rows(gr_lil_mat_t mat, slong * perm, slong r, slong s, gr_ctx_t ctx)
              int gr_lil_mat_permute_rows(gr_lil_mat_t mat, const slong * perm, gr_ctx_t ctx)
              int gr_lil_mat_invert_rows(gr_lil_mat_t mat, slong * perm, gr_ctx_t ctx)
              int gr_scr_mat_invert_rows(gr_scr_mat_t mat, slong * perm, gr_ctx_t ctx)

    Swap two rows in the matrix *mat*, permute the rows according to the given permutation,
    or invert all the rows. Note that the permutation *perm* is an input for the rows permutation
    function (required to be non-null), while for the other functions it may be optionally provided
    to keep track of the permutation(s) performed.

    Permuting the rows is TO BE IMPLEMENTED.

    Because of the nature of the sparse compressed row representation, swapping
    and permuting rows is an expensive operation and thus not provided.

Comparison
--------------------------------------------------------------------------------

.. function:: truth_t gr_csr_mat_is_zero(const gr_csr_mat_t mat, gr_ctx_t ctx) 
              truth_t gr_lil_mat_is_zero(const gr_lil_mat_t mat, gr_ctx_t ctx) 
              truth_t gr_coo_mat_is_zero(gr_coo_mat_t mat, gr_ctx_t ctx) 

    Return ``T_TRUE`` if *mat* has no nonzeroes, ``T_FALSE`` if has any element
    known to be nonzero, and ``T_UNKNOWN`` otherwise. For a coordinate-form matrix,
    if it is not in canonical form, it is canonicalized before returning.

.. function:: truth_t gr_csr_mat_is_one(const gr_csr_mat_t mat, gr_ctx_t ctx) 
              truth_t gr_lil_mat_is_one(const gr_lil_mat_t mat, gr_ctx_t ctx) 
              truth_t gr_coo_mat_is_one(const gr_coo_mat_t mat, gr_ctx_t ctx) 
              truth_t gr_csr_mat_is_neg_one(const gr_csr_mat_t mat, gr_ctx_t ctx) 
              truth_t gr_lil_mat_is_neg_one(const gr_lil_mat_t mat, gr_ctx_t ctx) 
              truth_t gr_coo_mat_is_neg_one(const gr_coo_mat_t mat, gr_ctx_t ctx) 

    Return ``T_TRUE`` if *mat* is square, has one nonzero element in each row (at the
    corresponding column), and that element is known to be equal to (negative) one; ``T_FALSE`` if
    the matrix is not square, has any off-diagonal element known to be nonzero, or has
    any diagonal known to not be (negative) one; and ``T_UNKNOWN`` otherwise.

    Functions for coordinate-form matrices are TO BE IMPLEMENTED.

.. function:: truth_t gr_csr_mat_equal(const gr_csr_mat_t mat1, const gr_csr_mat_t mat2, gr_ctx_t ctx)
              truth_t gr_lil_mat_equal(const gr_lil_mat_t mat1, const gr_lil_mat_t mat2, gr_ctx_t ctx)
              truth_t gr_coo_mat_equal(const gr_coo_mat_t mat1, const gr_coo_mat_t mat2, gr_ctx_t ctx)
              truth_t gr_csr_mat_equal_lil_mat(const gr_csr_mat_t mat1, const gr_lil_mat_t mat2, gr_ctx_t ctx)

    Returns ``T_TRUE`` if *mat1* and *mat2* represent the same matrix and ``T_FALSE`` otherwise.


Output
--------------------------------------------------------------------------------

.. function:: int gr_csr_mat_write_nz(gr_stream_t out, const gr_csr_mat_t mat, gr_ctx_t ctx)
              int gr_lil_mat_write_nz(gr_stream_t out, const gr_lil_mat_t mat, gr_ctx_t ctx)
              int gr_coo_mat_write_nz(gr_stream_t out, const gr_coo_mat_t mat, gr_ctx_t ctx)

    Write the nonzeros of *mat* to the stream *out*.  Using ``gr_mat_write`` together with
    ``gr_mat_set_csr_mat``, ``gr_mat_set_lil_mat``, or ``gr_mat_set_lil_mat``,
    if it is desired to write out the entire matrix, zeros and all.

.. function:: int gr_csr_mat_print_nz(const gr_csr_mat_t mat, gr_ctx_t ctx)
              int gr_lil_mat_print_nz(const gr_lil_mat_t mat, gr_ctx_t ctx)
              int gr_coo_mat_print_nz(const gr_coo_mat_t mat, gr_ctx_t ctx)

    Write the nonzeros of *mat* to *stdout*.  Using ``gr_mat_write`` together with
    ``gr_mat_set_csr_mat``, ``gr_mat_set_lil_mat``, or ``gr_mat_set_lil_mat``,
    if it is desired to write out the entire matrix, zeros and all.


Arithmetic
--------------------------------------------------------------------------------

.. function:: int gr_csr_mat_neg(gr_csr_mat_t res, const gr_csr_mat_t src, gr_ctx_t ctx)
              int gr_lil_mat_neg(gr_lil_mat_t res, const gr_lil_mat_t src, gr_ctx_t ctx)
              int gr_coo_mat_neg(gr_coo_mat_t res, const gr_coo_mat_t src, gr_ctx_t ctx)

    Set *res* to -*src*.

.. function:: int gr_lil_mat_add(gr_lil_mat_t res, const gr_lil_mat_t src1, const gr_lil_mat_t src2, slong len, gr_ctx_t ctx)
              int gr_lil_mat_sub(gr_lil_mat_t res, const gr_lil_mat_t src1, const gr_lil_mat_t src2, slong len, gr_ctx_t ctx)
              int gr_lil_mat_mul(gr_lil_mat_t res, const gr_lil_mat_t src1, const gr_lil_mat_t src2, slong len, gr_ctx_t ctx)
    
    Componentwise operations.  (We do not provide analogous division or exponentiation
    routines due since sparse inputs to these operations would be undefined or
    fully dense.)

.. function:: int gr_lil_mat_addmul_scalar(gr_lil_mat_t res, const gr_lil_mat_t src, gr_srcptr c, gr_ctx_t ctx)
              int gr_lil_mat_submul_scalar(gr_lil_mat_t res, const gr_lil_mat_t src, gr_srcptr c, gr_ctx_t ctx)
    
    Componentwise add and sub mul, with different options for the scalar.



Componentwise multiplication and division
--------------------------------------------------------------------------------

.. function:: int gr_csr_mat_mul_scalar(gr_csr_mat_t dst, const gr_csr_mat_t src, gr_srcptr c, gr_ctx_t ctx)
              int gr_csr_mat_mul_scalar_si(gr_csr_mat_t dst, const gr_csr_mat_t src, slong c, gr_ctx_t ctx)
              int gr_csr_mat_mul_scalar_ui(gr_csr_mat_t dst, const gr_csr_mat_t src, ulong c, gr_ctx_t ctx)
              int gr_csr_mat_mul_scalar_fmpz(gr_csr_mat_t dst, const gr_csr_mat_t src, const fmpz_t c, gr_ctx_t ctx)
              int gr_csr_mat_mul_scalar_fmpq(gr_csr_mat_t dst, const gr_csr_mat_t src, const fmpq_t c, gr_ctx_t ctx)
              int gr_csr_mat_mul_scalar_2exp_si(gr_csr_mat_t dst, const gr_csr_mat_t src, slong c, gr_ctx_t ctx)
              int gr_csr_mat_div_scalar(gr_csr_mat_t dst, const gr_csr_mat_t src, gr_srcptr c, gr_ctx_t ctx)
              int gr_csr_mat_div_scalar_si(gr_csr_mat_t dst, const gr_csr_mat_t src, slong c, gr_ctx_t ctx)
              int gr_csr_mat_div_scalar_ui(gr_csr_mat_t dst, const gr_csr_mat_t src, ulong c, gr_ctx_t ctx)
              int gr_csr_mat_div_scalar_fmpz(gr_csr_mat_t dst, const gr_csr_mat_t src, const fmpz_t c, gr_ctx_t ctx)
              int gr_csr_mat_div_scalar_fmpq(gr_csr_mat_t dst, const gr_csr_mat_t src, const fmpq_t c, gr_ctx_t ctx)
              int gr_csr_mat_divexact_scalar(gr_csr_mat_t dst, const gr_csr_mat_t src, gr_srcptr c, gr_ctx_t ctx)
              int gr_csr_mat_divexact_scalar_si(gr_csr_mat_t dst, const gr_csr_mat_t src, slong c, gr_ctx_t ctx)
              int gr_csr_mat_divexact_scalar_ui(gr_csr_mat_t dst, const gr_csr_mat_t src, ulong c, gr_ctx_t ctx)
              int gr_csr_mat_divexact_scalar_fmpz(gr_csr_mat_t dst, const gr_csr_mat_t src, const fmpz_t c, gr_ctx_t ctx)
              int gr_csr_mat_divexact_scalar_fmpq(gr_csr_mat_t dst, const gr_csr_mat_t src, const fmpq_t c, gr_ctx_t ctx)
              
              int gr_lil_mat_mul_scalar(gr_lil_mat_t dst, const gr_lil_mat_t src, gr_srcptr c, gr_ctx_t ctx)
              int gr_lil_mat_mul_scalar_si(gr_lil_mat_t dst, const gr_lil_mat_t src, slong c, gr_ctx_t ctx)
              int gr_lil_mat_mul_scalar_ui(gr_lil_mat_t dst, const gr_lil_mat_t src, ulong c, gr_ctx_t ctx)
              int gr_lil_mat_mul_scalar_fmpz(gr_lil_mat_t dst, const gr_lil_mat_t src, const fmpz_t c, gr_ctx_t ctx)
              int gr_lil_mat_mul_scalar_fmpq(gr_lil_mat_t dst, const gr_lil_mat_t src, const fmpq_t c, gr_ctx_t ctx)
              int gr_lil_mat_mul_scalar_2exp_si(gr_lil_mat_t dst, const gr_lil_mat_t src, slong c, gr_ctx_t ctx)
              int gr_lil_mat_div_scalar(gr_lil_mat_t dst, const gr_lil_mat_t src, gr_srcptr c, gr_ctx_t ctx)
              int gr_lil_mat_div_scalar_si(gr_lil_mat_t dst, const gr_lil_mat_t src, slong c, gr_ctx_t ctx)
              int gr_lil_mat_div_scalar_ui(gr_lil_mat_t dst, const gr_lil_mat_t src, ulong c, gr_ctx_t ctx)
              int gr_lil_mat_div_scalar_fmpz(gr_lil_mat_t dst, const gr_lil_mat_t src, const fmpz_t c, gr_ctx_t ctx)
              int gr_lil_mat_div_scalar_fmpq(gr_lil_mat_t dst, const gr_lil_mat_t src, const fmpq_t c, gr_ctx_t ctx)
              int gr_lil_mat_divexact_scalar(gr_lil_mat_t dst, const gr_lil_mat_t src, gr_srcptr c, gr_ctx_t ctx)
              int gr_lil_mat_divexact_scalar_si(gr_lil_mat_t dst, const gr_lil_mat_t src, slong c, gr_ctx_t ctx)
              int gr_lil_mat_divexact_scalar_ui(gr_lil_mat_t dst, const gr_lil_mat_t src, ulong c, gr_ctx_t ctx)
              int gr_lil_mat_divexact_scalar_fmpz(gr_lil_mat_t dst, const gr_lil_mat_t src, const fmpz_t c, gr_ctx_t ctx)
              int gr_lil_mat_divexact_scalar_fmpq(gr_lil_mat_t dst, const gr_lil_mat_t src, const fmpq_t c, gr_ctx_t ctx)
              
              int gr_coo_mat_mul_scalar(gr_coo_mat_t dst, const gr_coo_mat_t src, gr_srcptr c, gr_ctx_t ctx)
              int gr_coo_mat_mul_scalar_si(gr_coo_mat_t dst, const gr_coo_mat_t src, slong c, gr_ctx_t ctx)
              int gr_coo_mat_mul_scalar_ui(gr_coo_mat_t dst, const gr_coo_mat_t src, ulong c, gr_ctx_t ctx)
              int gr_coo_mat_mul_scalar_fmpz(gr_coo_mat_t dst, const gr_coo_mat_t src, const fmpz_t c, gr_ctx_t ctx)
              int gr_coo_mat_mul_scalar_fmpq(gr_coo_mat_t dst, const gr_coo_mat_t src, const fmpq_t c, gr_ctx_t ctx)
              int gr_coo_mat_mul_scalar_2exp_si(gr_coo_mat_t dst, const gr_coo_mat_t src, slong c, gr_ctx_t ctx)
              int gr_coo_mat_div_scalar(gr_coo_mat_t dst, const gr_coo_mat_t src, gr_srcptr c, gr_ctx_t ctx)
              int gr_coo_mat_div_scalar_si(gr_coo_mat_t dst, const gr_coo_mat_t src, slong c, gr_ctx_t ctx)
              int gr_coo_mat_div_scalar_ui(gr_coo_mat_t dst, const gr_coo_mat_t src, ulong c, gr_ctx_t ctx)
              int gr_coo_mat_div_scalar_fmpz(gr_coo_mat_t dst, const gr_coo_mat_t src, const fmpz_t c, gr_ctx_t ctx)
              int gr_coo_mat_div_scalar_fmpq(gr_coo_mat_t dst, const gr_coo_mat_t src, const fmpq_t c, gr_ctx_t ctx)
              int gr_coo_mat_divexact_scalar(gr_coo_mat_t dst, const gr_coo_mat_t src, gr_srcptr c, gr_ctx_t ctx)
              int gr_coo_mat_divexact_scalar_si(gr_coo_mat_t dst, const gr_coo_mat_t src, slong c, gr_ctx_t ctx)
              int gr_coo_mat_divexact_scalar_ui(gr_coo_mat_t dst, const gr_coo_mat_t src, ulong c, gr_ctx_t ctx)
              int gr_coo_mat_divexact_scalar_fmpz(gr_coo_mat_t dst, const gr_coo_mat_t src, const fmpz_t c, gr_ctx_t ctx)
              int gr_coo_mat_divexact_scalar_fmpq(gr_coo_mat_t dst, const gr_coo_mat_t src, const fmpq_t c, gr_ctx_t ctx)

    Set *dst* to be *src* multiplied or divided by *c*.
    (Addition and subtraction are not provided because they would create
    dense output.)

Arithmetic into dense matrices
--------------------------------------------------------------------------------

.. function:: int gr_mat_update_lil_mat_nz(gr_ptr dres, const gr_lil_mat_t src, gr_ctx_t ctx)
              int gr_mat_add_lil_mat(gr_ptr dres, gr_srcptr dvec1, const gr_lil_mat_t svec2, gr_ctx_t ctx)
              int gr_mat_sub_lil_mat(gr_ptr dres, gr_srcptr dvec1, const gr_lil_mat_t svec2, gr_ctx_t ctx)
              int gr_mat_mul_lil_mat_nz(gr_ptr dres, gr_srcptr dvec1, const gr_lil_mat_t svec2, gr_ctx_t ctx)
              int gr_mat_div_lil_mat_nz(gr_ptr dres, gr_srcptr dvec1, const gr_lil_mat_t svec2, gr_ctx_t ctx)
              int gr_mat_addmul_lil_mat_scalar(gr_ptr dres, const gr_lil_mat_t svec, gr_srcptr c, gr_ctx_t ctx)
              int gr_mat_submul_lil_mat_scalar(gr_ptr dres, const gr_lil_mat_t svec, gr_srcptr c, gr_ctx_t ctx)
              int gr_mat_addmul_lil_mat_scalar_si(gr_ptr dres, const gr_lil_mat_t svec, slong c, gr_ctx_t ctx)
              int gr_mat_submul_lil_mat_scalar_si(gr_ptr dres, const gr_lil_mat_t svec, slong c, gr_ctx_t ctx)
              int gr_mat_addmul_lil_mat_scalar_fmpz(gr_ptr dres, const gr_lil_mat_t svec, const fmpz_t c, gr_ctx_t ctx)
              int gr_mat_submul_lil_mat_scalar_fmpz(gr_ptr dres, const gr_lil_mat_t svec, const fmpz_t c, gr_ctx_t ctx)
    
    These functions facilitate accumulating a sparse matrix into a dense
    target.  They have one dense input, one sparse input, and a dense output
    (where the dense input and output are the same for the fused operations).
    For all functions, it is assumed that *dres* and *dvec1* have the same
    shape as *svec* or *svec2*, as appropriate.  All functions only modify
    the locations in *dres* at which the sparse matrix has a nonzero value:
    in particular, the functions *gr_mat_mul_lil_mat_nz* and
    *gr_mat_div_lil_mat_nz* behave very differently from their dense counterparts.

Sum and product
--------------------------------------------------------------------------------

.. function:: int gr_csr_mat_sum(gr_ptr res, const gr_csr_mat_t mat, gr_ctx_t ctx)
              int gr_lil_mat_sum(gr_ptr res, const gr_lil_mat_t mat, gr_ctx_t ctx)
              int gr_coo_mat_sum(gr_ptr res, const gr_coo_mat_t mat, gr_ctx_t ctx)

    Set *res* to the sum of the entries in *mat*.

.. function:: int gr_csr_mat_nz_product(gr_ptr res, const gr_csr_mat_t mat, gr_ctx_t ctx)
              int gr_lil_mat_nz_product(gr_ptr res, const gr_lil_mat_t mat, gr_ctx_t ctx)
              int gr_coo_mat_nz_product(gr_ptr res, const gr_coo_mat_t mat, gr_ctx_t ctx)

    Set *res* to the product of the nonzero entries in *mat*.


Matrix multiplication
--------------------------------------------------------------------------------

.. function:: int gr_csr_mat_mul_vec(gr_vec_t v, const gr_csr_mat_t A, const gr_vec_t u, gr_ctx_t ctx)
              int gr_lil_mat_mul_vec(gr_vec_t v, const gr_lil_mat_t A, const gr_vec_t u, gr_ctx_t ctx)

    Set *v* equal to `A \cdot u`, i.e., right multiplication by *u*.

.. function:: int gr_csr_mat_mul_mat_transpose(gr_mat_t Ct, const gr_csr_mat_t A, const gr_mat_t Bt, gr_ctx_t ctx)
              int gr_lil_mat_mul_mat_transpose(gr_mat_t Ct, const gr_lil_mat_t A, const gr_mat_t Bt, gr_ctx_t ctx)

    Set *C^T* equal to `A \cdot B^T`, i.e., right multiplication with *B* and *C* considered to
    be column matrices.

.. function:: int gr_csr_mat_mul_mat(gr_mat_t C, const gr_csr_mat_t A, const gr_mat_t B, gr_ctx_t ctx)
              int gr_lil_mat_mul_mat(gr_mat_t C, const gr_lil_mat_t A, const gr_mat_t B, gr_ctx_t ctx)

    Set *C* equal to `A \cdot B`, i.e., perform right multiplication by *B*.


Solving, nullvector, and nullspace computation
--------------------------------------------------------------------------------

.. function:: int gr_lil_mat_solve_lanczos(gr_ptr x, const gr_lil_mat_t M, gr_srcptr b, flint_rand_t state, gr_ctx_t ctx)
              int gr_lil_mat_solve_block_lanczos(gr_ptr x, const gr_lil_mat_t M, gr_srcptr b, slong block_size, flint_rand_t state, gr_ctx_t ctx)
              int gr_lil_mat_solve_wiedemann(gr_ptr x, const gr_lil_mat_t M, gr_srcptr b, gr_ctx_t ctx)
              int gr_lil_mat_solve_block_wiedemann(gr_ptr x, const gr_lil_mat_t M, gr_srcptr b, slong block_size, flint_rand_t state, gr_ctx_t ctx)

    Solve `Mx = b` for a sparse matrix `M` and *M->c* long column vector `b`, using either Lanczos' or Wiedemann's algorithm.
    Both are randomized algorithms which use a given random state machine, and thus may fail without provably no solution
    (returning `GR_UNABLE`). Both have block variants which are better for large matrices, and take an extra parameter of
    `block_size`. The (block) Wiedemann algorithm requires the given matrix to be square, but not symmetric: the Lanczos
    algorithm requires both, but assumes the given matrix may be neither, so instead solves `M^TMx = M^Tb`. Thus, it may
    return a *pseudo-solution*, which solves the latter but not the former.

.. function:: int gr_lil_mat_nullvector_lanczos(gr_ptr x, const gr_lil_mat_t M, flint_rand_t state, gr_ctx_t ctx)
              int gr_lil_mat_nullvector_block_lanczos(gr_ptr x, const gr_lil_mat_t M, slong block_size, flint_rand_t state, gr_ctx_t ctx)
              int gr_lil_mat_nullvector_wiedemann(gr_ptr x, const gr_lil_mat_t M, flint_rand_t state, gr_ctx_t ctx)
              int gr_lil_mat_nullvector_block_wiedemann(gr_ptr x, const gr_lil_mat_t M, slong block_size, flint_rand_t state, gr_ctx_t ctx)

    Find a nullvector for the sparse matrix `M`, using (block) Lanczos or Wiedemann.

.. function:: int gr_lil_mat_nullspace(gr_mat_t X, const gr_lil_mat_t M, flint_rand_t state, slong max_iters, const char *algorithm, slong block_size, gr_ctx_t ctx)
              int gr_lil_mat_nullspace_lanczos(gr_mat_t X, const gr_lil_mat_t M, flint_rand_t state, slong max_iters, gr_ctx_t ctx)
              int gr_lil_mat_nullspace_wiedemann(gr_mat_t X, const gr_lil_mat_t M, flint_rand_t state, slong max_iters, gr_ctx_t ctx)
              gr_lil_mat_nullspace_block_lanczos(gr_mat_t X, const gr_lil_mat_t M, flint_rand_t state, slong max_iters, slong block_size, gr_ctx_t ctx)
              int gr_lil_mat_nullspace_block_wiedemann(gr_mat_t X, const gr_lil_mat_t M, flint_rand_t state, slong max_iters, slong block_size, gr_ctx_t ctx)

    Find the nullspace for the sparse matrix `M`, using (block) Lanczos or Wiedemann.

.. raw:: latexint

    \newpage
