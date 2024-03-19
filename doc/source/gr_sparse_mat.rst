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
    * A :type:`gr_ptr` array of of length ``nnz`` providing nonzero values.

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
    be passed by reference. Note that ``gr_less_mat_t`` is the less efficient
    way to represent a sparse matrix, but is more convenient for modification.

.. function:: void gr_sparse_mat_nrows(gr_csr_mat_t mat, gr_ctx_t ctx)
              void gr_sparse_mat_ncols(gr_csr_mat_t mat, gr_ctx_t ctx)
              void gr_sparse_mat_nnz(gr_csr_mat_t mat, gr_ctx_t ctx)
              void gr_sparse_mat_nrows(gr_csr_mat_t mat, gr_ctx_t ctx)
              void gr_sparse_mat_ncols(gr_csr_mat_t mat, gr_ctx_t ctx)
              void gr_sparse_mat_nnz(gr_csr_mat_t mat, gr_ctx_t ctx)

    Get the number of rows/columns/nonzeroes of the given sparse matrix,
    in either representation.

.. macro:: GR_CSR_MAT_COL(mat, i, j)
           GR_LIL_MAT_COL(mat, i, j)

    Get a pointer to the column of the *j*-th nonzero in the *i*-th row of *mat*.
    There is no bounds checking.

.. macro:: GR_CSR_MAT_ENTRY(mat, i, j, sz)
           GR_LIL_MAT_ENTRY(mat, i, j, sz)

    Get a pointer to the value of the *j*-th nonzero entry in the *i*-th row of *mat*,
    given the size of each entry (typically obtained as ``ctx->sizeof_elem``).
    There is no bounds checking.

.. function:: ulong * gr_csr_mat_col_ptr(gr_csr_mat_t mat, slong i, slong j)
              const ulong * gr_csr_mat_col_srcptr(const gr_csr_mat_t mat, slong i, slong j)
              ulong * gr_lil_mat_col_ptr(gr_lil_mat_t mat, slong i, slong j)
              const ulong * gr_lil_mat_col_srcptr(const gr_lil_mat_t mat, slong i, slong j)

    Get a (const) pointer to the column of the *j*-th nonzero in the *i*-th row of *mat*.
    There is no bounds checking.

.. function:: gr_ptr gr_csr_mat_entry_ptr(gr_csr_mat_t mat, slong i, slong j, gr_ctx_t ctx)
              gr_srcptr gr_csr_mat_entry_srcptr(gr_csr_mat_t mat, slong i, slong j, gr_ctx_t ctx)
              gr_ptr gr_lil_mat_entry_ptr(gr_lil_mat_t mat, slong i, slong j, gr_ctx_t ctx)
              gr_srcptr gr_lil_mat_entry_srcptr(gr_lil_mat_t mat, slong i, slong j, gr_ctx_t ctx)

    Get a (const) pointer to the *j*-th nonzero entry in the *i*-th row of *mat*.
    There is no bounds checking.

.. function:: void gr_csr_mat_init(gr_csr_mat_t mat, slong rows, slong cols, gr_ctx_t ctx)
              void gr_lil_mat_init(gr_lil_mat_t mat, slong rows, slong cols, gr_ctx_t ctx)

    Initializes *mat* to a *rows* x *cols* matrix with no nonzeros.

.. function:: void gr_csr_mat_clear(gr_csr_mat_t vec, gr_ctx_t ctx)
              void gr_lil_mat_clear(gr_lil_mat_t vec, gr_ctx_t ctx)

    Clears the matrix *mat*.

.. function:: gr_csr_mat_swap(gr_csr_mat_t mat1, gr_csr_mat_t mat2, gr_ctx_t ctx)
              gr_lil_mat_swap(gr_lil_mat_t mat1, gr_lil_mat_t mat2, gr_ctx_t ctx)

    Swap data underlying two sparse matrices (no allocation or copying).

.. function:: void gr_csr_mat_fit_nnz(gr_csr_mat_t mat, slong nnz, gr_ctx_t ctx)

    Ensure that *mat* has enough storage to hold at least *nnz* nonzeros.  This does
    not change the dimensions of the matrix or the number of nonzeros stored.

    Note that, for matrices in *lil* form, one must perform the analogous operation
    on each row.

.. function:: void gr_csr_mat_shrink_to_nnz(gr_csr_mat_t mat, gr_ctx_t ctx)

    Reallocate the storage in *mat* down the current number of nonzeros.

    Note that, for matrices in *lil* form, the analogous operation is performed
    on each row.

.. function:: void gr_csr_mat_set_cols(gr_csr_mat_t mat, slong cols, gr_ctx_t ctx)
              void gr_lil_mat_set_cols(gr_llil_mat_t mat, slong cols, gr_ctx_t ctx)

    Set the nominal number of columns of the matrix *mat* to *cols*.  If *cols* is smaller than
    the current number of columns of *vec*, any entries whose columns are at least *cols*
    are truncated. That is, the number of nonzeros can change.



Getting, setting and conversion
--------------------------------------------------------------------------------

.. function:: int gr_csr_mat_get_entry(gr_ptr res, gr_csr_mat_t vec, slong row, slong col, gr_ctx_t ctx)
              int gr_lil_mat_get_entry(gr_ptr res, gr_lil_mat_t vec, slong row, slong col, gr_ctx_t ctx)

    Set *res* to be the entry at position (*row*, *col*).
    
    Because of the way sparse matrices are represented, this is logarithmic time in the number of nonzeros
    in the specified row.

.. function:: int gr_csr_mat_set_entry(gr_csr_mat_t vec, slong row, slong col, gr_srcptr entry, gr_ctx_t ctx)
              int gr_lil_mat_set_entry(gr_lil_mat_t vec, slong row, slong col, gr_srcptr entry, gr_ctx_t ctx)

    Set the the entry at location (*row*, *col*) to be *entry*.
    
    Because of the way sparse vectors are represented, it is not efficient to call this function
    repeatedly: for the list of lists representation, it is linear time in the number of nonzeros
    in the associated row; for the sparse compressed row representation, it is linear time in the
    overall number of nonzeroes. If possible, the entries to update should be batched up and
    given using `gr_lil_mat_update`, `gr_csr_mat_set_from_coo`,  `gr_lil_mat_set_from_coo`,
    `gr_csr_mat_set_from_coo_sorted_deduped`, or `gr_lil_mat_set_from_coo_sorted_deduped`.

.. function:: int gr_csr_mat_zero(gr_csr_mat_t mat, gr_ctx_t ctx)
              int gr_lil_mat_zero(gr_lil_mat_t mat, gr_ctx_t ctx)

    Set *mat* to the zero matrix.

.. function:: int gr_csr_mat_set(gr_csr_mat_t res, const gr_csr_mat_t mat, gr_ctx_t ctx)
              int gr_lil_mat_set(gr_lil_mat_t res, const gr_lil_mat_t mat, gr_ctx_t ctx)
              int gr_csr_mat_set_lil_mat(gr_csr_mat_t res, const gr_lil_mat_t mat, gr_ctx_t ctx)
              int gr_lil_mat_set_csr_mat(gr_lil_mat_t res, const gr_csr_mat_t mat, gr_ctx_t ctx)

    Set *res* to a copy of *mat*, possibly changing the type of representation.

.. function:: int gr_csr_mat_set_mat(gr_csr_mat_t res, gr_mat_t mat, gr_ctx_t ctx)
              int gr_lil_mat_set_mat(gr_lil_mat_t res, gr_mat_t mat, gr_ctx_t ctx)

    Set *res* from the (nominally) dense matrix *mat*.
    
.. function:: int gr_mat_set_csr_mat(gr_mat_t res, const gr_csr_mat_t mat, gr_ctx_t ctx)
              int gr_mat_set_lil_mat(gr_mat_t res, const gr_lil_mat_t mat, gr_ctx_t ctx)

    Set a dense matrix *res* from the sparse matrix *mat*, for either representation.

.. function:: int gr_csr_mat_init_set(gr_csr_mat_t res, const gr_csr_mat_t mat, gr_ctx_t ctx)
              int gr_lil_mat_init_set(gr_lil_mat_t res, const gr_lil_mat_t mat, gr_ctx_t ctx)
              int gr_csr_mat_init_set_lil_mat(gr_csr_mat_t res, const gr_lil_mat_t mat, gr_ctx_t ctx)
              int gr_lil_mat_init_set_csr_mat(gr_lil_mat_t res, const gr_csr_mat_t mat, gr_ctx_t ctx)
              int gr_csr_mat_init_set_mat(gr_csr_mat_t res, const gr_mat_t mat, gr_ctx_t ctx)
              int gr_lil_mat_init_set_mat(gr_lil_mat_t res, const gr_mat_t mat, gr_ctx_t ctx)
              int gr_mat_init_set_csr_mat(gr_mat_t res, const gr_csr_mat_t mat, gr_ctx_t ctx)
              int gr_mat_init_set_lil_mat(gr_mat_t res, const gr_lil_mat_t mat, gr_ctx_t ctx)

    Simultaneous initialize and setting.

.. function:: int gr_lil_mat_update(gr_lil_mat_t res, const gr_lil_mat_t_t src, gr_ctx_t ctx)

    Update *res* with the nonzeros in *src*.  That is, any (row, col) indices in *res* which also appear
    in *src* are overwritten with their values in *src*.  Any indices in *res* which do
    not appear in *src* are left unchanged.

.. function:: int gr_csr_mat_set_from_coo(gr_csr_mat_t mat, ulong * rows, ulong * cols, gr_srcptr entries, slong nnz, gr_ctx_t ctx)
              int gr_lil_mat_set_from_coo(gr_lil_mat_t vec, ulong * rows, ulong * cols, gr_srcptr entries, slong nnz, gr_ctx_t ctx)

    Set *mat* from the sparse data given by arrays *rows*, *cols*, and *entries*.
    
    These functions allow one to construct a matrix from nonzero data in coo(rdinate) form, i.e.,
    a list of triples of (row, col, entry) (provided as parallel arrays *rows*, *cols*, *entries*).
    No assumption is made that the indices are sorted nor that the entries are nonzero.  The values associated
    with duplicate indices are added together.

.. function:: int gr_csr_mat_set_from_coo_sorted_deduped(gr_csr_mat_t mat, ulong * rows, ulong * cols, gr_srcptr entries, slong nnz, gr_ctx_t ctx)
              int gr_csr_mat_set_from_coo_sorted_deduped(gr_csr_mat_t mat, ulong * rows, ulong * cols, gr_srcptr entries, slong nnz, gr_ctx_t ctx)

    Set *mat* from the sparse data given by arrays *rows*, *cols*, and *entries*.
    
    Similar to the set from coo functions, save that the entries are required to be nonzero and the
    (row, col) pairs distinct and sorted (by row then col).

.. function:: int gr_csr_mat_randtest(gr_csr_mat_t mat, double row_density, slong nrows, slong ncols, flint_rand_t state, gr_ctx_t ctx)
              int gr_lil_mat_randtest(gr_lil_mat_t mat, double row_density, slong nrows, slong ncols, flint_rand_t state, gr_ctx_t ctx)

    Initialize *mat* to a random *nrows* x *ncols* matrix with row density (fraction of nonzeros)
    *density*. The algorithm is suitable when *density* is small.
    Specifically, indices are generated randomly and deduped.  So if the
    density is larger than ``1/sqrt(ncols)``, the true density of the returned vector
    is likely to be lower than *density*.

.. function:: void gr_lil_mat_window_init(gr_lil_mat_t window, const gr_lil_mat_t mat, slong r1, slong c1, slong r2, slong c2, gr_ctx_t ctx)
              void gr_lil_mat_window_clear(gr_lil_mat_t window, gr_ctx_t ctx)

    A window is a view on a submatrix with a given interval of rows and columns, and is provided
    by using pointer offsets into the given matrix (so no copying is performed). The window
    produced is read-only.

.. function:: int gr_scr_mat_permute_cols(gr_scr_mat_t mat, slong * perm, gr_ctx_t ctx)
              int gr_lil_mat_permute_cols(gr_lil_mat_t mat, slong * perm, gr_ctx_t ctx)

    Permute the columns in *mat* according to the given permutation, i.e., ``mat[r][perm[i]] = mat[i]``.

.. function:: int gr_lil_mat_swap_rows(gr_lil_mat_t mat, slong * perm, slong r, slong s, gr_ctx_t ctx)
              int gr_lil_mat_permute_rows(gr_lil_mat_t mat, const slong * perm, gr_ctx_t ctx)
              int gr_lil_mat_invert_rows(gr_lil_mat_t mat, slong * perm, gr_ctx_t ctx)
              int gr_scr_mat_invert_rows(gr_scr_mat_t mat, slong * perm, gr_ctx_t ctx)

    Swap two rows in the matrix *mat*, permute the rows according to the given permutation,
    or invert all the rows. Note that the permutation *perm* is an input for the rows permutation
    function (required to be non-null), while for the other functions it may be optionally provided
    to keep track of the permutation(s) performed.

    Because of the nature of the sparse compressed row representation, swapping
    and permuting rows is an expensive operation and thus not provided.

Comparison
--------------------------------------------------------------------------------

.. function:: truth_t gr_csr_mat_is_zero(const gr_csr_mat_t mat, gr_ctx_t ctx) 
              truth_t gr_lil_mat_is_zero(const gr_lil_mat_t mat, gr_ctx_t ctx) 

    Return ``T_TRUE`` if *mat* has no nonzeroes, ``T_FALSE`` if has any element
    known to be nonzero, and ``T_UNKNOWN`` otherwise.

.. function:: truth_t gr_csr_mat_is_one(const gr_csr_mat_t mat, gr_ctx_t ctx) 
              truth_t gr_lil_mat_is_one(const gr_lil_mat_t mat, gr_ctx_t ctx) 
              truth_t gr_csr_mat_is_neg_one(const gr_csr_mat_t mat, gr_ctx_t ctx) 
              truth_t gr_lil_mat_is_neg_one(const gr_lil_mat_t mat, gr_ctx_t ctx) 

    Return ``T_TRUE`` if *mat* is square, has one nonzero element in each row (at the
    corresponding column), and that element is known to be equal to (negative) one; ``T_FALSE`` if
    the matrix is not square, has any off-diagonal element known to be nonzero, or has
    any diagonal known to not be (negative) one; and ``T_UNKNOWN`` otherwise.

.. function:: truth_t gr_csr_mat_equal(const gr_csr_mat_t mat1, const gr_csr_mat_t mat2, gr_ctx_t ctx)
              truth_t gr_lil_mat_equal(const gr_lil_mat_t mat1, const gr_lil_mat_t mat2, gr_ctx_t ctx)
              truth_t gr_csr_mat_equal_lil_mat(const gr_csr_mat_t mat1, const gr_lil_mat_t mat2, gr_ctx_t ctx)

    Returns ``T_TRUE`` if *mat1* and *mat2* represent the same matrix and ``T_FALSE`` otherwise.


Output
--------------------------------------------------------------------------------

.. function:: int gr_csr_mat_write_nz(gr_stream_t out, const gr_csr_mat_t vec, gr_ctx_t ctx)
              int gr_lil_mat_write_nz(gr_stream_t out, const gr_lil_mat_t vec, gr_ctx_t ctx)

    Write the nonzeros of *mat* to the stream *out*.  See ``gr_mat_set_csr_mat`` and
    ``gr_mat_set_lil_mat`` if it is desired to write out the entire matrix, zeros and all.

.. function:: int gr_csr_mat_print_nz(const gr_csr_mat_t mat, gr_ctx_t ctx)
              int gr_lil_mat_print_nz(const gr_lil_mat_t mat, gr_ctx_t ctx)

    Print the nonzeros of *mat* to ``stdout``.   See ``gr_mat_set_csr_mat`` and
    ``gr_mat_set_lil_mat`` if it is desired to print out the entire matrix, zeros and all.


Arithmetic
--------------------------------------------------------------------------------

.. function:: int gr_csr_mat_neg(gr_csr_mat_t res, const gr_csr_mat_t src, gr_ctx_t ctx)
              int gr_lil_mat_neg(gr_lil_mat_t res, const gr_csr_mat_t src, gr_ctx_t ctx)

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


Arithmetic into dense vectors
--------------------------------------------------------------------------------

.. function:: int gr_vec_update_sparse_vec_nz(gr_ptr dres, const gr_sparse_vec_t src, gr_ctx_t ctx)
              int gr_vec_add_sparse_vec(gr_ptr dres, gr_srcptr dvec1, const gr_sparse_vec_t svec2, gr_ctx_t ctx)
              int gr_vec_sub_sparse_vec(gr_ptr dres, gr_srcptr dvec1, const gr_sparse_vec_t svec2, gr_ctx_t ctx)
              int gr_vec_mul_sparse_vec_nz(gr_ptr dres, gr_srcptr dvec1, const gr_sparse_vec_t svec2, gr_ctx_t ctx)
              int gr_vec_div_sparse_vec_nz(gr_ptr dres, gr_srcptr dvec1, const gr_sparse_vec_t svec2, gr_ctx_t ctx)
              int gr_vec_addmul_sparse_vec_scalar(gr_ptr dres, const gr_sparse_vec_t svec, gr_srcptr c, gr_ctx_t ctx)
              int gr_vec_submul_sparse_vec_scalar(gr_ptr dres, const gr_sparse_vec_t svec, gr_srcptr c, gr_ctx_t ctx)
              int gr_vec_addmul_sparse_vec_scalar_si(gr_ptr dres, const gr_sparse_vec_t svec, slong c, gr_ctx_t ctx)
              int gr_vec_submul_sparse_vec_scalar_si(gr_ptr dres, const gr_sparse_vec_t svec, slong c, gr_ctx_t ctx)
              int gr_vec_addmul_sparse_vec_scalar_fmpz(gr_ptr dres, const gr_sparse_vec_t svec, const fmpz_t c, gr_ctx_t ctx)
              int gr_vec_submul_sparse_vec_scalar_fmpz(gr_ptr dres, const gr_sparse_vec_t svec, const fmpz_t c, gr_ctx_t ctx)
    
    These functions facilitate accumulating a sparse vector into a dense
    target.  They have one dense input, one sparse input, and a dense output
    (where the dense input and output are the same for the fused operations).
    For all functions, it is assumed that *dres* and *dvec1* have the same
    length as *svec* or *svec2*, as appropriate.  All functions only modify
    the locations in *dres* at which the sparse vector has a nonzero value:
    in particular, the functions *gr_vec_mul_sparse_vec_nz* and
    *gr_vec_div_sparse_vec_nz* behave very differently from their dense counterparts.


Componentwise multiplication and division
--------------------------------------------------------------------------------

.. function:: int gr_csr_mat_mul_scalar(gr_csr_mat_t res, const gr_csr_mat_t src, gr_srcptr c, gr_ctx_t ctx)
              int gr_csr_mat_mul_scalar_si(gr_csr_mat_t res, const gr_csr_mat_t src, slong c, gr_ctx_t ctx)
              int gr_csr_mat_mul_scalar_ui(gr_csr_mat_t res, const gr_csr_mat_t src, ulong c, gr_ctx_t ctx)
              int gr_csr_mat_mul_scalar_fmpz(gr_csr_mat_t res, const gr_csr_mat_t src, const fmpz_t c, gr_ctx_t ctx)
              int gr_csr_mat_mul_scalar_fmpq(gr_csr_mat_t res, const gr_csr_mat_t src, const fmpq_t c, gr_ctx_t ctx)
              int gr_csr_mat_mul_scalar_2exp_si(gr_csr_mat_t res, const gr_csr_mat_t src, slong c, gr_ctx_t ctx)
              int gr_csr_mat_div_scalar(gr_csr_mat_t res, const gr_csr_mat_t src, gr_srcptr c, gr_ctx_t ctx)
              int gr_csr_mat_div_scalar_si(gr_csr_mat_t res, const gr_csr_mat_t src, slong c, gr_ctx_t ctx)
              int gr_csr_mat_div_scalar_ui(gr_csr_mat_t res, const gr_csr_mat_t src, ulong c, gr_ctx_t ctx)
              int gr_csr_mat_div_scalar_fmpz(gr_csr_mat_t res, const gr_csr_mat_t src, const fmpz_t c, gr_ctx_t ctx)
              int gr_csr_mat_div_scalar_fmpq(gr_csr_mat_t res, const gr_csr_mat_t src, const fmpq_t c, gr_ctx_t ctx)
              int gr_csr_mat_divexact_scalar(gr_csr_mat_t res, const gr_csr_mat_t src, gr_srcptr c, gr_ctx_t ctx)
              int gr_csr_mat_divexact_scalar_si(gr_csr_mat_t res, const gr_csr_mat_t src, slong c, gr_ctx_t ctx)
              int gr_csr_mat_divexact_scalar_ui(gr_csr_mat_t res, const gr_csr_mat_t src, ulong c, gr_ctx_t ctx)
              int gr_csr_mat_divexact_scalar_fmpz(gr_csr_mat_t res, const gr_csr_mat_t src, const fmpz_t c, gr_ctx_t ctx)
              int gr_csr_mat_divexact_scalar_fmpq(gr_csr_mat_t res, const gr_csr_mat_t src, const fmpq_t c, gr_ctx_t ctx)
              int gr_lil_mat_mul_scalar(gr_lil_mat_t res, const gr_lil_mat_t src, gr_srcptr c, gr_ctx_t ctx)
              int gr_lil_mat_mul_scalar_si(gr_lil_mat_t res, const gr_lil_mat_t src, slong c, gr_ctx_t ctx)
              int gr_lil_mat_mul_scalar_ui(gr_lil_mat_t res, const gr_lil_mat_t src, ulong c, gr_ctx_t ctx)
              int gr_lil_mat_mul_scalar_fmpz(gr_lil_mat_t res, const gr_lil_mat_t src, const fmpz_t c, gr_ctx_t ctx)
              int gr_lil_mat_mul_scalar_fmpq(gr_lil_mat_t res, const gr_lil_mat_t src, const fmpq_t c, gr_ctx_t ctx)
              int gr_lil_mat_mul_scalar_2exp_si(gr_lil_mat_t res, const gr_lil_mat_t src, slong c, gr_ctx_t ctx)
              int gr_lil_mat_div_scalar(gr_lil_mat_t res, const gr_lil_mat_t src, gr_srcptr c, gr_ctx_t ctx)
              int gr_lil_mat_div_scalar_si(gr_lil_mat_t res, const gr_lil_mat_t src, slong c, gr_ctx_t ctx)
              int gr_lil_mat_div_scalar_ui(gr_lil_mat_t res, const gr_lil_mat_t src, ulong c, gr_ctx_t ctx)
              int gr_lil_mat_div_scalar_fmpz(gr_lil_mat_t res, const gr_lil_mat_t src, const fmpz_t c, gr_ctx_t ctx)
              int gr_lil_mat_div_scalar_fmpq(gr_lil_mat_t res, const gr_lil_mat_t src, const fmpq_t c, gr_ctx_t ctx)
              int gr_lil_mat_divexact_scalar(gr_lil_mat_t res, const gr_lil_mat_t src, gr_srcptr c, gr_ctx_t ctx)
              int gr_lil_mat_divexact_scalar_si(gr_lil_mat_t res, const gr_lil_mat_t src, slong c, gr_ctx_t ctx)
              int gr_lil_mat_divexact_scalar_ui(gr_lil_mat_t res, const gr_lil_mat_t src, ulong c, gr_ctx_t ctx)
              int gr_lil_mat_divexact_scalar_fmpz(gr_lil_mat_t res, const gr_lil_mat_t src, const fmpz_t c, gr_ctx_t ctx)
              int gr_lil_mat_divexact_scalar_fmpq(gr_lil_mat_t res, const gr_lil_mat_t src, const fmpq_t c, gr_ctx_t ctx)

    Set *res* to be *src* multiplied or divided by *c*.
    (Addition and subtraction are not provided because they would create
    dense output.)

Sum and product
--------------------------------------------------------------------------------

.. function:: int gr_csr_mat_sum(gr_ptr res, const gr_csr_mat_t mat, gr_ctx_t ctx)
              int gr_lil_mat_sum(gr_ptr res, const gr_lil_mat_t mat, gr_ctx_t ctx)

    Set *res* to the sum of the entries in *mat*.

.. function:: int gr_csr_mat_product(gr_ptr res, const gr_csr_mat_t mat, gr_ctx_t ctx)
              int gr_lil_mat_product(gr_ptr res, const gr_lil_mat_t mat, gr_ctx_t ctx)

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

    This is the standard procedure use by solving, and is maximally cache friendly.

.. function:: int gr_csr_mat_mul_mat(gr_mat_t C, const gr_csr_mat_t A, const gr_mat_t B, gr_ctx_t ctx)
              int gr_lil_mat_mul_mat(gr_mat_t C, const gr_lil_mat_t A, const gr_mat_t B, gr_ctx_t ctx)

    Set *C* equal to `A \cdot B`, i.e., perform right multiplication by *B*.


.. raw:: latex

    \newpage
