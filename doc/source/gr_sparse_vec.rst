.. _gr-sparse-vec:

**gr_sparse_vec.h** -- sparse vectors over generic rings
===============================================================================

A :type:`gr_sparse_vec_t` is a represention of a vector over a generic
ring which is optimized for the situation where there are many entries which
are zero (the vector is *sparse*, or has *low density*).  Internally, the
vector records only the positions of the nonzeros and their values.
Although technically the type and all its functions will operate correctly
with any density, a regular :type:`gr_vec_t` may be more efficient in the number
of nonzeros is not small.

Types and basic access
--------------------------------------------------------------------------------

.. type:: gr_sparse_vec_struct

.. type:: gr_sparse_vec_t
        
    This struct contains:

    * an `slong` value ``length``, the nominal length of the vector,
    * an `slong` value ``nnz``, the number of nonzero elements in vector,
    * an `slong` value ``alloc``, the maximum number of nonzero element currently storable in the vector,
    * a pointer to an `slong` array of indices (``inds``) of nonzeros.
    * a :type:`gr_ptr` array of values of nonzeros (``nzs``).

    The methods for sparse vectors maintain the following properties:

    * the indices ``inds`` are unique and sorted into strictly increasing order;
    * the ``nzs`` are nonzero (``gr_is_zero(entry, ctx)`` does not return ``T_TRUE``).
    * ``nnz <= alloc <= length``.
    
    A ``gr_sparse_vec_t`` is defined as an array of length one of type
    ``gr_sparse_vec_struct``, permitting a ``gr_sparse_vec_t`` to
    be passed by reference.

.. function:: void gr_sparse_vec_init(gr_sparse_vec_t vec, slong len, gr_ctx_t ctx)

    Initializes *vec* to an empty vector of length *len* in the ring *ctx*. The
    length must be nonnegative. No storage for nonzeroes is allocated.

.. function:: void gr_sparse_vec_clear(gr_sparse_vec_t vec, gr_ctx_t ctx)

    Clears the vector *vec*.

.. function:: slong gr_sparse_vec_length(const gr_sparse_vec_t vec)

    Return the nominal length of the vector (note: not the number of nonzeros).

.. function:: void gr_sparse_vec_set_length(gr_sparse_vec_t vec, slong len, gr_ctx_t ctx)

    Set the nominal length of the vector *vec* to *len*.  If *len* is smaller than
    the current length of *vec*, any entries whose indices are at least *len*
    are truncated.  That is, the number of nonzeros can change.

.. function:: slong gr_sparse_vec_nnz(const gr_sparse_vec_t vec)

    Get the number of nonzeros in *vec*.

.. function:: void gr_sparse_vec_fit_nnz(gr_sparse_vec_t vec, slong nnz, gr_ctx_t ctx)

    Ensure that *vec* has enough storage to hold at least *nnz* nonzeros.  This does
    not change the length of the vector or the number of nonzeros stored.

.. function:: void gr_sparse_vec_shrink_to_nnz(gr_sparse_vec_t vec, gr_ctx_t ctx)

    Reallocate the storage in *vec* down the current number of nonzeros.

.. function:: int gr_sparse_vec_from_entries(gr_sparse_vec_t vec, ulong * inds, gr_srcptr entries, slong nnz, int is_canonical, gr_ctx_t ctx)

    Construct *vec* from the sparse data given by *inds* and *entries* of length *nnz*.
    If ``is_canonical`` is true, the indices are assumed to be sorted and unique, and
    the entries are assumed to be nonzero. Otherwise, the function will sort and compress
    the entries to leave the vector in sorted, unique, nonzero form.

.. function:: int gr_sparse_vec_randtest(gr_sparse_vec_t vec, slong nnz, int replacement, flint_rand_t state, gr_ctx_t ctx)

    Initialize *vec* to a random vector with ``nnz`` nonzero columns, sampled
    either with or without replacement. For small values of ``nnz`` (below 
    ``1/sqrt(length)``), sampling with replacement is faster and likely to
    give the right number of nonzeroes. For larger values of ``nnz``, this
    function uses reservoir sampling to get the specified number of nonzeroes
    with uniform sampling. The entry at each column so sampled
    is generated using ``gr_randtest_nonzero``.

.. function:: int gr_sparse_vec_randtest(gr_sparse_vec_t vec, double prob, int replacement, flint_rand_t state, gr_ctx_t ctx)

    Initialize *vec* to a random vector, in which each column is independently
    sampled with probability ``prob``. The entry at each column so sampled
    is generated using ``gr_randtest_nonzero``.


Getting, setting and conversion
--------------------------------------------------------------------------------

.. macro:: GR_SPARSE_VEC_IND(vec, i)

    Macro to access the index of the *i*-th nonzero.
    The index must be in bounds (between 0 and ``vec->nnz``).

.. macro:: GR_SPARSE_VEC_ENTRY(vec, i, sz)

    Access the value of the *i*-th nonzero.
    The index must be in bounds (between 0 and ``vec->nnz``).
    (See also `gr_sparse_vec_get_entry`)

.. function:: ulong * gr_sparse_vec_ind_ptr(gr_sparse_vec_t vec, slong i, gr_ctx_t ctx)
              const ulong * gr_sparse_vec_ind_srcptr(const gr_sparse_vec_t vec, slong i, gr_ctx_t ctx)
              gr_ptr gr_sparse_vec_entry_ptr(gr_sparse_vec_t vec, slong i, gr_ctx_t ctx)
              gr_srcptr gr_sparse_vec_entry_srcptr(const gr_sparse_vec_t vec, slong i, gr_ctx_t ctx)

    Get pointers to the index or value of the *i*th nonzero. Unlike the above
    macros, these perform bounds checking, and return NULL if out of bounds.

.. function:: gr_ptr gr_sparse_vec_find_entry(gr_sparse_vec_t vec, slong ind, gr_ctx_t ctx)

    Return a pointer to the entry at index *ind*. If *ind* is not a index
    in which *vec* contains a nonzero (including if it is out of bounds),
    the function returns NULL. This is performed via binary search, so
    is worst-case logarithmic time.

.. function:: int gr_sparse_vec_get_entry(gr_ptr res, gr_sparse_vec_t vec, slong ind, gr_ctx_t ctx)

    Set *res* to be a copy of the entry at index *ind*.  If *ind* is not a index
    in which *vec* contains a nonzero, *res* is set to zero. This is performed via binary search, so
    is worst-case logarithmic time.

.. function:: int gr_sparse_vec_set_entry(gr_sparse_vec_t vec, slong ind, gr_srcptr entry, gr_ctx_t ctx)

    Set the the value at index *ind* to be *entry*.  Because of the way sparse
    vectors are represented, it is not efficient to call this function
    repeatedly (it is worst-case linear time in the number of nonzeros in *vec*). 
    In such cases, it is preferrable to construct a new sparse vector from the entries
    (using `gr_sparse_vec_set_from_entries` or `gr_sparse_vec_set_from_entries_sorted_deduped`)
    passed to `gr_sparse_vec_update`.

.. function:: int gr_sparse_vec_set(gr_sparse_vec_t dst, const gr_sparse_vec_t src, gr_ctx_t ctx)

    Set *dst* to a copy of *src*.

.. function:: int gr_sparse_vec_slice(gr_sparse_vec_t dst, const gr_sparse_vec_t src, slong ind_start, slong ind_end, gr_ctx_t ctx)

    Set *dst* to a copy of the slice of *src* given by any entries whose
    indices lie in the half open interval ``[ind_start, ind_end)``.
    Column indices are shifted by *ind_start* (a index of ``ind_start``
    would become ``0``).

.. function:: int gr_sparse_vec_set_vec(gr_sparse_vec_t dst, gr_srcptr src, slong len, gr_ctx_t ctx)

    Set *dst* from the (raw) dense vector *src* of length *len*.
    
.. function:: int gr_vec_set_sparse_vec(gr_ptr dst, gr_sparse_vec_t src, gr_ctx_t ctx)

    Set the raw vector *dst* from sparse vector *src*. *dst* must have
    sufficient space (i.e., ``vec->length``).

.. function:: void gr_sparse_vec_swap(gr_sparse_vec_t vec1, gr_sparse_vec_t vec2, gr_ctx_t ctx)

    Swap the sparse vectors *vec1* and *vec2*. The entries are swapped as pointers, so
    no additional allociation or copying is performed.

.. function:: int gr_sparse_vec_zero(gr_sparse_vec_t vec, gr_ctx_t ctx)

    Set *vec* to the zero vector.

.. function:: int gr_sparse_vec_one(gr_sparse_vec_t vec, slong ind, gr_ctx_t ctx)

    Set *vec* to the elementary vector with a one at index *ind*.

.. function:: int gr_sparse_vec_permute_inds(gr_sparse_vec_t vec, const gr_sparse_vec_t src, slong * p, gr_ctx_t ctx)

    Set *vec* to a copy of *src* with the indices permuted.  The
    indices are shifted as: ``vec[p[i]] = src[i]``.

Comparison
--------------------------------------------------------------------------------

.. function:: truth_t gr_sparse_vec_equal(const gr_sparse_vec_t vec1, const gr_sparse_vec_t vec2, gr_ctx_t ctx)

    Returns ``T_TRUE`` if *vec1* and *vec2* represent the same vector and ``T_FALSE`` otherwise.

.. function:: truth_t gr_sparse_vec_is_zero(const gr_sparse_vec_t vec, gr_ctx_t ctx) 

    Return ``T_TRUE`` if *vec* represents the zero vector and ``T_FALSE`` otherwise.


Output
--------------------------------------------------------------------------------

.. function:: int gr_sparse_vec_write_nz(gr_stream_t out, const gr_sparse_vec_t vec, gr_ctx_t ctx)

    Write the nonzeros of *vec* to the stream *out*.  See ``gr_vec_set_sparse_vec``
    if it is desired to print out the entire vector, zeros and all.

.. function:: int gr_sparse_vec_print_nz(const gr_sparse_vec_t vec, gr_ctx_t ctx)

    Print the nonzeros of *vec* to ``stdout``.  See ``gr_vec_set_sparse_vec``
    if it is desired to print out the entire vector, zeros and all.


Arithmetic
--------------------------------------------------------------------------------

.. function:: int gr_sparse_vec_neg(gr_sparse_vec_t dst, const gr_sparse_vec_t src, gr_ctx_t ctx)

    Set *dst* to -*src*

.. function:: int gr_sparse_vec_update(gr_sparse_vec_t dst, const gr_sparse_vec_t src, gr_ctx_t ctx)

    Update *dst* with the nonzeros in *src*.  That is, any index in *dst* which also appear
    in *src* are overwritten with their values in *src*.  Any indices in *dst* which do
    not appear in *src* are left unchanged.

.. function:: int gr_sparse_vec_add(gr_sparse_vec_t dst, const gr_sparse_vec_t src1, const gr_sparse_vec_t src2, slong len, gr_ctx_t ctx)
              int gr_sparse_vec_sub(gr_sparse_vec_t dst, const gr_sparse_vec_t src1, const gr_sparse_vec_t src2, slong len, gr_ctx_t ctx)
              int gr_sparse_vec_mul(gr_sparse_vec_t dst, const gr_sparse_vec_t src1, const gr_sparse_vec_t src2, slong len, gr_ctx_t ctx)
    
    Componentwise operations.  (We do not provide analogous division or exponentiation
    routines due since sparse inputs to these operations would be undefined or
    fully dense.)

.. function:: int gr_sparse_vec_add_other(gr_sparse_vec_t dst, const gr_sparse_vec_t src1, const gr_sparse_vec_t src2, gr_ctx_t ctx2, gr_ctx_t ctx)
              int gr_sparse_vec_sub_other(gr_sparse_vec_t dst, const gr_sparse_vec_t src1, const gr_sparse_vec_t src2, gr_ctx_t ctx2, gr_ctx_t ctx)
              int gr_sparse_vec_mul_other(gr_sparse_vec_t dst, const gr_sparse_vec_t src1, const gr_sparse_vec_t src2, gr_ctx_t ctx2, gr_ctx_t ctx)
    
    Componentwise operations where the second input is allowed to have a different ring.

.. function:: int gr_other_add_sparse_vec(gr_sparse_vec_t dst, const gr_sparse_vec_t src1, gr_ctx_t ctx1, const gr_sparse_vec_t src2, gr_ctx_t ctx)
              int gr_other_sub_sparse_vec(gr_sparse_vec_t dst, const gr_sparse_vec_t src1, gr_ctx_t ctx1, const gr_sparse_vec_t src2, gr_ctx_t ctx)
              int gr_other_mul_sparse_vec(gr_sparse_vec_t dst, const gr_sparse_vec_t src1, gr_ctx_t ctx1, const gr_sparse_vec_t src2, gr_ctx_t ctx)
    
    Componentwise operations where the first input is allowed to have a different ring.

.. function:: int gr_sparse_vec_addmul_scalar(gr_sparse_vec_t dst, const gr_sparse_vec_t src, gr_srcptr c, gr_ctx_t ctx)
              int gr_sparse_vec_submul_scalar(gr_sparse_vec_t dst, const gr_sparse_vec_t src, gr_srcptr c, gr_ctx_t ctx)
              int gr_sparse_vec_addmul_scalar_si(gr_sparse_vec_t dst, const gr_sparse_vec_t src, slong c, gr_ctx_t ctx)
              int gr_sparse_vec_submul_scalar_si(gr_sparse_vec_t dst, const gr_sparse_vec_t src, slong c, gr_ctx_t ctx)
    
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


Scalar multiplication and division
--------------------------------------------------------------------------------

.. function:: int gr_sparse_vec_mul_scalar(gr_sparse_vec_t dst, const gr_sparse_vec_t src, gr_srcptr c, gr_ctx_t ctx)
              int gr_sparse_vec_mul_scalar_si(gr_sparse_vec_t dst, const gr_sparse_vec_t src, slong c, gr_ctx_t ctx)
              int gr_sparse_vec_mul_scalar_ui(gr_sparse_vec_t dst, const gr_sparse_vec_t src, ulong c, gr_ctx_t ctx)
              int gr_sparse_vec_mul_scalar_fmpz(gr_sparse_vec_t dst, const gr_sparse_vec_t src, const fmpz_t c, gr_ctx_t ctx)
              int gr_sparse_vec_mul_scalar_fmpq(gr_sparse_vec_t dst, const gr_sparse_vec_t src, const fmpq_t c, gr_ctx_t ctx)
              int gr_sparse_vec_mul_scalar_2exp_si(gr_sparse_vec_t dst, const gr_sparse_vec_t src, slong c, gr_ctx_t ctx)
              int gr_sparse_vec_div_scalar(gr_sparse_vec_t dst, const gr_sparse_vec_t src, gr_srcptr c, gr_ctx_t ctx)
              int gr_sparse_vec_div_scalar_si(gr_sparse_vec_t dst, const gr_sparse_vec_t src, slong c, gr_ctx_t ctx)
              int gr_sparse_vec_div_scalar_ui(gr_sparse_vec_t dst, const gr_sparse_vec_t src, ulong c, gr_ctx_t ctx)
              int gr_sparse_vec_div_scalar_fmpz(gr_sparse_vec_t dst, const gr_sparse_vec_t src, const fmpz_t c, gr_ctx_t ctx)
              int gr_sparse_vec_div_scalar_fmpq(gr_sparse_vec_t dst, const gr_sparse_vec_t src, const fmpq_t c, gr_ctx_t ctx)
              int gr_sparse_vec_divexact_scalar(gr_sparse_vec_t dst, const gr_sparse_vec_t src, gr_srcptr c, gr_ctx_t ctx)
              int gr_sparse_vec_divexact_scalar_si(gr_sparse_vec_t dst, const gr_sparse_vec_t src, slong c, gr_ctx_t ctx)
              int gr_sparse_vec_divexact_scalar_ui(gr_sparse_vec_t dst, const gr_sparse_vec_t src, ulong c, gr_ctx_t ctx)
              int gr_sparse_vec_divexact_scalar_fmpz(gr_sparse_vec_t dst, const gr_sparse_vec_t src, const fmpz_t c, gr_ctx_t ctx)
              int gr_sparse_vec_divexact_scalar_fmpq(gr_sparse_vec_t dst, const gr_sparse_vec_t src, const fmpq_t c, gr_ctx_t ctx)

    Set *dst* to be *src* multiplied or divided by *c*.
    (Addition and subtraction are not provided because they would create
    dense output.)

Sum and product
--------------------------------------------------------------------------------

.. function:: int gr_sparse_vec_sum(gr_ptr res, const gr_sparse_vec_t vec, gr_ctx_t ctx)

    Set *res* to the sum of the entries in *vec*.

.. function:: int gr_sparse_vec_nz_product(gr_ptr res, const gr_sparse_vec_t vec, gr_ctx_t ctx)

    Set *res* to the product of the nonzero entries in *vec*.


Dot product
--------------------------------------------------------------------------------

.. function:: int gr_sparse_vec_dot(gr_ptr res, gr_srcptr c, int subtract, const gr_sparse_vec_t x, const gr_sparse_vec_t y, gr_ctx_t ctx)

    Set *res* equal to `c \pm x \cdot y`.

.. function:: int gr_sparse_vec_dot_vec(gr_ptr res, gr_srcptr c, int subtract, const gr_sparse_vec_t x, const gr_vec_t y, gr_ctx_t ctx)

    Set *res* equal to `c \pm x \cdot y`.




.. raw:: latex

    \newpage
