.. _fq-zech-sparse-vec:

**fq_zech_sparse_vec.h** -- sparse vectors over finite fields (Zech logarithm representation)
===============================================================================

Description.

Types, macros and constants
-------------------------------------------------------------------------------

.. type:: fq_zech_sparse_entry_t

    Holds a pair (ind, val), where ind (of type slong) is an index into the
    vector and val (of type fq_zech_t) is the value at that index

.. type:: fq_zech_sparse_vec_t

    Holds an array of nonzero entries (of type sparse_entry_struct *) of 
    specified length nnz, sorted by ind

Memory management
--------------------------------------------------------------------------------


.. function:: void fq_zech_sparse_vec_init(fq_zech_sparse_vec_t vec, const fq_zech_ctx_t ctx)

    Initializes an empty sparse vector (no allocation)

.. function:: void fq_zech_sparse_vec_clear(fq_zech_sparse_vec_t vec, slong len, const fq_zech_ctx_t ctx)

    Clears the entries of vec and frees the space allocated for it

.. function:: void fq_zech_sparse_vec_swap(fq_zech_sparse_vec_t vec1, fq_zech_sparse_vec_t vec2, const fq_zech_ctx_t ctx)

    Swaps two vectors (no reallocaton)

Instantiation
--------------------------------------------------------------------------------

.. function:: void fq_zech_sparse_vec_zero(fq_zech_sparse_vec_t vec, const fq_zech_ctx_t ctx)

    Sets vec to zero (the empty sparse vector)

.. function:: void fq_zech_sparse_vec_one(fq_zech_sparse_vec_t vec, slong ind, const fq_zech_ctx_t ctx)

    Sets vec to the ind-th basis vector (single one in position ind)

.. function:: void fq_zech_sparse_vec_set(fq_zech_sparse_vec_t dst, fq_zech_sparse_vec_t src, const fq_zech_ctx_t ctx)

    Makes dst a (deep) copy of src

.. function:: void fq_zech_sparse_vec_set_entry(fq_zech_sparse_vec_t vec, slong ind, const fq_zech_t val, const fq_zech_ctx_t ctx)

    Sets the value at index ind to val, either replacing an existing value or extending
    the array of entries

.. function:: void _fq_zech_sparse_vec_from_entries(fq_zech_sparse_vec_t vec, slong *inds, fq_zech_struct *vals, slong nnz, const fq_zech_ctx_t ctx)

    Constructs vec from a given sequence of indices and associated values, both of length nnz.
    Assumes no duplicate indices

Comparison
--------------------------------------------------------------------------------

.. function:: void fq_zech_sparse_is_zero(fq_zech_sparse_vec_t vec, const fq_zech_ctx_t ctx)

    Checks if the given vector is trivial (empty), returning `1` if so and `0` 
    otherwise

.. function:: void fq_zech_sparse_vec_equal(const fq_zech_sparse_vec_t vec1, const fq_zech_sparse_vec_t vec2, slong ioff, const fq_zech_ctx_t ctx)

    Checks if vec1 equals vec2 (with s specified column offset ioff), returning
    `1` if so and `0` otherwise

Indexing
--------------------------------------------------------------------------------

.. function:: fq_zech_t * fq_zech_sparse_vec_at(fq_zech_sparse_vec_t vec, slong ind, const fq_zech_ctx_t ctx)

    Returns a pointer to the value at the index ind (or NULL if index not found)


Conversion to/from dense vector
--------------------------------------------------------------------------------

.. function:: void fq_zech_sparse_vec_from_dense(fq_zech_sparse_vec_t dst, const fq_zech_struct *src, slong len, const fq_zech_ctx_t ctx)

    Converts the dense vector src of length len to a sparse vector

.. function:: void fq_zech_sparse_vec_to_dense(fq_zech_struct *dst, const fq_zech_sparse_vec_t src, slong len, const fq_zech_ctx_t ctx)

    Converts the sparse vector src to a dense vector of length len


Windows, concatenation, and splitting
--------------------------------------------------------------------------------

.. function:: void fq_zech_sparse_vec_window_init(fq_zech_sparse_vec_t window, const fq_zech_sparse_vec_t vec, slong i1, slong i2, const fq_zech_ctx_t ctx)

    Constructs a window on a the sparse vector vec between indices i1 and i2
    Note that window is only valid as long as original vector remains uzechified

.. function:: void fq_zech_sparse_vec_window_clear(fq_zech_sparse_vec_t window, const fq_zech_ctx_t ctx)

    Clears a window (for safety only)

.. function:: void fq_zech_sparse_vec_concat(fq_zech_sparse_vec_t res, const fq_zech_sparse_vec_t vec1, const fq_zech_sparse_vec_t vec2, slong len1, const fq_zech_ctx_t ctx)

    Concatenates two vectors vec1 and vec2 into res, with indices of vec2 
    offset by len1

.. function:: void fq_zech_sparse_vec_split(fq_zech_sparse_vec_t res1, fq_zech_sparse_vec_t res1, const fq_zech_sparse_vec_t vec, slong ind, const fq_zech_ctx_t ctx)

    Splits vec into two vectors res1 and res2, with res1 containing all entries 
    below index ind and res2 containing the rest

Permutation
--------------------------------------------------------------------------------

.. function:: void fq_zech_sparse_vec_permute_inds(fq_zech_sparse_vec_t vec, slong *P, const fq_zech_ctx_t ctx)

    Permutes the indices of vec according to P, and resorts


Randomization
--------------------------------------------------------------------------------


.. function:: void fq_zech_sparse_vec_randtest(fq_zech_sparse_vec_t vec, flint_rand_t state, slong nnz, slong len, const fq_zech_ctx_t ctx)

    Makes vec a sparse vector with nnz nonzero entries uniformly distributed
    between 0 and len - 1, with individual entries generated by fq_zech_randtest


Output
--------------------------------------------------------------------------------

.. function:: void fq_zech_sparse_vec_print_pretty(const fq_zech_sparse_vec_t vec, slong ioff, slong maxi, const fq_zech_ctx_t ctx)

    Prints the vector of given length to ``stdout`` in a human-readable format


Arithmetic
--------------------------------------------------------------------------------

.. function:: void fq_zech_sparse_vec_neg(fq_zech_sparse_vec_t v, const fq_zech_sparse_vec_t u, const fq_zech_ctx_t ctx)

    Sets ``v`` to the negation of ``u``

.. function:: void fq_zech_sparse_vec_scalar_mul_fq_zech(fq_zech_sparse_vec_t v, const fq_zech_sparse_vec_t u, const fq_zech_t c, const fq_zech_ctx_t ctx)

    Sets ``v`` to the scalar multiple of ``u`` by ``c``

.. function:: void fq_zech_sparse_vec_add(fq_zech_sparse_vec_t w, const fq_zech_sparse_vec_t u, const fq_zech_sparse_vec_t v, const fq_zech_ctx_t ctx)

    Sets ``w`` to the sum of ``u`` and ``v``

.. function:: void fq_zech_sparse_vec_sub(fq_zech_sparse_vec_t w, const fq_zech_sparse_vec_t u, const fq_zech_sparse_vec_t v, const fq_zech_ctx_t ctx)

    Sets ``w`` to the difference of ``u`` and ``v``

.. function:: void fq_zech_sparse_vec_scalar_addmul_fq_zech(fq_zech_sparse_vec_t w, const fq_zech_sparse_vec_t u, const fq_zech_sparse_vec_t v, const fq_zech_t c, const fq_zech_ctx_t ctx)

    Sets ``w`` to the sum of ``u`` and ``c` times ``v``

.. function:: void fq_zech_sparse_vec_scalar_addmul_fq_zech(fq_zech_sparse_vec_t w, const fq_zech_sparse_vec_t u, const fq_zech_sparse_vec_t v, const fq_zech_t c, const fq_zech_ctx_t ctx)

    Sets ``w`` to the difference of ``u`` and ``c` times ``v``

.. function:: void fq_zech_sparse_vec_dot(fq_zech_t ret, const fq_zech_sparse_vec_t u, const fq_zech_sparse_vec_t v, const fq_zech_ctx_t ctx)

    Sets ``ret`` to the dot product of ``u`` and ``v``

.. function:: void fq_zech_sparse_vec_dot_dense(fq_zech_t ret, const fq_zech_sparse_vec_t u, const fq_zech_struct * v, const fq_zech_ctx_t ctx)

    Sets ``ret`` to the dot product of (``u``, ``v``)
