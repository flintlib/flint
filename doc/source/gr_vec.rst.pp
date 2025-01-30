include(`macros.m4')dnl
.. _gr-vec:

**gr_vec.h** -- vectors over generic rings
===============================================================================

Types and basic operations
--------------------------------------------------------------------------------

.. type:: gr_vec_struct

.. type:: gr_vec_t

.. function:: void gr_vec_init(gr_vec_t vec, slong len, gr_ctx_t ctx)

    Initializes *vec* to a vector of length *len* with elements
    in the ring *ctx*. The length must be nonnegative.
    All entries are set to zero.

.. function:: void gr_vec_clear(gr_vec_t vec, gr_ctx_t ctx)

    Clears the vector *vec*.

.. macro:: GR_VEC_ENTRY(vec, i, sz)

    Macro to access the *i*-th element in the vector *vec*,
    indexed from zero, assuming that entries have size ``sz``.
    The index must be in bounds.

.. function:: gr_ptr gr_vec_entry_ptr(gr_vec_t vec, slong i, gr_ctx_t ctx)

    Returns a pointer to the *i*-th element in the vector *vec*,
    indexed from zero. The index must be in bounds.

.. function:: slong gr_vec_length(const gr_vec_t vec, gr_ctx_t ctx)

    Returns the length of the vector *vec*.

.. function:: void gr_vec_fit_length(gr_vec_t vec, slong len, gr_ctx_t ctx)

    Allocates space for at least *len* elements in the vector *vec*.
    This does not change the size of the vector.

.. function:: void gr_vec_set_length(gr_vec_t vec, slong len, gr_ctx_t ctx)

    Resizes the vector to length *len*, which must be nonnegative.
    The vector will be extended with zeros.

.. function:: int gr_vec_set(gr_vec_t res, const gr_vec_t src, gr_ctx_t ctx)

    Sets *res* to a copy of the vector *src*.

.. function:: int gr_vec_append(gr_vec_t vec, gr_srcptr x, gr_ctx_t ctx)

    Appends the element *x* to the end of vector *vec*.

.. function:: int _gr_vec_write(gr_stream_t out, gr_srcptr vec, slong len, gr_ctx_t ctx)
              int gr_vec_write(gr_stream_t out, const gr_vec_t vec, gr_ctx_t ctx)
              int gr_vec_print(const gr_vec_t vec, gr_ctx_t ctx)

.. macro:: GR_ENTRY(vec, i, size)

    Macro to access the *i*-th entry of a ``gr_ptr`` or ``gr_srcptr``
    vector *vec*, where each element is ``size`` bytes.

.. function:: void _gr_vec_init(gr_ptr vec, slong len, gr_ctx_t ctx)

    Initialize *len* elements of *vec* to the value 0.
    The pointer *vec* must already refer to allocated memory.

.. function:: void _gr_vec_clear(gr_ptr vec, slong len, gr_ctx_t ctx)

    Clears *len* elements of *vec*.
    This frees memory allocated by individual elements, but
    does not free the memory allocated by *vec* itself.

.. function:: void _gr_vec_swap(gr_ptr vec1, gr_ptr vec2, slong len, gr_ctx_t ctx)

    Swap the entries of *vec1* and *vec2*.

.. function:: int _gr_vec_randtest(gr_ptr res, flint_rand_t state, slong len, gr_ctx_t ctx)

.. function:: int _gr_vec_set(gr_ptr res, gr_srcptr src, slong len, gr_ctx_t ctx)

.. function:: truth_t _gr_vec_equal(gr_srcptr vec1, gr_srcptr vec2, slong len, gr_ctx_t ctx)

.. function:: int _gr_vec_zero(gr_ptr vec, slong len, gr_ctx_t ctx)

.. function:: truth_t _gr_vec_is_zero(gr_srcptr vec, slong len, gr_ctx_t ctx)

.. function:: int _gr_vec_normalise(slong * res, gr_srcptr vec, slong len, gr_ctx_t ctx)

.. function:: slong _gr_vec_normalise_weak(gr_srcptr vec, slong len, gr_ctx_t ctx)


Arithmetic
--------------------------------------------------------------------------------

.. function:: int _gr_vec_neg(gr_ptr res, gr_srcptr src, slong len, gr_ctx_t ctx)

.. function:: int _gr_vec_add(gr_ptr res, gr_srcptr src1, gr_srcptr src2, slong len, gr_ctx_t ctx)
              int _gr_vec_sub(gr_ptr res, gr_srcptr src1, gr_srcptr src2, slong len, gr_ctx_t ctx)
              int _gr_vec_mul(gr_ptr res, gr_srcptr src1, gr_srcptr src2, slong len, gr_ctx_t ctx)
              int _gr_vec_div(gr_ptr res, gr_srcptr src1, gr_srcptr src2, slong len, gr_ctx_t ctx)
              int _gr_vec_divexact(gr_ptr res, gr_srcptr src1, gr_srcptr src2, slong len, gr_ctx_t ctx)
              int _gr_vec_pow(gr_ptr res, gr_srcptr src1, gr_srcptr src2, slong len, gr_ctx_t ctx)

    Binary operations applied elementwise.

.. function:: int _gr_vec_add_scalar(gr_ptr vec1, gr_srcptr vec2, slong len, gr_srcptr c, gr_ctx_t ctx)
              int _gr_vec_sub_scalar(gr_ptr vec1, gr_srcptr vec2, slong len, gr_srcptr c, gr_ctx_t ctx)
              int _gr_vec_mul_scalar(gr_ptr vec1, gr_srcptr vec2, slong len, gr_srcptr c, gr_ctx_t ctx)
              int _gr_vec_div_scalar(gr_ptr vec1, gr_srcptr vec2, slong len, gr_srcptr c, gr_ctx_t ctx)
              int _gr_vec_divexact_scalar(gr_ptr vec1, gr_srcptr vec2, slong len, gr_srcptr c, gr_ctx_t ctx)
              int _gr_vec_pow_scalar(gr_ptr vec1, gr_srcptr vec2, slong len, gr_srcptr c, gr_ctx_t ctx)
              int _gr_scalar_add_vec(gr_ptr vec1, gr_srcptr c, gr_srcptr vec2, slong len, gr_ctx_t ctx)
              int _gr_scalar_sub_vec(gr_ptr vec1, gr_srcptr c, gr_srcptr vec2, slong len, gr_ctx_t ctx)
              int _gr_scalar_mul_vec(gr_ptr vec1, gr_srcptr c, gr_srcptr vec2, slong len, gr_ctx_t ctx)
              int _gr_scalar_div_vec(gr_ptr vec1, gr_srcptr c, gr_srcptr vec2, slong len, gr_ctx_t ctx)
              int _gr_scalar_divexact_vec(gr_ptr vec1, gr_srcptr c, gr_srcptr vec2, slong len, gr_ctx_t ctx)
              int _gr_scalar_pow_vec(gr_ptr vec1, gr_srcptr c, gr_srcptr vec2, slong len, gr_ctx_t ctx)

    Binary operations applied elementwise with a fixed scalar operand.

.. function:: int _gr_vec_add_other(gr_ptr vec1, gr_srcptr vec2, gr_srcptr vec3, gr_ctx_t ctx3, slong len, gr_ctx_t ctx)
              int _gr_vec_sub_other(gr_ptr vec1, gr_srcptr vec2, gr_srcptr vec3, gr_ctx_t ctx3, slong len, gr_ctx_t ctx)
              int _gr_vec_mul_other(gr_ptr vec1, gr_srcptr vec2, gr_srcptr vec3, gr_ctx_t ctx3, slong len, gr_ctx_t ctx)
              int _gr_vec_div_other(gr_ptr vec1, gr_srcptr vec2, gr_srcptr vec3, gr_ctx_t ctx3, slong len, gr_ctx_t ctx)
              int _gr_vec_divexact_other(gr_ptr vec1, gr_srcptr vec2, gr_srcptr vec3, gr_ctx_t ctx3, slong len, gr_ctx_t ctx)
              int _gr_vec_pow_other(gr_ptr vec1, gr_srcptr vec2, gr_srcptr vec3, gr_ctx_t ctx3, slong len, gr_ctx_t ctx)
              int _gr_other_add_vec(gr_ptr vec1, gr_srcptr vec2, gr_ctx_t ctx2, gr_srcptr vec3, slong len, gr_ctx_t ctx)
              int _gr_other_sub_vec(gr_ptr vec1, gr_srcptr vec2, gr_ctx_t ctx2, gr_srcptr vec3, slong len, gr_ctx_t ctx)
              int _gr_other_mul_vec(gr_ptr vec1, gr_srcptr vec2, gr_ctx_t ctx2, gr_srcptr vec3, slong len, gr_ctx_t ctx)
              int _gr_other_div_vec(gr_ptr vec1, gr_srcptr vec2, gr_ctx_t ctx2, gr_srcptr vec3, slong len, gr_ctx_t ctx)
              int _gr_other_divexact_vec(gr_ptr vec1, gr_srcptr vec2, gr_ctx_t ctx2, gr_srcptr vec3, slong len, gr_ctx_t ctx)
              int _gr_other_pow_vec(gr_ptr vec1, gr_srcptr vec2, gr_ctx_t ctx2, gr_srcptr vec3, slong len, gr_ctx_t ctx)

    Binary operations applied elementwise, allowing a different type for one of the vectors.

.. function:: int _gr_vec_add_scalar_other(gr_ptr vec1, gr_srcptr vec2, slong len, gr_srcptr c, gr_ctx_t cctx, gr_ctx_t ctx)
              int _gr_vec_sub_scalar_other(gr_ptr vec1, gr_srcptr vec2, slong len, gr_srcptr c, gr_ctx_t cctx, gr_ctx_t ctx)
              int _gr_vec_mul_scalar_other(gr_ptr vec1, gr_srcptr vec2, slong len, gr_srcptr c, gr_ctx_t cctx, gr_ctx_t ctx)
              int _gr_vec_div_scalar_other(gr_ptr vec1, gr_srcptr vec2, slong len, gr_srcptr c, gr_ctx_t cctx, gr_ctx_t ctx)
              int _gr_vec_divexact_scalar_other(gr_ptr vec1, gr_srcptr vec2, slong len, gr_srcptr c, gr_ctx_t cctx, gr_ctx_t ctx)
              int _gr_vec_pow_scalar_other(gr_ptr vec1, gr_srcptr vec2, slong len, gr_srcptr c, gr_ctx_t cctx, gr_ctx_t ctx)
              int _gr_scalar_other_add_vec(gr_ptr vec1, gr_srcptr c, gr_ctx_t cctx, gr_srcptr vec2, slong len, gr_ctx_t ctx)
              int _gr_scalar_other_sub_vec(gr_ptr vec1, gr_srcptr c, gr_ctx_t cctx, gr_srcptr vec2, slong len, gr_ctx_t ctx)
              int _gr_scalar_other_mul_vec(gr_ptr vec1, gr_srcptr c, gr_ctx_t cctx, gr_srcptr vec2, slong len, gr_ctx_t ctx)
              int _gr_scalar_other_div_vec(gr_ptr vec1, gr_srcptr c, gr_ctx_t cctx, gr_srcptr vec2, slong len, gr_ctx_t ctx)
              int _gr_scalar_other_divexact_vec(gr_ptr vec1, gr_srcptr c, gr_ctx_t cctx, gr_srcptr vec2, slong len, gr_ctx_t ctx)
              int _gr_scalar_other_pow_vec(gr_ptr vec1, gr_srcptr c, gr_ctx_t cctx, gr_srcptr vec2, slong len, gr_ctx_t ctx)
              int _gr_vec_add_scalar_si(gr_ptr vec1, gr_srcptr vec2, slong len, slong c, gr_ctx_t ctx)
              int _gr_vec_sub_scalar_si(gr_ptr vec1, gr_srcptr vec2, slong len, slong c, gr_ctx_t ctx)
              int _gr_vec_mul_scalar_si(gr_ptr vec1, gr_srcptr vec2, slong len, slong c, gr_ctx_t ctx)
              int _gr_vec_div_scalar_si(gr_ptr vec1, gr_srcptr vec2, slong len, slong c, gr_ctx_t ctx)
              int _gr_vec_divexact_scalar_si(gr_ptr vec1, gr_srcptr vec2, slong len, slong c, gr_ctx_t ctx)
              int _gr_vec_pow_scalar_si(gr_ptr vec1, gr_srcptr vec2, slong len, slong c, gr_ctx_t ctx)
              int _gr_vec_add_scalar_ui(gr_ptr vec1, gr_srcptr vec2, slong len, ulong c, gr_ctx_t ctx)
              int _gr_vec_sub_scalar_ui(gr_ptr vec1, gr_srcptr vec2, slong len, ulong c, gr_ctx_t ctx)
              int _gr_vec_mul_scalar_ui(gr_ptr vec1, gr_srcptr vec2, slong len, ulong c, gr_ctx_t ctx)
              int _gr_vec_div_scalar_ui(gr_ptr vec1, gr_srcptr vec2, slong len, ulong c, gr_ctx_t ctx)
              int _gr_vec_divexact_scalar_ui(gr_ptr vec1, gr_srcptr vec2, slong len, ulong c, gr_ctx_t ctx)
              int _gr_vec_pow_scalar_ui(gr_ptr vec1, gr_srcptr vec2, slong len, ulong c, gr_ctx_t ctx)
              int _gr_vec_add_scalar_fmpz(gr_ptr vec1, gr_srcptr vec2, slong len, const fmpz_t c, gr_ctx_t ctx)
              int _gr_vec_sub_scalar_fmpz(gr_ptr vec1, gr_srcptr vec2, slong len, const fmpz_t c, gr_ctx_t ctx)
              int _gr_vec_mul_scalar_fmpz(gr_ptr vec1, gr_srcptr vec2, slong len, const fmpz_t c, gr_ctx_t ctx)
              int _gr_vec_div_scalar_fmpz(gr_ptr vec1, gr_srcptr vec2, slong len, const fmpz_t c, gr_ctx_t ctx)
              int _gr_vec_divexact_scalar_fmpz(gr_ptr vec1, gr_srcptr vec2, slong len, const fmpz_t c, gr_ctx_t ctx)
              int _gr_vec_pow_scalar_fmpz(gr_ptr vec1, gr_srcptr vec2, slong len, const fmpz_t c, gr_ctx_t ctx)
              int _gr_vec_add_scalar_fmpq(gr_ptr vec1, gr_srcptr vec2, slong len, const fmpq_t c, gr_ctx_t ctx)
              int _gr_vec_sub_scalar_fmpq(gr_ptr vec1, gr_srcptr vec2, slong len, const fmpq_t c, gr_ctx_t ctx)
              int _gr_vec_mul_scalar_fmpq(gr_ptr vec1, gr_srcptr vec2, slong len, const fmpq_t c, gr_ctx_t ctx)
              int _gr_vec_div_scalar_fmpq(gr_ptr vec1, gr_srcptr vec2, slong len, const fmpq_t c, gr_ctx_t ctx)
              int _gr_vec_divexact_scalar_fmpq(gr_ptr vec1, gr_srcptr vec2, slong len, const fmpq_t c, gr_ctx_t ctx)
              int _gr_vec_pow_scalar_fmpq(gr_ptr vec1, gr_srcptr vec2, slong len, const fmpq_t c, gr_ctx_t ctx)

    Binary operations applied elementwise with a fixed scalar operand, allowing a different type
    for the scalar.

.. function:: int _gr_vec_addmul_scalar(gr_ptr vec1, gr_srcptr vec2, slong len, gr_srcptr c, gr_ctx_t ctx)
              int _gr_vec_submul_scalar(gr_ptr vec1, gr_srcptr vec2, slong len, gr_srcptr c, gr_ctx_t ctx)
              int _gr_vec_addmul_scalar_si(gr_ptr vec1, gr_srcptr vec2, slong len, slong c, gr_ctx_t ctx)
              int _gr_vec_submul_scalar_si(gr_ptr vec1, gr_srcptr vec2, slong len, slong c, gr_ctx_t ctx)

.. function:: int _gr_vec_mul_scalar_2exp_si(gr_ptr res, gr_srcptr vec, slong len, slong c, gr_ctx_t ctx)

Sums and products
--------------------------------------------------------------------------------

.. function:: int _gr_vec_sum(gr_ptr res, gr_srcptr vec, slong len, gr_ctx_t ctx)

.. function:: int _gr_vec_product(gr_ptr res, gr_srcptr vec, slong len, gr_ctx_t ctx)

Dot products
--------------------------------------------------------------------------------

.. function:: int _gr_vec_dot(gr_ptr res, gr_srcptr initial, int subtract, gr_srcptr vec1, gr_srcptr vec2, slong len, gr_ctx_t ctx)
              int _gr_vec_dot_si(gr_ptr res, gr_srcptr initial, int subtract, gr_srcptr vec1, const slong * vec2, slong len, gr_ctx_t ctx)
              int _gr_vec_dot_ui(gr_ptr res, gr_srcptr initial, int subtract, gr_srcptr vec1, const ulong * vec2, slong len, gr_ctx_t ctx)
              int _gr_vec_dot_fmpz(gr_ptr res, gr_srcptr initial, int subtract, gr_srcptr vec1, const fmpz * vec2, slong len, gr_ctx_t ctx)

    Sets *res* to `c \pm \sum_{i=0}^{n-1} a_i b_i`.

.. function:: int _gr_vec_dot_rev(gr_ptr res, gr_srcptr initial, int subtract, gr_srcptr vec1, gr_srcptr vec2, slong len, gr_ctx_t ctx)

    Sets *res* to `c \pm \sum_{i=0}^{n-1} a_i b_{n-1-i}`.

Other functions
--------------------------------------------------------------------------------

.. function:: int _gr_vec_step(gr_ptr vec, gr_srcptr start, gr_srcptr step, slong len, gr_ctx_t ctx)

.. function:: int _gr_vec_reciprocals(gr_ptr res, slong len, gr_ctx_t ctx)

    Sets *res* to the vector of reciprocals of the positive integers 1, 2, ... up to *len* inclusive.

.. function:: int _gr_vec_set_powers(gr_ptr res, gr_srcptr x, slong len, gr_ctx_t ctx)


.. raw:: latex

    \newpage
