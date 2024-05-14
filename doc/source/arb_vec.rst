.. _arb-vec:

**arb_vec.h** -- vectors of real numbers
==================================================================================================

Memory management
--------------------------------------------------------------------------------

.. function:: arb_ptr _arb_vec_init(slong n)

    Returns a pointer to an array of *n* initialized :type:`arb_struct`
    entries.

.. function:: void _arb_vec_clear(arb_ptr v, slong n)

    Clears an array of *n* initialized :type:`arb_struct` entries.

.. function:: void _arb_vec_swap(arb_ptr vec1, arb_ptr vec2, slong len)

    Swaps the entries of *vec1* and *vec2*.

.. function:: slong _arb_vec_allocated_bytes(arb_srcptr vec, slong len)

    Returns the total number of bytes allocated for this vector, i.e. the
    space taken up by the vector itself plus the sum of the internal heap
    allocation sizes for all its member elements.

.. function:: double _arb_vec_estimate_allocated_bytes(slong len, slong prec)

    Estimates the number of bytes that need to be allocated for a vector of
    *len* elements with *prec* bits of precision, including the space for
    internal limb data.
    This function returns a *double* to avoid overflow issues when both
    *len* and *prec* are large.

    This is only an approximation of the physical memory that will be used
    by an actual vector. In practice, the space varies with the content
    of the numbers; for example, zeros and small integers require no
    internal heap allocation even if the precision is huge.
    The estimate assumes that exponents will not be bignums.
    The actual amount may also be higher or lower due to overhead in the
    memory allocator or overcommitment by the operating system.

.. function:: void _arb_vec_trim(arb_ptr res, arb_srcptr vec, slong len)

    Applies :func:`arb_trim` elementwise.

.. function:: slong _arb_vec_bits(arb_srcptr x, slong len)

    Returns the maximum of :func:`arb_bits` for all entries in *vec*.

Comparisons
--------------------------------------------------------------------------------

.. function:: int _arb_vec_is_zero(arb_srcptr vec, slong len)

    Returns nonzero iff all entries in *x* are zero.

.. function:: int _arb_vec_is_finite(arb_srcptr x, slong len)

    Returns nonzero iff all entries in *x* certainly are finite.

.. function:: int _arb_vec_equal(arb_srcptr vec1, arb_srcptr vec2, slong len)

    Returns nonzero iff *vec1* and *vec2* are equal in the sense of
    :func:`arb_equal`, i.e. have both the same midpoint and radius element-wise.

.. function:: int _arb_vec_overlaps(arb_srcptr vec1, arb_srcptr vec2, slong len)

    Returns nonzero iff *vec1* overlaps *vec2* element-wise.

.. function:: int _arb_vec_contains(arb_srcptr vec1, arb_srcptr vec2, slong len)

    Returns nonzero iff *vec1* contains *vec2* element-wise.

Assignments and conversions
--------------------------------------------------------------------------------

.. function:: void _arb_vec_set(arb_ptr res, arb_srcptr vec, slong len)

    Sets *res* to a copy of *vec*.

.. function:: void _arb_vec_set_round(arb_ptr res, arb_srcptr vec, slong len, slong prec)

    Sets *res* to a copy of *vec*, rounding each entry to *prec* bits.

.. function:: void _arb_vec_zero(arb_ptr vec, slong n)

    Sets all entries in *vec* to zero.

.. function:: int _arb_vec_get_unique_fmpz_vec(fmpz * res,  arb_srcptr vec, slong len)

    Calls :func:`arb_get_unique_fmpz` element-wise and returns nonzero if
    all entries can be rounded uniquely to integers. If any entry in *vec*
    cannot be rounded uniquely to an integer, returns zero.

Arithmetic
--------------------------------------------------------------------------------

.. function:: void _arb_vec_add(arb_ptr C, arb_srcptr A, arb_srcptr B, slong n, slong prec)
              void _arb_vec_sub(arb_ptr C, arb_srcptr A, arb_srcptr B, slong n, slong prec)

    Performs `C = A \pm B` with *prec* bits of precision.

.. function:: void _arb_vec_neg(arb_ptr B, arb_srcptr A, slong n)

    Performs `B = -A`. Precision is preserved.

Scalar arithmetic
--------------------------------------------------------------------------------

.. function:: void _arb_vec_scalar_mul(arb_ptr res, arb_srcptr vec, slong len, const arb_t c, slong prec)
              void _arb_vec_scalar_div(arb_ptr res, arb_srcptr vec, slong len, const arb_t c, slong prec)
              void _arb_vec_scalar_mul_fmpz(arb_ptr res, arb_srcptr vec, slong len, const fmpz_t c, slong prec)
              void _arb_vec_scalar_mul_2exp_si(arb_ptr res, arb_srcptr src, slong len, slong c)
              void _arb_vec_scalar_addmul(arb_ptr res, arb_srcptr vec, slong len, const arb_t c, slong prec)

   Performs the respective scalar operation element-wise.

Error arithmetic
--------------------------------------------------------------------------------

.. function:: void _arb_vec_get_mag(mag_t bound, arb_srcptr vec, slong len)

    Sets *bound* to an upper bound for the entries in *vec*.

.. function:: void _arb_vec_add_error_arf_vec(arb_ptr res, arf_srcptr err, slong len)
              void _arb_vec_add_error_mag_vec(arb_ptr res, mag_srcptr err, slong len)

    Adds the magnitude of each entry in *err* to the radius of the corresponding
    entry in *res*.

.. function:: void _arb_vec_indeterminate(arb_ptr vec, slong len)

    Applies :func:`arb_indeterminate` element-wise.

Miscellaneous
--------------------------------------------------------------------------------

.. function:: void _arb_vec_set_powers(arb_ptr xs, const arb_t x, slong len, slong prec)

    Sets *xs* to the powers `1, x, x^2, \ldots, x^{len-1}`.

Input and output
--------------------------------------------------------------------------------

.. function:: void _arb_vec_printn(arb_srcptr vec, slong len, slong digits, ulong flags)
              void _arb_vec_printd(arb_srcptr vec, slong len, slong ndigits)

    Prints *vec* in decimal using :func:`arb_printn` or :func:`arb_printd` on each entry.
