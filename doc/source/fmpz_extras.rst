.. _fmpz_extras:

**fmpz_extras.h** -- extra methods for FLINT integers
===============================================================================

This module implements a few utility methods for the FLINT
multiprecision integer type (*fmpz_t*). It is mainly intended for internal use.

Memory-related methods
-------------------------------------------------------------------------------

.. function:: slong fmpz_allocated_bytes(const fmpz_t x)

    Returns the total number of bytes heap-allocated internally by this object.
    The count excludes the size of the structure itself. Add
    ``sizeof(fmpz)`` to get the size of the object as a whole.

Convenience methods
-------------------------------------------------------------------------------

.. function:: void fmpz_add_si(fmpz_t z, const fmpz_t x, slong y)

.. function:: void fmpz_sub_si(fmpz_t z, const fmpz_t x, slong y)

    Sets *z* to the sum (respectively difference) of *x* and *y*.

.. function:: void fmpz_adiv_q_2exp(fmpz_t z, const fmpz_t x, flint_bitcnt_t exp)

    Sets *z* to `x / 2^{exp}`, rounded away from zero.

.. function:: void fmpz_ui_mul_ui(fmpz_t x, ulong a, ulong b)

    Sets *x* to *a* times *b*.

.. function:: void fmpz_ui_pow_ui(fmpz_t x, ulong b, ulong e)

    Sets *x* to *b* raised to the power *e*.

.. function:: void fmpz_max(fmpz_t z, const fmpz_t x, const fmpz_t y)

.. function:: void fmpz_min(fmpz_t z, const fmpz_t x, const fmpz_t y)

    Sets *z* to the maximum (respectively minimum) of *x* and *y*.

Inlined arithmetic
-------------------------------------------------------------------------------

The *fmpz_t* bignum type uses an immediate representation for small
integers, specifically when the absolute value is at most `2^{62}-1` (on
64-bit machines) or `2^{30}-1` (on 32-bit machines).
The following methods completely inline the case
where all operands (and possibly some intermediate values in the calculation)
are known to be small.
This is faster in code where all values *almost certainly will be much
smaller than a full word*. In particular, these methods are used within
Arb for manipulating exponents of floating-point numbers.
Inlining slows down the general case, and increases code size,
so these methods should not be used gratuitously.

.. function:: void fmpz_add_inline(fmpz_t z, const fmpz_t x, const fmpz_t y)

.. function:: void fmpz_add_si_inline(fmpz_t z, const fmpz_t x, slong y)

.. function:: void fmpz_add_ui_inline(fmpz_t z, const fmpz_t x, ulong y)

    Sets *z* to the sum of *x* and *y*.

.. function:: void fmpz_sub_si_inline(fmpz_t z, const fmpz_t x, slong y)

    Sets *z* to the difference of *x* and *y*.

.. function:: void fmpz_add2_fmpz_si_inline(fmpz_t z, const fmpz_t x, const fmpz_t y, slong c)

    Sets *z* to the sum of *x*, *y*, and *c*.

.. function:: mp_size_t _fmpz_size(const fmpz_t x)

    Returns the number of limbs required to represent *x*.

.. function:: slong _fmpz_sub_small(const fmpz_t x, const fmpz_t y)

    Computes the difference of *x* and *y* and returns the result as
    an *slong*. The result is clamped between -*WORD_MAX* and *WORD_MAX*,
    i.e. between `\pm (2^{63}-1)` inclusive on a 64-bit machine.

.. function:: slong _fmpz_set_si_small(fmpz_t x, slong v)

    Sets *x* to the integer *v* which is required to be a value
    between *COEFF_MIN* and *COEFF_MAX* so that promotion to
    a bignum cannot occur.

Low-level conversions
-------------------------------------------------------------------------------

.. function:: void fmpz_set_mpn_large(fmpz_t z, mp_srcptr src, mp_size_t n, int negative)

    Sets *z* to the integer represented by the *n* limbs in the array *src*,
    or minus this value if *negative* is 1.
    Requires `n \ge 2` and that the top limb of *src* is nonzero.
    Note that *fmpz_set_ui*, *fmpz_neg_ui* can be used for single-limb integers.

.. macro:: FMPZ_GET_MPN_READONLY(zsign, zn, zptr, ztmp, zv)

    Given an *fmpz_t* *zv*, this macro sets *zptr* to a pointer to the limbs of *zv*,
    *zn* to the number of limbs, and *zsign* to a sign bit (0 if nonnegative,
    1 if negative). The variable *ztmp* must be a single *mp_limb_t*, which is
    used as a buffer. If *zv* is a small value, *zv* itself contains no limb
    array that *zptr* could point to, so the single limb is copied to *ztmp*
    and *zptr* is set to point to *ztmp*. The case where *zv*
    is zero is not handled specially, and *zn* is set to 1.

.. function:: void fmpz_lshift_mpn(fmpz_t z, mp_srcptr src, mp_size_t n, int negative, flint_bitcnt_t shift)

    Sets *z* to the integer represented by the *n* limbs in the array *src*,
    or minus this value if *negative* is 1, shifted left by *shift* bits.
    Requires `n \ge 1` and that the top limb of *src* is nonzero.

