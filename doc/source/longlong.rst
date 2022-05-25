.. _longlong:

**longlong.h** -- support functions for multi-word arithmetic
===============================================================================


Auxiliary asm macros
--------------------------------------------------------------------------------


.. macro:: umul_ppmm(high_prod, low_prod, multipler, multiplicand)

    Multiplies two single limb integers ``MULTIPLER`` and 
    ``MULTIPLICAND``, and generates a two limb product in 
    ``HIGH_PROD`` and ``LOW_PROD``.

.. macro:: smul_ppmm(high_prod, low_prod, multipler, multiplicand)

    As for ``umul_ppmm()`` but the numbers are signed.

.. macro:: udiv_qrnnd(quotient, remainder, high_numerator, low_numerator, denominator)

    Divides an unsigned integer, composed by the limb integers 
    ``HIGH_NUMERATOR`` and\\ ``LOW_NUMERATOR``, by ``DENOMINATOR`` 
    and places the quotient in ``QUOTIENT`` and the remainder in 
    ``REMAINDER``.  ``HIGH_NUMERATOR`` must be less than 
    ``DENOMINATOR`` for correct operation. 

.. macro:: sdiv_qrnnd(quotient, remainder, high_numerator, low_numerator, denominator)

    As for ``udiv_qrnnd()`` but the numbers are signed.  The quotient is 
    rounded towards `0`. Note that as the quotient is signed it must lie in 
    the range `[-2^63, 2^63)`.

.. macro:: count_leading_zeros(count, x)

    Counts the number of zero-bits from the msb to the first non-zero bit 
    in the limb ``x``.  This is the number of steps ``x`` needs to 
    be shifted left to set the msb. If ``x`` is `0` then count is 
    undefined.

.. macro:: count_trailing_zeros(count, x)

    As for ``count_leading_zeros()``, but counts from the least 
    significant end. If ``x`` is zero then count is undefined.

.. macro:: add_ssaaaa(high_sum, low_sum, high_addend_1, low_addend_1, high_addend_2, low_addend_2)

    Adds two limb integers, composed by ``HIGH_ADDEND_1`` and 
    ``LOW_ADDEND_1``, and\\ ``HIGH_ADDEND_2`` and ``LOW_ADDEND_2``, 
    respectively.  The result is placed in ``HIGH_SUM`` and 
    ``LOW_SUM``.  Overflow, i.e.\ carry out, is not stored anywhere, 
    and is lost.

.. macro:: add_sssaaaaaa(high_sum, mid_sum, low_sum, high_addend_1, mid_addend_1, low_addend_1, high_addend_2, mid_addend_2, low_addend_2)

    Adds two three limb integers. Carry out is lost.

.. macro:: sub_ddmmss(high_difference, low_difference, high_minuend, low_minuend, high_subtrahend, low_subtrahend)

    Subtracts two limb integers, composed by ``HIGH_MINUEND_1`` and 
    ``LOW_MINUEND_1``, and ``HIGH_SUBTRAHEND_2`` and 
    ``LOW_SUBTRAHEND_2``, respectively.  The result is placed in\\ 
    ``HIGH_DIFFERENCE`` and ``LOW_DIFFERENCE``.  Overflow, i.e.\ 
    carry out is not stored anywhere, and is lost.

.. macro:: sub_dddmmmsss(high_diff, mid_diff, low_diff, high_minuend_1, mid_minuend_1, low_minuend_1, high_subtrahend_2, mid_subtrahend_2, low_subtrahend_2)

    Subtracts two three limb integers. Borrow out is lost.

.. macro:: byte_swap(x)

    Swap the order of the bytes in the word `x`, i.e. most significant byte
    becomes least significant byte, etc.

.. macro:: invert_limb(invxl, xl)

    Deprecated: see :func:`n_preinvert_limb_prenorm`.

.. macro:: udiv_qrnnd_preinv(q, r, nh, nl, d, di)

    As for ``udiv_qrnnd()`` but takes a precomputed inverse ``di`` as 
    computed by ``invert_limb()``. The algorithm, in terms of the theorem 
    above, is::

        nadj = n1*(d-B/2) + n0
        xh, xl = (n2+n1)*(m-B)
        xh, xl += nadj + n2*B ( xh, xl = n2*B + (n2+n1)*(m-B) + n1*(d-B/2) + n0 )
        _q1 = B - xh - 1
        xh, xl = _q1*d + nh, nl - B*d = nh, nl - q1*d - d so that xh = 0 or -1
        r = xl + xh*d where xh is 0 if q1 is off by 1, otherwise -1
        q = xh - _q1 = xh + 1 + n2

