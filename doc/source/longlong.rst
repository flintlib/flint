.. _longlong:

**longlong.h** -- support functions for multi-word arithmetic
===============================================================================

Bit manipulation
-------------------------------------------------------------------------------

.. macro:: flint_clz(x)

    Returns the number of zero-bits from the msb to the first non-zero bit in
    the limb `x`.  This is the number of steps `x` needs to be shifted left to
    set the most significant bit in `x`. If `x` is zero then the return value is
    undefined.

.. macro:: flint_ctz(x)

    As for ``flint_clz()``, but counts from the least significant end. If `x` is
    zero then the return value is undefined.

.. function:: flint_bitcnt_t FLINT_BIT_COUNT(ulong x)

    Returns the number of binary bits required to represent *x*. If *x* is zero
    it returns *0*. This is an inline-function only.

.. macro:: FLINT_FLOG2(x)
           FLINT_CLOG2(x)

    For `x \ge 1`, it returns `\lfloor \log_2 x \rfloor`
    and `\lceil \log_2 x \rceil`, respectively.

Addition and subtraction
-------------------------------------------------------------------------------

.. note::

    When aliasing inputs with outputs in these addition and subtraction macros,
    make sure to have `s_{i}` aliased with `a_{i}` for addition macros, and
    `d_{i}` aliased with `m_{i}` for optimal performance. Moreover, keep
    immediates (in other words, constants known to the compiler) in the `b_{i}`
    variables for addition and `s_{i}` for subtraction.

.. macro:: add_ssaaaa(s1, s0, a1, a0, b1, b0)

    Sets `s_1` and `s_0` according
    to `c B^2 + s_1 B + s_0 = (a_1 B + a_0) + (b_1 B + b_0)`,
    where `B = 2^{\mathtt{FLINT\_BITS}}` is the base, and `c` is the carry from
    the addition which is not stored anywhere.

.. macro:: add_sssaaaaaa(s2, s1, s0, a2, a1, a0, b2, b1, b0)

    Works like ``add_ssaaaa``, but for two three-limbed integers. Carry is lost.

.. macro:: sub_ddmmss(d1, d0, m1, m0, s1, s0)

    Sets `d_1` and `d_0` to the difference between the two-limbed
    integers `m_1 B + m_0` and `s_1 B + s_0`,
    where `B = 2^{\mathtt{FLINT\_BITS}}`. Borrow from the subtraction is not
    stored anywhere.

.. macro:: sub_dddmmmsss(d2, d1, d0, m2, m1, m0, s2, s1, s0)

    Works like ``sub_dddmmmsss``, but for two three-limbed integers. Borrow is
    lost.

Multiplication
-------------------------------------------------------------------------------

.. macro:: umul_ppmm(p1, p0, u, v)

    Computes `p_1 B + p0 = u v`, where `B = 2^{\mathtt{FLINT\_BITS}}`.

.. macro:: smul_ppmm(p1, p0, u, v)

    Works like ``umul_ppmm`` but for signed numbers.

Division
-------------------------------------------------------------------------------

.. macro:: udiv_qrnnd(q, r, n1, n0, d)

    Computes the non-negative integers `q` and `r` in `d q + r = n_1 B + n_0`,
    where `B = 2^{\mathtt{FLINT\_BITS}}`. Assumes that `n_1 < d`.

.. macro:: sdiv_qrnnd(quotient, remainder, high_numerator, low_numerator, denominator)

    Works like ``udiv_qrnnd``, but for signed numbers.

.. macro:: udiv_qrnnd_preinv(q, r, n1, n0, d, di)

    Works like ``udiv_qrnnd``, but takes a precomputed inverse ``di`` as 
    computed by :func:`n_preinvert_limb_prenorm`. This function assumes ``d``
    is normalised, i.e. with most significant bit set.
