/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "longlong.h"
#include "mpn_extras.h"
#include "mp_real.h"
#include "fp_small.h"

/*
    Floating-point quotients with a guard limb (see mpn_extras.h).

    With A = {a, an} / B^an and Y = {b, bn} / B^bn in [1/2, 1), the
    quotient A / Y lies in (1/2, 2); the shift s = [A >= Y] brings
    2^-s A / Y into [1/2, 1), and T = 2^-s (A / Y) B^(n+1) is what the
    n + 1 limbs approximate.

    The operands are truncated to m = n + 2 limbs (A' = A - da, Y' =
    Y - dy with 0 <= da, dy < B^-m): A'/Y' - A/Y = (dy A - da Y) / (Y Y')
    is below 2 B^-m / (1/4) = 8 B^-m in absolute value, i.e. 8 / B units
    of the guard limb. The shift is taken from the truncated operands,
    so that for the exact quotient 2^-s A / Y may fall slightly outside
    [1/2, 1) (by less than 8 B^-m) when the comparison of the truncated
    operands differs from that of the full ones; the result is clamped
    to [B^(n+1) / 2, B^(n+1)), which keeps the error bound.

    The quotient itself is then computed as
      * n = 1: two 3/2 divisions on the top two limbs of each operand
        (truncation to 2 limbs: 8 units), exact floor;
      * a single-limb divisor: exact floor by a chain of 2/1 divisions;
      * short divisors (FLINT_MPN_DIV_SMALL_SHAPE): exact floor by the
        register-based division;
      * otherwise the schoolbook approximate division
        (_flint_mpn_divapprox_basecase_preinv1) up to DIVHIGH_BASECASE_MAX
        limbs, flint_mpn_divapprox beyond: the floor or one more;
      * beyond DIVHIGH_BASECASE_MAX limbs, a numerator of at most n / 8
        limbs multiplies the approximate reciprocal (flint_mpn_invapprox)
        of the divisor: within 2 units.
    So |Q - T| < 2 + 8 B^-m B^(n+1) / ... <= 10 in all cases, within
    FLINT_MPN_DIVHIGH_GUARD_ERR_BOUND = 16.
*/

/* [A >= Y] for fractions of an resp. bn limbs (zero-extended) */
FLINT_FORCE_INLINE ulong
_fraction_ge(nn_srcptr a, slong an, nn_srcptr b, slong bn)
{
    slong m;
    int c;

    if (a[an - 1] != b[bn - 1])
        return a[an - 1] > b[bn - 1];

    m = FLINT_MIN(an, bn);
    c = mpn_cmp(a + an - m, b + bn - m, m);
    if (c != 0)
        return c > 0;
    if (bn > an)
        return flint_mpn_zero_p(b, bn - an);
    return 1;
}

/* the n + 1 limbs of q normalised into [B^(n+1)/2, B^(n+1)): clamped
   above when the quotient reached B^(n+1) (top limb q[n+1] != 0), and
   below when the top bit is clear (both only by the truncation of the
   operands or the error of the approximation) */
FLINT_FORCE_INLINE void
_divhigh_clamp(nn_ptr q, ulong top, slong n)
{
    if (top != 0)
    {
        flint_mpn_store(q, n + 1, UWORD_MAX);
    }
    else if (!(q[n] >> (FLINT_BITS - 1)))
    {
        flint_mpn_zero(q, n);
        q[n] = UWORD(1) << (FLINT_BITS - 1);
    }
}

/* (t, nn = bn + n + 1) = floor(A' B^(bn + n + 1) / 2^s) for the fraction
   A' of an <= n + 2 limbs */
FLINT_FORCE_INLINE void
_divhigh_numerator(nn_ptr t, nn_srcptr a, slong an, slong bn, slong n, ulong s)
{
    slong nn = bn + n + 1, off = nn - an, i;

    FLINT_ASSERT(off >= 0);
    _flint_mpn_zero_small(t, off);
    if (s == 0)
        _flint_mpn_copy_small(t + off, a, an);
    else
    {
        /* a >> 1 with the bit shifted out at the top of the padding */
        ulong lo = 0;
        for (i = an - 1; i >= 0; i--)
        {
            t[off + i] = (a[i] >> 1) | lo;
            lo = a[i] << (FLINT_BITS - 1);
        }
        if (off > 0)
            t[off - 1] = lo;
    }
}

/* up to this many quotient limbs, the schoolbook approximate division
   with the 3/2 inverse is used for divisors of three or more limbs (for
   these sizes it beats the divide-and-conquer and short divisions of
   flint_mpn_divapprox: measured crossover 56-80 on Skylake, 80 on Zen 3);
   the stack holds the scratch */
#define DIVHIGH_BASECASE_MAX 80
#define DIVHIGH_STACK_LIMBS (3 * DIVHIGH_BASECASE_MAX + 8)

/* numerators of at most n / 8 limbs (1 / y in particular) go through the
   approximate reciprocal of the divisor (measured: 0.75-0.96 of the time
   of the division for 1 / y from 100 to 15000 limbs) */
#define DIVHIGH_SHORT_NUMERATOR(an, n) (8 * (an) <= (n))

ulong
flint_mpn_divhigh(nn_ptr q, nn_srcptr a, slong an, nn_srcptr b, slong bn, slong n)
{
    slong m = n + 2, nn;
    ulong s, stack[DIVHIGH_STACK_LIMBS];
    nn_ptr t, qq;

    FLINT_ASSERT(n >= 1 && an >= 1 && bn >= 1);
    FLINT_ASSERT(a[an - 1] >> (FLINT_BITS - 1));
    FLINT_ASSERT(b[bn - 1] >> (FLINT_BITS - 1));

    /* truncate the operands to m limbs */
    if (an > m)
    {
        a += an - m;
        an = m;
    }
    if (bn > m)
    {
        b += bn - m;
        bn = m;
    }

    if (n == 1)
    {
        ulong a1 = a[an - 1], a0 = (an >= 2) ? a[an - 2] : 0;
        ulong b1 = b[bn - 1], b0 = (bn >= 2) ? b[bn - 2] : 0;
        ulong q1, q0;

        _flint_mpn_fdiv_2_kernel(&q1, &q0, &s, a1, a0, b1, b0);
        q[1] = q1;
        q[0] = q0;
        /* two 3/2 divisions of a normalised numerator: in range already */
        return s;
    }

    s = _fraction_ge(a, an, b, bn);
    nn = bn + n + 1;

    if (n > DIVHIGH_BASECASE_MAX && DIVHIGH_SHORT_NUMERATOR(an, n) && bn != 1 &&
        !FLINT_MPN_DIV_SMALL_SHAPE(nn, bn, 0))
    {
        /* a short numerator (1 / y in particular): the approximate
           reciprocal I of b times a, which spares the multiplication of
           the full-length quotient that a division performs (and uses
           the Newton inversion with its lower cutoffs). With X = B^N / b,
           N = n + 1 + bn, I is floor(X) or floor(X) + 1 (n + 3 limbs,
           below 2 B^(n+1) + 2); T' = 2^-s a B^N / b / B^an for the
           truncated operands, and Q = floor(a I / (2^s B^an)) differs from
           T' by less than a / B^an < 1 plus the floor: |Q - T| < 8 + 2. */
        nn_ptr U, I;
        slong in = n + 3, alloc;

        /* one buffer: the reciprocal (in limbs, at U or U + 2, with the
           scratch of the Newton iteration), then the product a I (in + an
           limbs) in place for a single limb (the reciprocal and small
           integers), whose top in limbs from limb an are the quotient */
        alloc = (an == 1) ? in + 3 : in + 2 + an + in;
        U = (alloc <= DIVHIGH_STACK_LIMBS) ? stack : flint_malloc(alloc * sizeof(ulong));
        I = _flint_mpn_invapprox_buf(U, b, bn, n + 1 + bn);
        if (an == 1 && a[0] == (UWORD(1) << (FLINT_BITS - 1)))
        {
            /* a I / B = I / 2 */
            mpn_rshift(I, I, in, 1 + s);
            qq = I;
        }
        else
        {
            if (an == 1)
            {
                I[in] = mpn_mul_1(I, I, in, a[0]);
                qq = I + 1;
            }
            else
            {
                flint_mpn_mul(U + in + 2, I, in, a, an);
                qq = U + in + 2 + an;
            }
            if (s)
                mpn_rshift(qq, qq, in, 1);
        }
        _divhigh_clamp(qq, qq[n + 1] | qq[n + 2], n);
        _flint_mpn_copy_small(q, qq, n + 1);
        if (U != stack)
            flint_free(U);
        return s;
    }

    /* the quotient has n + 2 limbs (the top one zero, or one when the
       truncated comparison went the other way than the exact one): the
       schoolbook division writes the low n + 1 to q directly and returns
       the top one, the other divisions need a buffer */
    if (nn + n + 2 <= DIVHIGH_STACK_LIMBS)
        t = stack;
    else
        t = flint_malloc((nn + n + 2) * sizeof(ulong));
    _divhigh_numerator(t, a, an, bn, n, s);

    if (FLINT_MPN_DIV_SMALL_SHAPE(nn, bn, 0))
    {
        /* (one-limb divisors included: the register division chooses
           between the hardware chain and mpn_divrem_1 by length) */
        qq = t + nn;
        _flint_mpn_tdiv_qr_small(qq, NULL, t, nn, b, bn);
        _divhigh_clamp(qq, qq[n + 1], n);
        _flint_mpn_copy_small(q, qq, n + 1);
    }
    else if (n <= DIVHIGH_BASECASE_MAX)
    {
        ulong dinv = flint_mpn_preinv1(b[bn - 1], b[bn - 2]), top;
        top = _flint_mpn_divapprox_basecase_preinv1(q, t, nn, b, bn, dinv);
        _divhigh_clamp(q, top, n);
    }
    else
    {
        qq = t + nn;
        flint_mpn_divapprox(qq, t, nn, b, bn);
        _divhigh_clamp(qq, qq[n + 1], n);
        flint_mpn_copyi(q, qq, n + 1);
    }

    if (t != stack)
        flint_free(t);
    return s;
}

ulong
flint_mpn_invhigh(nn_ptr q, nn_srcptr b, slong bn, slong n)
{
    ulong half = UWORD(1) << (FLINT_BITS - 1);
    return flint_mpn_divhigh(q, &half, 1, b, bn, n);
}

/*
    The exact floor of T = 2^-s (A / Y) B^(n+1) with s = [A >= Y] from the
    full operands, and whether T is not an integer. The numerator
    A B^k / 2^s (k = bn - an + n + 1) is formed exactly when k >= 0; for
    k < 0 the divisor is extended by -k zero limbs instead, so that the
    quotient of the full A by Y B^-k is floor(T) and its remainder, with
    the bit shifted out of A when s = 1, decides exactness.
*/
int
flint_mpn_divhigh_exact(nn_ptr q, ulong * sp, nn_srcptr a, slong an, nn_srcptr b, slong bn, slong n)
{
    slong k, nn, dn;
    ulong s, odd;
    nn_ptr t, d, qq, r;
    int inexact;
    TMP_INIT;

    FLINT_ASSERT(n >= 1 && an >= 1 && bn >= 1);

    s = _fraction_ge(a, an, b, bn);
    k = bn - an + n + 1;
    nn = an + FLINT_MAX(k, 0);
    dn = bn + FLINT_MAX(-k, 0);
    FLINT_ASSERT(nn - dn + 1 == n + 2);

    TMP_START;
    t = TMP_ALLOC((nn + dn + (n + 2) + dn) * sizeof(ulong));
    d = t + nn;
    qq = d + dn;
    r = qq + (n + 2);

    flint_mpn_zero(t, nn - an);
    odd = 0;
    if (s == 0)
        flint_mpn_copyi(t + nn - an, a, an);
    else
    {
        odd = mpn_rshift(t + nn - an, a, an, 1);
        if (nn > an)
        {
            /* the bit shifted out lands in the padding: nothing lost */
            t[nn - an - 1] = odd;
            odd = 0;
        }
        else
            odd = odd >> (FLINT_BITS - 1);
    }

    if (dn > bn)
    {
        flint_mpn_zero(d, dn - bn);
        flint_mpn_copyi(d + dn - bn, b, bn);
    }
    else
        d = (nn_ptr) b;

    flint_mpn_tdiv_qr(qq, r, t, nn, d, dn);
    FLINT_ASSERT(qq[n + 1] == 0);

    inexact = (odd != 0) || !flint_mpn_zero_p(r, dn);
    flint_mpn_copyi(q, qq, n + 1);
    *sp = s;
    TMP_END;
    return inexact;
}
