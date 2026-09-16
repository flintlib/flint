/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <math.h>
#include "mpn_extras.h"

/*
    Upper bound on the number of limbs of x^e (x with xn limbs, top limb
    nonzero; e >= 0). Uses e log2(x) computed in double precision from the
    top two limbs, rounded up, which is essentially tight also for small
    bases (where e bits(x) would overestimate by up to a factor two). The
    bound includes the one limb of slack needed by the intermediate
    products of flint_mpn_pow.
*/
mp_size_t
flint_mpn_pow_bound_limbs(mp_srcptr x, mp_size_t xn, ulong e)
{
    double l, bits;

    FLINT_ASSERT(xn >= 1);
    FLINT_ASSERT(x[xn - 1] != 0);

    if (e == 0 || (xn == 1 && x[0] == 1))
        return 1;

    /* small results: the crude bound e bits(x) is cheap and tight enough */
    if (xn == 1 && e <= 4096)
    {
        ulong b = e * FLINT_BIT_COUNT(x[0]);
        if (b <= 4096)
            return b / FLINT_BITS + 2;
    }

    /* log2(x) <= (xn - 1) 64 + log2(top + (next + 1) 2^-64) + rounding */
    l = (double) x[xn - 1];
    if (xn > 1)
        l += ((double) x[xn - 2] + 1.0) * ldexp(1.0, -FLINT_BITS);
    l = log2(l) * (1.0 + 1e-14) + 1e-9 + (double) (xn - 1) * FLINT_BITS;

    bits = (double) e * l * (1.0 + 1e-14) + 2.0;

    if (bits > 1e18)
        return WORD_MAX / FLINT_BITS;

    /* one extra limb: the multiplications write the full product length
       even when the top limb turns out to be zero */
    return (mp_size_t) (bits / FLINT_BITS) + 2;
}

/* r = a * b with an >= 1, bn >= 1 (an + bn limbs written) */
static void
_mul(mp_ptr r, mp_srcptr a, mp_size_t an, mp_srcptr b, mp_size_t bn)
{
    if (an >= bn)
        flint_mpn_mul(r, a, an, b, bn);
    else
        flint_mpn_mul(r, b, bn, a, an);
}

/*
    res = x^e by left-to-right binary exponentiation; returns the number of
    limbs of the result. res must have room for
    flint_mpn_pow_bound_limbs(x, xn, e) limbs. Requires xn >= 1 and
    x[xn-1] != 0.

    Low zero limbs of x are stripped and accounted for by an offset in the
    output; trailing zero bits are shifted out only when the base is short
    (the shifted copy is then negligible), the whole limbs of the resulting
    shift of x^e also becoming an output offset and the remaining bits
    being restored by one in-place shift of the result. As long as the running power fits in one or two limbs, the
    squarings and multiplications use fixed-size limb arithmetic. After
    that the running power alternates between res and a scratch buffer,
    the parity of the remaining out-of-place operations being chosen so
    that the last one lands in res; multiplications by a single-limb base
    are done in place with mpn_mul_1 and do not flip the buffers.
*/
mp_size_t
flint_mpn_pow(mp_ptr res, mp_srcptr x, mp_size_t xn, ulong e)
{
    mp_size_t k = 0, rn, bound;
    ulong shift = 0;
    unsigned int v;
    slong i, j;
    mp_ptr scratch, dst, other;
    mp_srcptr src = NULL;
    mp_limb_t xs[8];
    slong nflips;
    TMP_INIT;

    FLINT_ASSERT(xn >= 1);
    FLINT_ASSERT(x[xn - 1] != 0);

    if (e == 0)
    {
        res[0] = 1;
        return 1;
    }

    /* low zero limbs: x = B^k x', x^e = B^(k e) x'^e */
    while (x[k] == 0)
        k++;
    x += k;
    xn -= k;
    k *= e;
    flint_mpn_zero(res, k);
    res += k;

    /* trailing zero bits, only for short bases: x = 2^v x', so
       x^e = x'^e 2^(v e); the whole-limb part of that shift is applied
       up front (zeroing the low limbs and moving the output pointer, like
       the zero limbs above) and only the remaining bit shift is done at
       the end */
    v = flint_ctz(x[0]);
    if (v != 0 && xn <= 8)
    {
        mp_size_t sl;

        mpn_rshift(xs, x, xn, v);
        if (xs[xn - 1] == 0)
            xn--;
        x = xs;
        shift = (ulong) v * e;
        sl = shift / FLINT_BITS;
        shift %= FLINT_BITS;
        flint_mpn_zero(res, sl);
        res += sl;
        k += sl;
    }

    if (e == 1 || (xn == 1 && x[0] == 1))
    {
        flint_mpn_copyi(res, x, xn);
        rn = xn;
        goto finish;
    }

    i = FLINT_BIT_COUNT(e) - 2;      /* next exponent bit to process */

    if (xn == 1)
    {
        mp_limb_t x0 = x[0], r0 = x0, r1, r2, hi, lo;

        /* running power in one limb */
        for ( ; i >= 0; i--)
        {
            umul_ppmm(hi, lo, r0, r0);
            r0 = lo;
            if (hi != 0)
            {
                r1 = hi;
                if ((e >> i) & 1)
                {
                    /* (the macro declares its own t1: use other names) */
                    mp_limb_t s1 = r1, s0 = r0;
                    FLINT_MPN_MUL_2X1(r2, r1, r0, s1, s0, x0);
                    if (r2 != 0)
                    {
                        res[0] = r0; res[1] = r1; res[2] = r2;
                        rn = 3;
                        i--;
                        goto generic;
                    }
                }
                i--;
                goto two_limbs;
            }
            if ((e >> i) & 1)
            {
                umul_ppmm(hi, lo, r0, x0);
                r0 = lo;
                if (hi != 0)
                {
                    r1 = hi;
                    i--;
                    goto two_limbs;
                }
            }
        }
        res[0] = r0;
        rn = 1;
        goto finish;

    two_limbs:
        /* (r1, r0) with r1 != 0: if bits remain, the next squaring gives
           at least three limbs and is written straight into the buffer the
           generic pass expects (res if the number of remaining squarings
           is odd, else scratch), so nothing is ever copied */
        if (i < 0)
        {
            res[0] = r0; res[1] = r1;
            rn = 2;
            goto finish;
        }

        TMP_START;
        nflips = i + 1;
        bound = flint_mpn_pow_bound_limbs(x, xn, e);
        scratch = TMP_ALLOC(bound * sizeof(mp_limb_t));
        dst = (nflips & 1) ? res : scratch;
        other = (nflips & 1) ? scratch : res;

        FLINT_MPN_SQR_2X2(dst[3], dst[2], dst[1], dst[0], r1, r0);
        rn = (dst[3] != 0) ? 4 : 3;
        if ((e >> i) & 1)
        {
            mp_limb_t c = mpn_mul_1(dst, dst, rn, x0);
            if (c != 0)
                dst[rn++] = c;
        }
        i--;
        src = dst;
        { mp_ptr t = dst; dst = other; other = t; }
        goto loop;
    }

    /* multi-limb base: the running power is x itself (read in place) */
    rn = xn;
    src = x;

generic:
    /* remaining out-of-place operations: one squaring per bit, plus one
       multiplication per set bit unless the base has a single limb; the
       last one must land in res, so the first goes to res iff their number
       is odd (a single-limb base arrives here with its running power in
       res, from which it is moved only if res is the next destination) */
    nflips = 0;
    for (j = i; j >= 0; j--)
        nflips += 1 + (xn > 1 && ((e >> j) & 1));

    if (nflips == 0)
        goto finish;

    TMP_START;
    bound = flint_mpn_pow_bound_limbs(x, xn, e);
    scratch = TMP_ALLOC(bound * sizeof(mp_limb_t));
    dst = (nflips & 1) ? res : scratch;
    other = (nflips & 1) ? scratch : res;

    if (xn == 1)
    {
        /* three limbs from the fixed-size stage, sitting in res */
        src = res;
        if (dst == res)
        {
            flint_mpn_copyi(scratch, res, rn);
            src = scratch;
        }
    }

loop:
    for ( ; i >= 0; i--)
    {
        flint_mpn_sqr(dst, src, rn);
        rn = 2 * rn - (dst[2 * rn - 1] == 0);
        src = dst;
        { mp_ptr t = dst; dst = other; other = t; }

        if ((e >> i) & 1)
        {
            if (xn == 1)
            {
                mp_limb_t c = mpn_mul_1((mp_ptr) src, src, rn, x[0]);
                if (c != 0)
                    ((mp_ptr) src)[rn++] = c;
            }
            else
            {
                _mul(dst, src, rn, x, xn);
                rn = rn + xn - (dst[rn + xn - 1] == 0);
                src = dst;
                { mp_ptr t = dst; dst = other; other = t; }
            }
        }
    }

    FLINT_ASSERT(src == res);
    TMP_END;

finish:
    if (shift != 0)
    {
        /* the bit part of the shift, in place */
        mp_limb_t c = mpn_lshift(res, res, rn, shift);
        if (c != 0)
            res[rn++] = c;
    }

    return rn + k;
}
