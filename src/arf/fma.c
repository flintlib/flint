/*
    Copyright (C) 2014, 2021, 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arf.h"
#include "mpn_extras.h"

/*
    The fused multiply-add res = z + (-1)^tsgnbit |x y| of finite nonzero
    operands (the core of arf_fma, arf_addmul and arf_submul).

    The exact product followed by one exact addition and rounding is the
    standard method. For long operands at high precision the product is
    truncated instead: flint_mpn_mulhigh gives its top an limbs plus a
    guard limb, a lower bound of the product with the true product less
    than FLINT_MPN_MULHIGH_GUARD_ERR_BOUND(an) units of the guard limb
    above it. The exact sum of z and that bound is formed, and rounded
    when its bits between the rounding point and the error band certify
    the rounding (_arf_rounding_certain: neither all zero nor all ones,
    so that the true sum rounds the same way in every mode and is
    inexact). Otherwise (the rounding point lies in the error band,
    cancellation shortened the sum, or the result may be exact) the exact
    product decides.
*/

/* the truncated product from this many limbs of output, for operands at
   least this long (below, the full product with one addition is as
   fast: measured) */
#ifndef ARF_FMA_HIGH_MIN_LIMBS
#define ARF_FMA_HIGH_MIN_LIMBS 24
#endif

FLINT_FORCE_INLINE int
_arf_fma_exact(arf_ptr res, arf_srcptr x, arf_srcptr y, arf_srcptr z,
    int tsgnbit, slong prec, arf_rnd_t rnd)
{
    slong xn, yn, zn, tn, alloc;
    nn_srcptr xptr, yptr, zptr;
    nn_ptr tptr, tptr2;
    fmpz_t texp;
    slong shift;
    int inexact;
    ARF_MUL_TMP_DECL

    ARF_GET_MPN_READONLY(xptr, xn, x);
    ARF_GET_MPN_READONLY(yptr, yn, y);
    ARF_GET_MPN_READONLY(zptr, zn, z);

    fmpz_init(texp);

    _fmpz_add2_fast(texp, ARF_EXPREF(x), ARF_EXPREF(y), 0);
    shift = _fmpz_sub_small(ARF_EXPREF(z), texp);

    alloc = tn = xn + yn;
    ARF_MUL_TMP_ALLOC(tptr2, alloc)
    tptr = tptr2;

    FLINT_MPN_MUL_WITH_SPECIAL_CASES(tptr, xptr, xn, yptr, yn);

    tn -= (tptr[0] == 0);
    tptr += (tptr[0] == 0);

    if (shift >= 0)
        inexact = _arf_add_mpn(res, zptr, zn, ARF_SGNBIT(z), ARF_EXPREF(z),
            tptr, tn, tsgnbit, shift, prec, rnd);
    else
        inexact = _arf_add_mpn(res, tptr, tn, tsgnbit, texp,
            zptr, zn, ARF_SGNBIT(z), -shift, prec, rnd);

    ARF_MUL_TMP_FREE(tptr2, alloc)
    fmpz_clear(texp);

    return inexact;
}

/* the truncated product; returns -1 (res untouched) when the exact
   product must decide */
static int
_arf_fma_high(arf_ptr res, arf_srcptr x, arf_srcptr y, arf_srcptr z,
    int tsgnbit, slong prec, arf_rnd_t rnd)
{
    slong xn, yn, zn, an, rn, shift, e, tn, alloc;
    nn_srcptr xptr, yptr, zptr, tp;
    nn_ptr t;
    fmpz_t texp;
    arf_t s;
    int inexact;
    ARF_MUL_TMP_DECL

    ARF_GET_MPN_READONLY(xptr, xn, x);
    ARF_GET_MPN_READONLY(yptr, yn, y);
    ARF_GET_MPN_READONLY(zptr, zn, z);

    fmpz_init(texp);
    _fmpz_add2_fast(texp, ARF_EXPREF(x), ARF_EXPREF(y), 0);
    shift = _fmpz_sub_small(ARF_EXPREF(z), texp);
    rn = (prec + FLINT_BITS - 1) / FLINT_BITS;

    /* the limbs of the product above the rounding point, two guard limbs
       and one for the normalisation; when the product lies entirely below
       the guard limbs of z, the addition only needs its sign */
    if (shift >= 0)
    {
        if (shift > (zn + 1) * FLINT_BITS && shift > prec + FLINT_BITS + 1)
        {
            ulong one = 1;
            inexact = _arf_add_mpn(res, zptr, zn, ARF_SGNBIT(z), ARF_EXPREF(z),
                &one, 1, tsgnbit, shift, prec, arf_rnd_relaxed_to_strict(rnd));
            fmpz_clear(texp);
            return inexact;
        }
        an = (prec - shift + FLINT_BITS - 1) / FLINT_BITS + 3;
        if (an < 3)
            an = 3;
    }
    else
        an = rn + 3;

    /* cancellation of comparable operands, or a product short enough to
       be formed in full: the exact path */
    if ((tsgnbit != ARF_SGNBIT(z) && shift >= -2 && shift <= 2) ||
        an + 1 >= xn + yn || _flint_mpn_mulhigh_prefer_full(xn, yn, an))
    {
        fmpz_clear(texp);
        return -1;
    }

    arf_init(s);

    if (shift < 0 && -shift > (an + 2) * FLINT_BITS)
    {
        /* z lies below two limbs appended under the guard limb: the sum is
           the product's limbs with 1 added to, or subtracted from, the
           lowest of them (an error of less than one unit of that limb) */
        alloc = an + 3;
        ARF_MUL_TMP_ALLOC(t, alloc)
        if (xptr == yptr && xn == yn)
            t[2] = flint_mpn_sqrhigh(t + 3, xptr, xn, an);
        else
            t[2] = flint_mpn_mulhigh(t + 3, xptr, xn, yptr, yn, an);
        t[1] = 0;
        t[0] = 0;
        if (tsgnbit == ARF_SGNBIT(z))
            t[0] = 1;
        else
            mpn_sub_1(t, t, an + 3, 1);
        _arf_set_round_mpn(s, &e, t, an + 3, tsgnbit, ARF_PREC_EXACT, ARF_RND_DOWN);
        _fmpz_add_fast(ARF_EXPREF(s), texp, e);
        ARF_MUL_TMP_FREE(t, alloc)
    }
    else
    {
        alloc = an + 1;
        ARF_MUL_TMP_ALLOC(t, alloc)
        if (xptr == yptr && xn == yn)
            t[0] = flint_mpn_sqrhigh(t + 1, xptr, xn, an);
        else
            t[0] = flint_mpn_mulhigh(t + 1, xptr, xn, yptr, yn, an);

        /* the exact sum of z and the lower bound (the guard limb as data) */
        if (shift >= 0)
            _arf_add_mpn(s, zptr, zn, ARF_SGNBIT(z), ARF_EXPREF(z),
                t, an + 1, tsgnbit, shift, ARF_PREC_EXACT, ARF_RND_DOWN);
        else
            _arf_add_mpn(s, t, an + 1, tsgnbit, texp,
                zptr, zn, ARF_SGNBIT(z), -shift, ARF_PREC_EXACT, ARF_RND_DOWN);
        ARF_MUL_TMP_FREE(t, alloc)
    }

    inexact = -1;
    if (!arf_is_zero(s))
    {
        /* the true sum lies within bound units of the guard limb above
           or below s: the position of that limb's unit relative to the
           lowest bit of s (below it when trailing zero limbs of the sum
           were stripped: the band then lies under the lowest bit), plus
           the bits of the bound */
        ARF_GET_MPN_READONLY(tp, tn, s);
        e = (_fmpz_sub_small(texp, ARF_EXPREF(s)) + tn * FLINT_BITS - (an + 1) * FLINT_BITS)
            + FLINT_BIT_COUNT(FLINT_MPN_MULHIGH_GUARD_ERR_BOUND(an));
        if (e < 0)
            e = 0;

        if (_arf_rounding_certain(tp, tn, e, prec))
        {
            arf_set_round(res, s, prec, arf_rnd_relaxed_to_strict(rnd));
            inexact = 1;
        }
    }

    arf_clear(s);
    fmpz_clear(texp);
    return inexact;
}

/* res = z + (-1)^tsgnbit |x y| for finite nonzero operands */
FLINT_FORCE_INLINE int
_arf_fma_finite(arf_ptr res, arf_srcptr x, arf_srcptr y, arf_srcptr z,
    int tsgnbit, slong prec, arf_rnd_t rnd)
{
    if (prec != ARF_PREC_EXACT && prec >= ARF_FMA_HIGH_MIN_LIMBS * FLINT_BITS &&
        ARF_XSIZE(x) >= 2 * ARF_FMA_HIGH_MIN_LIMBS && ARF_XSIZE(y) >= 2 * ARF_FMA_HIGH_MIN_LIMBS)
    {
        int r = _arf_fma_high(res, x, y, z, tsgnbit, prec, rnd);
        if (r >= 0)
            return r;
    }

    if (rnd >= ARF_RND_FAST)
        rnd = arf_rnd_relaxed_to_strict(rnd);

    return _arf_fma_exact(res, x, y, z, tsgnbit, prec, rnd);
}

int
arf_fma(arf_ptr res, arf_srcptr x, arf_srcptr y, arf_srcptr z, slong prec, arf_rnd_t rnd)
{
    if (arf_is_special(x) || arf_is_special(y) || arf_is_special(z))
    {
        if (arf_is_zero(z))
        {
            return arf_mul(res, x, y, prec, rnd);
        }
        else if (arf_is_finite(x) && arf_is_finite(y))
        {
            return arf_set_round(res, z, prec, rnd);
        }
        else
        {
            /* todo: speed up */
            arf_t t;
            int inexact;
            arf_init(t);
            arf_mul(t, x, y, ARF_PREC_EXACT, ARF_RND_DOWN);
            inexact = arf_add(res, z, t, prec, rnd);
            arf_clear(t);
            return inexact;
        }
    }

    return _arf_fma_finite(res, x, y, z, ARF_SGNBIT(x) ^ ARF_SGNBIT(y), prec, rnd);
}

int
arf_addmul(arf_ptr z, arf_srcptr x, arf_srcptr y, slong prec, arf_rnd_t rnd)
{
    if (arf_is_special(x) || arf_is_special(y) || arf_is_special(z))
    {
        if (arf_is_zero(z))
        {
            return arf_mul(z, x, y, prec, rnd);
        }
        else if (arf_is_finite(x) && arf_is_finite(y))
        {
            return arf_set_round(z, z, prec, rnd);
        }
        else
        {
            /* todo: speed up */
            arf_t t;
            int inexact;
            arf_init(t);
            arf_mul(t, x, y, ARF_PREC_EXACT, ARF_RND_DOWN);
            inexact = arf_add(z, z, t, prec, rnd);
            arf_clear(t);
            return inexact;
        }
    }

    return _arf_fma_finite(z, x, y, z, ARF_SGNBIT(x) ^ ARF_SGNBIT(y), prec, rnd);
}

int
arf_submul(arf_ptr z, arf_srcptr x, arf_srcptr y, slong prec, arf_rnd_t rnd)
{
    if (arf_is_special(x) || arf_is_special(y) || arf_is_special(z))
    {
        if (arf_is_zero(z))
        {
            return arf_neg_mul(z, x, y, prec, rnd);
        }
        else if (arf_is_finite(x) && arf_is_finite(y))
        {
            return arf_set_round(z, z, prec, rnd);
        }
        else
        {
            /* todo: speed up */
            arf_t t;
            int inexact;
            arf_init(t);
            arf_mul(t, x, y, ARF_PREC_EXACT, ARF_RND_DOWN);
            inexact = arf_sub(z, z, t, prec, rnd);
            arf_clear(t);
            return inexact;
        }
    }

    return _arf_fma_finite(z, x, y, z, ARF_SGNBIT(x) ^ ARF_SGNBIT(y) ^ 1, prec, rnd);
}
