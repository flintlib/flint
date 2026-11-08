/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "elem.h"

/*
    Arctangents: the kernels _mp_real_atan_fp and _mp_real_atan2_fp
    (mp_real.h) on the mantissa, after the tiny arguments and the
    in-register fast paths of arf/elem.h at one and two limbs.
*/

/* +-c pi/4 for c in {1, 2, 4} at the working precision (c ulps below) */
static int
_nfloat_pi4_mul(nfloat_ptr res, int c, int sgnbit, gr_ctx_t ctx)
{
    ulong t[NFLOAT_ELEM_MAX_LIMBS + 2];
    slong w = _nfloat_elem_limbs(ctx, 0);

    t[w] = mpn_mul_1(t, _mp_real_const_ptr(MP_REAL_CONST_ID_PI4, w), w, c);
    return _nfloat_set_mpn_err(res, t, w + 1, -FLINT_BITS * w, c, c, sgnbit, ctx);
}

int
nfloat_atan(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx)
{
    ulong y[NFLOAT_ELEM_MAX_LIMBS + 2];
    slong n, e, w;
    ulong err;
    int sgnbit;

    if (NFLOAT_IS_SPECIAL(x))
    {
        if (NFLOAT_IS_ZERO(x))
            return nfloat_zero(res, ctx);
        if (NFLOAT_IS_NAN(x))
            return nfloat_nan(res, ctx);
        return _nfloat_pi4_mul(res, 2, NFLOAT_IS_NEG_INF(x), ctx);
    }

    n = NFLOAT_CTX_NLIMBS(ctx);
    e = NFLOAT_EXP(x);
    sgnbit = NFLOAT_SGNBIT(x);

    /* atan(x) = x (1 - delta), 0 < delta < x^2 / 3 */
    if (2 * e <= -FLINT_BITS * n)
        return _nfloat_set_tiny_rel_err(res, x, -1, ctx);

#if FLINT_BITS == 64
    if (n == 1)
    {
        ulong t2 = 0, t1 = 0, t0 = 0;
        if (_arf_small_atan_1(&t2, &t1, &t0, &err, NFLOAT_D(x)[0], e) == 0)
            return _nfloat_small_set_1(res, t2, t1, t0, 0, err, sgnbit, ctx);
    }
    else if (n == 2)
    {
        ulong t3 = 0, t2 = 0, t1 = 0, t0 = 0;
        if (_arf_small_atan_2(&t3, &t2, &t1, &t0, &err, NFLOAT_D(x)[1], NFLOAT_D(x)[0], e) == 0)
            return _nfloat_small_set_2(res, t3, t2, t1, t0, 0, err, sgnbit, ctx);
    }
#endif

    _mp_real_atan_fp(y, &w, &err, NFLOAT_D(x), n, e, _nfloat_elem_bits(ctx, 0), NFLOAT_ELEM_MAX_LIMBS);
    return _nfloat_set_mpn_err(res, y, w + 1, -FLINT_BITS * w, err, err, sgnbit, ctx);
}

/* the magnitude comparison |a| < |b| for normal nfloats */
static int
_nfloat_cmpabs_normal(nfloat_srcptr a, nfloat_srcptr b, slong n)
{
    if (NFLOAT_EXP(a) != NFLOAT_EXP(b))
        return (NFLOAT_EXP(a) < NFLOAT_EXP(b)) ? -1 : 1;
    return mpn_cmp(NFLOAT_D(a), NFLOAT_D(b), n);
}

int
nfloat_atan2(nfloat_ptr res, nfloat_srcptr y, nfloat_srcptr x, gr_ctx_t ctx)
{
    ulong t[NFLOAT_ELEM_MAX_LIMBS + 3];
    nfloat_srcptr num, den;
    slong n, d, len, ue;
    ulong elo, ehi;
    int c, sigma, ysgn, cmp;

    if (NFLOAT_IS_SPECIAL(x) || NFLOAT_IS_SPECIAL(y))
    {
        if (NFLOAT_IS_NAN(x) || NFLOAT_IS_NAN(y))
            return nfloat_nan(res, ctx);

        if (NFLOAT_IS_ZERO(y))
        {
            /* 0 for x >= 0, pi for x < 0 */
            if (NFLOAT_IS_ZERO(x) || NFLOAT_IS_POS_INF(x) ||
                (!NFLOAT_IS_SPECIAL(x) && !NFLOAT_SGNBIT(x)))
                return nfloat_zero(res, ctx);
            return _nfloat_pi4_mul(res, 4, 0, ctx);
        }

        if (NFLOAT_IS_ZERO(x) && !NFLOAT_IS_SPECIAL(y))
            return _nfloat_pi4_mul(res, 2, NFLOAT_SGNBIT(y), ctx);

        /* infinities */
        return nfloat_nan(res, ctx);
    }

    n = NFLOAT_CTX_NLIMBS(ctx);
    ysgn = NFLOAT_SGNBIT(y);
    cmp = _nfloat_cmpabs_normal(y, x, n);

    /* |y| = |x|: pi/4 or 3 pi/4 */
    if (cmp == 0)
        return _nfloat_pi4_mul(res, NFLOAT_SGNBIT(x) ? 3 : 1, ysgn, ctx);

    /* theta = c pi/4 + sigma atan(t), t = num / den < 1 */
    if (cmp < 0)
    {
        num = y;
        den = x;
        c = NFLOAT_SGNBIT(x) ? 4 : 0;
        sigma = NFLOAT_SGNBIT(x) ? -1 : 1;
    }
    else
    {
        num = x;
        den = y;
        c = 2;
        sigma = NFLOAT_SGNBIT(x) ? 1 : -1;
    }

    /* t = (Mnum / Mden) 2^-d, d >= 0 */
    d = NFLOAT_EXP(den) - NFLOAT_EXP(num);

#if FLINT_BITS == 64
    if (n == 1)
    {
        if (_arf_small_atan2_1(t + 2, t + 1, t, &elo, NFLOAT_D(num)[0], NFLOAT_D(den)[0], d, c, sigma) == 0)
            return _nfloat_small_set_1(res, t[2], t[1], t[0], 0, elo, ysgn, ctx);
    }
    else if (n == 2)
    {
        if (_arf_small_atan2_2(t + 3, t + 2, t + 1, t, &elo, NFLOAT_D(num)[1], NFLOAT_D(num)[0],
                NFLOAT_D(den)[1], NFLOAT_D(den)[0], d, c, sigma) == 0)
            return _nfloat_small_set_2(res, t[3], t[2], t[1], t[0], 0, elo, ysgn, ctx);
    }
#endif

    _mp_real_atan2_fp(t, &len, &ue, &elo, &ehi, NFLOAT_D(num), n, NFLOAT_D(den), n, d, c, sigma, _nfloat_elem_bits(ctx, 0), NFLOAT_ELEM_MAX_LIMBS + 1);
    return _nfloat_set_mpn_err(res, t, len, ue, elo, ehi, ysgn, ctx);
}
