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
    Powers x^y = exp(y log |x|), with the sign (-1)^y for x < 0 and integer
    y (NaN otherwise). With |t| < 2^tb for t = y log |x| (which overflows
    the result unless |t| < 2^(FLINT_BITS - 3)), the logarithm is computed in fixed point
    to a relative accuracy of about 2^-(p + tb + 12) and multiplied exactly
    by the mantissa of y, so that t is known to an absolute accuracy far
    below 2^-p, and the exponential of t is evaluated with the same
    reduction as exp. Simple
    cases (y = +-1, 2, +-1/2, x = 1, powers of two to integer powers) are
    handled directly.
*/

/* 1 if y is an integer, with *odd set to its parity */
static int
_nfloat_is_int(int * odd, nfloat_srcptr y, slong n)
{
    slong e = NFLOAT_EXP(y), b, i;

    if (e <= 0)
        return 0;

    if (e >= FLINT_BITS * n + 1)
    {
        *odd = 0;
        return 1;
    }

    /* the bit of weight 1 is at position b = FLINT_BITS n - e of the mantissa */
    b = FLINT_BITS * n - e;
    for (i = 0; i < b / FLINT_BITS; i++)
        if (NFLOAT_D(y)[i] != 0)
            return 0;
    if (b % FLINT_BITS != 0 && (NFLOAT_D(y)[b / FLINT_BITS] << (FLINT_BITS - b % FLINT_BITS)) != 0)
        return 0;

    *odd = (NFLOAT_D(y)[b / FLINT_BITS] >> (b % FLINT_BITS)) & 1;
    return 1;
}

/* 1 if the mantissa is a power of two */
static int
_nfloat_is_pow2(nfloat_srcptr x, slong n)
{
    return NFLOAT_D(x)[n - 1] == (UWORD(1) << (FLINT_BITS - 1)) && flint_mpn_zero_p(NFLOAT_D(x), n - 1);
}

#if FLINT_BITS == 64
/*
    In-register fast paths at one and two limbs (64-bit only): log |x| =
    Lambda at three fraction limbs (_nfloat_small_log_3 on the mantissa,
    within NFLOAT_SMALL_LOG_ERR_3 = 26 units of B^-3, also near |x| = 1
    since only the absolute error of t matters; at one limb and |y| <
    2^40 the two-limb version within 18 units of B^-2 suffices, 2^-82
    relative to the result), t = y Lambda exactly
    from the mantissas and truncated to three resp. four fraction limbs,
    then the reduction and kernel of the exp fast paths. The error of t,
    below (26 |y| + 1) 2^-192, adds twice as many units of the result
    (exp(t +- delta) = exp(t)(1 +- 2 delta)): within 28 units of B^-2 for
    |y| < 2^63 at one limb, and below 2^-146 relative for |y| < 2^40 at
    two limbs. Other y, and |t| >= 2^61 (beyond the exponent range), take
    the general code.
*/

#define NFLOAT_POW_FAST_DECLINE (-1)

/* (r[0], ..., r[len - 1]) = (a[0], ..., a[len - 1]) >> s for 1 <= s <
   2 FLINT_BITS, truncated */
#define NFLOAT_POW_RSHIFT(r, a, len, s) \
    do { \
        slong __i; \
        ulong __s = (s); \
        if (__s >= FLINT_BITS) \
        { \
            for (__i = 0; __i < (len) - 1; __i++) \
                (r)[__i] = (a)[__i + 1]; \
            (r)[(len) - 1] = 0; \
            __s -= FLINT_BITS; \
        } \
        else \
        { \
            for (__i = 0; __i < (len); __i++) \
                (r)[__i] = (a)[__i]; \
        } \
        if (__s != 0) \
        { \
            for (__i = 0; __i < (len) - 1; __i++) \
                (r)[__i] = ((r)[__i] >> __s) | ((r)[__i + 1] << (FLINT_BITS - __s)); \
            (r)[(len) - 1] >>= __s; \
        } \
    } while (0)

static int
_nfloat_pow_fast_1(nfloat_ptr res, ulong dx, slong ex, ulong dy, slong ey, int ysgn, int rsgn, gr_ctx_t ctx)
{
    ulong l3, l2, l1, l0, h, l, q, t2, t1, t0, y2, y1, y0, eprop;
    ulong P[5], Z[5];
    slong k;
    int neg;

    if (ey > FLINT_BITS - 1 || ey < -(FLINT_BITS - 1))
        return NFLOAT_POW_FAST_DECLINE;

    if (ey <= 40)
    {
        /* Lambda at two fraction limbs suffices: (l3, l2, l1, 0) */
        neg = _nfloat_small_log_2(&l3, &l2, &l1, dx, 0, ex);
        l0 = 0;
    }
    else
        neg = _nfloat_small_log_3(&l3, &l2, &l1, &l0, dx, 0, 0, ex);

    /* P = Lambda my exactly, the value P / B^4 */
    umul_ppmm(P[1], P[0], l0, dy);
    umul_ppmm(h, l, l1, dy);
    add_ssaaaa(P[2], P[1], h, l, UWORD(0), P[1]);
    umul_ppmm(h, l, l2, dy);
    add_ssaaaa(P[3], P[2], h, l, UWORD(0), P[2]);
    umul_ppmm(h, l, l3, dy);
    add_ssaaaa(P[4], P[3], h, l, UWORD(0), P[3]);

    /* Z = floor(|t| B^3) = P >> (FLINT_BITS - ey) */
    NFLOAT_POW_RSHIFT(Z, P, 5, FLINT_BITS - ey);
    if (Z[4] != 0 || Z[3] >= (UWORD(1) << 61))
        return NFLOAT_POW_FAST_DECLINE;

    q = _nfloat_small_exp_reduce_3(&t2, &t1, &t0, Z[3], Z[2], Z[1], Z[0]);
    if (neg ^ ysgn)
    {
        nn_srcptr L = _mp_real_const_ptr(MP_REAL_CONST_ID_LOG2, 3);
        sub_dddmmmsss(t2, t1, t0, L[2], L[1], L[0], t2, t1, t0);
        k = -(slong) (q + 1);
    }
    else
        k = (slong) q;

    _mp_real_small_exp_2(&y2, &y1, &y0, t2, t1);

    /* 2 (18 2^ey 2^64 + 1) 2^-64 units with Lambda at two limbs,
       2 (26 2^ey + 1) 2^-64 units at three */
    if (ey <= 40)
        eprop = (ey >= 0) ? (UWORD(36) << ey) + 2 : 38;
    else
        eprop = (UWORD(52) >> (FLINT_BITS - ey)) + 2;
    return _nfloat_small_set_1(res, y2, y1, y0, k, NFLOAT_SMALL_EXP_ERR_2 + eprop, rsgn, ctx);
}

static int
_nfloat_pow_fast_2(nfloat_ptr res, ulong dx1, ulong dx0, slong ex, ulong dy1, ulong dy0,
    slong ey, int ysgn, int rsgn, gr_ctx_t ctx)
{
    ulong l3, l2, l1, l0, h, l, q, t3, t2, t1, t0, y3, y2, y1, y0, eprop;
    ulong a4, a3, a2, a1, a0, b4, b3, b2, b1, b0;
    ulong P[6], Z[6];
    slong k;
    int neg;

    if (ey > 40 || ey < -(FLINT_BITS - 1))
        return NFLOAT_POW_FAST_DECLINE;

    neg = _nfloat_small_log_3(&l3, &l2, &l1, &l0, dx1, dx0, 0, ex);

    /* P = Lambda (my1, my0) exactly, the value P / B^5 */
    umul_ppmm(a1, a0, l0, dy0);
    umul_ppmm(h, l, l1, dy0);
    add_ssaaaa(a2, a1, h, l, UWORD(0), a1);
    umul_ppmm(h, l, l2, dy0);
    add_ssaaaa(a3, a2, h, l, UWORD(0), a2);
    umul_ppmm(h, l, l3, dy0);
    add_ssaaaa(a4, a3, h, l, UWORD(0), a3);

    umul_ppmm(b1, b0, l0, dy1);
    umul_ppmm(h, l, l1, dy1);
    add_ssaaaa(b2, b1, h, l, UWORD(0), b1);
    umul_ppmm(h, l, l2, dy1);
    add_ssaaaa(b3, b2, h, l, UWORD(0), b2);
    umul_ppmm(h, l, l3, dy1);
    add_ssaaaa(b4, b3, h, l, UWORD(0), b3);

    P[0] = a0;
    add_sssssaaaaaaaaaa(P[5], P[4], P[3], P[2], P[1], UWORD(0), a4, a3, a2, a1,
        b4, b3, b2, b1, b0);

    /* Z = floor(|t| B^4) = P >> (FLINT_BITS - ey) */
    NFLOAT_POW_RSHIFT(Z, P, 6, FLINT_BITS - ey);
    if (Z[5] != 0 || Z[4] >= (UWORD(1) << 61))
        return NFLOAT_POW_FAST_DECLINE;

    q = _nfloat_small_exp_reduce_4(&t3, &t2, &t1, &t0, Z[4], Z[3], Z[2], Z[1], Z[0]);
    if (neg ^ ysgn)
    {
        nn_srcptr L = _mp_real_const_ptr(MP_REAL_CONST_ID_LOG2, 4);
        sub_ddddmmmmssss(t3, t2, t1, t0, L[3], L[2], L[1], L[0], t3, t2, t1, t0);
        k = -(slong) (q + 1);
    }
    else
        k = (slong) q;

    _mp_real_small_exp_3(&y3, &y2, &y1, &y0, t3, t2, t1);

    /* 2 (26 2^ey + 1) units of B^-3 */
    eprop = (ey >= 0) ? (UWORD(52) << ey) + 2 : 54;
    return _nfloat_small_set_2(res, y3, y2, y1, y0, k, NFLOAT_SMALL_EXP_ERR_3 + eprop, rsgn, ctx);
}
#endif

int
nfloat_pow(nfloat_ptr res, nfloat_srcptr x, nfloat_srcptr y, gr_ctx_t ctx)
{
    _nfloat_fix_t L;
    ulong out[NFLOAT_ELEM_MAX_LIMBS + 2];
    slong n, w, nk, k, tb;
    ulong err;
    int odd = 0, isint, rsgn, tsgn;

    if (NFLOAT_IS_SPECIAL(x) || NFLOAT_IS_SPECIAL(y))
    {
        if (NFLOAT_IS_ZERO(y))
            return nfloat_one(res, ctx);
        if (NFLOAT_IS_ZERO(x) && !NFLOAT_IS_SPECIAL(y) && !NFLOAT_SGNBIT(y))
            return nfloat_zero(res, ctx);
        return nfloat_nan(res, ctx);
    }

    n = NFLOAT_CTX_NLIMBS(ctx);
    isint = _nfloat_is_int(&odd, y, n);

    if (NFLOAT_SGNBIT(x) && !isint)
        return nfloat_nan(res, ctx);

    rsgn = NFLOAT_SGNBIT(x) && odd;

    /* |x| = 1 */
    if (NFLOAT_EXP(x) == 1 && _nfloat_is_pow2(x, n))
    {
        nfloat_one(res, ctx);
        NFLOAT_SGNBIT(res) = rsgn;
        return GR_SUCCESS;
    }

    /* y = +-1, 2, +-1/2 */
    if (_nfloat_is_pow2(y, n) && NFLOAT_EXP(y) >= 0 && NFLOAT_EXP(y) <= 2)
    {
        if (NFLOAT_EXP(y) == 1)
            return NFLOAT_SGNBIT(y) ? nfloat_inv(res, x, ctx) : nfloat_set(res, x, ctx);
        if (NFLOAT_EXP(y) == 2 && !NFLOAT_SGNBIT(y))
            return nfloat_sqr(res, x, ctx);
        if (NFLOAT_EXP(y) == 0)
            return NFLOAT_SGNBIT(y) ? nfloat_rsqrt(res, x, ctx) : nfloat_sqrt(res, x, ctx);
    }

    /* (2^e)^y = 2^(e y) for integer y, when it fits: |e|, |y| <
       2^(FLINT_BITS/2 - 2), so that e y + 1 fits in a word */
    if (isint && _nfloat_is_pow2(x, n) && NFLOAT_EXP(y) <= FLINT_BITS / 2 - 2)
    {
        slong e = NFLOAT_EXP(x) - 1;
        slong yi = (slong) (NFLOAT_D(y)[n - 1] >> (FLINT_BITS - NFLOAT_EXP(y)));

        if (NFLOAT_SGNBIT(y))
            yi = -yi;

        if (FLINT_ABS(e) < (WORD(1) << (FLINT_BITS / 2 - 2)))
        {
            nfloat_one(res, ctx);
            NFLOAT_SGNBIT(res) = rsgn;
            NFLOAT_EXP(res) = e * yi + 1;
            NFLOAT_HANDLE_UNDERFLOW_OVERFLOW(res, ctx);
            return GR_SUCCESS;
        }
    }

#if FLINT_BITS == 64
    if (n <= 2)
    {
        int status = (n == 1)
            ? _nfloat_pow_fast_1(res, NFLOAT_D(x)[0], NFLOAT_EXP(x), NFLOAT_D(y)[0],
                NFLOAT_EXP(y), NFLOAT_SGNBIT(y), rsgn, ctx)
            : _nfloat_pow_fast_2(res, NFLOAT_D(x)[1], NFLOAT_D(x)[0], NFLOAT_EXP(x),
                NFLOAT_D(y)[1], NFLOAT_D(y)[0], NFLOAT_EXP(y), NFLOAT_SGNBIT(y), rsgn, ctx);
        if (status != NFLOAT_POW_FAST_DECLINE)
            return status;
    }
#endif

    /* |t| < 2^tb with tb = e_y + bits(|E| + 1) (|log |x|| < |E| + 1); log |x|
       to tb + 4 more bits of relative accuracy */
    {
        slong E = NFLOAT_EXP(x);
        slong ey = NFLOAT_EXP(y);
        ulong aE = (E < 0) ? -(ulong) E : (ulong) E;

        tb = FLINT_MAX(ey, -FLINT_BITS * n - 8) + FLINT_BIT_COUNT(aE + 1);
        if (tb > NFLOAT_EXP_MAX_ARG_EXP + 2)
            tb = NFLOAT_EXP_MAX_ARG_EXP + 2;
        tb = FLINT_MAX(tb, 0);
    }

    /* t at w >= nk limbs, tb + 12 more bits than the result: the
       propagated error below takes a right shift */
    nk = _nfloat_elem_limbs(ctx, 0);
    w = nk + (tb + 12 + FLINT_BITS - 1) / FLINT_BITS;

    /* log |x| = (L->d, L->len) 2^(L->e) within L->err units, to tb + 12
       more bits of relative accuracy than the function precision: an
       absolute accuracy for t = y log |x| of about 2^-(p + 12) */
    _nfloat_log_fix(L, NFLOAT_D(x), n, NFLOAT_EXP(x), tb + 12, ctx);
    if (L->len == 0)
        return GR_UNABLE;

    {
        ulong P[2 * NFLOAT_ELEM_MAX_LIMBS + 2];
        ulong Tm[NFLOAT_ELEM_MAX_LIMBS + 2];
        slong pl, te, texp, sh, dd;
        ulong et, eprop;

        /* t = y L exactly from the mantissas: P 2^(L->e + e_y - FLINT_BITS n),
           within L->err |Y| 2^(L->e + e_y - FLINT_BITS n) < L->err 2^(L->e + e_y) */
        if (L->len >= n)
            flint_mpn_mul(P, L->d, L->len, NFLOAT_D(y), n);
        else
            flint_mpn_mul(P, NFLOAT_D(y), n, L->d, L->len);
        pl = L->len + n;
        while (P[pl - 1] == 0)
            pl--;

        /* the top w limbs, normalized: t = Tm 2^te within the units for
           the dropped bits plus the error of L, L->err 2^(L->e + e_y - te) */
        et = _nfloat_elem_normalize(Tm, &te, P, pl,
            L->e + NFLOAT_EXP(y) - FLINT_BITS * n, 0, w);

        texp = te + FLINT_BITS * w;
        tsgn = L->sgnbit ^ NFLOAT_SGNBIT(y);

        /* |t| >= 2^(FLINT_BITS - 3) overflows or underflows */
        if (texp > NFLOAT_EXP_MAX_ARG_EXP)
            return tsgn ? _nfloat_underflow(res, rsgn, ctx) : _nfloat_overflow(res, rsgn, ctx);

        /* tiny t: 1 + delta, |delta| < 2|t| */
        if (texp < -FLINT_BITS * n - 4)
            return _nfloat_one_plus_tiny(res, rsgn, tsgn, -FLINT_BITS * n, ctx);

        /* exp(t +- delta) = exp(t)(1 +- 2 delta) and exp(t) < 2^(FLINT_BITS nk + 1)
           units: delta 2^(FLINT_BITS nk + 2) units, for the truncation
           delta = et 2^te and the error of L, L->err 2^(L->e + e_y) */
        k = _nfloat_exp_mpn(out, &err, Tm, w, texp, tsgn, nk);

        sh = te + FLINT_BITS * nk + 2;
        dd = L->e + NFLOAT_EXP(y) + FLINT_BITS * nk + 2;
        eprop = 0;

        if (sh <= -FLINT_BITS)
            eprop += 1;
        else if (sh <= 0)
            eprop += (et >> (-sh)) + 1;
        else
            return GR_UNABLE;

        if (dd <= -FLINT_BITS)
            eprop += 1;
        else if (dd <= 0)
            eprop += (L->err >> (-dd)) + 1;
        else if (dd <= FLINT_BITS - 4 && L->err < (UWORD(1) << (FLINT_BITS - 4 - dd)))
            eprop += L->err << dd;     /* (below 2^(FLINT_BITS - 4): kernel units are finer than p) */
        else
            return GR_UNABLE;

        err += eprop;
    }

    return _nfloat_set_mpn_err(res, out, nk + 1, k - FLINT_BITS * nk, err, err, rsgn, ctx);
}
