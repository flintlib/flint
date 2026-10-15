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
    Exponentials.

    exp(x) for |x| < 2^(FLINT_BITS - 3) (larger arguments overflow or underflow the
    exponent range) is reduced as |x| = q log 2 + t, t in [0, log 2), in
    fixed point with one integral and nk + 1 fraction limbs, nk the
    kernel limbs: q is a lower approximation of |x| / log 2 from two
    limbs of 1/log 2 (exact or one low, then corrected; the reduction
    step is that of mp_real/exp.c), and log 2 is
    the floor L of log 2 B^(nk + 1), so that the reduction errs by
    (q + 1)(log 2 - L) < 2^(FLINT_BITS - 2) B^-(nk + 1), a quarter ulp of B^-nk. For
    x < 0, exp(x) = 2^-(q+1) exp(L - t'), t' = |x| - q log 2. Then
    exp(x) = 2^k exp(t) with the kernel on t at nk fraction limbs.
    Positive x < 1 go to the kernel directly, and small negative x to
    the alternating series of exp(-|x|).
*/

/* exp(x) for the normal nfloat x = (-1)^sgnbit (d, n) 2^(e - FLINT_BITS n),
   e <= NFLOAT_EXP_MAX_ARG_EXP, with nk kernel limbs: sets (y, nk + 1)
   (y needs nk + 2 limbs) and *err, and returns k such that
   exp(x) = Y 2^(k - FLINT_BITS nk) within *err units, Y = (y, nk + 1).
   Y is in [B^nk, 2 B^nk) except on the paths with k = 0 for x < 1. */
slong
_nfloat_exp_mpn(nn_ptr y, ulong * err, nn_srcptr d, slong n, slong e, int sgnbit, slong nk)
{
    ulong X[NFLOAT_ELEM_MAX_LIMBS + 3];
    nn_srcptr L;
    ulong q;
    slong N, k;
    int trunc;

    if (e <= 0)
    {
        slong z = -e;

        if (!sgnbit)
        {
            /* x in (0, 1): truncating x costs one ulp, times exp' < e */
            trunc = _nfloat_get_fixed(X, nk, nk, d, n, e);
            _mp_real_exp_kernel(y, err, X, nk);
            *err += 3 * trunc;
            return 0;
        }

        /* the series of exp(-|x|) where it is faster (exp' <= 1) */
        (void) z;
        trunc = _nfloat_get_fixed(X, nk, nk, d, n, e);
        if (_mp_real_exp_neg_series(y, err, X, nk))
        {
            *err += trunc;
            return 0;
        }
    }

    /* X = |x| with one integral and N fraction limbs (truncated) */
    N = nk + 1;
    _nfloat_get_fixed(X, N + 1, N, d, n, e);
    L = _mp_real_const_ptr(MP_REAL_CONST_ID_LOG2, N);

    /* q <= |x| / log 2, exact or one less, and t = X - q L resp.
       (q + 1) L - |x| from the two's complement (mp_real/impl.h) */
    q = _mp_real_exp_quotient(X[N], X[N - 1]);
    if (sgnbit)
        mpn_neg(X, X, N + 1);
    k = _mp_real_exp_reduce_step(X, L, N, q, sgnbit);

    FLINT_ASSERT(X[N] == 0);

    /* t within two ulps of B^-nk (the reduction, the truncation of x and
       of the low fraction limb), times exp' < 2 */
    _mp_real_exp_kernel(y, err, X + 1, nk);
    *err += 4;
    return k;
}

#if FLINT_BITS == 64
/*
    In-register fast paths at one and two limbs for all |x| below the
    overflow threshold (e <= NFLOAT_EXP_MAX_ARG_EXP): the reduction
    |x| = q log 2 + t (_nfloat_small_exp_reduce_3 resp. _4, with one
    integral and three resp. four fraction limbs and L, the floor of
    log 2, at as many limbs: the reduction errs by (q + 1) 2^-192 resp.
    (q + 1) 2^-256, below 0.36 ulp of B^-2 resp. B^-3 for q < 2^61.6),
    then the table-driven kernel of mp_real at two resp. three limbs on
    the top fraction limbs of t (or of L - t for x < 0, exp(x) =
    2^-(q+1) exp(L - t)): one more ulp for the dropped limb, twice
    through exp' < 2, within NFLOAT_SMALL_EXP_ERR_2 = ERR + 3 ulps in
    all. Results beyond the exponent range over- or underflow in the
    final normalization.

    expm1 subtracts 1 in the same fixed point (2^k Y - 1 = 2^k (Y - 2^-k)
    for x > 0, 1 - 2^-(q+1) Y' for x < 0), for |x| >= 2^-41, where the
    absolute error stays below 2^-81 resp. 2^-145 relative to the result.
*/

static int
_nfloat_exp_fast_1(nfloat_ptr res, ulong d, slong e, int sgnbit, gr_ctx_t ctx)
{
    ulong y2, y1, y0;
    slong k = _nfloat_small_exp_1(&y2, &y1, &y0, d, e, sgnbit);
    return _nfloat_small_set_1(res, y2, y1, y0, k, NFLOAT_SMALL_EXP_ERR_2, 0, ctx);
}

static int
_nfloat_exp_fast_2(nfloat_ptr res, ulong d1, ulong d0, slong e, int sgnbit, gr_ctx_t ctx)
{
    ulong y3, y2, y1, y0;
    slong k = _nfloat_small_exp_2(&y3, &y2, &y1, &y0, d1, d0, e, sgnbit);
    return _nfloat_small_set_2(res, y3, y2, y1, y0, k, NFLOAT_SMALL_EXP_ERR_3, 0, ctx);
}

/* expm1 for -40 <= e <= NFLOAT_EXP_MAX_ARG_EXP */
static int
_nfloat_expm1_fast_1(nfloat_ptr res, ulong d, slong e, int sgnbit, gr_ctx_t ctx)
{
    ulong y2, y1, y0, err = NFLOAT_SMALL_EXP_ERR_2 + 1;
    slong k = _nfloat_small_exp_1(&y2, &y1, &y0, d, e, sgnbit);

    if (!sgnbit)
    {
        /* Y - 2^-k, k >= 0: the bit 2 FLINT_BITS - k of Y (or below an
           ulp) */
        if (k <= 2 * FLINT_BITS)
        {
            slong b = 2 * FLINT_BITS - k;
            ulong b2 = 0, b1 = 0, b0 = 0;
            if (b >= 2 * FLINT_BITS)
                b2 = UWORD(1);
            else if (b >= FLINT_BITS)
                b1 = UWORD(1) << (b - FLINT_BITS);
            else
                b0 = UWORD(1) << b;
            sub_dddmmmsss(y2, y1, y0, y2, y1, y0, b2, b1, b0);
        }
        return _nfloat_small_set_1(res, y2, y1, y0, k, err, 0, ctx);
    }
    else
    {
        /* 1 - Y 2^(k - 2 FLINT_BITS), k <= -1 (one more unit for the
           shifted-out bits) */
        _nfloat_small_rshift_3(&y2, &y1, &y0, y2, y1, y0, -k);
        sub_dddmmmsss(y2, y1, y0, UWORD(1), UWORD(0), UWORD(0), y2, y1, y0);
        return _nfloat_small_set_1(res, y2, y1, y0, 0, err, 1, ctx);
    }
}

static int
_nfloat_expm1_fast_2(nfloat_ptr res, ulong d1, ulong d0, slong e, int sgnbit, gr_ctx_t ctx)
{
    ulong y3, y2, y1, y0, err = NFLOAT_SMALL_EXP_ERR_3 + 1;
    slong k = _nfloat_small_exp_2(&y3, &y2, &y1, &y0, d1, d0, e, sgnbit);

    if (!sgnbit)
    {
        if (k <= 3 * FLINT_BITS)
        {
            slong b = 3 * FLINT_BITS - k;
            ulong b3 = 0, b2 = 0, b1 = 0, b0 = 0;
            if (b >= 3 * FLINT_BITS)
                b3 = UWORD(1);
            else if (b >= 2 * FLINT_BITS)
                b2 = UWORD(1) << (b - 2 * FLINT_BITS);
            else if (b >= FLINT_BITS)
                b1 = UWORD(1) << (b - FLINT_BITS);
            else
                b0 = UWORD(1) << b;
            sub_ddddmmmmssss(y3, y2, y1, y0, y3, y2, y1, y0, b3, b2, b1, b0);
        }
        return _nfloat_small_set_2(res, y3, y2, y1, y0, k, err, 0, ctx);
    }
    else
    {
        _nfloat_small_rshift_4(&y3, &y2, &y1, &y0, y3, y2, y1, y0, -k);
        sub_ddddmmmmssss(y3, y2, y1, y0, UWORD(1), UWORD(0), UWORD(0), UWORD(0), y3, y2, y1, y0);
        return _nfloat_small_set_2(res, y3, y2, y1, y0, 0, err, 1, ctx);
    }
}
#endif

int
nfloat_exp(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx)
{
    ulong y[NFLOAT_ELEM_MAX_LIMBS + 2];
    ulong err;
    slong n, e, nk, k;
    int sgnbit;

    if (NFLOAT_IS_SPECIAL(x))
    {
        if (NFLOAT_IS_ZERO(x))
            return nfloat_one(res, ctx);
        if (NFLOAT_IS_POS_INF(x))
            return nfloat_pos_inf(res, ctx);
        if (NFLOAT_IS_NEG_INF(x))
            return nfloat_zero(res, ctx);
        return nfloat_nan(res, ctx);
    }

    n = NFLOAT_CTX_NLIMBS(ctx);
    e = NFLOAT_EXP(x);
    sgnbit = NFLOAT_SGNBIT(x);

    if (e > NFLOAT_EXP_MAX_ARG_EXP)
        return sgnbit ? _nfloat_underflow(res, 0, ctx) : _nfloat_overflow(res, 0, ctx);

    /* |x| < 2^e: exp(x) = 1 + delta, |delta| < 2^(e + 1) <= 2^(-FLINT_BITS n) */
    if (e < -FLINT_BITS * n
        || (e < -NFLOAT_CTX_FUNC_PREC(ctx) && !NFLOAT_CTX_HAS_DIRECTED_ROUNDING(ctx)))
        return _nfloat_one_plus_tiny(res, 0, sgnbit, -FLINT_BITS * n, ctx);

#if FLINT_BITS == 64
    if (n == 1)
        return _nfloat_exp_fast_1(res, NFLOAT_D(x)[0], e, sgnbit, ctx);
    if (n == 2)
        return _nfloat_exp_fast_2(res, NFLOAT_D(x)[1], NFLOAT_D(x)[0], e, sgnbit, ctx);
#endif

    nk = _nfloat_elem_limbs(ctx, 0);
    k = _nfloat_exp_mpn(y, &err, NFLOAT_D(x), n, e, sgnbit, nk);
    return _nfloat_set_mpn_err(res, y, nk + 1, k - FLINT_BITS * nk, err, err, 0, ctx);
}

int
nfloat_expm1(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx)
{
    ulong y[NFLOAT_ELEM_MAX_LIMBS + 2];
    ulong err;
    slong n, e, nk, k;
    int sgnbit;

    if (NFLOAT_IS_SPECIAL(x))
    {
        if (NFLOAT_IS_ZERO(x))
            return nfloat_zero(res, ctx);
        if (NFLOAT_IS_POS_INF(x))
            return nfloat_pos_inf(res, ctx);
        if (NFLOAT_IS_NEG_INF(x))
            return nfloat_neg_one(res, ctx);
        return nfloat_nan(res, ctx);
    }

    n = NFLOAT_CTX_NLIMBS(ctx);
    e = NFLOAT_EXP(x);
    sgnbit = NFLOAT_SGNBIT(x);

    if (e > NFLOAT_EXP_MAX_ARG_EXP)
    {
        if (!sgnbit)
            return _nfloat_overflow(res, 0, ctx);
        /* -1 + exp(x), exp(x) < 2^(-2^(FLINT_BITS - 3)) */
        return _nfloat_one_plus_tiny(res, 1, 1, -FLINT_BITS * n, ctx);
    }

    /* expm1(x) = x (1 + delta), delta = x/2 + x^2/6 + ..., |delta| < |x| */
    if (e < -FLINT_BITS * n
        || (e < -NFLOAT_CTX_FUNC_PREC(ctx) && !NFLOAT_CTX_HAS_DIRECTED_ROUNDING(ctx)))
        return _nfloat_set_tiny_rel_err(res, x, sgnbit ? -1 : 1, ctx);

#if FLINT_BITS == 64
    if (n == 1 && e >= -40)
        return _nfloat_expm1_fast_1(res, NFLOAT_D(x)[0], e, sgnbit, ctx);
    if (n == 2 && e >= -40)
        return _nfloat_expm1_fast_2(res, NFLOAT_D(x)[1], NFLOAT_D(x)[0], e, sgnbit, ctx);
#endif

    if (e <= 0)
    {
        /* |x| < 2^-z: an absolute accuracy of 2^-(p + z) */
        nk = _nfloat_elem_limbs(ctx, -e);
        k = _nfloat_exp_mpn(y, &err, NFLOAT_D(x), n, e, sgnbit, nk);

        if (!sgnbit)
        {
            /* exp(x) - 1 in [0, e - 1) */
            FLINT_ASSERT(k == 0);
            y[nk] -= 1;
            return _nfloat_set_mpn_err(res, y, nk + 1, -FLINT_BITS * nk, err, err, 0, ctx);
        }
        else
        {
            /* 1 - Y 2^k B^-nk, k in {0, -1, -2}: shift Y down by -k bits
               (one more unit for the dropped bits) */
            if (k != 0)
            {
                int dropped;
                FLINT_ASSERT(k >= -2);
                dropped = (y[0] & ((UWORD(1) << (-k)) - 1)) != 0;
                mpn_rshift(y, y, nk + 1, (unsigned int) (-k));
                err = (err >> (-k)) + 1 + dropped;
            }
            mpn_neg(y, y, nk + 1);
            y[nk] += 1;
            return _nfloat_set_mpn_err(res, y, nk + 1, -FLINT_BITS * nk, err, err, 1, ctx);
        }
    }
    else
    {
        nk = _nfloat_elem_limbs(ctx, 0);
        k = _nfloat_exp_mpn(y, &err, NFLOAT_D(x), n, e, sgnbit, nk);

        if (!sgnbit)
        {
            /* Y 2^(k - FLINT_BITS nk) - 1 with k >= 1: subtract the bit
               FLINT_BITS nk - k of Y, or a fraction of a unit */
            ulong elo = err;
            slong b = FLINT_BITS * nk - k;

            if (b >= 0)
                mpn_sub_1(y + b / FLINT_BITS, y + b / FLINT_BITS,
                    nk + 1 - b / FLINT_BITS, UWORD(1) << (b % FLINT_BITS));
            else
                elo++;

            return _nfloat_set_mpn_err(res, y, nk + 1, k - FLINT_BITS * nk, elo, err, 0, ctx);
        }
        else
        {
            /* 1 - exp(x), exp(x) < 1/2 (k <= -1): Y 2^k at nk fraction
               limbs (one more unit for the dropped bits) */
            slong s = -k;

            FLINT_ASSERT(k <= -1);

            if (s >= FLINT_BITS * (nk + 1))
            {
                _nfloat_zero_limbs(y, nk + 1);
                err = 1;
            }
            else
            {
                slong ls = s / FLINT_BITS, i;
                unsigned int bs = s % FLINT_BITS;
                int dropped = 0;

                for (i = 0; i < ls; i++)
                    dropped |= (y[i] != 0);
                if (bs != 0)
                {
                    dropped |= ((y[ls] << (FLINT_BITS - bs)) != 0);
                    mpn_rshift(y, y + ls, nk + 1 - ls, bs);
                }
                else
                    flint_mpn_copyi(y, y + ls, nk + 1 - ls);
                _nfloat_zero_limbs(y + nk + 1 - ls, ls);
                err = (s >= FLINT_BITS) ? 1 : (err >> s) + 1;
                err += dropped;
            }

            mpn_neg(y, y, nk + 1);
            y[nk] += 1;
            return _nfloat_set_mpn_err(res, y, nk + 1, -FLINT_BITS * nk, err, err, 1, ctx);
        }
    }
}

/*
    exp2(x) = 2^k exp(f log 2) with k = floor(x), f = x - k in [0, 1):
    the reduction is exact apart from the truncation of f (at nk + 1
    fraction limbs) and the product f L, L the floor of log 2 B^(nk + 1).
*/
int
nfloat_exp2(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx)
{
    ulong X[NFLOAT_ELEM_MAX_LIMBS + 3];
    ulong T[NFLOAT_ELEM_MAX_LIMBS + 2];
    ulong y[NFLOAT_ELEM_MAX_LIMBS + 2];
    nn_srcptr L;
    ulong err;
    slong n, e, nk, N, k;
    int sgnbit;

    if (NFLOAT_IS_SPECIAL(x))
        return nfloat_exp(res, x, ctx);

    n = NFLOAT_CTX_NLIMBS(ctx);
    e = NFLOAT_EXP(x);
    sgnbit = NFLOAT_SGNBIT(x);

    if (e > NFLOAT_EXP_MAX_ARG_EXP)
        return sgnbit ? _nfloat_underflow(res, 0, ctx) : _nfloat_overflow(res, 0, ctx);

    /* 2^x = 1 + delta, |delta| < 2^e */
    if (e < -FLINT_BITS * n
        || (e < -NFLOAT_CTX_FUNC_PREC(ctx) && !NFLOAT_CTX_HAS_DIRECTED_ROUNDING(ctx)))
        return _nfloat_one_plus_tiny(res, 0, sgnbit, -FLINT_BITS * n, ctx);

    nk = _nfloat_elem_limbs(ctx, 0);
    N = nk + 1;

    /* X = x with one integral limb in two's complement and N fraction
       limbs, truncated towards -inf */
    {
        int trunc;

        /* (written first only to keep gcc from warning) */
        X[0] = X[N] = 0;
        trunc = _nfloat_get_fixed(X, N + 1, N, NFLOAT_D(x), n, e);

        if (sgnbit)
        {
            mpn_neg(X, X, N + 1);
            if (trunc)
                mpn_sub_1(X, X, N + 1, 1);
        }

        k = (slong) X[N];

        /* an integer: 2^k exactly */
        if (!trunc && flint_mpn_zero_p(X, N))
        {
            nfloat_one(res, ctx);
            NFLOAT_EXP(res) = k + 1;
            NFLOAT_HANDLE_UNDERFLOW_OVERFLOW(res, ctx);
            return GR_SUCCESS;
        }

        /* t = f L in [0, log 2): the product truncated (2 ulps of B^-N),
           f's truncation (1 ulp) times log 2 */
        L = _mp_real_const_ptr(MP_REAL_CONST_ID_LOG2, N);
        flint_mpn_mulhigh_n(T, X, L, N);
    }

    /* t within two ulps of B^-nk, times exp' < 2 */
    _mp_real_exp_kernel(y, &err, T + 1, nk);
    err += 4;
    return _nfloat_set_mpn_err(res, y, nk + 1, k - FLINT_BITS * nk, err, err, 0, ctx);
}
