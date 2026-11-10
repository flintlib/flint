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
    Logarithms: the kernels live in mp_real/log_fp.c and the register fast
    paths in arf/elem.h; this file assembles their results into nfloats.
*/

/* res = log(u) resp. log2(u) for u = (d, dn) 2^(E - FLINT_BITS dn)
   normalized, by the kernel at extra more bits of precision */
static int
_nfloat_log_mpn(nfloat_ptr res, nn_srcptr d, slong dn, slong E, int base2, slong extra, gr_ctx_t ctx)
{
    nn_ptr y;
    ulong err;
    slong len, ue;
    int neg, status;
    TMP_INIT;

    TMP_START;
    y = TMP_ALLOC(_mp_real_log_fp_ylen(_nfloat_elem_bits(ctx, extra), dn) * sizeof(ulong));
    neg = _mp_real_log_fp(y, &len, &ue, &err, d, dn, E, base2, _nfloat_elem_bits(ctx, extra));

    if (len == 0)
        status = nfloat_zero(res, ctx);
    else
        status = _nfloat_set_mpn_err(res, y, len, ue, err, err, neg, ctx);

    TMP_END;
    return status;
}

#if FLINT_BITS == 64
static int
_nfloat_log_fast_1(nfloat_ptr res, ulong d, slong E, gr_ctx_t ctx)
{
    ulong r2 = 0, r1 = 0, r0 = 0;
    int neg = _arf_small_log_2(&r2, &r1, &r0, d, 0, E);
    return _nfloat_small_set_1(res, r2, r1, r0, 0, ARF_SMALL_LOG_ERR_2, neg, ctx);
}

static int
_nfloat_log_fast_2(nfloat_ptr res, ulong d1, ulong d0, slong E, gr_ctx_t ctx)
{
    ulong r3 = 0, r2 = 0, r1 = 0, r0 = 0;
    int neg = _arf_small_log_3(&r3, &r2, &r1, &r0, d1, d0, 0, E);
    return _nfloat_small_set_2(res, r3, r2, r1, r0, 0, ARF_SMALL_LOG_ERR_3, neg, ctx);
}
#endif

static int
_nfloat_log(nfloat_ptr res, nfloat_srcptr x, int base2, gr_ctx_t ctx)
{
    slong n = NFLOAT_CTX_NLIMBS(ctx), E;
    nn_srcptr d = NFLOAT_D(x);

    if (NFLOAT_IS_SPECIAL(x) || NFLOAT_SGNBIT(x))
    {
        if (NFLOAT_IS_ZERO(x))
            return nfloat_neg_inf(res, ctx);
        if (NFLOAT_IS_POS_INF(x))
            return nfloat_pos_inf(res, ctx);
        return nfloat_nan(res, ctx);
    }

    E = NFLOAT_EXP(x);

    if (base2 && d[n - 1] == (UWORD(1) << (FLINT_BITS - 1)) && flint_mpn_zero_p(d, n - 1))
        return nfloat_set_si(res, E - 1, ctx);

#if FLINT_BITS == 64
    if (!base2)
    {
        if (n == 1 && (ARF_SMALL_LOG_OK_1(d[0], E)))
            return _nfloat_log_fast_1(res, d[0], E, ctx);
        else if (n == 2 && (ARF_SMALL_LOG_OK_2(d[1], E)))
            return _nfloat_log_fast_2(res, d[1], d[0], E, ctx);
    }
#endif

    /* near 1: log(1 + D) = D (1 + delta), |delta| < |D| < 2^-z (in the
       default rounding mode, also when D is below the function
       precision) */
    if (!base2 && (E == 0 || E == 1))
    {
        slong z = _arf_log_lzb(d, n, E);

        if (z == WORD_MAX)
            return nfloat_zero(res, ctx);

        if (z >= FLINT_BITS * n || (z > NFLOAT_CTX_FUNC_PREC(ctx)
                && !NFLOAT_CTX_HAS_DIRECTED_ROUNDING(ctx)))
        {
            ulong D[NFLOAT_MAX_ALLOC];
            ulong t[NFLOAT_MAX_LIMBS + 1];
            int status;

            /* D exactly as an nfloat (at most FLINT_BITS n bits) */
            if (E == 1)
            {
                _mp_real_fp_get_fixed(t, n + 1, n, d, n, 1);
                status = nfloat_set_mpn_2exp(D, t, n, 0, 0, ctx);
            }
            else
            {
                mpn_neg(t, d, n);
                status = nfloat_set_mpn_2exp(D, t, n, 0, 1, ctx);
            }

            FLINT_ASSERT(status == GR_SUCCESS);
            (void) status;
            /* log(1 + D) < D for D > 0; |log(1 + D)| > |D| for D < 0 */
            return _nfloat_set_tiny_rel_err(res, D, (E == 0) ? 1 : -1, ctx);
        }
    }

    return _nfloat_log_mpn(res, d, n, E, base2, 0, ctx);
}

int
nfloat_log(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx)
{
    return _nfloat_log(res, x, 0, ctx);
}

int
nfloat_log2(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx)
{
    return _nfloat_log(res, x, 1, ctx);
}

int
nfloat_log1p(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx)
{
    slong n, e;
    int sgnbit;

    if (NFLOAT_IS_SPECIAL(x))
    {
        if (NFLOAT_IS_ZERO(x))
            return nfloat_zero(res, ctx);
        if (NFLOAT_IS_POS_INF(x))
            return nfloat_pos_inf(res, ctx);
        return nfloat_nan(res, ctx);
    }

    n = NFLOAT_CTX_NLIMBS(ctx);
    e = NFLOAT_EXP(x);
    sgnbit = NFLOAT_SGNBIT(x);

    /* log1p(x) = x (1 + delta), |delta| < |x| */
    if (e < -FLINT_BITS * n)
        return _nfloat_set_tiny_rel_err(res, x, sgnbit ? 1 : -1, ctx);

    if (sgnbit && e >= 1)
    {
        /* x <= -1 */
        if (e == 1 && NFLOAT_D(x)[n - 1] == (UWORD(1) << (FLINT_BITS - 1))
            && flint_mpn_zero_p(NFLOAT_D(x), n - 1))
            return nfloat_neg_inf(res, ctx);
        return nfloat_nan(res, ctx);
    }

#if FLINT_BITS == 64
    if (n == 1)
    {
        ulong r2 = 0, r1 = 0, r0 = 0, err = 0;
        int neg;
        if (_arf_small_log1p_1(&r2, &r1, &r0, &err, &neg, NFLOAT_D(x)[0], e, sgnbit) == 0)
            return _nfloat_small_set_1(res, r2, r1, r0, 0, err, neg, ctx);
    }
    else if (n == 2)
    {
        ulong r3 = 0, r2 = 0, r1 = 0, r0 = 0, err = 0;
        int neg;
        if (_arf_small_log1p_2(&r3, &r2, &r1, &r0, &err, &neg, NFLOAT_D(x)[1], NFLOAT_D(x)[0], e, sgnbit) == 0)
            return _nfloat_small_set_2(res, r3, r2, r1, r0, 0, err, neg, ctx);
    }
#endif

    {
        nn_ptr y;
        ulong err;
        slong len, ue;
        int neg, status;
        TMP_INIT;

        TMP_START;
        y = TMP_ALLOC(_mp_real_log1p_fp_ylen(_nfloat_elem_bits(ctx, 0), n, e) * sizeof(ulong));
        neg = _mp_real_log1p_fp(y, &len, &ue, &err, NFLOAT_D(x), n, e, sgnbit, _nfloat_elem_bits(ctx, 0));
        FLINT_ASSERT(len != 0);
        status = _nfloat_set_mpn_err(res, y, len, ue, err, err, neg, ctx);
        TMP_END;
        return status;
    }
}
