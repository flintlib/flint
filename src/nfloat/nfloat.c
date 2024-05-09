/*
    Copyright (C) 2024 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <limits.h>
#include <mpfr.h>
#include "long_extras.h"
#include "fmpz.h"
#include "arf.h"
#include "arb.h"
#include "nfloat.h"
#include "gr_generic.h"
#include "gr_special.h"

int
nfloat_write(gr_stream_t out, nfloat_srcptr x, gr_ctx_t ctx)
{
    gr_ctx_t arf_ctx;
    arf_t t;
    int status;

    gr_ctx_init_real_float_arf(arf_ctx, NFLOAT_CTX_PREC(ctx));
    arf_init(t);
    nfloat_get_arf(t, x, ctx);
    status = gr_write(out, t, arf_ctx);
    arf_clear(t);
    return status;
    gr_ctx_clear(arf_ctx);
}

int
nfloat_randtest(nfloat_ptr res, flint_rand_t state, gr_ctx_t ctx)
{
    arf_t t;
    int status;

    arf_init(t);
    arf_randtest(t, state, NFLOAT_CTX_PREC(ctx), 10);
    status = nfloat_set_arf(res, t, ctx);
    arf_clear(t);
    return status;
}

void
nfloat_swap(nfloat_ptr x, nfloat_ptr y, gr_ctx_t ctx)
{
    slong i, n = NFLOAT_CTX_DATA_NLIMBS(ctx);

    for (i = 0; i < n; i++)
        FLINT_SWAP(mp_limb_t, NFLOAT_DATA(x)[i], NFLOAT_DATA(y)[i]);
}

int
nfloat_set(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx)
{
    slong i, n = NFLOAT_CTX_DATA_NLIMBS(ctx);

    for (i = 0; i < n; i++)
        NFLOAT_DATA(res)[i] = NFLOAT_DATA(x)[i];

    return GR_SUCCESS;
}

truth_t
nfloat_equal(nfloat_srcptr x, nfloat_srcptr y, gr_ctx_t ctx)
{
    slong i, n = NFLOAT_CTX_NLIMBS(ctx);

    if (NFLOAT_IS_SPECIAL(x) || NFLOAT_IS_SPECIAL(y))
    {
        return (NFLOAT_EXP(x) == NFLOAT_EXP(y) &&
                NFLOAT_SGNBIT(x) == NFLOAT_SGNBIT(y)) ? T_TRUE : T_FALSE;
    }

    if (NFLOAT_EXP(x) != NFLOAT_EXP(y))
        return T_FALSE;

    if (NFLOAT_SGNBIT(x) != NFLOAT_SGNBIT(y))
        return T_FALSE;

    for (i = 0; i < n; i++)
        if (NFLOAT_D(x)[i] != NFLOAT_D(y)[i])
            return T_FALSE;

    return T_TRUE;
}

int
nfloat_one(nfloat_ptr res, gr_ctx_t ctx)
{
    slong n = NFLOAT_CTX_NLIMBS(ctx);
    NFLOAT_EXP(res) = 1;
    NFLOAT_SGNBIT(res) = 0;
    flint_mpn_zero(NFLOAT_D(res), n - 1);
    NFLOAT_D(res)[n - 1] = UWORD(1) << (FLINT_BITS - 1);
    return GR_SUCCESS;
}

int
nfloat_neg_one(nfloat_ptr res, gr_ctx_t ctx)
{
    slong n = NFLOAT_CTX_NLIMBS(ctx);
    NFLOAT_EXP(res) = 1;
    NFLOAT_SGNBIT(res) = 1;
    flint_mpn_zero(NFLOAT_D(res), n - 1);
    NFLOAT_D(res)[n - 1] = UWORD(1) << (FLINT_BITS - 1);
    return GR_SUCCESS;
}

truth_t
nfloat_is_one(nfloat_srcptr x, gr_ctx_t ctx)
{
    slong n = NFLOAT_CTX_NLIMBS(ctx);

    if (NFLOAT_IS_NAN(x))
        return T_UNKNOWN;

    return ((NFLOAT_EXP(x) == 1) && (NFLOAT_SGNBIT(x) == 0) &&
           flint_mpn_zero_p(NFLOAT_D(x), n - 1) &&
            (NFLOAT_D(x)[n - 1] == (UWORD(1) << (FLINT_BITS - 1)))) ? T_TRUE : T_FALSE;
}

truth_t
nfloat_is_neg_one(nfloat_srcptr x, gr_ctx_t ctx)
{
    slong n = NFLOAT_CTX_NLIMBS(ctx);

    if (NFLOAT_IS_NAN(x))
        return T_UNKNOWN;

    return (NFLOAT_EXP(x) == 1) && (NFLOAT_SGNBIT(x) == 1) &&
           flint_mpn_zero_p(NFLOAT_D(x), n - 1) &&
            (NFLOAT_D(x)[n - 1] == (UWORD(1) << (FLINT_BITS - 1))) ? T_TRUE : T_FALSE;
}

int
nfloat_pos_inf(nfloat_ptr res, gr_ctx_t ctx)
{
    if (!(NFLOAT_CTX_FLAGS(ctx) & NFLOAT_ALLOW_INF))
        return GR_DOMAIN;

    NFLOAT_EXP(res) = NFLOAT_EXP_POS_INF;
    NFLOAT_SGNBIT(res) = 0;
    return GR_SUCCESS;
}

int
nfloat_neg_inf(nfloat_ptr res, gr_ctx_t ctx)
{
    if (!(NFLOAT_CTX_FLAGS(ctx) & NFLOAT_ALLOW_INF))
        return GR_DOMAIN;

    NFLOAT_EXP(res) = NFLOAT_EXP_NEG_INF;
    NFLOAT_SGNBIT(res) = 0;
    return GR_SUCCESS;
}

int
nfloat_nan(nfloat_ptr res, gr_ctx_t ctx)
{
    if (!(NFLOAT_CTX_FLAGS(ctx) & NFLOAT_ALLOW_NAN))
        return GR_UNABLE;

    NFLOAT_EXP(res) = NFLOAT_EXP_NAN;
    NFLOAT_SGNBIT(res) = 0;
    return GR_SUCCESS;
}

int
_nfloat_underflow(nfloat_ptr res, int sgnbit, gr_ctx_t ctx)
{
    if (!(NFLOAT_CTX_FLAGS(ctx) & NFLOAT_ALLOW_UNDERFLOW))
        return GR_UNABLE;

    return nfloat_zero(res, ctx);
}

int
_nfloat_overflow(nfloat_ptr res, int sgnbit, gr_ctx_t ctx)
{
    if (!(NFLOAT_CTX_FLAGS(ctx) & NFLOAT_ALLOW_INF))
        return GR_UNABLE;

    return sgnbit ? nfloat_neg_inf(res, ctx) : nfloat_pos_inf(res, ctx);
}

int
nfloat_set_ui(nfloat_ptr res, ulong x, gr_ctx_t ctx)
{
    if (x == 0)
    {
        return nfloat_zero(res, ctx);
    }
    else
    {
        slong n = NFLOAT_CTX_NLIMBS(ctx);
        slong norm = flint_clz(x);

        NFLOAT_EXP(res) = FLINT_BITS - norm;
        NFLOAT_SGNBIT(res) = 0;
        flint_mpn_zero(NFLOAT_D(res), n - 1);
        NFLOAT_D(res)[n - 1] = x << norm;
        return GR_SUCCESS;
    }
}

int
nfloat_set_si(nfloat_ptr res, slong x, gr_ctx_t ctx)
{
    if (x == 0)
    {
        return nfloat_zero(res, ctx);
    }
    else
    {
        ulong ux = (x > 0) ? x : -x;
        slong n = NFLOAT_CTX_NLIMBS(ctx);
        slong norm = flint_clz(ux);

        NFLOAT_EXP(res) = FLINT_BITS - norm;
        NFLOAT_SGNBIT(res) = (x < 0);
        flint_mpn_zero(NFLOAT_D(res), n - 1);
        NFLOAT_D(res)[n - 1] = ux << norm;
        return GR_SUCCESS;
    }
}

int
nfloat_set_fmpz(nfloat_ptr res, const fmpz_t x, gr_ctx_t ctx)
{
    if (!COEFF_IS_MPZ(*x))
    {
        return nfloat_set_si(res, *x, ctx);
    }
    else
    {
        mpz_ptr z = COEFF_TO_PTR(*x);
        slong zn = z->_mp_size;

        if (zn > 0)
            return _nfloat_set_mpn_2exp(res, z->_mp_d, zn, zn * FLINT_BITS, 0, ctx);
        else
            return _nfloat_set_mpn_2exp(res, z->_mp_d, -zn, -zn * FLINT_BITS, 1, ctx);
    }
}

/* todo: fast code */
int
nfloat_set_fmpq(nfloat_ptr res, const fmpq_t v, gr_ctx_t ctx)
{
    arf_t t;
    int status;
    arf_init(t);
    arf_set_fmpq(t, v, NFLOAT_CTX_PREC(ctx), ARF_RND_DOWN);
    status = nfloat_set_arf(res, t, ctx);
    arf_clear(t);
    return status;
}

/* todo: fast code */
int
nfloat_set_d(nfloat_ptr res, double x, gr_ctx_t ctx)
{
    arf_t t;
    int status;
    arf_init(t);
    arf_set_d(t, x);
    status = nfloat_set_arf(res, t, ctx);
    arf_clear(t);
    return status;
}

int
nfloat_set_str(nfloat_ptr res, const char * x, gr_ctx_t ctx)
{
    int status;

    arb_t t;
    arb_init(t);

    if (!arb_set_str(t, x, NFLOAT_CTX_PREC(ctx) + 20))
    {
        arf_set_round(arb_midref(t), arb_midref(t), NFLOAT_CTX_PREC(ctx), ARF_RND_NEAR);
        status = nfloat_set_arf(res, arb_midref(t), ctx);
    }
    else
    {
        status = gr_generic_set_str_ring_exponents(res, x, ctx);
    }

    arb_clear(t);
    return status;
}

int
nfloat_set_other(nfloat_ptr res, gr_srcptr x, gr_ctx_t x_ctx, gr_ctx_t ctx)
{
    switch (x_ctx->which_ring)
    {
        case GR_CTX_NFLOAT:
            {
                slong nlimbs, x_nlimbs;

                if (NFLOAT_IS_SPECIAL(x))
                {
                    if (NFLOAT_IS_ZERO(x))
                        return nfloat_zero(res, ctx);
                    if (NFLOAT_IS_POS_INF(x))
                        return nfloat_pos_inf(res, ctx);
                    if (NFLOAT_IS_NEG_INF(x))
                        return nfloat_neg_inf(res, ctx);
                    return nfloat_nan(res, ctx);
                }

                nlimbs = NFLOAT_CTX_NLIMBS(ctx);
                x_nlimbs = NFLOAT_CTX_NLIMBS(x_ctx);

                NFLOAT_EXP(res) = NFLOAT_EXP(x);
                NFLOAT_SGNBIT(res) = NFLOAT_SGNBIT(x);

                if (nlimbs <= x_nlimbs)
                {
                    flint_mpn_copyi(NFLOAT_D(res), NFLOAT_D(x) + x_nlimbs - nlimbs, nlimbs);
                }
                else
                {
                    flint_mpn_zero(NFLOAT_D(res), nlimbs - x_nlimbs);
                    flint_mpn_copyi(NFLOAT_D(res) + nlimbs - x_nlimbs, NFLOAT_D(x), x_nlimbs);
                }

                return GR_SUCCESS;
            }

        case GR_CTX_FMPZ:
            return nfloat_set_fmpz(res, x, ctx);

        case GR_CTX_FMPQ:
            return nfloat_set_fmpq(res, x, ctx);

        case GR_CTX_REAL_FLOAT_ARF:
            return nfloat_set_arf(res, x, ctx);

        case GR_CTX_RR_ARB:
            return nfloat_set_arf(res, arb_midref((arb_srcptr) x), ctx);

        default:
            {
                int status;
                arf_t t;

                gr_ctx_t arf_ctx;
                arf_init(t);

                gr_ctx_init_real_float_arf(arf_ctx, NFLOAT_CTX_PREC(ctx));
                status = gr_set_other(t, x, x_ctx, arf_ctx);
                if (status == GR_SUCCESS)
                    status = nfloat_set_arf(res, t, ctx);

                arf_clear(t);
                gr_ctx_clear(arf_ctx);
                return status;
            }
    }
}

int
nfloat_set_arf(nfloat_ptr res, const arf_t x, gr_ctx_t ctx)
{
    if (arf_is_special(x))
    {
        if (arf_is_zero(x))
            return nfloat_zero(res, ctx);
        else if (arf_is_pos_inf(x))
            return nfloat_pos_inf(res, ctx);
        else if (arf_is_neg_inf(x))
            return nfloat_neg_inf(res, ctx);
        else
            return nfloat_nan(res, ctx);
    }
    else
    {
        slong exp, xn;
        mp_srcptr xp;
        int sgnbit;

        ARF_GET_MPN_READONLY(xp, xn, x);

        exp = ARF_EXP(x);
        sgnbit = ARF_SGNBIT(x);

        if (COEFF_IS_MPZ(exp) || exp < NFLOAT_MIN_EXP || exp > NFLOAT_MAX_EXP)
        {
            if (fmpz_sgn(ARF_EXPREF(x)) < 0)
                return _nfloat_underflow(res, sgnbit, ctx);
            else
                return _nfloat_overflow(res, sgnbit, ctx);
        }

        _nfloat_set_mpn_2exp(res, xp, xn, exp, sgnbit, ctx);
    }

    return GR_SUCCESS;
}

int
nfloat_get_arf(arf_t res, nfloat_srcptr x, gr_ctx_t ctx)
{
    if (NFLOAT_IS_SPECIAL(x))
    {
        if (NFLOAT_IS_ZERO(x))
            arf_zero(res);
        else if (NFLOAT_IS_POS_INF(x))
            arf_pos_inf(res);
        else if (NFLOAT_IS_NEG_INF(x))
            arf_neg_inf(res);
        else
            arf_nan(res);
    }
    else
    {
        arf_set_mpn(res, NFLOAT_D(x), NFLOAT_CTX_NLIMBS(ctx), NFLOAT_SGNBIT(x));
        arf_mul_2exp_si(res, res, NFLOAT_EXP(x) - FLINT_BITS * NFLOAT_CTX_NLIMBS(ctx));
    }

    return GR_SUCCESS;
}

int
nfloat_neg(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx)
{
    if (res != x)
    {
        slong i, n = NFLOAT_CTX_DATA_NLIMBS(ctx);

        for (i = 0; i < n; i++)
            NFLOAT_DATA(res)[i] = NFLOAT_DATA(x)[i];
    }

    if (NFLOAT_IS_SPECIAL(res))
    {
        if (NFLOAT_IS_POS_INF(res))
            NFLOAT_EXP(res) = NFLOAT_EXP_NEG_INF;
        else if (NFLOAT_IS_NEG_INF(res))
            NFLOAT_EXP(res) = NFLOAT_EXP_POS_INF;
    }
    else
    {
        NFLOAT_SGNBIT(res) ^= 1;
    }

    return GR_SUCCESS;
}

int
nfloat_abs(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx)
{
    if (res != x)
    {
        slong i, n = NFLOAT_CTX_DATA_NLIMBS(ctx);

        for (i = 0; i < n; i++)
            NFLOAT_DATA(res)[i] = NFLOAT_DATA(x)[i];
    }

    if (NFLOAT_IS_SPECIAL(res))
    {
        if (NFLOAT_IS_NEG_INF(res))
            NFLOAT_EXP(res) = NFLOAT_EXP_POS_INF;
    }
    else
    {
        NFLOAT_SGNBIT(res) = 0;
    }

    return GR_SUCCESS;
}

int
_nfloat_cmp(nfloat_srcptr x, nfloat_srcptr y, gr_ctx_t ctx)
{
    int xs, ys, mc;
    slong xe, ye;

    if (NFLOAT_IS_SPECIAL(x) || NFLOAT_IS_SPECIAL(y))
    {
        if (NFLOAT_CTX_HAS_INF_NAN(ctx))
        {
            if (NFLOAT_IS_NAN(x) || NFLOAT_IS_NAN(y))
                return 0;
            if (NFLOAT_IS_POS_INF(x)) return NFLOAT_IS_POS_INF(y) ? 0 : 1;
            if (NFLOAT_IS_NEG_INF(x)) return NFLOAT_IS_NEG_INF(y) ? 0 : -1;
            if (NFLOAT_IS_POS_INF(y)) return -1;
            if (NFLOAT_IS_NEG_INF(y)) return 1;
        }

        if (NFLOAT_IS_ZERO(x) && NFLOAT_IS_ZERO(y)) return 0;
        if (NFLOAT_IS_ZERO(x)) return NFLOAT_SGNBIT(y) ? 1 : -1;
        if (NFLOAT_IS_ZERO(y)) return NFLOAT_SGNBIT(x) ? -1 : 1;
    }

    xs = NFLOAT_SGNBIT(x);
    ys = NFLOAT_SGNBIT(y);

    if (xs != ys)
        return xs ? -1 : 1;

    xe = NFLOAT_EXP(x);
    ye = NFLOAT_EXP(y);

    if (xe != ye)
        return ((xe < ye) ^ xs) ? -1 : 1;

    mc = mpn_cmp(NFLOAT_D(x), NFLOAT_D(y), NFLOAT_CTX_NLIMBS(ctx));

    if (mc != 0)
        return ((mc < 0) ^ xs) ? -1 : 1;

    return 0;
}

int
_nfloat_cmpabs(nfloat_srcptr x, nfloat_srcptr y, gr_ctx_t ctx)
{
    int mc;
    slong xe, ye;

    if (NFLOAT_IS_SPECIAL(x) || NFLOAT_IS_SPECIAL(y))
    {
        if (NFLOAT_CTX_HAS_INF_NAN(ctx))
        {
            if (NFLOAT_IS_NAN(x) || NFLOAT_IS_NAN(y))
                return 0;
            if (NFLOAT_IS_INF(x)) return NFLOAT_IS_INF(y) ? 0 : 1;
            if (NFLOAT_IS_INF(y)) return -1;
        }

        if (NFLOAT_IS_ZERO(x) && NFLOAT_IS_ZERO(y)) return 0;
        if (NFLOAT_IS_ZERO(x)) return -1;
        if (NFLOAT_IS_ZERO(y)) return 1;
    }

    xe = NFLOAT_EXP(x);
    ye = NFLOAT_EXP(y);

    if (xe != ye)
        return (xe < ye) ? -1 : 1;

    mc = mpn_cmp(NFLOAT_D(x), NFLOAT_D(y), NFLOAT_CTX_NLIMBS(ctx));

    if (mc != 0)
        return (mc < 0) ? -1 : 1;

    return 0;
}

int
nfloat_cmp(int * res, nfloat_srcptr x, nfloat_srcptr y, gr_ctx_t ctx)
{
    if (NFLOAT_IS_NAN(x) || NFLOAT_IS_NAN(y))
        return GR_UNABLE;

    *res = _nfloat_cmp(x, y, ctx);
    return GR_SUCCESS;
}

int
nfloat_cmpabs(int * res, nfloat_srcptr x, nfloat_srcptr y, gr_ctx_t ctx)
{
    if (NFLOAT_IS_NAN(x) || NFLOAT_IS_NAN(y))
        return GR_UNABLE;

    *res = _nfloat_cmpabs(x, y, ctx);
    return GR_SUCCESS;
}

int
nfloat_sgn(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx)
{
    if (NFLOAT_IS_SPECIAL(x))
    {
        if (NFLOAT_IS_ZERO(x))
            return nfloat_zero(res, ctx);

        if (NFLOAT_IS_POS_INF(x))
            return nfloat_one(res, ctx);

        if (NFLOAT_IS_NEG_INF(x))
            return nfloat_neg_one(res, ctx);

        return nfloat_nan(res, ctx);
    }
    else
    {
        return NFLOAT_SGNBIT(x) ? nfloat_neg_one(res, ctx) : nfloat_one(res, ctx);
    }
}

int
nfloat_im(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx)
{
    return nfloat_zero(res, ctx);
}

static mp_limb_pair_t
_flint_mpn_mulhigh_normalised2(mp_ptr rp, mp_srcptr xp, mp_srcptr yp, mp_size_t n)
{
    mp_limb_pair_t ret;

    FLINT_ASSERT(n >= 1);

    if (rp == xp || rp == yp)
    {
        mp_ptr t;
        TMP_INIT;
        TMP_START;
        t = TMP_ALLOC(sizeof(mp_limb_t) * n);
        ret = _flint_mpn_mulhigh_normalised2(t, xp, yp, n);
        flint_mpn_copyi(rp, t, n);
        TMP_END;
        return ret;
    }

    if (xp == yp)
        ret.m1 = flint_mpn_sqrhigh(rp, xp, n);
    else
        ret.m1 = flint_mpn_mulhigh_n(rp, xp, yp, n);

    if (LIMB_MSB_IS_SET(rp[n - 1]))
    {
        ret.m2 = 0;
    }
    else
    {
        ret.m2 = 1;
        mpn_lshift(rp, rp, n, 1);
        rp[0] |= (ret.m1 >> (FLINT_BITS - 1));
        ret.m1 <<= 1;
    }

    return ret;
}

/* handles aliasing */
FLINT_FORCE_INLINE
mp_limb_pair_t flint_mpn_mulhigh_normalised2(mp_ptr rp, mp_srcptr xp, mp_srcptr yp, mp_size_t n)
{
    FLINT_ASSERT(n >= 1);

    if (FLINT_HAVE_MULHIGH_NORMALISED_FUNC(n))
        return flint_mpn_mulhigh_normalised_func_tab[n](rp, xp, yp);
    else
        return _flint_mpn_mulhigh_normalised2(rp, xp, yp, n);
}

FLINT_FORCE_INLINE int
n_signed_sub(ulong * r, ulong x, ulong y)
{
    if (x >= y)
    {
        *r = x - y;
        return 0;
    }
    else
    {
        *r = y - x;
        return 1;
    }
}

int
_nfloat_add_1(nfloat_ptr res, mp_limb_t x0, slong xexp, int xsgnbit, mp_limb_t y0, slong delta, gr_ctx_t ctx)
{
    mp_limb_t hi, lo;

    NFLOAT_SGNBIT(res) = xsgnbit;

    if (delta < FLINT_BITS)
    {
        add_ssaaaa(hi, lo, 0, x0, 0, y0 >> delta);

        if (hi == 0)
        {
            NFLOAT_D(res)[0] = lo;
            NFLOAT_EXP(res) = xexp;
        }
        else
        {
            NFLOAT_D(res)[0] = (lo >> 1) | (UWORD(1) << (FLINT_BITS - 1));
            NFLOAT_EXP(res) = xexp + 1;
            NFLOAT_HANDLE_OVERFLOW(res, ctx);
        }
    }
    else
    {
        NFLOAT_D(res)[0] = x0;
        NFLOAT_EXP(res) = xexp;
        NFLOAT_SGNBIT(res) = xsgnbit;
    }

    return GR_SUCCESS;
}

int
_nfloat_sub_1(nfloat_ptr res, mp_limb_t x0, slong xexp, int xsgnbit, mp_limb_t y0, slong delta, gr_ctx_t ctx)
{
    mp_limb_t u;
    slong norm;

    if (delta == 0)
    {
        NFLOAT_SGNBIT(res) = n_signed_sub(&u, x0, y0) ^ xsgnbit;

        if (u == 0)
            return nfloat_zero(res, ctx);
    }
    else if (delta < FLINT_BITS)
    {
        u = x0 - (y0 >> delta);
        NFLOAT_SGNBIT(res) = xsgnbit;
    }
    else
    {
        NFLOAT_D(res)[0] = x0;
        NFLOAT_EXP(res) = xexp;
        NFLOAT_SGNBIT(res) = xsgnbit;
        return GR_SUCCESS;
    }

    norm = flint_clz(u);
    NFLOAT_D(res)[0] = u << norm;
    NFLOAT_EXP(res) = xexp - norm;
    NFLOAT_HANDLE_UNDERFLOW(res, ctx);
    return GR_SUCCESS;
}

int
_nfloat_add_n(nfloat_ptr res, mp_srcptr xd, slong xexp, int xsgnbit, mp_srcptr yd, slong delta, slong nlimbs, gr_ctx_t ctx)
{
    slong shift_limbs, shift_bits;
    mp_limb_t cy;
    mp_limb_t t[NFLOAT_MAX_LIMBS];

    NFLOAT_SGNBIT(res) = xsgnbit;

    if (delta == 0)
    {
        cy = mpn_add_n(NFLOAT_D(res), xd, yd, nlimbs);
    }
    else if (delta < FLINT_BITS)
    {
        mpn_rshift(t, yd, nlimbs, delta);
        cy = mpn_add_n(NFLOAT_D(res), xd, t, nlimbs);
    }
    else if (delta < nlimbs * FLINT_BITS)
    {
        shift_limbs = delta / FLINT_BITS;
        shift_bits = delta % FLINT_BITS;
        if (shift_bits == 0)
            flint_mpn_copyi(t, yd + shift_limbs, nlimbs - shift_limbs);
        else
            mpn_rshift(t, yd + shift_limbs, nlimbs - shift_limbs, shift_bits);
        cy = mpn_add(NFLOAT_D(res), xd, nlimbs, t, nlimbs - shift_limbs);
    }
    else
    {
        flint_mpn_copyi(NFLOAT_D(res), xd, nlimbs);
        cy = 0;
    }

    if (cy == 0)
    {
        NFLOAT_EXP(res) = xexp;
    }
    else
    {
        mpn_rshift(NFLOAT_D(res), NFLOAT_D(res), nlimbs, 1);
        NFLOAT_D(res)[nlimbs - 1] |= (UWORD(1) << (FLINT_BITS - 1));
        NFLOAT_EXP(res) = xexp + 1;
        NFLOAT_HANDLE_OVERFLOW(res, ctx);
    }

    return GR_SUCCESS;
}

int
_nfloat_sub_n(nfloat_ptr res, mp_srcptr xd, slong xexp, int xsgnbit, mp_srcptr yd, slong delta, slong nlimbs, gr_ctx_t ctx)
{
    slong shift_limbs, shift_bits, n, norm;
    mp_limb_t t[NFLOAT_MAX_LIMBS];

    if (delta == 0)
    {
        NFLOAT_SGNBIT(res) = flint_mpn_signed_sub_n(NFLOAT_D(res), xd, yd, nlimbs) ^ xsgnbit;
    }
    else
    {
        NFLOAT_SGNBIT(res) = xsgnbit;

        if (delta < FLINT_BITS)
        {
            mpn_rshift(t, yd, nlimbs, delta);
            mpn_sub_n(NFLOAT_D(res), xd, t, nlimbs);
            NFLOAT_SGNBIT(res) = xsgnbit;
        }
        else if (delta < nlimbs * FLINT_BITS)
        {
            shift_limbs = delta / FLINT_BITS;
            shift_bits = delta % FLINT_BITS;
            if (shift_bits == 0)
                flint_mpn_copyi(t, yd + shift_limbs, nlimbs - shift_limbs);
            else
                mpn_rshift(t, yd + shift_limbs, nlimbs - shift_limbs, shift_bits);
            mpn_sub(NFLOAT_D(res), xd, nlimbs, t, nlimbs - shift_limbs);
            NFLOAT_SGNBIT(res) = xsgnbit;
        }
        else
        {
            flint_mpn_copyi(NFLOAT_D(res), xd, nlimbs);
            NFLOAT_EXP(res) = xexp;
            NFLOAT_SGNBIT(res) = xsgnbit;
            return GR_SUCCESS;
        }
    }

    n = nlimbs;
    MPN_NORM(NFLOAT_D(res), n);

    if (n != nlimbs)
    {
        if (n == 0)
            return nfloat_zero(res, ctx);

        norm = flint_clz(NFLOAT_D(res)[n - 1]);
        if (norm)
            mpn_lshift(NFLOAT_D(res) + nlimbs - n, NFLOAT_D(res), n, norm);
        else
            flint_mpn_copyd(NFLOAT_D(res) + nlimbs - n, NFLOAT_D(res), n);

        xexp -= (nlimbs - n) * FLINT_BITS + norm;
    }
    else
    {
        norm = flint_clz(NFLOAT_D(res)[nlimbs - 1]);
        if (norm)
            mpn_lshift(NFLOAT_D(res), NFLOAT_D(res), nlimbs, norm);
        xexp -= norm;
    }

    NFLOAT_EXP(res) = xexp;
    NFLOAT_HANDLE_UNDERFLOW(res, ctx);

    return GR_SUCCESS;
}

int
nfloat_add(nfloat_ptr res, nfloat_srcptr x, nfloat_srcptr y, gr_ctx_t ctx)
{
    slong xexp, yexp, delta, nlimbs;
    int xsgnbit, ysgnbit;

    if (NFLOAT_IS_SPECIAL(x) || NFLOAT_IS_SPECIAL(y))
    {
        if (NFLOAT_IS_ZERO(x))
            return nfloat_set(res, y, ctx);
        if (NFLOAT_IS_ZERO(y))
            return nfloat_set(res, x, ctx);
        if (NFLOAT_IS_POS_INF(x) && NFLOAT_IS_POS_INF(y))
            return nfloat_pos_inf(res, ctx);
        if (NFLOAT_IS_NEG_INF(x) && NFLOAT_IS_NEG_INF(y))
            return nfloat_neg_inf(res, ctx);
        return nfloat_nan(res, ctx);
    }

    nlimbs = NFLOAT_CTX_NLIMBS(ctx);

    xexp = NFLOAT_EXP(x);
    yexp = NFLOAT_EXP(y);

    if (xexp < yexp)
    {
        FLINT_SWAP(nfloat_srcptr, x, y);
        FLINT_SWAP(slong, xexp, yexp);
    }

    xsgnbit = NFLOAT_SGNBIT(x);
    ysgnbit = NFLOAT_SGNBIT(y);

    delta = xexp - yexp;

    if (xsgnbit == ysgnbit)
    {
        if (nlimbs == 1)
            return _nfloat_add_1(res, NFLOAT_D(x)[0], xexp, xsgnbit, NFLOAT_D(y)[0], delta, ctx);
        else
            return _nfloat_add_n(res, NFLOAT_D(x), xexp, xsgnbit, NFLOAT_D(y), delta, nlimbs, ctx);
    }
    else
    {
        if (nlimbs == 1)
            return _nfloat_sub_1(res, NFLOAT_D(x)[0], xexp, xsgnbit, NFLOAT_D(y)[0], delta, ctx);
        else
            return _nfloat_sub_n(res, NFLOAT_D(x), xexp, xsgnbit, NFLOAT_D(y), delta, nlimbs, ctx);
    }
}

int
nfloat_sub(nfloat_ptr res, nfloat_srcptr x, nfloat_srcptr y, gr_ctx_t ctx)
{
    slong xexp, yexp, delta, nlimbs;
    int xsgnbit, ysgnbit;

    if (NFLOAT_IS_SPECIAL(x) || NFLOAT_IS_SPECIAL(y))
    {
        if (NFLOAT_IS_ZERO(x))
            return nfloat_neg(res, y, ctx);
        if (NFLOAT_IS_ZERO(y))
            return nfloat_set(res, x, ctx);
        if (NFLOAT_IS_POS_INF(x) && NFLOAT_IS_NEG_INF(y))
            return nfloat_pos_inf(res, ctx);
        if (NFLOAT_IS_NEG_INF(x) && NFLOAT_IS_POS_INF(y))
            return nfloat_neg_inf(res, ctx);
        return nfloat_nan(res, ctx);
    }

    nlimbs = NFLOAT_CTX_NLIMBS(ctx);

    xexp = NFLOAT_EXP(x);
    yexp = NFLOAT_EXP(y);

    xsgnbit = NFLOAT_SGNBIT(x);
    ysgnbit = NFLOAT_SGNBIT(y) ^ 1;

    if (xexp < yexp)
    {
        FLINT_SWAP(nfloat_srcptr, x, y);
        FLINT_SWAP(slong, xexp, yexp);
        FLINT_SWAP(int, xsgnbit, ysgnbit);
    }

    delta = xexp - yexp;

    if (xsgnbit == ysgnbit)
    {
        if (nlimbs == 1)
            return _nfloat_add_1(res, NFLOAT_D(x)[0], xexp, xsgnbit, NFLOAT_D(y)[0], delta, ctx);
        else
            return _nfloat_add_n(res, NFLOAT_D(x), xexp, xsgnbit, NFLOAT_D(y), delta, nlimbs, ctx);
    }
    else
    {
        if (nlimbs == 1)
            return _nfloat_sub_1(res, NFLOAT_D(x)[0], xexp, xsgnbit, NFLOAT_D(y)[0], delta, ctx);
        else
            return _nfloat_sub_n(res, NFLOAT_D(x), xexp, xsgnbit, NFLOAT_D(y), delta, nlimbs, ctx);
    }
}

int
nfloat_mul(nfloat_ptr res, nfloat_srcptr x, nfloat_srcptr y, gr_ctx_t ctx)
{
    mp_limb_pair_t mul_res;
    slong nlimbs;

    if (NFLOAT_IS_SPECIAL(x) || NFLOAT_IS_SPECIAL(y))
    {
        if (NFLOAT_CTX_FLAGS(ctx) & (NFLOAT_ALLOW_INF | NFLOAT_ALLOW_NAN))
        {
            if (NFLOAT_IS_ZERO(x) && !NFLOAT_IS_SPECIAL(y))
                return nfloat_zero(res, ctx);
            if (NFLOAT_IS_ZERO(y) && !NFLOAT_IS_SPECIAL(x))
                return nfloat_zero(res, ctx);
            return GR_UNABLE;
        }
        else
        {
            return nfloat_zero(res, ctx);
        }
    }

    nlimbs = NFLOAT_CTX_NLIMBS(ctx);

    if (nlimbs == 1)
    {
        mp_limb_t hi, lo;

        umul_ppmm(hi, lo, NFLOAT_D(x)[0], NFLOAT_D(y)[0]);

        if (LIMB_MSB_IS_SET(hi))
        {
            NFLOAT_D(res)[0] = hi;
            NFLOAT_EXP(res) = NFLOAT_EXP(x) + NFLOAT_EXP(y);
        }
        else
        {
            NFLOAT_D(res)[0] = (hi << 1) | (lo >> (FLINT_BITS - 1));
            NFLOAT_EXP(res) = NFLOAT_EXP(x) + NFLOAT_EXP(y) - 1;
        }
    }
    else if (nlimbs == 2)
    {
        mp_limb_t r3, r2, r1, FLINT_SET_BUT_UNUSED(r0);

        FLINT_MPN_MUL_2X2(r3, r2, r1, r0, NFLOAT_D(x)[1], NFLOAT_D(x)[0], NFLOAT_D(y)[1], NFLOAT_D(y)[0]);

        if (LIMB_MSB_IS_SET(r3))
        {
            NFLOAT_D(res)[0] = r2;
            NFLOAT_D(res)[1] = r3;
            NFLOAT_EXP(res) = NFLOAT_EXP(x) + NFLOAT_EXP(y);
        }
        else
        {
            NFLOAT_D(res)[0] = (r2 << 1) | (r1 >> (FLINT_BITS - 1));
            NFLOAT_D(res)[1] = (r3 << 1) | (r2 >> (FLINT_BITS - 1));
            NFLOAT_EXP(res) = NFLOAT_EXP(x) + NFLOAT_EXP(y) - 1;
        }
    }
    else
    {
        mul_res = flint_mpn_mulhigh_normalised2(NFLOAT_D(res), NFLOAT_D(x), NFLOAT_D(y), NFLOAT_CTX_NLIMBS(ctx));
        NFLOAT_EXP(res) = NFLOAT_EXP(x) + NFLOAT_EXP(y) - mul_res.m2;
    }

    NFLOAT_SGNBIT(res) = NFLOAT_SGNBIT(x) ^ NFLOAT_SGNBIT(y);

    NFLOAT_HANDLE_UNDERFLOW_OVERFLOW(res, ctx);
    return GR_SUCCESS;
}

int
nfloat_mul_2exp_si(nfloat_ptr res, nfloat_srcptr x, slong y, gr_ctx_t ctx)
{
    if (NFLOAT_IS_SPECIAL(x))
    {
        return nfloat_set(res, x, ctx);
    }
    else
    {
        /* todo */
        if (y < NFLOAT_MIN_EXP || y > NFLOAT_MAX_EXP)
            return GR_UNABLE;

        nfloat_set(res, x, ctx);
        NFLOAT_EXP(res) += y;
        NFLOAT_HANDLE_UNDERFLOW_OVERFLOW(res, ctx);
        return GR_SUCCESS;
    }
}


int
nfloat_inv(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx)
{
    mpfr_t rf, xf;
    slong prec = NFLOAT_CTX_PREC(ctx);
    slong nlimbs = NFLOAT_CTX_NLIMBS(ctx);

    if (NFLOAT_IS_SPECIAL(x))
    {
        if (NFLOAT_IS_ZERO(x))
            return nfloat_nan(res, ctx);
        else
            return nfloat_nan(res, ctx); /* todo */
    }

    if (NFLOAT_D(x)[nlimbs - 1] == (UWORD(1) << (FLINT_BITS - 1)) &&
        flint_mpn_zero_p(NFLOAT_D(x), nlimbs - 1))
    {
        nfloat_set(res, x, ctx);
        NFLOAT_EXP(res) = 2 - NFLOAT_EXP(res);
        NFLOAT_HANDLE_UNDERFLOW_OVERFLOW(res, ctx);
        return GR_SUCCESS;
    }

    /* todo: make sure aliasing is correct */

    rf->_mpfr_d = NFLOAT_D(res);
    rf->_mpfr_prec = prec;
    rf->_mpfr_sign = 1;
    rf->_mpfr_exp = 0;

    xf->_mpfr_d = NFLOAT_D(x);
    xf->_mpfr_prec = prec;
    xf->_mpfr_sign = 1;
    xf->_mpfr_exp = 0;

    mpfr_ui_div(rf, 1, xf, MPFR_RNDZ);

    NFLOAT_EXP(res) = -NFLOAT_EXP(x) + rf->_mpfr_exp;
    NFLOAT_SGNBIT(res) = NFLOAT_SGNBIT(x);

    NFLOAT_HANDLE_UNDERFLOW_OVERFLOW(res, ctx);
    return GR_SUCCESS;
}

int
nfloat_div(nfloat_ptr res, nfloat_srcptr x, nfloat_srcptr y, gr_ctx_t ctx)
{
    mpfr_t rf, xf, yf;
    slong prec = NFLOAT_CTX_PREC(ctx);

    if (NFLOAT_IS_SPECIAL(x) || NFLOAT_IS_SPECIAL(y))
    {
        if (NFLOAT_IS_ZERO(x) && !NFLOAT_IS_SPECIAL(y))
            return nfloat_zero(res, ctx);
        else
            return nfloat_nan(res, ctx); /* todo */
    }

    /* todo: make sure aliasing is correct */
    rf->_mpfr_d = NFLOAT_D(res);
    rf->_mpfr_prec = prec;
    rf->_mpfr_sign = 1;
    rf->_mpfr_exp = 0;

    xf->_mpfr_d = NFLOAT_D(x);
    xf->_mpfr_prec = prec;
    xf->_mpfr_sign = 1;
    xf->_mpfr_exp = 0;

    yf->_mpfr_d = NFLOAT_D(y);
    yf->_mpfr_prec = prec;
    yf->_mpfr_sign = 1;
    yf->_mpfr_exp = 0;

    mpfr_div(rf, xf, yf, MPFR_RNDZ);

    NFLOAT_EXP(res) = NFLOAT_EXP(x) - NFLOAT_EXP(y) + rf->_mpfr_exp;
    NFLOAT_SGNBIT(res) = NFLOAT_SGNBIT(x) ^ NFLOAT_SGNBIT(y);

    NFLOAT_HANDLE_UNDERFLOW_OVERFLOW(res, ctx);
    return GR_SUCCESS;
}

/* check if mpfr_div_ui argument is ulong */
#if ULONG_MAX == UWORD_MAX

int
nfloat_div_ui(nfloat_ptr res, nfloat_srcptr x, ulong y, gr_ctx_t ctx)
{
    mpfr_t rf, xf;
    slong prec = NFLOAT_CTX_PREC(ctx);

    if (NFLOAT_IS_SPECIAL(x) || y == 0)
    {
        if (NFLOAT_IS_ZERO(x) && y != 0)
            return nfloat_zero(res, ctx);
        else
            return nfloat_nan(res, ctx); /* todo */
    }

    /* todo: make sure aliasing is correct */
    rf->_mpfr_d = NFLOAT_D(res);
    rf->_mpfr_prec = prec;
    rf->_mpfr_sign = 1;
    rf->_mpfr_exp = 0;

    xf->_mpfr_d = NFLOAT_D(x);
    xf->_mpfr_prec = prec;
    xf->_mpfr_sign = 1;
    xf->_mpfr_exp = 0;

    mpfr_div_ui(rf, xf, y, MPFR_RNDZ);

    NFLOAT_EXP(res) = NFLOAT_EXP(x) + rf->_mpfr_exp;
    NFLOAT_SGNBIT(res) = NFLOAT_SGNBIT(x);

    NFLOAT_HANDLE_UNDERFLOW_OVERFLOW(res, ctx);
    return GR_SUCCESS;
}

int
nfloat_div_si(nfloat_ptr res, nfloat_srcptr x, slong y, gr_ctx_t ctx)
{
    mpfr_t rf, xf;
    slong prec = NFLOAT_CTX_PREC(ctx);

    if (NFLOAT_IS_SPECIAL(x) || y == 0)
    {
        if (NFLOAT_IS_ZERO(x) && y != 0)
            return nfloat_zero(res, ctx);
        else
            return nfloat_nan(res, ctx); /* todo */
    }

    /* todo: make sure aliasing is correct */
    rf->_mpfr_d = NFLOAT_D(res);
    rf->_mpfr_prec = prec;
    rf->_mpfr_sign = 1;
    rf->_mpfr_exp = 0;

    xf->_mpfr_d = NFLOAT_D(x);
    xf->_mpfr_prec = prec;
    xf->_mpfr_sign = 1;
    xf->_mpfr_exp = 0;

    mpfr_div_si(rf, xf, y, MPFR_RNDZ);

    NFLOAT_EXP(res) = NFLOAT_EXP(x) + rf->_mpfr_exp;
    NFLOAT_SGNBIT(res) = NFLOAT_SGNBIT(x) ^ (y < 0);

    NFLOAT_HANDLE_UNDERFLOW_OVERFLOW(res, ctx);
    return GR_SUCCESS;
}

#else

int
nfloat_div_ui(nfloat_ptr res, nfloat_srcptr x, ulong y, gr_ctx_t ctx)
{
    arf_t t;
    int status;

    arf_init(t);
    nfloat_get_arf(t, x, ctx);
    arf_div_ui(t, t, y, NFLOAT_CTX_PREC(ctx), ARF_RND_DOWN);
    status = nfloat_set_arf(res, t, ctx);
    arf_clear(t);

    NFLOAT_HANDLE_UNDERFLOW_OVERFLOW(res, ctx);
    return status;
}

int
nfloat_div_si(nfloat_ptr res, nfloat_srcptr x, slong y, gr_ctx_t ctx)
{
    arf_t t;
    int status;

    arf_init(t);
    nfloat_get_arf(t, x, ctx);
    arf_div_si(t, t, y, NFLOAT_CTX_PREC(ctx), ARF_RND_DOWN);
    status = nfloat_set_arf(res, t, ctx);
    arf_clear(t);

    NFLOAT_HANDLE_UNDERFLOW_OVERFLOW(res, ctx);
    return status;
}

#endif

int
nfloat_sqrt(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx)
{
    mpfr_t xf, zf;
    int odd_exp;
    slong nlimbs = NFLOAT_CTX_NLIMBS(ctx);
    slong prec = NFLOAT_CTX_PREC(ctx);

    if (NFLOAT_IS_SPECIAL(x))
    {
        if (NFLOAT_IS_NEG_INF(x))
            return nfloat_nan(res, ctx);
        else
            return nfloat_set(res, x, ctx);
    }

    if (NFLOAT_SGNBIT(x))
        return nfloat_nan(res, ctx);

    odd_exp = NFLOAT_EXP(x) & 1;

    /* Powers of two */
    if (odd_exp &&
        NFLOAT_D(x)[nlimbs - 1] == (UWORD(1) << (FLINT_BITS - 1)) &&
        flint_mpn_zero_p(NFLOAT_D(x), nlimbs - 1))
    {
        nfloat_set(res, x, ctx);
        NFLOAT_EXP(res) = (NFLOAT_EXP(res) + 1) / 2;
        return GR_SUCCESS;
    }

    if (res == x)
    {
        zf->_mpfr_d = NFLOAT_D(x);
        zf->_mpfr_prec = prec;
        zf->_mpfr_sign = 1;
        zf->_mpfr_exp = odd_exp;

        mpfr_sqrt(zf, zf, MPFR_RNDZ);
    }
    else
    {
        xf->_mpfr_d = NFLOAT_D(x);
        xf->_mpfr_prec = prec;
        xf->_mpfr_sign = 1;
        xf->_mpfr_exp = odd_exp;

        zf->_mpfr_d = NFLOAT_D(res);
        zf->_mpfr_prec = prec;
        zf->_mpfr_sign = 1;
        zf->_mpfr_exp = 0;

        mpfr_sqrt(zf, xf, MPFR_RNDZ);
    }

    /* floor division */
    NFLOAT_EXP(res) = (NFLOAT_EXP(x) - (NFLOAT_EXP(x) < 0)) / 2 + zf->_mpfr_exp;

    NFLOAT_SGNBIT(res) = 0;

    return GR_SUCCESS;
}

int
nfloat_rsqrt(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx)
{
    mpfr_t xf, zf;
    int odd_exp;
    slong nlimbs = NFLOAT_CTX_NLIMBS(ctx);
    slong prec = NFLOAT_CTX_PREC(ctx);

    if (NFLOAT_IS_SPECIAL(x))
    {
        if (NFLOAT_IS_ZERO(x))
            return nfloat_pos_inf(res, ctx);
        else if (NFLOAT_IS_POS_INF(x))
            return nfloat_zero(res, ctx);
        else
            return nfloat_nan(res, ctx);
    }

    if (NFLOAT_SGNBIT(x))
        return nfloat_nan(res, ctx);

    odd_exp = NFLOAT_EXP(x) & 1;

    /* Powers of two */
    if (odd_exp &&
        NFLOAT_D(x)[nlimbs - 1] == (UWORD(1) << (FLINT_BITS - 1)) &&
        flint_mpn_zero_p(NFLOAT_D(x), nlimbs - 1))
    {
        nfloat_set(res, x, ctx);
        NFLOAT_EXP(res) = ((-NFLOAT_EXP(res) + 1)) / 2 + 1;
        return GR_SUCCESS;
    }

    if (res == x)
    {
        zf->_mpfr_d = NFLOAT_D(x);
        zf->_mpfr_prec = prec;
        zf->_mpfr_sign = 1;
        zf->_mpfr_exp = odd_exp;

        mpfr_rec_sqrt(zf, zf, MPFR_RNDZ);
    }
    else
    {
        xf->_mpfr_d = NFLOAT_D(x);
        xf->_mpfr_prec = prec;
        xf->_mpfr_sign = 1;
        xf->_mpfr_exp = odd_exp;

        zf->_mpfr_d = NFLOAT_D(res);
        zf->_mpfr_prec = prec;
        zf->_mpfr_sign = 1;
        zf->_mpfr_exp = 0;

        mpfr_rec_sqrt(zf, xf, MPFR_RNDZ);
    }

    /* floor division */
    NFLOAT_EXP(res) = -((NFLOAT_EXP(x) - (NFLOAT_EXP(x) < 0)) / 2) + zf->_mpfr_exp;

    NFLOAT_SGNBIT(res) = 0;

    return GR_SUCCESS;
}

static int
_nfloat_func1_via_arf(gr_method_unary_op func, nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx)
{
    gr_ctx_t arf_ctx;
    arf_t t;
    int status;

    gr_ctx_init_real_float_arf(arf_ctx, NFLOAT_CTX_PREC(ctx));
    arf_init(t);
    nfloat_get_arf(t, x, ctx);
    status = func(t, t, arf_ctx);
    if (status == GR_SUCCESS)
        status = nfloat_set_arf(res, t, ctx);
    arf_clear(t);
    gr_ctx_clear(arf_ctx);

    return status;
}

static int
_nfloat_func2_via_arf(gr_method_binary_op func, nfloat_ptr res, nfloat_srcptr x, nfloat_srcptr y, gr_ctx_t ctx)
{
    gr_ctx_t arf_ctx;
    arf_t t, u;
    int status;

    gr_ctx_init_real_float_arf(arf_ctx, NFLOAT_CTX_PREC(ctx));
    arf_init(t);
    arf_init(u);
    nfloat_get_arf(t, x, ctx);
    nfloat_get_arf(u, y, ctx);
    status = func(t, t, u, arf_ctx);
    if (status == GR_SUCCESS)
        status = nfloat_set_arf(res, t, ctx);
    arf_clear(t);
    arf_clear(u);
    gr_ctx_clear(arf_ctx);

    return status;
}

/* todo: fast code */
int nfloat_floor(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx) { return _nfloat_func1_via_arf((gr_method_unary_op) gr_floor, res, x, ctx); }
int nfloat_ceil(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx) { return _nfloat_func1_via_arf((gr_method_unary_op) gr_ceil, res, x, ctx); }
int nfloat_trunc(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx) { return _nfloat_func1_via_arf((gr_method_unary_op) gr_trunc, res, x, ctx); }
int nfloat_nint(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx) { return _nfloat_func1_via_arf((gr_method_unary_op) gr_nint, res, x, ctx); }

int nfloat_pi(nfloat_ptr res, gr_ctx_t ctx)
{
    slong nlimbs = NFLOAT_CTX_NLIMBS(ctx);

    FLINT_ASSERT(nlimbs <= ARB_PI4_TAB_LIMBS);

    NFLOAT_EXP(res) = 2;
    NFLOAT_SGNBIT(res) = 0;
    flint_mpn_copyi(NFLOAT_D(res), arb_pi4_tab + ARB_PI4_TAB_LIMBS - nlimbs, nlimbs);

    return GR_SUCCESS;
}

int nfloat_pow(nfloat_ptr res, nfloat_srcptr x, nfloat_srcptr y, gr_ctx_t ctx) { return _nfloat_func2_via_arf((gr_method_binary_op) gr_pow, res, x, y, ctx); }
int nfloat_exp(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx) { return _nfloat_func1_via_arf((gr_method_unary_op) gr_exp, res, x, ctx); }
int nfloat_expm1(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx) { return _nfloat_func1_via_arf((gr_method_unary_op) gr_expm1, res, x, ctx); }
int nfloat_log(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx) { return _nfloat_func1_via_arf((gr_method_unary_op) gr_log, res, x, ctx); }
int nfloat_log1p(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx) { return _nfloat_func1_via_arf((gr_method_unary_op) gr_log1p, res, x, ctx); }
int nfloat_sin(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx) { return _nfloat_func1_via_arf((gr_method_unary_op) gr_sin, res, x, ctx); }
int nfloat_cos(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx) { return _nfloat_func1_via_arf((gr_method_unary_op) gr_cos, res, x, ctx); }
int nfloat_tan(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx) { return _nfloat_func1_via_arf((gr_method_unary_op) gr_tan, res, x, ctx); }
int nfloat_sinh(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx) { return _nfloat_func1_via_arf((gr_method_unary_op) gr_sinh, res, x, ctx); }
int nfloat_cosh(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx) { return _nfloat_func1_via_arf((gr_method_unary_op) gr_cosh, res, x, ctx); }
int nfloat_tanh(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx) { return _nfloat_func1_via_arf((gr_method_unary_op) gr_tanh, res, x, ctx); }
int nfloat_atan(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx) { return _nfloat_func1_via_arf((gr_method_unary_op) gr_atan, res, x, ctx); }
int nfloat_gamma(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx) { return _nfloat_func1_via_arf((gr_method_unary_op) gr_gamma, res, x, ctx); }
int nfloat_zeta(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx) { return _nfloat_func1_via_arf((gr_method_unary_op) gr_zeta, res, x, ctx); }
