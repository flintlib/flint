/*
    Copyright (C) 2024 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "mpn_extras.h"
#include "gr.h"
#include "gr_mat.h"
#include "gr_generic.h"
#include "acf.h"
#include "acb.h"
#include "mag.h"
#include "nfloat.h"
#include "impl.h"
#include "elem.h"

static int
_flint_mpn_signed_add_n(nn_ptr res, nn_srcptr x, int xsgnbit, nn_srcptr y, int ysgnbit, mp_size_t n)
{
    if (xsgnbit == ysgnbit)
        mpn_add_n(res, x, y, n);
    else
    {
        if (mpn_cmp(x, y, n) >= 0)
            mpn_sub_n(res, x, y, n);
        else
        {
            mpn_sub_n(res, y, x, n);
            xsgnbit = !xsgnbit;
        }
    }

    return xsgnbit;
}

int
nfloat_complex_get_acf(acf_t res, nfloat_complex_srcptr x, gr_ctx_t ctx)
{
    int status;
    status = nfloat_get_arf(acf_realref(res), NFLOAT_COMPLEX_RE(x, ctx), ctx);
    status |= nfloat_get_arf(acf_imagref(res), NFLOAT_COMPLEX_IM(x, ctx), ctx);
    return status;
}

int
nfloat_complex_set_acf(nfloat_complex_ptr res, const acf_t x, gr_ctx_t ctx)
{
    int status;
    status = nfloat_set_arf(NFLOAT_COMPLEX_RE(res, ctx), acf_realref(x), ctx);
    status |= nfloat_set_arf(NFLOAT_COMPLEX_IM(res, ctx), acf_imagref(x), ctx);
    return status;
}

int
nfloat_complex_get_acb(acb_t res, nfloat_complex_srcptr x, gr_ctx_t ctx)
{
    int status;
    status = nfloat_get_arf(arb_midref(acb_realref(res)), NFLOAT_COMPLEX_RE(x, ctx), ctx);
    mag_zero(arb_radref(acb_realref(res)));
    status |= nfloat_get_arf(arb_midref(acb_imagref(res)), NFLOAT_COMPLEX_IM(x, ctx), ctx);
    mag_zero(arb_radref(acb_imagref(res)));
    return status;
}

int
nfloat_complex_set_acb(nfloat_complex_ptr res, const acb_t x, gr_ctx_t ctx)
{
    int status;
    status = nfloat_set_arf(NFLOAT_COMPLEX_RE(res, ctx), arb_midref(acb_realref(x)), ctx);
    status |= nfloat_set_arf(NFLOAT_COMPLEX_IM(res, ctx), arb_midref(acb_imagref(x)), ctx);
    return status;
}

int
nfloat_complex_write(gr_stream_t out, nfloat_complex_srcptr x, gr_ctx_t ctx)
{
    gr_ctx_t acf_ctx;
    acf_t t;
    int status;

    gr_ctx_init_complex_float_acf(acf_ctx, NFLOAT_CTX_PREC(ctx));
    acf_init(t);
    nfloat_get_arf(acf_realref(t), NFLOAT_COMPLEX_RE(x, ctx), ctx);
    nfloat_get_arf(acf_imagref(t), NFLOAT_COMPLEX_IM(x, ctx), ctx);
    status = gr_write(out, t, acf_ctx);
    acf_clear(t);
    return status;
    gr_ctx_clear(acf_ctx);
}

int
nfloat_complex_randtest(nfloat_complex_ptr res, flint_rand_t state, gr_ctx_t ctx)
{
    int status;
    status = nfloat_randtest(NFLOAT_COMPLEX_RE(res, ctx), state, ctx);
    status |= nfloat_randtest(NFLOAT_COMPLEX_IM(res, ctx), state, ctx);
    return status;
}

void
nfloat_complex_swap(nfloat_complex_ptr x, nfloat_complex_ptr y, gr_ctx_t ctx)
{
    slong i, n = NFLOAT_COMPLEX_CTX_DATA_NLIMBS(ctx);

    for (i = 0; i < n; i++)
        FLINT_SWAP(ulong, NFLOAT_DATA(x)[i], NFLOAT_DATA(y)[i]);
}

int
nfloat_complex_set(nfloat_complex_ptr res, nfloat_complex_ptr x, gr_ctx_t ctx)
{
    slong i, n = NFLOAT_COMPLEX_CTX_DATA_NLIMBS(ctx);

    for (i = 0; i < n; i++)
        NFLOAT_DATA(res)[i] = NFLOAT_DATA(x)[i];

    return GR_SUCCESS;
}

int
nfloat_complex_set_other(nfloat_complex_ptr res, gr_srcptr x, gr_ctx_t x_ctx, gr_ctx_t ctx)
{
    int status = GR_SUCCESS;

    switch (x_ctx->which_ring)
    {
        case GR_CTX_NFLOAT:
            status = nfloat_set_other(res, x, x_ctx, ctx);
            nfloat_zero(NFLOAT_COMPLEX_IM(res, ctx), ctx);
            return status;

        case GR_CTX_FMPZ:
            status = nfloat_set_fmpz(res, x, ctx);
            nfloat_zero(NFLOAT_COMPLEX_IM(res, ctx), ctx);
            return status;

        case GR_CTX_FMPQ:
            status = nfloat_set_fmpq(res, x, ctx);
            nfloat_zero(NFLOAT_COMPLEX_IM(res, ctx), ctx);
            return status;

        case GR_CTX_REAL_FLOAT_ARF:
            status = nfloat_set_arf(res, x, ctx);
            nfloat_zero(NFLOAT_COMPLEX_IM(res, ctx), ctx);
            return status;

        case GR_CTX_RR_ARB:
            status = nfloat_set_arf(res, arb_midref((arb_srcptr) x), ctx);
            nfloat_zero(NFLOAT_COMPLEX_IM(res, ctx), ctx);
            return status;

        case GR_CTX_COMPLEX_FLOAT_ACF:
            status = nfloat_complex_set_acf(res, x, ctx);
            return status;

        case GR_CTX_CC_ACB:
            status = nfloat_complex_set_acb(res, x, ctx);
            return status;

        default:
            {
                int status;
                acf_t t;

                gr_ctx_t acf_ctx;
                acf_init(t);

                gr_ctx_init_complex_float_acf(acf_ctx, NFLOAT_CTX_PREC(ctx));
                status = gr_set_other(t, x, x_ctx, acf_ctx);
                if (status == GR_SUCCESS)
                    status = nfloat_complex_set_acf(res, t, ctx);

                acf_clear(t);
                gr_ctx_clear(acf_ctx);
                return status;
            }
    }
}

int
nfloat_complex_one(nfloat_complex_ptr res, gr_ctx_t ctx)
{
    nfloat_one(NFLOAT_COMPLEX_RE(res, ctx), ctx);
    nfloat_zero(NFLOAT_COMPLEX_IM(res, ctx), ctx);
    return GR_SUCCESS;
}

int
nfloat_complex_neg_one(nfloat_complex_ptr res, gr_ctx_t ctx)
{
    nfloat_neg_one(NFLOAT_COMPLEX_RE(res, ctx), ctx);
    nfloat_zero(NFLOAT_COMPLEX_IM(res, ctx), ctx);
    return GR_SUCCESS;
}

truth_t
nfloat_complex_is_zero(nfloat_complex_srcptr x, gr_ctx_t ctx)
{
    return truth_and(nfloat_is_zero(NFLOAT_COMPLEX_RE(x, ctx), ctx),
                     nfloat_is_zero(NFLOAT_COMPLEX_IM(x, ctx), ctx));
}

truth_t
nfloat_complex_is_one(nfloat_complex_srcptr x, gr_ctx_t ctx)
{
    return truth_and(nfloat_is_one(NFLOAT_COMPLEX_RE(x, ctx), ctx),
                     nfloat_is_zero(NFLOAT_COMPLEX_IM(x, ctx), ctx));
}

truth_t
nfloat_complex_is_neg_one(nfloat_complex_srcptr x, gr_ctx_t ctx)
{
    return truth_and(nfloat_is_neg_one(NFLOAT_COMPLEX_RE(x, ctx), ctx),
                     nfloat_is_zero(NFLOAT_COMPLEX_IM(x, ctx), ctx));
}

int
nfloat_complex_i(nfloat_complex_ptr res, gr_ctx_t ctx)
{
    nfloat_zero(NFLOAT_COMPLEX_RE(res, ctx), ctx);
    nfloat_one(NFLOAT_COMPLEX_IM(res, ctx), ctx);
    return GR_SUCCESS;
}

int
nfloat_complex_pi(nfloat_complex_ptr res, gr_ctx_t ctx)
{
    nfloat_pi(NFLOAT_COMPLEX_RE(res, ctx), ctx);
    nfloat_zero(NFLOAT_COMPLEX_IM(res, ctx), ctx);
    return GR_SUCCESS;
}

/* todo: be smart when in-place */
int
nfloat_complex_conj(nfloat_complex_ptr res, nfloat_complex_srcptr x, gr_ctx_t ctx)
{
    nfloat_set(NFLOAT_COMPLEX_RE(res, ctx), NFLOAT_COMPLEX_RE(x, ctx), ctx);
    nfloat_neg(NFLOAT_COMPLEX_IM(res, ctx), NFLOAT_COMPLEX_IM(x, ctx), ctx);
    return GR_SUCCESS;
}

int
nfloat_complex_re(nfloat_complex_ptr res, nfloat_complex_srcptr x, gr_ctx_t ctx)
{
    nfloat_set(NFLOAT_COMPLEX_RE(res, ctx), NFLOAT_COMPLEX_RE(x, ctx), ctx);
    nfloat_zero(NFLOAT_COMPLEX_IM(res, ctx), ctx);
    return GR_SUCCESS;
}

int
nfloat_complex_im(nfloat_complex_ptr res, nfloat_complex_srcptr x, gr_ctx_t ctx)
{
    nfloat_set(NFLOAT_COMPLEX_RE(res, ctx), NFLOAT_COMPLEX_IM(x, ctx), ctx);
    nfloat_zero(NFLOAT_COMPLEX_IM(res, ctx), ctx);
    return GR_SUCCESS;
}

truth_t
nfloat_complex_equal(nfloat_complex_srcptr x, nfloat_complex_srcptr y, gr_ctx_t ctx)
{
    return truth_and(nfloat_equal(NFLOAT_COMPLEX_RE(x, ctx), NFLOAT_COMPLEX_RE(y, ctx), ctx),
                     nfloat_equal(NFLOAT_COMPLEX_IM(x, ctx), NFLOAT_COMPLEX_IM(y, ctx), ctx));
}

int
nfloat_complex_set_si(nfloat_complex_ptr res, slong x, gr_ctx_t ctx)
{
    nfloat_zero(NFLOAT_COMPLEX_IM(res, ctx), ctx);
    return nfloat_set_si(NFLOAT_COMPLEX_RE(res, ctx), x, ctx);
}

int
nfloat_complex_set_ui(nfloat_complex_ptr res, ulong x, gr_ctx_t ctx)
{
    nfloat_zero(NFLOAT_COMPLEX_IM(res, ctx), ctx);
    return nfloat_set_ui(NFLOAT_COMPLEX_RE(res, ctx), x, ctx);
}

int
nfloat_complex_set_fmpz(nfloat_complex_ptr res, const fmpz_t x, gr_ctx_t ctx)
{
    nfloat_zero(NFLOAT_COMPLEX_IM(res, ctx), ctx);
    return nfloat_set_fmpz(NFLOAT_COMPLEX_RE(res, ctx), x, ctx);
}

int
nfloat_complex_set_fmpq(nfloat_complex_ptr res, const fmpq_t x, gr_ctx_t ctx)
{
    nfloat_zero(NFLOAT_COMPLEX_IM(res, ctx), ctx);
    return nfloat_set_fmpq(NFLOAT_COMPLEX_RE(res, ctx), x, ctx);
}

int
nfloat_complex_set_d(nfloat_complex_ptr res, double x, gr_ctx_t ctx)
{
    nfloat_zero(NFLOAT_COMPLEX_IM(res, ctx), ctx);
    return nfloat_set_d(NFLOAT_COMPLEX_RE(res, ctx), x, ctx);
}

int
nfloat_complex_neg(nfloat_complex_ptr res, nfloat_complex_srcptr x, gr_ctx_t ctx)
{
    nfloat_neg(NFLOAT_COMPLEX_RE(res, ctx), NFLOAT_COMPLEX_RE(x, ctx), ctx);
    nfloat_neg(NFLOAT_COMPLEX_IM(res, ctx), NFLOAT_COMPLEX_IM(x, ctx), ctx);
    return GR_SUCCESS;
}

int
nfloat_complex_add(nfloat_complex_ptr res, nfloat_complex_srcptr x, nfloat_complex_srcptr y, gr_ctx_t ctx)
{
    int status = GR_SUCCESS;
    status |= nfloat_add(NFLOAT_COMPLEX_RE(res, ctx), NFLOAT_COMPLEX_RE(x, ctx), NFLOAT_COMPLEX_RE(y, ctx), ctx);
    status |= nfloat_add(NFLOAT_COMPLEX_IM(res, ctx), NFLOAT_COMPLEX_IM(x, ctx), NFLOAT_COMPLEX_IM(y, ctx), ctx);
    return status;
}

int
nfloat_complex_sub(nfloat_complex_ptr res, nfloat_complex_srcptr x, nfloat_complex_srcptr y, gr_ctx_t ctx)
{
    int status = GR_SUCCESS;
    status |= nfloat_sub(NFLOAT_COMPLEX_RE(res, ctx), NFLOAT_COMPLEX_RE(x, ctx), NFLOAT_COMPLEX_RE(y, ctx), ctx);
    status |= nfloat_sub(NFLOAT_COMPLEX_IM(res, ctx), NFLOAT_COMPLEX_IM(x, ctx), NFLOAT_COMPLEX_IM(y, ctx), ctx);
    return status;
}

static inline int
_flint_mpn_cmp_2(ulong a1, ulong a0, ulong b1, ulong b0)
{
    if (a1 != b1) return (a1 < b1) ? -1 : 1;
    if (a0 != b0) return (a0 < b0) ? -1 : 1;
    return 0;
}

static inline int
_flint_mpn_cmp_3(ulong a2, ulong a1, ulong a0, ulong b2, ulong b1, ulong b0)
{
    if (a2 != b2) return (a2 < b2) ? -1 : 1;
    if (a1 != b1) return (a1 < b1) ? -1 : 1;
    if (a0 != b0) return (a0 < b0) ? -1 : 1;
    return 0;
}

int
_nfloat_complex_sqr_naive(nfloat_ptr res1, nfloat_ptr res2, nfloat_srcptr a, nfloat_srcptr b, gr_ctx_t ctx)
{
    ulong a2[NFLOAT_MAX_ALLOC];
    ulong b2[NFLOAT_MAX_ALLOC];
    int status = GR_SUCCESS;

    status |= nfloat_mul(a2, a, a, ctx);
    status |= nfloat_mul(b2, b, b, ctx);
    status |= nfloat_mul(res2, a, b, ctx);
    status |= nfloat_mul_2exp_si(res2, res2, 1, ctx);
    status |= nfloat_sub(res1, a2, b2, ctx);

    return status;
}

int
_nfloat_complex_sqr_standard(nfloat_ptr res1, nfloat_ptr res2, nfloat_srcptr a, nfloat_srcptr b, gr_ctx_t ctx)
{
    int status = GR_SUCCESS;
    ulong a2[NFLOAT_MAX_LIMBS + 1];
    ulong b2[NFLOAT_MAX_LIMBS + 1];
    ulong hi, lo;
    int ssgnbit;
    slong a2exp, b2exp, a2b2exp, delta;
    slong n, norm;

    FLINT_ASSERT(!NFLOAT_IS_SPECIAL(a));
    FLINT_ASSERT(!NFLOAT_IS_SPECIAL(b));

    n = NFLOAT_CTX_NLIMBS(ctx);

    a2exp = 2 * NFLOAT_EXP(a);
    b2exp = 2 * NFLOAT_EXP(b);

    a2b2exp = FLINT_MAX(a2exp, b2exp);
    delta = a2b2exp - FLINT_MIN(a2exp, b2exp);

    /* TODO: this case is rare, but we could optimize it too */
    if (delta >= FLINT_BITS)
        return _nfloat_complex_sqr_naive(res1, res2, a, b, ctx);

    if (n == 1)
    {
        umul_ppmm(a2[1], a2[0], NFLOAT_D(a)[0], NFLOAT_D(a)[0]);
        umul_ppmm(b2[1], b2[0], NFLOAT_D(b)[0], NFLOAT_D(b)[0]);

        if (a2exp == b2exp)
        {
            if (_flint_mpn_cmp_2(a2[1], a2[0], b2[1], b2[0]) >= 0)
            {
                sub_ddmmss(a2[1], a2[0], a2[1], a2[0], b2[1], b2[0]);
                ssgnbit = 0;
            }
            else
            {
                sub_ddmmss(a2[1], a2[0], b2[1], b2[0], a2[1], a2[0]);
                ssgnbit = 1;
            }
        }
        else if (a2exp > b2exp)
        {
            b2[0] = (b2[0] >> delta) | (b2[1] << (FLINT_BITS - delta));
            b2[1] = (b2[1] >> delta);
            sub_ddmmss(a2[1], a2[0], a2[1], a2[0], b2[1], b2[0]);
            ssgnbit = 0;
        }
        else
        {
            a2[0] = (a2[0] >> delta) | (a2[1] << (FLINT_BITS - delta));
            a2[1] = (a2[1] >> delta);
            sub_ddmmss(a2[1], a2[0], b2[1], b2[0], a2[1], a2[0]);
            ssgnbit = 1;
        }

        umul_ppmm(hi, lo, NFLOAT_D(a)[0], NFLOAT_D(b)[0]);

        if (LIMB_MSB_IS_SET(hi))
        {
            NFLOAT_D(res2)[0] = hi;
            NFLOAT_EXP(res2) = NFLOAT_EXP(a) + NFLOAT_EXP(b) + 1;
        }
        else
        {
            NFLOAT_D(res2)[0] = (hi << 1) | (lo >> (FLINT_BITS - 1));
            NFLOAT_EXP(res2) = NFLOAT_EXP(a) + NFLOAT_EXP(b);
        }

        NFLOAT_SGNBIT(res2) = NFLOAT_SGNBIT(a) ^ NFLOAT_SGNBIT(b);
        NFLOAT_HANDLE_UNDERFLOW_OVERFLOW(res2, ctx);

        if (a2[1] == 0)
        {
            if (a2[0] == 0)
            {
                status |= nfloat_zero(res1, ctx);
                return status;
            }
            else
            {
                norm = flint_clz(a2[0]);
                a2[1] = a2[0] << norm;
                a2b2exp -= FLINT_BITS + norm;
            }
        }
        else
        {
            norm = flint_clz(a2[1]);

            if (norm != 0)
            {
                a2[1] = (a2[1] << norm) | (a2[0] >> (FLINT_BITS - norm));
                a2b2exp -= norm;
            }
        }

        NFLOAT_SGNBIT(res1) = ssgnbit;
        NFLOAT_EXP(res1) = a2b2exp;
        NFLOAT_D(res1)[0] = a2[1];
        NFLOAT_HANDLE_UNDERFLOW_OVERFLOW(res1, ctx);
        return status;
    }

    if (n == 2)
    {
        ulong FLINT_SET_BUT_UNUSED(uu);

        FLINT_MPN_SQR_2X2(a2[2], a2[1], a2[0], uu, NFLOAT_D(a)[1], NFLOAT_D(a)[0]);
        FLINT_MPN_SQR_2X2(b2[2], b2[1], b2[0], uu, NFLOAT_D(b)[1], NFLOAT_D(b)[0]);

        if (a2exp == b2exp)
        {
            a2b2exp = a2exp;

            if (_flint_mpn_cmp_3(a2[2], a2[1], a2[0], b2[2], b2[1], b2[0]) >= 0)
            {
                sub_dddmmmsss(a2[2], a2[1], a2[0], a2[2], a2[1], a2[0], b2[2], b2[1], b2[0]);
                ssgnbit = 0;
            }
            else
            {
                sub_dddmmmsss(a2[2], a2[1], a2[0], b2[2], b2[1], b2[0], a2[2], a2[1], a2[0]);
                ssgnbit = 1;
            }
        }
        else if (a2exp > b2exp)
        {
            b2[0] = (b2[0] >> delta) | (b2[1] << (FLINT_BITS - delta));
            b2[1] = (b2[1] >> delta) | (b2[2] << (FLINT_BITS - delta));
            b2[2] = (b2[2] >> delta);

            sub_dddmmmsss(a2[2], a2[1], a2[0], a2[2], a2[1], a2[0], b2[2], b2[1], b2[0]);
            ssgnbit = 0;
        }
        else
        {
            a2[0] = (a2[0] >> delta) | (a2[1] << (FLINT_BITS - delta));
            a2[1] = (a2[1] >> delta) | (a2[2] << (FLINT_BITS - delta));
            a2[2] = (a2[2] >> delta);

            sub_dddmmmsss(a2[2], a2[1], a2[0], b2[2], b2[1], b2[0], a2[2], a2[1], a2[0]);
            ssgnbit = 1;
        }

        status |= nfloat_mul(res2, a, b, ctx);
        if (NFLOAT_EXP(res2) >= NFLOAT_MIN_EXP && NFLOAT_EXP(res2) < NFLOAT_MAX_EXP)
            NFLOAT_EXP(res2)++;
        else
            status |= nfloat_mul_2exp_si(res2, res2, 1, ctx);

        status |= nfloat_2_set_3_2exp(res1, a2[2], a2[1], a2[0], a2b2exp, ssgnbit, ctx);
        return status;
    }
    else
    {
        a2[0] = flint_mpn_sqrhigh_n(a2 + 1, NFLOAT_D(a), n);
        b2[0] = flint_mpn_sqrhigh_n(b2 + 1, NFLOAT_D(b), n);

        if (n == 3)
        {
            if (a2exp == b2exp)
            {
                ssgnbit = (mpn_cmp(a2, b2, n + 1) < 0);
            }
            else if (a2exp > b2exp)
            {
                b2[0] = (b2[0] >> delta) | (b2[1] << (FLINT_BITS - delta));
                b2[1] = (b2[1] >> delta) | (b2[2] << (FLINT_BITS - delta));
                b2[2] = (b2[2] >> delta) | (b2[3] << (FLINT_BITS - delta));
                b2[3] = (b2[3] >> delta);
                ssgnbit = 0;
            }
            else
            {
                a2[0] = (a2[0] >> delta) | (a2[1] << (FLINT_BITS - delta));
                a2[1] = (a2[1] >> delta) | (a2[2] << (FLINT_BITS - delta));
                a2[2] = (a2[2] >> delta) | (a2[3] << (FLINT_BITS - delta));
                a2[3] = (a2[3] >> delta);
                ssgnbit = 1;
            }

            if (ssgnbit == 0)
                sub_ddddmmmmssss(a2[3], a2[2], a2[1], a2[0], a2[3], a2[2], a2[1], a2[0], b2[3], b2[2], b2[1], b2[0]);
            else
                sub_ddddmmmmssss(a2[3], a2[2], a2[1], a2[0], b2[3], b2[2], b2[1], b2[0], a2[3], a2[2], a2[1], a2[0]);
        }
        else if (n == 4)
        {
            if (a2exp == b2exp)
            {
                ssgnbit = (mpn_cmp(a2, b2, n + 1) < 0);
            }
            else if (a2exp > b2exp)
            {
                b2[0] = (b2[0] >> delta) | (b2[1] << (FLINT_BITS - delta));
                b2[1] = (b2[1] >> delta) | (b2[2] << (FLINT_BITS - delta));
                b2[2] = (b2[2] >> delta) | (b2[3] << (FLINT_BITS - delta));
                b2[3] = (b2[3] >> delta) | (b2[4] << (FLINT_BITS - delta));
                b2[4] = (b2[4] >> delta);
                ssgnbit = 0;
            }
            else
            {
                a2[0] = (a2[0] >> delta) | (a2[1] << (FLINT_BITS - delta));
                a2[1] = (a2[1] >> delta) | (a2[2] << (FLINT_BITS - delta));
                a2[2] = (a2[2] >> delta) | (a2[3] << (FLINT_BITS - delta));
                a2[3] = (a2[3] >> delta) | (a2[4] << (FLINT_BITS - delta));
                a2[4] = (a2[4] >> delta);
                ssgnbit = 1;
            }

            if (ssgnbit == 0)
                sub_dddddmmmmmsssss(a2[4], a2[3], a2[2], a2[1], a2[0], a2[4], a2[3], a2[2], a2[1], a2[0], b2[4], b2[3], b2[2], b2[1], b2[0]);
            else
                sub_dddddmmmmmsssss(a2[4], a2[3], a2[2], a2[1], a2[0], b2[4], b2[3], b2[2], b2[1], b2[0], a2[4], a2[3], a2[2], a2[1], a2[0]);
        }
        else
        {
            if (a2exp == b2exp)
            {
                ssgnbit = flint_mpn_signed_sub_n(a2, a2, b2, n + 1);
            }
            else if (a2exp > b2exp)
            {
                mpn_rshift(b2, b2, n + 1, a2exp - b2exp);
                mpn_sub_n(a2, a2, b2, n + 1);
                ssgnbit = 0;
            }
            else
            {
                mpn_rshift(a2, a2, n + 1, b2exp - a2exp);
                mpn_sub_n(a2, b2, a2, n + 1);
                ssgnbit = 1;
            }
        }
    }

    status |= nfloat_mul(res2, a, b, ctx);
    if (NFLOAT_EXP(res2) >= NFLOAT_MIN_EXP && NFLOAT_EXP(res2) < NFLOAT_MAX_EXP)
        NFLOAT_EXP(res2)++;
    else
        status |= nfloat_mul_2exp_si(res2, res2, 1, ctx);

    status |= nfloat_set_mpn_2exp(res1, a2, n + 1, a2b2exp, ssgnbit, ctx);

    return status;
}

int
_nfloat_complex_sqr_karatsuba(nfloat_ptr res1, nfloat_ptr res2, nfloat_srcptr a, nfloat_srcptr b, gr_ctx_t ctx)
{
    slong aexp, bexp, abexp, adelta, bdelta;
    int asgnbit, bsgnbit, ssgnbit, tsgnbit;
    int status;
    slong n;
    ulong aa[NFLOAT_MAX_LIMBS + 1];
    ulong bb[NFLOAT_MAX_LIMBS + 1];
    ulong s[NFLOAT_MAX_LIMBS + 1];
    ulong t[NFLOAT_MAX_LIMBS + 1];
    ulong u[NFLOAT_MAX_LIMBS + 1];
    ulong v[NFLOAT_MAX_LIMBS + 1];

    FLINT_ASSERT(!NFLOAT_IS_SPECIAL(a));
    FLINT_ASSERT(!NFLOAT_IS_SPECIAL(b));

    n = NFLOAT_CTX_NLIMBS(ctx);

    aexp = NFLOAT_EXP(a);
    bexp = NFLOAT_EXP(b);

    asgnbit = NFLOAT_SGNBIT(a);
    bsgnbit = NFLOAT_SGNBIT(b);

    /* Two extra bits allows adding without overflow. */
    abexp = FLINT_MAX(aexp, bexp) + 2;
    adelta = abexp - aexp;
    bdelta = abexp - bexp;

    /* We use one guard limb. There is cancellation of about
       max(adelta,bdelta) bits; abandon if we are running out
       of guard bits. */
    if (adelta >= FLINT_BITS - 4 || bdelta >= FLINT_BITS - 4)
        return _nfloat_complex_sqr_standard(res1, res2, a, b, ctx);

    aa[0] = mpn_rshift(aa + 1, NFLOAT_D(a), n, adelta);
    bb[0] = mpn_rshift(bb + 1, NFLOAT_D(b), n, bdelta);
    _flint_mpn_signed_add_n(v, aa, asgnbit, bb, bsgnbit, n + 1);
    flint_mpn_sqrhigh_n(s, v, n + 1);
    flint_mpn_sqrhigh_n(t, aa, n + 1);
    flint_mpn_sqrhigh_n(u, bb, n + 1);
    mpn_add_n(v, t, u, n + 1);
    ssgnbit = flint_mpn_signed_sub_n(s, s, v, n + 1);
    tsgnbit = flint_mpn_signed_sub_n(t, t, u, n + 1);

    status = GR_SUCCESS;
    status |= nfloat_set_mpn_2exp(res1, t, n + 1, 2 * abexp, tsgnbit, ctx);
    status |= nfloat_set_mpn_2exp(res2, s, n + 1, 2 * abexp, ssgnbit, ctx);
    return status;
}

int
_nfloat_complex_sqr(nfloat_ptr res1, nfloat_ptr res2, nfloat_srcptr a, nfloat_srcptr b, gr_ctx_t ctx)
{
    int status;

    if (NFLOAT_CTX_HAS_INF_NAN(ctx))
        return _nfloat_complex_sqr_naive(res1, res2, a, b, ctx);

    if (NFLOAT_IS_ZERO(b))
    {
        status = nfloat_sqr(res1, a, ctx);
        status |= nfloat_zero(res2, ctx);
        return status;
    }

    if (NFLOAT_IS_ZERO(a))
    {
        status = nfloat_sqr(res1, b, ctx);
        status |= nfloat_neg(res1, res1, ctx);
        status |= nfloat_zero(res2, ctx);
        return status;
    }

    if (NFLOAT_CTX_NLIMBS(ctx) < 20)
        return _nfloat_complex_sqr_standard(res1, res2, a, b, ctx);
    else
        return _nfloat_complex_sqr_karatsuba(res1, res2, a, b, ctx);
}

int
nfloat_complex_sqr(nfloat_complex_ptr res, nfloat_complex_srcptr x, gr_ctx_t ctx)
{
    return _nfloat_complex_sqr(NFLOAT_COMPLEX_RE(res, ctx),
                               NFLOAT_COMPLEX_IM(res, ctx),
                               NFLOAT_COMPLEX_RE(x, ctx),
                               NFLOAT_COMPLEX_IM(x, ctx), ctx);
}

int
_nfloat_complex_mul_naive(nfloat_ptr res1, nfloat_ptr res2, nfloat_srcptr a, nfloat_srcptr b, nfloat_srcptr c, nfloat_srcptr d, gr_ctx_t ctx)
{
    ulong ac[NFLOAT_MAX_ALLOC];
    ulong bd[NFLOAT_MAX_ALLOC];
    ulong ad[NFLOAT_MAX_ALLOC];
    ulong bc[NFLOAT_MAX_ALLOC];
    int status = GR_SUCCESS;

    status |= nfloat_mul(ac, a, c, ctx);
    status |= nfloat_mul(bd, b, d, ctx);
    status |= nfloat_mul(ad, a, d, ctx);
    status |= nfloat_mul(bc, b, c, ctx);
    status |= nfloat_sub(res1, ac, bd, ctx);
    status |= nfloat_add(res2, ad, bc, ctx);

    return status;
}

int
_nfloat_complex_mul_standard(nfloat_ptr res1, nfloat_ptr res2, nfloat_srcptr a, nfloat_srcptr b, nfloat_srcptr c, nfloat_srcptr d, gr_ctx_t ctx)
{
    slong aexp, bexp, cexp, dexp, acexp, bdexp, adexp, bcexp;
    slong sexp, texp, sdelta, tdelta;
    int asgnbit, bsgnbit, csgnbit, dsgnbit, usgnbit, vsgnbit, ssgnbit, tsgnbit;
    int status;
    slong n;
    ulong u[NFLOAT_MAX_LIMBS + 1];
    ulong v[NFLOAT_MAX_LIMBS + 1];
    ulong s[NFLOAT_MAX_LIMBS + 2];
    ulong t[NFLOAT_MAX_LIMBS + 2];

    FLINT_ASSERT(!NFLOAT_IS_SPECIAL(a));
    FLINT_ASSERT(!NFLOAT_IS_SPECIAL(b));
    FLINT_ASSERT(!NFLOAT_IS_SPECIAL(c));
    FLINT_ASSERT(!NFLOAT_IS_SPECIAL(d));

    n = NFLOAT_CTX_NLIMBS(ctx);

    aexp = NFLOAT_EXP(a);
    bexp = NFLOAT_EXP(b);
    cexp = NFLOAT_EXP(c);
    dexp = NFLOAT_EXP(d);

    asgnbit = NFLOAT_SGNBIT(a);
    bsgnbit = NFLOAT_SGNBIT(b);
    csgnbit = NFLOAT_SGNBIT(c);
    dsgnbit = NFLOAT_SGNBIT(d);

    /* ac - bd */
    acexp = aexp + cexp;
    bdexp = bexp + dexp;
    sexp = FLINT_MAX(acexp, bdexp);
    sdelta = sexp - FLINT_MIN(acexp, bdexp);

    /* ad + bc */
    adexp = aexp + dexp;
    bcexp = bexp + cexp;
    texp = FLINT_MAX(adexp, bcexp);
    tdelta = texp - FLINT_MIN(adexp, bcexp);

    /* todo */
    if (sdelta >= FLINT_BITS || tdelta >= FLINT_BITS)
        return _nfloat_complex_mul_naive(res1, res2, a, b, c, d, ctx);

    if (n == 1)
    {
        umul_ppmm(u[1], u[0], NFLOAT_D(a)[0], NFLOAT_D(c)[0]);
        umul_ppmm(v[1], v[0], NFLOAT_D(b)[0], NFLOAT_D(d)[0]);
        usgnbit = asgnbit ^ csgnbit;
        vsgnbit = !(bsgnbit ^ dsgnbit);

        if (sdelta != 0)
        {
            if (acexp > bdexp)
            {
                v[0] = (v[0] >> sdelta) | (v[1] << (FLINT_BITS - sdelta));
                v[1] >>= sdelta;
            }
            else
            {
                u[0] = (u[0] >> sdelta) | (u[1] << (FLINT_BITS - sdelta));
                u[1] >>= sdelta;
            }
        }

        if (usgnbit == vsgnbit)
        {
            add_sssaaaaaa(s[2], s[1], s[0], 0, u[1], u[0], 0, v[1], v[0]);
            ssgnbit = usgnbit;
        }
        else
        {
            if (_flint_mpn_cmp_2(u[1], u[0], v[1], v[0]) >= 0)
            {
                sub_ddmmss(s[1], s[0], u[1], u[0], v[1], v[0]);
                ssgnbit = usgnbit;
            }
            else
            {
                sub_ddmmss(s[1], s[0], v[1], v[0], u[1], u[0]);
                ssgnbit = !usgnbit;
            }

            s[2] = 0;
        }

        umul_ppmm(u[1], u[0], NFLOAT_D(a)[0], NFLOAT_D(d)[0]);
        umul_ppmm(v[1], v[0], NFLOAT_D(b)[0], NFLOAT_D(c)[0]);
        usgnbit = asgnbit ^ dsgnbit;
        vsgnbit = bsgnbit ^ csgnbit;

        if (tdelta != 0)
        {
            if (adexp > bcexp)
            {
                v[0] = (v[0] >> tdelta) | (v[1] << (FLINT_BITS - tdelta));
                v[1] >>= tdelta;
            }
            else
            {
                u[0] = (u[0] >> tdelta) | (u[1] << (FLINT_BITS - tdelta));
                u[1] >>= tdelta;
            }
        }

        if (usgnbit == vsgnbit)
        {
            add_sssaaaaaa(t[2], t[1], t[0], 0, u[1], u[0], 0, v[1], v[0]);
            tsgnbit = usgnbit;
        }
        else
        {
            if (_flint_mpn_cmp_2(u[1], u[0], v[1], v[0]) >= 0)
            {
                sub_ddmmss(t[1], t[0], u[1], u[0], v[1], v[0]);
                tsgnbit = usgnbit;
            }
            else
            {
                sub_ddmmss(t[1], t[0], v[1], v[0], u[1], u[0]);
                tsgnbit = !usgnbit;
            }

            t[2] = 0;
        }

        status = GR_SUCCESS;
        status |= nfloat_1_set_3_2exp(res1, s[2], s[1], s[0], sexp + FLINT_BITS, ssgnbit, ctx);
        status |= nfloat_1_set_3_2exp(res2, t[2], t[1], t[0], texp + FLINT_BITS, tsgnbit, ctx);
        return status;
    }
    else if (n == 2)
    {
        ulong FLINT_SET_BUT_UNUSED(uu);

        FLINT_MPN_MUL_2X2(u[2], u[1], u[0], uu, NFLOAT_D(a)[1], NFLOAT_D(a)[0], NFLOAT_D(c)[1], NFLOAT_D(c)[0]);
        FLINT_MPN_MUL_2X2(v[2], v[1], v[0], uu, NFLOAT_D(b)[1], NFLOAT_D(b)[0], NFLOAT_D(d)[1], NFLOAT_D(d)[0]);
        usgnbit = asgnbit ^ csgnbit;
        vsgnbit = !(bsgnbit ^ dsgnbit);

        if (sdelta != 0)
        {
            if (acexp > bdexp)
            {
                v[0] = (v[0] >> sdelta) | (v[1] << (FLINT_BITS - sdelta));
                v[1] = (v[1] >> sdelta) | (v[2] << (FLINT_BITS - sdelta));
                v[2] >>= sdelta;
            }
            else
            {
                u[0] = (u[0] >> sdelta) | (u[1] << (FLINT_BITS - sdelta));
                u[1] = (u[1] >> sdelta) | (u[2] << (FLINT_BITS - sdelta));
                u[2] >>= sdelta;
            }
        }

        if (usgnbit == vsgnbit)
        {
            add_ssssaaaaaaaa(s[3], s[2], s[1], s[0], 0, u[2], u[1], u[0], 0, v[2], v[1], v[0]);
            ssgnbit = usgnbit;
        }
        else
        {
            if (_flint_mpn_cmp_3(u[2], u[1], u[0], v[2], v[1], v[0]) >= 0)
            {
                sub_dddmmmsss(s[2], s[1], s[0], u[2], u[1], u[0], v[2], v[1], v[0]);
                ssgnbit = usgnbit;
            }
            else
            {
                sub_dddmmmsss(s[2], s[1], s[0], v[2], v[1], v[0], u[2], u[1], u[0]);
                ssgnbit = !usgnbit;
            }

            s[3] = 0;
        }

        FLINT_MPN_MUL_2X2(u[2], u[1], u[0], uu, NFLOAT_D(a)[1], NFLOAT_D(a)[0], NFLOAT_D(d)[1], NFLOAT_D(d)[0]);
        FLINT_MPN_MUL_2X2(v[2], v[1], v[0], uu, NFLOAT_D(b)[1], NFLOAT_D(b)[0], NFLOAT_D(c)[1], NFLOAT_D(c)[0]);
        usgnbit = asgnbit ^ dsgnbit;
        vsgnbit = bsgnbit ^ csgnbit;

        if (tdelta != 0)
        {
            if (adexp > bcexp)
            {
                v[0] = (v[0] >> tdelta) | (v[1] << (FLINT_BITS - tdelta));
                v[1] = (v[1] >> tdelta) | (v[2] << (FLINT_BITS - tdelta));
                v[2] >>= tdelta;
            }
            else
            {
                u[0] = (u[0] >> tdelta) | (u[1] << (FLINT_BITS - tdelta));
                u[1] = (u[1] >> tdelta) | (u[2] << (FLINT_BITS - tdelta));
                u[2] >>= tdelta;
            }
        }

        if (usgnbit == vsgnbit)
        {
            add_ssssaaaaaaaa(t[3], t[2], t[1], t[0], 0, u[2], u[1], u[0], 0, v[2], v[1], v[0]);
            tsgnbit = usgnbit;
        }
        else
        {
            if (_flint_mpn_cmp_3(u[2], u[1], u[0], v[2], v[1], v[0]) >= 0)
            {
                sub_dddmmmsss(t[2], t[1], t[0], u[2], u[1], u[0], v[2], v[1], v[0]);
                tsgnbit = usgnbit;
            }
            else
            {
                sub_dddmmmsss(t[2], t[1], t[0], v[2], v[1], v[0], u[2], u[1], u[0]);
                tsgnbit = !usgnbit;
            }

            t[3] = 0;
        }

        status = GR_SUCCESS;
        status |= nfloat_2_set_4_2exp(res1, s[3], s[2], s[1], s[0], sexp + FLINT_BITS, ssgnbit, ctx);
        status |= nfloat_2_set_4_2exp(res2, t[3], t[2], t[1], t[0], texp + FLINT_BITS, tsgnbit, ctx);
        return status;
    }
    else if (n == 3)
    {
        u[0] = flint_mpn_mulhigh_n(u + 1, NFLOAT_D(a), NFLOAT_D(c), n);
        v[0] = flint_mpn_mulhigh_n(v + 1, NFLOAT_D(b), NFLOAT_D(d), n);
        usgnbit = asgnbit ^ csgnbit;
        vsgnbit = !(bsgnbit ^ dsgnbit);

        if (sdelta != 0)
        {
            if (acexp > bdexp)
            {
                v[0] = (v[0] >> sdelta) | (v[1] << (FLINT_BITS - sdelta));
                v[1] = (v[1] >> sdelta) | (v[2] << (FLINT_BITS - sdelta));
                v[2] = (v[2] >> sdelta) | (v[3] << (FLINT_BITS - sdelta));
                v[3] >>= sdelta;
            }
            else
            {
                u[0] = (u[0] >> sdelta) | (u[1] << (FLINT_BITS - sdelta));
                u[1] = (u[1] >> sdelta) | (u[2] << (FLINT_BITS - sdelta));
                u[2] = (u[2] >> sdelta) | (u[3] << (FLINT_BITS - sdelta));
                u[3] >>= sdelta;
            }
        }

        if (usgnbit == vsgnbit)
        {
            add_sssssaaaaaaaaaa(s[4], s[3], s[2], s[1], s[0], 0, u[3], u[2], u[1], u[0], 0, v[3], v[2], v[1], v[0]);
            ssgnbit = usgnbit;
        }
        else
        {
            if (mpn_cmp(u, v, 4) >= 0)
            {
                sub_ddddmmmmssss(s[3], s[2], s[1], s[0], u[3], u[2], u[1], u[0], v[3], v[2], v[1], v[0]);
                ssgnbit = usgnbit;
            }
            else
            {
                sub_ddddmmmmssss(s[3], s[2], s[1], s[0], v[3], v[2], v[1], v[0], u[3], u[2], u[1], u[0]);
                ssgnbit = !usgnbit;
            }

            s[4] = 0;
        }

        u[0] = flint_mpn_mulhigh_n(u + 1, NFLOAT_D(a), NFLOAT_D(d), n);
        v[0] = flint_mpn_mulhigh_n(v + 1, NFLOAT_D(b), NFLOAT_D(c), n);
        usgnbit = asgnbit ^ dsgnbit;
        vsgnbit = bsgnbit ^ csgnbit;

        if (tdelta != 0)
        {
            if (adexp > bcexp)
            {
                v[0] = (v[0] >> tdelta) | (v[1] << (FLINT_BITS - tdelta));
                v[1] = (v[1] >> tdelta) | (v[2] << (FLINT_BITS - tdelta));
                v[2] = (v[2] >> tdelta) | (v[3] << (FLINT_BITS - tdelta));
                v[3] >>= tdelta;
            }
            else
            {
                u[0] = (u[0] >> tdelta) | (u[1] << (FLINT_BITS - tdelta));
                u[1] = (u[1] >> tdelta) | (u[2] << (FLINT_BITS - tdelta));
                u[2] = (u[2] >> tdelta) | (u[3] << (FLINT_BITS - tdelta));
                u[3] >>= tdelta;
            }
        }

        if (usgnbit == vsgnbit)
        {
            add_sssssaaaaaaaaaa(t[4], t[3], t[2], t[1], t[0], 0, u[3], u[2], u[1], u[0], 0, v[3], v[2], v[1], v[0]);
            tsgnbit = usgnbit;
        }
        else
        {
            if (mpn_cmp(u, v, 4) >= 0)
            {
                sub_ddddmmmmssss(t[3], t[2], t[1], t[0], u[3], u[2], u[1], u[0], v[3], v[2], v[1], v[0]);
                tsgnbit = usgnbit;
            }
            else
            {
                sub_ddddmmmmssss(t[3], t[2], t[1], t[0], v[3], v[2], v[1], v[0], u[3], u[2], u[1], u[0]);
                tsgnbit = !usgnbit;
            }

            t[4] = 0;
        }
    }
    else if (n == 4)
    {
        u[0] = flint_mpn_mulhigh_n(u + 1, NFLOAT_D(a), NFLOAT_D(c), n);
        v[0] = flint_mpn_mulhigh_n(v + 1, NFLOAT_D(b), NFLOAT_D(d), n);
        usgnbit = asgnbit ^ csgnbit;
        vsgnbit = !(bsgnbit ^ dsgnbit);

        if (sdelta != 0)
        {
            if (acexp > bdexp)
            {
                v[0] = (v[0] >> sdelta) | (v[1] << (FLINT_BITS - sdelta));
                v[1] = (v[1] >> sdelta) | (v[2] << (FLINT_BITS - sdelta));
                v[2] = (v[2] >> sdelta) | (v[3] << (FLINT_BITS - sdelta));
                v[3] = (v[3] >> sdelta) | (v[4] << (FLINT_BITS - sdelta));
                v[4] >>= sdelta;
            }
            else
            {
                u[0] = (u[0] >> sdelta) | (u[1] << (FLINT_BITS - sdelta));
                u[1] = (u[1] >> sdelta) | (u[2] << (FLINT_BITS - sdelta));
                u[2] = (u[2] >> sdelta) | (u[3] << (FLINT_BITS - sdelta));
                u[3] = (u[3] >> sdelta) | (u[4] << (FLINT_BITS - sdelta));
                u[4] >>= sdelta;
            }
        }

        if (usgnbit == vsgnbit)
        {
            add_ssssssaaaaaaaaaaaa(s[5], s[4], s[3], s[2], s[1], s[0], 0, u[4], u[3], u[2], u[1], u[0], 0, v[4], v[3], v[2], v[1], v[0]);
            ssgnbit = usgnbit;
        }
        else
        {
            if (mpn_cmp(u, v, 5) >= 0)
            {
                sub_dddddmmmmmsssss(s[4], s[3], s[2], s[1], s[0], u[4], u[3], u[2], u[1], u[0], v[4], v[3], v[2], v[1], v[0]);
                ssgnbit = usgnbit;
            }
            else
            {
                sub_dddddmmmmmsssss(s[4], s[3], s[2], s[1], s[0], v[4], v[3], v[2], v[1], v[0], u[4], u[3], u[2], u[1], u[0]);
                ssgnbit = !usgnbit;
            }

            s[5] = 0;
        }

        u[0] = flint_mpn_mulhigh_n(u + 1, NFLOAT_D(a), NFLOAT_D(d), n);
        v[0] = flint_mpn_mulhigh_n(v + 1, NFLOAT_D(b), NFLOAT_D(c), n);
        usgnbit = asgnbit ^ dsgnbit;
        vsgnbit = bsgnbit ^ csgnbit;

        if (tdelta != 0)
        {
            if (adexp > bcexp)
            {
                v[0] = (v[0] >> tdelta) | (v[1] << (FLINT_BITS - tdelta));
                v[1] = (v[1] >> tdelta) | (v[2] << (FLINT_BITS - tdelta));
                v[2] = (v[2] >> tdelta) | (v[3] << (FLINT_BITS - tdelta));
                v[3] = (v[3] >> tdelta) | (v[4] << (FLINT_BITS - tdelta));
                v[4] >>= tdelta;
            }
            else
            {
                u[0] = (u[0] >> tdelta) | (u[1] << (FLINT_BITS - tdelta));
                u[1] = (u[1] >> tdelta) | (u[2] << (FLINT_BITS - tdelta));
                u[2] = (u[2] >> tdelta) | (u[3] << (FLINT_BITS - tdelta));
                u[3] = (u[3] >> tdelta) | (u[4] << (FLINT_BITS - tdelta));
                u[4] >>= tdelta;
            }
        }

        if (usgnbit == vsgnbit)
        {
            add_ssssssaaaaaaaaaaaa(t[5], t[4], t[3], t[2], t[1], t[0], 0, u[4], u[3], u[2], u[1], u[0], 0, v[4], v[3], v[2], v[1], v[0]);
            tsgnbit = usgnbit;
        }
        else
        {
            if (mpn_cmp(u, v, 5) >= 0)
            {
                sub_dddddmmmmmsssss(t[4], t[3], t[2], t[1], t[0], u[4], u[3], u[2], u[1], u[0], v[4], v[3], v[2], v[1], v[0]);
                tsgnbit = usgnbit;
            }
            else
            {
                sub_dddddmmmmmsssss(t[4], t[3], t[2], t[1], t[0], v[4], v[3], v[2], v[1], v[0], u[4], u[3], u[2], u[1], u[0]);
                tsgnbit = !usgnbit;
            }

            t[5] = 0;
        }
    }
    else
    {
        u[0] = flint_mpn_mulhigh_n(u + 1, NFLOAT_D(a), NFLOAT_D(c), n);
        v[0] = flint_mpn_mulhigh_n(v + 1, NFLOAT_D(b), NFLOAT_D(d), n);
        usgnbit = asgnbit ^ csgnbit;
        vsgnbit = !(bsgnbit ^ dsgnbit);

        if (sdelta != 0)
        {
            if (acexp > bdexp)
                mpn_rshift(v, v, n + 1, sdelta);
            else
                mpn_rshift(u, u, n + 1, sdelta);
        }

        if (usgnbit == vsgnbit)
        {
            s[n + 1] = mpn_add_n(s, u, v, n + 1);
            ssgnbit = usgnbit;
        }
        else
        {
            ssgnbit = usgnbit ^ flint_mpn_signed_sub_n(s, u, v, n + 1);
            s[n + 1] = 0;
        }

        u[0] = flint_mpn_mulhigh_n(u + 1, NFLOAT_D(a), NFLOAT_D(d), n);
        v[0] = flint_mpn_mulhigh_n(v + 1, NFLOAT_D(b), NFLOAT_D(c), n);
        usgnbit = asgnbit ^ dsgnbit;
        vsgnbit = bsgnbit ^ csgnbit;

        if (tdelta != 0)
        {
            if (adexp > bcexp)
                mpn_rshift(v, v, n + 1, tdelta);
            else
                mpn_rshift(u, u, n + 1, tdelta);
        }

        if (usgnbit == vsgnbit)
        {
            t[n + 1] = mpn_add_n(t, u, v, n + 1);
            tsgnbit = usgnbit;
        }
        else
        {
            tsgnbit = usgnbit ^ flint_mpn_signed_sub_n(t, u, v, n + 1);
            t[n + 1] = 0;
        }
    }

    status = GR_SUCCESS;
    status |= nfloat_set_mpn_2exp(res1, s, n + 2, sexp + FLINT_BITS, ssgnbit, ctx);
    status |= nfloat_set_mpn_2exp(res2, t, n + 2, texp + FLINT_BITS, tsgnbit, ctx);

    return status;
}

int
_nfloat_complex_mul_karatsuba(nfloat_ptr res1, nfloat_ptr res2, nfloat_srcptr a, nfloat_srcptr b, nfloat_srcptr c, nfloat_srcptr d, gr_ctx_t ctx)
{
    slong aexp, bexp, cexp, dexp, abexp, cdexp, adelta, bdelta, cdelta, ddelta;
    int asgnbit, bsgnbit, csgnbit, dsgnbit, ssgnbit, tsgnbit, usgnbit;
    int status;
    slong n;

    FLINT_ASSERT(!NFLOAT_IS_SPECIAL(a));
    FLINT_ASSERT(!NFLOAT_IS_SPECIAL(b));
    FLINT_ASSERT(!NFLOAT_IS_SPECIAL(c));
    FLINT_ASSERT(!NFLOAT_IS_SPECIAL(d));

    n = NFLOAT_CTX_NLIMBS(ctx);

    aexp = NFLOAT_EXP(a);
    bexp = NFLOAT_EXP(b);
    cexp = NFLOAT_EXP(c);
    dexp = NFLOAT_EXP(d);

    asgnbit = NFLOAT_SGNBIT(a);
    bsgnbit = NFLOAT_SGNBIT(b);
    csgnbit = NFLOAT_SGNBIT(c);
    dsgnbit = NFLOAT_SGNBIT(d);

    abexp = FLINT_MAX(aexp, bexp) + 2;
    cdexp = FLINT_MAX(cexp, dexp) + 2;

    adelta = abexp - aexp;
    bdelta = abexp - bexp;
    cdelta = cdexp - cexp;
    ddelta = cdexp - dexp;

    if (adelta < FLINT_BITS && bdelta < FLINT_BITS &&
        cdelta < FLINT_BITS && ddelta < FLINT_BITS)
    {
        ulong aa[NFLOAT_MAX_LIMBS + 1];
        ulong bb[NFLOAT_MAX_LIMBS + 1];
        ulong cc[NFLOAT_MAX_LIMBS + 1];
        ulong dd[NFLOAT_MAX_LIMBS + 1];
        ulong s[NFLOAT_MAX_LIMBS + 1];
        ulong t[NFLOAT_MAX_LIMBS + 1];
        ulong u[NFLOAT_MAX_LIMBS + 1];
        ulong v[NFLOAT_MAX_LIMBS + 1];

        /*
            s = c * (a + b)
            t = a * (d - c)
            u = b * (c + d)
            re = s - u
            im = s + t
        */

        aa[0] = mpn_rshift(aa + 1, NFLOAT_D(a), n, adelta);
        bb[0] = mpn_rshift(bb + 1, NFLOAT_D(b), n, bdelta);
        cc[0] = mpn_rshift(cc + 1, NFLOAT_D(c), n, cdelta);
        dd[0] = mpn_rshift(dd + 1, NFLOAT_D(d), n, ddelta);

        ssgnbit = csgnbit ^ _flint_mpn_signed_add_n(v, aa, asgnbit, bb, bsgnbit, n + 1);
        flint_mpn_mulhigh_n(s, cc, v, n + 1);

        tsgnbit = asgnbit ^ _flint_mpn_signed_add_n(v, dd, dsgnbit, cc, !csgnbit, n + 1);
        flint_mpn_mulhigh_n(t, aa, v, n + 1);

        usgnbit = bsgnbit ^ _flint_mpn_signed_add_n(v, cc, csgnbit, dd, dsgnbit, n + 1);
        flint_mpn_mulhigh_n(u, bb, v, n + 1);

        usgnbit = _flint_mpn_signed_add_n(u, s, ssgnbit, u, !usgnbit, n + 1);
        tsgnbit = _flint_mpn_signed_add_n(t, s, ssgnbit, t, tsgnbit, n + 1);

        status = GR_SUCCESS;
        status |= nfloat_set_mpn_2exp(res1, u, n + 1, abexp + cdexp, usgnbit, ctx);
        status |= nfloat_set_mpn_2exp(res2, t, n + 1, abexp + cdexp, tsgnbit, ctx);
        return status;
    }

    return _nfloat_complex_mul_naive(res1, res2, a, b, c, d, ctx);
}

int
nfloat_complex_mul(nfloat_complex_ptr res, nfloat_complex_srcptr x, nfloat_complex_srcptr y, gr_ctx_t ctx)
{
    nfloat_srcptr a, b, c, d;
    nfloat_ptr r, s;
    int status;

    if (x == y)
        return nfloat_complex_sqr(res, x, ctx);

    r = NFLOAT_COMPLEX_RE(res, ctx);
    s = NFLOAT_COMPLEX_IM(res, ctx);

    a = NFLOAT_COMPLEX_RE(x, ctx);
    b = NFLOAT_COMPLEX_IM(x, ctx);
    c = NFLOAT_COMPLEX_RE(y, ctx);
    d = NFLOAT_COMPLEX_IM(y, ctx);

    if (NFLOAT_CTX_HAS_INF_NAN(ctx))
        return _nfloat_complex_mul_naive(r, s, a, b, c, d, ctx);

    if (NFLOAT_IS_ZERO(d))
    {
        status = nfloat_mul(s, b, c, ctx);
        status |= nfloat_mul(r, a, c, ctx);
        return status;
    }

    if (NFLOAT_IS_ZERO(b))
    {
        status = nfloat_mul(s, a, d, ctx);
        status |= nfloat_mul(r, a, c, ctx);
        return status;
    }

    if (NFLOAT_IS_ZERO(c))
    {
        ulong t[NFLOAT_MAX_ALLOC];
        status = nfloat_mul(t, b, d, ctx);
        status |= nfloat_mul(s, a, d, ctx);
        status |= nfloat_neg(r, t, ctx);
        return status;
    }

    if (NFLOAT_IS_ZERO(a))
    {
        ulong t[NFLOAT_MAX_ALLOC];
        status = nfloat_mul(t, b, d, ctx);
        status |= nfloat_mul(s, b, c, ctx);
        status |= nfloat_neg(r, t, ctx);
        return status;
    }

    if (NFLOAT_CTX_NLIMBS(ctx) < 12)
        return _nfloat_complex_mul_standard(r, s, a, b, c, d, ctx);
    else
        return _nfloat_complex_mul_karatsuba(r, s, a, b, c, d, ctx);
}

/*
    Division, inverse, square root and reciprocal square root.

    These share one reciprocal (or reciprocal square root) between the two
    parts of the result. The errors are a few ulp relative to the magnitude
    of the result; the rounding mode is not respected.

    At one limb, the regular cases (no zero or special part, and parts whose
    exponents differ by at most FLINT_BITS / 2 + 2 where |y|^2 or |z| is
    needed) are computed in registers:

    * |y|^2 is the top limb t of S = c^2 + d^2 (both squares at two limbs,
      the smaller one shifted to the scale of the larger), truncated, with
      relative error below 2^(1 - FLINT_BITS); both parts of the
      numerator are then divided by t, the quotients exactly truncated
      (one inverse of t shared, or hardware division where it is fast).

    * |z| = sqrt(S) is the exact floor of the square root of the top two
      limbs of S (_nfloat_isqrt_2); then T = |z| + |a| (no cancellation)
      at one limb, w = 1 / sqrt(2T) (square root) or 1 / sqrt(2 T |z|^2)
      (reciprocal square root) from the one-limb kernel, and the parts are
      T w and |b| w with the signs of the formulas below.
*/

/* (h, l) >>= k for 0 <= k < 2 FLINT_BITS */
#define NFLOAT_RSHIFT_2(h, l, k) \
    do { \
        if ((k) >= FLINT_BITS) \
        { \
            (l) = (h) >> ((k) - FLINT_BITS); \
            (h) = 0; \
        } \
        else if ((k) != 0) \
        { \
            (l) = ((l) >> (k)) | ((h) << (FLINT_BITS - (k))); \
            (h) = (h) >> (k); \
        } \
    } while (0)

/* (m, e) = the truncated product of one-limb mantissas u 2^eu and v 2^ev */
#define NFLOAT_MUL_1(m, e, u, eu, v, ev) \
    do { \
        ulong __h, __l; \
        umul_ppmm(__h, __l, (u), (v)); \
        if (__h >> (FLINT_BITS - 1)) \
        { \
            (e) = (eu) + (ev); \
            (m) = __h; \
        } \
        else \
        { \
            (e) = (eu) + (ev) - 1; \
            (m) = (__h << 1) | (__l >> (FLINT_BITS - 1)); \
        } \
    } while (0)

FLINT_FORCE_INLINE int
_nfloat_set_1_mant(nfloat_ptr res, ulong m, slong exp, int sgnbit, gr_ctx_t ctx)
{
    NFLOAT_CHECK_EXP_RANGE(res, exp, sgnbit, ctx);
    NFLOAT_D(res)[0] = m;
    NFLOAT_EXP(res) = exp;
    NFLOAT_SGNBIT(res) = sgnbit;
    return GR_SUCCESS;
}

/* S = c^2 + d^2 2^-k (or c^2 2^-k + d^2) for nonzero normal c, d at one
   limb, as (s2, s1, s0) with the smaller square shifted by k = 2 |EXP(c) -
   EXP(d)| < 2 FLINT_BITS and truncated; returns e = max(EXP(c), EXP(d)),
   so that c^2 + d^2 = S 2^(2e - 2 FLINT_BITS) (up to the truncation) */
FLINT_FORCE_INLINE slong
_nfloat_complex_sum_sqr_1(ulong * s2, ulong * s1, ulong * s0, nfloat_srcptr c, nfloat_srcptr d)
{
    slong ec = NFLOAT_EXP(c), ed = NFLOAT_EXP(d), e = FLINT_MAX(ec, ed);
    ulong h1, l1, h2, l2;
    slong k;

    umul_ppmm(h1, l1, NFLOAT_D(c)[0], NFLOAT_D(c)[0]);
    umul_ppmm(h2, l2, NFLOAT_D(d)[0], NFLOAT_D(d)[0]);
    k = 2 * (e - ec);
    NFLOAT_RSHIFT_2(h1, l1, k);
    k = 2 * (e - ed);
    NFLOAT_RSHIFT_2(h2, l2, k);
    add_sssaaaaaa(*s2, *s1, *s0, UWORD(0), h1, l1, UWORD(0), h2, l2);
    return e;
}

/* (re, im) = (p (-1)^psgn, q (-1)^qsgn) / |y|^2 with |y|^2 = c^2 + d^2 as
   above, at one limb; p and q are normal or zero, and may alias re, im, c
   or d */
static int
_nfloat_complex_div_sum_sqr_1(nfloat_ptr re, nfloat_ptr im, nfloat_srcptr p, int psgn,
    nfloat_srcptr q, int qsgn, nfloat_srcptr c, nfloat_srcptr d, gr_ctx_t ctx)
{
    ulong s2, s1, s0, t, qp, qq, rp, rq, sp, sq;
    slong e, L, texp, ep, eq;
    int pzero, qzero, status;

    e = _nfloat_complex_sum_sqr_1(&s2, &s1, &s0, c, d);

    /* t = the top limb of S (S >= 2^(2 FLINT_BITS - 2): the larger square
       is not shifted), L = the bit length of S */
    if (s2 != 0)
    {
        t = (UWORD(1) << (FLINT_BITS - 1)) | (s1 >> 1);
        L = 2 * FLINT_BITS + 1;
    }
    else if (s1 >> (FLINT_BITS - 1))
    {
        t = s1;
        L = 2 * FLINT_BITS;
    }
    else
    {
        t = (s1 << 1) | (s0 >> (FLINT_BITS - 1));
        L = 2 * FLINT_BITS - 1;
    }

    /* |y|^2 = (t / B) 2^texp */
    texp = L + 2 * e - 2 * FLINT_BITS;

    pzero = NFLOAT_IS_ZERO(p);
    qzero = NFLOAT_IS_ZERO(q);
    qp = qq = sp = sq = 0;

#if FLINT_PREINVERT_LIMB_USE_NATIVE
    if (!pzero)
        qp = _nfloat_div_1_kernel(&rp, &sp, NFLOAT_D(p)[0], t);
    if (!qzero)
        qq = _nfloat_div_1_kernel(&rq, &sq, NFLOAT_D(q)[0], t);
#else
    {
        ulong tinv;
        NFLOAT_INVERT_LIMB(tinv, t);
        if (!pzero)
            qp = _nfloat_div_1_kernel_preinv(&rp, &sp, NFLOAT_D(p)[0], t, tinv);
        if (!qzero)
            qq = _nfloat_div_1_kernel_preinv(&rq, &sq, NFLOAT_D(q)[0], t, tinv);
    }
#endif
    (void) rp;
    (void) rq;

    ep = pzero ? 0 : NFLOAT_EXP(p) - texp + (slong) sp;
    eq = qzero ? 0 : NFLOAT_EXP(q) - texp + (slong) sq;

    status = pzero ? nfloat_zero(re, ctx) : _nfloat_set_1_mant(re, qp, ep, psgn, ctx);
    status |= qzero ? nfloat_zero(im, ctx) : _nfloat_set_1_mant(im, qq, eq, qsgn, ctx);
    return status;
}

#if FLINT_BITS == 64

/* sqrt(a + bi) (rsqrt = 0) or 1 / sqrt(a + bi) (rsqrt = 1) at one limb, for
   nonzero normal a, b with |EXP(a) - EXP(b)| <= FLINT_BITS / 2 + 2 */
static int
_nfloat_complex_sqrt_rsqrt_1(nfloat_ptr re, nfloat_ptr im, nfloat_srcptr a, nfloat_srcptr b, int rsqrt, gr_ctx_t ctx)
{
    ulong A = NFLOAT_D(a)[0], Bm = NFLOAT_D(b)[0];
    slong ea = NFLOAT_EXP(a), eb = NFLOAT_EXP(b), e, er, eT, eQ, ew, ebig, esmall;
    int asgn = NFLOAT_SGNBIT(a), bsgn = NFLOAT_SGNBIT(b), status;
    ulong s2, s1, s0, n1, n0, rm, a1, a0, c, t1, t0, T, Q, w, big, small, tq;
    unsigned int sh;

    /* |z| = (rm / B) 2^er: the floor of the square root of S, taken
       as (n1, n0) 2^(2j) with n1 >= B/4 */
    e = _nfloat_complex_sum_sqr_1(&s2, &s1, &s0, a, b);
    if (s2 != 0)
    {
        n1 = (s2 << (FLINT_BITS - 2)) | (s1 >> 2);
        n0 = (s1 << (FLINT_BITS - 2)) | (s0 >> 2);
        er = e + 1;
    }
    else
    {
        n1 = s1;
        n0 = s0;
        er = e;
    }
    _nfloat_isqrt_2(&rm, n1, n0);

    /* T = |z| + |a| = (T / B) 2^eT, truncated; |z| >= |a|, with the shift
       er - ea <= FLINT_BITS / 2 + 3 */
    sh = er - ea;
    a1 = A >> sh;
    a0 = (sh == 0) ? 0 : (A << (FLINT_BITS - sh));
    add_sssaaaaaa(c, t1, t0, UWORD(0), rm, UWORD(0), UWORD(0), a1, a0);
    (void) t0;
    if (c)
    {
        T = (UWORD(1) << (FLINT_BITS - 1)) | (t1 >> 1);
        eT = er + 1;
    }
    else
    {
        T = t1;
        eT = er;
    }

    /* w = 1 / sqrt(Q 2^eQ), Q 2^eQ = 2 T resp. 2 T |z|^2 */
    if (!rsqrt)
    {
        Q = T;
        eQ = eT + 1;
    }
    else
    {
        NFLOAT_MUL_1(Q, eQ, T, eT, rm, er);
        NFLOAT_MUL_1(Q, eQ, Q, eQ, rm, er);
        eQ += 1;
    }
    tq = eQ & 1;
    _nfloat_rsqrt_1_approx(&w, Q, tq, 0);
    ew = (1 - (slong) tq) - (eQ - (slong) tq) / 2;

    NFLOAT_MUL_1(big, ebig, T, eT, w, ew);
    NFLOAT_MUL_1(small, esmall, Bm, eb, w, ew);

    /* sqrt: (T w, b w) for a >= 0, (|b| w, sgn(b) T w) for a < 0;
       rsqrt: the conjugates divided by |z| */
    if (!asgn)
    {
        status = _nfloat_set_1_mant(re, big, ebig, 0, ctx);
        status |= _nfloat_set_1_mant(im, small, esmall, bsgn ^ rsqrt, ctx);
    }
    else
    {
        status = _nfloat_set_1_mant(re, small, esmall, 0, ctx);
        status |= _nfloat_set_1_mant(im, big, ebig, bsgn ^ rsqrt, ctx);
    }

    return status;
}

#endif

/* exponents for which squares, products of two and their reciprocals
   cannot overflow or underflow, so that no scaling is needed */
#define NFLOAT_COMPLEX_EXP_MODERATE(e) (FLINT_ABS(e) < NFLOAT_MAX_EXP / 8)

/* x *= 2^e in place, directly on the exponent where that is safe */
FLINT_FORCE_INLINE int
_nfloat_mul_2exp_si_ip(nfloat_ptr x, slong e, gr_ctx_t ctx)
{
    if (!NFLOAT_IS_SPECIAL(x) && NFLOAT_COMPLEX_EXP_MODERATE(NFLOAT_EXP(x)))
    {
        NFLOAT_EXP(x) += e;
        return GR_SUCCESS;
    }
    return nfloat_mul_2exp_si(x, x, e, ctx);
}

/* r = |a + bi| within about one ulp (not correctly rounded): the
   correctly rounded hypot at one and two limbs, where it is fastest, and
   for extreme exponents; otherwise the square root of a^2 + b^2 */
static int
_nfloat_complex_abs_approx(nfloat_ptr r, nfloat_srcptr a, nfloat_srcptr b, gr_ctx_t ctx)
{
    if (NFLOAT_CTX_NLIMBS(ctx) >= 3 &&
        NFLOAT_COMPLEX_EXP_MODERATE(NFLOAT_EXP(a)) && NFLOAT_COMPLEX_EXP_MODERATE(NFLOAT_EXP(b)))
    {
        ulong t[NFLOAT_MAX_ALLOC];
        int status;

        status = nfloat_sqr(r, a, ctx);
        status |= nfloat_sqr(t, b, ctx);
        status |= nfloat_add(r, r, t, ctx);
        status |= nfloat_sqrt(r, r, ctx);
        return status;
    }

    return nfloat_hypot(r, a, b, ctx);
}

/* scales x by 2^-e in place (exactly, barring underflow of a part) */
static int
_nfloat_complex_scale_parts(nfloat_ptr r, nfloat_ptr s, slong e, gr_ctx_t ctx)
{
    int status = GR_SUCCESS;
    if (e != 0)
    {
        status |= nfloat_mul_2exp_si(r, r, -e, ctx);
        status |= nfloat_mul_2exp_si(s, s, -e, ctx);
    }
    return status;
}

int
nfloat_complex_inv(nfloat_complex_ptr res, nfloat_complex_srcptr x, gr_ctx_t ctx)
{
    nfloat_srcptr a, b;
    nfloat_ptr r, s;
    int status;
    slong n = NFLOAT_CTX_NLIMBS(ctx);

    r = NFLOAT_COMPLEX_RE(res, ctx);
    s = NFLOAT_COMPLEX_IM(res, ctx);

    a = NFLOAT_COMPLEX_RE(x, ctx);
    b = NFLOAT_COMPLEX_IM(x, ctx);

    if (NFLOAT_IS_ZERO(b))
    {
        status = nfloat_inv(r, a, ctx);
        nfloat_zero(s, ctx);
        return status;
    }

    if (NFLOAT_IS_ZERO(a))
    {
        status = nfloat_inv(s, b, ctx);
        nfloat_neg(s, s, ctx);
        nfloat_zero(r, ctx);
        return status;
    }

    {
        ulong a2[NFLOAT_MAX_ALLOC];
        ulong b2[NFLOAT_MAX_ALLOC];
        ulong t[NFLOAT_MAX_ALLOC];
        slong ea = NFLOAT_EXP(a), eb = NFLOAT_EXP(b), e, lim;

        /* (a - bi) / (a^2 + b^2), avoiding spurious over- and underflow */
        lim = FLINT_BITS / 2 * n + 2;

        if (ea - eb > lim)
        {
            /* 1/a - (b/a)/a i, the relative error below one ulp */
            status = nfloat_div(t, b, a, ctx);
            status |= nfloat_div(s, t, a, ctx);
            status |= nfloat_inv(r, a, ctx);
            status |= nfloat_neg(s, s, ctx);
            return status;
        }

        if (eb - ea > lim)
        {
            /* (a/b)/b - 1/b i */
            status = nfloat_div(t, a, b, ctx);
            status |= nfloat_div(r, t, b, ctx);
            status |= nfloat_inv(s, b, ctx);
            status |= nfloat_neg(s, s, ctx);
            return status;
        }

        if (n == 1 && !NFLOAT_IS_SPECIAL(a) && !NFLOAT_IS_SPECIAL(b))
            return _nfloat_complex_div_sum_sqr_1(r, s, a, NFLOAT_SGNBIT(a), b, !NFLOAT_SGNBIT(b), a, b, ctx);

        if (NFLOAT_COMPLEX_EXP_MODERATE(ea) && NFLOAT_COMPLEX_EXP_MODERATE(eb))
        {
            /* (a - bi) / (a^2 + b^2) with one reciprocal (r, s may alias
               a, b: each product reads its own part) */
            status = nfloat_sqr(t, a, ctx);
            status |= nfloat_sqr(b2, b, ctx);
            status |= nfloat_add(t, t, b2, ctx);
            status |= nfloat_inv(t, t, ctx);
            status |= nfloat_mul(r, a, t, ctx);
            status |= nfloat_mul(s, b, t, ctx);
            status |= nfloat_neg(s, s, ctx);
            return status;
        }

        /* scale by 2^-e (exact), one reciprocal of a^2 + b^2 */
        e = FLINT_MAX(ea, eb);
        nfloat_set(a2, a, ctx);
        nfloat_set(b2, b, ctx);
        NFLOAT_EXP(a2) -= e;
        NFLOAT_EXP(b2) -= e;

        status = nfloat_sqr(t, a2, ctx);
        status |= nfloat_sqr(s, b2, ctx);
        status |= nfloat_add(t, t, s, ctx);
        status |= nfloat_inv(t, t, ctx);
        status |= nfloat_mul(r, a2, t, ctx);
        status |= nfloat_mul(s, b2, t, ctx);
        status |= nfloat_neg(s, s, ctx);
        status |= _nfloat_complex_scale_parts(r, s, e, ctx);
        return status;
    }
}

int
nfloat_complex_div(nfloat_complex_ptr res, nfloat_complex_srcptr x, nfloat_complex_srcptr y, gr_ctx_t ctx)
{
    nfloat_srcptr a, b, c, d;
    nfloat_ptr r, s;
    int status = GR_SUCCESS;
    slong n = NFLOAT_CTX_NLIMBS(ctx);

    r = NFLOAT_COMPLEX_RE(res, ctx);
    s = NFLOAT_COMPLEX_IM(res, ctx);

    a = NFLOAT_COMPLEX_RE(x, ctx);
    b = NFLOAT_COMPLEX_IM(x, ctx);
    c = NFLOAT_COMPLEX_RE(y, ctx);
    d = NFLOAT_COMPLEX_IM(y, ctx);

    /* todo: other special cases */

    if (NFLOAT_IS_ZERO(d))
    {
        if (NFLOAT_IS_ZERO(b))
        {
            status = nfloat_div(r, a, c, ctx);
            nfloat_zero(s, ctx);
        }
        else if (NFLOAT_IS_ZERO(a))
        {
            status = nfloat_div(s, b, c, ctx);
            nfloat_zero(r, ctx);
        }
        else
        {
            status = nfloat_div(s, b, c, ctx);
            status |= nfloat_div(r, a, c, ctx);
        }
    }
    else if (NFLOAT_IS_ZERO(c))
    {
        if (NFLOAT_IS_ZERO(b))
        {
            status = nfloat_div(s, a, d, ctx);
            nfloat_neg(s, s, ctx);
            nfloat_zero(r, ctx);
        }
        else if (NFLOAT_IS_ZERO(a))
        {
            status = nfloat_div(r, b, d, ctx);
            nfloat_zero(s, ctx);
        }
        else
        {
            status = nfloat_div(r, a, d, ctx);
            status |= nfloat_div(s, b, d, ctx);
            nfloat_swap(r, s, ctx);
            nfloat_neg(s, s, ctx);
        }
    }
    else
    {
        ulong u[2 * NFLOAT_MAX_ALLOC];
        ulong t[NFLOAT_MAX_ALLOC];
        slong ec = NFLOAT_EXP(c), ed = NFLOAT_EXP(d), e;
        slong lim = FLINT_BITS / 2 * n + 2;
        nfloat_ptr ur = NFLOAT_COMPLEX_RE(u, ctx), ui = NFLOAT_COMPLEX_IM(u, ctx);

        if (FLINT_ABS(ec - ed) > lim)
        {
            /* x / y = x (1 / y), with 1 / y as above */
            status = nfloat_complex_inv(u, y, ctx);
            status |= nfloat_complex_mul(res, x, u, ctx);
            return status;
        }

        if (n == 1 && !NFLOAT_IS_SPECIAL(c) && !NFLOAT_IS_SPECIAL(d))
        {
            /* x conj(y) / |y|^2 in registers, x conj(y) unscaled unless it
               overflows */
            int st;

            _nfloat_copy_limbs(ur, c, NFLOAT_CTX_DATA_NLIMBS(ctx));
            _nfloat_copy_limbs(ui, d, NFLOAT_CTX_DATA_NLIMBS(ctx));
            NFLOAT_SGNBIT(ui) ^= 1;

            st = nfloat_complex_mul(u, x, u, ctx);

            if (st == GR_SUCCESS &&
                (!NFLOAT_IS_SPECIAL(ur) || NFLOAT_IS_ZERO(ur)) &&
                (!NFLOAT_IS_SPECIAL(ui) || NFLOAT_IS_ZERO(ui)))
                return _nfloat_complex_div_sum_sqr_1(r, s, ur, NFLOAT_SGNBIT(ur), ui, NFLOAT_SGNBIT(ui), c, d, ctx);
        }

        if (NFLOAT_COMPLEX_EXP_MODERATE(ec) && NFLOAT_COMPLEX_EXP_MODERATE(ed) &&
            NFLOAT_COMPLEX_EXP_MODERATE(NFLOAT_EXP(a)) && NFLOAT_COMPLEX_EXP_MODERATE(NFLOAT_EXP(b)))
        {
            /* x conj(y) / |y|^2 with one reciprocal of |y|^2 */
            _nfloat_copy_limbs(ur, c, NFLOAT_CTX_DATA_NLIMBS(ctx));
            _nfloat_copy_limbs(ui, d, NFLOAT_CTX_DATA_NLIMBS(ctx));
            NFLOAT_SGNBIT(ui) ^= 1;
            status = nfloat_sqr(t, c, ctx);
            status |= nfloat_complex_mul(u, x, u, ctx);
            {
                ulong d2[NFLOAT_MAX_ALLOC];
                status |= nfloat_sqr(d2, d, ctx);
                status |= nfloat_add(t, t, d2, ctx);
            }
            status |= nfloat_inv(t, t, ctx);
            status |= nfloat_mul(r, ur, t, ctx);
            status |= nfloat_mul(s, ui, t, ctx);
            return status;
        }

        /* x conj(y) / |y|^2, with y scaled by 2^-e (exact) and one
           reciprocal of |y|^2 */
        e = FLINT_MAX(ec, ed);
        nfloat_set(ur, c, ctx);
        nfloat_neg(ui, d, ctx);
        NFLOAT_EXP(ur) -= e;
        NFLOAT_EXP(ui) -= e;

        status = nfloat_complex_mul(u, x, u, ctx);

        {
            ulong c2[NFLOAT_MAX_ALLOC];
            ulong d2[NFLOAT_MAX_ALLOC];

            nfloat_set(c2, c, ctx);
            nfloat_set(d2, d, ctx);
            NFLOAT_EXP(c2) -= e;
            NFLOAT_EXP(d2) -= e;
            status |= nfloat_sqr(t, c2, ctx);
            status |= nfloat_sqr(c2, d2, ctx);
            status |= nfloat_add(t, t, c2, ctx);
            status |= nfloat_inv(t, t, ctx);
            status |= nfloat_mul(r, ur, t, ctx);
            status |= nfloat_mul(s, ui, t, ctx);
            status |= _nfloat_complex_scale_parts(r, s, e, ctx);
        }
    }

    return status;
}

int
nfloat_complex_sqrt(nfloat_complex_ptr res, nfloat_complex_srcptr x, gr_ctx_t ctx)
{
    nfloat_srcptr a, b;
    nfloat_ptr c, d;

    ulong r[NFLOAT_MAX_ALLOC];
    ulong t[NFLOAT_MAX_ALLOC];

    int status = GR_SUCCESS;

    c = NFLOAT_COMPLEX_RE(res, ctx);
    d = NFLOAT_COMPLEX_IM(res, ctx);

    a = NFLOAT_COMPLEX_RE(x, ctx);
    b = NFLOAT_COMPLEX_IM(x, ctx);

    if (NFLOAT_IS_ZERO(b))
    {
        if (!NFLOAT_SGNBIT(a))
        {
            status |= nfloat_sqrt(c, a, ctx);
            status |= nfloat_zero(d, ctx);
        }
        else
        {
            status |= nfloat_neg(d, a, ctx);
            status |= nfloat_sqrt(d, d, ctx);
            status |= nfloat_zero(c, ctx);
        }

        return status;
    }

    if (NFLOAT_IS_ZERO(a))
    {
        if (!NFLOAT_SGNBIT(b))
        {
            status |= nfloat_mul_2exp_si(c, b, -1, ctx);
            status |= nfloat_sqrt(c, c, ctx);
            status |= nfloat_set(d, c, ctx);
        }
        else
        {
            status |= nfloat_mul_2exp_si(c, b, -1, ctx);
            status |= nfloat_neg(c, c, ctx);
            status |= nfloat_sqrt(c, c, ctx);
            status |= nfloat_neg(d, c, ctx);
        }

        return status;
    }

#if FLINT_BITS == 64
    if (NFLOAT_CTX_NLIMBS(ctx) == 1 && !NFLOAT_IS_SPECIAL(a) && !NFLOAT_IS_SPECIAL(b) &&
        FLINT_ABS(NFLOAT_EXP(a) - NFLOAT_EXP(b)) <= FLINT_BITS / 2 + 2)
        return _nfloat_complex_sqrt_rsqrt_1(c, d, a, b, 0, ctx);
#endif

    /* r = |z|, q = sqrt((r + |a|) / 2), sqrt(z) = (q, b / (2q)) for a >= 0,
       (|b| / (2q), sgn(b) q) for a < 0 (a single square root and division
       beat a reciprocal square root and two products at all precisions) */
    {
        int asgnbit = NFLOAT_SGNBIT(a), bsgnbit = NFLOAT_SGNBIT(b);

        status |= _nfloat_complex_abs_approx(r, a, b, ctx);
        if (asgnbit)
            status |= nfloat_sub(t, r, a, ctx);
        else
            status |= nfloat_add(t, r, a, ctx);
        status |= _nfloat_mul_2exp_si_ip(t, -1, ctx);

        if (!asgnbit)
        {
            /* c may alias a, which is no longer needed */
            status |= nfloat_sqrt(c, t, ctx);
            status |= nfloat_div(d, b, c, ctx);
            status |= _nfloat_mul_2exp_si_ip(d, -1, ctx);
        }
        else
        {
            status |= nfloat_sqrt(t, t, ctx);
            status |= nfloat_div(c, b, t, ctx);
            NFLOAT_SGNBIT(c) = 0;
            status |= _nfloat_mul_2exp_si_ip(c, -1, ctx);
            status |= nfloat_set(d, t, ctx);
            if (bsgnbit)
                status |= nfloat_neg(d, d, ctx);
        }
    }

    return status;
}

int
nfloat_complex_rsqrt(nfloat_complex_ptr res, nfloat_complex_srcptr x, gr_ctx_t ctx)
{
    nfloat_srcptr a, b;
    nfloat_ptr c, d;

    ulong r[NFLOAT_MAX_ALLOC];
    ulong t[NFLOAT_MAX_ALLOC];
    ulong v[NFLOAT_MAX_ALLOC];

    int status = GR_SUCCESS;

    c = NFLOAT_COMPLEX_RE(res, ctx);
    d = NFLOAT_COMPLEX_IM(res, ctx);

    a = NFLOAT_COMPLEX_RE(x, ctx);
    b = NFLOAT_COMPLEX_IM(x, ctx);

    if (NFLOAT_IS_ZERO(b))
    {
        if (!NFLOAT_SGNBIT(a))
        {
            status |= nfloat_rsqrt(c, a, ctx);
            status |= nfloat_zero(d, ctx);
        }
        else
        {
            status |= nfloat_neg(d, a, ctx);
            status |= nfloat_rsqrt(d, d, ctx);
            status |= nfloat_neg(d, d, ctx);
            status |= nfloat_zero(c, ctx);
        }

        return status;
    }

    if (NFLOAT_IS_ZERO(a))
    {
        if (!NFLOAT_SGNBIT(b))
        {
            status |= nfloat_mul_2exp_si(c, b, 1, ctx);
            status |= nfloat_rsqrt(c, c, ctx);
            status |= nfloat_neg(d, c, ctx);
        }
        else
        {
            status |= nfloat_mul_2exp_si(c, b, 1, ctx);
            status |= nfloat_neg(c, c, ctx);
            status |= nfloat_rsqrt(c, c, ctx);
            status |= nfloat_set(d, c, ctx);
        }

        return status;
    }

#if FLINT_BITS == 64
    if (NFLOAT_CTX_NLIMBS(ctx) == 1 && !NFLOAT_IS_SPECIAL(a) && !NFLOAT_IS_SPECIAL(b) &&
        FLINT_ABS(NFLOAT_EXP(a) - NFLOAT_EXP(b)) <= FLINT_BITS / 2 + 2)
        return _nfloat_complex_sqrt_rsqrt_1(c, d, a, b, 1, ctx);
#endif

    /* T = |z| + |a|, v = 1 / sqrt(2 T |z|^2): 1 / sqrt(z) = (T v, -b v) for
       a >= 0, (|b| v, -sgn(b) T v) for a < 0 */
    {
        int asgnbit = NFLOAT_SGNBIT(a), bsgnbit = NFLOAT_SGNBIT(b);

        if (NFLOAT_CTX_NLIMBS(ctx) >= 3 &&
            NFLOAT_COMPLEX_EXP_MODERATE(NFLOAT_EXP(a)) && NFLOAT_COMPLEX_EXP_MODERATE(NFLOAT_EXP(b)))
        {
            /* as in _nfloat_complex_abs_approx, keeping u = b^2 for
               2 r^2 T = r (b^2 + T^2) (squarings are cheaper than
               products at high precision) */
            ulong u[NFLOAT_MAX_ALLOC];

            status |= nfloat_sqr(u, b, ctx);
            status |= nfloat_sqr(r, a, ctx);
            status |= nfloat_add(r, r, u, ctx);
            status |= nfloat_sqrt(r, r, ctx);
            if (asgnbit)
                status |= nfloat_sub(t, r, a, ctx);
            else
                status |= nfloat_add(t, r, a, ctx);
            status |= nfloat_sqr(v, t, ctx);
            status |= nfloat_add(v, v, u, ctx);
            status |= nfloat_mul(v, v, r, ctx);
        }
        else
        {
            status |= nfloat_hypot(r, a, b, ctx);
            if (asgnbit)
                status |= nfloat_sub(t, r, a, ctx);
            else
                status |= nfloat_add(t, r, a, ctx);
            status |= nfloat_mul(v, t, r, ctx);
            status |= nfloat_mul(v, v, r, ctx);
            status |= _nfloat_mul_2exp_si_ip(v, 1, ctx);
        }
        status |= nfloat_rsqrt(v, v, ctx);

        if (!asgnbit)
        {
            status |= nfloat_mul(d, b, v, ctx);
            status |= nfloat_mul(c, t, v, ctx);
            status |= nfloat_neg(d, d, ctx);
        }
        else
        {
            status |= nfloat_mul(c, b, v, ctx);
            status |= nfloat_mul(d, t, v, ctx);
            NFLOAT_SGNBIT(c) = 0;
            if (!bsgnbit)
                status |= nfloat_neg(d, d, ctx);
        }
    }

    return status;
}

int
nfloat_complex_mul_2exp_si(nfloat_complex_ptr res, nfloat_complex_srcptr x, slong y, gr_ctx_t ctx)
{
    int status = GR_SUCCESS;
    status |= nfloat_mul_2exp_si(NFLOAT_COMPLEX_RE(res, ctx), NFLOAT_COMPLEX_RE(x, ctx), y, ctx);
    status |= nfloat_mul_2exp_si(NFLOAT_COMPLEX_IM(res, ctx), NFLOAT_COMPLEX_IM(x, ctx), y, ctx);
    return status;
}

int
nfloat_complex_cmp(int * res, nfloat_complex_srcptr x, nfloat_complex_srcptr y, gr_ctx_t ctx)
{
    if (NFLOAT_CTX_HAS_INF_NAN(ctx))
        return GR_UNABLE;

    if (!NFLOAT_IS_ZERO(NFLOAT_COMPLEX_IM(x, ctx)) ||
        !NFLOAT_IS_ZERO(NFLOAT_COMPLEX_IM(y, ctx)))
        return GR_DOMAIN;

    return nfloat_cmp(res, NFLOAT_COMPLEX_RE(x, ctx), NFLOAT_COMPLEX_RE(y, ctx), ctx);
}

#include "double_extras.h"

int
nfloat_complex_cmpabs(int * res, nfloat_complex_srcptr x, nfloat_complex_srcptr y, gr_ctx_t ctx)
{
    nfloat_srcptr a, b, c, d;
    slong aexp, bexp, cexp, dexp, xexp, yexp, exp;
    slong xn = NFLOAT_CTX_NLIMBS(ctx);

    if (NFLOAT_CTX_HAS_INF_NAN(ctx))
        return GR_UNABLE;

    a = NFLOAT_COMPLEX_RE(x, ctx);
    b = NFLOAT_COMPLEX_IM(x, ctx);
    c = NFLOAT_COMPLEX_RE(y, ctx);
    d = NFLOAT_COMPLEX_IM(y, ctx);

    if (NFLOAT_IS_ZERO(b))
    {
        if (NFLOAT_IS_ZERO(d))
            return nfloat_cmpabs(res, a, c, ctx);
        if (NFLOAT_IS_ZERO(c))
            return nfloat_cmpabs(res, a, d, ctx);
        if (NFLOAT_IS_ZERO(a))
        {
            *res = -1;
            return GR_SUCCESS;
        }
    }

    if (NFLOAT_IS_ZERO(a))
    {
        if (NFLOAT_IS_ZERO(d))
            return nfloat_cmpabs(res, b, c, ctx);
        if (NFLOAT_IS_ZERO(c))
            return nfloat_cmpabs(res, b, d, ctx);
    }

    if (NFLOAT_IS_ZERO(c) && NFLOAT_IS_ZERO(d))
    {
        *res = 1;
        return GR_SUCCESS;
    }

    aexp = NFLOAT_EXP(a);
    bexp = NFLOAT_EXP(b);
    cexp = NFLOAT_EXP(c);
    dexp = NFLOAT_EXP(d);

    /* 0.5 * 2^xexp <= |x| < sqrt(2) * 2^xexp */
    xexp = FLINT_MAX(aexp, bexp);
    /* 0.5 * 2^yexp <= |y| < sqrt(2) * 2^yexp */
    yexp = FLINT_MAX(cexp, dexp);

    if (xexp + 2 < yexp)
    {
        *res = -1;
        return GR_SUCCESS;
    }

    if (xexp > yexp + 2)
    {
        *res = 1;
        return GR_SUCCESS;
    }

    exp = FLINT_MAX(xexp, yexp);

    double tt, xx = 0.0, yy = 0.0;

    if (aexp >= exp - 53)
    {
        tt = d_mul_2exp_inrange(NFLOAT_D(a)[xn - 1], aexp - exp - FLINT_BITS);
        xx += tt * tt;
    }

    if (bexp >= exp - 53)
    {
        tt = d_mul_2exp_inrange(NFLOAT_D(b)[xn - 1], bexp - exp - FLINT_BITS);
        xx += tt * tt;
    }

    if (cexp >= exp - 53)
    {
        tt = d_mul_2exp_inrange(NFLOAT_D(c)[xn - 1], cexp - exp - FLINT_BITS);
        yy += tt * tt;
    }

    if (dexp >= exp - 53)
    {
        tt = d_mul_2exp_inrange(NFLOAT_D(d)[xn - 1], dexp - exp - FLINT_BITS);
        yy += tt * tt;
    }

    if (xx < yy * 0.999999)
    {
        *res = -1;
        return GR_SUCCESS;
    }

    if (xx * 0.999999 > yy)
    {
        *res = 1;
        return GR_SUCCESS;
    }

    arf_struct s[5];

    arf_init(s + 0);
    arf_init(s + 1);
    arf_init(s + 2);
    arf_init(s + 3);
    arf_init(s + 4);

    nfloat_get_arf(s + 0, a, ctx);
    nfloat_get_arf(s + 1, b, ctx);
    nfloat_get_arf(s + 2, c, ctx);
    nfloat_get_arf(s + 3, d, ctx);

    arf_mul(s + 0, s + 0, s + 0, ARF_PREC_EXACT, ARF_RND_DOWN);
    arf_mul(s + 1, s + 1, s + 1, ARF_PREC_EXACT, ARF_RND_DOWN);
    arf_mul(s + 2, s + 2, s + 2, ARF_PREC_EXACT, ARF_RND_DOWN);
    arf_mul(s + 3, s + 3, s + 3, ARF_PREC_EXACT, ARF_RND_DOWN);
    arf_neg(s + 2, s + 2);
    arf_neg(s + 3, s + 3);
    arf_sum(s + 4, s, 4, 30, ARF_RND_DOWN);
    *res = arf_sgn(s + 4);

    arf_clear(s + 0);
    arf_clear(s + 1);
    arf_clear(s + 2);
    arf_clear(s + 3);
    arf_clear(s + 4);

    return GR_SUCCESS;
}

int
nfloat_complex_abs(nfloat_complex_ptr res, nfloat_complex_srcptr x, gr_ctx_t ctx)
{
    int status;

    if (NFLOAT_CTX_HAS_INF_NAN(ctx))
        return GR_UNABLE;

    /* (correctly rounded, also with directed rounding) */
    status = nfloat_hypot(NFLOAT_COMPLEX_RE(res, ctx), NFLOAT_COMPLEX_RE(x, ctx), NFLOAT_COMPLEX_IM(x, ctx), ctx);
    status |= nfloat_zero(NFLOAT_COMPLEX_IM(res, ctx), ctx);
    return status;
}

void
_nfloat_complex_vec_init(nfloat_complex_ptr res, slong len, gr_ctx_t ctx)
{
    _nfloat_vec_init(res, 2 * len, ctx);
}

void
_nfloat_complex_vec_clear(nfloat_complex_ptr FLINT_UNUSED(res), slong FLINT_UNUSED(len), gr_ctx_t FLINT_UNUSED(ctx))
{
    return;
}

int
_nfloat_complex_vec_zero(nfloat_complex_ptr res, slong len, gr_ctx_t ctx)
{
    return _nfloat_vec_zero(res, 2 * len, ctx);
}

int
_nfloat_complex_vec_set(nfloat_complex_ptr res, nfloat_complex_srcptr x, slong len, gr_ctx_t ctx)
{
    return _nfloat_vec_set(res, x, 2 * len, ctx);
}

int
_nfloat_complex_vec_add(nfloat_complex_ptr res, nfloat_complex_srcptr x, nfloat_complex_srcptr y, slong len, gr_ctx_t ctx)
{
    return _nfloat_vec_add(res, x, y, 2 * len, ctx);
}

int
_nfloat_complex_vec_sub(nfloat_complex_ptr res, nfloat_complex_srcptr x, nfloat_complex_srcptr y, slong len, gr_ctx_t ctx)
{
    return _nfloat_vec_sub(res, x, y, 2 * len, ctx);
}


int _nfloat_complex_methods_initialized = 0;

gr_static_method_table _nfloat_complex_methods;

gr_method_tab_input _nfloat_complex_methods_input[] =
{
    {GR_METHOD_CTX_WRITE,       (gr_funcptr) nfloat_ctx_write},
    {GR_METHOD_CTX_IS_RING,     (gr_funcptr) gr_generic_ctx_predicate_false},
    {GR_METHOD_CTX_IS_COMMUTATIVE_RING, (gr_funcptr) gr_generic_ctx_predicate_false},
    {GR_METHOD_CTX_IS_INTEGRAL_DOMAIN,  (gr_funcptr) gr_generic_ctx_predicate_false},
    {GR_METHOD_CTX_IS_FIELD,            (gr_funcptr) gr_generic_ctx_predicate_false},
    {GR_METHOD_CTX_IS_UNIQUE_FACTORIZATION_DOMAIN,
                                (gr_funcptr) gr_generic_ctx_predicate_false},
    {GR_METHOD_CTX_IS_FINITE,
                                (gr_funcptr) gr_generic_ctx_predicate_false},
    {GR_METHOD_CTX_IS_FINITE_CHARACTERISTIC,
                                (gr_funcptr) gr_generic_ctx_predicate_false},
    {GR_METHOD_CTX_IS_ALGEBRAICALLY_CLOSED,
                                (gr_funcptr) gr_generic_ctx_predicate_false},
    {GR_METHOD_CTX_IS_ORDERED_RING,
                                (gr_funcptr) gr_generic_ctx_predicate_false},
    {GR_METHOD_CTX_IS_APPROX_COMMUTATIVE_RING,     (gr_funcptr) gr_generic_ctx_predicate_true},
    {GR_METHOD_CTX_IS_EXACT,    (gr_funcptr) gr_generic_ctx_predicate_false},
    {GR_METHOD_CTX_IS_CANONICAL,
                                (gr_funcptr) gr_generic_ctx_predicate_false},

    {GR_METHOD_CTX_HAS_REAL_PREC, (gr_funcptr) gr_generic_ctx_predicate_true},
    {GR_METHOD_CTX_SET_REAL_PREC, (gr_funcptr) nfloat_ctx_set_real_prec},
    {GR_METHOD_CTX_GET_REAL_PREC, (gr_funcptr) nfloat_ctx_get_real_prec},

    {GR_METHOD_INIT,            (gr_funcptr) nfloat_complex_init},
    {GR_METHOD_CLEAR,           (gr_funcptr) nfloat_complex_clear},
    {GR_METHOD_SWAP,            (gr_funcptr) nfloat_complex_swap},
    {GR_METHOD_SET_SHALLOW,     (gr_funcptr) nfloat_complex_set},
    {GR_METHOD_RANDTEST,        (gr_funcptr) nfloat_complex_randtest},
    {GR_METHOD_WRITE,           (gr_funcptr) nfloat_complex_write},
    {GR_METHOD_ZERO,            (gr_funcptr) nfloat_complex_zero},
    {GR_METHOD_ONE,             (gr_funcptr) nfloat_complex_one},
    {GR_METHOD_NEG_ONE,         (gr_funcptr) nfloat_complex_neg_one},
    {GR_METHOD_IS_ZERO,         (gr_funcptr) nfloat_complex_is_zero},
    {GR_METHOD_IS_ONE,          (gr_funcptr) nfloat_complex_is_one},
    {GR_METHOD_IS_NEG_ONE,      (gr_funcptr) nfloat_complex_is_neg_one},
    {GR_METHOD_EQUAL,           (gr_funcptr) nfloat_complex_equal},
    {GR_METHOD_SET,             (gr_funcptr) nfloat_complex_set},
    {GR_METHOD_SET_SI,          (gr_funcptr) nfloat_complex_set_si},
    {GR_METHOD_SET_UI,          (gr_funcptr) nfloat_complex_set_ui},
    {GR_METHOD_SET_FMPZ,        (gr_funcptr) nfloat_complex_set_fmpz},
    {GR_METHOD_SET_FMPQ,        (gr_funcptr) nfloat_complex_set_fmpq},
    {GR_METHOD_SET_D,           (gr_funcptr) nfloat_complex_set_d},
    {GR_METHOD_SET_STR,         (gr_funcptr) gr_generic_set_str_ring_exponents},
    {GR_METHOD_SET_OTHER,       (gr_funcptr) nfloat_complex_set_other},
/*
    {GR_METHOD_GET_FMPZ,        (gr_funcptr) nfloat_complex_get_fmpz},
    {GR_METHOD_GET_FMPQ,        (gr_funcptr) nfloat_complex_get_fmpq},
    {GR_METHOD_GET_UI,          (gr_funcptr) nfloat_complex_get_ui},
    {GR_METHOD_GET_SI,          (gr_funcptr) nfloat_complex_get_si},
    {GR_METHOD_GET_D,           (gr_funcptr) nfloat_complex_get_d},
*/

    {GR_METHOD_NEG,             (gr_funcptr) nfloat_complex_neg},
    {GR_METHOD_ADD,             (gr_funcptr) nfloat_complex_add},
/*
    {GR_METHOD_ADD_UI,          (gr_funcptr) nfloat_complex_add_ui},
    {GR_METHOD_ADD_SI,          (gr_funcptr) nfloat_complex_add_si},
    {GR_METHOD_ADD_FMPZ,        (gr_funcptr) nfloat_complex_add_fmpz},
*/
    {GR_METHOD_SUB,             (gr_funcptr) nfloat_complex_sub},
/*
    {GR_METHOD_SUB_UI,          (gr_funcptr) nfloat_complex_sub_ui},
    {GR_METHOD_SUB_SI,          (gr_funcptr) nfloat_complex_sub_si},
    {GR_METHOD_SUB_FMPZ,        (gr_funcptr) nfloat_complex_sub_fmpz},
*/
    {GR_METHOD_MUL,             (gr_funcptr) nfloat_complex_mul},
/*
    {GR_METHOD_MUL_UI,          (gr_funcptr) nfloat_complex_mul_ui},
    {GR_METHOD_MUL_SI,          (gr_funcptr) nfloat_complex_mul_si},
    {GR_METHOD_MUL_FMPZ,        (gr_funcptr) nfloat_complex_mul_fmpz},
    {GR_METHOD_MUL_TWO,         (gr_funcptr) nfloat_complex_mul_two},
*/
/*
    {GR_METHOD_ADDMUL,          (gr_funcptr) nfloat_complex_addmul},
    {GR_METHOD_SUBMUL,          (gr_funcptr) nfloat_complex_submul},
*/
    {GR_METHOD_SQR,             (gr_funcptr) nfloat_complex_sqr},
    {GR_METHOD_DIV,             (gr_funcptr) nfloat_complex_div},
/*
    {GR_METHOD_DIV_UI,          (gr_funcptr) nfloat_complex_div_ui},
    {GR_METHOD_DIV_SI,          (gr_funcptr) nfloat_complex_div_si},
*/
/*
    {GR_METHOD_DIV_FMPZ,        (gr_funcptr) nfloat_complex_div_fmpz},
*/
    {GR_METHOD_INV,             (gr_funcptr) nfloat_complex_inv},
    {GR_METHOD_MUL_2EXP_SI,     (gr_funcptr) nfloat_complex_mul_2exp_si},
/*
    {GR_METHOD_MUL_2EXP_FMPZ,      (gr_funcptr) nfloat_complex_mul_2exp_fmpz},
    {GR_METHOD_SET_FMPZ_2EXP_FMPZ, (gr_funcptr) nfloat_complex_set_fmpz_2exp_fmpz},
    {GR_METHOD_GET_FMPZ_2EXP_FMPZ, (gr_funcptr) nfloat_complex_get_fmpz_2exp_fmpz},
*/

    {GR_METHOD_POW,             (gr_funcptr) nfloat_complex_pow},
/*
    {GR_METHOD_POW_UI,          (gr_funcptr) nfloat_complex_pow_ui},
    {GR_METHOD_POW_SI,          (gr_funcptr) nfloat_complex_pow_si},
    {GR_METHOD_POW_FMPZ,        (gr_funcptr) nfloat_complex_pow_fmpz},
    {GR_METHOD_POW_FMPQ,        (gr_funcptr) nfloat_complex_pow_fmpq},
*/
    {GR_METHOD_SQRT,            (gr_funcptr) nfloat_complex_sqrt},
    {GR_METHOD_RSQRT,           (gr_funcptr) nfloat_complex_rsqrt},
/*
    {GR_METHOD_POS_INF,         (gr_funcptr) nfloat_complex_pos_inf},
    {GR_METHOD_NEG_INF,         (gr_funcptr) nfloat_complex_neg_inf},
    {GR_METHOD_UINF,            (gr_funcptr) gr_not_in_domain},
    {GR_METHOD_UNDEFINED,       (gr_funcptr) nfloat_complex_nan},
    {GR_METHOD_UNKNOWN,         (gr_funcptr) nfloat_complex_nan},

    {GR_METHOD_FLOOR,           (gr_funcptr) nfloat_complex_floor},
    {GR_METHOD_CEIL,            (gr_funcptr) nfloat_complex_ceil},
    {GR_METHOD_TRUNC,           (gr_funcptr) nfloat_complex_trunc},
    {GR_METHOD_NINT,            (gr_funcptr) nfloat_complex_nint},
*/

    {GR_METHOD_ABS,             (gr_funcptr) nfloat_complex_abs},
    {GR_METHOD_CONJ,            (gr_funcptr) nfloat_complex_conj},
    {GR_METHOD_RE,              (gr_funcptr) nfloat_complex_re},
    {GR_METHOD_IM,              (gr_funcptr) nfloat_complex_im},
/*
    {GR_METHOD_SGN,             (gr_funcptr) nfloat_complex_sgn},
    {GR_METHOD_CSGN,            (gr_funcptr) nfloat_complex_sgn},
*/
    {GR_METHOD_CMP,             (gr_funcptr) nfloat_complex_cmp},
    {GR_METHOD_CMPABS,          (gr_funcptr) nfloat_complex_cmpabs},

    {GR_METHOD_I,               (gr_funcptr) nfloat_complex_i},
    {GR_METHOD_PI,              (gr_funcptr) nfloat_complex_pi},
    {GR_METHOD_EXP,             (gr_funcptr) nfloat_complex_exp},
    {GR_METHOD_EXPM1,           (gr_funcptr) nfloat_complex_expm1},
    {GR_METHOD_EXP_PI_I,        (gr_funcptr) nfloat_complex_exp_pi_i},
    {GR_METHOD_EXP2,            (gr_funcptr) nfloat_complex_exp2},
    {GR_METHOD_EXP10,           (gr_funcptr) nfloat_complex_exp10},
    {GR_METHOD_LOG,             (gr_funcptr) nfloat_complex_log},
    {GR_METHOD_LOG1P,           (gr_funcptr) nfloat_complex_log1p},
    {GR_METHOD_LOG_PI_I,        (gr_funcptr) nfloat_complex_log_pi_i},
    {GR_METHOD_LOG2,            (gr_funcptr) nfloat_complex_log2},
    {GR_METHOD_LOG10,           (gr_funcptr) nfloat_complex_log10},
    {GR_METHOD_SIN,             (gr_funcptr) nfloat_complex_sin},
    {GR_METHOD_COS,             (gr_funcptr) nfloat_complex_cos},
    {GR_METHOD_TAN,             (gr_funcptr) nfloat_complex_tan},
    {GR_METHOD_COT,             (gr_funcptr) nfloat_complex_cot},
    {GR_METHOD_SEC,             (gr_funcptr) nfloat_complex_sec},
    {GR_METHOD_CSC,             (gr_funcptr) nfloat_complex_csc},
    {GR_METHOD_SIN_PI,          (gr_funcptr) nfloat_complex_sin_pi},
    {GR_METHOD_COS_PI,          (gr_funcptr) nfloat_complex_cos_pi},
    {GR_METHOD_TAN_PI,          (gr_funcptr) nfloat_complex_tan_pi},
    {GR_METHOD_COT_PI,          (gr_funcptr) nfloat_complex_cot_pi},
    {GR_METHOD_SEC_PI,          (gr_funcptr) nfloat_complex_sec_pi},
    {GR_METHOD_CSC_PI,          (gr_funcptr) nfloat_complex_csc_pi},
    {GR_METHOD_SINC,            (gr_funcptr) nfloat_complex_sinc},
    {GR_METHOD_SINC_PI,         (gr_funcptr) nfloat_complex_sinc_pi},
    {GR_METHOD_SINH,            (gr_funcptr) nfloat_complex_sinh},
    {GR_METHOD_COSH,            (gr_funcptr) nfloat_complex_cosh},
    {GR_METHOD_TANH,            (gr_funcptr) nfloat_complex_tanh},
    {GR_METHOD_COTH,            (gr_funcptr) nfloat_complex_coth},
    {GR_METHOD_SECH,            (gr_funcptr) nfloat_complex_sech},
    {GR_METHOD_CSCH,            (gr_funcptr) nfloat_complex_csch},
    {GR_METHOD_ASIN,            (gr_funcptr) nfloat_complex_asin},
    {GR_METHOD_ACOS,            (gr_funcptr) nfloat_complex_acos},
    {GR_METHOD_ATAN,            (gr_funcptr) nfloat_complex_atan},
    {GR_METHOD_ACOT,            (gr_funcptr) nfloat_complex_acot},
    {GR_METHOD_ASEC,            (gr_funcptr) nfloat_complex_asec},
    {GR_METHOD_ACSC,            (gr_funcptr) nfloat_complex_acsc},
    {GR_METHOD_ASINH,           (gr_funcptr) nfloat_complex_asinh},
    {GR_METHOD_ACOSH,           (gr_funcptr) nfloat_complex_acosh},
    {GR_METHOD_ATANH,           (gr_funcptr) nfloat_complex_atanh},
    {GR_METHOD_ACOTH,           (gr_funcptr) nfloat_complex_acoth},
    {GR_METHOD_ASECH,           (gr_funcptr) nfloat_complex_asech},
    {GR_METHOD_ACSCH,           (gr_funcptr) nfloat_complex_acsch},
    {GR_METHOD_ASIN_PI,         (gr_funcptr) nfloat_complex_asin_pi},
    {GR_METHOD_ACOS_PI,         (gr_funcptr) nfloat_complex_acos_pi},
    {GR_METHOD_ATAN_PI,         (gr_funcptr) nfloat_complex_atan_pi},
    {GR_METHOD_ACOT_PI,         (gr_funcptr) nfloat_complex_acot_pi},
    {GR_METHOD_ASEC_PI,         (gr_funcptr) nfloat_complex_asec_pi},
    {GR_METHOD_ACSC_PI,         (gr_funcptr) nfloat_complex_acsc_pi},
    {GR_METHOD_SIN_COS,         (gr_funcptr) nfloat_complex_sin_cos},
    {GR_METHOD_SIN_COS_PI,      (gr_funcptr) nfloat_complex_sin_cos_pi},
    {GR_METHOD_SINH_COSH,       (gr_funcptr) nfloat_complex_sinh_cosh},
/*
    {GR_METHOD_GAMMA,            (gr_funcptr) nfloat_complex_gamma},
    {GR_METHOD_ZETA,             (gr_funcptr) nfloat_complex_zeta},
*/

    {GR_METHOD_VEC_INIT,        (gr_funcptr) _nfloat_complex_vec_init},
    {GR_METHOD_VEC_CLEAR,       (gr_funcptr) _nfloat_complex_vec_clear},
    {GR_METHOD_VEC_SET,         (gr_funcptr) _nfloat_complex_vec_set},
    {GR_METHOD_VEC_ZERO,        (gr_funcptr) _nfloat_complex_vec_zero},
    {GR_METHOD_VEC_ADD,                 (gr_funcptr) _nfloat_complex_vec_add},
    {GR_METHOD_VEC_SUB,                 (gr_funcptr) _nfloat_complex_vec_sub},
/*
    {GR_METHOD_VEC_MUL,                 (gr_funcptr) _nfloat_complex_vec_mul},
    {GR_METHOD_VEC_MUL_SCALAR,          (gr_funcptr) _nfloat_complex_vec_mul_scalar},
    {GR_METHOD_VEC_ADDMUL_SCALAR,          (gr_funcptr) _nfloat_complex_vec_addmul_scalar},
    {GR_METHOD_VEC_SUBMUL_SCALAR,          (gr_funcptr) _nfloat_complex_vec_submul_scalar},
*/
    {GR_METHOD_VEC_DOT,         (gr_funcptr) _nfloat_complex_vec_dot},
    {GR_METHOD_VEC_DOT_REV,     (gr_funcptr) _nfloat_complex_vec_dot_rev},
/*
    {GR_METHOD_POLY_MULLOW,     (gr_funcptr) nfloat_complex_poly_mullow},
    {GR_METHOD_POLY_ROOTS_OTHER,(gr_funcptr) nfloat_complex_poly_roots_other},
*/
    {GR_METHOD_MAT_MUL,         (gr_funcptr) nfloat_complex_mat_mul},
    {GR_METHOD_MAT_NONSINGULAR_SOLVE_TRIL,  (gr_funcptr) nfloat_complex_mat_nonsingular_solve_tril},
    {GR_METHOD_MAT_NONSINGULAR_SOLVE_TRIU,  (gr_funcptr) nfloat_complex_mat_nonsingular_solve_triu},
    {GR_METHOD_MAT_LU,                      (gr_funcptr) nfloat_complex_mat_lu},
    {GR_METHOD_MAT_DET,         (gr_funcptr) gr_mat_det_generic_field},
    {GR_METHOD_MAT_FIND_NONZERO_PIVOT,     (gr_funcptr) gr_mat_find_nonzero_pivot_large_abs},

    {0,                         (gr_funcptr) NULL},
};

int
nfloat_complex_ctx_init(gr_ctx_t ctx, slong prec, int flags)
{
    slong nlimbs;

    if (prec <= 0 || prec > NFLOAT_MAX_LIMBS * FLINT_BITS)
        return GR_UNABLE;

    nlimbs = (prec + FLINT_BITS - 1) / FLINT_BITS;

    ctx->which_ring = GR_CTX_NFLOAT_COMPLEX;
    ctx->sizeof_elem = 2 * sizeof(ulong) * (nlimbs + NFLOAT_HEADER_LIMBS);
    ctx->size_limit = WORD_MAX;

    NFLOAT_CTX_NLIMBS(ctx) = nlimbs;
    NFLOAT_CTX_FLAGS(ctx) = flags;
    NFLOAT_CTX_RND(ctx) = 0;
    NFLOAT_CTX_FUNC_PREC(ctx) = nlimbs * FLINT_BITS;

    ctx->methods = _nfloat_complex_methods;

    if (!_nfloat_complex_methods_initialized)
    {
        gr_method_tab_init(_nfloat_complex_methods, _nfloat_complex_methods_input);
        _nfloat_complex_methods_initialized = 1;
    }

    return GR_SUCCESS;
}
