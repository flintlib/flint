/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/*
    Complex decimal floating-point numbers: pairs of decfloats with
    componentwise correct rounding. The real part is rounded in the
    context rounding mode and the imaginary part in the context rounding
    mode for imaginary parts.
*/

#include <string.h>
#include "qqbar.h"
#include "decimal.h"
#include "fmpq.h"
#include "arf.h"
#include "arb.h"
#include "acf.h"
#include "acb.h"
#include "gr.h"
#include "gr_generic.h"
#include "gr_vec.h"
#include "gr_mat.h"

#define RND(ctx) DECIMAL_CTX_RND(ctx)
#define RND_IM(ctx) DECIMAL_CTX_RND_IM(ctx)
#define PREC(ctx) DECIMAL_CTX_PREC(ctx)
#define EXACT_RND (DECIMAL_RND_DOWN | DECIMAL_RND_NOLIMITS)

/* status for real-valued conversions of a nonreal number */
static int
_nonreal_status(const deccfloat_t x)
{
    return DECFLOAT_IS_NAN(&x->im) ? GR_UNABLE : GR_DOMAIN;
}

/* ------------------------------------------------------------------------- */
/*    Exact helpers                                                          */
/* ------------------------------------------------------------------------- */

/* res = x1 y1 +/- x2 y2 with a single rounding. A zero factor annihilates
   its term even if the other factor is infinite, so that zero components
   act as exact zeros (i * inf = inf * i). */
int
_decfloat_dot2(decfloat_t res, const decfloat_t x1, const decfloat_t y1, const decfloat_t x2, const decfloat_t y2, int subtract, slong prec, int rnd, decimal_rounding_info * info, decmag_ptr err, gr_ctx_t ctx)
{
    decfloat_t t1, t2;
    int status, z1, z2;

    z1 = DECFLOAT_IS_ZERO(x1) || DECFLOAT_IS_ZERO(y1);
    z2 = DECFLOAT_IS_ZERO(x2) || DECFLOAT_IS_ZERO(y2);

    if (z1 && z2)
    {
        if (info != NULL)
        {
            info->inexact = 0;
            info->increased = 0;
            info->underflow = 0;
            info->overflow = 0;
        }
        if (err != NULL)
            _decmag_zero(err, ctx);
        return decfloat_zero(res, ctx);
    }

    if (z2)
        return _decfloat_mul(res, x1, y1, prec, rnd, info, err, ctx);

    if (z1)
    {
        if (!subtract)
            return _decfloat_mul(res, x2, y2, prec, rnd, info, err, ctx);

        status = _decfloat_mul(res, x2, y2, prec, DECIMAL_RND_NEGATE(rnd), info, err, ctx);
        if (status == GR_SUCCESS)
        {
            if (DECFLOAT_IS_POS_INF(res)) _decfloat_neg_inf(res);
            else if (DECFLOAT_IS_NEG_INF(res)) _decfloat_pos_inf(res);
            else res->m.size = -res->m.size;
        }
        return status;
    }

    decfloat_init(t1, ctx);
    decfloat_init(t2, ctx);

    status = _decfloat_mul(t1, x1, y1, DECIMAL_PREC_EXACT, EXACT_RND, NULL, NULL, ctx);
    status |= _decfloat_mul(t2, x2, y2, DECIMAL_PREC_EXACT, EXACT_RND, NULL, NULL, ctx);

    if (status == GR_SUCCESS)
        status = _decfloat_add(res, t1, t2, subtract, prec, rnd, info, err, ctx);

    decfloat_clear(t1, ctx);
    decfloat_clear(t2, ctx);
    return status;
}

/* res = x / 2 exactly */
static int
_decfloat_half_exact(decfloat_t res, const decfloat_t x, gr_ctx_t ctx)
{
    decfloat_t five;
    int status;
    decfloat_init(five, ctx);
    GR_MUST_SUCCEED(decfloat_set_round_si(five, 5, DECIMAL_PREC_EXACT, EXACT_RND, ctx));
    status = _decfloat_mul(res, x, five, DECIMAL_PREC_EXACT, EXACT_RND, NULL, NULL, ctx);
    if (status == GR_SUCCESS)
        status = decfloat_mul_10exp_si_round(res, res, -1, DECIMAL_PREC_EXACT, EXACT_RND, ctx);
    decfloat_clear(five, ctx);
    return status;
}

/* exact arithmetic; GR_UNABLE for special values */
int
_deccfloat_mul_exact(deccfloat_t res, const deccfloat_t x, const deccfloat_t y, gr_ctx_t ctx)
{
    decfloat_t re, im;
    int status;

    if (!_deccfloat_is_finite(x) || !_deccfloat_is_finite(y))
        return GR_UNABLE;

    decfloat_init(re, ctx);
    decfloat_init(im, ctx);
    status = _decfloat_dot2(re, &x->re, &y->re, &x->im, &y->im, 1, DECIMAL_PREC_EXACT, EXACT_RND, NULL, NULL, ctx);
    status |= _decfloat_dot2(im, &x->re, &y->im, &x->im, &y->re, 0, DECIMAL_PREC_EXACT, EXACT_RND, NULL, NULL, ctx);
    decfloat_swap(&res->re, re, ctx);
    decfloat_swap(&res->im, im, ctx);
    decfloat_clear(re, ctx);
    decfloat_clear(im, ctx);
    return status;
}

int
_deccfloat_add_exact(deccfloat_t res, const deccfloat_t x, const deccfloat_t y, int subtract, gr_ctx_t ctx)
{
    int status;

    if (!_deccfloat_is_finite(x) || !_deccfloat_is_finite(y))
        return GR_UNABLE;

    status = _decfloat_add(&res->re, &x->re, &y->re, subtract, DECIMAL_PREC_EXACT, EXACT_RND, NULL, NULL, ctx);
    status |= _decfloat_add(&res->im, &x->im, &y->im, subtract, DECIMAL_PREC_EXACT, EXACT_RND, NULL, NULL, ctx);
    return status;
}

/* 1/x exactly if it is a Gaussian decimal fraction */
int
_deccfloat_inv_exact(deccfloat_t res, const deccfloat_t x, gr_ctx_t ctx)
{
    decfloat_t den, re, im;
    int status;

    if (!_deccfloat_is_finite(x) || (DECFLOAT_IS_ZERO(&x->re) && DECFLOAT_IS_ZERO(&x->im)))
        return GR_UNABLE;

    decfloat_init(den, ctx);
    decfloat_init(re, ctx);
    decfloat_init(im, ctx);

    status = _decfloat_dot2(den, &x->re, &x->re, &x->im, &x->im, 0, DECIMAL_PREC_EXACT, EXACT_RND, NULL, NULL, ctx);
    if (status == GR_SUCCESS)
        status = _decfloat_div(re, &x->re, den, DECIMAL_PREC_EXACT, EXACT_RND, NULL, NULL, ctx);
    if (status == GR_SUCCESS)
        status = _decfloat_div(im, &x->im, den, DECIMAL_PREC_EXACT, EXACT_RND, NULL, NULL, ctx);

    if (status == GR_SUCCESS)
    {
        im->m.size = -im->m.size;
        decfloat_swap(&res->re, re, ctx);
        decfloat_swap(&res->im, im, ctx);
    }

    decfloat_clear(den, ctx);
    decfloat_clear(re, ctx);
    decfloat_clear(im, ctx);
    return status;
}

int
_deccfloat_pow_ui_exact(deccfloat_t res, const deccfloat_t x, ulong n, gr_ctx_t ctx)
{
    deccfloat_t y, t;
    int status = GR_SUCCESS;

    if (n == 0)
        return deccfloat_one(res, ctx);

    if (!_deccfloat_is_finite(x))
        return GR_UNABLE;

    if (n == 1)
    {
        status = decfloat_set_round(&res->re, &x->re, DECIMAL_PREC_EXACT, EXACT_RND, ctx);
        status |= decfloat_set_round(&res->im, &x->im, DECIMAL_PREC_EXACT, EXACT_RND, ctx);
        return status;
    }

    deccfloat_init(y, ctx);
    deccfloat_init(t, ctx);

    /* binary powering from the top bit */
    {
        int bit = FLINT_BIT_COUNT(n) - 2;
        status = _deccfloat_pow_ui_exact(y, x, 1, ctx);

        for ( ; bit >= 0 && status == GR_SUCCESS; bit--)
        {
            status = _deccfloat_mul_exact(y, y, y, ctx);
            if (status == GR_SUCCESS && ((n >> bit) & 1))
                status = _deccfloat_mul_exact(y, y, x, ctx);
        }
    }

    if (status == GR_SUCCESS)
        deccfloat_swap(res, y, ctx);

    deccfloat_clear(y, ctx);
    deccfloat_clear(t, ctx);
    return status;
}

/* principal square root, exactly if it is a Gaussian decimal fraction;
   requires a nonreal argument */
int
_deccfloat_sqrt_exact(deccfloat_t res, const deccfloat_t x, gr_ctx_t ctx)
{
    decfloat_t r, p, q;
    int status;

    if (!_deccfloat_is_finite(x) || DECFLOAT_IS_ZERO(&x->im))
        return GR_UNABLE;

    decfloat_init(r, ctx);
    decfloat_init(p, ctx);
    decfloat_init(q, ctx);

    /* r = |x| */
    status = _decfloat_dot2(r, &x->re, &x->re, &x->im, &x->im, 0, DECIMAL_PREC_EXACT, EXACT_RND, NULL, NULL, ctx);
    if (status == GR_SUCCESS)
        status = _decfloat_sqrt(r, r, DECIMAL_PREC_EXACT, EXACT_RND, NULL, NULL, ctx);
    /* p = sqrt((r + a) / 2) */
    if (status == GR_SUCCESS)
        status = _decfloat_add(p, r, &x->re, 0, DECIMAL_PREC_EXACT, EXACT_RND, NULL, NULL, ctx);
    if (status == GR_SUCCESS)
        status = _decfloat_half_exact(p, p, ctx);
    if (status == GR_SUCCESS)
        status = _decfloat_sqrt(p, p, DECIMAL_PREC_EXACT, EXACT_RND, NULL, NULL, ctx);
    /* q = b / (2 p) */
    if (status == GR_SUCCESS)
        status = _decfloat_add(q, p, p, 0, DECIMAL_PREC_EXACT, EXACT_RND, NULL, NULL, ctx);
    if (status == GR_SUCCESS)
        status = _decfloat_div(q, &x->im, q, DECIMAL_PREC_EXACT, EXACT_RND, NULL, NULL, ctx);

    if (status == GR_SUCCESS)
    {
        decfloat_swap(&res->re, p, ctx);
        decfloat_swap(&res->im, q, ctx);
    }

    decfloat_clear(r, ctx);
    decfloat_clear(p, ctx);
    decfloat_clear(q, ctx);
    return status;
}

int
_deccfloat_finalize(deccfloat_t res, gr_ctx_t ctx)
{
    int status = _decfloat_finalize(&res->re, ctx);
    status |= _decfloat_finalize(&res->im, ctx);
    return status;
}

/* ------------------------------------------------------------------------- */
/*    Arithmetic                                                             */
/* ------------------------------------------------------------------------- */

int
deccfloat_neg(deccfloat_t res, const deccfloat_t x, gr_ctx_t ctx)
{
    int status;
    status = decfloat_neg_round(&res->re, &x->re, PREC(ctx), RND(ctx), ctx);
    status |= decfloat_neg_round(&res->im, &x->im, PREC(ctx), RND_IM(ctx), ctx);
    return status;
}

int
deccfloat_conj(deccfloat_t res, const deccfloat_t x, gr_ctx_t ctx)
{
    int status;
    status = decfloat_set_round(&res->re, &x->re, PREC(ctx), RND(ctx), ctx);
    status |= decfloat_neg_round(&res->im, &x->im, PREC(ctx), RND_IM(ctx), ctx);
    return status;
}

int
deccfloat_re(deccfloat_t res, const deccfloat_t x, gr_ctx_t ctx)
{
    int status = decfloat_set_round(&res->re, &x->re, PREC(ctx), RND(ctx), ctx);
    decfloat_zero(&res->im, ctx);
    return status;
}

int
deccfloat_im(deccfloat_t res, const deccfloat_t x, gr_ctx_t ctx)
{
    int status = decfloat_set_round(&res->re, &x->im, PREC(ctx), RND(ctx), ctx);
    decfloat_zero(&res->im, ctx);
    return status;
}

int
deccfloat_add(deccfloat_t res, const deccfloat_t x, const deccfloat_t y, gr_ctx_t ctx)
{
    int status;
    status = decfloat_add_round(&res->re, &x->re, &y->re, PREC(ctx), RND(ctx), ctx);
    status |= decfloat_add_round(&res->im, &x->im, &y->im, PREC(ctx), RND_IM(ctx), ctx);
    return status;
}

int
deccfloat_sub(deccfloat_t res, const deccfloat_t x, const deccfloat_t y, gr_ctx_t ctx)
{
    int status;
    status = decfloat_sub_round(&res->re, &x->re, &y->re, PREC(ctx), RND(ctx), ctx);
    status |= decfloat_sub_round(&res->im, &x->im, &y->im, PREC(ctx), RND_IM(ctx), ctx);
    return status;
}

int
deccfloat_mul(deccfloat_t res, const deccfloat_t x, const deccfloat_t y, gr_ctx_t ctx)
{
    decfloat_t re, im;
    int status;

    decfloat_init(re, ctx);
    decfloat_init(im, ctx);

    status = _decfloat_dot2(re, &x->re, &y->re, &x->im, &y->im, 1, PREC(ctx), RND(ctx), NULL, NULL, ctx);
    status |= _decfloat_dot2(im, &x->re, &y->im, &x->im, &y->re, 0, PREC(ctx), RND_IM(ctx), NULL, NULL, ctx);

    decfloat_swap(&res->re, re, ctx);
    decfloat_swap(&res->im, im, ctx);
    decfloat_clear(re, ctx);
    decfloat_clear(im, ctx);
    return status;
}

int
deccfloat_sqr(deccfloat_t res, const deccfloat_t x, gr_ctx_t ctx)
{
    decfloat_t re, im;
    int status;

    if (DECFLOAT_IS_ZERO(&x->im) && DECFLOAT_IS_FINITE(&x->re))
    {
        status = _decfloat_mul(&res->re, &x->re, &x->re, PREC(ctx), RND(ctx), NULL, NULL, ctx);
        decfloat_zero(&res->im, ctx);
        return status;
    }

    decfloat_init(re, ctx);
    decfloat_init(im, ctx);

    status = _decfloat_dot2(re, &x->re, &x->re, &x->im, &x->im, 1, PREC(ctx), RND(ctx), NULL, NULL, ctx);

    /* im = 2ab, rounded once */
    if (DECFLOAT_IS_ZERO(&x->re))
        decfloat_zero(im, ctx);
    else
    {
        status |= _decfloat_mul(im, &x->re, &x->im, DECIMAL_PREC_EXACT, EXACT_RND, NULL, NULL, ctx);
        status |= _decfloat_add(im, im, im, 0, PREC(ctx), RND_IM(ctx), NULL, NULL, ctx);
    }

    decfloat_swap(&res->re, re, ctx);
    decfloat_swap(&res->im, im, ctx);
    decfloat_clear(re, ctx);
    decfloat_clear(im, ctx);
    return status;
}

int
deccfloat_div(deccfloat_t res, const deccfloat_t x, const deccfloat_t y, gr_ctx_t ctx)
{
    decfloat_t re, im, den;
    int status;

    decfloat_init(re, ctx);
    decfloat_init(im, ctx);

    if (DECFLOAT_IS_ZERO(&y->im) && DECFLOAT_IS_FINITE(&y->re))
    {
        /* real y; a zero part of x stays zero (1 / 0 = inf, i / 0 = i inf) */
        if (DECFLOAT_IS_ZERO(&x->re) && !DECFLOAT_IS_ZERO(&x->im))
            status = decfloat_zero(re, ctx);
        else
            status = _decfloat_div(re, &x->re, &y->re, PREC(ctx), RND(ctx), NULL, NULL, ctx);
        if (DECFLOAT_IS_ZERO(&x->im) && !DECFLOAT_IS_ZERO(&x->re))
            status |= decfloat_zero(im, ctx);
        else
            status |= _decfloat_div(im, &x->im, &y->re, PREC(ctx), RND_IM(ctx), NULL, NULL, ctx);
    }
    else if (DECFLOAT_IS_ZERO(&y->re) && DECFLOAT_IS_FINITE(&y->im))
    {
        /* imaginary y: (a + bi) / (di) = b/d + a/(-d) i */
        decfloat_t t;
        decfloat_init(t, ctx);
        GR_MUST_SUCCEED(decfloat_neg_round(t, &y->im, DECIMAL_PREC_EXACT, EXACT_RND, ctx));
        status = _decfloat_div(re, &x->im, &y->im, PREC(ctx), RND(ctx), NULL, NULL, ctx);
        status |= _decfloat_div(im, &x->re, t, PREC(ctx), RND_IM(ctx), NULL, NULL, ctx);
        decfloat_clear(t, ctx);
    }
    else
    {
        /* (a + bi)(c - di) / (c^2 + d^2), numerators and denominator exact */
        decfloat_init(den, ctx);
        status = _decfloat_dot2(den, &y->re, &y->re, &y->im, &y->im, 0, DECIMAL_PREC_EXACT, EXACT_RND, NULL, NULL, ctx);
        status |= _decfloat_dot2(re, &x->re, &y->re, &x->im, &y->im, 0, DECIMAL_PREC_EXACT, EXACT_RND, NULL, NULL, ctx);
        status |= _decfloat_dot2(im, &x->im, &y->re, &x->re, &y->im, 1, DECIMAL_PREC_EXACT, EXACT_RND, NULL, NULL, ctx);
        if (status == GR_SUCCESS)
        {
            status = _decfloat_div(re, re, den, PREC(ctx), RND(ctx), NULL, NULL, ctx);
            status |= _decfloat_div(im, im, den, PREC(ctx), RND_IM(ctx), NULL, NULL, ctx);
        }
        decfloat_clear(den, ctx);
    }

    decfloat_swap(&res->re, re, ctx);
    decfloat_swap(&res->im, im, ctx);
    decfloat_clear(re, ctx);
    decfloat_clear(im, ctx);
    return status;
}

int
deccfloat_inv(deccfloat_t res, const deccfloat_t x, gr_ctx_t ctx)
{
    deccfloat_t one;
    int status;
    deccfloat_init(one, ctx);
    GR_MUST_SUCCEED(deccfloat_one(one, ctx));
    status = deccfloat_div(res, one, x, ctx);
    deccfloat_clear(one, ctx);
    return status;
}

int
deccfloat_sqrt(deccfloat_t res, const deccfloat_t x, gr_ctx_t ctx)
{
    int status;

    if (DECFLOAT_IS_ZERO(&x->im))
    {
        if (DECFLOAT_IS_NAN(&x->re) || _decfloat_sgn(&x->re, ctx) >= 0)
        {
            status = decfloat_sqrt_round(&res->re, &x->re, PREC(ctx), RND(ctx), ctx);
            decfloat_zero(&res->im, ctx);
        }
        else
        {
            decfloat_t t;
            decfloat_init(t, ctx);
            GR_MUST_SUCCEED(decfloat_neg_round(t, &x->re, DECIMAL_PREC_EXACT, EXACT_RND, ctx));
            status = decfloat_sqrt_round(&res->im, t, PREC(ctx), RND_IM(ctx), ctx);
            decfloat_zero(&res->re, ctx);
            decfloat_clear(t, ctx);
        }
        return status;
    }

    if (DECFLOAT_IS_ZERO(&x->re) && DECFLOAT_IS_FINITE(&x->im))
    {
        /* sqrt(bi) = sqrt(|b|/2) (1 +/- i) */
        decfloat_t t, re, im;
        int sgn = _decfloat_sgn(&x->im, ctx);

        decfloat_init(t, ctx);
        decfloat_init(re, ctx);
        decfloat_init(im, ctx);
        GR_MUST_SUCCEED(decfloat_abs_round(t, &x->im, DECIMAL_PREC_EXACT, EXACT_RND, ctx));
        GR_MUST_SUCCEED(_decfloat_half_exact(t, t, ctx));
        status = decfloat_sqrt_round(re, t, PREC(ctx), RND(ctx), ctx);
        status |= decfloat_sqrt_round(im, t, PREC(ctx), (sgn > 0) ? RND_IM(ctx) : DECIMAL_RND_NEGATE(RND_IM(ctx)), ctx);
        if (sgn < 0 && status == GR_SUCCESS)
            im->m.size = -im->m.size;
        decfloat_swap(&res->re, re, ctx);
        decfloat_swap(&res->im, im, ctx);
        decfloat_clear(t, ctx);
        decfloat_clear(re, ctx);
        decfloat_clear(im, ctx);
        return status;
    }

    if (_deccfloat_is_finite(x))
    {
        deccfloat_t t;
        deccfloat_init(t, ctx);
        status = _deccfloat_sqrt_exact(t, x, ctx);
        if (status == GR_SUCCESS)
            status = deccfloat_set(res, t, ctx);
        deccfloat_clear(t, ctx);
        if (status == GR_SUCCESS)
            return status;
    }

    {
        deccfloat_srcptr a[1];
        deccfloat_ptr r[1];
        a[0] = x;
        r[0] = res;
        return _deccfloat_acb_gr(r, 1, a, 1, 0, 0, GR_METHOD_SQRT, ctx);
    }
}

int
deccfloat_rsqrt(deccfloat_t res, const deccfloat_t x, gr_ctx_t ctx)
{
    int status;

    if (DECFLOAT_IS_ZERO(&x->im))
    {
        if (DECFLOAT_IS_NAN(&x->re) || _decfloat_sgn(&x->re, ctx) >= 0)
        {
            status = decfloat_rsqrt_round(&res->re, &x->re, PREC(ctx), RND(ctx), ctx);
            decfloat_zero(&res->im, ctx);
        }
        else
        {
            /* 1/sqrt(-t) = -i / sqrt(t) */
            decfloat_t t;
            decfloat_init(t, ctx);
            GR_MUST_SUCCEED(decfloat_neg_round(t, &x->re, DECIMAL_PREC_EXACT, EXACT_RND, ctx));
            status = decfloat_rsqrt_round(&res->im, t, PREC(ctx), DECIMAL_RND_NEGATE(RND_IM(ctx)), ctx);
            if (status == GR_SUCCESS)
                res->im.m.size = -res->im.m.size;
            decfloat_zero(&res->re, ctx);
            decfloat_clear(t, ctx);
        }
        return status;
    }

    if (_deccfloat_is_finite(x))
    {
        deccfloat_t t;
        deccfloat_init(t, ctx);
        status = _deccfloat_sqrt_exact(t, x, ctx);
        if (status == GR_SUCCESS)
            status = deccfloat_inv(res, t, ctx);
        deccfloat_clear(t, ctx);
        if (status == GR_SUCCESS)
            return status;
    }

    {
        deccfloat_srcptr a[1];
        deccfloat_ptr r[1];
        a[0] = x;
        r[0] = res;
        return _deccfloat_acb_gr(r, 1, a, 1, 0, 0, GR_METHOD_RSQRT, ctx);
    }
}

int
deccfloat_abs(deccfloat_t res, const deccfloat_t x, gr_ctx_t ctx)
{
    int status;

    if (DECFLOAT_IS_ZERO(&x->im))
        status = decfloat_abs_round(&res->re, &x->re, PREC(ctx), RND(ctx), ctx);
    else if (DECFLOAT_IS_ZERO(&x->re))
        status = decfloat_abs_round(&res->re, &x->im, PREC(ctx), RND(ctx), ctx);
    else
    {
        decfloat_t t;
        decfloat_init(t, ctx);
        status = _decfloat_dot2(t, &x->re, &x->re, &x->im, &x->im, 0, DECIMAL_PREC_EXACT, EXACT_RND, NULL, NULL, ctx);
        if (status == GR_SUCCESS)
            status = _decfloat_sqrt(&res->re, t, PREC(ctx), RND(ctx), NULL, NULL, ctx);
        decfloat_clear(t, ctx);
    }

    decfloat_zero(&res->im, ctx);
    return status;
}

int
deccfloat_arg(deccfloat_t res, const deccfloat_t x, gr_ctx_t ctx)
{
    decfloat_t t;
    int status;
    decfloat_init(t, ctx);
    status = decfloat_atan2(t, &x->im, &x->re, ctx);
    decfloat_swap(&res->re, t, ctx);
    decfloat_zero(&res->im, ctx);
    decfloat_clear(t, ctx);
    return status;
}

int
deccfloat_csgn(deccfloat_t res, const deccfloat_t x, gr_ctx_t ctx)
{
    int s;

    if (_deccfloat_is_nan(x))
        return deccfloat_nan(res, ctx);

    s = _decfloat_sgn(&x->re, ctx);
    if (s == 0)
        s = _decfloat_sgn(&x->im, ctx);

    return deccfloat_set_si(res, s, ctx);
}

int
deccfloat_sgn(deccfloat_t res, const deccfloat_t x, gr_ctx_t ctx)
{
    int status;

    if (_deccfloat_is_nan(x))
        return deccfloat_nan(res, ctx);

    if (DECFLOAT_IS_ZERO(&x->im))
        return deccfloat_set_si(res, _decfloat_sgn(&x->re, ctx), ctx);

    if (DECFLOAT_IS_ZERO(&x->re))
    {
        decfloat_zero(&res->re, ctx);
        return decfloat_set_si(&res->im, _decfloat_sgn(&x->im, ctx), ctx);
    }

    if (_deccfloat_is_finite(x))
    {
        /* exact modulus: correctly rounded quotients */
        decfloat_t r;
        decfloat_init(r, ctx);
        status = _decfloat_dot2(r, &x->re, &x->re, &x->im, &x->im, 0, DECIMAL_PREC_EXACT, EXACT_RND, NULL, NULL, ctx);
        if (status == GR_SUCCESS)
            status = _decfloat_sqrt(r, r, DECIMAL_PREC_EXACT, EXACT_RND, NULL, NULL, ctx);
        if (status == GR_SUCCESS)
        {
            decfloat_t re, im;
            decfloat_init(re, ctx);
            decfloat_init(im, ctx);
            status = _decfloat_div(re, &x->re, r, PREC(ctx), RND(ctx), NULL, NULL, ctx);
            status |= _decfloat_div(im, &x->im, r, PREC(ctx), RND_IM(ctx), NULL, NULL, ctx);
            decfloat_swap(&res->re, re, ctx);
            decfloat_swap(&res->im, im, ctx);
            decfloat_clear(re, ctx);
            decfloat_clear(im, ctx);
        }
        decfloat_clear(r, ctx);
        if (status == GR_SUCCESS)
            return status;
    }

    {
        deccfloat_srcptr a[1];
        deccfloat_ptr r[1];
        a[0] = x;
        r[0] = res;
        return _deccfloat_acb_gr(r, 1, a, 1, 0, 0, GR_METHOD_SGN, ctx);
    }
}

typedef int (*deccfloat_binary_fn)(deccfloat_t, const deccfloat_t, const deccfloat_t, gr_ctx_t);

/* op(x, y) for an integer y of type DECIMAL_SCALAR_* */
DECIMAL_DRIVER int
_deccfloat_scalar_op(deccfloat_t res, const deccfloat_t x, const void * y, int type, deccfloat_binary_fn op, gr_ctx_t ctx)
{
    deccfloat_t t;
    int status;
    deccfloat_init(t, ctx);
    status = _decfloat_set_scalar_exact(&t->re, y, type, ctx);
    if (status == GR_SUCCESS)
        status = op(res, x, t, ctx);
    deccfloat_clear(t, ctx);
    return status;
}

#define DEF_SCALAR(name, T, type, op) \
int name(deccfloat_t res, const deccfloat_t x, T y, gr_ctx_t ctx) \
{ \
    return _deccfloat_scalar_op(res, x, &y, type, op, ctx); \
}

#define DEF_SCALAR_FMPZ(name, op) \
int name(deccfloat_t res, const deccfloat_t x, const fmpz_t y, gr_ctx_t ctx) \
{ \
    return _deccfloat_scalar_op(res, x, y, 2, op, ctx); \
}

DEF_SCALAR(deccfloat_add_ui, ulong, 0, deccfloat_add)
DEF_SCALAR(deccfloat_add_si, slong, 1, deccfloat_add)
DEF_SCALAR_FMPZ(deccfloat_add_fmpz, deccfloat_add)
DEF_SCALAR(deccfloat_sub_ui, ulong, 0, deccfloat_sub)
DEF_SCALAR(deccfloat_sub_si, slong, 1, deccfloat_sub)
DEF_SCALAR_FMPZ(deccfloat_sub_fmpz, deccfloat_sub)
DEF_SCALAR(deccfloat_mul_ui, ulong, 0, deccfloat_mul)
DEF_SCALAR(deccfloat_mul_si, slong, 1, deccfloat_mul)
DEF_SCALAR_FMPZ(deccfloat_mul_fmpz, deccfloat_mul)
DEF_SCALAR(deccfloat_div_ui, ulong, 0, deccfloat_div)
DEF_SCALAR(deccfloat_div_si, slong, 1, deccfloat_div)
DEF_SCALAR_FMPZ(deccfloat_div_fmpz, deccfloat_div)

int
deccfloat_mul_decfloat(deccfloat_t res, const deccfloat_t x, const decfloat_t y, gr_ctx_t ctx)
{
    int status;
    status = _decfloat_mul(&res->re, &x->re, y, PREC(ctx), RND(ctx), NULL, NULL, ctx);
    status |= _decfloat_mul(&res->im, &x->im, y, PREC(ctx), RND_IM(ctx), NULL, NULL, ctx);
    return status;
}

int
deccfloat_div_decfloat(deccfloat_t res, const deccfloat_t x, const decfloat_t y, gr_ctx_t ctx)
{
    int status;
    status = _decfloat_div(&res->re, &x->re, y, PREC(ctx), RND(ctx), NULL, NULL, ctx);
    status |= _decfloat_div(&res->im, &x->im, y, PREC(ctx), RND_IM(ctx), NULL, NULL, ctx);
    return status;
}

int
deccfloat_mul_i(deccfloat_t res, const deccfloat_t x, gr_ctx_t ctx)
{
    decfloat_t t;
    int status;
    decfloat_init(t, ctx);
    /* (a + bi) i = -b + ai */
    status = decfloat_neg_round(t, &x->im, PREC(ctx), RND(ctx), ctx);
    status |= decfloat_set_round(&res->im, &x->re, PREC(ctx), RND_IM(ctx), ctx);
    decfloat_swap(&res->re, t, ctx);
    decfloat_clear(t, ctx);
    return status;
}

int
deccfloat_div_i(deccfloat_t res, const deccfloat_t x, gr_ctx_t ctx)
{
    decfloat_t t;
    int status;
    decfloat_init(t, ctx);
    /* (a + bi) / i = b - ai */
    status = decfloat_set_round(t, &x->im, PREC(ctx), RND(ctx), ctx);
    status |= decfloat_neg_round(&res->im, &x->re, PREC(ctx), RND_IM(ctx), ctx);
    decfloat_swap(&res->re, t, ctx);
    decfloat_clear(t, ctx);
    return status;
}

int
deccfloat_mul_two(deccfloat_t res, const deccfloat_t x, gr_ctx_t ctx)
{
    return deccfloat_add(res, x, x, ctx);
}

int
deccfloat_mul_10exp_fmpz(deccfloat_t res, const deccfloat_t x, const fmpz_t e, gr_ctx_t ctx)
{
    int status;
    status = decfloat_mul_10exp_fmpz_round(&res->re, &x->re, e, PREC(ctx), RND(ctx), ctx);
    status |= decfloat_mul_10exp_fmpz_round(&res->im, &x->im, e, PREC(ctx), RND_IM(ctx), ctx);
    return status;
}

int
deccfloat_mul_10exp_si(deccfloat_t res, const deccfloat_t x, slong e, gr_ctx_t ctx)
{
    int status;
    status = decfloat_mul_10exp_si_round(&res->re, &x->re, e, PREC(ctx), RND(ctx), ctx);
    status |= decfloat_mul_10exp_si_round(&res->im, &x->im, e, PREC(ctx), RND_IM(ctx), ctx);
    return status;
}

int
deccfloat_mul_2exp_fmpz(deccfloat_t res, const deccfloat_t x, const fmpz_t e, gr_ctx_t ctx)
{
    int status;
    status = decfloat_mul_2exp_fmpz_round(&res->re, &x->re, e, PREC(ctx), RND(ctx), ctx);
    status |= decfloat_mul_2exp_fmpz_round(&res->im, &x->im, e, PREC(ctx), RND_IM(ctx), ctx);
    return status;
}

int
deccfloat_mul_2exp_si(deccfloat_t res, const deccfloat_t x, slong e, gr_ctx_t ctx)
{
    fmpz_t t;
    int status;
    fmpz_init_set_si(t, e);
    status = deccfloat_mul_2exp_fmpz(res, x, t, ctx);
    fmpz_clear(t);
    return status;
}

#define DEF_ROUND_TO_INT(name, fn) \
int name(deccfloat_t res, const deccfloat_t x, gr_ctx_t ctx) \
{ \
    if (!_deccfloat_is_real(x)) \
        return _nonreal_status(x); \
    decfloat_zero(&res->im, ctx); \
    return fn(&res->re, &x->re, PREC(ctx), RND(ctx), ctx); \
}

DEF_ROUND_TO_INT(deccfloat_floor, decfloat_floor_round)
DEF_ROUND_TO_INT(deccfloat_ceil, decfloat_ceil_round)
DEF_ROUND_TO_INT(deccfloat_trunc, decfloat_trunc_round)
DEF_ROUND_TO_INT(deccfloat_nint, decfloat_nint_round)

int
deccfloat_cmp(int * res, const deccfloat_t x, const deccfloat_t y, gr_ctx_t ctx)
{
    if (!_deccfloat_is_real(x) || !_deccfloat_is_real(y))
    {
        *res = 0;
        return (_deccfloat_is_nan(x) || _deccfloat_is_nan(y)) ? GR_UNABLE : GR_DOMAIN;
    }

    return decfloat_cmp(res, &x->re, &y->re, ctx);
}

int
deccfloat_cmpabs(int * res, const deccfloat_t x, const deccfloat_t y, gr_ctx_t ctx)
{
    decfloat_t s, t;
    int status;

    if (_deccfloat_is_real(x) && _deccfloat_is_real(y))
        return decfloat_cmpabs(res, &x->re, &y->re, ctx);

    if (_deccfloat_is_nan(x) || _deccfloat_is_nan(y))
    {
        *res = 0;
        return GR_UNABLE;
    }

    decfloat_init(s, ctx);
    decfloat_init(t, ctx);
    status = _decfloat_dot2(s, &x->re, &x->re, &x->im, &x->im, 0, DECIMAL_PREC_EXACT, EXACT_RND, NULL, NULL, ctx);
    status |= _decfloat_dot2(t, &y->re, &y->re, &y->im, &y->im, 0, DECIMAL_PREC_EXACT, EXACT_RND, NULL, NULL, ctx);
    if (status == GR_SUCCESS)
        status = decfloat_cmp(res, s, t, ctx);
    else
        *res = 0;
    decfloat_clear(s, ctx);
    decfloat_clear(t, ctx);
    return status;
}

/* The remaining functions are not performance-critical. */
PUSH_OPTIONS
OPTIMIZE_OSIZE

/* ------------------------------------------------------------------------- */
/*    Basic functions                                                        */
/* ------------------------------------------------------------------------- */

int
deccfloat_one(deccfloat_t res, gr_ctx_t ctx)
{
    decfloat_zero(&res->im, ctx);
    return decfloat_one(&res->re, ctx);
}

int
deccfloat_neg_one(deccfloat_t res, gr_ctx_t ctx)
{
    decfloat_zero(&res->im, ctx);
    return decfloat_neg_one(&res->re, ctx);
}

int
deccfloat_i(deccfloat_t res, gr_ctx_t ctx)
{
    decfloat_zero(&res->re, ctx);
    return decfloat_one(&res->im, ctx);
}

int
deccfloat_nan(deccfloat_t res, gr_ctx_t ctx)
{
    int status = decfloat_nan(&res->re, ctx);
    status |= decfloat_nan(&res->im, ctx);
    return status;
}

int
deccfloat_pos_inf(deccfloat_t res, gr_ctx_t ctx)
{
    decfloat_zero(&res->im, ctx);
    return decfloat_pos_inf(&res->re, ctx);
}

int
deccfloat_neg_inf(deccfloat_t res, gr_ctx_t ctx)
{
    decfloat_zero(&res->im, ctx);
    return decfloat_neg_inf(&res->re, ctx);
}

truth_t
deccfloat_is_zero(const deccfloat_t x, gr_ctx_t ctx)
{
    truth_t a = decfloat_is_zero(&x->re, ctx);
    truth_t b = decfloat_is_zero(&x->im, ctx);
    return truth_and(a, b);
}

truth_t
deccfloat_is_one(const deccfloat_t x, gr_ctx_t ctx)
{
    return truth_and(decfloat_is_one(&x->re, ctx), decfloat_is_zero(&x->im, ctx));
}

truth_t
deccfloat_is_neg_one(const deccfloat_t x, gr_ctx_t ctx)
{
    return truth_and(decfloat_is_neg_one(&x->re, ctx), decfloat_is_zero(&x->im, ctx));
}

truth_t
deccfloat_equal(const deccfloat_t x, const deccfloat_t y, gr_ctx_t ctx)
{
    return truth_and(decfloat_equal(&x->re, &y->re, ctx), decfloat_equal(&x->im, &y->im, ctx));
}

truth_t
deccfloat_is_integer(const deccfloat_t x, gr_ctx_t ctx)
{
    return truth_and(decfloat_is_integer(&x->re, ctx), decfloat_is_zero(&x->im, ctx));
}

/* ------------------------------------------------------------------------- */
/*    Assignment and conversion                                              */
/* ------------------------------------------------------------------------- */

int
deccfloat_set_round(deccfloat_t res, const deccfloat_t x, slong prec, int rnd, int rnd_im, gr_ctx_t ctx)
{
    int status;
    status = decfloat_set_round(&res->re, &x->re, prec, rnd, ctx);
    status |= decfloat_set_round(&res->im, &x->im, prec, rnd_im, ctx);
    return status;
}

int
deccfloat_set(deccfloat_t res, const deccfloat_t x, gr_ctx_t ctx)
{
    return deccfloat_set_round(res, x, PREC(ctx), RND(ctx), RND_IM(ctx), ctx);
}

int
deccfloat_set_decfloat(deccfloat_t res, const decfloat_t x, gr_ctx_t ctx)
{
    decfloat_zero(&res->im, ctx);
    return decfloat_set_round(&res->re, x, PREC(ctx), RND(ctx), ctx);
}

int
deccfloat_set_decfloat_decfloat(deccfloat_t res, const decfloat_t re, const decfloat_t im, gr_ctx_t ctx)
{
    int status;
    status = decfloat_set_round(&res->re, re, PREC(ctx), RND(ctx), ctx);
    status |= decfloat_set_round(&res->im, im, PREC(ctx), RND_IM(ctx), ctx);
    return status;
}

#define DEF_SET_REAL(name, T, fn) \
int name(deccfloat_t res, T x, gr_ctx_t ctx) \
{ \
    decfloat_zero(&res->im, ctx); \
    return fn(&res->re, x, ctx); \
}

DEF_SET_REAL(deccfloat_set_si, slong, decfloat_set_si)
DEF_SET_REAL(deccfloat_set_ui, ulong, decfloat_set_ui)
DEF_SET_REAL(deccfloat_set_fmpz, const fmpz_t, decfloat_set_fmpz)
DEF_SET_REAL(deccfloat_set_fmpq, const fmpq_t, decfloat_set_fmpq)
DEF_SET_REAL(deccfloat_set_d, double, decfloat_set_d)

int
deccfloat_set_fmpz_10exp_fmpz(deccfloat_t res, const fmpz_t m, const fmpz_t e, gr_ctx_t ctx)
{
    decfloat_zero(&res->im, ctx);
    return decfloat_set_fmpz_10exp_fmpz(&res->re, m, e, ctx);
}

int
deccfloat_set_fmpz_2exp_fmpz(deccfloat_t res, const fmpz_t m, const fmpz_t e, gr_ctx_t ctx)
{
    decfloat_zero(&res->im, ctx);
    return decfloat_set_fmpz_2exp_fmpz(&res->re, m, e, ctx);
}

int
deccfloat_set_str(deccfloat_t res, const char * s, gr_ctx_t ctx)
{
    int status;

    /* plain real literals; everything else goes through the expression
       parser, which understands I */
    status = _decfloat_set_str_literal(&res->re, s, PREC(ctx), RND(ctx), ctx);

    if (status == GR_SUCCESS)
    {
        decfloat_zero(&res->im, ctx);
        return GR_SUCCESS;
    }

    if (status == GR_DOMAIN)
        return gr_generic_set_str_ring_exponents(res, s, ctx);

    return status;
}

int
deccfloat_set_acf(deccfloat_t res, const acf_t x, gr_ctx_t ctx)
{
    int status;
    status = decfloat_set_round_arf(&res->re, acf_realref(x), PREC(ctx), RND(ctx), ctx);
    status |= decfloat_set_round_arf(&res->im, acf_imagref(x), PREC(ctx), RND_IM(ctx), ctx);
    return status;
}

int
deccfloat_set_acb(deccfloat_t res, const acb_t x, gr_ctx_t ctx)
{
    acf_t t;

    if (DECIMAL_CTX_IS_EXACT(ctx) && !acb_is_exact(x))
        return GR_UNABLE;

    *acf_realref(t) = *arb_midref(acb_realref(x));
    *acf_imagref(t) = *arb_midref(acb_imagref(x));
    return deccfloat_set_acf(res, t, ctx);
}

int
deccfloat_get_acf(acf_t res, const deccfloat_t x, slong prec_bits, int rnd, gr_ctx_t ctx)
{
    int status;
    status = decfloat_get_arf(acf_realref(res), &x->re, prec_bits, rnd, ctx);
    status |= decfloat_get_arf(acf_imagref(res), &x->im, prec_bits, rnd, ctx);
    return status;
}

int
deccfloat_get_acb(acb_t res, const deccfloat_t x, slong prec_bits, gr_ctx_t ctx)
{
    int status;
    status = decfloat_get_arb(acb_realref(res), &x->re, prec_bits, ctx);
    status |= decfloat_get_arb(acb_imagref(res), &x->im, prec_bits, ctx);
    return status;
}

int
deccfloat_get_fmpz(fmpz_t res, const deccfloat_t x, gr_ctx_t ctx)
{
    if (!_deccfloat_is_real(x))
        return _nonreal_status(x);
    return decfloat_get_fmpz(res, &x->re, ctx);
}

int
deccfloat_get_fmpq(fmpq_t res, const deccfloat_t x, gr_ctx_t ctx)
{
    if (!_deccfloat_is_real(x))
        return _nonreal_status(x);
    return decfloat_get_fmpq(res, &x->re, ctx);
}

int
deccfloat_get_si(slong * res, const deccfloat_t x, gr_ctx_t ctx)
{
    if (!_deccfloat_is_real(x))
        return _nonreal_status(x);
    return decfloat_get_si(res, &x->re, ctx);
}

int
deccfloat_get_ui(ulong * res, const deccfloat_t x, gr_ctx_t ctx)
{
    if (!_deccfloat_is_real(x))
        return _nonreal_status(x);
    return decfloat_get_ui(res, &x->re, ctx);
}

int
deccfloat_get_d(double * res, const deccfloat_t x, gr_ctx_t ctx)
{
    if (!_deccfloat_is_real(x))
        return _nonreal_status(x);
    return decfloat_get_d(res, &x->re, ctx);
}

int
deccfloat_get_re(decfloat_t res, const deccfloat_t x, gr_ctx_t ctx)
{
    return decfloat_set_round(res, &x->re, DECIMAL_PREC_EXACT, EXACT_RND, ctx);
}

int
deccfloat_get_im(decfloat_t res, const deccfloat_t x, gr_ctx_t ctx)
{
    return decfloat_set_round(res, &x->im, DECIMAL_PREC_EXACT, EXACT_RND, ctx);
}

/* exact conversion of a decfloat from another decimal context */
static int
_set_decfloat_other_exact(decfloat_t res, const decfloat_t x, gr_ctx_t x_ctx, gr_ctx_t ctx)
{
    slong saved = PREC(ctx);
    int status;
    PREC(ctx) = DECIMAL_PREC_EXACT;
    status = _decfloat_set_decfloat_other(res, x, x_ctx, ctx);
    PREC(ctx) = saved;
    return status;
}

int
deccfloat_set_other(deccfloat_t res, gr_srcptr x, gr_ctx_t x_ctx, gr_ctx_t ctx)
{
    switch (x_ctx->which_ring)
    {
        case GR_CTX_FMPZ:
            return deccfloat_set_fmpz(res, x, ctx);

        case GR_CTX_FMPQ:
            return deccfloat_set_fmpq(res, x, ctx);

        case GR_CTX_REAL_FLOAT_ARF:
            decfloat_zero(&res->im, ctx);
            return decfloat_set_arf(&res->re, x, ctx);

        case GR_CTX_COMPLEX_FLOAT_ACF:
            return deccfloat_set_acf(res, x, ctx);

        case GR_CTX_RR_ARB:
            if (DECIMAL_CTX_IS_EXACT(ctx) && !arb_is_exact((arb_srcptr) x))
                return GR_UNABLE;
            decfloat_zero(&res->im, ctx);
            return decfloat_set_arf(&res->re, arb_midref((arb_srcptr) x), ctx);

        case GR_CTX_CC_ACB:
            return deccfloat_set_acb(res, x, ctx);

        case GR_CTX_DECFLOAT:
            decfloat_zero(&res->im, ctx);
            return decfloat_set_other(&res->re, x, x_ctx, ctx);

        case GR_CTX_DECBALL:
            decfloat_zero(&res->im, ctx);
            return decfloat_set_other(&res->re, x, x_ctx, ctx);

        case GR_CTX_DECCFLOAT:
        case GR_CTX_DECCBALL:
            {
                deccfloat_t t;
                decfloat_srcptr re, im;
                int status;

                if (x_ctx->which_ring == GR_CTX_DECCFLOAT)
                {
                    re = DECCFLOAT_REALREF((deccfloat_srcptr) x);
                    im = DECCFLOAT_IMAGREF((deccfloat_srcptr) x);
                }
                else
                {
                    if (DECIMAL_CTX_IS_EXACT(ctx) && !_deccball_is_exact((deccball_srcptr) x, x_ctx))
                        return GR_UNABLE;
                    re = DECBALL_MIDREF(DECCBALL_REALREF((deccball_srcptr) x));
                    im = DECBALL_MIDREF(DECCBALL_IMAGREF((deccball_srcptr) x));
                }

                deccfloat_init(t, ctx);
                status = _set_decfloat_other_exact(&t->re, re, x_ctx, ctx);
                status |= _set_decfloat_other_exact(&t->im, im, x_ctx, ctx);
                if (status == GR_SUCCESS)
                    status = deccfloat_set(res, t, ctx);
                deccfloat_clear(t, ctx);
                return status;
            }

        case GR_CTX_REAL_ALGEBRAIC_QQBAR:
        case GR_CTX_COMPLEX_ALGEBRAIC_QQBAR:
            return deccfloat_set_qqbar(res, x, ctx);

        default:
            {
                gr_ctx_t cctx;
                acb_t z;
                int status;

                gr_ctx_init_complex_acb(cctx, 20 + (DECIMAL_CTX_IS_EXACT(ctx) ? 64 : 4 * PREC(ctx)));
                acb_init(z);

                status = gr_set_other(z, x, x_ctx, cctx);

                if (status == GR_SUCCESS)
                    status = deccfloat_set_acb(res, z, ctx);

                acb_clear(z);
                gr_ctx_clear(cctx);

                return status;
            }
    }
}

/* ------------------------------------------------------------------------- */
/*    Output, random                                                         */
/* ------------------------------------------------------------------------- */

char *
deccfloat_get_str(const deccfloat_t x, gr_ctx_t ctx)
{
    char * a, * b, * s;

    if (DECFLOAT_IS_ZERO(&x->im))
        return decfloat_get_str(&x->re, ctx);

    if (DECFLOAT_IS_ZERO(&x->re))
    {
        b = decfloat_get_str(&x->im, ctx);
        s = flint_malloc(strlen(b) + 3);
        strcpy(s, b);
        strcat(s, "*I");
        flint_free(b);
        return s;
    }

    a = decfloat_get_str(&x->re, ctx);

    if (DECFLOAT_SGNBIT(&x->im) || DECFLOAT_IS_NEG_INF(&x->im))
    {
        decfloat_t t;
        decfloat_init(t, ctx);
        GR_MUST_SUCCEED(decfloat_neg_round(t, &x->im, DECIMAL_PREC_EXACT, EXACT_RND, ctx));
        b = decfloat_get_str(t, ctx);
        decfloat_clear(t, ctx);
        s = flint_malloc(strlen(a) + strlen(b) + 8);
        strcpy(s, "(");
        strcat(s, a);
        strcat(s, " - ");
    }
    else
    {
        b = decfloat_get_str(&x->im, ctx);
        s = flint_malloc(strlen(a) + strlen(b) + 8);
        strcpy(s, "(");
        strcat(s, a);
        strcat(s, " + ");
    }

    strcat(s, b);
    strcat(s, "*I)");
    flint_free(a);
    flint_free(b);
    return s;
}

int
deccfloat_write(gr_stream_t out, const deccfloat_t x, gr_ctx_t ctx)
{
    return gr_stream_write_free(out, deccfloat_get_str(x, ctx));
}

int
deccfloat_randtest(deccfloat_t res, flint_rand_t state, gr_ctx_t ctx)
{
    int status;

    status = decfloat_randtest(&res->re, state, ctx);
    status |= decfloat_randtest(&res->im, state, ctx);

    /* exercise the real and imaginary special cases */
    switch (n_randint(state, 8))
    {
        case 0:
        case 1:
            decfloat_zero(&res->im, ctx);
            break;
        case 2:
            decfloat_zero(&res->re, ctx);
            break;
        default:
            break;
    }

    return status;
}

int
deccfloat_randtest_special(deccfloat_t res, flint_rand_t state, gr_ctx_t ctx)
{
    int status;

    status = decfloat_randtest_special(&res->re, state, ctx);
    status |= decfloat_randtest_special(&res->im, state, ctx);

    switch (n_randint(state, 8))
    {
        case 0:
        case 1:
            decfloat_zero(&res->im, ctx);
            break;
        case 2:
            decfloat_zero(&res->re, ctx);
            break;
        default:
            break;
    }

    return status;
}

/* ------------------------------------------------------------------------- */
/*    Powers                                                                 */
/* ------------------------------------------------------------------------- */

/* whether x^n is cheap enough to compute exactly */
static int
_pow_exact_feasible(const deccfloat_t x, const fmpz_t n, gr_ctx_t ctx)
{
    slong D;
    ulong an;

    /* the bound below requires |n| <= 100000 anyway */
    if (fmpz_bits(n) > 20)
        return 0;

    an = FLINT_ABS(fmpz_get_si(n));
    /* size of the parts as exact fractions */
    D = decfloat_digits(&x->re, ctx) + decfloat_digits(&x->im, ctx) + 2;
    if (!DECFLOAT_IS_SPECIAL(&x->re)) D += FLINT_MIN(FLINT_ABS(_decfloat_val10_clamped(&x->re, ctx)), 1000000);
    if (!DECFLOAT_IS_SPECIAL(&x->im)) D += FLINT_MIN(FLINT_ABS(_decfloat_val10_clamped(&x->im, ctx)), 1000000);

    return (ulong) D * an <= 100000;
}

int
deccfloat_pow_fmpz(deccfloat_t res, const deccfloat_t x, const fmpz_t n, gr_ctx_t ctx)
{
    int status;

    if (fmpz_is_zero(n))
        return deccfloat_one(res, ctx);

    if (_deccfloat_is_nan(x))
        return deccfloat_nan(res, ctx);

    if (fmpz_is_one(n))
        return deccfloat_set(res, x, ctx);

    if (DECFLOAT_IS_ZERO(&x->im))
    {
        decfloat_zero(&res->im, ctx);
        return decfloat_pow_fmpz(&res->re, &x->re, n, ctx);
    }

    if (DECFLOAT_IS_ZERO(&x->re) && DECFLOAT_IS_FINITE(&x->im))
    {
        /* (bi)^n = b^n i^n */
        decfloat_t t;
        int k = fmpz_fdiv_ui(n, 4), negative = (k >= 2), rnd;

        decfloat_init(t, ctx);
        rnd = (k & 1) ? RND_IM(ctx) : RND(ctx);
        if (negative)
            rnd = DECIMAL_RND_NEGATE(rnd);

        {
            int saved = RND(ctx);
            RND(ctx) = rnd;
            status = decfloat_pow_fmpz(t, &x->im, n, ctx);
            RND(ctx) = saved;
        }

        if (status == GR_SUCCESS)
        {
            if (negative)
            {
                if (DECFLOAT_IS_POS_INF(t)) _decfloat_neg_inf(t);
                else if (DECFLOAT_IS_NEG_INF(t)) _decfloat_pos_inf(t);
                else t->m.size = -t->m.size;
            }

            if (k & 1)
            {
                decfloat_zero(&res->re, ctx);
                decfloat_swap(&res->im, t, ctx);
            }
            else
            {
                decfloat_zero(&res->im, ctx);
                decfloat_swap(&res->re, t, ctx);
            }
        }

        decfloat_clear(t, ctx);
        return status;
    }

    if (_deccfloat_is_finite(x) && _pow_exact_feasible(x, n, ctx))
    {
        deccfloat_t t;
        deccfloat_init(t, ctx);
        status = _deccfloat_pow_ui_exact(t, x, FLINT_ABS(fmpz_get_si(n)), ctx);
        if (status == GR_SUCCESS)
        {
            if (fmpz_sgn(n) > 0)
                status = deccfloat_set(res, t, ctx);
            else
                status = deccfloat_inv(res, t, ctx);
        }
        deccfloat_clear(t, ctx);
        return status;
    }

    if (fmpz_equal_si(n, -1))
        return deccfloat_inv(res, x, ctx);

    {
        gr_ctx_t actx;
        acb_t a, r;
        slong prec = PREC(ctx), wp, wpbits;
        int rr = 0;

        if (prec == DECIMAL_PREC_EXACT)
            return GR_UNABLE;

        acb_init(a);
        acb_init(r);
        status = GR_UNABLE;

        for (wp = prec + 10; wp <= _decimal_ziv_wp_max(prec, FLINT_MAX(decfloat_digits(&x->re, ctx), decfloat_digits(&x->im, ctx))); wp *= 2)
        {
            deccball_t Y;
            gr_ctx_t bctx;

            wpbits = _decimal_digits_to_bits(wp);
            gr_ctx_init_complex_acb(actx, wpbits);
            status = deccfloat_get_acb(a, x, wpbits, ctx);
            if (status == GR_SUCCESS)
                status = gr_pow_fmpz(r, a, n, actx);
            gr_ctx_clear(actx);

            if (status != GR_SUCCESS)
            {
                if (status & GR_DOMAIN)
                    break;
                continue;
            }

            _gr_ctx_init_decimal(bctx, DECIMAL_CTX_CBALL, DECIMAL_CTX_E(ctx), wp, DECIMAL_RND_DOWN, 0);
            decimal_ctx_set_rad_prec(bctx, DECMAG_MAX_PREC);
            deccball_init(Y, bctx);

            if (acb_is_finite(r) && deccball_set_acb(Y, r, bctx) == GR_SUCCESS)
            {
                decfloat_t tre, tim;
                decfloat_init(tre, ctx);
                decfloat_init(tim, ctx);
                rr = _decfloat_round_ball(tre, &Y->re, prec, RND(ctx), bctx, ctx);
                if (rr == 1)
                    rr = _decfloat_round_ball(tim, &Y->im, prec, RND_IM(ctx), bctx, ctx);
                if (rr == 1)
                {
                    decfloat_swap(&res->re, tre, ctx);
                    decfloat_swap(&res->im, tim, ctx);
                    status = _deccfloat_finalize(res, ctx);
                }
                decfloat_clear(tre, ctx);
                decfloat_clear(tim, ctx);
            }
            else
                rr = 0;

            deccball_clear(Y, bctx);
            gr_ctx_clear(bctx);

            if (rr == 1)
                break;
            status = GR_UNABLE;
            if (rr == -1)
                break;
        }

        acb_clear(a);
        acb_clear(r);
        return status;
    }
}

int
deccfloat_pow_ui(deccfloat_t res, const deccfloat_t x, ulong n, gr_ctx_t ctx)
{
    fmpz_t t;
    int status;
    fmpz_init_set_ui(t, n);
    status = deccfloat_pow_fmpz(res, x, t, ctx);
    fmpz_clear(t);
    return status;
}

int
deccfloat_pow_si(deccfloat_t res, const deccfloat_t x, slong n, gr_ctx_t ctx)
{
    fmpz_t t;
    int status;
    fmpz_init_set_si(t, n);
    status = deccfloat_pow_fmpz(res, x, t, ctx);
    fmpz_clear(t);
    return status;
}

int
deccfloat_pow(deccfloat_t res, const deccfloat_t x, const deccfloat_t y, gr_ctx_t ctx)
{
    int status;

    if (_deccfloat_is_nan(x) || _deccfloat_is_nan(y))
        return deccfloat_nan(res, ctx);

    if (DECFLOAT_IS_ZERO(&y->im))
    {
        /* integer exponent */
        if (_decfloat_is_int(&y->re, ctx) && (DECFLOAT_IS_ZERO(&y->re) || _decfloat_sci_exp_clamped(&y->re, ctx) < DECIMAL_WORD_DIGITS - 1))
        {
            fmpz_t n;
            fmpz_init(n);
            if (decfloat_get_fmpz(n, &y->re, ctx) == GR_SUCCESS && fmpz_bits(n) <= 62)
            {
                status = deccfloat_pow_fmpz(res, x, n, ctx);
                fmpz_clear(n);
                return status;
            }
            fmpz_clear(n);
        }

        /* real base */
        if (DECFLOAT_IS_ZERO(&x->im))
        {
            decfloat_t t;
            decfloat_init(t, ctx);
            status = decfloat_pow(t, &x->re, &y->re, ctx);

            if (status == GR_DOMAIN && DECFLOAT_IS_FINITE(&x->re) && DECFLOAT_IS_FINITE(&y->re) && _decfloat_sgn(&x->re, ctx) < 0)
            {
                /* (-t)^(p/2) = t^(p/2) i^p: y = +/- (j + 1/2) has a single
                   digit 5 at 10^-1, and p mod 4 is given by the parity of
                   the digit of |y| at 10^0 */
                if (_decfloat_val10_clamped(&y->re, ctx) == -1 && decfloat_get_digit_si(&y->re, -1, ctx) == 5)
                {
                    decfloat_t u;
                    int j = decfloat_get_digit_si(&y->re, 0, ctx) % 2;
                    int k = ((j == 0) == (_decfloat_sgn(&y->re, ctx) > 0)) ? 1 : 3;   /* p mod 4 */
                    int rnd = (k == 1) ? RND_IM(ctx) : DECIMAL_RND_NEGATE(RND_IM(ctx));
                    int saved = RND(ctx);

                    decfloat_init(u, ctx);
                    GR_MUST_SUCCEED(decfloat_neg_round(u, &x->re, DECIMAL_PREC_EXACT, EXACT_RND, ctx));
                    RND(ctx) = rnd;
                    status = decfloat_pow(t, u, &y->re, ctx);
                    RND(ctx) = saved;

                    if (status == GR_SUCCESS)
                    {
                        if (k == 3)
                        {
                            if (DECFLOAT_IS_POS_INF(t)) _decfloat_neg_inf(t);
                            else if (DECFLOAT_IS_NEG_INF(t)) _decfloat_pos_inf(t);
                            else t->m.size = -t->m.size;
                        }
                        decfloat_zero(&res->re, ctx);
                        decfloat_swap(&res->im, t, ctx);
                    }
                    decfloat_clear(u, ctx);
                    decfloat_clear(t, ctx);
                    if (status == GR_SUCCESS)
                        return status;
                    goto ziv;
                }
            }

            if (status == GR_SUCCESS)
            {
                decfloat_swap(&res->re, t, ctx);
                decfloat_zero(&res->im, ctx);
            }
            decfloat_clear(t, ctx);
            if (status != GR_DOMAIN)
                return status;
            goto ziv;
        }

        /* y = +/- 1/2: a single digit 5 at 10^-1 */
        if (!DECFLOAT_IS_SPECIAL(&y->re) && _decfloat_sci_exp_clamped(&y->re, ctx) == -1
            && _decfloat_val10_clamped(&y->re, ctx) == -1 && decfloat_get_digit_si(&y->re, -1, ctx) == 5)
        {
            if (_decfloat_sgn(&y->re, ctx) > 0)
                return deccfloat_sqrt(res, x, ctx);
            else
                return deccfloat_rsqrt(res, x, ctx);
        }
    }

    /* zero base */
    if (DECFLOAT_IS_ZERO(&x->re) && DECFLOAT_IS_ZERO(&x->im))
    {
        if (_decfloat_sgn(&y->re, ctx) > 0)
            return deccfloat_zero(res, ctx);
        if (DECFLOAT_IS_ZERO(&y->re) && DECFLOAT_IS_ZERO(&y->im))
            return deccfloat_one(res, ctx);
        return GR_DOMAIN;
    }

ziv:
    {
        deccfloat_srcptr a[2];
        deccfloat_ptr r[1];
        a[0] = x;
        a[1] = y;
        r[0] = res;
        return _deccfloat_acb_gr(r, 1, a, 2, 0, 0, GR_METHOD_POW, ctx);
    }
}

/* ------------------------------------------------------------------------- */
/*    Vectors                                                                */
/* ------------------------------------------------------------------------- */

int
deccfloat_vec_dot(deccfloat_t res, const deccfloat_t initial, int subtract, deccfloat_srcptr vec1, deccfloat_srcptr vec2, slong len, gr_ctx_t ctx)
{
    return gr_generic_vec_dot(res, initial, subtract, vec1, vec2, len, ctx);
}

int
deccfloat_vec_dot_rev(deccfloat_t res, const deccfloat_t initial, int subtract, deccfloat_srcptr vec1, deccfloat_srcptr vec2, slong len, gr_ctx_t ctx)
{
    return gr_generic_vec_dot_rev(res, initial, subtract, vec1, vec2, len, ctx);
}

/* ------------------------------------------------------------------------- */
/*    Method table                                                           */
/* ------------------------------------------------------------------------- */

static truth_t
_deccfloat_ctx_is_ring(gr_ctx_t ctx)
{
    /* the exact ring Z[1/10][i] */
    return (DECIMAL_CTX_IS_EXACT(ctx) && !(DECIMAL_CTX_FLAGS(ctx) & (DECIMAL_ALLOW_INF | DECIMAL_ALLOW_NAN | DECIMAL_ALLOW_UNDERFLOW))) ? T_TRUE : T_FALSE;
}

static truth_t
_deccfloat_ctx_is_approx(gr_ctx_t ctx)
{
    return (_deccfloat_ctx_is_ring(ctx) == T_TRUE) ? T_FALSE : T_TRUE;
}

static truth_t
_deccfloat_ctx_is_exact(gr_ctx_t ctx)
{
    return DECIMAL_CTX_IS_EXACT(ctx) ? T_TRUE : T_FALSE;
}

static truth_t
_deccfloat_ctx_is_canonical(gr_ctx_t ctx)
{
    return (DECIMAL_CTX_FLAGS(ctx) & DECIMAL_ALLOW_NAN) ? T_FALSE : T_TRUE;
}

int _deccfloat_methods_initialized = 0;
gr_static_method_table _deccfloat_methods;

gr_method_tab_input _deccfloat_methods_input[] =
{
    {GR_METHOD_CTX_CLEAR,       (gr_funcptr) decimal_ctx_clear},
    {GR_METHOD_CTX_WRITE,       (gr_funcptr) decimal_ctx_write},
    {GR_METHOD_CTX_IS_RING,     (gr_funcptr) _deccfloat_ctx_is_ring},
    {GR_METHOD_CTX_IS_COMMUTATIVE_RING, (gr_funcptr) _deccfloat_ctx_is_ring},
    {GR_METHOD_CTX_IS_INTEGRAL_DOMAIN,  (gr_funcptr) _deccfloat_ctx_is_ring},
    {GR_METHOD_CTX_IS_FIELD,            (gr_funcptr) gr_generic_ctx_predicate_false},
    {GR_METHOD_CTX_IS_UNIQUE_FACTORIZATION_DOMAIN, (gr_funcptr) _deccfloat_ctx_is_ring},
    {GR_METHOD_CTX_IS_FINITE,   (gr_funcptr) gr_generic_ctx_predicate_false},
    {GR_METHOD_CTX_IS_FINITE_CHARACTERISTIC, (gr_funcptr) gr_generic_ctx_predicate_false},
    {GR_METHOD_CTX_IS_ALGEBRAICALLY_CLOSED, (gr_funcptr) gr_generic_ctx_predicate_false},
    {GR_METHOD_CTX_IS_ORDERED_RING, (gr_funcptr) gr_generic_ctx_predicate_false},
    {GR_METHOD_CTX_IS_APPROX_COMMUTATIVE_RING, (gr_funcptr) _deccfloat_ctx_is_approx},
    {GR_METHOD_CTX_IS_RATIONAL_VECTOR_SPACE, (gr_funcptr) _decimal_ctx_is_vector_space},
    {GR_METHOD_CTX_IS_REAL_VECTOR_SPACE, (gr_funcptr) _decimal_ctx_is_vector_space},
    {GR_METHOD_CTX_IS_COMPLEX_VECTOR_SPACE, (gr_funcptr) _decimal_ctx_is_vector_space},
    {GR_METHOD_CTX_IS_EXACT,    (gr_funcptr) _deccfloat_ctx_is_exact},
    {GR_METHOD_CTX_IS_CANONICAL, (gr_funcptr) _deccfloat_ctx_is_canonical},
    {GR_METHOD_CTX_HAS_REAL_PREC, (gr_funcptr) gr_generic_ctx_predicate_true},
    {GR_METHOD_CTX_SET_REAL_PREC, (gr_funcptr) _decimal_ctx_set_real_prec},
    {GR_METHOD_CTX_GET_REAL_PREC, (gr_funcptr) _decimal_ctx_get_real_prec},

    {GR_METHOD_INIT,            (gr_funcptr) deccfloat_init},
    {GR_METHOD_CLEAR,           (gr_funcptr) deccfloat_clear},
    {GR_METHOD_SWAP,            (gr_funcptr) deccfloat_swap},
    {GR_METHOD_SET_SHALLOW,     (gr_funcptr) deccfloat_set_shallow},
    {GR_METHOD_RANDTEST,        (gr_funcptr) deccfloat_randtest},
    {GR_METHOD_WRITE,           (gr_funcptr) deccfloat_write},
    {GR_METHOD_ZERO,            (gr_funcptr) deccfloat_zero},
    {GR_METHOD_ONE,             (gr_funcptr) deccfloat_one},
    {GR_METHOD_NEG_ONE,         (gr_funcptr) deccfloat_neg_one},
    {GR_METHOD_IS_ZERO,         (gr_funcptr) deccfloat_is_zero},
    {GR_METHOD_IS_ONE,          (gr_funcptr) deccfloat_is_one},
    {GR_METHOD_IS_NEG_ONE,      (gr_funcptr) deccfloat_is_neg_one},
    {GR_METHOD_EQUAL,           (gr_funcptr) deccfloat_equal},
    {GR_METHOD_SET,             (gr_funcptr) deccfloat_set},
    {GR_METHOD_SET_SI,          (gr_funcptr) deccfloat_set_si},
    {GR_METHOD_SET_UI,          (gr_funcptr) deccfloat_set_ui},
    {GR_METHOD_SET_FMPZ,        (gr_funcptr) deccfloat_set_fmpz},
    {GR_METHOD_SET_FMPQ,        (gr_funcptr) deccfloat_set_fmpq},
    {GR_METHOD_SET_D,           (gr_funcptr) deccfloat_set_d},
    {GR_METHOD_SET_STR,         (gr_funcptr) deccfloat_set_str},
    {GR_METHOD_SET_FMPZ_10EXP_FMPZ, (gr_funcptr) deccfloat_set_fmpz_10exp_fmpz},
    {GR_METHOD_SET_FMPZ_2EXP_FMPZ, (gr_funcptr) deccfloat_set_fmpz_2exp_fmpz},
    {GR_METHOD_SET_OTHER,       (gr_funcptr) deccfloat_set_other},
    {GR_METHOD_GET_FMPZ,        (gr_funcptr) deccfloat_get_fmpz},
    {GR_METHOD_GET_FMPQ,        (gr_funcptr) deccfloat_get_fmpq},
    {GR_METHOD_GET_UI,          (gr_funcptr) deccfloat_get_ui},
    {GR_METHOD_GET_SI,          (gr_funcptr) deccfloat_get_si},
    {GR_METHOD_GET_D,           (gr_funcptr) deccfloat_get_d},

    {GR_METHOD_NEG,             (gr_funcptr) deccfloat_neg},
    {GR_METHOD_ADD,             (gr_funcptr) deccfloat_add},
    {GR_METHOD_ADD_UI,          (gr_funcptr) deccfloat_add_ui},
    {GR_METHOD_ADD_SI,          (gr_funcptr) deccfloat_add_si},
    {GR_METHOD_ADD_FMPZ,        (gr_funcptr) deccfloat_add_fmpz},
    {GR_METHOD_SUB,             (gr_funcptr) deccfloat_sub},
    {GR_METHOD_SUB_UI,          (gr_funcptr) deccfloat_sub_ui},
    {GR_METHOD_SUB_SI,          (gr_funcptr) deccfloat_sub_si},
    {GR_METHOD_SUB_FMPZ,        (gr_funcptr) deccfloat_sub_fmpz},
    {GR_METHOD_MUL,             (gr_funcptr) deccfloat_mul},
    {GR_METHOD_MUL_UI,          (gr_funcptr) deccfloat_mul_ui},
    {GR_METHOD_MUL_SI,          (gr_funcptr) deccfloat_mul_si},
    {GR_METHOD_MUL_FMPZ,        (gr_funcptr) deccfloat_mul_fmpz},
    {GR_METHOD_MUL_TWO,         (gr_funcptr) deccfloat_mul_two},
    {GR_METHOD_SQR,             (gr_funcptr) deccfloat_sqr},
    {GR_METHOD_DIV,             (gr_funcptr) deccfloat_div},
    {GR_METHOD_DIV_UI,          (gr_funcptr) deccfloat_div_ui},
    {GR_METHOD_DIV_SI,          (gr_funcptr) deccfloat_div_si},
    {GR_METHOD_DIV_FMPZ,        (gr_funcptr) deccfloat_div_fmpz},
    {GR_METHOD_INV,             (gr_funcptr) deccfloat_inv},
    {GR_METHOD_MUL_2EXP_SI,     (gr_funcptr) deccfloat_mul_2exp_si},
    {GR_METHOD_MUL_2EXP_FMPZ,   (gr_funcptr) deccfloat_mul_2exp_fmpz},
    {GR_METHOD_POW_UI,          (gr_funcptr) deccfloat_pow_ui},
    {GR_METHOD_POW_SI,          (gr_funcptr) deccfloat_pow_si},
    {GR_METHOD_POW_FMPZ,        (gr_funcptr) deccfloat_pow_fmpz},
    {GR_METHOD_POW,             (gr_funcptr) deccfloat_pow},
    {GR_METHOD_SQRT,            (gr_funcptr) deccfloat_sqrt},
    {GR_METHOD_RSQRT,           (gr_funcptr) deccfloat_rsqrt},
    {GR_METHOD_POS_INF,         (gr_funcptr) deccfloat_pos_inf},
    {GR_METHOD_NEG_INF,         (gr_funcptr) deccfloat_neg_inf},
    {GR_METHOD_UINF,            (gr_funcptr) gr_not_in_domain},
    {GR_METHOD_UNDEFINED,       (gr_funcptr) deccfloat_nan},
    {GR_METHOD_UNKNOWN,         (gr_funcptr) deccfloat_nan},
    {GR_METHOD_FLOOR,           (gr_funcptr) deccfloat_floor},
    {GR_METHOD_CEIL,            (gr_funcptr) deccfloat_ceil},
    {GR_METHOD_TRUNC,           (gr_funcptr) deccfloat_trunc},
    {GR_METHOD_NINT,            (gr_funcptr) deccfloat_nint},
    {GR_METHOD_I,               (gr_funcptr) deccfloat_i},
    {GR_METHOD_ABS,             (gr_funcptr) deccfloat_abs},
    {GR_METHOD_CONJ,            (gr_funcptr) deccfloat_conj},
    {GR_METHOD_RE,              (gr_funcptr) deccfloat_re},
    {GR_METHOD_IM,              (gr_funcptr) deccfloat_im},
    {GR_METHOD_SGN,             (gr_funcptr) deccfloat_sgn},
    {GR_METHOD_CSGN,            (gr_funcptr) deccfloat_csgn},
    {GR_METHOD_ARG,             (gr_funcptr) deccfloat_arg},
    {GR_METHOD_CMP,             (gr_funcptr) deccfloat_cmp},
    {GR_METHOD_CMPABS,          (gr_funcptr) deccfloat_cmpabs},
    {GR_METHOD_VEC_DOT,         (gr_funcptr) deccfloat_vec_dot},
    {GR_METHOD_VEC_DOT_REV,     (gr_funcptr) deccfloat_vec_dot_rev},
    {GR_METHOD_MAT_DET,         (gr_funcptr) gr_mat_det_generic_field},
    {GR_METHOD_MAT_FIND_NONZERO_PIVOT, (gr_funcptr) gr_mat_find_nonzero_pivot_large_abs},
    {0,                         (gr_funcptr) NULL},
};

POP_OPTIONS
