/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <string.h>
#include "qqbar.h"
#include "decimal.h"
#include "mag.h"
#include "fmpq.h"
#include "fmpz_extras.h"
#include "arf.h"
#include "arb.h"
#include "acb.h"
#include "gr.h"
#include "gr_generic.h"
#include "gr_mat.h"

/* sign tests: lo = mid - rad, hi = mid + rad */
static void
_decball_bounds(decfloat_t lo, decfloat_t hi, const decball_t x, gr_ctx_t ctx)
{
    decfloat_t r;
    decfloat_init(r, ctx);
    GR_MUST_SUCCEED(_decmag_get_decfloat(r, &x->rad, ctx));
    /* rounding toward the outside keeps bounds valid; we use a modest precision */
    GR_MUST_SUCCEED(_decfloat_add(lo, &x->mid, r, 1, DECIMAL_CTX_RAD_PREC(ctx) + 3, DECIMAL_RND_FLOOR | DECIMAL_RND_NOLIMITS, NULL, NULL, ctx));
    GR_MUST_SUCCEED(_decfloat_add(hi, &x->mid, r, 0, DECIMAL_CTX_RAD_PREC(ctx) + 3, DECIMAL_RND_CEIL | DECIMAL_RND_NOLIMITS, NULL, NULL, ctx));
    decfloat_clear(r, ctx);
}

/* ------------------------------------------------------------------------- */
/*    Constants and assignment                                               */
/* ------------------------------------------------------------------------- */

int
decball_zero(decball_t res, gr_ctx_t ctx)
{
    _decmag_zero(&res->rad, ctx);
    return decfloat_zero(&res->mid, ctx);
}

int
decball_one(decball_t res, gr_ctx_t ctx)
{
    _decmag_zero(&res->rad, ctx);
    return decfloat_one(&res->mid, ctx);
}

int
decball_neg_one(decball_t res, gr_ctx_t ctx)
{
    _decmag_zero(&res->rad, ctx);
    return decfloat_neg_one(&res->mid, ctx);
}

int
decball_zero_pm_inf(decball_t res, gr_ctx_t ctx)
{
    decfloat_zero(&res->mid, ctx);
    _decmag_inf(&res->rad, ctx);
    return GR_SUCCESS;
}

/* Add the rounding error of the midpoint operation to rad. */
void
_decball_add_rounding_error(decball_t res, const decimal_rounding_info * info, const decmag_t precise_err, slong prec, gr_ctx_t ctx)
{
    /* exponent underflow: the true value lies in [0 +/- 10^emin] */
    if (info->underflow)
    {
        decmag_t t;
        _decmag_init(t, ctx);
        _decmag_set_10exp_si(t, DECIMAL_CTX_EMIN(ctx), ctx);
        _decmag_add(&res->rad, &res->rad, t, ctx);
        _decmag_clear(t, ctx);
        return;
    }

    /* exponent overflow: unknown magnitude */
    if (info->overflow)
    {
        _decmag_inf(&res->rad, ctx);
        return;
    }

    if (!info->inexact)
        return;

    if (DECIMAL_CTX_PRECISE_RADIUS(ctx))
    {
        _decmag_add(&res->rad, &res->rad, precise_err, ctx);
    }
    else
    {
        decmag_t ulp;
        int rnd = DECIMAL_CTX_RND(ctx);
        _decmag_init(ulp, ctx);
        _decmag_set_ulp(ulp, &res->mid, prec, ctx);
        if ((rnd & DECIMAL_RND_MASK) >= DECIMAL_RND_NEAR)
        {
            /* half an ulp: 5 * 10^(E - prec) */
            fmpz_t t;
            fmpz_init(t);
            decfloat_get_sci_exp(t, &res->mid, ctx);
            fmpz_sub_ui(t, t, prec);
            _decmag_set_ui_10exp_fmpz(ulp, 5, t, ctx);
            fmpz_clear(t);
        }
        _decmag_add(&res->rad, &res->rad, ulp, ctx);
        _decmag_clear(ulp, ctx);
    }
}

/* Rounds the midpoint of a ball to prec digits (context rounding mode),
   adding the rounding error to the radius. */
int
decball_set_round(decball_t res, const decball_t x, slong prec, gr_ctx_t ctx)
{
    decimal_rounding_info info;
    decmag_t err;
    int status;

    _decmag_init(err, ctx);
    status = decfloat_set_round_info(&res->mid, &x->mid, prec, DECIMAL_CTX_RND(ctx),
        &info, DECIMAL_CTX_PRECISE_RADIUS(ctx) ? err : NULL, ctx);

    if (status == GR_SUCCESS)
    {
        _decmag_set(&res->rad, &x->rad, ctx);
        _decball_add_rounding_error(res, &info, err, prec, ctx);
    }

    _decmag_clear(err, ctx);
    return status;
}

int
decball_set(decball_t res, const decball_t x, gr_ctx_t ctx)
{
    return decball_set_round(res, x, DECIMAL_CTX_PREC(ctx), ctx);
}

int
decball_set_round2(decball_t res, const decball_t x, slong prec, slong rad_prec, gr_ctx_t ctx)
{
    int status = decball_set_round(res, x, prec, ctx);
    if (status == GR_SUCCESS)
        _decmag_set_round(&res->rad, &res->rad, rad_prec, ctx);
    return status;
}

int
decball_set_decfloat(decball_t res, const decfloat_t x, gr_ctx_t ctx)
{
    decimal_rounding_info info;
    decmag_t err;
    int status;

    _decmag_init(err, ctx);
    status = decfloat_set_round_info(&res->mid, x, DECIMAL_CTX_PREC(ctx), DECIMAL_CTX_RND(ctx),
        &info, DECIMAL_CTX_PRECISE_RADIUS(ctx) ? err : NULL, ctx);
    _decmag_zero(&res->rad, ctx);
    if (status == GR_SUCCESS)
        _decball_add_rounding_error(res, &info, err, DECIMAL_CTX_PREC(ctx), ctx);
    _decmag_clear(err, ctx);
    return status;
}

DECIMAL_DRIVER int
_decball_set_scalar(decball_t res, const void * y, int type, gr_ctx_t ctx)
{
    decfloat_t t;
    int status;
    decfloat_init(t, ctx);
    status = _decfloat_set_scalar_exact(t, y, type, ctx);
    if (status == GR_SUCCESS)
        status = decball_set_decfloat(res, t, ctx);
    decfloat_clear(t, ctx);
    return status;
}

int decball_set_si(decball_t res, slong x, gr_ctx_t ctx) { return _decball_set_scalar(res, &x, DECIMAL_SCALAR_SI, ctx); }
int decball_set_ui(decball_t res, ulong x, gr_ctx_t ctx) { return _decball_set_scalar(res, &x, DECIMAL_SCALAR_UI, ctx); }
int decball_set_fmpz(decball_t res, const fmpz_t x, gr_ctx_t ctx) { return _decball_set_scalar(res, x, DECIMAL_SCALAR_FMPZ, ctx); }
int decball_set_d(decball_t res, double x, gr_ctx_t ctx) { return _decball_set_scalar(res, &x, DECIMAL_SCALAR_D, ctx); }

/* ------------------------------------------------------------------------- */
/*    Arithmetic                                                             */
/* ------------------------------------------------------------------------- */

int
decball_neg(decball_t res, const decball_t x, gr_ctx_t ctx)
{
    decimal_rounding_info info;
    decmag_t err;
    int status;

    _decmag_init(err, ctx);
    status = decfloat_set_round_info(&res->mid, &x->mid, DECIMAL_CTX_PREC(ctx),
        DECIMAL_RND_NEGATE(DECIMAL_CTX_RND(ctx)),
        &info, DECIMAL_CTX_PRECISE_RADIUS(ctx) ? err : NULL, ctx);

    if (status == GR_SUCCESS)
    {
        res->mid.m.size = -res->mid.m.size;

        _decmag_set(&res->rad, &x->rad, ctx);
        _decball_add_rounding_error(res, &info, err, DECIMAL_CTX_PREC(ctx), ctx);
    }

    _decmag_clear(err, ctx);
    return status;
}

int
decball_abs(decball_t res, const decball_t x, gr_ctx_t ctx)
{
    if (_decfloat_sgn(&x->mid, ctx) < 0)
        return decball_neg(res, x, ctx);
    return decball_set(res, x, ctx);
}

int
decball_add_round(decball_t res, const decball_t x, const decball_t y, slong prec, gr_ctx_t ctx)
{
    decimal_rounding_info info;
    decmag_t err;
    int status;

    _decmag_init(err, ctx);
    status = _decfloat_add(&res->mid, &x->mid, &y->mid, 0, prec, DECIMAL_CTX_RND(ctx),
        &info, DECIMAL_CTX_PRECISE_RADIUS(ctx) ? err : NULL, ctx);

    if (status == GR_SUCCESS)
    {
        _decmag_add(&res->rad, &x->rad, &y->rad, ctx);
        _decball_add_rounding_error(res, &info, err, prec, ctx);
    }

    _decmag_clear(err, ctx);
    return status;
}

int
decball_sub_round(decball_t res, const decball_t x, const decball_t y, slong prec, gr_ctx_t ctx)
{
    decimal_rounding_info info;
    decmag_t err;
    int status;

    _decmag_init(err, ctx);
    status = _decfloat_add(&res->mid, &x->mid, &y->mid, 1, prec, DECIMAL_CTX_RND(ctx),
        &info, DECIMAL_CTX_PRECISE_RADIUS(ctx) ? err : NULL, ctx);

    if (status == GR_SUCCESS)
    {
        _decmag_add(&res->rad, &x->rad, &y->rad, ctx);
        _decball_add_rounding_error(res, &info, err, prec, ctx);
    }

    _decmag_clear(err, ctx);
    return status;
}

int
decball_mul_round(decball_t res, const decball_t x, const decball_t y, slong prec, gr_ctx_t ctx)
{
    decimal_rounding_info info;
    decmag_t err, rad, t, u;
    int status;

    _decmag_init(err, ctx);
    _decmag_init(rad, ctx);
    _decmag_init(t, ctx);
    _decmag_init(u, ctx);

    /* radius first (midpoints are read before being overwritten) */
    if (DECMAG_IS_ZERO(&x->rad) && DECMAG_IS_ZERO(&y->rad))
    {
        _decmag_zero(rad, ctx);
    }
    else
    {
        _decmag_set_decfloat(t, &x->mid, ctx);
        _decmag_set_decfloat(u, &y->mid, ctx);
        _decmag_mul(rad, t, &y->rad, ctx);          /* |xm| yr */
        _decmag_addmul(rad, u, &x->rad, ctx);       /* |ym| xr */
        _decmag_addmul(rad, &x->rad, &y->rad, ctx); /* xr yr */
    }

    status = _decfloat_mul(&res->mid, &x->mid, &y->mid, prec, DECIMAL_CTX_RND(ctx),
        &info, DECIMAL_CTX_PRECISE_RADIUS(ctx) ? err : NULL, ctx);

    if (status == GR_SUCCESS)
    {
        _decmag_swap(&res->rad, rad, ctx);
        _decball_add_rounding_error(res, &info, err, prec, ctx);
    }

    _decmag_clear(err, ctx);
    _decmag_clear(rad, ctx);
    _decmag_clear(t, ctx);
    _decmag_clear(u, ctx);
    return status;
}

int
decball_sqr(decball_t res, const decball_t x, gr_ctx_t ctx)
{
    return decball_mul_round(res, x, x, DECIMAL_CTX_PREC(ctx), ctx);
}

int
decball_div_round(decball_t res, const decball_t x, const decball_t y, slong prec, gr_ctx_t ctx)
{
    decimal_rounding_info info;
    decmag_t err, rad, t, u, v;
    int status;

    if (DECFLOAT_IS_ZERO(&y->mid) && DECMAG_IS_ZERO(&y->rad))
        return GR_DOMAIN;

    if (DECMAG_IS_INF(&y->rad))
        return GR_UNABLE;

    if (_decball_contains_zero(y, ctx))
        return GR_UNABLE;

    _decmag_init(err, ctx);
    _decmag_init(rad, ctx);
    _decmag_init(t, ctx);
    _decmag_init(u, ctx);
    _decmag_init(v, ctx);

    if (DECMAG_IS_ZERO(&x->rad) && DECMAG_IS_ZERO(&y->rad))
    {
        _decmag_zero(rad, ctx);
    }
    else
    {
        /* (|xm| yr + |ym| xr) / (|ym| (|ym| - yr)) */
        _decmag_set_decfloat(t, &x->mid, ctx);
        _decmag_set_decfloat(u, &y->mid, ctx);
        _decmag_mul(rad, t, &y->rad, ctx);
        _decmag_addmul(rad, u, &x->rad, ctx);

        _decmag_set_decfloat_lower(u, &y->mid, ctx);
        _decmag_sub_lower(v, u, &y->rad, ctx);
        _decmag_mul_lower(v, v, u, ctx);
        _decmag_div(rad, rad, v, ctx);
    }

    status = _decfloat_div(&res->mid, &x->mid, &y->mid, prec, DECIMAL_CTX_RND(ctx),
        &info, DECIMAL_CTX_PRECISE_RADIUS(ctx) ? err : NULL, ctx);

    if (status == GR_SUCCESS)
    {
        _decmag_swap(&res->rad, rad, ctx);
        _decball_add_rounding_error(res, &info, err, prec, ctx);
    }

    _decmag_clear(err, ctx);
    _decmag_clear(rad, ctx);
    _decmag_clear(t, ctx);
    _decmag_clear(u, ctx);
    _decmag_clear(v, ctx);
    return status;
}

int
decball_inv_round(decball_t res, const decball_t x, slong prec, gr_ctx_t ctx)
{
    decball_t one;
    int status;
    decball_init(one, ctx);
    GR_MUST_SUCCEED(decball_one(one, ctx));
    status = decball_div_round(res, one, x, prec, ctx);
    decball_clear(one, ctx);
    return status;
}

int
decball_sqrt_round(decball_t res, const decball_t x, slong prec, gr_ctx_t ctx)
{
    decimal_rounding_info info;
    decmag_t err, rad, lo, t;
    int status;

    if (DECMAG_IS_INF(&x->rad))
        return GR_UNABLE;

    if (_decball_is_negative(x, ctx))
        return GR_DOMAIN;

    if (!_decball_is_nonnegative(x, ctx))
        return GR_UNABLE;

    _decmag_init(err, ctx);
    _decmag_init(rad, ctx);
    _decmag_init(lo, ctx);
    _decmag_init(t, ctx);

    if (DECMAG_IS_ZERO(&x->rad))
    {
        _decmag_zero(rad, ctx);
    }
    else
    {
        _decmag_set_decfloat_lower(lo, &x->mid, ctx);
        _decmag_sub_lower(lo, lo, &x->rad, ctx);

        if (DECMAG_IS_ZERO(lo))
        {
            /* |sqrt(t) - sqrt(m)| <= sqrt(m + r) */
            _decmag_set_decfloat(t, &x->mid, ctx);
            _decmag_add(t, t, &x->rad, ctx);
            _decmag_sqrt(rad, t, ctx);
        }
        else
        {
            /* r / (2 sqrt(lo)) */
            _decmag_sqrt_lower(t, lo, ctx);
            _decmag_mul_ui(t, t, 2, ctx);
            /* mul_ui rounds up; we need a lower bound of 2 sqrt(lo): recompute */
            _decmag_sqrt_lower(t, lo, ctx);
            {
                decmag_t two;
                _decmag_init(two, ctx);
                _decmag_set_ui_lower(two, 2, ctx);
                _decmag_mul_lower(t, t, two, ctx);
                _decmag_clear(two, ctx);
            }
            _decmag_div(rad, &x->rad, t, ctx);
        }
    }

    status = _decfloat_sqrt(&res->mid, &x->mid, prec, DECIMAL_CTX_RND(ctx),
        &info, DECIMAL_CTX_PRECISE_RADIUS(ctx) ? err : NULL, ctx);

    if (status == GR_SUCCESS)
    {
        _decmag_swap(&res->rad, rad, ctx);
        _decball_add_rounding_error(res, &info, err, prec, ctx);
    }

    _decmag_clear(err, ctx);
    _decmag_clear(rad, ctx);
    _decmag_clear(lo, ctx);
    _decmag_clear(t, ctx);
    return status;
}

int
decball_rsqrt(decball_t res, const decball_t x, gr_ctx_t ctx)
{
    decball_t t;
    int status;

    if (DECMAG_IS_INF(&x->rad))
        return GR_UNABLE;
    if (_decball_is_nonpositive(x, ctx))
        return GR_DOMAIN;
    if (!_decball_is_positive(x, ctx))
        return GR_UNABLE;

    decball_init(t, ctx);
    status = decball_sqrt_round(t, x, DECIMAL_CTX_PREC(ctx) + 5, ctx);
    if (status == GR_SUCCESS)
        status = decball_inv_round(res, t, DECIMAL_CTX_PREC(ctx), ctx);
    decball_clear(t, ctx);
    return status;
}

int decball_add(decball_t res, const decball_t x, const decball_t y, gr_ctx_t ctx) { return decball_add_round(res, x, y, DECIMAL_CTX_PREC(ctx), ctx); }
int decball_sub(decball_t res, const decball_t x, const decball_t y, gr_ctx_t ctx) { return decball_sub_round(res, x, y, DECIMAL_CTX_PREC(ctx), ctx); }
int decball_mul(decball_t res, const decball_t x, const decball_t y, gr_ctx_t ctx) { return decball_mul_round(res, x, y, DECIMAL_CTX_PREC(ctx), ctx); }
int decball_div(decball_t res, const decball_t x, const decball_t y, gr_ctx_t ctx) { return decball_div_round(res, x, y, DECIMAL_CTX_PREC(ctx), ctx); }
int decball_inv(decball_t res, const decball_t x, gr_ctx_t ctx) { return decball_inv_round(res, x, DECIMAL_CTX_PREC(ctx), ctx); }
int decball_sqrt(decball_t res, const decball_t x, gr_ctx_t ctx) { return decball_sqrt_round(res, x, DECIMAL_CTX_PREC(ctx), ctx); }

typedef int (*decball_round_op)(decball_t, const decball_t, const decball_t, slong, gr_ctx_t);

DECIMAL_DRIVER int
_decball_scalar_op(decball_t res, const decball_t x, const void * y, int type, decball_round_op op, gr_ctx_t ctx)
{
    decball_t t;
    int status;
    decball_init(t, ctx);
    status = _decfloat_set_scalar_exact(&t->mid, y, type, ctx);
    if (status == GR_SUCCESS)
        status = op(res, x, t, DECIMAL_CTX_PREC(ctx), ctx);
    decball_clear(t, ctx);
    return status;
}

int decball_add_ui(decball_t res, const decball_t x, ulong y, gr_ctx_t ctx) { return _decball_scalar_op(res, x, &y, DECIMAL_SCALAR_UI, decball_add_round, ctx); }
int decball_add_si(decball_t res, const decball_t x, slong y, gr_ctx_t ctx) { return _decball_scalar_op(res, x, &y, DECIMAL_SCALAR_SI, decball_add_round, ctx); }
int decball_add_fmpz(decball_t res, const decball_t x, const fmpz_t y, gr_ctx_t ctx) { return _decball_scalar_op(res, x, y, DECIMAL_SCALAR_FMPZ, decball_add_round, ctx); }
int decball_sub_ui(decball_t res, const decball_t x, ulong y, gr_ctx_t ctx) { return _decball_scalar_op(res, x, &y, DECIMAL_SCALAR_UI, decball_sub_round, ctx); }
int decball_sub_si(decball_t res, const decball_t x, slong y, gr_ctx_t ctx) { return _decball_scalar_op(res, x, &y, DECIMAL_SCALAR_SI, decball_sub_round, ctx); }
int decball_sub_fmpz(decball_t res, const decball_t x, const fmpz_t y, gr_ctx_t ctx) { return _decball_scalar_op(res, x, y, DECIMAL_SCALAR_FMPZ, decball_sub_round, ctx); }
int decball_mul_ui(decball_t res, const decball_t x, ulong y, gr_ctx_t ctx) { return _decball_scalar_op(res, x, &y, DECIMAL_SCALAR_UI, decball_mul_round, ctx); }
int decball_mul_si(decball_t res, const decball_t x, slong y, gr_ctx_t ctx) { return _decball_scalar_op(res, x, &y, DECIMAL_SCALAR_SI, decball_mul_round, ctx); }
int decball_mul_fmpz(decball_t res, const decball_t x, const fmpz_t y, gr_ctx_t ctx) { return _decball_scalar_op(res, x, y, DECIMAL_SCALAR_FMPZ, decball_mul_round, ctx); }
int decball_div_ui(decball_t res, const decball_t x, ulong y, gr_ctx_t ctx) { return _decball_scalar_op(res, x, &y, DECIMAL_SCALAR_UI, decball_div_round, ctx); }
int decball_div_si(decball_t res, const decball_t x, slong y, gr_ctx_t ctx) { return _decball_scalar_op(res, x, &y, DECIMAL_SCALAR_SI, decball_div_round, ctx); }
int decball_div_fmpz(decball_t res, const decball_t x, const fmpz_t y, gr_ctx_t ctx) { return _decball_scalar_op(res, x, y, DECIMAL_SCALAR_FMPZ, decball_div_round, ctx); }

int
decball_mul_two(decball_t res, const decball_t x, gr_ctx_t ctx)
{
    return decball_add_round(res, x, x, DECIMAL_CTX_PREC(ctx), ctx);
}

int
decball_mul_10exp_fmpz(decball_t res, const decball_t x, const fmpz_t e, gr_ctx_t ctx)
{
    decimal_rounding_info info;
    decmag_t err;
    int status;

    _decmag_init(err, ctx);
    status = _decfloat_mul_10exp(&res->mid, &x->mid, e, DECIMAL_CTX_PREC(ctx), DECIMAL_CTX_RND(ctx),
        &info, DECIMAL_CTX_PRECISE_RADIUS(ctx) ? err : NULL, ctx);
    if (status == GR_SUCCESS)
    {
        _decmag_mul_10exp_fmpz(&res->rad, &x->rad, e, ctx);
        _decball_add_rounding_error(res, &info, err, DECIMAL_CTX_PREC(ctx), ctx);
    }
    _decmag_clear(err, ctx);
    return status;
}

int
decball_mul_10exp_si(decball_t res, const decball_t x, slong e, gr_ctx_t ctx)
{
    fmpz_t t;
    int status;
    fmpz_init_set_si(t, e);
    status = decball_mul_10exp_fmpz(res, x, t, ctx);
    fmpz_clear(t);
    return status;
}

int
decball_mul_2exp_fmpz(decball_t res, const decball_t x, const fmpz_t e, gr_ctx_t ctx)
{
    decball_t p;
    fmpz_t one;
    int status;

    decball_init(p, ctx);
    fmpz_init_set_ui(one, 1);
    _decmag_zero(&p->rad, ctx);
    if (fmpz_bits(e) <= 16)
    {
        status = decfloat_set_round_fmpz_2exp_fmpz(&p->mid, one, e, DECIMAL_PREC_EXACT, DECIMAL_RND_DOWN | DECIMAL_RND_NOLIMITS, ctx);
    }
    else
    {
        /* an enclosure of 2^e (the exact power has too many digits),
           computed without exponent limits */
        gr_ctx_t bctx;
        arb_t t;
        _gr_ctx_init_decimal(bctx, DECIMAL_CTX_BALL, DECIMAL_CTX_E(ctx), DECIMAL_CTX_PREC(ctx) + 5, DECIMAL_RND_DOWN, 0);
        decimal_ctx_set_rad_prec(bctx, DECIMAL_CTX_RAD_PREC(ctx));
        arb_init(t);
        arb_one(t);
        arb_mul_2exp_fmpz(t, t, e);
        status = decball_set_arb(p, t, bctx);
        arb_clear(t);
        gr_ctx_clear(bctx);
    }
    if (status == GR_SUCCESS)
        status = decball_mul(res, x, p, ctx);
    decball_clear(p, ctx);
    fmpz_clear(one);
    return status;
}

int
decball_mul_2exp_si(decball_t res, const decball_t x, slong e, gr_ctx_t ctx)
{
    fmpz_t t;
    int status;
    fmpz_init_set_si(t, e);
    status = decball_mul_2exp_fmpz(res, x, t, ctx);
    fmpz_clear(t);
    return status;
}

/* The remaining functions are not performance-critical. */
PUSH_OPTIONS
OPTIMIZE_OSIZE

/* ------------------------------------------------------------------------- */
/*    Conversion, components and bounds                                      */
/* ------------------------------------------------------------------------- */

int
decball_set_fmpz_10exp_fmpz(decball_t res, const fmpz_t m, const fmpz_t e, gr_ctx_t ctx)
{
    decfloat_t t;
    int status;
    decfloat_init(t, ctx);
    status = decfloat_set_round_fmpz_10exp_fmpz(t, m, e, DECIMAL_PREC_EXACT, DECIMAL_RND_DOWN | DECIMAL_RND_NOLIMITS, ctx);
    if (status == GR_SUCCESS)
        status = decball_set_decfloat(res, t, ctx);
    decfloat_clear(t, ctx);
    return status;
}

int
decball_set_fmpq(decball_t res, const fmpq_t x, gr_ctx_t ctx)
{
    decball_t a, b;
    int status;

    if (fmpz_is_one(fmpq_denref(x)))
        return decball_set_fmpz(res, fmpq_numref(x), ctx);

    /* exact operands, which may be outside the exponent range */
    decball_init(a, ctx);
    decball_init(b, ctx);
    status = decfloat_set_round_fmpz(&a->mid, fmpq_numref(x), DECIMAL_PREC_EXACT, DECIMAL_RND_DOWN | DECIMAL_RND_NOLIMITS, ctx);
    status |= decfloat_set_round_fmpz(&b->mid, fmpq_denref(x), DECIMAL_PREC_EXACT, DECIMAL_RND_DOWN | DECIMAL_RND_NOLIMITS, ctx);
    if (status == GR_SUCCESS)
        status = decball_div(res, a, b, ctx);
    decball_clear(a, ctx);
    decball_clear(b, ctx);
    return status;
}

int
decball_set_interval_mid_rad(decball_t res, const decball_t m, const decball_t r, gr_ctx_t ctx)
{
    decmag_t t;
    int status;

    _decmag_init(t, ctx);
    decball_get_abs_ubound(t, r, ctx);
    status = decball_set(res, m, ctx);
    _decmag_add(&res->rad, &res->rad, t, ctx);
    _decmag_clear(t, ctx);
    return status;
}

int
decball_set_str(decball_t res, const char * s, gr_ctx_t ctx)
{
    decfloat_t t;
    int status;

    /* fast path for plain literals; everything else (including
       "mid +/- rad" and "[mid +/- rad]") goes through the generic
       expression parser, which calls decball_set_interval_mid_rad */
    decfloat_init(t, ctx);
    status = _decfloat_set_str_literal(t, s, DECIMAL_PREC_EXACT, DECIMAL_RND_DOWN | DECIMAL_RND_NOLIMITS, ctx);

    if (status == GR_SUCCESS)
        status = decball_set_decfloat(res, t, ctx);
    else if (status == GR_DOMAIN)
        status = gr_generic_set_str_ring_exponents(res, s, ctx);

    decfloat_clear(t, ctx);
    return status;
}

/* conversion from a decfloat in another decimal context */
int
_decball_set_decfloat_other(decball_t res, const decfloat_t x, gr_ctx_t x_ctx, gr_ctx_t ctx)
{
    decfloat_t t;
    int status;
    decfloat_init(t, ctx);
    /* exact conversion through the float ring with exact precision */
    {
        slong saved = DECIMAL_CTX_PREC(ctx);
        DECIMAL_CTX_PREC(ctx) = DECIMAL_PREC_EXACT;
        status = _decfloat_set_decfloat_other(t, x, x_ctx, ctx);
        DECIMAL_CTX_PREC(ctx) = saved;
    }
    if (status == GR_SUCCESS)
        status = decball_set_decfloat(res, t, ctx);
    decfloat_clear(t, ctx);
    return status;
}

/* conversion from a decball in another decimal context */
int
_decball_set_decball_other(decball_t res, const decball_t y, gr_ctx_t x_ctx, gr_ctx_t ctx)
{
    int status;

    if (DECIMAL_CTX_E(x_ctx) == DECIMAL_CTX_E(ctx))
        return decball_set(res, y, ctx);

    status = _decball_set_decfloat_other(res, &y->mid, x_ctx, ctx);

    if (status == GR_SUCCESS && !DECMAG_IS_ZERO(&y->rad))
        _decmag_add(&res->rad, &res->rad, &y->rad, ctx);

    return status;
}

int
decball_set_other(decball_t res, gr_srcptr x, gr_ctx_t x_ctx, gr_ctx_t ctx)
{
    switch (x_ctx->which_ring)
    {
        case GR_CTX_FMPZ:
            return decball_set_fmpz(res, x, ctx);

        case GR_CTX_FMPQ:
            return decball_set_fmpq(res, x, ctx);

        case GR_CTX_REAL_FLOAT_ARF:
            {
                decfloat_t t;
                int status;
                decfloat_init(t, ctx);
                status = decfloat_set_round_arf(t, x, DECIMAL_PREC_EXACT, DECIMAL_RND_DOWN | DECIMAL_RND_NOLIMITS, ctx);
                if (status == GR_SUCCESS)
                    status = decball_set_decfloat(res, t, ctx);
                decfloat_clear(t, ctx);
                return status;
            }

        case GR_CTX_RR_ARB:
            return decball_set_arb(res, x, ctx);

        case GR_CTX_DECFLOAT:
            return _decball_set_decfloat_other(res, x, x_ctx, ctx);

        case GR_CTX_DECBALL:
            return _decball_set_decball_other(res, x, x_ctx, ctx);

        case GR_CTX_DECCFLOAT:
            if (!_deccfloat_is_real((deccfloat_srcptr) x))
                return _deccfloat_is_nan((deccfloat_srcptr) x) ? GR_UNABLE : GR_DOMAIN;
            return _decball_set_decfloat_other(res, DECCFLOAT_REALREF((deccfloat_srcptr) x), x_ctx, ctx);

        case GR_CTX_DECCBALL:
            if (!_deccball_is_real((deccball_srcptr) x, x_ctx))
                return _decball_contains_zero(DECCBALL_IMAGREF((deccball_srcptr) x), x_ctx) ? GR_UNABLE : GR_DOMAIN;
            return _decball_set_decball_other(res, DECCBALL_REALREF((deccball_srcptr) x), x_ctx, ctx);

        case GR_CTX_REAL_ALGEBRAIC_QQBAR:
        case GR_CTX_COMPLEX_ALGEBRAIC_QQBAR:
            return decball_set_qqbar(res, x, ctx);

        default:
            {
                gr_ctx_t cctx;
                acb_t z;
                int status;

                gr_ctx_init_complex_acb(cctx, 20 + 4 * DECIMAL_CTX_PREC(ctx));
                acb_init(z);

                status = gr_set_other(z, x, x_ctx, cctx);

                if (status == GR_SUCCESS)
                {
                    if (acb_is_real(z))
                        status = decball_set_arb(res, acb_realref(z), ctx);
                    else
                        status = GR_DOMAIN;
                }

                acb_clear(z);
                gr_ctx_clear(cctx);

                return status;
            }
    }
}

int
decball_get_mid(decfloat_t res, const decball_t x, gr_ctx_t ctx)
{
    return decfloat_set_round(res, &x->mid, DECIMAL_PREC_EXACT, DECIMAL_RND_DOWN | DECIMAL_RND_NOLIMITS, ctx);
}

void
decball_get_rad(decmag_t res, const decball_t x, gr_ctx_t ctx)
{
    _decmag_set(res, &x->rad, ctx);
}

int
decball_mid(decball_t res, const decball_t x, gr_ctx_t ctx)
{
    _decmag_zero(&res->rad, ctx);
    return decball_get_mid(&res->mid, x, ctx);
}

int
decball_rad(decball_t res, const decball_t x, gr_ctx_t ctx)
{
    _decmag_zero(&res->rad, ctx);
    return _decmag_get_decfloat(&res->mid, &x->rad, ctx);
}

int
decball_shell(decball_t res, const decball_t x, gr_ctx_t ctx)
{
    _decmag_set(&res->rad, &x->rad, ctx);
    return decfloat_zero(&res->mid, ctx);
}

/* lower (upper) endpoint of x, rounded toward -inf (+inf) to prec digits;
   GR_DOMAIN for an infinite radius (the endpoint is not a real number) */
static int
_decball_endpoint(decfloat_t res, const decball_t x, int upper, slong prec, gr_ctx_t ctx)
{
    decfloat_t r;
    int status;

    if (DECMAG_IS_INF(&x->rad))
        return GR_DOMAIN;

    decfloat_init(r, ctx);
    GR_MUST_SUCCEED(_decmag_get_decfloat(r, &x->rad, ctx));
    status = _decfloat_add(res, &x->mid, r, !upper, prec,
        upper ? DECIMAL_RND_CEIL : DECIMAL_RND_FLOOR, NULL, NULL, ctx);
    decfloat_clear(r, ctx);
    return status;
}

int
decball_lower(decball_t res, const decball_t x, gr_ctx_t ctx)
{
    _decmag_zero(&res->rad, ctx);
    return _decball_endpoint(&res->mid, x, 0, DECIMAL_CTX_PREC(ctx), ctx);
}

int
decball_upper(decball_t res, const decball_t x, gr_ctx_t ctx)
{
    _decmag_zero(&res->rad, ctx);
    return _decball_endpoint(&res->mid, x, 1, DECIMAL_CTX_PREC(ctx), ctx);
}

int
decball_abs_upper(decball_t res, const decball_t x, gr_ctx_t ctx)
{
    decfloat_t r;
    int status;

    _decmag_zero(&res->rad, ctx);

    if (DECMAG_IS_INF(&x->rad))
        return GR_DOMAIN;

    decfloat_init(r, ctx);
    GR_MUST_SUCCEED(_decmag_get_decfloat(r, &x->rad, ctx));
    /* |mid| + r rounded up */
    status = _decfloat_add(&res->mid, &x->mid, r, DECFLOAT_SGNBIT(&x->mid), DECIMAL_CTX_PREC(ctx),
        DECIMAL_RND_UP, NULL, NULL, ctx);
    if (status == GR_SUCCESS)
        res->mid.m.size = FLINT_ABS(res->mid.m.size);
    decfloat_clear(r, ctx);
    return status;
}

int
decball_abs_lower(decball_t res, const decball_t x, gr_ctx_t ctx)
{
    decfloat_t r;
    int status;

    _decmag_zero(&res->rad, ctx);

    if (DECMAG_IS_INF(&x->rad))
        return decfloat_zero(&res->mid, ctx);

    decfloat_init(r, ctx);
    GR_MUST_SUCCEED(_decmag_get_decfloat(r, &x->rad, ctx));
    /* |mid| - r rounded toward zero, or zero if negative */
    status = _decfloat_add(&res->mid, &x->mid, r, !DECFLOAT_SGNBIT(&x->mid), DECIMAL_CTX_PREC(ctx),
        DECIMAL_RND_DOWN, NULL, NULL, ctx);
    if (status == GR_SUCCESS)
    {
        if (DECFLOAT_SGNBIT(&x->mid) != DECFLOAT_SGNBIT(&res->mid))
            decfloat_zero(&res->mid, ctx);
        else
            res->mid.m.size = FLINT_ABS(res->mid.m.size);
    }
    decfloat_clear(r, ctx);
    return status;
}

int
decball_add_rad(decball_t res, const decball_t x, const decball_t r, gr_ctx_t ctx)
{
    return decball_set_interval_mid_rad(res, x, r, ctx);
}

int
decball_set_interval_mid_inf(decball_t res, const decball_t m, gr_ctx_t ctx)
{
    int status = decball_set(res, m, ctx);
    _decmag_inf(&res->rad, ctx);
    return status;
}

/* the smallest ball containing lo and hi */
int
decball_set_interval(decball_t res, const decball_t lo, const decball_t hi, gr_ctx_t ctx)
{
    decfloat_t d, m, A, B;
    decmag_t r1, r2;
    slong wp;
    int status;

    if (DECMAG_IS_INF(&lo->rad) || DECMAG_IS_INF(&hi->rad))
        return decball_zero_pm_inf(res, ctx);

    decfloat_init(d, ctx);
    decfloat_init(m, ctx);
    decfloat_init(A, ctx);
    decfloat_init(B, ctx);
    _decmag_init(r1, ctx);
    _decmag_init(r2, ctx);

    /* endpoints A <= B of the hull, rounded outward at a precision high
       enough to make the widening negligible */
    wp = DECIMAL_CTX_IS_EXACT(ctx) ? DECIMAL_PREC_EXACT : DECIMAL_CTX_PREC(ctx) + DECMAG_MAX_PREC + 2;
    status = _decball_endpoint(A, lo, 0, wp, ctx);
    status |= _decball_endpoint(d, hi, 0, wp, ctx);
    if (status == GR_SUCCESS && _decfloat_cmp(d, A, ctx) < 0)
        decfloat_swap(A, d, ctx);
    status |= _decball_endpoint(B, hi, 1, wp, ctx);
    status |= _decball_endpoint(d, lo, 1, wp, ctx);
    if (status == GR_SUCCESS && _decfloat_cmp(d, B, ctx) > 0)
        decfloat_swap(B, d, ctx);

    /* midpoint (A + B) / 2 rounded to the context precision; the
       distances to the endpoints are then bounded from above */
    if (status == GR_SUCCESS)
        status = decfloat_add(m, A, B, ctx);
    if (status == GR_SUCCESS)
        status = decfloat_div_ui(m, m, 2, ctx);

    if (status == GR_SUCCESS)
    {
        GR_MUST_SUCCEED(_decfloat_add(d, m, A, 1, DECIMAL_CTX_RAD_PREC(ctx) + 3, DECIMAL_RND_UP | DECIMAL_RND_NOLIMITS, NULL, NULL, ctx));
        _decmag_set_decfloat(r1, d, ctx);
        GR_MUST_SUCCEED(_decfloat_add(d, B, m, 1, DECIMAL_CTX_RAD_PREC(ctx) + 3, DECIMAL_RND_UP | DECIMAL_RND_NOLIMITS, NULL, NULL, ctx));
        _decmag_set_decfloat(r2, d, ctx);
        _decmag_max(&res->rad, r1, r2, ctx);
        decfloat_swap(&res->mid, m, ctx);
    }

    decfloat_clear(d, ctx);
    decfloat_clear(m, ctx);
    decfloat_clear(A, ctx);
    decfloat_clear(B, ctx);
    _decmag_clear(r1, ctx);
    _decmag_clear(r2, ctx);
    return status;
}

void
decball_get_abs_ubound(decmag_t res, const decball_t x, gr_ctx_t ctx)
{
    decmag_t t;
    _decmag_init(t, ctx);
    _decmag_set_decfloat(t, &x->mid, ctx);
    _decmag_add(res, t, &x->rad, ctx);
    _decmag_clear(t, ctx);
}

void
decball_get_abs_lbound(decmag_t res, const decball_t x, gr_ctx_t ctx)
{
    decmag_t t;
    _decmag_init(t, ctx);
    _decmag_set_decfloat_lower(t, &x->mid, ctx);
    _decmag_sub_lower(res, t, &x->rad, ctx);
    _decmag_clear(t, ctx);
}

/* ------------------------------------------------------------------------- */
/*    Rounding to integers, powers and radius operations                     */
/* ------------------------------------------------------------------------- */

static int
_decball_round_to_int(decball_t res, const decball_t x, int int_rnd, gr_ctx_t ctx)
{
    int status;

    if (DECMAG_IS_ZERO(&x->rad))
    {
        status = _decfloat_round_to_int(&res->mid, &x->mid, int_rnd, DECIMAL_PREC_EXACT, DECIMAL_RND_DOWN, NULL, NULL, ctx);
        _decmag_zero(&res->rad, ctx);
        if (status == GR_SUCCESS)
            status = decball_set_round(res, res, DECIMAL_CTX_PREC(ctx), ctx);
        return status;
    }
    else if (DECMAG_IS_INF(&x->rad))
    {
        return decball_zero_pm_inf(res, ctx);
    }
    else
    {
        /* if the ball contains no integer, floor/ceil are constant on it */
        decfloat_t lo, hi, flo, fhi;
        decfloat_init(lo, ctx);
        decfloat_init(hi, ctx);
        decfloat_init(flo, ctx);
        decfloat_init(fhi, ctx);

        _decball_bounds(lo, hi, x, ctx);
        status = _decfloat_round_to_int(flo, lo, int_rnd, DECIMAL_PREC_EXACT, DECIMAL_RND_DOWN | DECIMAL_RND_NOLIMITS, NULL, NULL, ctx);
        status |= _decfloat_round_to_int(fhi, hi, int_rnd, DECIMAL_PREC_EXACT, DECIMAL_RND_DOWN | DECIMAL_RND_NOLIMITS, NULL, NULL, ctx);

        if (status == GR_SUCCESS && decfloat_equal(flo, fhi, ctx) == T_TRUE)
        {
            decfloat_swap(&res->mid, flo, ctx);
            _decmag_zero(&res->rad, ctx);
            status = decball_set_round(res, res, DECIMAL_CTX_PREC(ctx), ctx);
        }
        else if (status == GR_SUCCESS)
        {
            /* [round(mid) +/- (rad + 1)] */
            decmag_t r;
            _decmag_init(r, ctx);
            _decmag_set(r, &x->rad, ctx);
            status = _decfloat_round_to_int(&res->mid, &x->mid, int_rnd, DECIMAL_PREC_EXACT, DECIMAL_RND_DOWN, NULL, NULL, ctx);
            if (status == GR_SUCCESS)
            {
                decmag_t one;
                _decmag_init(one, ctx);
                _decmag_one(one, ctx);
                _decmag_add(&res->rad, r, one, ctx);
                _decmag_clear(one, ctx);
                status = decball_set_round(res, res, DECIMAL_CTX_PREC(ctx), ctx);
            }
            _decmag_clear(r, ctx);
        }

        decfloat_clear(lo, ctx);
        decfloat_clear(hi, ctx);
        decfloat_clear(flo, ctx);
        decfloat_clear(fhi, ctx);
        return status;
    }
}

int decball_floor(decball_t res, const decball_t x, gr_ctx_t ctx) { return _decball_round_to_int(res, x, DECIMAL_RND_FLOOR, ctx); }
int decball_ceil(decball_t res, const decball_t x, gr_ctx_t ctx) { return _decball_round_to_int(res, x, DECIMAL_RND_CEIL, ctx); }
int decball_trunc(decball_t res, const decball_t x, gr_ctx_t ctx) { return _decball_round_to_int(res, x, DECIMAL_RND_DOWN, ctx); }
int decball_nint(decball_t res, const decball_t x, gr_ctx_t ctx) { return _decball_round_to_int(res, x, DECIMAL_RND_NEAR, ctx); }

/* binary powering; does not depend on the method table of ctx, so that
   it can be called with a complex context */
int
_decball_pow_fmpz_binexp(decball_t res, const decball_t x, const fmpz_t y, gr_ctx_t ctx)
{
    decball_t t;
    slong i, bits;
    int status = GR_SUCCESS;

    if (fmpz_is_zero(y))
        return decball_one(res, ctx);

    if (fmpz_is_one(y))
        return decball_set(res, x, ctx);

    decball_init(t, ctx);
    status = decball_set(t, x, ctx);

    /* bits of |y| */
    {
        fmpz_t a;
        fmpz_init(a);
        fmpz_abs(a, y);
        bits = fmpz_bits(a);
        for (i = bits - 2; i >= 0 && status == GR_SUCCESS; i--)
        {
            status = decball_mul(t, t, t, ctx);
            if (status == GR_SUCCESS && fmpz_tstbit(a, i))
                status = decball_mul(t, t, x, ctx);
        }
        fmpz_clear(a);
    }

    if (status == GR_SUCCESS)
    {
        if (fmpz_sgn(y) < 0)
            status = decball_inv(res, t, ctx);
        else
            decball_swap(res, t, ctx);
    }

    decball_clear(t, ctx);
    return status;
}

int
decball_pow_fmpz(decball_t res, const decball_t x, const fmpz_t y, gr_ctx_t ctx)
{
    /* exact powers (in particular powers of ten) followed by a single rounding */
    if (DECMAG_IS_ZERO(&x->rad) && DECFLOAT_IS_FINITE(&x->mid))
    {
        decfloat_t t;
        int r;
        decfloat_init(t, ctx);
        r = _decfloat_pow_int_exact(t, &x->mid, y, DECIMAL_PREC_EXACT, DECIMAL_RND_DOWN | DECIMAL_RND_NOLIMITS, ctx);
        if (r == 1)
            r = decball_set_decfloat(res, t, ctx);
        else
            r = -1;
        decfloat_clear(t, ctx);
        if (r != -1)
            return r;
    }

    if (fmpz_bits(y) < 40)
        return _decball_pow_fmpz_binexp(res, x, y, ctx);

    {
        decball_t t;
        int status;
        decball_init(t, ctx);
        status = decball_set_fmpz(t, y, ctx);
        if (status == GR_SUCCESS)
            status = decball_pow(res, x, t, ctx);
        decball_clear(t, ctx);
        return status;
    }
}

int
decball_pow_ui(decball_t res, const decball_t x, ulong y, gr_ctx_t ctx)
{
    fmpz_t t;
    int status;
    fmpz_init_set_ui(t, y);
    status = decball_pow_fmpz(res, x, t, ctx);
    fmpz_clear(t);
    return status;
}

int
decball_pow_si(decball_t res, const decball_t x, slong y, gr_ctx_t ctx)
{
    fmpz_t t;
    int status;
    fmpz_init_set_si(t, y);
    status = decball_pow_fmpz(res, x, t, ctx);
    fmpz_clear(t);
    return status;
}

int
decball_add_error_decmag(decball_t res, const decmag_t err, gr_ctx_t ctx)
{
    _decmag_add(&res->rad, &res->rad, err, ctx);
    return GR_SUCCESS;
}

int
decball_add_error_10exp_si(decball_t res, slong e, gr_ctx_t ctx)
{
    decmag_t t;
    _decmag_init(t, ctx);
    _decmag_set_10exp_si(t, e, ctx);
    _decmag_add(&res->rad, &res->rad, t, ctx);
    _decmag_clear(t, ctx);
    return GR_SUCCESS;
}

int
decball_add_error_decfloat(decball_t res, const decfloat_t err, gr_ctx_t ctx)
{
    decmag_t t;
    _decmag_init(t, ctx);
    _decmag_set_decfloat(t, err, ctx);
    _decmag_add(&res->rad, &res->rad, t, ctx);
    _decmag_clear(t, ctx);
    return GR_SUCCESS;
}

slong
decball_rel_accuracy_digits(const decball_t x, gr_ctx_t ctx)
{
    fmpz_t Em, Er;
    slong res;

    if (DECMAG_IS_ZERO(&x->rad))
        return DECFLOAT_IS_FINITE(&x->mid) ? DECIMAL_PREC_EXACT : -DECIMAL_PREC_EXACT;

    if (!DECFLOAT_IS_FINITE(&x->mid) || DECMAG_IS_INF(&x->rad))
        return -DECIMAL_PREC_EXACT;

    fmpz_init(Em);
    fmpz_init(Er);

    _decmag_get_sci_exp(Er, &x->rad);

    if (!DECFLOAT_IS_ZERO(&x->mid))
        decfloat_get_sci_exp(Em, &x->mid, ctx);

    fmpz_sub(Em, Em, Er);

    if (fmpz_fits_si(Em))
    {
        res = fmpz_get_si(Em);
        res = FLINT_MIN(res, DECIMAL_PREC_EXACT);
        res = FLINT_MAX(res, -DECIMAL_PREC_EXACT);
    }
    else
        res = (fmpz_sgn(Em) > 0) ? DECIMAL_PREC_EXACT : -DECIMAL_PREC_EXACT;

    fmpz_clear(Em);
    fmpz_clear(Er);
    return res;
}

/* Round the midpoint to the number of digits justified by the radius. */
int
decball_trim(decball_t res, const decball_t x, gr_ctx_t ctx)
{
    slong Em, Er, prec;
    fmpz_t t;

    if (DECMAG_IS_ZERO(&x->rad) || DECFLOAT_IS_ZERO(&x->mid) || DECMAG_IS_INF(&x->rad))
        return decball_set(res, x, ctx);

    /* rad ~ 10^Er, mid ~ 10^Em: keep Em - Er + rad_prec + 1 digits */
    if (!decfloat_get_sci_exp_si(&Em, &x->mid, ctx))
        return decball_set(res, x, ctx);

    fmpz_init(t);
    _decmag_get_sci_exp(t, &x->rad);
    if (!fmpz_fits_si(t))
    {
        fmpz_clear(t);
        return decball_set(res, x, ctx);
    }
    Er = fmpz_get_si(t);
    fmpz_clear(t);

    prec = Em - Er + DECIMAL_CTX_RAD_PREC(ctx) + 1;
    if (prec < 1)
        prec = 1;
    if (prec >= DECIMAL_CTX_PREC(ctx))
        return decball_set(res, x, ctx);

    return decball_set_round(res, x, prec, ctx);
}

/* ------------------------------------------------------------------------- */
/*    Predicates                                                             */
/* ------------------------------------------------------------------------- */

/* |mid| <= rad ? */
int
_decball_contains_zero(const decball_t x, gr_ctx_t ctx)
{
    decfloat_t t;
    int c;

    if (DECFLOAT_IS_ZERO(&x->mid))
        return 1;
    if (DECMAG_IS_ZERO(&x->rad))
        return 0;
    if (DECMAG_IS_INF(&x->rad))
        return 1;

    decfloat_init(t, ctx);
    GR_MUST_SUCCEED(_decmag_get_decfloat(t, &x->rad, ctx));
    c = _decfloat_cmpabs(&x->mid, t, ctx);
    decfloat_clear(t, ctx);
    return c <= 0;
}

truth_t
decball_is_zero(const decball_t x, gr_ctx_t ctx)
{
    if (DECFLOAT_IS_ZERO(&x->mid) && DECMAG_IS_ZERO(&x->rad))
        return T_TRUE;
    if (_decball_contains_zero(x, ctx))
        return T_UNKNOWN;
    return T_FALSE;
}

int
_decball_contains_decfloat(const decball_t x, const decfloat_t y, gr_ctx_t ctx)
{
    decfloat_t t, r;
    int c, status;

    if (DECFLOAT_IS_NAN(y) || DECFLOAT_IS_INF(y))
        return 0;
    if (DECMAG_IS_INF(&x->rad))
        return 1;
    if (DECMAG_IS_ZERO(&x->rad))
        return decfloat_equal(&x->mid, y, ctx) == T_TRUE;

    /* |mid - y| <= rad, exactly */
    decfloat_init(t, ctx);
    decfloat_init(r, ctx);
    status = _decfloat_add(t, &x->mid, y, 1, DECIMAL_PREC_EXACT, DECIMAL_RND_DOWN | DECIMAL_RND_NOLIMITS, NULL, NULL, ctx);
    if (status != GR_SUCCESS)
    {
        /* exponents too far apart to subtract exactly: then the difference
           is dominated by the larger operand; compare its magnitude */
        const decfloat_struct * big = (_decfloat_cmpabs(&x->mid, y, ctx) >= 0) ? &x->mid : y;
        decmag_t bm;
        _decmag_init(bm, ctx);
        _decmag_set_decfloat_lower(bm, big, ctx);
        /* |mid - y| >= |big| - |small| > |big| / 2 roughly; use lower bound |big|/10 */
        _decmag_div_ui(bm, bm, 10, ctx);
        c = (_decmag_cmp(bm, &x->rad, ctx) <= 0);
        _decmag_clear(bm, ctx);
        decfloat_clear(t, ctx);
        decfloat_clear(r, ctx);
        return c;
    }
    GR_MUST_SUCCEED(_decmag_get_decfloat(r, &x->rad, ctx));
    c = _decfloat_cmpabs(t, r, ctx);
    decfloat_clear(t, ctx);
    decfloat_clear(r, ctx);
    return c <= 0;
}

int
_decball_contains_fmpz(const decball_t x, const fmpz_t y, gr_ctx_t ctx)
{
    decfloat_t t;
    int c;
    decfloat_init(t, ctx);
    GR_MUST_SUCCEED(decfloat_set_round_fmpz(t, y, DECIMAL_PREC_EXACT, DECIMAL_RND_DOWN | DECIMAL_RND_NOLIMITS, ctx));
    c = _decball_contains_decfloat(x, t, ctx);
    decfloat_clear(t, ctx);
    return c;
}

int
_decball_contains_si(const decball_t x, slong y, gr_ctx_t ctx)
{
    fmpz_t t;
    int c;
    fmpz_init_set_si(t, y);
    c = _decball_contains_fmpz(x, t, ctx);
    fmpz_clear(t);
    return c;
}

int
_decball_contains_fmpq(const decball_t x, const fmpq_t y, gr_ctx_t ctx)
{
    /* mid - rad <= y <= mid + rad, decided exactly with rationals */
    fmpq_t m, r, lo, hi;
    int result;

    if (fmpz_is_one(fmpq_denref(y)))
        return _decball_contains_fmpz(x, fmpq_numref(y), ctx);

    if (DECMAG_IS_INF(&x->rad))
        return 1;

    fmpq_init(m);
    fmpq_init(r);
    fmpq_init(lo);
    fmpq_init(hi);

    if (decfloat_get_fmpq(m, &x->mid, ctx) != GR_SUCCESS)
    {
        /* enormous exponent: fall back to a float comparison */
        decfloat_t t;
        decfloat_init(t, ctx);
        GR_MUST_SUCCEED(decfloat_set_round_fmpq(t, y, DECIMAL_PREC_EXACT, DECIMAL_RND_DOWN | DECIMAL_RND_NOLIMITS, ctx));
        result = _decball_contains_decfloat(x, t, ctx);
        decfloat_clear(t, ctx);
    }
    else
    {
        if (!fmpz_fits_si(&x->rad.exp))
        {
            result = fmpz_sgn(&x->rad.exp) > 0;   /* huge radius: contains; tiny: treat as not (mid is a decimal, y is not) */
            if (!result)
                result = fmpq_equal(m, y);
        }
        else
        {
            _decmag_get_fmpq(r, &x->rad, ctx);
            fmpq_sub(lo, m, r);
            fmpq_add(hi, m, r);
            result = (fmpq_cmp(lo, y) <= 0) && (fmpq_cmp(y, hi) <= 0);
        }
    }

    fmpq_clear(m);
    fmpq_clear(r);
    fmpq_clear(lo);
    fmpq_clear(hi);
    return result;
}

/* lower bound for |xm - ym| (exact difference rounded down at rp+3 digits),
   and upper bound in *hi */
static void
_decball_middist(decmag_t lo, decmag_t hi, const decball_t x, const decball_t y, gr_ctx_t ctx)
{
    decfloat_t t;
    decimal_rounding_info info;
    slong prec = DECIMAL_CTX_RAD_PREC(ctx) + 3;

    decfloat_init(t, ctx);
    GR_MUST_SUCCEED(_decfloat_add(t, &x->mid, &y->mid, 1, prec, DECIMAL_RND_DOWN | DECIMAL_RND_NOLIMITS, &info, NULL, ctx));
    _decmag_set_decfloat_lower(lo, t, ctx);
    _decmag_set_decfloat(hi, t, ctx);
    if (info.inexact)
    {
        decmag_t u;
        _decmag_init(u, ctx);
        _decmag_set_ulp(u, t, prec, ctx);
        _decmag_add(hi, hi, u, ctx);
        _decmag_clear(u, ctx);
    }
    decfloat_clear(t, ctx);
}

/* u = 10^E where E = sci_exp(x) - prec + 1, i.e. one unit in the last
   place of x at precision prec. x must be finite and nonzero. */
static void
_decfloat_ulp(decfloat_t u, const decfloat_t x, slong prec, gr_ctx_t ctx)
{
    fmpz_t E, one;
    fmpz_init(E);
    fmpz_init_set_ui(one, 1);
    decfloat_get_sci_exp(E, x, ctx);
    fmpz_sub_ui(E, E, prec - 1);
    GR_MUST_SUCCEED(decfloat_set_round_fmpz_10exp_fmpz(u, one, E, DECIMAL_PREC_EXACT, DECIMAL_RND_DOWN | DECIMAL_RND_NOLIMITS, ctx));
    fmpz_clear(E);
    fmpz_clear(one);
}

/*
    Decide exactly whether |xm - ym| + B <= C, where B and C are decimal
    floating-point numbers with few digits (radii; C >= 0). Returns 1 if
    true, 0 if false, and -1 if the decision would need an unreasonable
    amount of work (astronomically large exponent differences).

    We compute D = |xm - ym| and S = D + B with directed rounding at
    increasing precision p; inexact flags give the enclosures
    D <= |xm - ym| < D + ulp(D) and S <= D + B < S + ulp(S), which decide
    the comparison unless C lies in the narrow interval [S, S + 2 ulp).
    Since C has at most DECMAG_MAX_PREC < p digits, C is then S or S + ulp,
    which are decided by the exactness flags or by a higher precision.
*/
static int
_decball_decide_leq(const decfloat_t xm, const decfloat_t ym, const decfloat_t B, const decfloat_t C, gr_ctx_t ctx)
{
    decfloat_t D, S, U, u;
    decimal_rounding_info infoD, infoS;
    slong p;
    int result = -1;
    int status;

    decfloat_init(D, ctx);
    decfloat_init(S, ctx);
    decfloat_init(U, ctx);
    decfloat_init(u, ctx);

    for (p = DECMAG_MAX_PREC + 4; ; p *= 4)
    {
        if (p > (WORD(1) << 16))
        {
            /* last resort: exact computation */
            status = _decfloat_add(D, xm, ym, 1, DECIMAL_PREC_EXACT, DECIMAL_RND_DOWN | DECIMAL_RND_NOLIMITS, NULL, NULL, ctx);
            if (status == GR_SUCCESS)
            {
                D->m.size = FLINT_ABS(D->m.size);
                status = _decfloat_add(S, D, B, 0, DECIMAL_PREC_EXACT, DECIMAL_RND_DOWN | DECIMAL_RND_NOLIMITS, NULL, NULL, ctx);
            }
            if (status == GR_SUCCESS)
                result = (_decfloat_cmp(S, C, ctx) <= 0);
            break;
        }

        /* D <= |xm - ym| < D + ulp(D) */
        GR_MUST_SUCCEED(_decfloat_add(D, xm, ym, 1, p, DECIMAL_RND_DOWN | DECIMAL_RND_NOLIMITS, &infoD, NULL, ctx));
        D->m.size = FLINT_ABS(D->m.size);

        /* S <= D + B < S + ulp(S) */
        GR_MUST_SUCCEED(_decfloat_add(S, D, B, 0, p, DECIMAL_RND_FLOOR | DECIMAL_RND_NOLIMITS, &infoS, NULL, ctx));

        if (_decfloat_cmp(S, C, ctx) > 0)
        {
            result = 0;
            break;
        }

        /* upper bound U = S + ulp(S) + ulp(D) for the inexact terms */
        decfloat_set(U, S, ctx);
        if (infoS.inexact)
        {
            _decfloat_ulp(u, S, p, ctx);
            GR_MUST_SUCCEED(_decfloat_add(U, U, u, 0, DECIMAL_PREC_EXACT, DECIMAL_RND_DOWN | DECIMAL_RND_NOLIMITS, NULL, NULL, ctx));
        }
        if (infoD.inexact)
        {
            _decfloat_ulp(u, D, p, ctx);
            GR_MUST_SUCCEED(_decfloat_add(U, U, u, 0, DECIMAL_PREC_EXACT, DECIMAL_RND_DOWN | DECIMAL_RND_NOLIMITS, NULL, NULL, ctx));
        }

        if (_decfloat_cmp(U, C, ctx) <= 0)
        {
            result = 1;
            break;
        }

        /* now S <= C < U, so something was inexact and |xm - ym| + B > S */
        if (decfloat_equal(S, C, ctx) == T_TRUE)
        {
            result = 0;
            break;
        }

        if (!infoD.inexact)
        {
            /* D + B in (S, S + ulp(S)) and C in that interval too: impossible
               as C has fewer digits than p, so C >= S + ulp(S) = U */
            result = 1;
            break;
        }
    }

    decfloat_clear(D, ctx);
    decfloat_clear(S, ctx);
    decfloat_clear(U, ctx);
    decfloat_clear(u, ctx);

    return result;
}

int
_decball_overlaps(const decball_t x, const decball_t y, gr_ctx_t ctx)
{
    decmag_t dlo, dhi, slo, shi;
    int result;

    if (DECMAG_IS_INF(&x->rad) || DECMAG_IS_INF(&y->rad))
        return 1;

    if (decfloat_equal(&x->mid, &y->mid, ctx) == T_TRUE)
        return 1;

    if (DECMAG_IS_ZERO(&x->rad) && DECMAG_IS_ZERO(&y->rad))
        return 0;

    _decmag_init(dlo, ctx);
    _decmag_init(dhi, ctx);
    _decmag_init(slo, ctx);
    _decmag_init(shi, ctx);

    _decball_middist(dlo, dhi, x, y, ctx);
    _decmag_add(shi, &x->rad, &y->rad, ctx);
    _decmag_add_lower(slo, &x->rad, &y->rad, ctx);

    if (_decmag_cmp(dlo, shi, ctx) > 0)
        result = 0;
    else if (_decmag_cmp(dhi, slo, ctx) <= 0)
        result = 1;
    else
    {
        /* |xm - ym| <= xr + yr  <=>  |xm - ym| - yr <= xr */
        decfloat_t B, C;
        decfloat_init(B, ctx);
        decfloat_init(C, ctx);
        GR_MUST_SUCCEED(_decmag_get_decfloat(B, &y->rad, ctx));
        B->m.size = -B->m.size;
        GR_MUST_SUCCEED(_decmag_get_decfloat(C, &x->rad, ctx));
        result = _decball_decide_leq(&x->mid, &y->mid, B, C, ctx);
        if (result < 0)
            result = 1;   /* conservative */
        decfloat_clear(B, ctx);
        decfloat_clear(C, ctx);
    }

    _decmag_clear(dlo, ctx);
    _decmag_clear(dhi, ctx);
    _decmag_clear(slo, ctx);
    _decmag_clear(shi, ctx);
    return result;
}

int
_decball_contains(const decball_t x, const decball_t y, gr_ctx_t ctx)
{
    /* |xm - ym| + yr <= xr */
    decmag_t dlo, dhi, t;
    int result;

    if (DECMAG_IS_INF(&x->rad))
        return 1;
    if (DECMAG_IS_INF(&y->rad))
        return 0;

    if (DECMAG_IS_ZERO(&y->rad))
        return _decball_contains_decfloat(x, &y->mid, ctx);

    _decmag_init(dlo, ctx);
    _decmag_init(dhi, ctx);
    _decmag_init(t, ctx);

    _decball_middist(dlo, dhi, x, y, ctx);

    _decmag_add(t, dhi, &y->rad, ctx);
    if (_decmag_cmp(t, &x->rad, ctx) <= 0)
        result = 1;
    else
    {
        _decmag_add_lower(t, dlo, &y->rad, ctx);
        if (_decmag_cmp(t, &x->rad, ctx) > 0)
            result = 0;
        else
        {
            /* |xm - ym| + yr <= xr */
            decfloat_t B, C;
            decfloat_init(B, ctx);
            decfloat_init(C, ctx);
            GR_MUST_SUCCEED(_decmag_get_decfloat(B, &y->rad, ctx));
            GR_MUST_SUCCEED(_decmag_get_decfloat(C, &x->rad, ctx));
            result = _decball_decide_leq(&x->mid, &y->mid, B, C, ctx);
            if (result < 0)
                result = 0;   /* conservative */
            decfloat_clear(B, ctx);
            decfloat_clear(C, ctx);
        }
    }

    _decmag_clear(dlo, ctx);
    _decmag_clear(dhi, ctx);
    _decmag_clear(t, ctx);
    return result;
}

truth_t
decball_equal(const decball_t x, const decball_t y, gr_ctx_t ctx)
{
    if (DECMAG_IS_ZERO(&x->rad) && DECMAG_IS_ZERO(&y->rad))
        return decfloat_equal(&x->mid, &y->mid, ctx);

    if (_decball_overlaps(x, y, ctx))
        return T_UNKNOWN;

    return T_FALSE;
}

truth_t
decball_is_one(const decball_t x, gr_ctx_t ctx)
{
    if (DECMAG_IS_ZERO(&x->rad))
        return decfloat_is_one(&x->mid, ctx);
    if (_decball_contains_si(x, 1, ctx))
        return T_UNKNOWN;
    return T_FALSE;
}

truth_t
decball_is_neg_one(const decball_t x, gr_ctx_t ctx)
{
    if (DECMAG_IS_ZERO(&x->rad))
        return decfloat_is_neg_one(&x->mid, ctx);
    if (_decball_contains_si(x, -1, ctx))
        return T_UNKNOWN;
    return T_FALSE;
}

int
_decball_is_positive(const decball_t x, gr_ctx_t ctx)
{
    /* mid > rad */
    decfloat_t t;
    int c;
    if (DECMAG_IS_INF(&x->rad)) return 0;
    if (_decfloat_sgn(&x->mid, ctx) <= 0) return 0;
    if (DECMAG_IS_ZERO(&x->rad)) return 1;
    decfloat_init(t, ctx);
    GR_MUST_SUCCEED(_decmag_get_decfloat(t, &x->rad, ctx));
    c = _decfloat_cmp(&x->mid, t, ctx);
    decfloat_clear(t, ctx);
    return c > 0;
}

int
_decball_is_negative(const decball_t x, gr_ctx_t ctx)
{
    decball_t t;
    int r;
    *t = *x;
    t->mid.m.size = -t->mid.m.size;
    r = _decball_is_positive(t, ctx);
    return r;
}

int
_decball_is_nonnegative(const decball_t x, gr_ctx_t ctx)
{
    /* mid >= rad */
    decfloat_t t;
    int c;
    if (DECMAG_IS_INF(&x->rad)) return 0;
    if (_decfloat_sgn(&x->mid, ctx) < 0) return 0;
    if (DECMAG_IS_ZERO(&x->rad)) return 1;
    decfloat_init(t, ctx);
    GR_MUST_SUCCEED(_decmag_get_decfloat(t, &x->rad, ctx));
    c = _decfloat_cmp(&x->mid, t, ctx);
    decfloat_clear(t, ctx);
    return c >= 0;
}

int
_decball_is_nonpositive(const decball_t x, gr_ctx_t ctx)
{
    decball_t t;
    *t = *x;
    t->mid.m.size = -t->mid.m.size;
    return _decball_is_nonnegative(t, ctx);
}

int _decball_contains_negative(const decball_t x, gr_ctx_t ctx) { return !_decball_is_nonnegative(x, ctx); }
int _decball_contains_positive(const decball_t x, gr_ctx_t ctx) { return !_decball_is_nonpositive(x, ctx); }
int _decball_contains_nonnegative(const decball_t x, gr_ctx_t ctx) { return !_decball_is_negative(x, ctx); }
int _decball_contains_nonpositive(const decball_t x, gr_ctx_t ctx) { return !_decball_is_positive(x, ctx); }

truth_t
decball_is_integer(const decball_t x, gr_ctx_t ctx)
{
    decfloat_t lo, hi, flo;
    int status, contains_int;

    if (DECMAG_IS_ZERO(&x->rad))
        return _decfloat_is_int(&x->mid, ctx) ? T_TRUE : T_FALSE;
    if (DECMAG_IS_INF(&x->rad))
        return T_UNKNOWN;

    /* does [mid - rad, mid + rad] contain an integer? */
    decfloat_init(lo, ctx);
    decfloat_init(hi, ctx);
    decfloat_init(flo, ctx);

    {
        decfloat_t r;
        decfloat_init(r, ctx);
        GR_MUST_SUCCEED(_decmag_get_decfloat(r, &x->rad, ctx));
        status = _decfloat_add(lo, &x->mid, r, 1, DECIMAL_PREC_EXACT, DECIMAL_RND_DOWN | DECIMAL_RND_NOLIMITS, NULL, NULL, ctx);
        status |= _decfloat_add(hi, &x->mid, r, 0, DECIMAL_PREC_EXACT, DECIMAL_RND_DOWN | DECIMAL_RND_NOLIMITS, NULL, NULL, ctx);
        decfloat_clear(r, ctx);
    }

    if (status != GR_SUCCESS)
    {
        contains_int = 1;
    }
    else
    {
        /* ceil(lo) <= hi */
        GR_MUST_SUCCEED(_decfloat_round_to_int(flo, lo, DECIMAL_RND_CEIL, DECIMAL_PREC_EXACT, DECIMAL_RND_DOWN | DECIMAL_RND_NOLIMITS, NULL, NULL, ctx));
        contains_int = (_decfloat_cmp(flo, hi, ctx) <= 0);
    }

    decfloat_clear(lo, ctx);
    decfloat_clear(hi, ctx);
    decfloat_clear(flo, ctx);

    return contains_int ? T_UNKNOWN : T_FALSE;
}

int
decball_get_fmpz(fmpz_t res, const decball_t x, gr_ctx_t ctx)
{
    truth_t t = decball_is_integer(x, ctx);
    if (t == T_TRUE)
        return decfloat_get_fmpz(res, &x->mid, ctx);
    if (t == T_FALSE)
        return GR_DOMAIN;
    return GR_UNABLE;
}

int
decball_get_fmpq(fmpq_t res, const decball_t x, gr_ctx_t ctx)
{
    if (DECMAG_IS_ZERO(&x->rad))
        return decfloat_get_fmpq(res, &x->mid, ctx);
    return GR_UNABLE;
}

int
decball_get_d(double * res, const decball_t x, gr_ctx_t ctx)
{
    return decfloat_get_d(res, &x->mid, ctx);
}

int
decball_cmp(int * res, const decball_t x, const decball_t y, gr_ctx_t ctx)
{
    if ((DECMAG_IS_ZERO(&x->rad) && DECMAG_IS_ZERO(&y->rad)) || !_decball_overlaps(x, y, ctx))
    {
        *res = _decfloat_cmp(&x->mid, &y->mid, ctx);
        return GR_SUCCESS;
    }

    *res = 0;
    return GR_UNABLE;
}

int
decball_cmpabs(int * res, const decball_t x, const decball_t y, gr_ctx_t ctx)
{
    decball_t t, u;
    *t = *x;
    *u = *y;
    if (t->mid.m.size < 0) t->mid.m.size = -t->mid.m.size;
    if (u->mid.m.size < 0) u->mid.m.size = -u->mid.m.size;
    return decball_cmp(res, t, u, ctx);
}

int
decball_sgn(decball_t res, const decball_t x, gr_ctx_t ctx)
{
    if (_decball_is_positive(x, ctx))
        return decball_one(res, ctx);
    if (_decball_is_negative(x, ctx))
        return decball_neg_one(res, ctx);
    if (DECFLOAT_IS_ZERO(&x->mid) && DECMAG_IS_ZERO(&x->rad))
        return decball_zero(res, ctx);

    /* sign in [-1, 1] */
    decfloat_zero(&res->mid, ctx);
    _decmag_one(&res->rad, ctx);
    return GR_SUCCESS;
}

/* ------------------------------------------------------------------------- */
/*    Output, random                                                         */
/* ------------------------------------------------------------------------- */

char *
decball_get_str(const decball_t x, gr_ctx_t ctx)
{
    char stack_buf[128];
    char rad_digits[24];
    char * digits;
    char * s;
    slong L, Lr, n, need, pos, bound;
    fmpz_t t, tr;
    int sci = (DECIMAL_CTX_FLAGS(ctx) & DECIMAL_WRITE_SCIENTIFIC) != 0;

    if (DECMAG_IS_ZERO(&x->rad) || DECFLOAT_IS_ZERO(&x->mid) || DECMAG_IS_INF(&x->rad))
    {
        char * m = decfloat_get_str(&x->mid, ctx);
        char * r;

        if (DECMAG_IS_ZERO(&x->rad))
            return m;

        r = _decmag_get_str(&x->rad, ctx);
        s = flint_malloc(strlen(m) + strlen(r) + 8);
        strcpy(s, "[");
        strcat(s, m);
        strcat(s, " +/- ");
        strcat(s, r);
        strcat(s, "]");
        flint_free(m);
        flint_free(r);
        return s;
    }

    n = FLINT_ABS(x->mid.m.size);
    need = n * DECIMAL_CTX_E(ctx) + 1;
    digits = (need <= (slong) sizeof(stack_buf)) ? stack_buf : flint_malloc(need);

    fmpz_init(t);
    fmpz_init(tr);
    L = _decfloat_get_digits(digits, t, &x->mid, ctx);
    Lr = _decmag_get_digits(rad_digits, tr, &x->rad, ctx);

    bound = _decimal_write_number_bound(L, t) + _decimal_write_number_bound(Lr, tr) + 8;
    s = flint_malloc(bound);
    pos = 0;
    s[pos++] = '[';
    pos += _decimal_write_number(s + pos, x->mid.m.size < 0, digits, L, t, sci ? 1 : 0);
    memcpy(s + pos, " +/- ", 5);
    pos += 5;
    pos += _decimal_write_number(s + pos, 0, rad_digits, Lr, tr, sci ? 1 : 2);
    s[pos++] = ']';
    s[pos] = '\0';

    if (digits != stack_buf)
        flint_free(digits);
    fmpz_clear(t);
    fmpz_clear(tr);
    return s;
}

int
decball_write(gr_stream_t out, const decball_t x, gr_ctx_t ctx)
{
    return gr_stream_write_free(out, decball_get_str(x, ctx));
}

int
decball_randtest(decball_t res, flint_rand_t state, gr_ctx_t ctx)
{
    int status;

    /* like arb_randtest, only finite balls are generated (special values
       break the ring axioms tested by the generic test suite) */
    do
    {
        status = decfloat_randtest(&res->mid, state, ctx);
    }
    while (!DECFLOAT_IS_FINITE(&res->mid));

    if (n_randint(state, 2) || DECFLOAT_IS_ZERO(&res->mid))
    {
        _decmag_zero(&res->rad, ctx);
    }
    else
    {
        /* radius relative to the midpoint */
        slong E;
        if (decfloat_get_sci_exp_si(&E, &res->mid, ctx))
        {
            slong k = E - DECIMAL_CTX_PREC(ctx) + n_randint(state, 4) - 1;
            ulong m = 1 + n_randint(state, DECIMAL_CTX_RAD_POW(ctx) - 1);
            fmpz_t t;
            fmpz_init_set_si(t, k);
            if (n_randint(state, 10) == 0)
                fmpz_add_si(t, t, (slong) n_randint(state, 2 * DECIMAL_CTX_PREC(ctx) + 4) - DECIMAL_CTX_PREC(ctx));
            _decmag_set_ui_10exp_fmpz(&res->rad, m, t, ctx);
            fmpz_clear(t);
        }
        else
        {
            _decmag_zero(&res->rad, ctx);
        }
    }

    return status;
}

/* ------------------------------------------------------------------------- */
/*    Method table                                                           */
/* ------------------------------------------------------------------------- */

static truth_t
_decball_ctx_is_canonical(gr_ctx_t ctx)
{
    return T_FALSE;
}

int _decball_methods_initialized = 0;
gr_static_method_table _decball_methods;

gr_method_tab_input _decball_methods_input[] =
{
    {GR_METHOD_CTX_CLEAR,       (gr_funcptr) decimal_ctx_clear},
    {GR_METHOD_CTX_WRITE,       (gr_funcptr) decimal_ctx_write},
    {GR_METHOD_CTX_IS_RING,     (gr_funcptr) gr_generic_ctx_predicate_true},
    {GR_METHOD_CTX_IS_COMMUTATIVE_RING, (gr_funcptr) gr_generic_ctx_predicate_true},
    {GR_METHOD_CTX_IS_INTEGRAL_DOMAIN,  (gr_funcptr) gr_generic_ctx_predicate_true},
    {GR_METHOD_CTX_IS_FIELD,            (gr_funcptr) gr_generic_ctx_predicate_true},
    {GR_METHOD_CTX_IS_UNIQUE_FACTORIZATION_DOMAIN, (gr_funcptr) gr_generic_ctx_predicate_true},
    {GR_METHOD_CTX_IS_FINITE,   (gr_funcptr) gr_generic_ctx_predicate_false},
    {GR_METHOD_CTX_IS_FINITE_CHARACTERISTIC, (gr_funcptr) gr_generic_ctx_predicate_false},
    {GR_METHOD_CTX_IS_ALGEBRAICALLY_CLOSED, (gr_funcptr) gr_generic_ctx_predicate_false},
    {GR_METHOD_CTX_IS_ORDERED_RING, (gr_funcptr) gr_generic_ctx_predicate_true},
    {GR_METHOD_CTX_IS_RATIONAL_VECTOR_SPACE, (gr_funcptr) _decimal_ctx_is_vector_space},
    {GR_METHOD_CTX_IS_REAL_VECTOR_SPACE, (gr_funcptr) _decimal_ctx_is_vector_space},
    {GR_METHOD_CTX_IS_COMPLEX_VECTOR_SPACE, (gr_funcptr) gr_generic_ctx_predicate_false},
    {GR_METHOD_CTX_IS_EXACT,    (gr_funcptr) gr_generic_ctx_predicate_false},
    {GR_METHOD_CTX_IS_CANONICAL, (gr_funcptr) _decball_ctx_is_canonical},
    {GR_METHOD_CTX_HAS_REAL_PREC, (gr_funcptr) gr_generic_ctx_predicate_true},
    {GR_METHOD_CTX_SET_REAL_PREC, (gr_funcptr) _decimal_ctx_set_real_prec},
    {GR_METHOD_CTX_GET_REAL_PREC, (gr_funcptr) _decimal_ctx_get_real_prec},

    {GR_METHOD_INIT,            (gr_funcptr) decball_init},
    {GR_METHOD_CLEAR,           (gr_funcptr) decball_clear},
    {GR_METHOD_SWAP,            (gr_funcptr) decball_swap},
    {GR_METHOD_SET_SHALLOW,     (gr_funcptr) decball_set_shallow},
    {GR_METHOD_RANDTEST,        (gr_funcptr) decball_randtest},
    {GR_METHOD_WRITE,           (gr_funcptr) decball_write},
    {GR_METHOD_ZERO,            (gr_funcptr) decball_zero},
    {GR_METHOD_ONE,             (gr_funcptr) decball_one},
    {GR_METHOD_NEG_ONE,         (gr_funcptr) decball_neg_one},
    {GR_METHOD_IS_ZERO,         (gr_funcptr) decball_is_zero},
    {GR_METHOD_IS_ONE,          (gr_funcptr) decball_is_one},
    {GR_METHOD_IS_NEG_ONE,      (gr_funcptr) decball_is_neg_one},
    {GR_METHOD_EQUAL,           (gr_funcptr) decball_equal},
    {GR_METHOD_SET,             (gr_funcptr) decball_set},
    {GR_METHOD_SET_SI,          (gr_funcptr) decball_set_si},
    {GR_METHOD_SET_UI,          (gr_funcptr) decball_set_ui},
    {GR_METHOD_SET_FMPZ,        (gr_funcptr) decball_set_fmpz},
    {GR_METHOD_SET_FMPQ,        (gr_funcptr) decball_set_fmpq},
    {GR_METHOD_SET_D,           (gr_funcptr) decball_set_d},
    {GR_METHOD_SET_STR,         (gr_funcptr) decball_set_str},
    {GR_METHOD_SET_FMPZ_10EXP_FMPZ, (gr_funcptr) decball_set_fmpz_10exp_fmpz},
    {GR_METHOD_SET_INTERVAL_MID_RAD, (gr_funcptr) decball_set_interval_mid_rad},
    {GR_METHOD_SET_INTERVAL_MID_INF, (gr_funcptr) decball_set_interval_mid_inf},
    {GR_METHOD_SET_OTHER,       (gr_funcptr) decball_set_other},
    {GR_METHOD_GET_FMPZ,        (gr_funcptr) decball_get_fmpz},
    {GR_METHOD_GET_FMPQ,        (gr_funcptr) decball_get_fmpq},
    {GR_METHOD_GET_D,           (gr_funcptr) decball_get_d},

    {GR_METHOD_NEG,             (gr_funcptr) decball_neg},
    {GR_METHOD_ADD,             (gr_funcptr) decball_add},
    {GR_METHOD_ADD_UI,          (gr_funcptr) decball_add_ui},
    {GR_METHOD_ADD_SI,          (gr_funcptr) decball_add_si},
    {GR_METHOD_ADD_FMPZ,        (gr_funcptr) decball_add_fmpz},
    {GR_METHOD_SUB,             (gr_funcptr) decball_sub},
    {GR_METHOD_SUB_UI,          (gr_funcptr) decball_sub_ui},
    {GR_METHOD_SUB_SI,          (gr_funcptr) decball_sub_si},
    {GR_METHOD_SUB_FMPZ,        (gr_funcptr) decball_sub_fmpz},
    {GR_METHOD_MUL,             (gr_funcptr) decball_mul},
    {GR_METHOD_MUL_UI,          (gr_funcptr) decball_mul_ui},
    {GR_METHOD_MUL_SI,          (gr_funcptr) decball_mul_si},
    {GR_METHOD_MUL_FMPZ,        (gr_funcptr) decball_mul_fmpz},
    {GR_METHOD_MUL_TWO,         (gr_funcptr) decball_mul_two},
    {GR_METHOD_SQR,             (gr_funcptr) decball_sqr},
    {GR_METHOD_DIV,             (gr_funcptr) decball_div},
    {GR_METHOD_DIV_UI,          (gr_funcptr) decball_div_ui},
    {GR_METHOD_DIV_SI,          (gr_funcptr) decball_div_si},
    {GR_METHOD_DIV_FMPZ,        (gr_funcptr) decball_div_fmpz},
    {GR_METHOD_INV,             (gr_funcptr) decball_inv},
    {GR_METHOD_MUL_2EXP_SI,     (gr_funcptr) decball_mul_2exp_si},
    {GR_METHOD_MUL_2EXP_FMPZ,   (gr_funcptr) decball_mul_2exp_fmpz},
    {GR_METHOD_POW_UI,          (gr_funcptr) decball_pow_ui},
    {GR_METHOD_POW_SI,          (gr_funcptr) decball_pow_si},
    {GR_METHOD_POW_FMPZ,        (gr_funcptr) decball_pow_fmpz},
    {GR_METHOD_SQRT,            (gr_funcptr) decball_sqrt},
    {GR_METHOD_RSQRT,           (gr_funcptr) decball_rsqrt},
    {GR_METHOD_UINF,            (gr_funcptr) gr_not_in_domain},
    {GR_METHOD_FLOOR,           (gr_funcptr) decball_floor},
    {GR_METHOD_CEIL,            (gr_funcptr) decball_ceil},
    {GR_METHOD_TRUNC,           (gr_funcptr) decball_trunc},
    {GR_METHOD_NINT,            (gr_funcptr) decball_nint},
    {GR_METHOD_ABS,             (gr_funcptr) decball_abs},
    {GR_METHOD_CONJ,            (gr_funcptr) decball_set},
    {GR_METHOD_RE,              (gr_funcptr) decball_set},
    {GR_METHOD_IM,              (gr_funcptr) decball_zero},
    {GR_METHOD_SGN,             (gr_funcptr) decball_sgn},
    {GR_METHOD_CSGN,            (gr_funcptr) decball_sgn},
    {GR_METHOD_CMP,             (gr_funcptr) decball_cmp},
    {GR_METHOD_CMPABS,          (gr_funcptr) decball_cmpabs},
    {GR_METHOD_I,               (gr_funcptr) gr_not_in_domain},
    /* balls represent real numbers: no infinities or undefined values */
    {GR_METHOD_POS_INF,         (gr_funcptr) gr_not_in_domain},
    {GR_METHOD_NEG_INF,         (gr_funcptr) gr_not_in_domain},
    {GR_METHOD_UINF,            (gr_funcptr) gr_not_in_domain},
    {GR_METHOD_UNDEFINED,       (gr_funcptr) gr_not_in_domain},
    {GR_METHOD_UNKNOWN,         (gr_funcptr) gr_not_in_domain},
    {GR_METHOD_ATAN2, (gr_funcptr) decball_atan2},
    {GR_METHOD_POW,             (gr_funcptr) decball_pow},
    {GR_METHOD_MAT_DET,         (gr_funcptr) gr_mat_det_generic_field},
    {GR_METHOD_MAT_FIND_NONZERO_PIVOT, (gr_funcptr) gr_mat_find_nonzero_pivot_large_abs},
    {0,                         (gr_funcptr) NULL},
};

POP_OPTIONS
