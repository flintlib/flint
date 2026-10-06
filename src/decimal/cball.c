/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/*
    Complex decimal balls: pairs of real decimal balls (rectangular
    enclosures), like acb.
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

#define PREC(ctx) DECIMAL_CTX_PREC(ctx)
#define EXACT_RND (DECIMAL_RND_DOWN | DECIMAL_RND_NOLIMITS)

/* status for real-valued conversions of a nonreal ball */
static int
_nonreal_status(const deccball_t x, gr_ctx_t ctx)
{
    return _decball_contains_zero(&x->im, ctx) ? GR_UNABLE : GR_DOMAIN;
}

/* ------------------------------------------------------------------------- */
/*    Arithmetic                                                             */
/* ------------------------------------------------------------------------- */

/* res = x1 y1 +/- x2 y2 with the midpoint products computed exactly and
   rounded once */
int
_decball_dot2(decball_t res, const decball_t x1, const decball_t y1, const decball_t x2, const decball_t y2, int subtract, slong prec, gr_ctx_t ctx)
{
    decimal_rounding_info info;
    decmag_t err, rad, t, u;
    int status;

    if (!_decball_is_finite(x1, ctx) || !_decball_is_finite(y1, ctx) || !_decball_is_finite(x2, ctx) || !_decball_is_finite(y2, ctx))
    {
        decball_t a, b;
        decball_init(a, ctx);
        decball_init(b, ctx);
        status = decball_mul_round(a, x1, y1, prec, ctx);
        status |= decball_mul_round(b, x2, y2, prec, ctx);
        if (status == GR_SUCCESS)
        {
            if (subtract)
                status = decball_sub_round(res, a, b, prec, ctx);
            else
                status = decball_add_round(res, a, b, prec, ctx);
        }
        decball_clear(a, ctx);
        decball_clear(b, ctx);
        return status;
    }

    _decmag_init(err, ctx);
    _decmag_init(rad, ctx);
    _decmag_init(t, ctx);
    _decmag_init(u, ctx);

    /* radius bound first (midpoints may be overwritten) */
    if (!DECMAG_IS_ZERO(&x1->rad) || !DECMAG_IS_ZERO(&y1->rad))
    {
        _decmag_set_decfloat(t, &x1->mid, ctx);
        _decmag_set_decfloat(u, &y1->mid, ctx);
        _decmag_mul(rad, t, &y1->rad, ctx);
        _decmag_addmul(rad, u, &x1->rad, ctx);
        _decmag_addmul(rad, &x1->rad, &y1->rad, ctx);
    }

    if (!DECMAG_IS_ZERO(&x2->rad) || !DECMAG_IS_ZERO(&y2->rad))
    {
        _decmag_set_decfloat(t, &x2->mid, ctx);
        _decmag_set_decfloat(u, &y2->mid, ctx);
        _decmag_addmul(rad, t, &y2->rad, ctx);
        _decmag_addmul(rad, u, &x2->rad, ctx);
        _decmag_addmul(rad, &x2->rad, &y2->rad, ctx);
    }

    status = _decfloat_dot2(&res->mid, &x1->mid, &y1->mid, &x2->mid, &y2->mid, subtract, prec, DECIMAL_CTX_RND(ctx),
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
deccball_neg(deccball_t res, const deccball_t x, gr_ctx_t ctx)
{
    int status;
    status = decball_neg(&res->re, &x->re, ctx);
    status |= decball_neg(&res->im, &x->im, ctx);
    return status;
}

int
deccball_conj(deccball_t res, const deccball_t x, gr_ctx_t ctx)
{
    int status;
    status = decball_set(&res->re, &x->re, ctx);
    status |= decball_neg(&res->im, &x->im, ctx);
    return status;
}

int
deccball_re(deccball_t res, const deccball_t x, gr_ctx_t ctx)
{
    int status = decball_set(&res->re, &x->re, ctx);
    decball_zero(&res->im, ctx);
    return status;
}

int
deccball_im(deccball_t res, const deccball_t x, gr_ctx_t ctx)
{
    int status = decball_set(&res->re, &x->im, ctx);
    decball_zero(&res->im, ctx);
    return status;
}

int
deccball_add(deccball_t res, const deccball_t x, const deccball_t y, gr_ctx_t ctx)
{
    int status;
    status = decball_add_round(&res->re, &x->re, &y->re, PREC(ctx), ctx);
    status |= decball_add_round(&res->im, &x->im, &y->im, PREC(ctx), ctx);
    return status;
}

int
deccball_sub(deccball_t res, const deccball_t x, const deccball_t y, gr_ctx_t ctx)
{
    int status;
    status = decball_sub_round(&res->re, &x->re, &y->re, PREC(ctx), ctx);
    status |= decball_sub_round(&res->im, &x->im, &y->im, PREC(ctx), ctx);
    return status;
}

#define IS_ZERO_BALL(b) (DECFLOAT_IS_ZERO(&(b)->mid) && DECMAG_IS_ZERO(&(b)->rad))

int
deccball_mul(deccball_t res, const deccball_t x, const deccball_t y, gr_ctx_t ctx)
{
    decball_t re, im;
    int status;

    decball_init(re, ctx);
    decball_init(im, ctx);

    if (IS_ZERO_BALL(&x->im))
    {
        status = decball_mul_round(re, &x->re, &y->re, PREC(ctx), ctx);
        if (IS_ZERO_BALL(&y->im))
            decball_zero(im, ctx);
        else
            status |= decball_mul_round(im, &x->re, &y->im, PREC(ctx), ctx);
    }
    else if (IS_ZERO_BALL(&y->im))
    {
        status = decball_mul_round(re, &x->re, &y->re, PREC(ctx), ctx);
        status |= decball_mul_round(im, &x->im, &y->re, PREC(ctx), ctx);
    }
    else
    {
        status = _decball_dot2(re, &x->re, &y->re, &x->im, &y->im, 1, PREC(ctx), ctx);
        status |= _decball_dot2(im, &x->re, &y->im, &x->im, &y->re, 0, PREC(ctx), ctx);
    }

    decball_swap(&res->re, re, ctx);
    decball_swap(&res->im, im, ctx);
    decball_clear(re, ctx);
    decball_clear(im, ctx);
    return status;
}

int
deccball_sqr(deccball_t res, const deccball_t x, gr_ctx_t ctx)
{
    return deccball_mul(res, x, x, ctx);
}

DECIMAL_DRIVER int
_deccball_acb_unary(deccball_t res, const deccball_t x, int method, gr_ctx_t ctx)
{
    deccball_srcptr a[1];
    deccball_ptr r[1];
    a[0] = x;
    r[0] = res;
    return _deccball_acb_gr(r, 1, a, 1, 0, 0, method, ctx);
}

int
deccball_div(deccball_t res, const deccball_t x, const deccball_t y, gr_ctx_t ctx)
{
    decball_t re, im, den;
    int status;

    decball_init(re, ctx);
    decball_init(im, ctx);

    if (IS_ZERO_BALL(&y->im))
    {
        status = decball_div_round(re, &x->re, &y->re, PREC(ctx), ctx);
        status |= decball_div_round(im, &x->im, &y->re, PREC(ctx), ctx);
    }
    else if (IS_ZERO_BALL(&y->re))
    {
        /* (a + bi) / (di) = b/d - (a/d) i */
        status = decball_div_round(re, &x->im, &y->im, PREC(ctx), ctx);
        status |= decball_div_round(im, &x->re, &y->im, PREC(ctx), ctx);
        if (status == GR_SUCCESS)
            status = decball_neg(im, im, ctx);
    }
    else
    {
        decball_init(den, ctx);
        status = _decball_dot2(den, &y->re, &y->re, &y->im, &y->im, 0, PREC(ctx) + 5, ctx);
        status |= _decball_dot2(re, &x->re, &y->re, &x->im, &y->im, 0, PREC(ctx) + 5, ctx);
        status |= _decball_dot2(im, &x->im, &y->re, &x->re, &y->im, 1, PREC(ctx) + 5, ctx);
        if (status == GR_SUCCESS)
        {
            status = decball_div_round(re, re, den, PREC(ctx), ctx);
            status |= decball_div_round(im, im, den, PREC(ctx), ctx);
        }
        decball_clear(den, ctx);

        /* the rectangular bound for |y|^2 may include zero even when y
           excludes it: fall back to acb */
        if (status == GR_UNABLE && !_deccball_contains_zero(y, ctx))
        {
            deccball_srcptr a[2];
            deccball_ptr r[1];
            decball_clear(re, ctx);
            decball_clear(im, ctx);
            a[0] = x;
            a[1] = y;
            r[0] = res;
            return _deccball_acb_gr(r, 1, a, 2, 0, 0, GR_METHOD_DIV, ctx);
        }
    }

    decball_swap(&res->re, re, ctx);
    decball_swap(&res->im, im, ctx);
    decball_clear(re, ctx);
    decball_clear(im, ctx);
    return status;
}

int
deccball_inv(deccball_t res, const deccball_t x, gr_ctx_t ctx)
{
    deccball_t one;
    int status;
    deccball_init(one, ctx);
    GR_MUST_SUCCEED(deccball_one(one, ctx));
    status = deccball_div(res, one, x, ctx);
    deccball_clear(one, ctx);
    return status;
}

int
deccball_sqrt(deccball_t res, const deccball_t x, gr_ctx_t ctx)
{
    if (_deccball_is_real(x, ctx))
    {
        int status = decball_sqrt(&res->re, &x->re, ctx);
        if (status == GR_SUCCESS)
        {
            decball_zero(&res->im, ctx);
            return status;
        }
        if (status == GR_DOMAIN)
        {
            /* sqrt(-t) = i sqrt(t) */
            decball_t t;
            decball_init(t, ctx);
            GR_MUST_SUCCEED(decball_neg(t, &x->re, ctx));
            status = decball_sqrt(&res->im, t, ctx);
            decball_zero(&res->re, ctx);
            decball_clear(t, ctx);
            return status;
        }
    }

    return _deccball_acb_unary(res, x, GR_METHOD_SQRT, ctx);
}

int
deccball_rsqrt(deccball_t res, const deccball_t x, gr_ctx_t ctx)
{
    if (_deccball_is_real(x, ctx))
    {
        int status = decball_rsqrt(&res->re, &x->re, ctx);
        if (status == GR_SUCCESS)
        {
            decball_zero(&res->im, ctx);
            return status;
        }
    }

    return _deccball_acb_unary(res, x, GR_METHOD_RSQRT, ctx);
}

int
deccball_abs(deccball_t res, const deccball_t x, gr_ctx_t ctx)
{
    int status;

    if (_deccball_is_real(x, ctx))
        status = decball_abs(&res->re, &x->re, ctx);
    else if (IS_ZERO_BALL(&x->re))
        status = decball_abs(&res->re, &x->im, ctx);
    else
        return _deccball_acb_unary(res, x, GR_METHOD_ABS, ctx);

    decball_zero(&res->im, ctx);
    return status;
}

int
deccball_arg(deccball_t res, const deccball_t x, gr_ctx_t ctx)
{
    return _deccball_acb_unary(res, x, GR_METHOD_ARG, ctx);
}

int
deccball_sgn(deccball_t res, const deccball_t x, gr_ctx_t ctx)
{
    if (_deccball_is_real(x, ctx))
    {
        int status = decball_sgn(&res->re, &x->re, ctx);
        decball_zero(&res->im, ctx);
        return status;
    }

    return _deccball_acb_unary(res, x, GR_METHOD_SGN, ctx);
}

int
deccball_csgn(deccball_t res, const deccball_t x, gr_ctx_t ctx)
{
    return _deccball_acb_unary(res, x, GR_METHOD_CSGN, ctx);
}

typedef int (*deccball_binary_fn)(deccball_t, const deccball_t, const deccball_t, gr_ctx_t);

/* op(x, y) for an integer y of type DECIMAL_SCALAR_* */
DECIMAL_DRIVER int
_deccball_scalar_op(deccball_t res, const deccball_t x, const void * y, int type, deccball_binary_fn op, gr_ctx_t ctx)
{
    deccball_t t;
    int status;
    deccball_init(t, ctx);
    status = _decfloat_set_scalar_exact(&t->re.mid, y, type, ctx);
    if (status == GR_SUCCESS)
        status = op(res, x, t, ctx);
    deccball_clear(t, ctx);
    return status;
}

#define DEF_SCALAR(name, T, type, op) \
int name(deccball_t res, const deccball_t x, T y, gr_ctx_t ctx) \
{ \
    return _deccball_scalar_op(res, x, &y, type, op, ctx); \
}

#define DEF_SCALAR_FMPZ(name, op) \
int name(deccball_t res, const deccball_t x, const fmpz_t y, gr_ctx_t ctx) \
{ \
    return _deccball_scalar_op(res, x, y, 2, op, ctx); \
}

DEF_SCALAR(deccball_add_ui, ulong, 0, deccball_add)
DEF_SCALAR(deccball_add_si, slong, 1, deccball_add)
DEF_SCALAR_FMPZ(deccball_add_fmpz, deccball_add)
DEF_SCALAR(deccball_sub_ui, ulong, 0, deccball_sub)
DEF_SCALAR(deccball_sub_si, slong, 1, deccball_sub)
DEF_SCALAR_FMPZ(deccball_sub_fmpz, deccball_sub)
DEF_SCALAR(deccball_mul_ui, ulong, 0, deccball_mul)
DEF_SCALAR(deccball_mul_si, slong, 1, deccball_mul)
DEF_SCALAR_FMPZ(deccball_mul_fmpz, deccball_mul)
DEF_SCALAR(deccball_div_ui, ulong, 0, deccball_div)
DEF_SCALAR(deccball_div_si, slong, 1, deccball_div)
DEF_SCALAR_FMPZ(deccball_div_fmpz, deccball_div)

int
deccball_mul_decball(deccball_t res, const deccball_t x, const decball_t y, gr_ctx_t ctx)
{
    int status;
    status = decball_mul_round(&res->re, &x->re, y, PREC(ctx), ctx);
    status |= decball_mul_round(&res->im, &x->im, y, PREC(ctx), ctx);
    return status;
}

int
deccball_div_decball(deccball_t res, const deccball_t x, const decball_t y, gr_ctx_t ctx)
{
    int status;
    status = decball_div_round(&res->re, &x->re, y, PREC(ctx), ctx);
    status |= decball_div_round(&res->im, &x->im, y, PREC(ctx), ctx);
    return status;
}

int
deccball_mul_i(deccball_t res, const deccball_t x, gr_ctx_t ctx)
{
    decball_t t;
    int status;
    decball_init(t, ctx);
    status = decball_neg(t, &x->im, ctx);
    status |= decball_set(&res->im, &x->re, ctx);
    decball_swap(&res->re, t, ctx);
    decball_clear(t, ctx);
    return status;
}

int
deccball_div_i(deccball_t res, const deccball_t x, gr_ctx_t ctx)
{
    decball_t t;
    int status;
    decball_init(t, ctx);
    status = decball_set(t, &x->im, ctx);
    status |= decball_neg(&res->im, &x->re, ctx);
    decball_swap(&res->re, t, ctx);
    decball_clear(t, ctx);
    return status;
}

int
deccball_mul_two(deccball_t res, const deccball_t x, gr_ctx_t ctx)
{
    return deccball_add(res, x, x, ctx);
}

int
deccball_mul_10exp_fmpz(deccball_t res, const deccball_t x, const fmpz_t e, gr_ctx_t ctx)
{
    int status;
    status = decball_mul_10exp_fmpz(&res->re, &x->re, e, ctx);
    status |= decball_mul_10exp_fmpz(&res->im, &x->im, e, ctx);
    return status;
}

int
deccball_mul_10exp_si(deccball_t res, const deccball_t x, slong e, gr_ctx_t ctx)
{
    int status;
    status = decball_mul_10exp_si(&res->re, &x->re, e, ctx);
    status |= decball_mul_10exp_si(&res->im, &x->im, e, ctx);
    return status;
}

int
deccball_mul_2exp_fmpz(deccball_t res, const deccball_t x, const fmpz_t e, gr_ctx_t ctx)
{
    int status;
    status = decball_mul_2exp_fmpz(&res->re, &x->re, e, ctx);
    status |= decball_mul_2exp_fmpz(&res->im, &x->im, e, ctx);
    return status;
}

int
deccball_mul_2exp_si(deccball_t res, const deccball_t x, slong e, gr_ctx_t ctx)
{
    int status;
    status = decball_mul_2exp_si(&res->re, &x->re, e, ctx);
    status |= decball_mul_2exp_si(&res->im, &x->im, e, ctx);
    return status;
}

#define DEF_ROUND_TO_INT(name, fn) \
int name(deccball_t res, const deccball_t x, gr_ctx_t ctx) \
{ \
    if (!_deccball_is_real(x, ctx)) \
        return _nonreal_status(x, ctx); \
    decball_zero(&res->im, ctx); \
    return fn(&res->re, &x->re, ctx); \
}

DEF_ROUND_TO_INT(deccball_floor, decball_floor)
DEF_ROUND_TO_INT(deccball_ceil, decball_ceil)
DEF_ROUND_TO_INT(deccball_trunc, decball_trunc)
DEF_ROUND_TO_INT(deccball_nint, decball_nint)

int
deccball_cmp(int * res, const deccball_t x, const deccball_t y, gr_ctx_t ctx)
{
    if (!_deccball_is_real(x, ctx) || !_deccball_is_real(y, ctx))
    {
        *res = 0;
        return (_decball_contains_zero(&x->im, ctx) && _decball_contains_zero(&y->im, ctx)) ? GR_UNABLE : GR_DOMAIN;
    }

    return decball_cmp(res, &x->re, &y->re, ctx);
}

int
deccball_cmpabs(int * res, const deccball_t x, const deccball_t y, gr_ctx_t ctx)
{
    decball_t s, t;
    int status;

    if (_deccball_is_real(x, ctx) && _deccball_is_real(y, ctx))
        return decball_cmpabs(res, &x->re, &y->re, ctx);

    decball_init(s, ctx);
    decball_init(t, ctx);
    status = _decball_dot2(s, &x->re, &x->re, &x->im, &x->im, 0, PREC(ctx) + 5, ctx);
    status |= _decball_dot2(t, &y->re, &y->re, &y->im, &y->im, 0, PREC(ctx) + 5, ctx);
    if (status == GR_SUCCESS)
        status = decball_cmp(res, s, t, ctx);
    else
        *res = 0;
    decball_clear(s, ctx);
    decball_clear(t, ctx);
    return status;
}

/* The remaining functions are not performance-critical. */
PUSH_OPTIONS
OPTIMIZE_OSIZE

/* ------------------------------------------------------------------------- */
/*    Basic functions                                                        */
/* ------------------------------------------------------------------------- */

int
deccball_zero(deccball_t res, gr_ctx_t ctx)
{
    decball_zero(&res->re, ctx);
    decball_zero(&res->im, ctx);
    return GR_SUCCESS;
}

int
deccball_one(deccball_t res, gr_ctx_t ctx)
{
    decball_zero(&res->im, ctx);
    return decball_one(&res->re, ctx);
}

int
deccball_neg_one(deccball_t res, gr_ctx_t ctx)
{
    decball_zero(&res->im, ctx);
    return decball_neg_one(&res->re, ctx);
}

int
deccball_i(deccball_t res, gr_ctx_t ctx)
{
    decball_zero(&res->re, ctx);
    return decball_one(&res->im, ctx);
}

truth_t
deccball_is_zero(const deccball_t x, gr_ctx_t ctx)
{
    return truth_and(decball_is_zero(&x->re, ctx), decball_is_zero(&x->im, ctx));
}

truth_t
deccball_is_one(const deccball_t x, gr_ctx_t ctx)
{
    return truth_and(decball_is_one(&x->re, ctx), decball_is_zero(&x->im, ctx));
}

truth_t
deccball_is_neg_one(const deccball_t x, gr_ctx_t ctx)
{
    return truth_and(decball_is_neg_one(&x->re, ctx), decball_is_zero(&x->im, ctx));
}

truth_t
deccball_equal(const deccball_t x, const deccball_t y, gr_ctx_t ctx)
{
    return truth_and(decball_equal(&x->re, &y->re, ctx), decball_equal(&x->im, &y->im, ctx));
}

truth_t
deccball_is_integer(const deccball_t x, gr_ctx_t ctx)
{
    return truth_and(decball_is_integer(&x->re, ctx), decball_is_zero(&x->im, ctx));
}

int
_deccball_contains_zero(const deccball_t x, gr_ctx_t ctx)
{
    return _decball_contains_zero(&x->re, ctx) && _decball_contains_zero(&x->im, ctx);
}

int
_deccball_contains(const deccball_t x, const deccball_t y, gr_ctx_t ctx)
{
    return _decball_contains(&x->re, &y->re, ctx) && _decball_contains(&x->im, &y->im, ctx);
}

int
_deccball_overlaps(const deccball_t x, const deccball_t y, gr_ctx_t ctx)
{
    return _decball_overlaps(&x->re, &y->re, ctx) && _decball_overlaps(&x->im, &y->im, ctx);
}

int
_deccball_contains_deccfloat(const deccball_t x, const deccfloat_t y, gr_ctx_t ctx)
{
    return _decball_contains_decfloat(&x->re, &y->re, ctx) && _decball_contains_decfloat(&x->im, &y->im, ctx);
}

/* ------------------------------------------------------------------------- */
/*    Assignment and conversion                                              */
/* ------------------------------------------------------------------------- */

int
deccball_set_round(deccball_t res, const deccball_t x, slong prec, gr_ctx_t ctx)
{
    int status;
    status = decball_set_round(&res->re, &x->re, prec, ctx);
    status |= decball_set_round(&res->im, &x->im, prec, ctx);
    return status;
}

int
deccball_set(deccball_t res, const deccball_t x, gr_ctx_t ctx)
{
    return deccball_set_round(res, x, PREC(ctx), ctx);
}

int
deccball_set_round2(deccball_t res, const deccball_t x, slong prec, slong rad_prec, gr_ctx_t ctx)
{
    int status;
    status = decball_set_round2(&res->re, &x->re, prec, rad_prec, ctx);
    status |= decball_set_round2(&res->im, &x->im, prec, rad_prec, ctx);
    return status;
}

int
deccball_mid(deccball_t res, const deccball_t x, gr_ctx_t ctx)
{
    int status;
    status = decball_mid(&res->re, &x->re, ctx);
    status |= decball_mid(&res->im, &x->im, ctx);
    return status;
}

int
deccball_shell(deccball_t res, const deccball_t x, gr_ctx_t ctx)
{
    int status;
    status = decball_shell(&res->re, &x->re, ctx);
    status |= decball_shell(&res->im, &x->im, ctx);
    return status;
}

int
deccball_add_rad(deccball_t res, const deccball_t x, const deccball_t r, gr_ctx_t ctx)
{
    return deccball_set_interval_mid_rad(res, x, r, ctx);
}

/* an infinite real radius (the radius +inf is real), as in acb */
int
deccball_set_interval_mid_inf(deccball_t res, const deccball_t m, gr_ctx_t ctx)
{
    int status = deccball_set(res, m, ctx);
    _decmag_inf(&res->re.rad, ctx);
    return status;
}

int
deccball_set_decball(deccball_t res, const decball_t x, gr_ctx_t ctx)
{
    decball_zero(&res->im, ctx);
    return decball_set(&res->re, x, ctx);
}

int
deccball_set_decball_decball(deccball_t res, const decball_t re, const decball_t im, gr_ctx_t ctx)
{
    int status;
    status = decball_set(&res->re, re, ctx);
    status |= decball_set(&res->im, im, ctx);
    return status;
}

int
deccball_set_decfloat(deccball_t res, const decfloat_t x, gr_ctx_t ctx)
{
    decball_zero(&res->im, ctx);
    return decball_set_decfloat(&res->re, x, ctx);
}

int
deccball_set_deccfloat(deccball_t res, const deccfloat_t x, gr_ctx_t ctx)
{
    int status;
    status = decball_set_decfloat(&res->re, &x->re, ctx);
    status |= decball_set_decfloat(&res->im, &x->im, ctx);
    return status;
}

#define DEF_SET_REAL(name, T, fn) \
int name(deccball_t res, T x, gr_ctx_t ctx) \
{ \
    decball_zero(&res->im, ctx); \
    return fn(&res->re, x, ctx); \
}

DEF_SET_REAL(deccball_set_si, slong, decball_set_si)
DEF_SET_REAL(deccball_set_ui, ulong, decball_set_ui)
DEF_SET_REAL(deccball_set_fmpz, const fmpz_t, decball_set_fmpz)
DEF_SET_REAL(deccball_set_fmpq, const fmpq_t, decball_set_fmpq)
DEF_SET_REAL(deccball_set_d, double, decball_set_d)

int
deccball_set_fmpz_10exp_fmpz(deccball_t res, const fmpz_t m, const fmpz_t e, gr_ctx_t ctx)
{
    decball_zero(&res->im, ctx);
    return decball_set_fmpz_10exp_fmpz(&res->re, m, e, ctx);
}

int
deccball_set_fmpz_2exp_fmpz(deccball_t res, const fmpz_t m, const fmpz_t e, gr_ctx_t ctx)
{
    decfloat_t t;
    int status;
    decfloat_init(t, ctx);
    status = decfloat_set_round_fmpz_2exp_fmpz(t, m, e, DECIMAL_PREC_EXACT, EXACT_RND, ctx);
    if (status == GR_SUCCESS)
        status = deccball_set_decfloat(res, t, ctx);
    decfloat_clear(t, ctx);
    return status;
}

int
deccball_set_str(deccball_t res, const char * s, gr_ctx_t ctx)
{
    decfloat_t t;
    int status;

    decfloat_init(t, ctx);
    status = _decfloat_set_str_literal(t, s, DECIMAL_PREC_EXACT, EXACT_RND, ctx);

    if (status == GR_SUCCESS)
        status = deccball_set_decfloat(res, t, ctx);
    else if (status == GR_DOMAIN)
        status = gr_generic_set_str_ring_exponents(res, s, ctx);

    decfloat_clear(t, ctx);
    return status;
}

int
deccball_set_interval_mid_rad(deccball_t res, const deccball_t m, const deccball_t r, gr_ctx_t ctx)
{
    decmag_t t, u;
    int status;

    _decmag_init(t, ctx);
    _decmag_init(u, ctx);
    decball_get_abs_ubound(t, &r->re, ctx);
    decball_get_abs_ubound(u, &r->im, ctx);
    status = deccball_set(res, m, ctx);
    _decmag_add(&res->re.rad, &res->re.rad, t, ctx);
    _decmag_add(&res->im.rad, &res->im.rad, u, ctx);
    _decmag_clear(t, ctx);
    _decmag_clear(u, ctx);
    return status;
}

int
deccball_set_acb(deccball_t res, const acb_t x, gr_ctx_t ctx)
{
    int status;
    status = decball_set_arb(&res->re, acb_realref(x), ctx);
    status |= decball_set_arb(&res->im, acb_imagref(x), ctx);
    return status;
}

int
deccball_get_acb(acb_t res, const deccball_t x, slong prec_bits, gr_ctx_t ctx)
{
    int status;
    status = decball_get_arb(acb_realref(res), &x->re, prec_bits, ctx);
    status |= decball_get_arb(acb_imagref(res), &x->im, prec_bits, ctx);
    return status;
}

int
deccball_get_fmpz(fmpz_t res, const deccball_t x, gr_ctx_t ctx)
{
    if (!_deccball_is_real(x, ctx))
        return _nonreal_status(x, ctx);
    return decball_get_fmpz(res, &x->re, ctx);
}

int
deccball_get_fmpq(fmpq_t res, const deccball_t x, gr_ctx_t ctx)
{
    if (!_deccball_is_real(x, ctx))
        return _nonreal_status(x, ctx);
    return decball_get_fmpq(res, &x->re, ctx);
}

int
deccball_get_si(slong * res, const deccball_t x, gr_ctx_t ctx)
{
    fmpz_t t;
    int status;
    fmpz_init(t);
    status = deccball_get_fmpz(t, x, ctx);
    if (status == GR_SUCCESS)
    {
        if (fmpz_fits_si(t))
            *res = fmpz_get_si(t);
        else
            status = GR_DOMAIN;
    }
    fmpz_clear(t);
    return status;
}

int
deccball_get_ui(ulong * res, const deccball_t x, gr_ctx_t ctx)
{
    fmpz_t t;
    int status;
    fmpz_init(t);
    status = deccball_get_fmpz(t, x, ctx);
    if (status == GR_SUCCESS)
    {
        if (fmpz_sgn(t) >= 0 && fmpz_abs_fits_ui(t))
            *res = fmpz_get_ui(t);
        else
            status = GR_DOMAIN;
    }
    fmpz_clear(t);
    return status;
}

int
deccball_get_d(double * res, const deccball_t x, gr_ctx_t ctx)
{
    if (!_deccball_is_real(x, ctx))
        return _nonreal_status(x, ctx);
    return decball_get_d(res, &x->re, ctx);
}

int
deccball_get_mid(deccfloat_t res, const deccball_t x, gr_ctx_t ctx)
{
    int status;
    status = decball_get_mid(&res->re, &x->re, ctx);
    status |= decball_get_mid(&res->im, &x->im, ctx);
    return status;
}

int
deccball_set_other(deccball_t res, gr_srcptr x, gr_ctx_t x_ctx, gr_ctx_t ctx)
{
    switch (x_ctx->which_ring)
    {
        case GR_CTX_FMPZ:
            return deccball_set_fmpz(res, x, ctx);

        case GR_CTX_FMPQ:
            return deccball_set_fmpq(res, x, ctx);

        case GR_CTX_REAL_FLOAT_ARF:
        case GR_CTX_RR_ARB:
        case GR_CTX_DECFLOAT:
        case GR_CTX_DECBALL:
            decball_zero(&res->im, ctx);
            return decball_set_other(&res->re, x, x_ctx, ctx);

        case GR_CTX_COMPLEX_FLOAT_ACF:
            {
                deccfloat_t t;
                int status;
                deccfloat_init(t, ctx);
                status = decfloat_set_round_arf(&t->re, acf_realref((acf_srcptr) x), DECIMAL_PREC_EXACT, EXACT_RND, ctx);
                status |= decfloat_set_round_arf(&t->im, acf_imagref((acf_srcptr) x), DECIMAL_PREC_EXACT, EXACT_RND, ctx);
                if (status == GR_SUCCESS)
                    status = deccball_set_deccfloat(res, t, ctx);
                deccfloat_clear(t, ctx);
                return status;
            }

        case GR_CTX_CC_ACB:
            return deccball_set_acb(res, x, ctx);

        case GR_CTX_DECCFLOAT:
            {
                int status;
                status = _decball_set_decfloat_other(&res->re, DECCFLOAT_REALREF((deccfloat_srcptr) x), x_ctx, ctx);
                status |= _decball_set_decfloat_other(&res->im, DECCFLOAT_IMAGREF((deccfloat_srcptr) x), x_ctx, ctx);
                return status;
            }

        case GR_CTX_DECCBALL:
            {
                int status;
                status = _decball_set_decball_other(&res->re, DECCBALL_REALREF((deccball_srcptr) x), x_ctx, ctx);
                status |= _decball_set_decball_other(&res->im, DECCBALL_IMAGREF((deccball_srcptr) x), x_ctx, ctx);
                return status;
            }

        case GR_CTX_REAL_ALGEBRAIC_QQBAR:
        case GR_CTX_COMPLEX_ALGEBRAIC_QQBAR:
            return deccball_set_qqbar(res, x, ctx);

        default:
            {
                gr_ctx_t cctx;
                acb_t z;
                int status;

                gr_ctx_init_complex_acb(cctx, 20 + 4 * PREC(ctx));
                acb_init(z);

                status = gr_set_other(z, x, x_ctx, cctx);

                if (status == GR_SUCCESS)
                    status = deccball_set_acb(res, z, ctx);

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
deccball_get_str(const deccball_t x, gr_ctx_t ctx)
{
    char * a, * b, * s;

    if (_deccball_is_real(x, ctx))
        return decball_get_str(&x->re, ctx);

    if (DECFLOAT_IS_ZERO(&x->re.mid) && DECMAG_IS_ZERO(&x->re.rad))
    {
        b = decball_get_str(&x->im, ctx);
        s = flint_malloc(strlen(b) + 3);
        strcpy(s, b);
        strcat(s, "*I");
        flint_free(b);
        return s;
    }

    a = decball_get_str(&x->re, ctx);

    if (DECMAG_IS_ZERO(&x->im.rad) && DECFLOAT_SGNBIT(&x->im.mid))
    {
        decball_t t;
        decball_init(t, ctx);
        GR_MUST_SUCCEED(decfloat_neg_round(&t->mid, &x->im.mid, DECIMAL_PREC_EXACT, EXACT_RND, ctx));
        b = decball_get_str(t, ctx);
        decball_clear(t, ctx);
        s = flint_malloc(strlen(a) + strlen(b) + 8);
        strcpy(s, "(");
        strcat(s, a);
        strcat(s, " - ");
    }
    else
    {
        b = decball_get_str(&x->im, ctx);
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
deccball_write(gr_stream_t out, const deccball_t x, gr_ctx_t ctx)
{
    return gr_stream_write_free(out, deccball_get_str(x, ctx));
}

int
deccball_randtest(deccball_t res, flint_rand_t state, gr_ctx_t ctx)
{
    int status;

    status = decball_randtest(&res->re, state, ctx);
    status |= decball_randtest(&res->im, state, ctx);

    switch (n_randint(state, 8))
    {
        case 0:
        case 1:
            decball_zero(&res->im, ctx);
            break;
        case 2:
            decball_zero(&res->re, ctx);
            break;
        default:
            break;
    }

    return status;
}

/* ------------------------------------------------------------------------- */
/*    Powers                                                                 */
/* ------------------------------------------------------------------------- */

int
deccball_pow_fmpz(deccball_t res, const deccball_t x, const fmpz_t n, gr_ctx_t ctx)
{
    if (_deccball_is_real(x, ctx))
    {
        decball_zero(&res->im, ctx);
        return decball_pow_fmpz(&res->re, &x->re, n, ctx);
    }

    if (fmpz_bits(n) < 40)
        return gr_generic_pow_fmpz(res, x, n, ctx);

    {
        deccball_t t;
        int status;
        deccball_init(t, ctx);
        status = deccball_set_fmpz(t, n, ctx);
        if (status == GR_SUCCESS)
            status = deccball_pow(res, x, t, ctx);
        deccball_clear(t, ctx);
        return status;
    }
}

int
deccball_pow_ui(deccball_t res, const deccball_t x, ulong n, gr_ctx_t ctx)
{
    fmpz_t t;
    int status;
    fmpz_init_set_ui(t, n);
    status = deccball_pow_fmpz(res, x, t, ctx);
    fmpz_clear(t);
    return status;
}

int
deccball_pow_si(deccball_t res, const deccball_t x, slong n, gr_ctx_t ctx)
{
    fmpz_t t;
    int status;
    fmpz_init_set_si(t, n);
    status = deccball_pow_fmpz(res, x, t, ctx);
    fmpz_clear(t);
    return status;
}

int
deccball_pow(deccball_t res, const deccball_t x, const deccball_t y, gr_ctx_t ctx)
{
    deccball_srcptr a[2];
    deccball_ptr r[1];

    if (_deccball_is_exact(y, ctx) && _deccball_is_real(y, ctx) && _decfloat_is_int(&y->re.mid, ctx)
        && (DECFLOAT_IS_ZERO(&y->re.mid) || _decfloat_sci_exp_clamped(&y->re.mid, ctx) < 12))
    {
        fmpz_t n;
        fmpz_init(n);
        if (decfloat_get_fmpz(n, &y->re.mid, ctx) == GR_SUCCESS && fmpz_bits(n) < 40)
        {
            int status = gr_generic_pow_fmpz(res, x, n, ctx);
            fmpz_clear(n);
            return status;
        }
        fmpz_clear(n);
    }

    if (_deccball_is_real(x, ctx) && _deccball_is_real(y, ctx))
    {
        int status = decball_pow(&res->re, &x->re, &y->re, ctx);
        if (status == GR_SUCCESS)
        {
            decball_zero(&res->im, ctx);
            return status;
        }
    }

    a[0] = x;
    a[1] = y;
    r[0] = res;
    return _deccball_acb_gr(r, 1, a, 2, 0, 0, GR_METHOD_POW, ctx);
}

/* ------------------------------------------------------------------------- */
/*    Radius operations                                                      */
/* ------------------------------------------------------------------------- */

int
deccball_add_error_decmag(deccball_t res, const decmag_t err, gr_ctx_t ctx)
{
    int status;
    status = decball_add_error_decmag(&res->re, err, ctx);
    status |= decball_add_error_decmag(&res->im, err, ctx);
    return status;
}

int
deccball_trim(deccball_t res, const deccball_t x, gr_ctx_t ctx)
{
    int status;
    status = decball_trim(&res->re, &x->re, ctx);
    status |= decball_trim(&res->im, &x->im, ctx);
    return status;
}

/* like acb_rel_accuracy_bits: the larger midpoint against the larger
   radius */
slong
deccball_rel_accuracy_digits(const deccball_t x, gr_ctx_t ctx)
{
    fmpz_t Em, Er, t;
    slong res;

    if (DECMAG_IS_ZERO(&x->re.rad) && DECMAG_IS_ZERO(&x->im.rad))
        return _deccball_is_finite(x, ctx) ? DECIMAL_PREC_EXACT : -DECIMAL_PREC_EXACT;

    if (!_deccball_is_finite(x, ctx))
        return -DECIMAL_PREC_EXACT;

    fmpz_init(Em);
    fmpz_init(Er);
    fmpz_init(t);

    /* larger radius exponent */
    if (_decmag_cmp(&x->re.rad, &x->im.rad, ctx) >= 0)
        _decmag_get_sci_exp(Er, &x->re.rad);
    else
        _decmag_get_sci_exp(Er, &x->im.rad);

    /* larger midpoint exponent */
    if (!DECFLOAT_IS_ZERO(&x->re.mid))
        decfloat_get_sci_exp(Em, &x->re.mid, ctx);
    if (!DECFLOAT_IS_ZERO(&x->im.mid))
    {
        decfloat_get_sci_exp(t, &x->im.mid, ctx);
        if (DECFLOAT_IS_ZERO(&x->re.mid) || fmpz_cmp(t, Em) > 0)
            fmpz_swap(Em, t);
    }

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
    fmpz_clear(t);
    return res;
}

/* ------------------------------------------------------------------------- */
/*    Method table                                                           */
/* ------------------------------------------------------------------------- */

static truth_t
_deccball_ctx_is_canonical(gr_ctx_t ctx)
{
    return T_FALSE;
}

int _deccball_methods_initialized = 0;
gr_static_method_table _deccball_methods;

gr_method_tab_input _deccball_methods_input[] =
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
    {GR_METHOD_CTX_IS_ALGEBRAICALLY_CLOSED, (gr_funcptr) gr_generic_ctx_predicate_true},
    {GR_METHOD_CTX_IS_ORDERED_RING, (gr_funcptr) gr_generic_ctx_predicate_false},
    {GR_METHOD_CTX_IS_RATIONAL_VECTOR_SPACE, (gr_funcptr) _decimal_ctx_is_vector_space},
    {GR_METHOD_CTX_IS_REAL_VECTOR_SPACE, (gr_funcptr) _decimal_ctx_is_vector_space},
    {GR_METHOD_CTX_IS_COMPLEX_VECTOR_SPACE, (gr_funcptr) _decimal_ctx_is_vector_space},
    {GR_METHOD_CTX_IS_EXACT,    (gr_funcptr) gr_generic_ctx_predicate_false},
    {GR_METHOD_CTX_IS_CANONICAL, (gr_funcptr) _deccball_ctx_is_canonical},
    {GR_METHOD_CTX_HAS_REAL_PREC, (gr_funcptr) gr_generic_ctx_predicate_true},
    {GR_METHOD_CTX_SET_REAL_PREC, (gr_funcptr) _decimal_ctx_set_real_prec},
    {GR_METHOD_CTX_GET_REAL_PREC, (gr_funcptr) _decimal_ctx_get_real_prec},

    {GR_METHOD_INIT,            (gr_funcptr) deccball_init},
    {GR_METHOD_CLEAR,           (gr_funcptr) deccball_clear},
    {GR_METHOD_SWAP,            (gr_funcptr) deccball_swap},
    {GR_METHOD_SET_SHALLOW,     (gr_funcptr) deccball_set_shallow},
    {GR_METHOD_RANDTEST,        (gr_funcptr) deccball_randtest},
    {GR_METHOD_WRITE,           (gr_funcptr) deccball_write},
    {GR_METHOD_ZERO,            (gr_funcptr) deccball_zero},
    {GR_METHOD_ONE,             (gr_funcptr) deccball_one},
    {GR_METHOD_NEG_ONE,         (gr_funcptr) deccball_neg_one},
    {GR_METHOD_IS_ZERO,         (gr_funcptr) deccball_is_zero},
    {GR_METHOD_IS_ONE,          (gr_funcptr) deccball_is_one},
    {GR_METHOD_IS_NEG_ONE,      (gr_funcptr) deccball_is_neg_one},
    {GR_METHOD_EQUAL,           (gr_funcptr) deccball_equal},
    {GR_METHOD_SET,             (gr_funcptr) deccball_set},
    {GR_METHOD_SET_SI,          (gr_funcptr) deccball_set_si},
    {GR_METHOD_SET_UI,          (gr_funcptr) deccball_set_ui},
    {GR_METHOD_SET_FMPZ,        (gr_funcptr) deccball_set_fmpz},
    {GR_METHOD_SET_FMPQ,        (gr_funcptr) deccball_set_fmpq},
    {GR_METHOD_SET_D,           (gr_funcptr) deccball_set_d},
    {GR_METHOD_SET_STR,         (gr_funcptr) deccball_set_str},
    {GR_METHOD_SET_FMPZ_10EXP_FMPZ, (gr_funcptr) deccball_set_fmpz_10exp_fmpz},
    {GR_METHOD_SET_FMPZ_2EXP_FMPZ, (gr_funcptr) deccball_set_fmpz_2exp_fmpz},
    {GR_METHOD_SET_INTERVAL_MID_RAD, (gr_funcptr) deccball_set_interval_mid_rad},
    {GR_METHOD_SET_INTERVAL_MID_INF, (gr_funcptr) deccball_set_interval_mid_inf},
    {GR_METHOD_SET_OTHER,       (gr_funcptr) deccball_set_other},
    {GR_METHOD_GET_FMPZ,        (gr_funcptr) deccball_get_fmpz},
    {GR_METHOD_GET_FMPQ,        (gr_funcptr) deccball_get_fmpq},
    {GR_METHOD_GET_SI,          (gr_funcptr) deccball_get_si},
    {GR_METHOD_GET_UI,          (gr_funcptr) deccball_get_ui},
    {GR_METHOD_GET_D,           (gr_funcptr) deccball_get_d},

    {GR_METHOD_NEG,             (gr_funcptr) deccball_neg},
    {GR_METHOD_ADD,             (gr_funcptr) deccball_add},
    {GR_METHOD_ADD_UI,          (gr_funcptr) deccball_add_ui},
    {GR_METHOD_ADD_SI,          (gr_funcptr) deccball_add_si},
    {GR_METHOD_ADD_FMPZ,        (gr_funcptr) deccball_add_fmpz},
    {GR_METHOD_SUB,             (gr_funcptr) deccball_sub},
    {GR_METHOD_SUB_UI,          (gr_funcptr) deccball_sub_ui},
    {GR_METHOD_SUB_SI,          (gr_funcptr) deccball_sub_si},
    {GR_METHOD_SUB_FMPZ,        (gr_funcptr) deccball_sub_fmpz},
    {GR_METHOD_MUL,             (gr_funcptr) deccball_mul},
    {GR_METHOD_MUL_UI,          (gr_funcptr) deccball_mul_ui},
    {GR_METHOD_MUL_SI,          (gr_funcptr) deccball_mul_si},
    {GR_METHOD_MUL_FMPZ,        (gr_funcptr) deccball_mul_fmpz},
    {GR_METHOD_MUL_TWO,         (gr_funcptr) deccball_mul_two},
    {GR_METHOD_SQR,             (gr_funcptr) deccball_sqr},
    {GR_METHOD_DIV,             (gr_funcptr) deccball_div},
    {GR_METHOD_DIV_UI,          (gr_funcptr) deccball_div_ui},
    {GR_METHOD_DIV_SI,          (gr_funcptr) deccball_div_si},
    {GR_METHOD_DIV_FMPZ,        (gr_funcptr) deccball_div_fmpz},
    {GR_METHOD_INV,             (gr_funcptr) deccball_inv},
    {GR_METHOD_MUL_2EXP_SI,     (gr_funcptr) deccball_mul_2exp_si},
    {GR_METHOD_MUL_2EXP_FMPZ,   (gr_funcptr) deccball_mul_2exp_fmpz},
    {GR_METHOD_POW_UI,          (gr_funcptr) deccball_pow_ui},
    {GR_METHOD_POW_SI,          (gr_funcptr) deccball_pow_si},
    {GR_METHOD_POW_FMPZ,        (gr_funcptr) deccball_pow_fmpz},
    {GR_METHOD_POW,             (gr_funcptr) deccball_pow},
    {GR_METHOD_SQRT,            (gr_funcptr) deccball_sqrt},
    {GR_METHOD_RSQRT,           (gr_funcptr) deccball_rsqrt},
    {GR_METHOD_UINF,            (gr_funcptr) gr_not_in_domain},
    {GR_METHOD_FLOOR,           (gr_funcptr) deccball_floor},
    {GR_METHOD_CEIL,            (gr_funcptr) deccball_ceil},
    {GR_METHOD_TRUNC,           (gr_funcptr) deccball_trunc},
    {GR_METHOD_NINT,            (gr_funcptr) deccball_nint},
    {GR_METHOD_I,               (gr_funcptr) deccball_i},
    {GR_METHOD_ABS,             (gr_funcptr) deccball_abs},
    {GR_METHOD_CONJ,            (gr_funcptr) deccball_conj},
    {GR_METHOD_RE,              (gr_funcptr) deccball_re},
    {GR_METHOD_IM,              (gr_funcptr) deccball_im},
    {GR_METHOD_SGN,             (gr_funcptr) deccball_sgn},
    {GR_METHOD_CSGN,            (gr_funcptr) deccball_csgn},
    /* balls represent complex numbers: no infinities or undefined values */
    {GR_METHOD_POS_INF,         (gr_funcptr) gr_not_in_domain},
    {GR_METHOD_NEG_INF,         (gr_funcptr) gr_not_in_domain},
    {GR_METHOD_UINF,            (gr_funcptr) gr_not_in_domain},
    {GR_METHOD_UNDEFINED,       (gr_funcptr) gr_not_in_domain},
    {GR_METHOD_UNKNOWN,         (gr_funcptr) gr_not_in_domain},
    {GR_METHOD_ARG,             (gr_funcptr) deccball_arg},
    {GR_METHOD_CMP,             (gr_funcptr) deccball_cmp},
    {GR_METHOD_CMPABS,          (gr_funcptr) deccball_cmpabs},
    {GR_METHOD_MAT_DET,         (gr_funcptr) gr_mat_det_generic_field},
    {GR_METHOD_MAT_FIND_NONZERO_PIVOT, (gr_funcptr) gr_mat_find_nonzero_pivot_large_abs},
    {0,                         (gr_funcptr) NULL},
};

POP_OPTIONS
