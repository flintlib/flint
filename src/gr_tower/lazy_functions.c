/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/*
    Elementary functions of the lazy tower field, expressed through
    exp, log, sqrt, pi and i (which the field represents exactly), and
    real-number operations (floor, comparisons) decided by combining
    enclosures with exact zero tests. Everything here goes through the
    public gr interface of the context.

    Inverse functions are computed by the logarithmic formulas defining
    their principal branches (with the same conventions on the branch
    cuts as Arb), and the result is checked numerically against the
    principal value as a safeguard (GR_UNABLE if the two disagree).
*/

#include "fmpq.h"
#include "acb.h"
#include "gr.h"
#include "gr_special.h"
#include "gr_tower.h"
#include "gr_tower_lazy.h"
#include "gr_tower/impl.h"

PUSH_OPTIONS
OPTIMIZE_OSIZE

#define LAZY_CHECK_PREC 64

/* Whether the enclosures of x and y (at a modest precision) overlap. */
static int
_gr_tower_lazy_overlaps(const gr_tower_lazy_elem_t x, const acb_t y, gr_ctx_t ctx)
{
    acb_t z;
    int res;
    acb_init(z);
    res = (gr_tower_lazy_get_acb(z, x, LAZY_CHECK_PREC, ctx) == GR_SUCCESS) && acb_overlaps(z, y);
    acb_clear(z);
    return res;
}

/* Multiplies by i, divides by i, and the constants pi/2 and 2i. */
static int _gr_tower_lazy_mul_i(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t x, gr_ctx_t ctx)
{
    gr_ptr i;
    int status;
    GR_TMP_INIT(i, ctx);
    status = gr_i(i, ctx);
    status |= gr_mul(res, x, i, ctx);
    GR_TMP_CLEAR(i, ctx);
    return status;
}

static int _gr_tower_lazy_div_i(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t x, gr_ctx_t ctx)
{
    int status = _gr_tower_lazy_mul_i(res, x, ctx);
    return status | gr_neg(res, res, ctx);
}

/* -------------------------------------------------------------------- */
/* trigonometric and hyperbolic functions                                */
/* -------------------------------------------------------------------- */

/* e = exp(x), f = exp(-x) = 1/e */
static int
_exp_pair(gr_ptr e, gr_ptr f, gr_srcptr x, gr_ctx_t ctx)
{
    int status = gr_exp(e, x, ctx);
    if (status == GR_SUCCESS)
        status = gr_inv(f, e, ctx);
    return status;
}

static int
_gr_tower_lazy_sinh_impl(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t x, gr_ctx_t ctx)
{
    gr_ptr e, f;
    int status;
    GR_TMP_INIT2(e, f, ctx);
    status = _exp_pair(e, f, x, ctx);
    status |= gr_sub(res, e, f, ctx);
    status |= gr_div_ui(res, res, 2, ctx);
    GR_TMP_CLEAR2(e, f, ctx);
    return status;
}

static int
_gr_tower_lazy_cosh_impl(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t x, gr_ctx_t ctx)
{
    gr_ptr e, f;
    int status;
    GR_TMP_INIT2(e, f, ctx);
    status = _exp_pair(e, f, x, ctx);
    status |= gr_add(res, e, f, ctx);
    status |= gr_div_ui(res, res, 2, ctx);
    GR_TMP_CLEAR2(e, f, ctx);
    return status;
}

static int
_gr_tower_lazy_tanh_impl(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t x, gr_ctx_t ctx)
{
    gr_ptr e, f, s;
    int status;
    GR_TMP_INIT3(e, f, s, ctx);
    status = _exp_pair(e, f, x, ctx);
    status |= gr_sub(s, e, f, ctx);
    status |= gr_add(e, e, f, ctx);
    if (status == GR_SUCCESS)
        status = gr_div(res, s, e, ctx);
    GR_TMP_CLEAR3(e, f, s, ctx);
    return status;
}

static int
_gr_tower_lazy_sin_impl(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t x, gr_ctx_t ctx)
{
    /* sin(x) = sinh(i x) / i */
    gr_ptr t;
    int status;
    GR_TMP_INIT(t, ctx);
    status = _gr_tower_lazy_mul_i(t, x, ctx);
    status |= gr_tower_lazy_sinh(t, t, ctx);
    status |= _gr_tower_lazy_div_i(res, t, ctx);
    GR_TMP_CLEAR(t, ctx);
    return status;
}

static int
_gr_tower_lazy_cos_impl(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t x, gr_ctx_t ctx)
{
    gr_ptr t;
    int status;
    GR_TMP_INIT(t, ctx);
    status = _gr_tower_lazy_mul_i(t, x, ctx);
    status |= gr_tower_lazy_cosh(res, t, ctx);
    GR_TMP_CLEAR(t, ctx);
    return status;
}

static int
_gr_tower_lazy_tan_impl(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t x, gr_ctx_t ctx)
{
    /* tan(x) = tanh(i x) / i */
    gr_ptr t;
    int status;
    GR_TMP_INIT(t, ctx);
    status = _gr_tower_lazy_mul_i(t, x, ctx);
    status |= gr_tower_lazy_tanh(t, t, ctx);
    status |= _gr_tower_lazy_div_i(res, t, ctx);
    GR_TMP_CLEAR(t, ctx);
    return status;
}

/* -------------------------------------------------------------------- */
/* inverse functions                                                     */
/* -------------------------------------------------------------------- */

/* Accepts res if it overlaps the enclosure of the principal value. */
static int
_check_principal(gr_ptr res, const acb_t principal, gr_ctx_t ctx)
{
    if (!acb_is_finite(principal))
        return GR_UNABLE;
    return _gr_tower_lazy_overlaps(res, principal, ctx) ? GR_SUCCESS : GR_UNABLE;
}

typedef void (*acb_unary_func)(acb_t, const acb_t, slong);

static int
_principal_value(acb_t p, gr_srcptr x, acb_unary_func f, gr_ctx_t ctx)
{
    acb_t z;
    int status;
    acb_init(z);
    status = gr_tower_lazy_get_acb(z, x, LAZY_CHECK_PREC, ctx);
    if (status == GR_SUCCESS)
        f(p, z, LAZY_CHECK_PREC);
    acb_clear(z);
    return status;
}

static int
_gr_tower_lazy_atanh_impl(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t x, gr_ctx_t ctx)
{
    /* (log(1+x) - log(1-x)) / 2 */
    gr_ptr a, b;
    acb_t p;
    int status;

    GR_TMP_INIT2(a, b, ctx);
    acb_init(p);
    status = _principal_value(p, x, acb_atanh, ctx);
    status |= gr_add_ui(a, x, 1, ctx);
    status |= gr_log(a, a, ctx);
    status |= gr_neg(b, x, ctx);
    status |= gr_add_ui(b, b, 1, ctx);
    status |= gr_log(b, b, ctx);
    status |= gr_sub(res, a, b, ctx);
    status |= gr_div_ui(res, res, 2, ctx);
    if (status == GR_SUCCESS)
        status = _check_principal(res, p, ctx);
    acb_clear(p);
    GR_TMP_CLEAR2(a, b, ctx);
    return status;
}

static int
_gr_tower_lazy_atan_impl(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t x, gr_ctx_t ctx)
{
    /* atan(x) = atanh(i x) / i */
    gr_ptr t;
    acb_t p;
    int status;

    GR_TMP_INIT(t, ctx);
    acb_init(p);
    status = _principal_value(p, x, acb_atan, ctx);
    status |= _gr_tower_lazy_mul_i(t, x, ctx);
    if (status == GR_SUCCESS)
    {
        /* the formula directly, the branch being checked here */
        gr_ptr a, b;
        GR_TMP_INIT2(a, b, ctx);
        status |= gr_add_ui(a, t, 1, ctx);
        status |= gr_log(a, a, ctx);
        status |= gr_neg(b, t, ctx);
        status |= gr_add_ui(b, b, 1, ctx);
        status |= gr_log(b, b, ctx);
        status |= gr_sub(res, a, b, ctx);
        status |= gr_div_ui(res, res, 2, ctx);
        status |= _gr_tower_lazy_div_i(res, res, ctx);
        GR_TMP_CLEAR2(a, b, ctx);
    }
    if (status == GR_SUCCESS)
        status = _check_principal(res, p, ctx);
    acb_clear(p);
    GR_TMP_CLEAR(t, ctx);
    return status;
}

static int
_gr_tower_lazy_asinh_impl(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t x, gr_ctx_t ctx)
{
    /* log(x + sqrt(x^2 + 1)) */
    gr_ptr t;
    acb_t p;
    int status;

    GR_TMP_INIT(t, ctx);
    acb_init(p);
    status = _principal_value(p, x, acb_asinh, ctx);
    status |= gr_sqr(t, x, ctx);
    status |= gr_add_ui(t, t, 1, ctx);
    status |= gr_sqrt(t, t, ctx);
    status |= gr_add(t, t, x, ctx);
    status |= gr_log(res, t, ctx);
    if (status == GR_SUCCESS)
        status = _check_principal(res, p, ctx);
    acb_clear(p);
    GR_TMP_CLEAR(t, ctx);
    return status;
}

static int
_gr_tower_lazy_asin_impl(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t x, gr_ctx_t ctx)
{
    /* asin(x) = asinh(i x) / i */
    gr_ptr t;
    acb_t p;
    int status;

    GR_TMP_INIT(t, ctx);
    acb_init(p);
    status = _principal_value(p, x, acb_asin, ctx);
    status |= _gr_tower_lazy_mul_i(t, x, ctx);
    if (status == GR_SUCCESS)
    {
        gr_ptr u;
        GR_TMP_INIT(u, ctx);
        status |= gr_sqr(u, t, ctx);
        status |= gr_add_ui(u, u, 1, ctx);
        status |= gr_sqrt(u, u, ctx);
        status |= gr_add(u, u, t, ctx);
        status |= gr_log(res, u, ctx);
        status |= _gr_tower_lazy_div_i(res, res, ctx);
        GR_TMP_CLEAR(u, ctx);
    }
    if (status == GR_SUCCESS)
        status = _check_principal(res, p, ctx);
    acb_clear(p);
    GR_TMP_CLEAR(t, ctx);
    return status;
}

static int
_gr_tower_lazy_acos_impl(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t x, gr_ctx_t ctx)
{
    /* pi/2 - asin(x) */
    gr_ptr t;
    acb_t p;
    int status;

    GR_TMP_INIT(t, ctx);
    acb_init(p);
    status = _principal_value(p, x, acb_acos, ctx);
    status |= gr_tower_lazy_asin(t, x, ctx);
    status |= gr_pi(res, ctx);
    status |= gr_div_ui(res, res, 2, ctx);
    status |= gr_sub(res, res, t, ctx);
    if (status == GR_SUCCESS)
        status = _check_principal(res, p, ctx);
    acb_clear(p);
    GR_TMP_CLEAR(t, ctx);
    return status;
}

static int
_gr_tower_lazy_acosh_impl(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t x, gr_ctx_t ctx)
{
    /* log(x + sqrt(x+1) sqrt(x-1)) */
    gr_ptr t, u;
    acb_t p;
    int status;

    GR_TMP_INIT2(t, u, ctx);
    acb_init(p);
    status = _principal_value(p, x, acb_acosh, ctx);
    status |= gr_add_ui(t, x, 1, ctx);
    status |= gr_sqrt(t, t, ctx);
    status |= gr_sub_ui(u, x, 1, ctx);
    status |= gr_sqrt(u, u, ctx);
    status |= gr_mul(t, t, u, ctx);
    status |= gr_add(t, t, x, ctx);
    status |= gr_log(res, t, ctx);
    if (status == GR_SUCCESS)
        status = _check_principal(res, p, ctx);
    acb_clear(p);
    GR_TMP_CLEAR2(t, u, ctx);
    return status;
}

/* -------------------------------------------------------------------- */
/* real and imaginary parts, argument, sign                              */
/* -------------------------------------------------------------------- */

static int
_gr_tower_lazy_re_impl(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t x, gr_ctx_t ctx)
{
    gr_ptr c;
    int status;
    GR_TMP_INIT(c, ctx);
    status = gr_conj(c, x, ctx);
    status |= gr_add(res, x, c, ctx);
    status |= gr_div_ui(res, res, 2, ctx);
    GR_TMP_CLEAR(c, ctx);
    return status;
}

static int
_gr_tower_lazy_im_impl(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t x, gr_ctx_t ctx)
{
    gr_ptr c;
    int status;
    GR_TMP_INIT(c, ctx);
    status = gr_conj(c, x, ctx);
    status |= gr_sub(res, x, c, ctx);
    status |= gr_div_ui(res, res, 2, ctx);
    status |= _gr_tower_lazy_div_i(res, res, ctx);
    GR_TMP_CLEAR(c, ctx);
    return status;
}

static int
_gr_tower_lazy_arg_impl(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t x, gr_ctx_t ctx)
{
    /* arg(x) = im(log(x)) for x != 0; arg(0) = 0 */
    truth_t z = gr_is_zero(x, ctx);
    gr_ptr t;
    int status;

    if (z == T_TRUE)
        return gr_zero(res, ctx);
    if (z == T_UNKNOWN)
        return GR_UNABLE;

    GR_TMP_INIT(t, ctx);
    status = gr_log(t, x, ctx);
    if (status == GR_SUCCESS)
        status = gr_tower_lazy_im(res, t, ctx);
    GR_TMP_CLEAR(t, ctx);
    return status;
}

static int
_gr_tower_lazy_sgn_impl(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t x, gr_ctx_t ctx)
{
    /* x / |x| */
    truth_t z = gr_is_zero(x, ctx);
    gr_ptr t;
    int status;

    if (z == T_TRUE)
        return gr_zero(res, ctx);
    if (z == T_UNKNOWN)
        return GR_UNABLE;

    GR_TMP_INIT(t, ctx);
    status = gr_abs(t, x, ctx);
    if (status == GR_SUCCESS)
        status = gr_div(res, x, t, ctx);
    GR_TMP_CLEAR(t, ctx);
    return status;
}

/* -------------------------------------------------------------------- */
/* real numbers: comparisons and rounding                                */
/* -------------------------------------------------------------------- */

/*
    Whether x is real: its enclosure has an exactly zero imaginary part
    (which is the case for elements of real generators), or x equals its
    conjugate.
*/
static truth_t
_gr_tower_lazy_is_real_impl(const gr_tower_lazy_elem_t x, gr_ctx_t ctx)
{
    acb_t z;
    truth_t res = T_UNKNOWN;
    gr_ptr c;

    acb_init(z);
    if (gr_tower_lazy_get_acb(z, x, LAZY_CHECK_PREC, ctx) == GR_SUCCESS)
    {
        if (arb_is_zero(acb_imagref(z)))
            res = T_TRUE;
        else if (!arb_contains_zero(acb_imagref(z)))
            res = T_FALSE;
    }
    acb_clear(z);

    if (res != T_UNKNOWN)
        return res;

    GR_TMP_INIT(c, ctx);
    if (gr_conj(c, x, ctx) == GR_SUCCESS)
        res = gr_equal(c, x, ctx);
    GR_TMP_CLEAR(c, ctx);
    return res;
}

/*
    The sign of the real number x: -1, 0 or 1; GR_DOMAIN if x is not
    known to be real, GR_UNABLE if the sign cannot be decided.
*/
static int
_gr_tower_lazy_real_sign_impl(int * sign, const gr_tower_lazy_elem_t x, gr_ctx_t ctx)
{
    truth_t z;
    acb_t v;
    slong prec;
    int status = GR_UNABLE, tested = 0, unknown = 0;

    z = gr_tower_lazy_is_real(x, ctx);
    if (z != T_TRUE)
        return (z == T_FALSE) ? GR_DOMAIN : GR_UNABLE;

    /* numerically first (an exact zero test can be far more expensive
       than a few more bits), the exact test only when the enclosures
       keep containing zero */
    acb_init(v);
    /* (up to the numeric precision limit when the exact test is
       undecided; a provably nonzero x is separated eventually, with a
       generous safety cap) */
    for (prec = LAZY_CHECK_PREC; prec <= gr_tower_lazy_ctx_get_option(ctx, GR_TOWER_OPT_NUMERIC_PREC_LIMIT) * (unknown ? 1 : 16); prec *= 2)
    {
        if (gr_tower_lazy_get_acb(v, x, prec, ctx) != GR_SUCCESS)
            break;
        if (arb_is_positive(acb_realref(v)))
        {
            *sign = 1;
            status = GR_SUCCESS;
            break;
        }
        if (arb_is_negative(acb_realref(v)))
        {
            *sign = -1;
            status = GR_SUCCESS;
            break;
        }
        if (prec >= 8 * LAZY_CHECK_PREC && !tested)
        {
            tested = 1;
            z = gr_is_zero(x, ctx);
            /* not provably zero or nonzero (no structure theorem covers
               it): only more precision can still separate x from zero,
               as for 1 - erf(100) against exp(-10^4) (14500 bits) */
            if (z == T_UNKNOWN)
            {
                unknown = 1;
                continue;
            }
            if (z == T_TRUE)
            {
                *sign = 0;
                status = GR_SUCCESS;
                break;
            }
        }
    }
    acb_clear(v);
    return status;
}

static int
_gr_tower_lazy_cmp_impl(int * res, const gr_tower_lazy_elem_t x, const gr_tower_lazy_elem_t y, gr_ctx_t ctx)
{
    gr_ptr d;
    int status;

    /* two real numbers whose enclosures are disjoint: ordered without
       forming x - y (which merges their towers) */
    if (gr_tower_lazy_is_real(x, ctx) == T_TRUE && gr_tower_lazy_is_real(y, ctx) == T_TRUE)
    {
        acb_t a, b;
        slong prec;
        int done = 0;
        acb_init(a);
        acb_init(b);
        for (prec = LAZY_CHECK_PREC; prec <= gr_tower_lazy_ctx_get_option(ctx, GR_TOWER_OPT_PREC_LIMIT) && !done; prec *= 2)
        {
            if (gr_tower_lazy_get_acb(a, x, prec, ctx) != GR_SUCCESS || gr_tower_lazy_get_acb(b, y, prec, ctx) != GR_SUCCESS)
                break;
            if (arb_lt(acb_realref(a), acb_realref(b)))
                *res = -1, done = 1;
            else if (arb_gt(acb_realref(a), acb_realref(b)))
                *res = 1, done = 1;
        }
        acb_clear(a);
        acb_clear(b);
        if (done)
            return GR_SUCCESS;
    }

    GR_TMP_INIT(d, ctx);
    status = gr_sub(d, x, y, ctx);
    if (status == GR_SUCCESS)
        status = _gr_tower_lazy_real_sign_locked(res, d, ctx);
    GR_TMP_CLEAR(d, ctx);
    return status;
}

/* csgn(x): the sign of Re(x), or of Im(x) when Re(x) = 0 */
static int
_gr_tower_lazy_csgn_impl(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t x, gr_ctx_t ctx)
{
    gr_ptr t;
    int status, sign = 0;
    GR_TMP_INIT(t, ctx);
    status = gr_tower_lazy_re(t, x, ctx);
    if (status == GR_SUCCESS)
        status = _gr_tower_lazy_real_sign_locked(&sign, t, ctx);
    if (status == GR_SUCCESS && sign == 0)
    {
        status = gr_tower_lazy_im(t, x, ctx);
        if (status == GR_SUCCESS)
            status = _gr_tower_lazy_real_sign_locked(&sign, t, ctx);
    }
    if (status == GR_SUCCESS)
        status = gr_set_si(res, sign, ctx);
    GR_TMP_CLEAR(t, ctx);
    return status;
}

/* the sign of |x| - |y| = sign of x conj(x) - y conj(y) */
static int
_gr_tower_lazy_cmpabs_impl(int * res, const gr_tower_lazy_elem_t x, const gr_tower_lazy_elem_t y, gr_ctx_t ctx)
{
    gr_ptr d, t;
    int status;
    GR_TMP_INIT2(d, t, ctx);
    status = gr_conj(d, x, ctx);
    status |= gr_mul(d, d, x, ctx);
    status |= gr_conj(t, y, ctx);
    status |= gr_mul(t, t, y, ctx);
    status |= gr_sub(d, d, t, ctx);
    if (status == GR_SUCCESS)
        status = _gr_tower_lazy_real_sign_locked(res, d, ctx);
    GR_TMP_CLEAR2(d, t, ctx);
    return status;
}

/*
    floor(x) for real x: the candidate n = floor(midpoint) is verified by
    the exact sign of x - n and x - (n + 1).
*/
static int _gr_tower_lazy_floor_real(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t x, gr_ctx_t ctx);

/* the real part of x when x is not known to be real (the rounding
   functions act on the real part, as for qqbar and acb) */
static int
_gr_tower_lazy_real_part(gr_tower_lazy_elem_t res, int * is_real, const gr_tower_lazy_elem_t x, gr_ctx_t ctx)
{
    int status = GR_SUCCESS;
    *is_real = (gr_tower_lazy_is_real(x, ctx) == T_TRUE);
    if (!*is_real)
    {
        status = gr_tower_lazy_re(res, x, ctx);
        if (status == GR_SUCCESS && gr_tower_lazy_is_real(res, ctx) != T_TRUE)
            status = GR_UNABLE;
    }
    return status;
}

static int
_gr_tower_lazy_floor_impl(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t x, gr_ctx_t ctx)
{
    gr_ptr t;
    int status, real;

    GR_TMP_INIT(t, ctx);
    status = _gr_tower_lazy_real_part(t, &real, x, ctx);
    if (status == GR_SUCCESS)
        status = _gr_tower_lazy_floor_real(res, real ? x : t, ctx);
    GR_TMP_CLEAR(t, ctx);
    return status;
}

static int
_gr_tower_lazy_floor_real(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t x, gr_ctx_t ctx)
{
    acb_t v;
    fmpz_t n;
    gr_ptr t;
    slong prec;
    int status = GR_UNABLE, s;

    acb_init(v);
    fmpz_init(n);
    GR_TMP_INIT(t, ctx);

    for (prec = LAZY_CHECK_PREC; prec <= gr_tower_lazy_ctx_get_option(ctx, GR_TOWER_OPT_CERTIFY_PREC_LIMIT) * 4; prec *= 2)
    {
        arf_t m;
        if (gr_tower_lazy_get_acb(v, x, prec, ctx) != GR_SUCCESS)
            break;
        arf_init(m);
        arf_floor(m, arb_midref(acb_realref(v)));
        if (!arf_is_finite(m) || arf_cmpabs_2exp_si(m, 10000) > 0)
        {
            arf_clear(m);
            break;
        }
        arf_get_fmpz(n, m, ARF_RND_FLOOR);
        arf_clear(m);

        /* n <= x < n + 1 ? */
        status = gr_sub_fmpz(t, x, n, ctx);
        if (status != GR_SUCCESS)
            break;
        status = _gr_tower_lazy_real_sign_locked(&s, t, ctx);
        if (status != GR_SUCCESS)
            break;
        if (s < 0)
        {
            status = GR_UNABLE;
            continue;   /* try a better enclosure */
        }
        status = gr_sub_ui(t, t, 1, ctx);
        if (status != GR_SUCCESS)
            break;
        status = _gr_tower_lazy_real_sign_locked(&s, t, ctx);
        if (status != GR_SUCCESS)
            break;
        if (s < 0)
        {
            status = gr_set_fmpz(res, n, ctx);
            break;
        }
        status = GR_UNABLE;
    }

    acb_clear(v);
    fmpz_clear(n);
    GR_TMP_CLEAR(t, ctx);
    return status;
}

static int
_gr_tower_lazy_ceil_impl(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t x_in, gr_ctx_t ctx)
{
    /* ceil(x) = -floor(-x) (of the real part) */
    gr_ptr t;
    int status, real;
    GR_TMP_INIT(t, ctx);
    status = _gr_tower_lazy_real_part(t, &real, x_in, ctx);
    if (status == GR_SUCCESS)
    {
        status = gr_neg(t, real ? x_in : t, ctx);
        status |= _gr_tower_lazy_floor_real(t, t, ctx);
        status |= gr_neg(res, t, ctx);
    }
    GR_TMP_CLEAR(t, ctx);
    return status;
}

static int
_gr_tower_lazy_trunc_impl(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t x_in, gr_ctx_t ctx)
{
    gr_ptr t;
    int sign, status, real;
    GR_TMP_INIT(t, ctx);
    status = _gr_tower_lazy_real_part(t, &real, x_in, ctx);
    if (status == GR_SUCCESS)
        status = _gr_tower_lazy_real_sign_locked(&sign, real ? x_in : t, ctx);
    if (status == GR_SUCCESS)
    {
        if (sign >= 0)
            status = gr_tower_lazy_floor(res, real ? x_in : t, ctx);
        else
            status = gr_tower_lazy_ceil(res, real ? x_in : t, ctx);
    }
    GR_TMP_CLEAR(t, ctx);
    return status;
}

static int _gr_tower_lazy_nint_real(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t x, gr_ctx_t ctx);

static int
_gr_tower_lazy_nint_impl(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t x_in, gr_ctx_t ctx)
{
    gr_ptr t;
    int status, real;
    GR_TMP_INIT(t, ctx);
    status = _gr_tower_lazy_real_part(t, &real, x_in, ctx);
    if (status == GR_SUCCESS)
        status = _gr_tower_lazy_nint_real(res, real ? x_in : t, ctx);
    GR_TMP_CLEAR(t, ctx);
    return status;
}

static int
_gr_tower_lazy_nint_real(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t x, gr_ctx_t ctx)
{
    /* floor(x + 1/2), rounding half to even */
    gr_ptr t, u;
    int status, sign;
    fmpq_t half;

    GR_TMP_INIT2(t, u, ctx);
    fmpq_init(half);
    fmpq_set_si(half, 1, 2);
    status = gr_add_fmpq(t, x, half, ctx);
    status |= gr_tower_lazy_floor(u, t, ctx);
    if (status == GR_SUCCESS)
    {
        /* exactly halfway? then round to even */
        status = gr_sub(t, t, u, ctx);
        status |= _gr_tower_lazy_real_sign_locked(&sign, t, ctx);
        if (status == GR_SUCCESS && sign == 0)
        {
            fmpz_t n;
            fmpz_init(n);
            status = gr_get_fmpz(n, u, ctx);
            if (status == GR_SUCCESS && fmpz_is_odd(n))
                status = gr_sub_ui(u, u, 1, ctx);
            fmpz_clear(n);
        }
    }
    if (status == GR_SUCCESS)
        status = gr_set(res, u, ctx);
    fmpq_clear(half);
    GR_TMP_CLEAR2(t, u, ctx);
    return status;
}

/* A double approximation of a real element (the midpoint of an
   enclosure); GR_DOMAIN for elements not known to be real. */
static int
_gr_tower_lazy_get_d_impl(double * res, const gr_tower_lazy_elem_t x, gr_ctx_t ctx)
{
    acb_t z;
    int status;

    if (gr_tower_lazy_is_real(x, ctx) != T_TRUE)
        return GR_DOMAIN;

    acb_init(z);
    status = gr_tower_lazy_get_acb(z, x, 64, ctx);
    if (status == GR_SUCCESS)
        *res = arf_get_d(arb_midref(acb_realref(z)), ARF_RND_NEAR);
    acb_clear(z);
    return status;
}

/* locked entry points (see lazy.c), with the restrictions of the real
   and algebraic subfields: an elementary function of an algebraic
   number is transcendental except at its special values (Lindemann-
   Weierstrass, Gelfond-Schneider), and the inverse functions leave the
   real line outside their real domains */
#define LOCKED(call) int status; _gr_tower_lazy_lock(ctx); status = call; _gr_tower_lazy_unlock(ctx); return status;

/* in an algebraic context: f(x) with f(0) = value at 0 (0 or 1), else GR_DOMAIN */
static int
_alg_special_at_zero(gr_ptr res, gr_srcptr x, int one, gr_ctx_t ctx)
{
    truth_t z = gr_is_zero(x, ctx);
    if (z == T_TRUE)
        return one ? gr_one(res, ctx) : gr_zero(res, ctx);
    return (z == T_FALSE) ? GR_DOMAIN : GR_UNABLE;
}

/* in an algebraic context: f(x) with f(1) = 0, else GR_DOMAIN */
static int
_alg_special_at_one(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)
{
    truth_t o = gr_is_one(x, ctx);
    if (o == T_TRUE)
        return gr_zero(res, ctx);
    return (o == T_FALSE) ? GR_DOMAIN : GR_UNABLE;
}

/* in a real context: whether lo <= x <= hi (strictly if strict) */
static int
_real_in_interval(gr_srcptr x, slong lo, slong hi, int strict, gr_ctx_t ctx)
{
    gr_ptr t;
    int c, status;
    GR_TMP_INIT(t, ctx);
    status = gr_set_si(t, lo, ctx);
    status |= _gr_tower_lazy_cmp_impl(&c, x, t, ctx);
    if (status == GR_SUCCESS && (c < 0 || (strict && c == 0)))
        status = GR_DOMAIN;
    if (status == GR_SUCCESS)
    {
        status = gr_set_si(t, hi, ctx);
        status |= _gr_tower_lazy_cmp_impl(&c, x, t, ctx);
        if (status == GR_SUCCESS && (c > 0 || (strict && c == 0)))
            status = GR_DOMAIN;
    }
    GR_TMP_CLEAR(t, ctx);
    return status;
}

#define FLAGS(ctx) gr_tower_lazy_ctx_field_flags(ctx)

/* -------------------------------------------------------------------- */
/* real forms: tan and atan generators                                   */
/* -------------------------------------------------------------------- */

/* In a real field, whether x is real (then the real forms are used). */
static int
_real_arg(gr_srcptr x, gr_ctx_t ctx)
{
    return (FLAGS(ctx) & GR_TOWER_LAZY_REAL) && _gr_tower_lazy_is_real_exact(x, ctx) == T_TRUE;
}

/* Whether sin, cos and tan of x take the real forms: in a real field, or
   with the option GR_TOWER_OPT_TRIG_FORM = GR_TOWER_TRIG_TANGENT, for real x. */
static int
_trig_real_arg(gr_srcptr x, gr_ctx_t ctx)
{
    if (FLAGS(ctx) & GR_TOWER_LAZY_REAL)
        return _real_arg(x, ctx);
    return gr_tower_lazy_ctx_get_option(ctx, GR_TOWER_OPT_TRIG_FORM) == GR_TOWER_TRIG_TANGENT &&
           _gr_tower_lazy_is_real_exact(x, ctx) == T_TRUE;
}

/* Whether x = r pi with r rational. */
static int
_pi_multiple(fmpq_t r, gr_srcptr x, gr_ctx_t ctx)
{
    gr_ptr t;
    int found = 0;
    acb_t z;

    /* (quick numerical exclusion: x / pi far from a rational with a
       small denominator is not tested exactly -- the exact test needs
       only the division, but is not free) */
    acb_init(z);
    GR_TMP_INIT(t, ctx);
    if (gr_pi(t, ctx) == GR_SUCCESS && gr_div(t, x, t, ctx) == GR_SUCCESS &&
            gr_get_fmpq(r, t, ctx) == GR_SUCCESS)
        found = 1;
    GR_TMP_CLEAR(t, ctx);
    acb_clear(z);
    return found;
}

/* sin, cos or tan (which = 0, 1, 2) of a real x through u = tan(x/2)
   (tan(x) itself for tan); at rational multiples of pi (denominators up
   to the option GR_TOWER_OPT_TRIG_PI_LIMIT), the algebraic counterpart: polynomials in
   tan(pi/M) */

static int
_real_sin_cos_tan(gr_ptr res, gr_srcptr x, int which, gr_ctx_t ctx)
{
    gr_ptr u, v, w;
    fmpq_t r;
    int status;

    fmpq_init(r);
    if (_pi_multiple(r, x, ctx))
    {
        /* the tangent normal form (a polynomial in tan(pi/M)) */
        status = GR_UNABLE;
        if (fmpz_cmp_ui(fmpq_denref(r), gr_tower_lazy_ctx_get_option(ctx, GR_TOWER_OPT_TRIG_PI_LIMIT)) <= 0)
            status = _gr_tower_lazy_trig_pi_real_locked(res, r, which, ctx);
        fmpq_clear(r);
        return status;
    }

    /* x = r pi + y: the addition formulas, with the algebraic values at
       r pi (so that no generator has pi in its argument) */
    {
        gr_ptr y;
        GR_TMP_INIT(y, ctx);
        if (_gr_tower_lazy_pi_part_locked(r, y, x, ctx) && !fmpq_is_zero(r))
        {
            gr_ptr sr, cr, sy, cy;
            GR_TMP_INIT4(sr, cr, sy, cy, ctx);
            status = gr_pi(sr, ctx);
            status |= gr_mul_fmpq(sr, sr, r, ctx);
            status |= _gr_tower_lazy_cos_impl(cr, sr, ctx);
            status |= _gr_tower_lazy_sin_impl(sr, sr, ctx);
            if (which == 2)
            {
                /* (sin(r pi) + cos(r pi) tan y) / (cos(r pi) - sin(r pi) tan y) */
                status |= _real_sin_cos_tan(sy, y, 2, ctx);
                status |= gr_mul(cy, cr, sy, ctx);
                status |= gr_add(cy, cy, sr, ctx);
                status |= gr_mul(sy, sr, sy, ctx);
                status |= gr_sub(sy, cr, sy, ctx);
                if (status == GR_SUCCESS)
                    status = gr_div(res, cy, sy, ctx);
            }
            else
            {
                status |= _real_sin_cos_tan(sy, y, 0, ctx);
                status |= _real_sin_cos_tan(cy, y, 1, ctx);
                if (which == 0)
                {
                    /* sin(r pi) cos y + cos(r pi) sin y */
                    status |= gr_mul(cy, sr, cy, ctx);
                    status |= gr_mul(sy, cr, sy, ctx);
                    status |= gr_add(res, cy, sy, ctx);
                }
                else
                {
                    /* cos(r pi) cos y - sin(r pi) sin y */
                    status |= gr_mul(cy, cr, cy, ctx);
                    status |= gr_mul(sy, sr, sy, ctx);
                    status |= gr_sub(res, cy, sy, ctx);
                }
            }
            GR_TMP_CLEAR4(sr, cr, sy, cy, ctx);
            GR_TMP_CLEAR(y, ctx);
            fmpq_clear(r);
            return status;
        }
        GR_TMP_CLEAR(y, ctx);
    }
    fmpq_clear(r);

    if (which == 2)
        return _gr_tower_lazy_trig_locked(res, x, GR_TOWER_TAN, ctx);

    GR_TMP_INIT3(u, v, w, ctx);
    status = gr_div_ui(u, x, 2, ctx);
    status |= _gr_tower_lazy_trig_locked(u, u, GR_TOWER_TAN, ctx);
    status |= gr_sqr(v, u, ctx);
    status |= gr_add_ui(w, v, 1, ctx);           /* 1 + u^2 */
    if (which == 0)
    {
        status |= gr_mul_ui(v, u, 2, ctx);
        status |= gr_div(res, v, w, ctx);         /* 2u / (1 + u^2) */
    }
    else
    {
        status |= gr_neg(v, v, ctx);
        status |= gr_add_ui(v, v, 1, ctx);
        status |= gr_div(res, v, w, ctx);         /* (1 - u^2) / (1 + u^2) */
    }
    GR_TMP_CLEAR3(u, v, w, ctx);
    return status;
}

#define REAL_ANGLE_DENOMINATOR_LIMIT 360

/* atan(y) for a real y: r pi when atan(y) / pi = r is a rational with a
   small denominator (recognized numerically, verified exactly through
   tan(r pi) = y), otherwise the generator atan(y) */
static int
_real_atan(gr_ptr res, gr_srcptr y, gr_ctx_t ctx)
{
    acb_t z;
    arb_t t, pi;
    slong q;
    int status = GR_SUCCESS, found = 0;

    acb_init(z);
    arb_init(t);
    arb_init(pi);

    if (gr_tower_lazy_get_acb(z, y, 128, ctx) == GR_SUCCESS)
    {
        arb_atan(t, acb_realref(z), 128);
        arb_const_pi(pi, 128);
        arb_div(t, t, pi, 128);

        for (q = 1; q <= REAL_ANGLE_DENOMINATOR_LIMIT && !found; q++)
        {
            arb_t tq;
            fmpz_t p;
            arb_init(tq);
            fmpz_init(p);
            arb_mul_si(tq, t, q, 128);
            arf_get_fmpz(p, arb_midref(tq), ARF_RND_NEAR);
            arb_sub_fmpz(tq, tq, p, 128);
            if (arb_contains_zero(tq) && mag_cmp_2exp_si(arb_radref(tq), -60) < 0)
            {
                /* candidate r = p / q: tan(r pi) == y? */
                gr_ptr c;
                fmpq_t r;
                GR_TMP_INIT(c, ctx);
                fmpq_init(r);
                fmpz_set(fmpq_numref(r), p);
                fmpz_set_si(fmpq_denref(r), q);
                fmpq_canonicalise(r);
                if (fmpz_is_one(fmpq_denref(r)) || fmpz_cmp_ui(fmpq_denref(r), 2) != 0)
                {
                    if (gr_pi(c, ctx) == GR_SUCCESS && gr_mul_fmpq(c, c, r, ctx) == GR_SUCCESS &&
                        _gr_tower_lazy_tan_impl(c, c, ctx) == GR_SUCCESS && gr_equal(c, y, ctx) == T_TRUE)
                    {
                        status = gr_pi(res, ctx);
                        status |= gr_mul_fmpq(res, res, r, ctx);
                        found = 1;
                    }
                }
                fmpq_clear(r);
                GR_TMP_CLEAR(c, ctx);
                if (!found)
                    q = REAL_ANGLE_DENOMINATOR_LIMIT;   /* (the candidate failed: no other) */
            }
            arb_clear(tq);
            fmpz_clear(p);
        }
    }

    acb_clear(z);
    arb_clear(t);
    arb_clear(pi);

    if (found)
        return status;
    return _gr_tower_lazy_trig_locked(res, y, GR_TOWER_ATAN, ctx);
}

/* asin(y) for real y: atan(y / sqrt(1 - y^2)), +/- pi/2 (GR_DOMAIN
   outside [-1, 1], where the value is not real: checked here, since the
   real view's check applies only to its outermost operation) */
static int
_real_asin(gr_ptr res, gr_srcptr y, gr_ctx_t ctx)
{
    gr_ptr t;
    int status;
    truth_t one, mone;

    status = _real_in_interval(y, -1, 1, 0, ctx);
    if (status != GR_SUCCESS)
        return status;

    GR_TMP_INIT(t, ctx);
    one = gr_is_one(y, ctx);
    mone = gr_is_neg_one(y, ctx);
    if (one == T_TRUE || mone == T_TRUE)
    {
        status = gr_pi(res, ctx);
        status |= gr_div_si(res, res, (one == T_TRUE) ? 2 : -2, ctx);
    }
    else if (one == T_UNKNOWN || mone == T_UNKNOWN)
        status = GR_UNABLE;
    else
    {
        status = gr_sqr(t, y, ctx);
        status |= gr_neg(t, t, ctx);
        status |= gr_add_ui(t, t, 1, ctx);
        status |= gr_sqrt(t, t, ctx);
        status |= gr_div(t, y, t, ctx);
        if (status == GR_SUCCESS)
            status = _real_atan(res, t, ctx);
    }
    GR_TMP_CLEAR(t, ctx);
    return status;
}

static int _gr_tower_lazy_sin_impl(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t x, gr_ctx_t ctx);
static int _gr_tower_lazy_cos_impl(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t x, gr_ctx_t ctx);

/* the real forms where they apply, else the implementations through
   exp and log */
static int
_gr_tower_lazy_sin_real(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t x, gr_ctx_t ctx)
{
    if (_trig_real_arg(x, ctx))
    {
        int status = _real_sin_cos_tan(res, x, 0, ctx);
        if (status == GR_SUCCESS || status == GR_DOMAIN)
            return status;      /* (GR_DOMAIN: a pole of tan) */
    }
    return _gr_tower_lazy_sin_impl(res, x, ctx);
}

static int
_gr_tower_lazy_cos_real(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t x, gr_ctx_t ctx)
{
    if (_trig_real_arg(x, ctx))
    {
        int status = _real_sin_cos_tan(res, x, 1, ctx);
        if (status == GR_SUCCESS || status == GR_DOMAIN)
            return status;      /* (GR_DOMAIN: a pole of tan) */
    }
    return _gr_tower_lazy_cos_impl(res, x, ctx);
}

static int
_gr_tower_lazy_tan_real(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t x, gr_ctx_t ctx)
{
    if (_trig_real_arg(x, ctx))
    {
        int status = _real_sin_cos_tan(res, x, 2, ctx);
        if (status == GR_SUCCESS || status == GR_DOMAIN)
            return status;      /* (GR_DOMAIN: a pole of tan) */
    }
    return _gr_tower_lazy_tan_impl(res, x, ctx);
}

static int
_gr_tower_lazy_atan_real(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t x, gr_ctx_t ctx)
{
    if (_real_arg(x, ctx))
        return _real_atan(res, x, ctx);
    return _gr_tower_lazy_atan_impl(res, x, ctx);
}

#define _gr_tower_lazy_sinh_real _gr_tower_lazy_sinh_impl
#define _gr_tower_lazy_cosh_real _gr_tower_lazy_cosh_impl
#define _gr_tower_lazy_tanh_real _gr_tower_lazy_tanh_impl
#define _gr_tower_lazy_asinh_real _gr_tower_lazy_asinh_impl

#define ALG(ctx) ((FLAGS(ctx) & GR_TOWER_LAZY_ALGEBRAIC) && _gr_tower_lazy_outermost(ctx))
#define REAL(ctx) ((FLAGS(ctx) & GR_TOWER_LAZY_REAL) && _gr_tower_lazy_outermost(ctx))

#define ENTRY_ZERO(name, one) \
int gr_tower_lazy_##name(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t x, gr_ctx_t ctx) \
{ int status; _gr_tower_lazy_lock(ctx); \
  status = ALG(ctx) ? _alg_special_at_zero(res, x, one, ctx) : _gr_tower_lazy_##name##_real(res, x, ctx); \
  if (status == GR_SUCCESS && REAL(ctx)) status = _gr_tower_lazy_realify_locked(res, ctx); \
  _gr_tower_lazy_unlock(ctx); return status; }

ENTRY_ZERO(sinh, 0)
ENTRY_ZERO(cosh, 1)
ENTRY_ZERO(tanh, 0)
ENTRY_ZERO(sin, 0)
ENTRY_ZERO(cos, 1)
ENTRY_ZERO(tan, 0)
ENTRY_ZERO(atan, 0)
ENTRY_ZERO(asinh, 0)

int gr_tower_lazy_atanh(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t x, gr_ctx_t ctx)
{
    int status = GR_SUCCESS;
    _gr_tower_lazy_lock(ctx);
    if (ALG(ctx))
        status = _alg_special_at_zero(res, x, 0, ctx);
    else
    {
        if (REAL(ctx))
            status = _real_in_interval(x, -1, 1, 1, ctx);
        if (status == GR_SUCCESS)
            status = _gr_tower_lazy_atanh_impl(res, x, ctx);
    }
    _gr_tower_lazy_unlock(ctx);
    return status;
}

int gr_tower_lazy_asin(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t x, gr_ctx_t ctx)
{
    int status = GR_SUCCESS;
    _gr_tower_lazy_lock(ctx);
    if (ALG(ctx))
        status = _alg_special_at_zero(res, x, 0, ctx);
    else
    {
        if (REAL(ctx))
            status = _real_in_interval(x, -1, 1, 0, ctx);
        if (status == GR_SUCCESS)
            status = (REAL(ctx) || _real_arg(x, ctx)) && (FLAGS(ctx) & GR_TOWER_LAZY_REAL) ? _real_asin(res, x, ctx) : _gr_tower_lazy_asin_impl(res, x, ctx);
        if (status == GR_SUCCESS && REAL(ctx))
            status = _gr_tower_lazy_realify_locked(res, ctx);
    }
    _gr_tower_lazy_unlock(ctx);
    return status;
}

int gr_tower_lazy_acos(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t x, gr_ctx_t ctx)
{
    int status = GR_SUCCESS;
    _gr_tower_lazy_lock(ctx);
    if (ALG(ctx))
        status = _alg_special_at_one(res, x, ctx);
    else
    {
        if (REAL(ctx))
            status = _real_in_interval(x, -1, 1, 0, ctx);
        if (status == GR_SUCCESS && (FLAGS(ctx) & GR_TOWER_LAZY_REAL) && _real_arg(x, ctx))
        {
            /* pi/2 - asin(x) */
            gr_ptr t;
            GR_TMP_INIT(t, ctx);
            status = _real_asin(t, x, ctx);
            status |= gr_pi(res, ctx);
            status |= gr_div_ui(res, res, 2, ctx);
            status |= gr_sub(res, res, t, ctx);
            GR_TMP_CLEAR(t, ctx);
        }
        else if (status == GR_SUCCESS)
            status = _gr_tower_lazy_acos_impl(res, x, ctx);
        if (status == GR_SUCCESS && REAL(ctx))
            status = _gr_tower_lazy_realify_locked(res, ctx);
    }
    _gr_tower_lazy_unlock(ctx);
    return status;
}

int gr_tower_lazy_acosh(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t x, gr_ctx_t ctx)
{
    int status = GR_SUCCESS;
    _gr_tower_lazy_lock(ctx);
    if (ALG(ctx))
        status = _alg_special_at_one(res, x, ctx);
    else
    {
        if (REAL(ctx))
        {
            gr_ptr t;
            int c;
            GR_TMP_INIT(t, ctx);
            status = gr_one(t, ctx);
            status |= _gr_tower_lazy_cmp_impl(&c, x, t, ctx);
            if (status == GR_SUCCESS && c < 0)
                status = GR_DOMAIN;
            GR_TMP_CLEAR(t, ctx);
        }
        if (status == GR_SUCCESS)
            status = _gr_tower_lazy_acosh_impl(res, x, ctx);
    }
    _gr_tower_lazy_unlock(ctx);
    return status;
}

/* arg(x) is algebraic only when it is 0, for x > 0 */
int gr_tower_lazy_arg(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t x, gr_ctx_t ctx)
{
    int status;
    _gr_tower_lazy_lock(ctx);
    if (ALG(ctx))
    {
        /* arg(x) = 0 for x >= 0 (as in the complex field, arg(0) = 0);
           otherwise a nonzero angle, which is not algebraic */
        int sgn;
        status = _gr_tower_lazy_real_sign_impl(&sgn, x, ctx);
        if (status == GR_SUCCESS)
            status = (sgn >= 0) ? gr_zero(res, ctx) : GR_DOMAIN;
        else if (status == GR_DOMAIN)
            status = GR_DOMAIN;   /* (not real: a nonzero angle, and pi is not present) */
    }
    else
        status = _gr_tower_lazy_arg_impl(res, x, ctx);
    _gr_tower_lazy_unlock(ctx);
    return status;
}

int gr_tower_lazy_re(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t x, gr_ctx_t ctx) { LOCKED(_gr_tower_lazy_re_impl(res, x, ctx)) }
int gr_tower_lazy_im(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t x, gr_ctx_t ctx) { LOCKED(_gr_tower_lazy_im_impl(res, x, ctx)) }
int gr_tower_lazy_sgn(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t x, gr_ctx_t ctx) { LOCKED(_gr_tower_lazy_sgn_impl(res, x, ctx)) }
int gr_tower_lazy_csgn(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t x, gr_ctx_t ctx) { LOCKED(_gr_tower_lazy_csgn_impl(res, x, ctx)) }
int gr_tower_lazy_cmpabs(int * res, const gr_tower_lazy_elem_t x, const gr_tower_lazy_elem_t y, gr_ctx_t ctx) { LOCKED(_gr_tower_lazy_cmpabs_impl(res, x, y, ctx)) }
int gr_tower_lazy_floor(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t x, gr_ctx_t ctx) { LOCKED(_gr_tower_lazy_floor_impl(res, x, ctx)) }
int gr_tower_lazy_ceil(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t x, gr_ctx_t ctx) { LOCKED(_gr_tower_lazy_ceil_impl(res, x, ctx)) }
int gr_tower_lazy_trunc(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t x, gr_ctx_t ctx) { LOCKED(_gr_tower_lazy_trunc_impl(res, x, ctx)) }
int gr_tower_lazy_nint(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t x, gr_ctx_t ctx) { LOCKED(_gr_tower_lazy_nint_impl(res, x, ctx)) }
truth_t gr_tower_lazy_is_real(const gr_tower_lazy_elem_t x, gr_ctx_t ctx) { truth_t status; _gr_tower_lazy_lock(ctx); status = (FLAGS(ctx) & GR_TOWER_LAZY_REAL) && _gr_tower_lazy_outermost(ctx) ? T_TRUE : _gr_tower_lazy_is_real_impl(x, ctx); _gr_tower_lazy_unlock(ctx); return status; }
/* the actual test, for the checks of a real context */
truth_t _gr_tower_lazy_is_real_exact(const gr_tower_lazy_elem_t x, gr_ctx_t ctx) { truth_t status; _gr_tower_lazy_lock(ctx); status = _gr_tower_lazy_is_real_impl(x, ctx); _gr_tower_lazy_unlock(ctx); return status; }
int _gr_tower_lazy_real_sign_locked(int * sign, const gr_tower_lazy_elem_t x, gr_ctx_t ctx) { LOCKED(_gr_tower_lazy_real_sign_impl(sign, x, ctx)) }
int gr_tower_lazy_cmp(int * res, const gr_tower_lazy_elem_t x, const gr_tower_lazy_elem_t y, gr_ctx_t ctx) { LOCKED(_gr_tower_lazy_cmp_impl(res, x, y, ctx)) }
int gr_tower_lazy_get_d(double * res, const gr_tower_lazy_elem_t x, gr_ctx_t ctx) { LOCKED(_gr_tower_lazy_get_d_impl(res, x, ctx)) }

/* -------------------------------------------------------------------- */
/* helpers of the special functions (lazy_special.c, lazy_hypgeom.c,    */
/* lazy_modular.c, lazy_dirichlet.c)                                     */
/* -------------------------------------------------------------------- */

/* the sign of the real number x: numerically if the enclosure decides,
   else exactly */
int
_gr_tower_lazy_real_sign_fast(int * sgn, gr_srcptr x, gr_ctx_t ctx)
{
    acb_t z;
    int status;

    acb_init(z);
    status = gr_tower_lazy_get_acb(z, x, LAZY_CHECK_PREC, ctx);
    if (status == GR_SUCCESS && arb_is_positive(acb_realref(z)))
        *sgn = 1;
    else if (status == GR_SUCCESS && arb_is_negative(acb_realref(z)))
        *sgn = -1;
    else
        status = _gr_tower_lazy_real_sign_locked(sgn, x, ctx);
    acb_clear(z);
    return status;
}

/* the sign of Re(x) - c (c = NULL: 0), likewise */
int
_gr_tower_lazy_re_cmp(int * sgn, gr_srcptr x, const fmpq_t c, gr_ctx_t ctx)
{
    acb_t z;
    arb_t d;
    int status;

    acb_init(z);
    arb_init(d);
    status = gr_tower_lazy_get_acb(z, x, LAZY_CHECK_PREC, ctx);
    if (status == GR_SUCCESS)
    {
        if (c != NULL)
            arb_set_fmpq(d, c, LAZY_CHECK_PREC);
        arb_sub(d, acb_realref(z), d, LAZY_CHECK_PREC);
    }
    if (status == GR_SUCCESS && arb_is_positive(d))
        *sgn = 1;
    else if (status == GR_SUCCESS && arb_is_negative(d))
        *sgn = -1;
    else
    {
        gr_ptr t;
        GR_TMP_INIT(t, ctx);
        status = gr_re(t, x, ctx);
        if (c != NULL)
            status |= gr_sub_fmpq(t, t, c, ctx);
        if (status == GR_SUCCESS)
            status = _gr_tower_lazy_real_sign_locked(sgn, t, ctx);
        GR_TMP_CLEAR(t, ctx);
    }
    acb_clear(z);
    arb_clear(d);
    return status;
}

/* the sign of Im(x) */
int
_gr_tower_lazy_im_sign(int * sgn, gr_srcptr x, gr_ctx_t ctx)
{
    gr_ptr t;
    int status;
    GR_TMP_INIT(t, ctx);
    status = gr_im(t, x, ctx);
    if (status == GR_SUCCESS)
        status = _gr_tower_lazy_real_sign_fast(sgn, t, ctx);
    GR_TMP_CLEAR(t, ctx);
    return status;
}

/* The view restrictions on a result: in a real view, the result must be
   real; in an algebraic view, algebraic (special function values at
   algebraic points are usually transcendental, but this is not known
   in general: GR_UNABLE). */
int
_gr_tower_lazy_view_finish(int status, gr_ptr res, int real, int alg, gr_ctx_t ctx)
{
    if (status != GR_SUCCESS)
        return status;

    if (alg && _gr_tower_lazy_is_algebraic_repr_locked(res, ctx) != T_TRUE)
        return GR_UNABLE;

    if (real)
    {
        truth_t t = _gr_tower_lazy_is_real_exact(res, ctx);
        if (t == T_FALSE)
            return GR_DOMAIN;
        if (t == T_UNKNOWN)
            return GR_UNABLE;
        status = _gr_tower_lazy_realify_locked(res, ctx);
    }

    return status;
}

POP_OPTIONS
