/*
    Copyright (C) 2023 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpq.h"
#include "gr_vec.h"
#include "gr_poly.h"

#ifdef __GNUC__
# define strcmp __builtin_strcmp
#else
# include <string.h>
#endif

/* arb wrappers todo: elementary functions; pfq, coulomb, zeta/dirichlet deflated */

static const char * default_var = "x";

int
gr_poly_add_series(gr_poly_t res, const gr_poly_t poly1,
              const gr_poly_t poly2, slong n, gr_ctx_t ctx)
{
    int status = GR_SUCCESS;
    slong len1, len2, max = FLINT_MAX(poly1->length, poly2->length);

    if (n < 0)
       n = 0;

    max = FLINT_MIN(max, n);
    len1 = FLINT_MIN(poly1->length, max);
    len2 = FLINT_MIN(poly2->length, max);

    gr_poly_fit_length(res, max, ctx);
    status |= _gr_poly_add(res->coeffs, poly1->coeffs, len1, poly2->coeffs, len2, ctx);
    _gr_poly_set_length(res, max, ctx);
    _gr_poly_normalise(res, ctx);
    return status;
}

int
gr_poly_sub_series(gr_poly_t res, const gr_poly_t poly1,
              const gr_poly_t poly2, slong n, gr_ctx_t ctx)
{
    int status = GR_SUCCESS;
    slong len1, len2, max = FLINT_MAX(poly1->length, poly2->length);

    if (n < 0)
       n = 0;

    max = FLINT_MIN(max, n);
    len1 = FLINT_MIN(poly1->length, max);
    len2 = FLINT_MIN(poly2->length, max);

    gr_poly_fit_length(res, max, ctx);
    status |= _gr_poly_sub(res->coeffs, poly1->coeffs, len1, poly2->coeffs, len2, ctx);
    _gr_poly_set_length(res, max, ctx);
    _gr_poly_normalise(res, ctx);
    return status;
}

#define SERIES_ERR_EXACT WORD_MAX
#define SERIES_ERR_MAX WORD_MAX / 4

typedef struct
{
    gr_poly_struct poly;
    slong error;
}
gr_series_struct;

typedef gr_series_struct gr_series_t[1];

typedef struct
{
    slong mod;      /* exact truncation */
    slong prec;     /* default approximate truncation */
}
gr_series_ctx_struct;

typedef gr_series_ctx_struct gr_series_ctx_t[1];


void
gr_series_init(gr_series_t res, gr_series_ctx_t sctx, gr_ctx_t cctx)
{
    gr_poly_init(&res->poly, cctx);
    res->error = SERIES_ERR_EXACT;
}

void
gr_series_clear(gr_series_t res, gr_series_ctx_t sctx, gr_ctx_t cctx)
{
    gr_poly_clear(&res->poly, cctx);
}

int
gr_series_zero(gr_series_t res, gr_series_ctx_t sctx, gr_ctx_t cctx)
{
    res->error = SERIES_ERR_EXACT;
    return gr_poly_zero(&res->poly, cctx);
}

void
gr_series_swap(gr_series_t x, gr_series_t y, gr_series_ctx_t sctx, gr_ctx_t cctx)
{
    gr_series_t tmp;
    *tmp = *x;
    *x = *y;
    *y = *tmp;
}

int
gr_series_randtest(gr_series_t res, flint_rand_t state, slong len, gr_series_ctx_t sctx, gr_ctx_t cctx)
{
    int status = gr_poly_randtest(&res->poly, state, len, cctx);

    len = FLINT_MAX(len, 0);
    len = FLINT_MIN(len, SERIES_ERR_MAX);

    if (n_randint(state, 2))
        res->error = SERIES_ERR_EXACT;
    else
        res->error = n_randint(state, len + 1);

    return status;
}

int
gr_series_write(gr_stream_t out, const gr_series_t x, gr_series_ctx_t sctx, gr_ctx_t cctx)
{
    gr_poly_write(out, &x->poly, "x", cctx);

    if (x->error != SERIES_ERR_EXACT)
    {
        gr_stream_write(out, " + O(x^");
        gr_stream_write_si(out, x->error);
        gr_stream_write(out, ")");
    }

    if (sctx->mod != SERIES_ERR_EXACT)
    {
        gr_stream_write(out, " (mod x^");
        gr_stream_write_si(out, sctx->mod);
        gr_stream_write(out, ")");
    }

    return GR_SUCCESS;
}

int
gr_series_one(gr_series_t res, gr_series_ctx_t sctx, gr_ctx_t cctx)
{
    if (sctx->mod == 0)
    {
        res->error = SERIES_ERR_EXACT;
        return gr_poly_zero(&res->poly, cctx);
    }
    else if (sctx->prec == 0)
    {
        res->error = 0;
        return gr_poly_zero(&res->poly, cctx);
    }
    else
    {
        res->error = SERIES_ERR_EXACT;
        return gr_poly_one(&res->poly, cctx);
    }
}

/* todo: truncate without set */
int
gr_series_set(gr_series_t res, const gr_series_t x, gr_series_ctx_t sctx, gr_ctx_t cctx)
{
    int status = GR_SUCCESS;
    slong len, trunc;

    res->error = x->error;
    status |= gr_poly_set(&res->poly, &x->poly, cctx);

    trunc = FLINT_MIN(sctx->prec, sctx->mod);
    trunc = FLINT_MIN(trunc, res->error);
    len = res->poly.length;

    if (len > trunc)
    {
        if (len <= sctx->mod)
            res->error = SERIES_ERR_EXACT;
        if (len > sctx->prec)
            res->error = FLINT_MIN(res->error, sctx->prec);

        status |= gr_poly_truncate(&res->poly, &res->poly, trunc, cctx);
    }

    return status;
}

int
gr_series_gen(gr_series_t res, gr_series_ctx_t sctx, gr_ctx_t cctx)
{
    int status = GR_SUCCESS;

    status |= gr_poly_zero(&res->poly, cctx);
    status |= gr_poly_set_coeff_ui(&res->poly, 1, 1, cctx);
    res->error = SERIES_ERR_EXACT;

    /* force truncation if needed */
    status |= gr_series_set(res, res, sctx, cctx);
    return status;
}

/* todo: truncate before neg */
int
gr_series_neg(gr_series_t res, const gr_series_t x, gr_series_ctx_t sctx, gr_ctx_t cctx)
{
    int status = GR_SUCCESS;
    slong len, trunc;

    res->error = x->error;
    status |= gr_poly_neg(&res->poly, &x->poly, cctx);

    trunc = FLINT_MIN(sctx->prec, sctx->mod);
    trunc = FLINT_MIN(trunc, res->error);
    len = res->poly.length;

    if (len > trunc)
    {
        if (len <= sctx->mod)
            res->error = SERIES_ERR_EXACT;
        if (len > sctx->prec)
            res->error = FLINT_MIN(res->error, sctx->prec);

        status |= gr_poly_truncate(&res->poly, &res->poly, trunc, cctx);
    }

    return status;
}


int
gr_series_set_gr_poly(gr_series_t res, const gr_poly_t x, gr_series_ctx_t sctx, gr_ctx_t cctx)
{
    gr_series_t tmp;
    tmp->poly = *x;
    tmp->error = SERIES_ERR_EXACT;
    return gr_series_set(res, tmp, sctx, cctx);
}

int
gr_series_set_scalar(gr_series_t res, gr_srcptr x, gr_series_ctx_t sctx, gr_ctx_t cctx)
{
    if (gr_is_zero(x, cctx) == T_TRUE)
    {
        return gr_series_zero(res, sctx, cctx);
    }
    else
    {
        gr_series_t tmp;
        tmp->poly.coeffs = (gr_ptr) x;
        tmp->poly.length = 1;
        tmp->poly.alloc = 1;
        tmp->error = SERIES_ERR_EXACT;
        return gr_series_set(res, tmp, sctx, cctx);
    }
}

int
gr_series_set_si(gr_series_t res, slong c, gr_series_ctx_t sctx, gr_ctx_t cctx)
{
    if (c == 0)
    {
        return gr_series_zero(res, sctx, cctx);
    }
    else
    {
        gr_ptr t;
        int status = GR_SUCCESS;
        GR_TMP_INIT(t, cctx);
        status |= gr_set_si(t, c, cctx);
        status |= gr_series_set_scalar(res, t, sctx, cctx);
        GR_TMP_CLEAR(t, cctx);
        return status;
    }
}

int
gr_series_set_ui(gr_series_t res, ulong c, gr_series_ctx_t sctx, gr_ctx_t cctx)
{
    if (c == 0)
    {
        return gr_series_zero(res, sctx, cctx);
    }
    else
    {
        gr_ptr t;
        int status = GR_SUCCESS;
        GR_TMP_INIT(t, cctx);
        status |= gr_set_ui(t, c, cctx);
        status |= gr_series_set_scalar(res, t, sctx, cctx);
        GR_TMP_CLEAR(t, cctx);
        return status;
    }
}

int
gr_series_set_fmpz(gr_series_t res, const fmpz_t c, gr_series_ctx_t sctx, gr_ctx_t cctx)
{
    if (fmpz_is_zero(c))
    {
        return gr_series_zero(res, sctx, cctx);
    }
    else
    {
        gr_ptr t;
        int status = GR_SUCCESS;
        GR_TMP_INIT(t, cctx);
        status |= gr_set_fmpz(t, c, cctx);
        status |= gr_series_set_scalar(res, t, sctx, cctx);
        GR_TMP_CLEAR(t, cctx);
        return status;
    }
}

/* todo: division by zero valid mod x^0? */
int
gr_series_set_fmpq(gr_series_t res, const fmpq_t c, gr_series_ctx_t sctx, gr_ctx_t cctx)
{
    if (fmpq_is_zero(c))
    {
        return gr_series_zero(res, sctx, cctx);
    }
    else
    {
        gr_ptr t;
        int status = GR_SUCCESS;
        GR_TMP_INIT(t, cctx);
        status |= gr_set_fmpq(t, c, cctx);
        status |= gr_series_set_scalar(res, t, sctx, cctx);
        GR_TMP_CLEAR(t, cctx);
        return status;
    }
}


truth_t
gr_series_is_zero(const gr_series_t x, gr_series_ctx_t sctx, gr_ctx_t cctx)
{
    truth_t is_zero;

    is_zero = gr_poly_is_zero(&x->poly, cctx);

    if (is_zero == T_FALSE)
        return T_FALSE;

    if (x->error == SERIES_ERR_EXACT && is_zero == T_TRUE)
        return T_TRUE;

    return T_UNKNOWN;
}

/* todo: recursive version */
int
gr_series_make_exact(gr_series_t x, gr_series_ctx_t sctx, gr_ctx_t cctx)
{
    x->error = SERIES_ERR_EXACT;
    return GR_SUCCESS;
}



truth_t
_gr_poly_equal2(gr_srcptr poly1, slong len1, gr_srcptr poly2, slong len2, gr_ctx_t ctx)
{
    truth_t eq, eq2;
    slong sz = ctx->sizeof_elem;

    eq = _gr_vec_equal(poly1, poly2, len2, ctx);

    if (len1 == len2 || eq == T_FALSE)
        return eq;

    eq2 = _gr_vec_is_zero(GR_ENTRY(poly1, len2, sz), len1 - len2, ctx);

    return truth_and(eq, eq2);
}

truth_t
gr_series_is_one(const gr_series_t x, gr_series_ctx_t sctx, gr_ctx_t cctx)
{
    truth_t is_zero, is_one;
    slong xlen = x->poly.length;

    if (x->error <= 0)
        return T_UNKNOWN;

    if (sctx->mod == 0)
        return T_TRUE;

    if (xlen == 0)
    {
        /* deal with characteristic 0 */
        if (gr_ctx_is_finite_characteristic(cctx) != T_FALSE)
        {
            gr_ptr tmp;
            GR_TMP_INIT(tmp, cctx);
            is_one = gr_is_one(tmp, cctx);
            GR_TMP_CLEAR(tmp, cctx);
            return is_one;
        }

        return T_FALSE;
    }

    is_one = gr_is_one(x->poly.coeffs, cctx);

    if (is_one == T_FALSE)
        return T_FALSE;

    if (xlen >= 2)
    {
        is_zero = _gr_vec_is_zero(GR_ENTRY(x->poly.coeffs, 1, cctx->sizeof_elem), FLINT_MIN(xlen - 1, x->error - 1), cctx);
        if (is_zero == T_FALSE)
            return T_FALSE;
    }
    else
    {
        is_zero = T_TRUE;
    }

    if (x->error == SERIES_ERR_EXACT && is_zero == T_TRUE && is_one == T_TRUE)
        return T_TRUE;

    return T_UNKNOWN;
}


truth_t
gr_series_equal(const gr_series_t x, const gr_series_t y, gr_series_ctx_t sctx, gr_ctx_t cctx)
{
    truth_t equal;
    slong len, xlen, ylen, xerr, yerr, err;

    xlen = x->poly.length;
    ylen = y->poly.length;
    xerr = x->error;
    yerr = y->error;
    err = FLINT_MIN(xerr, yerr);

    len = FLINT_MAX(xlen, ylen);
    len = FLINT_MIN(len, sctx->mod);
    len = FLINT_MIN(len, err);

    /* compare exactly or within the error bound */
    if (xlen >= ylen)
        equal = _gr_poly_equal2(x->poly.coeffs, FLINT_MIN(xlen, len), y->poly.coeffs, FLINT_MIN(ylen, len), cctx);
    else
        equal = _gr_poly_equal2(y->poly.coeffs, FLINT_MIN(ylen, len), x->poly.coeffs, FLINT_MIN(xlen, len), cctx);

    if (equal == T_FALSE)
        return T_FALSE;

    /* terms >= x^mod are zero by definition */
    if (err >= sctx->mod && equal == T_TRUE)
        return T_TRUE;

    return T_UNKNOWN;
}

int
gr_series_add(gr_series_t res, const gr_series_t x, const gr_series_t y, gr_series_ctx_t sctx, gr_ctx_t cctx)
{
    slong len, xlen, ylen, xerr, yerr, err;
    int status = GR_SUCCESS;

    xlen = x->poly.length;
    ylen = y->poly.length;
    xerr = x->error;
    yerr = y->error;
    err = FLINT_MIN(xerr, yerr);

    /* length of the polynomial sum, ignoring errors */
    len = FLINT_MAX(xlen, ylen);

    /* the result will be truncated */
    if (len > sctx->prec)
        err = FLINT_MIN(err, sctx->prec);

    len = FLINT_MIN(len, sctx->mod);
    len = FLINT_MIN(len, sctx->prec);
    len = FLINT_MIN(len, err);

    /* terms >= x^mod are zero by definition */
    if (err >= sctx->mod)
        err = SERIES_ERR_EXACT;

    res->error = err;
    status |= gr_poly_add_series(&res->poly, &x->poly, &y->poly, len, cctx);
    return status;
}

int
gr_series_sub(gr_series_t res, const gr_series_t x, const gr_series_t y, gr_series_ctx_t sctx, gr_ctx_t cctx)
{
    slong len, xlen, ylen, xerr, yerr, err;
    int status = GR_SUCCESS;

    xlen = x->poly.length;
    ylen = y->poly.length;
    xerr = x->error;
    yerr = y->error;
    err = FLINT_MIN(xerr, yerr);

    /* length of the polynomial sum, ignoring errors */
    len = FLINT_MAX(xlen, ylen);

    /* the result will be truncated */
    if (len > sctx->prec)
        err = FLINT_MIN(err, sctx->prec);

    len = FLINT_MIN(len, sctx->mod);
    len = FLINT_MIN(len, sctx->prec);
    len = FLINT_MIN(len, err);

    /* terms >= x^mod are zero by definition */
    if (err >= sctx->mod)
        err = SERIES_ERR_EXACT;

    res->error = err;
    status |= gr_poly_sub_series(&res->poly, &x->poly, &y->poly, len, cctx);
    return status;
}

int
gr_series_mul(gr_series_t res, const gr_series_t x, const gr_series_t y, gr_series_ctx_t sctx, gr_ctx_t cctx)
{
    slong len, xlen, ylen, xerr, yerr, err;
    int status = GR_SUCCESS;

    xlen = x->poly.length;
    ylen = y->poly.length;
    xerr = x->error;
    yerr = y->error;
    err = FLINT_MIN(xerr, yerr);

    if (xlen == 0 && xerr == SERIES_ERR_EXACT)
        return gr_series_zero(res, sctx, cctx);

    if (ylen == 0 && yerr == SERIES_ERR_EXACT)
        return gr_series_zero(res, sctx, cctx);

    /* length of the polynomial product, ignoring errors */
    if (xlen == 0 || ylen == 0)
        len = 0;
    else
        len = xlen + ylen - 1;

    /* the result will be truncated */
    if (len > sctx->prec)
        err = FLINT_MIN(err, sctx->prec);

    len = FLINT_MIN(len, sctx->mod);
    len = FLINT_MIN(len, sctx->prec);
    len = FLINT_MIN(len, err);

    /* terms >= x^mod are zero by definition */
    if (err >= sctx->mod)
        err = SERIES_ERR_EXACT;

    res->error = err;
    status |= gr_poly_mullow(&res->poly, &x->poly, &y->poly, len, cctx);
    return status;
}

int
gr_series_inv(gr_series_t res, const gr_series_t x, gr_series_ctx_t sctx, gr_ctx_t cctx)
{
    slong len, xlen, xerr, err;
    int status = GR_SUCCESS;

    xlen = x->poly.length;
    xerr = x->error;
    err = xerr;

    if (xlen == 0 && sctx->mod == 0)
        return gr_series_zero(res, sctx, cctx);

    if (xlen == 0 && xerr == SERIES_ERR_EXACT)
    {
        truth_t zero = gr_ctx_is_zero_ring(cctx);

        if (zero == T_TRUE)
            return gr_series_zero(res, sctx, cctx);
        if (zero == T_UNKNOWN)
            return GR_UNABLE;

        return GR_DOMAIN;
    }

    if (xlen == 0 || xerr == 0)
        return GR_UNABLE;

    len = FLINT_MIN(sctx->mod, sctx->prec);
    len = FLINT_MIN(len, err);
    err = len;

    /* terms >= x^mod are zero by definition */
    if (err >= sctx->mod)
        err = SERIES_ERR_EXACT;

    res->error = err;
    status |= gr_poly_inv_series(&res->poly, &x->poly, len, cctx);
    return status;
}

/* todo: handle valuations */
int
gr_series_div(gr_series_t res, const gr_series_t x, const gr_series_t y, gr_series_ctx_t sctx, gr_ctx_t cctx)
{
    slong len, ylen, xerr, yerr, err;
    int status = GR_SUCCESS;

    ylen = y->poly.length;
    xerr = x->error;
    yerr = y->error;

    if (ylen == 0 && yerr == SERIES_ERR_EXACT)
        return GR_DOMAIN;

    if (ylen == 0 || yerr == 0)
        return GR_UNABLE;

    err = FLINT_MIN(xerr, yerr);
    err = FLINT_MIN(err, sctx->prec);

    len = FLINT_MIN(sctx->mod, sctx->prec);
    len = FLINT_MIN(len, err);

    /* terms >= x^mod are zero by definition */
    if (err >= sctx->mod)
        err = SERIES_ERR_EXACT;

    res->error = err;
    status |= gr_poly_div_series(&res->poly, &x->poly, &y->poly, len, cctx);
    return status;
}

int
gr_series_sqrt(gr_series_t res, const gr_series_t x, gr_series_ctx_t sctx, gr_ctx_t cctx)
{
    slong len, xlen, xerr, err;
    int status = GR_SUCCESS;

    xlen = x->poly.length;
    xerr = x->error;
    err = xerr;

    if (xlen == 0 && xerr == SERIES_ERR_EXACT)
        return gr_series_zero(res, sctx, cctx);

    if (xlen == 0 || xerr == 0)
        return GR_UNABLE;

    if (xlen > 1 && gr_is_zero(x->poly.coeffs, cctx) != T_FALSE)
        return GR_UNABLE;

    len = FLINT_MIN(sctx->mod, sctx->prec);
    len = FLINT_MIN(len, err);
    err = len;

    /* terms >= x^mod are zero by definition */
    if (err >= sctx->mod)
        err = SERIES_ERR_EXACT;

    res->error = err;
    status |= gr_poly_sqrt_series(&res->poly, &x->poly, len, cctx);
    return status;
}


#define UNARY_POLY_WRAPPER(func) \
int \
gr_series_ ## func(gr_series_t res, const gr_series_t x, gr_series_ctx_t sctx, gr_ctx_t cctx) \
{ \
    slong len, xlen, xerr, err; \
    int status = GR_SUCCESS; \
 \
    xlen = x->poly.length; \
    xerr = x->error; \
    err = xerr; \
 \
    len = FLINT_MIN(sctx->mod, sctx->prec); \
    len = FLINT_MIN(len, err); \
    err = len; \
 \
    if (xlen <= 1 && xerr == SERIES_ERR_EXACT) \
    { \
        len = FLINT_MIN(len, 1); \
        err = SERIES_ERR_EXACT; \
    } \
 \
    if (err >= sctx->mod) \
        err = SERIES_ERR_EXACT; \
 \
    res->error = err; \
    status |= gr_poly_ ## func ## _series(&res->poly, &x->poly, len, cctx); \
    return status; \
} \

UNARY_POLY_WRAPPER(exp)
UNARY_POLY_WRAPPER(log)
UNARY_POLY_WRAPPER(rsqrt)
UNARY_POLY_WRAPPER(tan)

UNARY_POLY_WRAPPER(asin)
UNARY_POLY_WRAPPER(acos)
UNARY_POLY_WRAPPER(atan)
UNARY_POLY_WRAPPER(asinh)
UNARY_POLY_WRAPPER(acosh)
UNARY_POLY_WRAPPER(atanh)


#include "arb_poly.h"
#include "acb_poly.h"
#include "arb_hypgeom.h"
#include "acb_hypgeom.h"
#include "acb_dirichlet.h"
#include "acb_elliptic.h"
#include "acb_modular.h"

#define ARB_WRAPPER(func, arb_func, acb_func) \
int \
gr_series_ ## func(gr_series_t res, const gr_series_t x, gr_series_ctx_t sctx, gr_ctx_t cctx) \
{ \
    slong xlen, len, xerr, err; \
    slong prec; \
    int status = GR_SUCCESS; \
 \
    if (cctx->which_ring != GR_CTX_RR_ARB && cctx->which_ring != GR_CTX_CC_ACB) \
        return GR_UNABLE; \
 \
    xlen = x->poly.length; \
    xerr = x->error; \
    err = xerr; \
 \
    len = FLINT_MIN(sctx->mod, sctx->prec); \
    len = FLINT_MIN(len, err); \
    err = len; \
 \
    if (err >= sctx->mod) \
        err = SERIES_ERR_EXACT; \
  \
    if (xlen <= 1 && xerr == SERIES_ERR_EXACT) \
    { \
        len = FLINT_MIN(len, 1); \
        err = SERIES_ERR_EXACT; \
    } \
 \
    res->error = err; \
 \
    prec = _gr_ctx_get_real_prec(cctx); \
 \
    if (cctx->which_ring == GR_CTX_RR_ARB) \
    { \
        arb_func((arb_poly_struct *) &res->poly, (arb_poly_struct *) &x->poly, len, prec); \
        if (!_arb_vec_is_finite((arb_ptr) res->poly.coeffs, res->poly.length)) \
            status = GR_UNABLE; \
    } \
    else \
    { \
        acb_func((acb_poly_struct *) &res->poly, (acb_poly_struct *) &x->poly, len, prec); \
        if (!_arb_vec_is_finite((arb_ptr) res->poly.coeffs, 2 * res->poly.length)) \
            status = GR_UNABLE; \
    } \
 \
    return status; \
} \

ARB_WRAPPER(gamma, arb_poly_gamma_series, acb_poly_gamma_series)
ARB_WRAPPER(rgamma, arb_poly_rgamma_series, acb_poly_rgamma_series)
ARB_WRAPPER(lgamma, arb_poly_lgamma_series, acb_poly_lgamma_series)
ARB_WRAPPER(digamma, arb_poly_digamma_series, acb_poly_digamma_series)
ARB_WRAPPER(erf, arb_hypgeom_erf_series, acb_hypgeom_erf_series)
ARB_WRAPPER(erfc, arb_hypgeom_erfc_series, acb_hypgeom_erfc_series)
ARB_WRAPPER(erfi, arb_hypgeom_erfi_series, acb_hypgeom_erfi_series)
ARB_WRAPPER(exp_integral_ei, arb_hypgeom_ei_series, acb_hypgeom_ei_series)
ARB_WRAPPER(cos_integral, arb_hypgeom_ci_series, acb_hypgeom_ci_series)
ARB_WRAPPER(cosh_integral, arb_hypgeom_chi_series, acb_hypgeom_chi_series)
ARB_WRAPPER(sin_integral, arb_hypgeom_si_series, acb_hypgeom_si_series)
ARB_WRAPPER(sinh_integral, arb_hypgeom_shi_series, acb_hypgeom_shi_series)

int
gr_series_fresnel(gr_series_t res1, gr_series_t res2, const gr_series_t x, int normalized, gr_series_ctx_t sctx, gr_ctx_t cctx)
{
    slong xlen, len, xerr, err;
    slong prec;
    int status = GR_SUCCESS;

    if (cctx->which_ring != GR_CTX_RR_ARB && cctx->which_ring != GR_CTX_CC_ACB)
        return GR_UNABLE;

    xlen = x->poly.length;
    xerr = x->error;
    err = xerr;

    len = FLINT_MIN(sctx->mod, sctx->prec);
    len = FLINT_MIN(len, err);
    err = len;

    if (err >= sctx->mod)
        err = SERIES_ERR_EXACT;

    if (xlen <= 1 && xerr == SERIES_ERR_EXACT)
    {
        len = FLINT_MIN(len, 1);
        err = SERIES_ERR_EXACT;
    }

    if (res1 != NULL) res1->error = err;
    if (res2 != NULL) res2->error = err;

    prec = _gr_ctx_get_real_prec(cctx);

    if (cctx->which_ring == GR_CTX_RR_ARB)
    {
        arb_hypgeom_fresnel_series(res1 ? (arb_poly_struct *) &res1->poly : NULL,
                                   res2 ? (arb_poly_struct *) &res2->poly : NULL,
                                    (arb_poly_struct *) &x->poly, normalized, len, prec);
        if (res1 && !_arb_vec_is_finite((arb_ptr) res1->poly.coeffs, res1->poly.length))
            status = GR_UNABLE;
        if (res2 && !_arb_vec_is_finite((arb_ptr) res2->poly.coeffs, res2->poly.length))
            status = GR_UNABLE;
    }
    else
    {
        acb_hypgeom_fresnel_series(res1 ? (acb_poly_struct *) &res1->poly : NULL,
                                   res2 ? (acb_poly_struct *) &res2->poly : NULL,
                                    (acb_poly_struct *) &x->poly, normalized, len, prec);
        if (res1 && !_arb_vec_is_finite((arb_ptr) res1->poly.coeffs, 2 * res1->poly.length))
            status = GR_UNABLE;
        if (res2 && !_arb_vec_is_finite((arb_ptr) res2->poly.coeffs, 2 * res2->poly.length))
            status = GR_UNABLE;
    }

    return status;
}

int
gr_series_fresnel_s(gr_series_t res, const gr_series_t x, int normalized, gr_series_ctx_t sctx, gr_ctx_t cctx)
{
    return gr_series_fresnel(res, NULL, x, normalized, sctx, cctx);
}

int
gr_series_fresnel_c(gr_series_t res, const gr_series_t x, int normalized, gr_series_ctx_t sctx, gr_ctx_t cctx)
{
    return gr_series_fresnel(NULL, res, x, normalized, sctx, cctx);
}

int
gr_series_airy(gr_series_t res1, gr_series_t res2, gr_series_t res3, gr_series_t res4, const gr_series_t x, gr_series_ctx_t sctx, gr_ctx_t cctx)
{
    slong xlen, len, xerr, err;
    slong prec;
    int status = GR_SUCCESS;

    if (cctx->which_ring != GR_CTX_RR_ARB && cctx->which_ring != GR_CTX_CC_ACB)
        return GR_UNABLE;

    xlen = x->poly.length;
    xerr = x->error;
    err = xerr;

    len = FLINT_MIN(sctx->mod, sctx->prec);
    len = FLINT_MIN(len, err);
    err = len;

    if (err >= sctx->mod)
        err = SERIES_ERR_EXACT;

    if (xlen <= 1 && xerr == SERIES_ERR_EXACT)
    {
        len = FLINT_MIN(len, 1);
        err = SERIES_ERR_EXACT;
    }

    if (res1 != NULL) res1->error = err;
    if (res2 != NULL) res2->error = err;
    if (res3 != NULL) res3->error = err;
    if (res4 != NULL) res4->error = err;

    prec = _gr_ctx_get_real_prec(cctx);

    if (cctx->which_ring == GR_CTX_RR_ARB)
    {
        arb_hypgeom_airy_series(res1 ? (arb_poly_struct *) &res1->poly : NULL,
                                res2 ? (arb_poly_struct *) &res2->poly : NULL,
                                res3 ? (arb_poly_struct *) &res3->poly : NULL,
                                res4 ? (arb_poly_struct *) &res4->poly : NULL,
                                (arb_poly_struct *) &x->poly, len, prec);
        if (res1 && !_arb_vec_is_finite((arb_ptr) res1->poly.coeffs, res1->poly.length))
            status = GR_UNABLE;
        if (res2 && !_arb_vec_is_finite((arb_ptr) res2->poly.coeffs, res2->poly.length))
            status = GR_UNABLE;
        if (res3 && !_arb_vec_is_finite((arb_ptr) res3->poly.coeffs, res3->poly.length))
            status = GR_UNABLE;
        if (res4 && !_arb_vec_is_finite((arb_ptr) res4->poly.coeffs, res4->poly.length))
            status = GR_UNABLE;
    }
    else
    {
        acb_hypgeom_airy_series(res1 ? (acb_poly_struct *) &res1->poly : NULL,
                                res2 ? (acb_poly_struct *) &res2->poly : NULL,
                                res3 ? (acb_poly_struct *) &res3->poly : NULL,
                                res4 ? (acb_poly_struct *) &res4->poly : NULL,
                                (acb_poly_struct *) &x->poly, len, prec);
        if (res1 && !_arb_vec_is_finite((arb_ptr) res1->poly.coeffs, 2 * res1->poly.length))
            status = GR_UNABLE;
        if (res2 && !_arb_vec_is_finite((arb_ptr) res2->poly.coeffs, 2 * res2->poly.length))
            status = GR_UNABLE;
        if (res3 && !_arb_vec_is_finite((arb_ptr) res3->poly.coeffs, 2 * res3->poly.length))
            status = GR_UNABLE;
        if (res4 && !_arb_vec_is_finite((arb_ptr) res4->poly.coeffs, 2 * res4->poly.length))
            status = GR_UNABLE;
    }

    return status;
}

int
gr_series_airy_ai(gr_series_t res, const gr_series_t x, gr_series_ctx_t sctx, gr_ctx_t cctx)
{
    return gr_series_airy(res, NULL, NULL, NULL, x, sctx, cctx);
}

int
gr_series_airy_ai_prime(gr_series_t res, const gr_series_t x, gr_series_ctx_t sctx, gr_ctx_t cctx)
{
    return gr_series_airy(NULL, res, NULL, NULL, x, sctx, cctx);
}

int
gr_series_airy_bi(gr_series_t res, const gr_series_t x, gr_series_ctx_t sctx, gr_ctx_t cctx)
{
    return gr_series_airy(NULL, NULL, res, NULL, x, sctx, cctx);
}

int
gr_series_airy_bi_prime(gr_series_t res, const gr_series_t x, gr_series_ctx_t sctx, gr_ctx_t cctx)
{
    return gr_series_airy(NULL, NULL, NULL, res, x, sctx, cctx);
}


int
gr_series_log_integral(gr_series_t res, const gr_series_t x, int offset, gr_series_ctx_t sctx, gr_ctx_t cctx)
{
    slong xlen, len, xerr, err;
    slong prec;
    int status = GR_SUCCESS;

    if (cctx->which_ring != GR_CTX_RR_ARB && cctx->which_ring != GR_CTX_CC_ACB)
        return GR_UNABLE;

    xlen = x->poly.length;
    xerr = x->error;
    err = xerr;

    len = FLINT_MIN(sctx->mod, sctx->prec);
    len = FLINT_MIN(len, err);
    err = len;

    if (err >= sctx->mod)
        err = SERIES_ERR_EXACT;

    if (xlen <= 1 && xerr == SERIES_ERR_EXACT)
    {
        len = FLINT_MIN(len, 1);
        err = SERIES_ERR_EXACT;
    }

    res->error = err;

    prec = _gr_ctx_get_real_prec(cctx);

    if (cctx->which_ring == GR_CTX_RR_ARB)
    {
        arb_hypgeom_li_series((arb_poly_struct *) &res->poly, (arb_poly_struct *) &x->poly, offset, len, prec);
        if (!_arb_vec_is_finite((arb_ptr) res->poly.coeffs, res->poly.length))
            status = GR_UNABLE; \
    }
    else
    {
        acb_hypgeom_li_series((acb_poly_struct *) &res->poly, (acb_poly_struct *) &x->poly, offset, len, prec);
        if (!_arb_vec_is_finite((arb_ptr) res->poly.coeffs, 2 * res->poly.length))
            status = GR_UNABLE;
    }

    return status;
}

int
gr_series_gamma_upper(gr_series_t res, const gr_series_t s, const gr_series_t x, int regularized, gr_series_ctx_t sctx, gr_ctx_t cctx)
{
    slong xlen, len, xerr, err;
    slong prec;
    int status = GR_SUCCESS;

    if (cctx->which_ring != GR_CTX_RR_ARB && cctx->which_ring != GR_CTX_CC_ACB)
        return GR_UNABLE;

    xlen = x->poly.length;
    xerr = x->error;
    err = xerr;

    len = FLINT_MIN(sctx->mod, sctx->prec);
    len = FLINT_MIN(len, err);
    err = len;

    if (err >= sctx->mod)
        err = SERIES_ERR_EXACT;

    if (xlen <= 1 && xerr == SERIES_ERR_EXACT)
    {
        len = FLINT_MIN(len, 1);
        err = SERIES_ERR_EXACT;
    }

    /* we only handle constant s */
    if (len >= 2 && s->poly.length >= 2)
        return GR_UNABLE;

    res->error = err;

    prec = _gr_ctx_get_real_prec(cctx);

    if (cctx->which_ring == GR_CTX_RR_ARB)
    {
        arb_t t;
        arb_init(t);
        arb_poly_get_coeff_arb(t, (arb_poly_struct *) &s->poly, 0);
        arb_hypgeom_gamma_upper_series((arb_poly_struct *) &res->poly, t, (arb_poly_struct *) &x->poly, regularized, len, prec);
        if (!_arb_vec_is_finite((arb_ptr) res->poly.coeffs, res->poly.length))
            status = GR_UNABLE;
        arb_clear(t);
    }
    else
    {
        acb_t t;
        acb_init(t);
        acb_poly_get_coeff_acb(t, (acb_poly_struct *) &s->poly, 0);
        acb_hypgeom_gamma_upper_series((acb_poly_struct *) &res->poly, t, (acb_poly_struct *) &x->poly, regularized, len, prec);
        if (!_arb_vec_is_finite((arb_ptr) res->poly.coeffs, 2 * res->poly.length))
            status = GR_UNABLE;
        acb_clear(t);
    }

    return status;
}

int
gr_series_gamma_lower(gr_series_t res, const gr_series_t s, const gr_series_t x, int regularized, gr_series_ctx_t sctx, gr_ctx_t cctx)
{
    slong xlen, len, xerr, err;
    slong prec;
    int status = GR_SUCCESS;

    if (cctx->which_ring != GR_CTX_RR_ARB && cctx->which_ring != GR_CTX_CC_ACB)
        return GR_UNABLE;

    xlen = x->poly.length;
    xerr = x->error;
    err = xerr;

    len = FLINT_MIN(sctx->mod, sctx->prec);
    len = FLINT_MIN(len, err);
    err = len;

    if (err >= sctx->mod)
        err = SERIES_ERR_EXACT;

    if (xlen <= 1 && xerr == SERIES_ERR_EXACT)
    {
        len = FLINT_MIN(len, 1);
        err = SERIES_ERR_EXACT;
    }

    /* we only handle constant s */
    if (len >= 2 && s->poly.length >= 2)
        return GR_UNABLE;

    res->error = err;

    prec = _gr_ctx_get_real_prec(cctx);

    if (cctx->which_ring == GR_CTX_RR_ARB)
    {
        arb_t t;
        arb_init(t);
        arb_poly_get_coeff_arb(t, (arb_poly_struct *) &s->poly, 0);
        arb_hypgeom_gamma_lower_series((arb_poly_struct *) &res->poly, t, (arb_poly_struct *) &x->poly, regularized, len, prec);
        if (!_arb_vec_is_finite((arb_ptr) res->poly.coeffs, res->poly.length))
            status = GR_UNABLE;
        arb_clear(t);
    }
    else
    {
        acb_t t;
        acb_init(t);
        acb_poly_get_coeff_acb(t, (acb_poly_struct *) &s->poly, 0);
        acb_hypgeom_gamma_lower_series((acb_poly_struct *) &res->poly, t, (acb_poly_struct *) &x->poly, regularized, len, prec);
        if (!_arb_vec_is_finite((arb_ptr) res->poly.coeffs, 2 * res->poly.length))
            status = GR_UNABLE;
        acb_clear(t);
    }

    return status;
}

int
gr_series_beta_lower(gr_series_t res, const gr_series_t a, const gr_series_t b, const gr_series_t x, int regularized, gr_series_ctx_t sctx, gr_ctx_t cctx)
{
    slong xlen, len, xerr, err;
    slong prec;
    int status = GR_SUCCESS;

    if (cctx->which_ring != GR_CTX_RR_ARB && cctx->which_ring != GR_CTX_CC_ACB)
        return GR_UNABLE;

    xlen = x->poly.length;
    xerr = x->error;
    err = xerr;

    len = FLINT_MIN(sctx->mod, sctx->prec);
    len = FLINT_MIN(len, err);
    err = len;

    if (err >= sctx->mod)
        err = SERIES_ERR_EXACT;

    if (xlen <= 1 && xerr == SERIES_ERR_EXACT)
    {
        len = FLINT_MIN(len, 1);
        err = SERIES_ERR_EXACT;
    }

    /* we only handle constant a, b */
    if (len >= 2 && a->poly.length >= 2)
        return GR_UNABLE;
    if (len >= 2 && b->poly.length >= 2)
        return GR_UNABLE;

    res->error = err;

    prec = _gr_ctx_get_real_prec(cctx);

    if (cctx->which_ring == GR_CTX_RR_ARB)
    {
        arb_t t, u;
        arb_init(t);
        arb_init(u);
        arb_poly_get_coeff_arb(t, (arb_poly_struct *) &a->poly, 0);
        arb_poly_get_coeff_arb(u, (arb_poly_struct *) &b->poly, 0);
        arb_hypgeom_beta_lower_series((arb_poly_struct *) &res->poly, t, u, (arb_poly_struct *) &x->poly, regularized, len, prec);
        if (!_arb_vec_is_finite((arb_ptr) res->poly.coeffs, res->poly.length))
            status = GR_UNABLE;
        arb_clear(t);
        arb_clear(u);
    }
    else
    {
        acb_t t, u;
        acb_init(t);
        acb_init(u);
        acb_poly_get_coeff_acb(t, (acb_poly_struct *) &a->poly, 0);
        acb_poly_get_coeff_acb(u, (acb_poly_struct *) &b->poly, 0);
        acb_hypgeom_beta_lower_series((acb_poly_struct *) &res->poly, t, u, (acb_poly_struct *) &x->poly, regularized, len, prec);
        if (!_arb_vec_is_finite((arb_ptr) res->poly.coeffs, 2 * res->poly.length))
            status = GR_UNABLE;
        acb_clear(t);
        acb_clear(u);
    }

    return status;
}

int
gr_series_polylog(gr_series_t res, const gr_series_t s, const gr_series_t z, gr_series_ctx_t sctx, gr_ctx_t cctx)
{
    slong xlen, len, xerr, err;
    slong prec;
    int status = GR_SUCCESS;

    if (cctx->which_ring != GR_CTX_CC_ACB)
        return GR_UNABLE;

    xlen = s->poly.length;
    xerr = s->error;
    err = xerr;

    len = FLINT_MIN(sctx->mod, sctx->prec);
    len = FLINT_MIN(len, err);
    err = len;

    if (err >= sctx->mod)
        err = SERIES_ERR_EXACT;

    if (xlen <= 1 && xerr == SERIES_ERR_EXACT)
    {
        len = FLINT_MIN(len, 1);
        err = SERIES_ERR_EXACT;
    }

    /* we only handle constant z */
    if (len >= 2 && z->poly.length >= 2)
        return GR_UNABLE;

    res->error = err;

    prec = _gr_ctx_get_real_prec(cctx);

    {
        acb_t t;
        acb_init(t);
        acb_poly_get_coeff_acb(t, (acb_poly_struct *) &z->poly, 0);
        acb_poly_polylog_series((acb_poly_struct *) &res->poly, (acb_poly_struct *) &s->poly, t, len, prec);
        if (!_arb_vec_is_finite((arb_ptr) res->poly.coeffs, 2 * res->poly.length))
            status = GR_UNABLE;
        acb_clear(t);
    }

    return status;
}

int
gr_series_hurwitz_zeta(gr_series_t res, const gr_series_t s, const gr_series_t z, gr_series_ctx_t sctx, gr_ctx_t cctx)
{
    slong xlen, len, xerr, err;
    slong prec;
    int status = GR_SUCCESS;

    if (cctx->which_ring != GR_CTX_CC_ACB)
        return GR_UNABLE;

    xlen = s->poly.length;
    xerr = s->error;
    err = xerr;

    len = FLINT_MIN(sctx->mod, sctx->prec);
    len = FLINT_MIN(len, err);
    err = len;

    if (err >= sctx->mod)
        err = SERIES_ERR_EXACT;

    if (xlen <= 1 && xerr == SERIES_ERR_EXACT)
    {
        len = FLINT_MIN(len, 1);
        err = SERIES_ERR_EXACT;
    }

    /* we only handle constant z */
    if (len >= 2 && z->poly.length >= 2)
        return GR_UNABLE;

    res->error = err;

    prec = _gr_ctx_get_real_prec(cctx);

    {
        acb_t t;
        acb_init(t);
        acb_poly_get_coeff_acb(t, (acb_poly_struct *) &z->poly, 0);
        acb_poly_zeta_series((acb_poly_struct *) &res->poly, (acb_poly_struct *) &s->poly, t, 0, len, prec);
        if (!_arb_vec_is_finite((arb_ptr) res->poly.coeffs, 2 * res->poly.length))
            status = GR_UNABLE;
        acb_clear(t);
    }

    return status;
}

int
gr_series_dirichlet_l(gr_series_t res, const dirichlet_group_t G, const dirichlet_char_t chi, const gr_series_t x, gr_series_ctx_t sctx, gr_ctx_t cctx)
{
    slong xlen, len, xerr, err;
    slong prec;
    int status = GR_SUCCESS;

    if (cctx->which_ring != GR_CTX_CC_ACB)
        return GR_UNABLE;

    xlen = x->poly.length;
    xerr = x->error;
    err = xerr;

    len = FLINT_MIN(sctx->mod, sctx->prec);
    len = FLINT_MIN(len, err);
    err = len;

    if (err >= sctx->mod)
        err = SERIES_ERR_EXACT;

    if (xlen <= 1 && xerr == SERIES_ERR_EXACT)
    {
        len = FLINT_MIN(len, 1);
        err = SERIES_ERR_EXACT;
    }

    res->error = err;

    prec = _gr_ctx_get_real_prec(cctx);

    {
        acb_dirichlet_l_series((acb_poly_struct *) &res->poly, (acb_poly_struct *) &x->poly, G, chi, 0, len, prec);
        if (!_arb_vec_is_finite((arb_ptr) res->poly.coeffs, 2 * res->poly.length))
            status = GR_UNABLE;
    }

    return status;
}

int
gr_series_dirichlet_hardy_theta(gr_series_t res, const dirichlet_group_t G, const dirichlet_char_t chi, const gr_series_t x, gr_series_ctx_t sctx, gr_ctx_t cctx)
{
    slong xlen, len, xerr, err;
    slong prec;
    int status = GR_SUCCESS;

    if (cctx->which_ring != GR_CTX_CC_ACB)
        return GR_UNABLE;

    xlen = x->poly.length;
    xerr = x->error;
    err = xerr;

    len = FLINT_MIN(sctx->mod, sctx->prec);
    len = FLINT_MIN(len, err);
    err = len;

    if (err >= sctx->mod)
        err = SERIES_ERR_EXACT;

    if (xlen <= 1 && xerr == SERIES_ERR_EXACT)
    {
        len = FLINT_MIN(len, 1);
        err = SERIES_ERR_EXACT;
    }

    res->error = err;

    prec = _gr_ctx_get_real_prec(cctx);

    {
        acb_dirichlet_hardy_theta_series((acb_poly_struct *) &res->poly, (acb_poly_struct *) &x->poly, G, chi, len, prec);
        if (!_arb_vec_is_finite((arb_ptr) res->poly.coeffs, 2 * res->poly.length))
            status = GR_UNABLE;
    }

    return status;
}

int
gr_series_dirichlet_hardy_z(gr_series_t res, const dirichlet_group_t G, const dirichlet_char_t chi, const gr_series_t x, gr_series_ctx_t sctx, gr_ctx_t cctx)
{
    slong xlen, len, xerr, err;
    slong prec;
    int status = GR_SUCCESS;

    if (cctx->which_ring != GR_CTX_CC_ACB)
        return GR_UNABLE;

    xlen = x->poly.length;
    xerr = x->error;
    err = xerr;

    len = FLINT_MIN(sctx->mod, sctx->prec);
    len = FLINT_MIN(len, err);
    err = len;

    if (err >= sctx->mod)
        err = SERIES_ERR_EXACT;

    if (xlen <= 1 && xerr == SERIES_ERR_EXACT)
    {
        len = FLINT_MIN(len, 1);
        err = SERIES_ERR_EXACT;
    }

    res->error = err;

    prec = _gr_ctx_get_real_prec(cctx);

    {
        acb_dirichlet_hardy_z_series((acb_poly_struct *) &res->poly, (acb_poly_struct *) &x->poly, G, chi, len, prec);
        if (!_arb_vec_is_finite((arb_ptr) res->poly.coeffs, 2 * res->poly.length))
            status = GR_UNABLE;
    }

    return status;
}

int
gr_series_jacobi_theta(gr_series_t res1, gr_series_t res2, gr_series_t res3, gr_series_t res4, const gr_series_t x, const gr_series_t tau, gr_series_ctx_t sctx, gr_ctx_t cctx)
{
    slong xlen, len, xerr, err;
    slong prec;
    int status = GR_SUCCESS;

    if (cctx->which_ring != GR_CTX_CC_ACB)
        return GR_UNABLE;

    xlen = x->poly.length;
    xerr = x->error;
    err = xerr;

    len = FLINT_MIN(sctx->mod, sctx->prec);
    len = FLINT_MIN(len, err);
    err = len;

    if (err >= sctx->mod)
        err = SERIES_ERR_EXACT;

    if (xlen <= 1 && xerr == SERIES_ERR_EXACT)
    {
        len = FLINT_MIN(len, 1);
        err = SERIES_ERR_EXACT;
    }

    /* we only handle constant z */
    if (len >= 2 && tau->poly.length >= 2)
        return GR_UNABLE;

    if (res1 != NULL) res1->error = err;
    if (res2 != NULL) res2->error = err;
    if (res3 != NULL) res3->error = err;
    if (res4 != NULL) res4->error = err;

    prec = _gr_ctx_get_real_prec(cctx);

    {
        acb_t t;
        acb_init(t);
        acb_poly_get_coeff_acb(t, (acb_poly_struct *) &tau->poly, 0);
        acb_modular_theta_series(res1 ? (acb_poly_struct *) &res1->poly : NULL,
                                res2 ? (acb_poly_struct *) &res2->poly : NULL,
                                res3 ? (acb_poly_struct *) &res3->poly : NULL,
                                res4 ? (acb_poly_struct *) &res4->poly : NULL,
                                (acb_poly_struct *) &x->poly, t, len, prec);
        if (res1 && !_arb_vec_is_finite((arb_ptr) res1->poly.coeffs, 2 * res1->poly.length))
            status = GR_UNABLE;
        if (res2 && !_arb_vec_is_finite((arb_ptr) res2->poly.coeffs, 2 * res2->poly.length))
            status = GR_UNABLE;
        if (res3 && !_arb_vec_is_finite((arb_ptr) res3->poly.coeffs, 2 * res3->poly.length))
            status = GR_UNABLE;
        if (res4 && !_arb_vec_is_finite((arb_ptr) res4->poly.coeffs, 2 * res4->poly.length))
            status = GR_UNABLE;
        acb_clear(t);
    }

    return status;
}

int
gr_series_jacobi_theta_1(gr_series_t res, const gr_series_t x, const gr_series_t tau, gr_series_ctx_t sctx, gr_ctx_t cctx)
{
    return gr_series_jacobi_theta(res, NULL, NULL, NULL, x, tau, sctx, cctx);
}

int
gr_series_jacobi_theta_2(gr_series_t res, const gr_series_t x, const gr_series_t tau, gr_series_ctx_t sctx, gr_ctx_t cctx)
{
    return gr_series_jacobi_theta(NULL, res, NULL, NULL, x, tau, sctx, cctx);
}

int
gr_series_jacobi_theta_3(gr_series_t res, const gr_series_t x, const gr_series_t tau, gr_series_ctx_t sctx, gr_ctx_t cctx)
{
    return gr_series_jacobi_theta(NULL, NULL, res, NULL, x, tau, sctx, cctx);
}

int
gr_series_jacobi_theta_4(gr_series_t res, const gr_series_t x, const gr_series_t tau, gr_series_ctx_t sctx, gr_ctx_t cctx)
{
    return gr_series_jacobi_theta(NULL, NULL, NULL, res, x, tau, sctx, cctx);
}


int
gr_series_agm1(gr_series_t res, const gr_series_t x, gr_series_ctx_t sctx, gr_ctx_t cctx)
{
    slong xlen, len, xerr, err;
    slong prec;
    int status = GR_SUCCESS;

    if (cctx->which_ring != GR_CTX_CC_ACB)
        return GR_UNABLE;

    xlen = x->poly.length;
    xerr = x->error;
    err = xerr;

    len = FLINT_MIN(sctx->mod, sctx->prec);
    len = FLINT_MIN(len, err);
    err = len;

    if (err >= sctx->mod)
        err = SERIES_ERR_EXACT;

    if (xlen <= 1 && xerr == SERIES_ERR_EXACT)
    {
        len = FLINT_MIN(len, 1);
        err = SERIES_ERR_EXACT;
    }

    res->error = err;

    prec = _gr_ctx_get_real_prec(cctx);

    {
        acb_poly_agm1_series((acb_poly_struct *) &res->poly, (acb_poly_struct *) &x->poly, len, prec);
        if (!_arb_vec_is_finite((arb_ptr) res->poly.coeffs, 2 * res->poly.length))
            status = GR_UNABLE;
    }

    return status;
}

int
gr_series_elliptic_k(gr_series_t res, const gr_series_t x, gr_series_ctx_t sctx, gr_ctx_t cctx)
{
    slong xlen, len, xerr, err;
    slong prec;
    int status = GR_SUCCESS;

    if (cctx->which_ring != GR_CTX_CC_ACB)
        return GR_UNABLE;

    xlen = x->poly.length;
    xerr = x->error;
    err = xerr;

    len = FLINT_MIN(sctx->mod, sctx->prec);
    len = FLINT_MIN(len, err);
    err = len;

    if (err >= sctx->mod)
        err = SERIES_ERR_EXACT;

    if (xlen <= 1 && xerr == SERIES_ERR_EXACT)
    {
        len = FLINT_MIN(len, 1);
        err = SERIES_ERR_EXACT;
    }

    res->error = err;

    prec = _gr_ctx_get_real_prec(cctx);

    {
        acb_poly_elliptic_k_series((acb_poly_struct *) &res->poly, (acb_poly_struct *) &x->poly, len, prec);
        if (!_arb_vec_is_finite((arb_ptr) res->poly.coeffs, 2 * res->poly.length))
            status = GR_UNABLE;
    }

    return status;
}

int
gr_series_weierstrass_p(gr_series_t res, const gr_series_t x, const gr_series_t tau, gr_series_ctx_t sctx, gr_ctx_t cctx)
{
    slong xlen, len, xerr, err;
    slong prec;
    int status = GR_SUCCESS;

    if (cctx->which_ring != GR_CTX_CC_ACB)
        return GR_UNABLE;

    xlen = x->poly.length;
    xerr = x->error;
    err = xerr;

    len = FLINT_MIN(sctx->mod, sctx->prec);
    len = FLINT_MIN(len, err);
    err = len;

    if (err >= sctx->mod)
        err = SERIES_ERR_EXACT;

    if (xlen <= 1 && xerr == SERIES_ERR_EXACT)
    {
        len = FLINT_MIN(len, 1);
        err = SERIES_ERR_EXACT;
    }

    /* we only handle constant tau */
    if (len >= 2 && tau->poly.length >= 2)
        return GR_UNABLE;

    res->error = err;

    prec = _gr_ctx_get_real_prec(cctx);

    {
        acb_t t;
        acb_init(t);
        acb_poly_get_coeff_acb(t, (acb_poly_struct *) &tau->poly, 0);
        acb_elliptic_p_series((acb_poly_struct *) &res->poly, (acb_poly_struct *) &x->poly, t, len, prec);
        if (!_arb_vec_is_finite((arb_ptr) res->poly.coeffs, 2 * res->poly.length))
            status = GR_UNABLE;
        acb_clear(t);
    }

    return status;
}

typedef struct
{
    gr_series_struct * entries;
    slong alloc;
    slong length;
}
gr_series_vec_struct;

typedef gr_series_vec_struct gr_series_vec_t[1];

int
gr_series_hypgeom_pfq(gr_series_t res, const gr_series_vec_t a, const gr_series_vec_t b, const gr_series_t x, int regularized, gr_series_ctx_t sctx, gr_ctx_t cctx)
{
    slong len, err, i;
    slong prec;
    slong p, q, pp, qq;
    int have_one;
    int status = GR_SUCCESS;
    acb_poly_struct *aa, *bb;

    if (cctx->which_ring != GR_CTX_CC_ACB)
        return GR_UNABLE;

    p = a->length;
    q = b->length;

    err = x->error;
    for (i = 0; i < p; i++)
        err = FLINT_MIN(err, a->entries[i].error);
    for (i = 0; i < q; i++)
        err = FLINT_MIN(err, b->entries[i].error);

    len = FLINT_MIN(sctx->mod, sctx->prec);
    len = FLINT_MIN(len, err);
    err = len;

    if (err >= sctx->mod)
        err = SERIES_ERR_EXACT;

    aa = GR_TMP_ALLOC(sizeof(acb_poly_struct) * (p + q + 1));

    have_one = 0;

    for (i = 0; i < p; i++)
    {
        if (!have_one && acb_poly_is_one((acb_poly_struct *) (&a->entries[i].poly)))
        {
            have_one = 1;
            continue;
        }

        aa[i - have_one].coeffs = a->entries[i].poly.coeffs;
        aa[i - have_one].alloc = a->entries[i].poly.alloc;
        aa[i - have_one].length = a->entries[i].poly.length;
    }

    if (have_one)
    {
        pp = p - 1;
        qq = q;
        bb = aa + pp;
    }
    else
    {
        pp = p;
        qq = q + 1;
        bb = aa + pp;

        acb_poly_init(bb + qq - 1);
        acb_poly_one(bb + qq - 1);
    }

    for (i = 0; i < q; i++)
    {
        bb[i].coeffs = b->entries[i].poly.coeffs;
        bb[i].alloc = b->entries[i].poly.alloc;
        bb[i].length = b->entries[i].poly.length;
    }

    prec = _gr_ctx_get_real_prec(cctx);

    res->error = err;

    acb_hypgeom_pfq_series_direct((acb_poly_struct *) &res->poly, aa, pp, bb, qq, (const acb_poly_struct *) &x->poly, regularized, -1, len, prec);

    if (!_acb_vec_is_finite((acb_ptr) res->poly.coeffs, res->poly.length))
        status = GR_UNABLE;

    if (!have_one)
        acb_poly_clear(bb + qq - 1);

    GR_TMP_FREE(aa, sizeof(acb_poly_struct) * (p + q + 1));

    return status;
}

typedef struct
{
    gr_ctx_struct * base_ring;
    gr_series_ctx_struct sctx;
    char * var;
}
series_ctx_t;

#define SERIES_CTX(ring_ctx) ((series_ctx_t *)((ring_ctx)))
#define SERIES_ELEM_CTX(ring_ctx) (SERIES_CTX(ring_ctx)->base_ring)
#define SERIES_SCTX(ring_ctx) (&(((series_ctx_t *)((ring_ctx)))->sctx))



static int _gr_gr_series_ctx_write(gr_stream_t out, gr_ctx_t ctx)
{
    gr_stream_write(out, "Power series over ");
    gr_ctx_write(out, SERIES_ELEM_CTX(ctx));
    if (ctx->which_ring == GR_CTX_GR_SERIES_MOD)
    {
        gr_stream_write(out, " mod x^");
        gr_stream_write_si(out, SERIES_SCTX(ctx)->mod);
    }

    gr_stream_write(out, " with precision ");
    gr_stream_write(out, "O(x^");
    gr_stream_write_si(out, SERIES_SCTX(ctx)->prec);
    gr_stream_write(out, ")");

    return GR_SUCCESS;
}

static void _gr_gr_series_init(gr_series_t res, gr_ctx_t ctx) { gr_series_init(res, SERIES_SCTX(ctx), SERIES_ELEM_CTX(ctx)); }
static void _gr_gr_series_clear(gr_series_t res, gr_ctx_t ctx) { gr_series_clear(res, SERIES_SCTX(ctx), SERIES_ELEM_CTX(ctx)); }
static void _gr_gr_series_swap(gr_series_t x, gr_series_t y, gr_ctx_t ctx) { gr_series_swap(x, y, SERIES_SCTX(ctx), SERIES_ELEM_CTX(ctx)); }
static int _gr_gr_series_zero(gr_series_t res, gr_ctx_t ctx) { return gr_series_zero(res, SERIES_SCTX(ctx), SERIES_ELEM_CTX(ctx)); }
static int _gr_gr_series_randtest(gr_series_t res, flint_rand_t state, gr_ctx_t ctx) { return gr_series_randtest(res, state, 6, SERIES_SCTX(ctx), SERIES_ELEM_CTX(ctx)); }
static int _gr_gr_series_write(gr_stream_t out, const gr_series_t x, gr_ctx_t ctx) { return gr_series_write(out, x, SERIES_SCTX(ctx), SERIES_ELEM_CTX(ctx)); }
static int _gr_gr_series_one(gr_series_t res, gr_ctx_t ctx) { return gr_series_one(res, SERIES_SCTX(ctx), SERIES_ELEM_CTX(ctx)); }
static int _gr_gr_series_gen(gr_series_t res, gr_ctx_t ctx) { return gr_series_gen(res, SERIES_SCTX(ctx), SERIES_ELEM_CTX(ctx)); }
static int _gr_gr_series_set(gr_series_t res, const gr_series_t x, gr_ctx_t ctx) { return gr_series_set(res, x, SERIES_SCTX(ctx), SERIES_ELEM_CTX(ctx)); }
static void _gr_gr_series_set_shallow(gr_series_t res, const gr_series_t x, gr_ctx_t ctx) { *res = *x; }
static int _gr_gr_series_set_si(gr_series_t res, slong c, gr_ctx_t ctx) { return gr_series_set_si(res, c, SERIES_SCTX(ctx), SERIES_ELEM_CTX(ctx)); }
static int _gr_gr_series_set_ui(gr_series_t res, ulong c, gr_ctx_t ctx) { return gr_series_set_ui(res, c, SERIES_SCTX(ctx), SERIES_ELEM_CTX(ctx)); }
static int _gr_gr_series_set_fmpz(gr_series_t res, const fmpz_t c, gr_ctx_t ctx) { return gr_series_set_fmpz(res, c, SERIES_SCTX(ctx), SERIES_ELEM_CTX(ctx)); }
static int _gr_gr_series_set_fmpq(gr_series_t res, const fmpq_t c, gr_ctx_t ctx) { return gr_series_set_fmpq(res, c, SERIES_SCTX(ctx), SERIES_ELEM_CTX(ctx)); }
static truth_t _gr_gr_series_is_zero(const gr_series_t x, gr_ctx_t ctx) { return gr_series_is_zero(x, SERIES_SCTX(ctx), SERIES_ELEM_CTX(ctx)); }
static truth_t _gr_gr_series_is_one(const gr_series_t x, gr_ctx_t ctx) { return gr_series_is_one(x, SERIES_SCTX(ctx), SERIES_ELEM_CTX(ctx)); }
static truth_t _gr_gr_series_equal(const gr_series_t x, const gr_series_t y, gr_ctx_t ctx) { return gr_series_equal(x, y, SERIES_SCTX(ctx), SERIES_ELEM_CTX(ctx)); }
static int _gr_gr_series_neg(gr_series_t res, const gr_series_t x, gr_ctx_t ctx) { return gr_series_neg(res, x, SERIES_SCTX(ctx), SERIES_ELEM_CTX(ctx)); }
static int _gr_gr_series_add(gr_series_t res, const gr_series_t x, const gr_series_t y, gr_ctx_t ctx) { return gr_series_add(res, x, y, SERIES_SCTX(ctx), SERIES_ELEM_CTX(ctx)); }
static int _gr_gr_series_sub(gr_series_t res, const gr_series_t x, const gr_series_t y, gr_ctx_t ctx) { return gr_series_sub(res, x, y, SERIES_SCTX(ctx), SERIES_ELEM_CTX(ctx)); }
static int _gr_gr_series_mul(gr_series_t res, const gr_series_t x, const gr_series_t y, gr_ctx_t ctx) { return gr_series_mul(res, x, y, SERIES_SCTX(ctx), SERIES_ELEM_CTX(ctx)); }
static int _gr_gr_series_inv(gr_series_t res, const gr_series_t x, gr_ctx_t ctx) { return gr_series_inv(res, x, SERIES_SCTX(ctx), SERIES_ELEM_CTX(ctx)); }
static int _gr_gr_series_div(gr_series_t res, const gr_series_t x, const gr_series_t y, gr_ctx_t ctx) { return gr_series_div(res, x, y, SERIES_SCTX(ctx), SERIES_ELEM_CTX(ctx)); }
static int _gr_gr_series_sqrt(gr_series_t res, const gr_series_t x, gr_ctx_t ctx) { return gr_series_sqrt(res, x, SERIES_SCTX(ctx), SERIES_ELEM_CTX(ctx)); }
static int _gr_gr_series_rsqrt(gr_series_t res, const gr_series_t x, gr_ctx_t ctx) { return gr_series_rsqrt(res, x, SERIES_SCTX(ctx), SERIES_ELEM_CTX(ctx)); }
static int _gr_gr_series_exp(gr_series_t res, const gr_series_t x, gr_ctx_t ctx) { return gr_series_exp(res, x, SERIES_SCTX(ctx), SERIES_ELEM_CTX(ctx)); }
static int _gr_gr_series_log(gr_series_t res, const gr_series_t x, gr_ctx_t ctx) { return gr_series_log(res, x, SERIES_SCTX(ctx), SERIES_ELEM_CTX(ctx)); }
static int _gr_gr_series_tan(gr_series_t res, const gr_series_t x, gr_ctx_t ctx) { return gr_series_tan(res, x, SERIES_SCTX(ctx), SERIES_ELEM_CTX(ctx)); }
static int _gr_gr_series_asin(gr_series_t res, const gr_series_t x, gr_ctx_t ctx) { return gr_series_asin(res, x, SERIES_SCTX(ctx), SERIES_ELEM_CTX(ctx)); }
static int _gr_gr_series_acos(gr_series_t res, const gr_series_t x, gr_ctx_t ctx) { return gr_series_acos(res, x, SERIES_SCTX(ctx), SERIES_ELEM_CTX(ctx)); }
static int _gr_gr_series_atan(gr_series_t res, const gr_series_t x, gr_ctx_t ctx) { return gr_series_atan(res, x, SERIES_SCTX(ctx), SERIES_ELEM_CTX(ctx)); }
static int _gr_gr_series_asinh(gr_series_t res, const gr_series_t x, gr_ctx_t ctx) { return gr_series_asinh(res, x, SERIES_SCTX(ctx), SERIES_ELEM_CTX(ctx)); }
static int _gr_gr_series_acosh(gr_series_t res, const gr_series_t x, gr_ctx_t ctx) { return gr_series_acosh(res, x, SERIES_SCTX(ctx), SERIES_ELEM_CTX(ctx)); }
static int _gr_gr_series_atanh(gr_series_t res, const gr_series_t x, gr_ctx_t ctx) { return gr_series_atanh(res, x, SERIES_SCTX(ctx), SERIES_ELEM_CTX(ctx)); }

static int _gr_gr_series_gamma(gr_series_t res, const gr_series_t x, gr_ctx_t ctx) { return gr_series_gamma(res, x, SERIES_SCTX(ctx), SERIES_ELEM_CTX(ctx)); }
static int _gr_gr_series_rgamma(gr_series_t res, const gr_series_t x, gr_ctx_t ctx) { return gr_series_rgamma(res, x, SERIES_SCTX(ctx), SERIES_ELEM_CTX(ctx)); }
static int _gr_gr_series_lgamma(gr_series_t res, const gr_series_t x, gr_ctx_t ctx) { return gr_series_lgamma(res, x, SERIES_SCTX(ctx), SERIES_ELEM_CTX(ctx)); }
static int _gr_gr_series_digamma(gr_series_t res, const gr_series_t x, gr_ctx_t ctx) { return gr_series_digamma(res, x, SERIES_SCTX(ctx), SERIES_ELEM_CTX(ctx)); }
static int _gr_gr_series_erf(gr_series_t res, const gr_series_t x, gr_ctx_t ctx) { return gr_series_erf(res, x, SERIES_SCTX(ctx), SERIES_ELEM_CTX(ctx)); }
static int _gr_gr_series_erfc(gr_series_t res, const gr_series_t x, gr_ctx_t ctx) { return gr_series_erfc(res, x, SERIES_SCTX(ctx), SERIES_ELEM_CTX(ctx)); }
static int _gr_gr_series_erfi(gr_series_t res, const gr_series_t x, gr_ctx_t ctx) { return gr_series_erfi(res, x, SERIES_SCTX(ctx), SERIES_ELEM_CTX(ctx)); }

static int _gr_gr_series_fresnel(gr_series_t res1, gr_series_t res2, const gr_series_t x, int normalized, gr_ctx_t ctx) { return gr_series_fresnel(res1, res2, x, normalized, SERIES_SCTX(ctx), SERIES_ELEM_CTX(ctx)); }
static int _gr_gr_series_fresnel_s(gr_series_t res, const gr_series_t x, int normalized, gr_ctx_t ctx) { return gr_series_fresnel_s(res, x, normalized, SERIES_SCTX(ctx), SERIES_ELEM_CTX(ctx)); }
static int _gr_gr_series_fresnel_c(gr_series_t res, const gr_series_t x, int normalized, gr_ctx_t ctx) { return gr_series_fresnel_c(res, x, normalized, SERIES_SCTX(ctx), SERIES_ELEM_CTX(ctx)); }

static int _gr_gr_series_airy(gr_series_t res1, gr_series_t res2, gr_series_t res3, gr_series_t res4, const gr_series_t x, gr_ctx_t ctx) { return gr_series_airy(res1, res2, res3, res4, x, SERIES_SCTX(ctx), SERIES_ELEM_CTX(ctx)); }
static int _gr_gr_series_airy_ai(gr_series_t res, const gr_series_t x, gr_ctx_t ctx) { return gr_series_airy_ai(res, x, SERIES_SCTX(ctx), SERIES_ELEM_CTX(ctx)); }
static int _gr_gr_series_airy_ai_prime(gr_series_t res, const gr_series_t x, gr_ctx_t ctx) { return gr_series_airy_ai_prime(res, x, SERIES_SCTX(ctx), SERIES_ELEM_CTX(ctx)); }
static int _gr_gr_series_airy_bi(gr_series_t res, const gr_series_t x, gr_ctx_t ctx) { return gr_series_airy_bi(res, x, SERIES_SCTX(ctx), SERIES_ELEM_CTX(ctx)); }
static int _gr_gr_series_airy_bi_prime(gr_series_t res, const gr_series_t x, gr_ctx_t ctx) { return gr_series_airy_bi_prime(res, x, SERIES_SCTX(ctx), SERIES_ELEM_CTX(ctx)); }

static int _gr_gr_series_exp_integral_ei(gr_series_t res, const gr_series_t x, gr_ctx_t ctx) { return gr_series_erf(res, x, SERIES_SCTX(ctx), SERIES_ELEM_CTX(ctx)); }
static int _gr_gr_series_cos_integral(gr_series_t res, const gr_series_t x, gr_ctx_t ctx) { return gr_series_cos_integral(res, x, SERIES_SCTX(ctx), SERIES_ELEM_CTX(ctx)); }
static int _gr_gr_series_cosh_integral(gr_series_t res, const gr_series_t x, gr_ctx_t ctx) { return gr_series_cosh_integral(res, x, SERIES_SCTX(ctx), SERIES_ELEM_CTX(ctx)); }
static int _gr_gr_series_sin_integral(gr_series_t res, const gr_series_t x, gr_ctx_t ctx) { return gr_series_sin_integral(res, x, SERIES_SCTX(ctx), SERIES_ELEM_CTX(ctx)); }
static int _gr_gr_series_sinh_integral(gr_series_t res, const gr_series_t x, gr_ctx_t ctx) { return gr_series_sinh_integral(res, x, SERIES_SCTX(ctx), SERIES_ELEM_CTX(ctx)); }
static int _gr_gr_series_log_integral(gr_series_t res, const gr_series_t x, int offset, gr_ctx_t ctx) { return gr_series_log_integral(res, x, offset, SERIES_SCTX(ctx), SERIES_ELEM_CTX(ctx)); }
static int _gr_gr_series_gamma_upper(gr_series_t res, const gr_series_t x, const gr_series_t y, int regularized, gr_ctx_t ctx) { return gr_series_gamma_upper(res, x, y, regularized, SERIES_SCTX(ctx), SERIES_ELEM_CTX(ctx)); }
static int _gr_gr_series_gamma_lower(gr_series_t res, const gr_series_t x, const gr_series_t y, int regularized, gr_ctx_t ctx) { return gr_series_gamma_lower(res, x, y, regularized, SERIES_SCTX(ctx), SERIES_ELEM_CTX(ctx)); }
static int _gr_gr_series_beta_lower(gr_series_t res, const gr_series_t a, const gr_series_t b, const gr_series_t x, int regularized, gr_ctx_t ctx) { return gr_series_beta_lower(res, a, b, x, regularized, SERIES_SCTX(ctx), SERIES_ELEM_CTX(ctx)); }

static int _gr_gr_series_hypgeom_pfq(gr_series_t res, const gr_series_vec_t a, const gr_series_vec_t b, const gr_series_t x, int regularized, gr_ctx_t ctx) { return gr_series_hypgeom_pfq(res, a, b, x, regularized, SERIES_SCTX(ctx), SERIES_ELEM_CTX(ctx)); }

static int _gr_gr_series_polylog(gr_series_t res, const gr_series_t s, const gr_series_t z, gr_ctx_t ctx) { return gr_series_polylog(res, s, z, SERIES_SCTX(ctx), SERIES_ELEM_CTX(ctx)); }
static int _gr_gr_series_hurwitz_zeta(gr_series_t res, const gr_series_t s, const gr_series_t z, gr_ctx_t ctx) { return gr_series_hurwitz_zeta(res, s, z, SERIES_SCTX(ctx), SERIES_ELEM_CTX(ctx)); }

static int _gr_gr_series_dirichlet_l(gr_series_t res, const dirichlet_group_t G, const dirichlet_char_t chi, const gr_series_t s, const gr_series_t z, gr_ctx_t ctx) { return gr_series_dirichlet_l(res, G, chi, s, SERIES_SCTX(ctx), SERIES_ELEM_CTX(ctx)); }
static int _gr_gr_series_dirichlet_hardy_theta(gr_series_t res, const dirichlet_group_t G, const dirichlet_char_t chi, const gr_series_t s, const gr_series_t z, gr_ctx_t ctx) { return gr_series_dirichlet_hardy_theta(res, G, chi, s, SERIES_SCTX(ctx), SERIES_ELEM_CTX(ctx)); }
static int _gr_gr_series_dirichlet_hardy_z(gr_series_t res, const dirichlet_group_t G, const dirichlet_char_t chi, const gr_series_t s, const gr_series_t z, gr_ctx_t ctx) { return gr_series_dirichlet_hardy_z(res, G, chi, s, SERIES_SCTX(ctx), SERIES_ELEM_CTX(ctx)); }

static int _gr_gr_series_jacobi_theta(gr_series_t res1, gr_series_t res2, gr_series_t res3, gr_series_t res4, const gr_series_t x, const gr_series_t tau, gr_ctx_t ctx) { return gr_series_jacobi_theta(res1, res2, res3, res4, x, tau, SERIES_SCTX(ctx), SERIES_ELEM_CTX(ctx)); }
static int _gr_gr_series_jacobi_theta_1(gr_series_t res, const gr_series_t x, const gr_series_t tau, gr_ctx_t ctx) { return gr_series_jacobi_theta_1(res, x, tau, SERIES_SCTX(ctx), SERIES_ELEM_CTX(ctx)); }
static int _gr_gr_series_jacobi_theta_2(gr_series_t res, const gr_series_t x, const gr_series_t tau, gr_ctx_t ctx) { return gr_series_jacobi_theta_2(res, x, tau, SERIES_SCTX(ctx), SERIES_ELEM_CTX(ctx)); }
static int _gr_gr_series_jacobi_theta_3(gr_series_t res, const gr_series_t x, const gr_series_t tau, gr_ctx_t ctx) { return gr_series_jacobi_theta_3(res, x, tau, SERIES_SCTX(ctx), SERIES_ELEM_CTX(ctx)); }
static int _gr_gr_series_jacobi_theta_4(gr_series_t res, const gr_series_t x, const gr_series_t tau, gr_ctx_t ctx) { return gr_series_jacobi_theta_4(res, x, tau, SERIES_SCTX(ctx), SERIES_ELEM_CTX(ctx)); }

static int _gr_gr_series_agm1(gr_series_t res, const gr_series_t x, gr_ctx_t ctx) { return gr_series_agm1(res, x, SERIES_SCTX(ctx), SERIES_ELEM_CTX(ctx)); }
static int _gr_gr_series_elliptic_k(gr_series_t res, const gr_series_t x, gr_ctx_t ctx) { return gr_series_elliptic_k(res, x, SERIES_SCTX(ctx), SERIES_ELEM_CTX(ctx)); }
static int _gr_gr_series_weierstrass_p(gr_series_t res, const gr_series_t x, const gr_series_t tau, gr_ctx_t ctx) { return gr_series_weierstrass_p(res, x, tau, SERIES_SCTX(ctx), SERIES_ELEM_CTX(ctx)); }


static int
_gr_gr_series_set_other(gr_series_t res, gr_srcptr x, gr_ctx_t x_ctx, gr_ctx_t ctx)
{
    if (x_ctx == ctx)
    {
        return _gr_gr_series_set(res, x, ctx);
    }
    else if (x_ctx == SERIES_ELEM_CTX(ctx))
    {
        return gr_series_set_scalar(res, x, SERIES_SCTX(ctx), SERIES_ELEM_CTX(ctx));
    }
    else if ((x_ctx->which_ring == GR_CTX_GR_SERIES || x_ctx->which_ring == GR_CTX_GR_SERIES_MOD)
        && !strcmp(SERIES_CTX(x_ctx)->var, SERIES_CTX(ctx)->var))
    {
        int status = GR_SUCCESS;

        /* todo: check compatibility? */

        status |= gr_poly_set_gr_poly_other(&res->poly, &((gr_series_struct *) x)->poly, SERIES_ELEM_CTX(x_ctx), SERIES_ELEM_CTX(ctx));
        res->error = ((gr_series_struct *) x)->error;
        /* possible truncation */
        status |= gr_series_set(res, res, SERIES_SCTX(ctx), SERIES_ELEM_CTX(ctx));
        return status;
    }
    else if (x_ctx->which_ring == GR_CTX_GR_POLY && !strcmp(POLYNOMIAL_CTX(x_ctx)->var, SERIES_CTX(ctx)->var))
    {
        int status = GR_SUCCESS;

        /* todo: check compatibility? */

        status |= gr_poly_set_gr_poly_other(&res->poly, x, SERIES_ELEM_CTX(x_ctx), SERIES_ELEM_CTX(ctx));
        res->error = SERIES_ERR_EXACT;
        /* possible truncation */
        status |= gr_series_set(res, res, SERIES_SCTX(ctx), SERIES_ELEM_CTX(ctx));
        return status;
    }
    else
    {
        int status = GR_SUCCESS;

        gr_poly_fit_length(&res->poly, 1, SERIES_ELEM_CTX(ctx));
        status = gr_set_other(res->poly.coeffs, x, x_ctx, SERIES_ELEM_CTX(ctx));
        if (status == GR_SUCCESS)
        {
            _gr_poly_set_length(&res->poly, 1, SERIES_ELEM_CTX(ctx));
            _gr_poly_normalise(&res->poly, SERIES_ELEM_CTX(ctx));
        }
        else
            _gr_poly_set_length(&res->poly, 0, SERIES_ELEM_CTX(ctx));

        res->error = SERIES_ERR_EXACT;
        /* possible truncation */
        status |= gr_series_set(res, res, SERIES_SCTX(ctx), SERIES_ELEM_CTX(ctx));

        return status;
    }

    return GR_UNABLE;
}



int _gr_series_methods_initialized = 0;

gr_static_method_table _gr_series_methods;

gr_method_tab_input _gr_series_methods_input[] =
{
    {GR_METHOD_CTX_WRITE,   (gr_funcptr) _gr_gr_series_ctx_write},
    {GR_METHOD_INIT,        (gr_funcptr) _gr_gr_series_init},
    {GR_METHOD_CLEAR,       (gr_funcptr) _gr_gr_series_clear},
    {GR_METHOD_SWAP,        (gr_funcptr) _gr_gr_series_swap},
    {GR_METHOD_SET_SHALLOW, (gr_funcptr) _gr_gr_series_set_shallow},
    {GR_METHOD_RANDTEST,    (gr_funcptr) _gr_gr_series_randtest},
    {GR_METHOD_WRITE,       (gr_funcptr) _gr_gr_series_write},
    {GR_METHOD_ZERO,        (gr_funcptr) _gr_gr_series_zero},
    {GR_METHOD_ONE,         (gr_funcptr) _gr_gr_series_one},
    {GR_METHOD_IS_ZERO,     (gr_funcptr) _gr_gr_series_is_zero},
    {GR_METHOD_IS_ONE,      (gr_funcptr) _gr_gr_series_is_one},
    {GR_METHOD_EQUAL,       (gr_funcptr) _gr_gr_series_equal},
    {GR_METHOD_GEN,         (gr_funcptr) _gr_gr_series_gen},
    {GR_METHOD_SET,         (gr_funcptr) _gr_gr_series_set},
    {GR_METHOD_SET_UI,      (gr_funcptr) _gr_gr_series_set_ui},
    {GR_METHOD_SET_SI,      (gr_funcptr) _gr_gr_series_set_si},
    {GR_METHOD_SET_FMPZ,    (gr_funcptr) _gr_gr_series_set_fmpz},
    {GR_METHOD_SET_FMPQ,    (gr_funcptr) _gr_gr_series_set_fmpq},
    {GR_METHOD_SET_OTHER,   (gr_funcptr) _gr_gr_series_set_other},
    {GR_METHOD_NEG,         (gr_funcptr) _gr_gr_series_neg},
    {GR_METHOD_ADD,         (gr_funcptr) _gr_gr_series_add},
    {GR_METHOD_SUB,         (gr_funcptr) _gr_gr_series_sub},
    {GR_METHOD_MUL,         (gr_funcptr) _gr_gr_series_mul},
    {GR_METHOD_INV,         (gr_funcptr) _gr_gr_series_inv},
    {GR_METHOD_DIV,         (gr_funcptr) _gr_gr_series_div},
    {GR_METHOD_DIV,         (gr_funcptr) _gr_gr_series_div},
    {GR_METHOD_SQRT,        (gr_funcptr) _gr_gr_series_sqrt},
    {GR_METHOD_RSQRT,       (gr_funcptr) _gr_gr_series_rsqrt},
    {GR_METHOD_EXP,         (gr_funcptr) _gr_gr_series_exp},
    {GR_METHOD_LOG,         (gr_funcptr) _gr_gr_series_log},
    {GR_METHOD_TAN,         (gr_funcptr) _gr_gr_series_tan},
    {GR_METHOD_ASIN,        (gr_funcptr) _gr_gr_series_asin},
    {GR_METHOD_ACOS,        (gr_funcptr) _gr_gr_series_acos},
    {GR_METHOD_ATAN,        (gr_funcptr) _gr_gr_series_atan},
    {GR_METHOD_ASINH,       (gr_funcptr) _gr_gr_series_asinh},
    {GR_METHOD_ACOSH,       (gr_funcptr) _gr_gr_series_acosh},
    {GR_METHOD_ATANH,       (gr_funcptr) _gr_gr_series_atanh},
    {GR_METHOD_GAMMA,       (gr_funcptr) _gr_gr_series_gamma},
    {GR_METHOD_RGAMMA,      (gr_funcptr) _gr_gr_series_rgamma},
    {GR_METHOD_LGAMMA,      (gr_funcptr) _gr_gr_series_lgamma},
    {GR_METHOD_DIGAMMA,     (gr_funcptr) _gr_gr_series_digamma},
    {GR_METHOD_ERF,         (gr_funcptr) _gr_gr_series_erf},
    {GR_METHOD_ERFC,        (gr_funcptr) _gr_gr_series_erfc},
    {GR_METHOD_ERFI,        (gr_funcptr) _gr_gr_series_erfi},
    {GR_METHOD_FRESNEL,     (gr_funcptr) _gr_gr_series_fresnel},
    {GR_METHOD_FRESNEL_S,   (gr_funcptr) _gr_gr_series_fresnel_s},
    {GR_METHOD_FRESNEL_C,   (gr_funcptr) _gr_gr_series_fresnel_c},
    {GR_METHOD_AIRY,        (gr_funcptr) _gr_gr_series_airy},
    {GR_METHOD_AIRY_AI,        (gr_funcptr) _gr_gr_series_airy_ai},
    {GR_METHOD_AIRY_AI_PRIME,  (gr_funcptr) _gr_gr_series_airy_ai_prime},
    {GR_METHOD_AIRY_BI,        (gr_funcptr) _gr_gr_series_airy_bi},
    {GR_METHOD_AIRY_BI_PRIME,  (gr_funcptr) _gr_gr_series_airy_bi_prime},
    {GR_METHOD_EXP_INTEGRAL_EI,       (gr_funcptr) _gr_gr_series_exp_integral_ei},
    {GR_METHOD_COS_INTEGRAL,          (gr_funcptr) _gr_gr_series_cos_integral},
    {GR_METHOD_COSH_INTEGRAL,         (gr_funcptr) _gr_gr_series_cosh_integral},
    {GR_METHOD_SIN_INTEGRAL,          (gr_funcptr) _gr_gr_series_sin_integral},
    {GR_METHOD_SINH_INTEGRAL,         (gr_funcptr) _gr_gr_series_sinh_integral},
    {GR_METHOD_LOG_INTEGRAL,          (gr_funcptr) _gr_gr_series_log_integral},
    {GR_METHOD_GAMMA_UPPER,           (gr_funcptr) _gr_gr_series_gamma_upper},
    {GR_METHOD_GAMMA_LOWER,           (gr_funcptr) _gr_gr_series_gamma_lower},
    {GR_METHOD_BETA_LOWER,            (gr_funcptr) _gr_gr_series_beta_lower},
    {GR_METHOD_HYPGEOM_PFQ,           (gr_funcptr) _gr_gr_series_hypgeom_pfq},
    {GR_METHOD_HURWITZ_ZETA,          (gr_funcptr) _gr_gr_series_hurwitz_zeta},
    {GR_METHOD_POLYLOG,               (gr_funcptr) _gr_gr_series_polylog},
    {GR_METHOD_DIRICHLET_L,           (gr_funcptr) _gr_gr_series_dirichlet_l},
    {GR_METHOD_DIRICHLET_HARDY_Z,     (gr_funcptr) _gr_gr_series_dirichlet_hardy_z},
    {GR_METHOD_DIRICHLET_HARDY_THETA, (gr_funcptr) _gr_gr_series_dirichlet_hardy_theta},
    {GR_METHOD_JACOBI_THETA,          (gr_funcptr) _gr_gr_series_jacobi_theta},
    {GR_METHOD_JACOBI_THETA_1,          (gr_funcptr) _gr_gr_series_jacobi_theta_1},
    {GR_METHOD_JACOBI_THETA_2,          (gr_funcptr) _gr_gr_series_jacobi_theta_2},
    {GR_METHOD_JACOBI_THETA_3,          (gr_funcptr) _gr_gr_series_jacobi_theta_3},
    {GR_METHOD_JACOBI_THETA_4,          (gr_funcptr) _gr_gr_series_jacobi_theta_4},
    {GR_METHOD_AGM1,                   (gr_funcptr) _gr_gr_series_agm1},
    {GR_METHOD_ELLIPTIC_K,             (gr_funcptr) _gr_gr_series_elliptic_k},
    {GR_METHOD_WEIERSTRASS_P,          (gr_funcptr) _gr_gr_series_weierstrass_p},

    {0,                     (gr_funcptr) NULL},
};

void
gr_ctx_init_gr_series(gr_ctx_t ctx, gr_ctx_t base_ring, slong prec)
{
    ctx->which_ring = GR_CTX_GR_SERIES;
    ctx->sizeof_elem = sizeof(gr_series_struct);
    ctx->size_limit = WORD_MAX;

    SERIES_CTX(ctx)->base_ring = (gr_ctx_struct *) base_ring;
    SERIES_CTX(ctx)->var = (char *) default_var;
    SERIES_CTX(ctx)->sctx.mod = SERIES_ERR_EXACT;
    SERIES_CTX(ctx)->sctx.prec = FLINT_MIN(FLINT_MAX(0, prec), SERIES_ERR_MAX);

    ctx->methods = _gr_series_methods;

    if (!_gr_series_methods_initialized)
    {
        gr_method_tab_init(_gr_series_methods, _gr_series_methods_input);
        _gr_series_methods_initialized = 1;
    }
}

void
gr_ctx_init_gr_series_mod(gr_ctx_t ctx, gr_ctx_t base_ring, slong mod)
{
    ctx->which_ring = GR_CTX_GR_SERIES_MOD;
    ctx->sizeof_elem = sizeof(gr_series_struct);
    ctx->size_limit = WORD_MAX;

    if (mod >= SERIES_ERR_EXACT)
        flint_throw(FLINT_ERROR, "(%s)\n", __func__);

    SERIES_CTX(ctx)->base_ring = (gr_ctx_struct *) base_ring;
    SERIES_CTX(ctx)->var = (char *) default_var;
    SERIES_CTX(ctx)->sctx.mod = FLINT_MAX(0, mod);
    SERIES_CTX(ctx)->sctx.prec = mod;

    ctx->methods = _gr_series_methods;

    if (!_gr_series_methods_initialized)
    {
        gr_method_tab_init(_gr_series_methods, _gr_series_methods_input);
        _gr_series_methods_initialized = 1;
    }
}
