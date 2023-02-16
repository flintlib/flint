/*
    Copyright (C) 2023 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "gr_vec.h"
#include "gr_poly.h"

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
    gr_poly_write(out, &x->poly, cctx);

    if (x->error != SERIES_ERR_EXACT)
    {
        gr_stream_write(out, " + O(x^");
        gr_stream_write_si(out, x->error);
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

        status |= gr_poly_truncate(&res->poly, trunc, cctx);
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

        status |= gr_poly_truncate(&res->poly, trunc, cctx);
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

/*

from flint2 import *; R = PowerSeriesRing_gr_series(ZZ, 10); x = R.gen(); (1+x)**20


*/



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
static truth_t _gr_gr_series_equal(const gr_series_t x, const gr_series_t y, gr_ctx_t ctx) { return gr_series_equal(x, y, SERIES_SCTX(ctx), SERIES_ELEM_CTX(ctx)); }
static int _gr_gr_series_neg(gr_series_t res, const gr_series_t x, gr_ctx_t ctx) { return gr_series_neg(res, x, SERIES_SCTX(ctx), SERIES_ELEM_CTX(ctx)); }
static int _gr_gr_series_add(gr_series_t res, const gr_series_t x, const gr_series_t y, gr_ctx_t ctx) { return gr_series_add(res, x, y, SERIES_SCTX(ctx), SERIES_ELEM_CTX(ctx)); }
static int _gr_gr_series_sub(gr_series_t res, const gr_series_t x, const gr_series_t y, gr_ctx_t ctx) { return gr_series_sub(res, x, y, SERIES_SCTX(ctx), SERIES_ELEM_CTX(ctx)); }
static int _gr_gr_series_mul(gr_series_t res, const gr_series_t x, const gr_series_t y, gr_ctx_t ctx) { return gr_series_mul(res, x, y, SERIES_SCTX(ctx), SERIES_ELEM_CTX(ctx)); }


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
/*
    else if (x_ctx->which_ring == GR_CTX_GR_POLY && !strcmp(POLYNOMIAL_CTX(x_ctx)->var, POLYNOMIAL_CTX(ctx)->var))
    {
        return gr_poly_set_gr_poly_other(res, x, POLYNOMIAL_ELEM_CTX(x_ctx), POLYNOMIAL_ELEM_CTX(ctx));
    }
    else if (x_ctx->which_ring == GR_CTX_GR_VEC)
    {
        gr_poly_t tmp;
        tmp->coeffs = ((gr_vec_struct *) x)->entries;
        tmp->length = ((gr_vec_struct *) x)->length;

        return gr_poly_set_gr_poly_other(res, tmp, VECTOR_CTX(x_ctx)->base_ring, POLYNOMIAL_ELEM_CTX(ctx));
    }
    else
    {
        int status = GR_SUCCESS;

        gr_poly_fit_length(res, 1, POLYNOMIAL_ELEM_CTX(ctx));
        status = gr_set_other(res->coeffs, x, x_ctx, POLYNOMIAL_ELEM_CTX(ctx));
        if (status == GR_SUCCESS)
        {
            _gr_poly_set_length(res, 1, POLYNOMIAL_ELEM_CTX(ctx));
            _gr_poly_normalise(res, POLYNOMIAL_ELEM_CTX(ctx));
        }
        else
            _gr_poly_set_length(res, 0, POLYNOMIAL_ELEM_CTX(ctx));
        return status;
    }
*/

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
        flint_abort();

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
