/*
    Copyright (C) 2023, 2025 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <string.h>
#include "fmpq.h"
#include "gr_vec.h"
#include "gr_poly.h"
#include "gr_generic.h"
#include "gr_series.h"


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

void
gr_series_init(gr_series_t res, gr_ctx_t ctx)
{
    gr_poly_init(GR_SERIES_POLY(res), GR_SERIES_ELEM_CTX(ctx));
    GR_SERIES_ERROR(res) = GR_SERIES_ERR_EXACT;
}

void
gr_series_clear(gr_series_t res, gr_ctx_t ctx)
{
    gr_poly_clear(GR_SERIES_POLY(res), GR_SERIES_ELEM_CTX(ctx));
}

int
gr_series_zero(gr_series_t res, gr_ctx_t ctx)
{
    GR_SERIES_ERROR(res) = GR_SERIES_ERR_EXACT;
    return gr_poly_zero(GR_SERIES_POLY(res), GR_SERIES_ELEM_CTX(ctx));
}

void
gr_series_swap(gr_series_t x, gr_series_t y, gr_ctx_t ctx)
{
    gr_series_t tmp;
    *tmp = *x;
    *x = *y;
    *y = *tmp;
}

void
gr_series_set_shallow(gr_series_t res, const gr_series_t x, gr_ctx_t ctx)
{
    *res = *x;
}

slong
_gr_series_get_error(const gr_series_t x, gr_ctx_t ctx)
{
    return GR_SERIES_ERROR(x);
}

truth_t
_gr_series_is_exact(const gr_series_t x, gr_ctx_t ctx)
{
    return (GR_SERIES_ERROR(x) == GR_SERIES_ERR_EXACT) ? T_TRUE : T_FALSE;
}

void
_gr_series_set_error(gr_series_t x, slong err, gr_ctx_t ctx)
{
    err = FLINT_MAX(err, 0);
    err = FLINT_MIN(err, GR_SERIES_ERR_MAX);

    GR_SERIES_ERROR(x) = err;

    if (GR_SERIES_POLY(x)->length > err)
        GR_MUST_SUCCEED(gr_poly_truncate(GR_SERIES_POLY(x), GR_SERIES_POLY(x), err, GR_SERIES_ELEM_CTX(ctx)));
}

void
_gr_series_make_exact(gr_series_t x, gr_ctx_t ctx)
{
    GR_SERIES_ERROR(x) = GR_SERIES_ERR_EXACT;
}

int
gr_series_randtest(gr_series_t res, flint_rand_t state, gr_ctx_t ctx)
{
    slong len;

    len = n_randint(state, 5);
    len = FLINT_MAX(len, 0);
    len = FLINT_MIN(len, GR_SERIES_ERR_MAX);

    int status = gr_poly_randtest(GR_SERIES_POLY(res), state, len, GR_SERIES_ELEM_CTX(ctx));

    if (n_randint(state, 2))
        _gr_series_make_exact(res, ctx);
    else
        _gr_series_set_error(res, n_randint(state, len + 1), ctx);

    return status;
}

int
gr_series_write(gr_stream_t out, const gr_series_t x, gr_ctx_t ctx)
{
    const char * var = GR_SERIES_CTX(ctx)->var;

    gr_poly_write(out, GR_SERIES_POLY(x), var, GR_SERIES_ELEM_CTX(ctx));

    if (GR_SERIES_ERROR(x) != GR_SERIES_ERR_EXACT)
    {
        gr_stream_write(out, " + O(");
        gr_stream_write(out, var);
        gr_stream_write(out, "^");
        gr_stream_write_si(out, GR_SERIES_ERROR(x));
        gr_stream_write(out, ")");
    }

    return GR_SUCCESS;
}

int
gr_series_one(gr_series_t res, gr_ctx_t ctx)
{
    if (GR_SERIES_PREC(ctx) == 0)
    {
        GR_SERIES_ERROR(res) = 0;
        return gr_poly_zero(GR_SERIES_POLY(res), GR_SERIES_ELEM_CTX(ctx));
    }
    else
    {
        GR_SERIES_ERROR(res) = GR_SERIES_ERR_EXACT;
        return gr_poly_one(GR_SERIES_POLY(res), GR_SERIES_ELEM_CTX(ctx));
    }
}

/* todo: truncate without set */
int
gr_series_set(gr_series_t res, const gr_series_t x, gr_ctx_t ctx)
{
    int status = GR_SUCCESS;
    slong len, trunc;

    GR_SERIES_ERROR(res) = GR_SERIES_ERROR(x);
    status |= gr_poly_set(GR_SERIES_POLY(res), GR_SERIES_POLY(x), GR_SERIES_ELEM_CTX(ctx));

    trunc = FLINT_MIN(GR_SERIES_PREC(ctx), GR_SERIES_ERROR(res));
    len = GR_SERIES_POLY(res)->length;

    if (len > trunc)
    {
        if (len > GR_SERIES_PREC(ctx))
            GR_SERIES_ERROR(res) = FLINT_MIN(GR_SERIES_ERROR(res), GR_SERIES_PREC(ctx));

        status |= gr_poly_truncate(GR_SERIES_POLY(res), GR_SERIES_POLY(res), trunc, GR_SERIES_ELEM_CTX(ctx));
    }

    return status;
}

int
gr_series_gen(gr_series_t res, gr_ctx_t ctx)
{
    if (GR_SERIES_PREC(ctx) >= 2)
    {
        GR_SERIES_ERROR(res) = GR_SERIES_ERR_EXACT;
        return gr_poly_gen(GR_SERIES_POLY(res), GR_SERIES_ELEM_CTX(ctx));
    }
    else
    {
        GR_SERIES_ERROR(res) = GR_SERIES_PREC(ctx);
        return gr_poly_zero(GR_SERIES_POLY(res), GR_SERIES_ELEM_CTX(ctx));
    }
}

/* todo: truncate before neg */
int
gr_series_neg(gr_series_t res, const gr_series_t x, gr_ctx_t ctx)
{
    int status = GR_SUCCESS;
    slong len, trunc;

    GR_SERIES_ERROR(res) = GR_SERIES_ERROR(x);
    status |= gr_poly_neg(GR_SERIES_POLY(res), GR_SERIES_POLY(x), GR_SERIES_ELEM_CTX(ctx));

    trunc = FLINT_MIN(GR_SERIES_PREC(ctx), GR_SERIES_ERROR(res));
    len = GR_SERIES_POLY(res)->length;

    if (len > trunc)
    {
        if (len > GR_SERIES_PREC(ctx))
            GR_SERIES_ERROR(res) = FLINT_MIN(GR_SERIES_ERROR(res), GR_SERIES_PREC(ctx));

        status |= gr_poly_truncate(GR_SERIES_POLY(res), GR_SERIES_POLY(res), trunc, GR_SERIES_ELEM_CTX(ctx));
    }

    return status;
}


int
gr_series_set_gr_poly(gr_series_t res, const gr_poly_t x, gr_ctx_t ctx)
{
    gr_series_t tmp;
    tmp->poly = *x;
    tmp->error = GR_SERIES_ERR_EXACT;

    return gr_series_set(res, tmp, ctx);
}

int
gr_series_set_scalar(gr_series_t res, gr_srcptr x, gr_ctx_t ctx)
{
    if (gr_is_zero(x, GR_SERIES_ELEM_CTX(ctx)) == T_TRUE)
    {
        return gr_series_zero(res, ctx);
    }
    else
    {
        gr_series_t tmp;
        tmp->poly.coeffs = (gr_ptr) x;
        tmp->poly.length = 1;
        tmp->poly.alloc = 1;
        tmp->error = GR_SERIES_ERR_EXACT;
        return gr_series_set(res, tmp, ctx);
    }
}

/* todo: optimize */
int
gr_series_set_si(gr_series_t res, slong c, gr_ctx_t ctx)
{
    if (c == 0)
    {
        return gr_series_zero(res, ctx);
    }
    else
    {
        gr_ptr t;
        int status = GR_SUCCESS;
        GR_TMP_INIT(t, GR_SERIES_ELEM_CTX(ctx));
        status |= gr_set_si(t, c, GR_SERIES_ELEM_CTX(ctx));
        status |= gr_series_set_scalar(res, t, ctx);
        GR_TMP_CLEAR(t, GR_SERIES_ELEM_CTX(ctx));
        return status;
    }
}

/* todo: optimize */
int
gr_series_set_ui(gr_series_t res, ulong c, gr_ctx_t ctx)
{
    if (c == 0)
    {
        return gr_series_zero(res, ctx);
    }
    else
    {
        gr_ptr t;
        int status = GR_SUCCESS;
        GR_TMP_INIT(t, GR_SERIES_ELEM_CTX(ctx));
        status |= gr_set_ui(t, c, GR_SERIES_ELEM_CTX(ctx));
        status |= gr_series_set_scalar(res, t, ctx);
        GR_TMP_CLEAR(t, GR_SERIES_ELEM_CTX(ctx));
        return status;
    }
}

int
gr_series_set_fmpz(gr_series_t res, const fmpz_t c, gr_ctx_t ctx)
{
    if (fmpz_is_zero(c))
    {
        return gr_series_zero(res, ctx);
    }
    else
    {
        gr_ptr t;
        int status = GR_SUCCESS;
        GR_TMP_INIT(t, GR_SERIES_ELEM_CTX(ctx));
        status |= gr_set_fmpz(t, c, GR_SERIES_ELEM_CTX(ctx));
        status |= gr_series_set_scalar(res, t, ctx);
        GR_TMP_CLEAR(t, GR_SERIES_ELEM_CTX(ctx));
        return status;
    }
}

int
gr_series_set_fmpq(gr_series_t res, const fmpq_t c, gr_ctx_t ctx)
{
    if (fmpq_is_zero(c))
    {
        return gr_series_zero(res, ctx);
    }
    else
    {
        gr_ptr t;
        int status = GR_SUCCESS;
        GR_TMP_INIT(t, GR_SERIES_ELEM_CTX(ctx));
        status |= gr_set_fmpq(t, c, GR_SERIES_ELEM_CTX(ctx));
        status |= gr_series_set_scalar(res, t, ctx);
        GR_TMP_CLEAR(t, GR_SERIES_ELEM_CTX(ctx));
        return status;
    }
}


truth_t
gr_series_is_zero(const gr_series_t x, gr_ctx_t ctx)
{
    truth_t is_zero;
    slong err = GR_SERIES_ERROR(x);
    slong len = GR_SERIES_POLY(x)->length;

    if (err == GR_SERIES_ERR_EXACT)
        return gr_poly_is_zero(GR_SERIES_POLY(x), GR_SERIES_ELEM_CTX(ctx));

    is_zero = _gr_vec_is_zero(GR_SERIES_POLY(x)->coeffs, FLINT_MIN(err, len), GR_SERIES_ELEM_CTX(ctx));

    if (is_zero == T_FALSE)
        return T_FALSE;

    return T_UNKNOWN;
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
gr_series_is_one(const gr_series_t x, gr_ctx_t ctx)
{
    truth_t is_zero, is_one;
    slong xlen = GR_SERIES_POLY(x)->length;

    if (GR_SERIES_ERROR(x) <= 0)
        return T_UNKNOWN;

    if (xlen == 0)
        return gr_ctx_is_zero_ring(GR_SERIES_ELEM_CTX(ctx));

    is_one = gr_is_one(GR_SERIES_POLY(x)->coeffs, GR_SERIES_ELEM_CTX(ctx));

    if (is_one == T_FALSE)
        return T_FALSE;

    if (xlen >= 2)
    {
        is_zero = _gr_vec_is_zero(GR_ENTRY(GR_SERIES_POLY(x)->coeffs, 1, GR_SERIES_ELEM_CTX(ctx)->sizeof_elem), FLINT_MIN(xlen - 1, GR_SERIES_ERROR(x) - 1), GR_SERIES_ELEM_CTX(ctx));
        if (is_zero == T_FALSE)
            return T_FALSE;
    }
    else
    {
        is_zero = T_TRUE;
    }

    if (GR_SERIES_ERROR(x) == GR_SERIES_ERR_EXACT && is_zero == T_TRUE && is_one == T_TRUE)
        return T_TRUE;

    return T_UNKNOWN;
}


truth_t
gr_series_equal(const gr_series_t x, const gr_series_t y, gr_ctx_t ctx)
{
    truth_t equal;
    slong len, xlen, ylen, xerr, yerr, err;

    xlen = GR_SERIES_POLY(x)->length;
    ylen = GR_SERIES_POLY(y)->length;
    xerr = GR_SERIES_ERROR(x);
    yerr = GR_SERIES_ERROR(y);
    err = FLINT_MIN(xerr, yerr);

    len = FLINT_MAX(xlen, ylen);
    len = FLINT_MIN(len, err);

    /* compare exactly or within the error bound */
    if (xlen >= ylen)
        equal = _gr_poly_equal2(x->poly.coeffs, FLINT_MIN(xlen, len), y->poly.coeffs, FLINT_MIN(ylen, len), GR_SERIES_ELEM_CTX(ctx));
    else
        equal = _gr_poly_equal2(y->poly.coeffs, FLINT_MIN(ylen, len), x->poly.coeffs, FLINT_MIN(xlen, len), GR_SERIES_ELEM_CTX(ctx));

    if (err == GR_SERIES_ERR_EXACT && equal == T_TRUE)
        return T_TRUE;

    if (equal == T_FALSE)
        return T_FALSE;

    return T_UNKNOWN;
}

int
gr_series_add(gr_series_t res, const gr_series_t x, const gr_series_t y, gr_ctx_t ctx)
{
    slong len, xlen, ylen, xerr, yerr, err;
    int status = GR_SUCCESS;

    xlen = GR_SERIES_POLY(x)->length;
    ylen = GR_SERIES_POLY(y)->length;
    xerr = GR_SERIES_ERROR(x);
    yerr = GR_SERIES_ERROR(y);
    err = FLINT_MIN(xerr, yerr);

    /* length of the polynomial sum, ignoring errors */
    len = FLINT_MAX(xlen, ylen);

    /* the result will be truncated */
    if (len > GR_SERIES_PREC(ctx))
        err = FLINT_MIN(err, GR_SERIES_PREC(ctx));

    len = FLINT_MIN(len, GR_SERIES_PREC(ctx));
    len = FLINT_MIN(len, err);

    GR_SERIES_ERROR(res) = err;
    status |= gr_poly_add_series(GR_SERIES_POLY(res), GR_SERIES_POLY(x), GR_SERIES_POLY(y), len, GR_SERIES_ELEM_CTX(ctx));
    return status;
}

int
gr_series_sub(gr_series_t res, const gr_series_t x, const gr_series_t y, gr_ctx_t ctx)
{
    slong len, xlen, ylen, xerr, yerr, err;
    int status = GR_SUCCESS;

    xlen = GR_SERIES_POLY(x)->length;
    ylen = GR_SERIES_POLY(y)->length;
    xerr = GR_SERIES_ERROR(x);
    yerr = GR_SERIES_ERROR(y);
    err = FLINT_MIN(xerr, yerr);

    /* length of the polynomial sum, ignoring errors */
    len = FLINT_MAX(xlen, ylen);

    /* the result will be truncated */
    if (len > GR_SERIES_PREC(ctx))
        err = FLINT_MIN(err, GR_SERIES_PREC(ctx));

    len = FLINT_MIN(len, GR_SERIES_PREC(ctx));
    len = FLINT_MIN(len, err);

    GR_SERIES_ERROR(res) = err;
    status |= gr_poly_sub_series(GR_SERIES_POLY(res), GR_SERIES_POLY(x), GR_SERIES_POLY(y), len, GR_SERIES_ELEM_CTX(ctx));
    return status;
}

int
gr_series_mul(gr_series_t res, const gr_series_t x, const gr_series_t y, gr_ctx_t ctx)
{
    slong len, xlen, ylen, xerr, yerr, err;
    int status = GR_SUCCESS;

    xlen = GR_SERIES_POLY(x)->length;
    ylen = GR_SERIES_POLY(y)->length;
    xerr = GR_SERIES_ERROR(x);
    yerr = GR_SERIES_ERROR(y);
    err = FLINT_MIN(xerr, yerr);

    if (xlen == 0 && xerr == GR_SERIES_ERR_EXACT)
        return gr_series_zero(res, ctx);

    if (ylen == 0 && yerr == GR_SERIES_ERR_EXACT)
        return gr_series_zero(res, ctx);

    /* length of the polynomial product, ignoring errors */
    if (xlen == 0 || ylen == 0)
        len = 0;
    else
        len = xlen + ylen - 1;

    /* the result will be truncated */
    if (len > GR_SERIES_PREC(ctx))
        err = FLINT_MIN(err, GR_SERIES_PREC(ctx));

    len = FLINT_MIN(len, GR_SERIES_PREC(ctx));
    len = FLINT_MIN(len, err);

    GR_SERIES_ERROR(res) = err;
    status |= gr_poly_mullow(GR_SERIES_POLY(res), GR_SERIES_POLY(x), GR_SERIES_POLY(y), len, GR_SERIES_ELEM_CTX(ctx));
    return status;
}

int
gr_series_inv(gr_series_t res, const gr_series_t x, gr_ctx_t ctx)
{
    slong len, xlen, xerr, err;
    int status = GR_SUCCESS;

    xlen = GR_SERIES_POLY(x)->length;
    xerr = GR_SERIES_ERROR(x);
    err = xerr;

    if (xlen == 0 && xerr == GR_SERIES_ERR_EXACT)
    {
        truth_t zero = gr_ctx_is_zero_ring(GR_SERIES_ELEM_CTX(ctx));

        if (zero == T_TRUE)
            return gr_series_zero(res, ctx);
        if (zero == T_UNKNOWN)
            return GR_UNABLE;

        return GR_DOMAIN;
    }

    if (xlen == 0 || xerr == 0)
        return GR_UNABLE;

    len = FLINT_MIN(GR_SERIES_PREC(ctx), err);
    err = len;

    GR_SERIES_ERROR(res) = err;
    status |= gr_poly_inv_series(GR_SERIES_POLY(res), GR_SERIES_POLY(x), len, GR_SERIES_ELEM_CTX(ctx));
    return status;
}

truth_t
gr_series_coeff_is_zero(const gr_series_t x, slong i, gr_ctx_t ctx)
{
    if (i >= GR_SERIES_ERROR(x))
        return T_UNKNOWN;

    if (i >= x->poly.length)
        return T_TRUE;

    if (i < 0)
        return T_TRUE;

    return gr_is_zero(gr_poly_coeff_srcptr(GR_SERIES_POLY(x), i, GR_SERIES_ELEM_CTX(ctx)), GR_SERIES_ELEM_CTX(ctx));
}

/* Todo: optimizations for len == 1 denominator */
/* Todo: user gr_poly_series_divexact when there is one (currently only basecase) */
int
_gr_series_div(gr_series_t res, const gr_series_t x, const gr_series_t y, int divexact, gr_ctx_t ctx)
{
    slong len, xlen, ylen, xerr, yerr, err;
    int status = GR_SUCCESS;
    truth_t is_zero;
    slong val;

    xlen = GR_SERIES_POLY(x)->length;
    ylen = GR_SERIES_POLY(y)->length;
    xerr = GR_SERIES_ERROR(x);
    yerr = GR_SERIES_ERROR(y);

    /* Divide out common leading zeros. */
    val = 0;
    while (val < FLINT_MAX(xlen, ylen))
    {
        is_zero = gr_series_coeff_is_zero(y, val, ctx);
        if (is_zero == T_UNKNOWN)
            return GR_UNABLE;
        if (is_zero == T_FALSE)
            break;

        is_zero = gr_series_coeff_is_zero(x, val, ctx);
        if (is_zero == T_UNKNOWN)
            return GR_UNABLE;
        if (is_zero == T_FALSE)
            return GR_DOMAIN;

        val++;
    }

    xlen = FLINT_MAX(xlen - val, 0);
    ylen = FLINT_MAX(ylen - val, 0);

    if (xerr != GR_SERIES_ERR_EXACT) xerr -= val;
    if (yerr != GR_SERIES_ERR_EXACT) yerr -= val;

    /* x / 0 = GR_DOMAIN, unless we are in the zero ring. */
    if (ylen == 0 && yerr == GR_SERIES_ERR_EXACT)
    {
        is_zero = gr_ctx_is_zero_ring(GR_SERIES_ELEM_CTX(ctx));

        if (is_zero == T_FALSE)
            return GR_DOMAIN;
        else if (is_zero == T_UNKNOWN)
            return GR_UNABLE;
        else
            return GR_SUCCESS;
    }

    /* If y = 0 + O(t^n), we do not know anything. */
    if (ylen <= 0 || yerr == 0)
        return GR_UNABLE;

    /* If not over a field, we must have an invertible denominator to guarantee
       that the quotient exists, unless this is a polynomial division. */
    if (!divexact && gr_ctx_is_field(GR_SERIES_ELEM_CTX(ctx)) != T_TRUE)
    {
        gr_srcptr y0 = gr_poly_coeff_srcptr(GR_SERIES_POLY(y), val, GR_SERIES_ELEM_CTX(ctx));

        if (gr_is_invertible(y0, GR_SERIES_ELEM_CTX(ctx)) != T_TRUE)
        {
            if (xerr == GR_SERIES_ERR_EXACT && yerr == GR_SERIES_ERR_EXACT)
            {
                gr_poly_t xb, yb;
                gr_poly_t q, r;
                int status;

                xb->coeffs = (gr_ptr) gr_poly_coeff_srcptr(GR_SERIES_POLY(x), val, GR_SERIES_ELEM_CTX(ctx));
                xb->length = xlen;
                yb->coeffs = (gr_ptr) gr_poly_coeff_srcptr(GR_SERIES_POLY(y), val, GR_SERIES_ELEM_CTX(ctx));
                yb->length = ylen;

                gr_poly_init(q, GR_SERIES_ELEM_CTX(ctx));
                gr_poly_init(r, GR_SERIES_ELEM_CTX(ctx));

                status = gr_poly_divrem(q, r, xb, yb, GR_SERIES_ELEM_CTX(ctx));

                if (status == GR_SUCCESS && gr_poly_is_zero(r, GR_SERIES_ELEM_CTX(ctx)) == T_TRUE)
                {
                    gr_poly_swap(GR_SERIES_POLY(res), q, GR_SERIES_ELEM_CTX(ctx));
                    GR_SERIES_ERROR(res) = GR_SERIES_ERR_EXACT;
                }
                else if (xb->length == 1 && yb->length == 1 && gr_poly_is_zero(r, GR_SERIES_ELEM_CTX(ctx)) == T_FALSE)
                {
                    status = GR_DOMAIN;
                }
                else
                {
                    status = GR_UNABLE;
                }

                gr_poly_clear(q, GR_SERIES_ELEM_CTX(ctx));
                gr_poly_clear(r, GR_SERIES_ELEM_CTX(ctx));

                return status;
            }

            return GR_UNABLE;
        }
    }

    err = FLINT_MIN(xerr, yerr);
    err = FLINT_MIN(err, GR_SERIES_PREC(ctx));

    if (err <= 0)
    {
        status |= gr_poly_zero(GR_SERIES_POLY(res), GR_SERIES_ELEM_CTX(ctx));
        GR_SERIES_ERROR(res) = 0;
        return status;
    }

    len = FLINT_MIN(GR_SERIES_PREC(ctx), err);

    /* If we have an exact polynomial division, see if we have an exact result. */
    /* TODO: Special case for length 1 divisor. */
    /* TODO: If the polynomials are short compared to the precision, try
             a checked polynomial division before doing a power series division. */
    if (xerr == GR_SERIES_ERR_EXACT && yerr == GR_SERIES_ERR_EXACT)
    {
        gr_poly_t t, u;

        gr_poly_init(t, GR_SERIES_ELEM_CTX(ctx));
        gr_poly_init(u, GR_SERIES_ELEM_CTX(ctx));

        if (val == 0)
        {
            /* Compute one extra coefficient to quickly discard most inexact divisions */
            status |= gr_poly_div_series(t, GR_SERIES_POLY(x), GR_SERIES_POLY(y), len + 1, GR_SERIES_ELEM_CTX(ctx));

            if (t->length <= len)
            {
                if (status == GR_SUCCESS)
                {
                    status |= gr_poly_mul(u, t, GR_SERIES_POLY(y), GR_SERIES_ELEM_CTX(ctx));
                    if (gr_poly_equal(GR_SERIES_POLY(x), u, GR_SERIES_ELEM_CTX(ctx)) == T_TRUE)
                        err = GR_SERIES_ERR_EXACT;
                }
            }
        }
        else
        {
            gr_poly_t xb, yb;

            xb->coeffs = (gr_ptr) gr_poly_coeff_srcptr(GR_SERIES_POLY(x), val, GR_SERIES_ELEM_CTX(ctx));
            xb->length = xlen;
            yb->coeffs = (gr_ptr) gr_poly_coeff_srcptr(GR_SERIES_POLY(y), val, GR_SERIES_ELEM_CTX(ctx));
            yb->length = ylen;

            status |= gr_poly_div_series(t, xb, yb, len + 1, GR_SERIES_ELEM_CTX(ctx));

            if (t->length <= len)
            {
                if (status == GR_SUCCESS)
                {
                    status |= gr_poly_mul(u, t, yb, GR_SERIES_ELEM_CTX(ctx));
                    if (gr_poly_equal(xb, u, GR_SERIES_ELEM_CTX(ctx)) == T_TRUE)
                        err = GR_SERIES_ERR_EXACT;
                }
            }
        }

        status |= gr_poly_truncate(t, t, len, GR_SERIES_ELEM_CTX(ctx));

        gr_poly_swap(GR_SERIES_POLY(res), t, GR_SERIES_ELEM_CTX(ctx));
        GR_SERIES_ERROR(res) = err;

        gr_poly_clear(t, GR_SERIES_ELEM_CTX(ctx));
        gr_poly_clear(u, GR_SERIES_ELEM_CTX(ctx));
        return status;
    }

    if (val == 0)
    {
        GR_SERIES_ERROR(res) = err;
        status |= gr_poly_div_series(GR_SERIES_POLY(res), GR_SERIES_POLY(x), GR_SERIES_POLY(y), len, GR_SERIES_ELEM_CTX(ctx));
        return status;
    }
    else
    {
        gr_poly_t xb, yb;

        xb->coeffs = (gr_ptr) gr_poly_coeff_srcptr(GR_SERIES_POLY(x), val, GR_SERIES_ELEM_CTX(ctx));
        xb->length = xlen;
        yb->coeffs = (gr_ptr) gr_poly_coeff_srcptr(GR_SERIES_POLY(y), val, GR_SERIES_ELEM_CTX(ctx));
        yb->length = ylen;

        if (x == res || y == res)
        {
            gr_series_t t;
            gr_series_init(t, ctx);
            t->error = err;
            status |= gr_poly_div_series(GR_SERIES_POLY(t), xb, yb, len, GR_SERIES_ELEM_CTX(ctx));
            gr_series_swap(res, t, ctx);
            gr_series_clear(t, ctx);
            return status;
        }
        else
        {
            GR_SERIES_ERROR(res) = err;
            status |= gr_poly_div_series(GR_SERIES_POLY(res), xb, yb, len, GR_SERIES_ELEM_CTX(ctx));
            return status;
        }
    }
}

int
gr_series_div(gr_series_t res, const gr_series_t x, const gr_series_t y, gr_ctx_t ctx)
{
    return _gr_series_div(res, x, y, 0, ctx);
}

int
gr_series_divexact(gr_series_t res, const gr_series_t x, const gr_series_t y, gr_ctx_t ctx)
{
    return _gr_series_div(res, x, y, 1, ctx);
}

int
gr_series_sqrt(gr_series_t res, const gr_series_t x, gr_ctx_t ctx)
{
    slong len, xlen, xerr, err;
    int status = GR_SUCCESS;

    xlen = GR_SERIES_POLY(x)->length;
    xerr = GR_SERIES_ERROR(x);
    err = xerr;

    if (xlen == 0 && xerr == GR_SERIES_ERR_EXACT)
        return gr_series_zero(res, ctx);

    if (xlen == 0 || xerr == 0)
        return GR_UNABLE;

    if (xlen == 1 && xerr == GR_SERIES_ERR_EXACT)
    {
        gr_poly_fit_length(GR_SERIES_POLY(res), 1, GR_SERIES_ELEM_CTX(ctx));

        status = gr_sqrt(GR_SERIES_POLY(res)->coeffs, GR_SERIES_POLY(x)->coeffs, GR_SERIES_ELEM_CTX(ctx));

        if (status == GR_SUCCESS && GR_SERIES_PREC(ctx) != 0)
        {
            _gr_poly_set_length(GR_SERIES_POLY(res), 1, GR_SERIES_ELEM_CTX(ctx));
            _gr_poly_normalise(GR_SERIES_POLY(res), GR_SERIES_ELEM_CTX(ctx));
            GR_SERIES_ERROR(res) = GR_SERIES_ERR_EXACT;
        }
        else
        {
            _gr_poly_set_length(GR_SERIES_POLY(res), 0, GR_SERIES_ELEM_CTX(ctx));
            GR_SERIES_ERROR(res) = 0;
        }

        return status;
    }

    /* todo: handle even valuations; exact polynomial square roots */

    if (gr_ctx_is_field(GR_SERIES_ELEM_CTX(ctx)) != T_TRUE)
        return GR_UNABLE;

    /* todo: function to check for odd characteristic */
    if (gr_ctx_is_finite_characteristic(GR_SERIES_ELEM_CTX(ctx)) != T_FALSE)
    {
        gr_ptr t;
        GR_TMP_INIT(t, GR_SERIES_ELEM_CTX(ctx));

        status = gr_set_ui(t, 2, GR_SERIES_ELEM_CTX(ctx));
        if (status == GR_SUCCESS)
            status = (gr_is_invertible(t, GR_SERIES_ELEM_CTX(ctx)) == T_TRUE) ? GR_SUCCESS : GR_UNABLE;

        GR_TMP_CLEAR(t, GR_SERIES_ELEM_CTX(ctx));

        if (status != GR_SUCCESS)
            return GR_UNABLE;
    }

    if (xlen > 1 && gr_is_zero(GR_SERIES_POLY(x)->coeffs, GR_SERIES_ELEM_CTX(ctx)) != T_FALSE)
        return GR_UNABLE;

    len = FLINT_MIN(GR_SERIES_PREC(ctx), err);
    err = len;

    GR_SERIES_ERROR(res) = err;
    status |= gr_poly_sqrt_series(GR_SERIES_POLY(res), GR_SERIES_POLY(x), len, GR_SERIES_ELEM_CTX(ctx));
    return status;
}

int
gr_series_rsqrt(gr_series_t res, const gr_series_t x, gr_ctx_t ctx)
{
    slong len, xlen, xerr, err;
    int status = GR_SUCCESS;

    xlen = GR_SERIES_POLY(x)->length;
    xerr = GR_SERIES_ERROR(x);
    err = xerr;

    if (gr_ctx_is_zero_ring(GR_SERIES_ELEM_CTX(ctx)) == T_TRUE)
        return gr_series_zero(res, ctx);

    if (xerr == 0)
        return GR_UNABLE;

    if (xlen == 0)
        return GR_DOMAIN;

    if (xlen == 1 && xerr == GR_SERIES_ERR_EXACT)
    {
        gr_poly_fit_length(GR_SERIES_POLY(res), 1, GR_SERIES_ELEM_CTX(ctx));

        status = gr_rsqrt(GR_SERIES_POLY(res)->coeffs, GR_SERIES_POLY(x)->coeffs, GR_SERIES_ELEM_CTX(ctx));

        if (status == GR_SUCCESS && GR_SERIES_PREC(ctx) != 0)
        {
            _gr_poly_set_length(GR_SERIES_POLY(res), 1, GR_SERIES_ELEM_CTX(ctx));
            _gr_poly_normalise(GR_SERIES_POLY(res), GR_SERIES_ELEM_CTX(ctx));
            GR_SERIES_ERROR(res) = GR_SERIES_ERR_EXACT;
        }
        else
        {
            _gr_poly_set_length(GR_SERIES_POLY(res), 0, GR_SERIES_ELEM_CTX(ctx));
            GR_SERIES_ERROR(res) = 0;
        }

        return status;
    }

    /* todo: handle even valuations; exact polynomial square roots */

    if (gr_ctx_is_field(GR_SERIES_ELEM_CTX(ctx)) != T_TRUE)
        return GR_UNABLE;

    /* todo: function to check for odd characteristic */
    if (gr_ctx_is_finite_characteristic(GR_SERIES_ELEM_CTX(ctx)) != T_FALSE)
    {
        gr_ptr t;
        GR_TMP_INIT(t, GR_SERIES_ELEM_CTX(ctx));

        status = gr_set_ui(t, 2, GR_SERIES_ELEM_CTX(ctx));
        if (status == GR_SUCCESS)
            status = (gr_is_invertible(t, GR_SERIES_ELEM_CTX(ctx)) == T_TRUE) ? GR_SUCCESS : GR_UNABLE;

        GR_TMP_CLEAR(t, GR_SERIES_ELEM_CTX(ctx));

        if (status != GR_SUCCESS)
            return GR_UNABLE;
    }

    if (xlen > 1 && gr_is_zero(GR_SERIES_POLY(x)->coeffs, GR_SERIES_ELEM_CTX(ctx)) != T_FALSE)
        return GR_UNABLE;

    len = FLINT_MIN(GR_SERIES_PREC(ctx), err);
    err = len;

    GR_SERIES_ERROR(res) = err;
    status |= gr_poly_rsqrt_series(GR_SERIES_POLY(res), GR_SERIES_POLY(x), len, GR_SERIES_ELEM_CTX(ctx));
    return status;
}



#define UNARY_POLY_WRAPPER(func) \
int \
gr_series_ ## func(gr_series_t res, const gr_series_t x, gr_ctx_t ctx) \
{ \
    slong len, xlen, xerr, err; \
    int status = GR_SUCCESS; \
 \
    if (gr_ctx_is_rational_vector_space(GR_SERIES_ELEM_CTX(ctx)) != T_TRUE) \
        return GR_UNABLE; \
\
    xlen = GR_SERIES_POLY(x)->length; \
    xerr = GR_SERIES_ERROR(x); \
    err = xerr; \
 \
    len = FLINT_MIN(GR_SERIES_PREC(ctx), err); \
    err = len; \
 \
    if (xlen <= 1 && xerr == GR_SERIES_ERR_EXACT) \
    { \
        len = FLINT_MIN(len, 1); \
        err = GR_SERIES_ERR_EXACT; \
    } \
 \
    GR_SERIES_ERROR(res) = err; \
    status |= gr_poly_ ## func ## _series(GR_SERIES_POLY(res), GR_SERIES_POLY(x), len, GR_SERIES_ELEM_CTX(ctx)); \
    return status; \
} \

UNARY_POLY_WRAPPER(exp)
UNARY_POLY_WRAPPER(log)
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

static int
arb_poly_is_finite(const arb_poly_t x)
{
    return _arb_vec_is_finite(x->coeffs, x->length);
}

static int
acb_poly_is_finite(const acb_poly_t x)
{
    return _acb_vec_is_finite(x->coeffs, x->length);
}

#define ARB_WRAPPER(func, arb_func, acb_func) \
int \
gr_series_ ## func(gr_series_t res, const gr_series_t x, gr_ctx_t ctx) \
{ \
    slong xlen, len, xerr, err; \
    slong prec; \
    int status = GR_SUCCESS; \
 \
    if (GR_SERIES_ELEM_CTX(ctx)->which_ring != GR_CTX_RR_ARB && GR_SERIES_ELEM_CTX(ctx)->which_ring != GR_CTX_CC_ACB) \
        return GR_UNABLE; \
 \
    xlen = GR_SERIES_POLY(x)->length; \
    xerr = GR_SERIES_ERROR(x); \
    err = xerr; \
 \
    len = FLINT_MIN(GR_SERIES_PREC(ctx), err); \
    err = len; \
 \
    if (xlen <= 1 && xerr == GR_SERIES_ERR_EXACT) \
    { \
        len = FLINT_MIN(len, 1); \
        err = GR_SERIES_ERR_EXACT; \
    } \
 \
    GR_SERIES_ERROR(res) = err; \
 \
    prec = _gr_ctx_get_real_prec(GR_SERIES_ELEM_CTX(ctx)); \
 \
    if (GR_SERIES_ELEM_CTX(ctx)->which_ring == GR_CTX_RR_ARB) \
    { \
        arb_func((arb_poly_struct *) GR_SERIES_POLY(res), (arb_poly_struct *) GR_SERIES_POLY(x), len, prec); \
        if (!arb_poly_is_finite((arb_poly_struct *) GR_SERIES_POLY(res))) \
            status = GR_UNABLE; \
    } \
    else \
    { \
        acb_func((acb_poly_struct *) GR_SERIES_POLY(res), (acb_poly_struct *) GR_SERIES_POLY(x), len, prec); \
        if (!acb_poly_is_finite((acb_poly_struct *) GR_SERIES_POLY(res))) \
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
gr_series_fresnel(gr_series_t res1, gr_series_t res2, const gr_series_t x, int normalized, gr_ctx_t ctx)
{
    slong xlen, len, xerr, err;
    slong prec;
    int status = GR_SUCCESS;

    if (GR_SERIES_ELEM_CTX(ctx)->which_ring != GR_CTX_RR_ARB && GR_SERIES_ELEM_CTX(ctx)->which_ring != GR_CTX_CC_ACB)
        return GR_UNABLE;

    xlen = GR_SERIES_POLY(x)->length;
    xerr = GR_SERIES_ERROR(x);
    err = xerr;

    len = FLINT_MIN(GR_SERIES_PREC(ctx), err);
    err = len;

    if (xlen <= 1 && xerr == GR_SERIES_ERR_EXACT)
    {
        len = FLINT_MIN(len, 1);
        err = GR_SERIES_ERR_EXACT;
    }

    if (res1 != NULL) res1->error = err;
    if (res2 != NULL) res2->error = err;

    prec = _gr_ctx_get_real_prec(GR_SERIES_ELEM_CTX(ctx));

    if (GR_SERIES_ELEM_CTX(ctx)->which_ring == GR_CTX_RR_ARB)
    {
        arb_hypgeom_fresnel_series(res1 ? (arb_poly_struct *) GR_SERIES_POLY(res1) : NULL,
                                   res2 ? (arb_poly_struct *) GR_SERIES_POLY(res2) : NULL,
                                    (arb_poly_struct *) GR_SERIES_POLY(x), normalized, len, prec);
        if (res1 && !arb_poly_is_finite((arb_poly_struct *) GR_SERIES_POLY(res1)))
            status = GR_UNABLE;
        if (res2 && !arb_poly_is_finite((arb_poly_struct *) GR_SERIES_POLY(res2)))
            status = GR_UNABLE;
    }
    else
    {
        acb_hypgeom_fresnel_series(res1 ? (acb_poly_struct *) GR_SERIES_POLY(res1) : NULL,
                                   res2 ? (acb_poly_struct *) GR_SERIES_POLY(res2) : NULL,
                                    (acb_poly_struct *) GR_SERIES_POLY(x), normalized, len, prec);
        if (res1 && !acb_poly_is_finite((acb_poly_struct *) GR_SERIES_POLY(res1)))
            status = GR_UNABLE;
        if (res2 && !acb_poly_is_finite((acb_poly_struct *) GR_SERIES_POLY(res2)))
            status = GR_UNABLE;
    }

    return status;
}

int
gr_series_fresnel_s(gr_series_t res, const gr_series_t x, int normalized, gr_ctx_t ctx)
{
    return gr_series_fresnel(res, NULL, x, normalized, ctx);
}

int
gr_series_fresnel_c(gr_series_t res, const gr_series_t x, int normalized, gr_ctx_t ctx)
{
    return gr_series_fresnel(NULL, res, x, normalized, ctx);
}

int
gr_series_airy(gr_series_t res1, gr_series_t res2, gr_series_t res3, gr_series_t res4, const gr_series_t x, gr_ctx_t ctx)
{
    slong xlen, len, xerr, err;
    slong prec;
    int status = GR_SUCCESS;

    if (GR_SERIES_ELEM_CTX(ctx)->which_ring != GR_CTX_RR_ARB && GR_SERIES_ELEM_CTX(ctx)->which_ring != GR_CTX_CC_ACB)
        return GR_UNABLE;

    xlen = GR_SERIES_POLY(x)->length;
    xerr = GR_SERIES_ERROR(x);
    err = xerr;

    len = FLINT_MIN(GR_SERIES_PREC(ctx), err);
    err = len;

    if (xlen <= 1 && xerr == GR_SERIES_ERR_EXACT)
    {
        len = FLINT_MIN(len, 1);
        err = GR_SERIES_ERR_EXACT;
    }

    if (res1 != NULL) res1->error = err;
    if (res2 != NULL) res2->error = err;
    if (res3 != NULL) res3->error = err;
    if (res4 != NULL) res4->error = err;

    prec = _gr_ctx_get_real_prec(GR_SERIES_ELEM_CTX(ctx));

    if (GR_SERIES_ELEM_CTX(ctx)->which_ring == GR_CTX_RR_ARB)
    {
        arb_hypgeom_airy_series(res1 ? (arb_poly_struct *) GR_SERIES_POLY(res1) : NULL,
                                res2 ? (arb_poly_struct *) GR_SERIES_POLY(res2) : NULL,
                                res3 ? (arb_poly_struct *) GR_SERIES_POLY(res3) : NULL,
                                res4 ? (arb_poly_struct *) GR_SERIES_POLY(res4) : NULL,
                                (arb_poly_struct *) GR_SERIES_POLY(x), len, prec);
        if (res1 && !arb_poly_is_finite((arb_poly_struct *) GR_SERIES_POLY(res1)))
            status = GR_UNABLE;
        if (res2 && !arb_poly_is_finite((arb_poly_struct *) GR_SERIES_POLY(res2)))
            status = GR_UNABLE;
        if (res3 && !arb_poly_is_finite((arb_poly_struct *) GR_SERIES_POLY(res3)))
            status = GR_UNABLE;
        if (res4 && !arb_poly_is_finite((arb_poly_struct *) GR_SERIES_POLY(res4)))
            status = GR_UNABLE;
    }
    else
    {
        acb_hypgeom_airy_series(res1 ? (acb_poly_struct *) GR_SERIES_POLY(res1) : NULL,
                                res2 ? (acb_poly_struct *) GR_SERIES_POLY(res2) : NULL,
                                res3 ? (acb_poly_struct *) GR_SERIES_POLY(res3) : NULL,
                                res4 ? (acb_poly_struct *) GR_SERIES_POLY(res4) : NULL,
                                (acb_poly_struct *) GR_SERIES_POLY(x), len, prec);
        if (res1 && !acb_poly_is_finite((acb_poly_struct *) GR_SERIES_POLY(res1)))
            status = GR_UNABLE;
        if (res2 && !acb_poly_is_finite((acb_poly_struct *) GR_SERIES_POLY(res2)))
            status = GR_UNABLE;
        if (res3 && !acb_poly_is_finite((acb_poly_struct *) GR_SERIES_POLY(res3)))
            status = GR_UNABLE;
        if (res4 && !acb_poly_is_finite((acb_poly_struct *) GR_SERIES_POLY(res4)))
            status = GR_UNABLE;
    }

    return status;
}

int
gr_series_airy_ai(gr_series_t res, const gr_series_t x, gr_ctx_t ctx)
{
    return gr_series_airy(res, NULL, NULL, NULL, x, ctx);
}

int
gr_series_airy_ai_prime(gr_series_t res, const gr_series_t x, gr_ctx_t ctx)
{
    return gr_series_airy(NULL, res, NULL, NULL, x, ctx);
}

int
gr_series_airy_bi(gr_series_t res, const gr_series_t x, gr_ctx_t ctx)
{
    return gr_series_airy(NULL, NULL, res, NULL, x, ctx);
}

int
gr_series_airy_bi_prime(gr_series_t res, const gr_series_t x, gr_ctx_t ctx)
{
    return gr_series_airy(NULL, NULL, NULL, res, x, ctx);
}


int
gr_series_log_integral(gr_series_t res, const gr_series_t x, int offset, gr_ctx_t ctx)
{
    slong xlen, len, xerr, err;
    slong prec;
    int status = GR_SUCCESS;

    if (GR_SERIES_ELEM_CTX(ctx)->which_ring != GR_CTX_RR_ARB && GR_SERIES_ELEM_CTX(ctx)->which_ring != GR_CTX_CC_ACB)
        return GR_UNABLE;

    xlen = GR_SERIES_POLY(x)->length;
    xerr = GR_SERIES_ERROR(x);
    err = xerr;

    len = FLINT_MIN(GR_SERIES_PREC(ctx), err);
    err = len;

    if (xlen <= 1 && xerr == GR_SERIES_ERR_EXACT)
    {
        len = FLINT_MIN(len, 1);
        err = GR_SERIES_ERR_EXACT;
    }

    GR_SERIES_ERROR(res) = err;

    prec = _gr_ctx_get_real_prec(GR_SERIES_ELEM_CTX(ctx));

    if (GR_SERIES_ELEM_CTX(ctx)->which_ring == GR_CTX_RR_ARB)
    {
        arb_hypgeom_li_series((arb_poly_struct *) GR_SERIES_POLY(res), (arb_poly_struct *) GR_SERIES_POLY(x), offset, len, prec);
        if (!arb_poly_is_finite((arb_poly_struct *) GR_SERIES_POLY(res)))
            status = GR_UNABLE;
    }
    else
    {
        acb_hypgeom_li_series((acb_poly_struct *) GR_SERIES_POLY(res), (acb_poly_struct *) GR_SERIES_POLY(x), offset, len, prec);
        if (!acb_poly_is_finite((acb_poly_struct *) GR_SERIES_POLY(res)))
            status = GR_UNABLE;
    }

    return status;
}

int
gr_series_gamma_upper(gr_series_t res, const gr_series_t s, const gr_series_t x, int regularized, gr_ctx_t ctx)
{
    slong xlen, len, xerr, err;
    slong prec;
    int status = GR_SUCCESS;

    if (GR_SERIES_ELEM_CTX(ctx)->which_ring != GR_CTX_RR_ARB && GR_SERIES_ELEM_CTX(ctx)->which_ring != GR_CTX_CC_ACB)
        return GR_UNABLE;

    xlen = GR_SERIES_POLY(x)->length;
    xerr = GR_SERIES_ERROR(x);
    err = xerr;

    len = FLINT_MIN(GR_SERIES_PREC(ctx), err);
    err = len;

    if (xlen <= 1 && xerr == GR_SERIES_ERR_EXACT)
    {
        len = FLINT_MIN(len, 1);
        err = GR_SERIES_ERR_EXACT;
    }

    /* we only handle constant s */
    if (len >= 2 && GR_SERIES_POLY(s)->length >= 2)
        return GR_UNABLE;

    GR_SERIES_ERROR(res) = err;

    prec = _gr_ctx_get_real_prec(GR_SERIES_ELEM_CTX(ctx));

    if (GR_SERIES_ELEM_CTX(ctx)->which_ring == GR_CTX_RR_ARB)
    {
        arb_t t;
        arb_init(t);
        arb_poly_get_coeff_arb(t, (arb_poly_struct *) GR_SERIES_POLY(s), 0);
        arb_hypgeom_gamma_upper_series((arb_poly_struct *) GR_SERIES_POLY(res), t, (arb_poly_struct *) GR_SERIES_POLY(x), regularized, len, prec);
        if (!arb_poly_is_finite((arb_poly_struct *) GR_SERIES_POLY(res)))
            status = GR_UNABLE;
        arb_clear(t);
    }
    else
    {
        acb_t t;
        acb_init(t);
        acb_poly_get_coeff_acb(t, (acb_poly_struct *) GR_SERIES_POLY(s), 0);
        acb_hypgeom_gamma_upper_series((acb_poly_struct *) GR_SERIES_POLY(res), t, (acb_poly_struct *) GR_SERIES_POLY(x), regularized, len, prec);
        if (!acb_poly_is_finite((acb_poly_struct *) GR_SERIES_POLY(res)))
            status = GR_UNABLE;
        acb_clear(t);
    }

    return status;
}

int
gr_series_gamma_lower(gr_series_t res, const gr_series_t s, const gr_series_t x, int regularized, gr_ctx_t ctx)
{
    slong xlen, len, xerr, err;
    slong prec;
    int status = GR_SUCCESS;

    if (GR_SERIES_ELEM_CTX(ctx)->which_ring != GR_CTX_RR_ARB && GR_SERIES_ELEM_CTX(ctx)->which_ring != GR_CTX_CC_ACB)
        return GR_UNABLE;

    xlen = GR_SERIES_POLY(x)->length;
    xerr = GR_SERIES_ERROR(x);
    err = xerr;

    len = FLINT_MIN(GR_SERIES_PREC(ctx), err);
    err = len;

    if (xlen <= 1 && xerr == GR_SERIES_ERR_EXACT)
    {
        len = FLINT_MIN(len, 1);
        err = GR_SERIES_ERR_EXACT;
    }

    /* we only handle constant s */
    if (len >= 2 && GR_SERIES_POLY(s)->length >= 2)
        return GR_UNABLE;

    GR_SERIES_ERROR(res) = err;

    prec = _gr_ctx_get_real_prec(GR_SERIES_ELEM_CTX(ctx));

    if (GR_SERIES_ELEM_CTX(ctx)->which_ring == GR_CTX_RR_ARB)
    {
        arb_t t;
        arb_init(t);
        arb_poly_get_coeff_arb(t, (arb_poly_struct *) GR_SERIES_POLY(s), 0);
        arb_hypgeom_gamma_lower_series((arb_poly_struct *) GR_SERIES_POLY(res), t, (arb_poly_struct *) GR_SERIES_POLY(x), regularized, len, prec);
        if (!arb_poly_is_finite((arb_poly_struct *) GR_SERIES_POLY(res)))
            status = GR_UNABLE;
        arb_clear(t);
    }
    else
    {
        acb_t t;
        acb_init(t);
        acb_poly_get_coeff_acb(t, (acb_poly_struct *) GR_SERIES_POLY(s), 0);
        acb_hypgeom_gamma_lower_series((acb_poly_struct *) GR_SERIES_POLY(res), t, (acb_poly_struct *) GR_SERIES_POLY(x), regularized, len, prec);
        if (!acb_poly_is_finite((acb_poly_struct *) GR_SERIES_POLY(res)))
            status = GR_UNABLE;
        acb_clear(t);
    }

    return status;
}

int
gr_series_beta_lower(gr_series_t res, const gr_series_t a, const gr_series_t b, const gr_series_t x, int regularized, gr_ctx_t ctx)
{
    slong xlen, len, xerr, err;
    slong prec;
    int status = GR_SUCCESS;

    if (GR_SERIES_ELEM_CTX(ctx)->which_ring != GR_CTX_RR_ARB && GR_SERIES_ELEM_CTX(ctx)->which_ring != GR_CTX_CC_ACB)
        return GR_UNABLE;

    xlen = GR_SERIES_POLY(x)->length;
    xerr = GR_SERIES_ERROR(x);
    err = xerr;

    len = FLINT_MIN(GR_SERIES_PREC(ctx), err);
    err = len;

    if (xlen <= 1 && xerr == GR_SERIES_ERR_EXACT)
    {
        len = FLINT_MIN(len, 1);
        err = GR_SERIES_ERR_EXACT;
    }

    /* we only handle constant a, b */
    if (len >= 2 && GR_SERIES_POLY(a)->length >= 2)
        return GR_UNABLE;
    if (len >= 2 && GR_SERIES_POLY(b)->length >= 2)
        return GR_UNABLE;

    GR_SERIES_ERROR(res) = err;

    prec = _gr_ctx_get_real_prec(GR_SERIES_ELEM_CTX(ctx));

    if (GR_SERIES_ELEM_CTX(ctx)->which_ring == GR_CTX_RR_ARB)
    {
        arb_t t, u;
        arb_init(t);
        arb_init(u);
        arb_poly_get_coeff_arb(t, (arb_poly_struct *) GR_SERIES_POLY(a), 0);
        arb_poly_get_coeff_arb(u, (arb_poly_struct *) GR_SERIES_POLY(b), 0);
        arb_hypgeom_beta_lower_series((arb_poly_struct *) GR_SERIES_POLY(res), t, u, (arb_poly_struct *) GR_SERIES_POLY(x), regularized, len, prec);
        if (!arb_poly_is_finite((arb_poly_struct *) GR_SERIES_POLY(res)))
            status = GR_UNABLE;
        arb_clear(t);
        arb_clear(u);
    }
    else
    {
        acb_t t, u;
        acb_init(t);
        acb_init(u);
        acb_poly_get_coeff_acb(t, (acb_poly_struct *) GR_SERIES_POLY(a), 0);
        acb_poly_get_coeff_acb(u, (acb_poly_struct *) GR_SERIES_POLY(b), 0);
        acb_hypgeom_beta_lower_series((acb_poly_struct *) GR_SERIES_POLY(res), t, u, (acb_poly_struct *) GR_SERIES_POLY(x), regularized, len, prec);
        if (!acb_poly_is_finite((acb_poly_struct *) GR_SERIES_POLY(res)))
            status = GR_UNABLE;
        acb_clear(t);
        acb_clear(u);
    }

    return status;
}

int
gr_series_polylog(gr_series_t res, const gr_series_t s, const gr_series_t z, gr_ctx_t ctx)
{
    slong xlen, len, xerr, err;
    slong prec;
    int status = GR_SUCCESS;

    if (GR_SERIES_ELEM_CTX(ctx)->which_ring != GR_CTX_CC_ACB)
        return GR_UNABLE;

    xlen = GR_SERIES_POLY(s)->length;
    xerr = GR_SERIES_ERROR(s);
    err = xerr;

    len = FLINT_MIN(GR_SERIES_PREC(ctx), err);
    err = len;

    if (xlen <= 1 && xerr == GR_SERIES_ERR_EXACT)
    {
        len = FLINT_MIN(len, 1);
        err = GR_SERIES_ERR_EXACT;
    }

    /* we only handle constant z */
    if (len >= 2 && GR_SERIES_POLY(z)->length >= 2)
        return GR_UNABLE;

    GR_SERIES_ERROR(res) = err;

    prec = _gr_ctx_get_real_prec(GR_SERIES_ELEM_CTX(ctx));

    {
        acb_t t;
        acb_init(t);
        acb_poly_get_coeff_acb(t, (acb_poly_struct *) GR_SERIES_POLY(z), 0);
        acb_poly_polylog_series((acb_poly_struct *) GR_SERIES_POLY(res), (acb_poly_struct *) GR_SERIES_POLY(s), t, len, prec);
        if (!acb_poly_is_finite((acb_poly_struct *) GR_SERIES_POLY(res)))
            status = GR_UNABLE;
        acb_clear(t);
    }

    return status;
}

int
gr_series_hurwitz_zeta(gr_series_t res, const gr_series_t s, const gr_series_t z, gr_ctx_t ctx)
{
    slong xlen, len, xerr, err;
    slong prec;
    int status = GR_SUCCESS;

    if (GR_SERIES_ELEM_CTX(ctx)->which_ring != GR_CTX_CC_ACB)
        return GR_UNABLE;

    xlen = GR_SERIES_POLY(s)->length;
    xerr = GR_SERIES_ERROR(s);
    err = xerr;

    len = FLINT_MIN(GR_SERIES_PREC(ctx), err);
    err = len;

    if (xlen <= 1 && xerr == GR_SERIES_ERR_EXACT)
    {
        len = FLINT_MIN(len, 1);
        err = GR_SERIES_ERR_EXACT;
    }

    /* we only handle constant z */
    if (len >= 2 && GR_SERIES_POLY(z)->length >= 2)
        return GR_UNABLE;

    GR_SERIES_ERROR(res) = err;

    prec = _gr_ctx_get_real_prec(GR_SERIES_ELEM_CTX(ctx));

    {
        acb_t t;
        acb_init(t);
        acb_poly_get_coeff_acb(t, (acb_poly_struct *) GR_SERIES_POLY(z), 0);
        acb_poly_zeta_series((acb_poly_struct *) GR_SERIES_POLY(res), (acb_poly_struct *) GR_SERIES_POLY(s), t, 0, len, prec);
        if (!acb_poly_is_finite((acb_poly_struct *) GR_SERIES_POLY(res)))
            status = GR_UNABLE;
        acb_clear(t);
    }

    return status;
}

int
gr_series_dirichlet_l(gr_series_t res, const dirichlet_group_t G, const dirichlet_char_t chi, const gr_series_t x, gr_ctx_t ctx)
{
    slong xlen, len, xerr, err;
    slong prec;
    int status = GR_SUCCESS;

    if (GR_SERIES_ELEM_CTX(ctx)->which_ring != GR_CTX_CC_ACB)
        return GR_UNABLE;

    xlen = GR_SERIES_POLY(x)->length;
    xerr = GR_SERIES_ERROR(x);
    err = xerr;

    len = FLINT_MIN(GR_SERIES_PREC(ctx), err);
    err = len;

    if (xlen <= 1 && xerr == GR_SERIES_ERR_EXACT)
    {
        len = FLINT_MIN(len, 1);
        err = GR_SERIES_ERR_EXACT;
    }

    GR_SERIES_ERROR(res) = err;

    prec = _gr_ctx_get_real_prec(GR_SERIES_ELEM_CTX(ctx));

    {
        acb_dirichlet_l_series((acb_poly_struct *) GR_SERIES_POLY(res), (acb_poly_struct *) GR_SERIES_POLY(x), G, chi, 0, len, prec);
        if (!acb_poly_is_finite((acb_poly_struct *) GR_SERIES_POLY(res)))
            status = GR_UNABLE;
    }

    return status;
}

int
gr_series_dirichlet_hardy_theta(gr_series_t res, const dirichlet_group_t G, const dirichlet_char_t chi, const gr_series_t x, gr_ctx_t ctx)
{
    slong xlen, len, xerr, err;
    slong prec;
    int status = GR_SUCCESS;

    if (GR_SERIES_ELEM_CTX(ctx)->which_ring != GR_CTX_CC_ACB)
        return GR_UNABLE;

    xlen = GR_SERIES_POLY(x)->length;
    xerr = GR_SERIES_ERROR(x);
    err = xerr;

    len = FLINT_MIN(GR_SERIES_PREC(ctx), err);
    err = len;

    if (xlen <= 1 && xerr == GR_SERIES_ERR_EXACT)
    {
        len = FLINT_MIN(len, 1);
        err = GR_SERIES_ERR_EXACT;
    }

    GR_SERIES_ERROR(res) = err;

    prec = _gr_ctx_get_real_prec(GR_SERIES_ELEM_CTX(ctx));

    {
        acb_dirichlet_hardy_theta_series((acb_poly_struct *) GR_SERIES_POLY(res), (acb_poly_struct *) GR_SERIES_POLY(x), G, chi, len, prec);
        if (!acb_poly_is_finite((acb_poly_struct *) GR_SERIES_POLY(res)))
            status = GR_UNABLE;
    }

    return status;
}

int
gr_series_dirichlet_hardy_z(gr_series_t res, const dirichlet_group_t G, const dirichlet_char_t chi, const gr_series_t x, gr_ctx_t ctx)
{
    slong xlen, len, xerr, err;
    slong prec;
    int status = GR_SUCCESS;

    if (GR_SERIES_ELEM_CTX(ctx)->which_ring != GR_CTX_CC_ACB)
        return GR_UNABLE;

    xlen = GR_SERIES_POLY(x)->length;
    xerr = GR_SERIES_ERROR(x);
    err = xerr;

    len = FLINT_MIN(GR_SERIES_PREC(ctx), err);
    err = len;

    if (xlen <= 1 && xerr == GR_SERIES_ERR_EXACT)
    {
        len = FLINT_MIN(len, 1);
        err = GR_SERIES_ERR_EXACT;
    }

    GR_SERIES_ERROR(res) = err;

    prec = _gr_ctx_get_real_prec(GR_SERIES_ELEM_CTX(ctx));

    {
        acb_dirichlet_hardy_z_series((acb_poly_struct *) GR_SERIES_POLY(res), (acb_poly_struct *) GR_SERIES_POLY(x), G, chi, len, prec);
        if (!acb_poly_is_finite((acb_poly_struct *) GR_SERIES_POLY(res)))
            status = GR_UNABLE;
    }

    return status;
}

int
gr_series_jacobi_theta(gr_series_t res1, gr_series_t res2, gr_series_t res3, gr_series_t res4, const gr_series_t x, const gr_series_t tau, gr_ctx_t ctx)
{
    slong xlen, len, xerr, err;
    slong prec;
    int status = GR_SUCCESS;

    if (GR_SERIES_ELEM_CTX(ctx)->which_ring != GR_CTX_CC_ACB)
        return GR_UNABLE;

    xlen = GR_SERIES_POLY(x)->length;
    xerr = GR_SERIES_ERROR(x);
    err = xerr;

    len = FLINT_MIN(GR_SERIES_PREC(ctx), err);
    err = len;

    if (xlen <= 1 && xerr == GR_SERIES_ERR_EXACT)
    {
        len = FLINT_MIN(len, 1);
        err = GR_SERIES_ERR_EXACT;
    }

    /* we only handle constant z */
    if (len >= 2 && GR_SERIES_POLY(tau)->length >= 2)
        return GR_UNABLE;

    if (res1 != NULL) res1->error = err;
    if (res2 != NULL) res2->error = err;
    if (res3 != NULL) res3->error = err;
    if (res4 != NULL) res4->error = err;

    prec = _gr_ctx_get_real_prec(GR_SERIES_ELEM_CTX(ctx));

    {
        acb_t t;
        acb_init(t);
        acb_poly_get_coeff_acb(t, (acb_poly_struct *) GR_SERIES_POLY(tau), 0);
        acb_modular_theta_series(res1 ? (acb_poly_struct *) GR_SERIES_POLY(res1) : NULL,
                                res2 ? (acb_poly_struct *) GR_SERIES_POLY(res2) : NULL,
                                res3 ? (acb_poly_struct *) GR_SERIES_POLY(res3) : NULL,
                                res4 ? (acb_poly_struct *) GR_SERIES_POLY(res4) : NULL,
                                (acb_poly_struct *) GR_SERIES_POLY(x), t, len, prec);
        if (res1 && !acb_poly_is_finite((acb_poly_struct *) GR_SERIES_POLY(res1)))
            status = GR_UNABLE;
        if (res2 && !acb_poly_is_finite((acb_poly_struct *) GR_SERIES_POLY(res2)))
            status = GR_UNABLE;
        if (res3 && !acb_poly_is_finite((acb_poly_struct *) GR_SERIES_POLY(res3)))
            status = GR_UNABLE;
        if (res4 && !acb_poly_is_finite((acb_poly_struct *) GR_SERIES_POLY(res4)))
            status = GR_UNABLE;
        acb_clear(t);
    }

    return status;
}

int
gr_series_jacobi_theta_1(gr_series_t res, const gr_series_t x, const gr_series_t tau, gr_ctx_t ctx)
{
    return gr_series_jacobi_theta(res, NULL, NULL, NULL, x, tau, ctx);
}

int
gr_series_jacobi_theta_2(gr_series_t res, const gr_series_t x, const gr_series_t tau, gr_ctx_t ctx)
{
    return gr_series_jacobi_theta(NULL, res, NULL, NULL, x, tau, ctx);
}

int
gr_series_jacobi_theta_3(gr_series_t res, const gr_series_t x, const gr_series_t tau, gr_ctx_t ctx)
{
    return gr_series_jacobi_theta(NULL, NULL, res, NULL, x, tau, ctx);
}

int
gr_series_jacobi_theta_4(gr_series_t res, const gr_series_t x, const gr_series_t tau, gr_ctx_t ctx)
{
    return gr_series_jacobi_theta(NULL, NULL, NULL, res, x, tau, ctx);
}


int
gr_series_agm1(gr_series_t res, const gr_series_t x, gr_ctx_t ctx)
{
    slong xlen, len, xerr, err;
    slong prec;
    int status = GR_SUCCESS;

    if (GR_SERIES_ELEM_CTX(ctx)->which_ring != GR_CTX_CC_ACB)
        return GR_UNABLE;

    xlen = GR_SERIES_POLY(x)->length;
    xerr = GR_SERIES_ERROR(x);
    err = xerr;

    len = FLINT_MIN(GR_SERIES_PREC(ctx), err);
    err = len;

    if (xlen <= 1 && xerr == GR_SERIES_ERR_EXACT)
    {
        len = FLINT_MIN(len, 1);
        err = GR_SERIES_ERR_EXACT;
    }

    GR_SERIES_ERROR(res) = err;

    prec = _gr_ctx_get_real_prec(GR_SERIES_ELEM_CTX(ctx));

    {
        acb_poly_agm1_series((acb_poly_struct *) GR_SERIES_POLY(res), (acb_poly_struct *) GR_SERIES_POLY(x), len, prec);
        if (!acb_poly_is_finite((acb_poly_struct *) GR_SERIES_POLY(res)))
            status = GR_UNABLE;
    }

    return status;
}

int
gr_series_elliptic_k(gr_series_t res, const gr_series_t x, gr_ctx_t ctx)
{
    slong xlen, len, xerr, err;
    slong prec;
    int status = GR_SUCCESS;

    if (GR_SERIES_ELEM_CTX(ctx)->which_ring != GR_CTX_CC_ACB)
        return GR_UNABLE;

    xlen = GR_SERIES_POLY(x)->length;
    xerr = GR_SERIES_ERROR(x);
    err = xerr;

    len = FLINT_MIN(GR_SERIES_PREC(ctx), err);
    err = len;

    if (xlen <= 1 && xerr == GR_SERIES_ERR_EXACT)
    {
        len = FLINT_MIN(len, 1);
        err = GR_SERIES_ERR_EXACT;
    }

    GR_SERIES_ERROR(res) = err;

    prec = _gr_ctx_get_real_prec(GR_SERIES_ELEM_CTX(ctx));

    {
        acb_poly_elliptic_k_series((acb_poly_struct *) GR_SERIES_POLY(res), (acb_poly_struct *) GR_SERIES_POLY(x), len, prec);
        if (!acb_poly_is_finite((acb_poly_struct *) GR_SERIES_POLY(res)))
            status = GR_UNABLE;
    }

    return status;
}

int
gr_series_weierstrass_p(gr_series_t res, const gr_series_t x, const gr_series_t tau, gr_ctx_t ctx)
{
    slong xlen, len, xerr, err;
    slong prec;
    int status = GR_SUCCESS;

    if (GR_SERIES_ELEM_CTX(ctx)->which_ring != GR_CTX_CC_ACB)
        return GR_UNABLE;

    xlen = GR_SERIES_POLY(x)->length;
    xerr = GR_SERIES_ERROR(x);
    err = xerr;

    len = FLINT_MIN(GR_SERIES_PREC(ctx), err);
    err = len;

    if (xlen <= 1 && xerr == GR_SERIES_ERR_EXACT)
    {
        len = FLINT_MIN(len, 1);
        err = GR_SERIES_ERR_EXACT;
    }

    /* we only handle constant tau */
    if (len >= 2 && GR_SERIES_POLY(tau)->length >= 2)
        return GR_UNABLE;

    GR_SERIES_ERROR(res) = err;

    prec = _gr_ctx_get_real_prec(GR_SERIES_ELEM_CTX(ctx));

    {
        acb_t t;
        acb_init(t);
        acb_poly_get_coeff_acb(t, (acb_poly_struct *) GR_SERIES_POLY(tau), 0);
        acb_elliptic_p_series((acb_poly_struct *) GR_SERIES_POLY(res), (acb_poly_struct *) GR_SERIES_POLY(x), t, len, prec);
        if (!acb_poly_is_finite((acb_poly_struct *) GR_SERIES_POLY(res)))
            status = GR_UNABLE;
        acb_clear(t);
    }

    return status;
}

int
gr_series_hypgeom_pfq(gr_series_t res, const gr_series_vec_t a, const gr_series_vec_t b, const gr_series_t x, int regularized, gr_ctx_t ctx)
{
    slong len, err, i;
    slong prec;
    slong p, q, pp, qq;
    int have_one;
    int status = GR_SUCCESS;
    acb_poly_struct *aa, *bb;

    if (GR_SERIES_ELEM_CTX(ctx)->which_ring != GR_CTX_CC_ACB)
        return GR_UNABLE;

    p = a->length;
    q = b->length;

    err = GR_SERIES_ERROR(x);
    for (i = 0; i < p; i++)
        err = FLINT_MIN(err, a->entries[i].error);
    for (i = 0; i < q; i++)
        err = FLINT_MIN(err, b->entries[i].error);

    len = FLINT_MIN(GR_SERIES_PREC(ctx), err);
    err = len;

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

    prec = _gr_ctx_get_real_prec(GR_SERIES_ELEM_CTX(ctx));

    GR_SERIES_ERROR(res) = err;

    acb_hypgeom_pfq_series_direct((acb_poly_struct *) GR_SERIES_POLY(res), aa, pp, bb, qq, (const acb_poly_struct *) GR_SERIES_POLY(x), regularized, -1, len, prec);

    if (!acb_poly_is_finite((acb_poly_struct *) GR_SERIES_POLY(res)))
        status = GR_UNABLE;

    if (!have_one)
        acb_poly_clear(bb + qq - 1);

    GR_TMP_FREE(aa, sizeof(acb_poly_struct) * (p + q + 1));

    return status;
}

void gr_series_ctx_clear(gr_ctx_t ctx)
{
    if (GR_SERIES_CTX(ctx)->var != default_var)
    {
        flint_free(GR_SERIES_CTX(ctx)->var);
    }

}

int gr_series_ctx_write(gr_stream_t out, gr_ctx_t ctx)
{
    gr_stream_write(out, "Power series over ");
    gr_ctx_write(out, GR_SERIES_ELEM_CTX(ctx));
    gr_stream_write(out, " with precision ");
    gr_stream_write(out, "O(");
    gr_stream_write(out, GR_SERIES_CTX(ctx)->var);
    gr_stream_write(out, "^");
    gr_stream_write_si(out, GR_SERIES_PREC(ctx));
    gr_stream_write(out, ")");
    return GR_SUCCESS;
}

truth_t
gr_series_ctx_is_ring(gr_ctx_t ctx)
{
    return gr_ctx_is_ring(GR_SERIES_ELEM_CTX(ctx));
}

truth_t
gr_series_ctx_is_commutative_ring(gr_ctx_t ctx)
{
    return gr_ctx_is_commutative_ring(GR_SERIES_ELEM_CTX(ctx));
}

truth_t
gr_series_ctx_is_integral_domain(gr_ctx_t ctx)
{
    return gr_ctx_is_integral_domain(GR_SERIES_ELEM_CTX(ctx));
}

truth_t
gr_series_ctx_is_rational_vector_space(gr_ctx_t ctx)
{
    return gr_ctx_is_rational_vector_space(GR_SERIES_ELEM_CTX(ctx));
}

truth_t
gr_series_ctx_is_real_vector_space(gr_ctx_t ctx)
{
    return gr_ctx_is_real_vector_space(GR_SERIES_ELEM_CTX(ctx));
}

truth_t
gr_series_ctx_is_complex_vector_space(gr_ctx_t ctx)
{
    return gr_ctx_is_complex_vector_space(GR_SERIES_ELEM_CTX(ctx));
}

int gr_series_ctx_set_gen_name(gr_ctx_t ctx, const char * s)
{
    slong len;
    len = strlen(s);

    if (GR_SERIES_CTX(ctx)->var == default_var)
        GR_SERIES_CTX(ctx)->var = NULL;

    GR_SERIES_CTX(ctx)->var = flint_realloc(GR_SERIES_CTX(ctx)->var, len + 1);
    memcpy(GR_SERIES_CTX(ctx)->var, s, len + 1);
    return GR_SUCCESS;
}

int gr_series_ctx_set_gen_names(gr_ctx_t ctx, const char ** s)
{
    return gr_series_ctx_set_gen_name(ctx, s[0]);
}


int
gr_series_gens_recursive(gr_vec_t vec, gr_ctx_t ctx)
{
    int status;
    gr_vec_t vec1;
    slong i, n;

    /* Get generators of the element ring */
    gr_vec_init(vec1, 0, GR_SERIES_ELEM_CTX(ctx));
    status = gr_gens_recursive(vec1, GR_SERIES_ELEM_CTX(ctx));
    n = vec1->length;

    gr_vec_set_length(vec, n + 1, ctx);

    /* Promote to series */
    for (i = 0; i < n; i++)
        status |= gr_series_set_scalar(gr_vec_entry_ptr(vec, i, ctx),
                gr_vec_entry_srcptr(vec1, i, GR_SERIES_ELEM_CTX(ctx)), ctx);

    status |= gr_series_gen(gr_vec_entry_ptr(vec, n, ctx), ctx);

    gr_vec_clear(vec1, GR_SERIES_ELEM_CTX(ctx));

    return status;
}

int
gr_series_set_other(gr_series_t res, gr_srcptr x, gr_ctx_t x_ctx, gr_ctx_t ctx)
{
    if (x_ctx == ctx)
    {
        return gr_series_set(res, x, ctx);
    }
    else if (x_ctx == GR_SERIES_ELEM_CTX(ctx))
    {
        return gr_series_set_scalar(res, x, ctx);
    }
    else if (x_ctx->which_ring == GR_CTX_GR_SERIES && !strcmp(GR_SERIES_CTX(x_ctx)->var, GR_SERIES_CTX(ctx)->var))
    {
        int status = GR_SUCCESS;
        const gr_series_struct * xser = (const gr_series_struct *) x;

        /* todo: check compatibility? */
        status |= gr_poly_set_gr_poly_other(GR_SERIES_POLY(res), GR_SERIES_POLY(xser), GR_SERIES_ELEM_CTX(x_ctx), GR_SERIES_ELEM_CTX(ctx));
        GR_SERIES_ERROR(res) = GR_SERIES_ERROR(xser);
        /* possible truncation */
        status |= gr_series_set(res, res, ctx);
        return status;
    }
    else if (x_ctx->which_ring == GR_CTX_GR_POLY && !strcmp(POLYNOMIAL_CTX(x_ctx)->var, GR_SERIES_CTX(ctx)->var))
    {
        int status = GR_SUCCESS;

        /* todo: check compatibility? */
        status |= gr_poly_set_gr_poly_other(GR_SERIES_POLY(res), x, POLYNOMIAL_ELEM_CTX(x_ctx), GR_SERIES_ELEM_CTX(ctx));
        GR_SERIES_ERROR(res) = GR_SERIES_ERR_EXACT;
        /* possible truncation */
        status |= gr_series_set(res, res, ctx);
        return status;
    }
    else
    {
        int status = GR_SUCCESS;

        gr_poly_fit_length(GR_SERIES_POLY(res), 1, GR_SERIES_ELEM_CTX(ctx));
        status = gr_set_other(GR_SERIES_POLY(res)->coeffs, x, x_ctx, GR_SERIES_ELEM_CTX(ctx));
        if (status == GR_SUCCESS)
        {
            _gr_poly_set_length(GR_SERIES_POLY(res), 1, GR_SERIES_ELEM_CTX(ctx));
            _gr_poly_normalise(GR_SERIES_POLY(res), GR_SERIES_ELEM_CTX(ctx));
        }
        else
            _gr_poly_set_length(GR_SERIES_POLY(res), 0, GR_SERIES_ELEM_CTX(ctx));

        GR_SERIES_ERROR(res) = GR_SERIES_ERR_EXACT;
        /* possible truncation */
        status |= gr_series_set(res, res, ctx);

        return status;
    }

    return GR_UNABLE;
}



int _gr_series_methods_initialized = 0;

gr_static_method_table _gr_series_methods;

gr_method_tab_input _gr_series_methods_input[] =
{
    {GR_METHOD_CTX_CLEAR,   (gr_funcptr) gr_series_ctx_clear},
    {GR_METHOD_CTX_WRITE,   (gr_funcptr) gr_series_ctx_write},
    {GR_METHOD_CTX_SET_GEN_NAME, (gr_funcptr) gr_series_ctx_set_gen_name},
    {GR_METHOD_CTX_SET_GEN_NAMES, (gr_funcptr) gr_series_ctx_set_gen_names},
    {GR_METHOD_CTX_IS_RING, (gr_funcptr) gr_series_ctx_is_ring},
    {GR_METHOD_CTX_IS_COMMUTATIVE_RING, (gr_funcptr) gr_series_ctx_is_commutative_ring},
    {GR_METHOD_CTX_IS_INTEGRAL_DOMAIN, (gr_funcptr) gr_series_ctx_is_integral_domain},
    {GR_METHOD_CTX_IS_RATIONAL_VECTOR_SPACE, (gr_funcptr) gr_series_ctx_is_rational_vector_space},
    {GR_METHOD_CTX_IS_REAL_VECTOR_SPACE, (gr_funcptr) gr_series_ctx_is_real_vector_space},
    {GR_METHOD_CTX_IS_COMPLEX_VECTOR_SPACE, (gr_funcptr) gr_series_ctx_is_complex_vector_space},
    {GR_METHOD_INIT,        (gr_funcptr) gr_series_init},
    {GR_METHOD_CLEAR,       (gr_funcptr) gr_series_clear},
    {GR_METHOD_SWAP,        (gr_funcptr) gr_series_swap},
    {GR_METHOD_SET_SHALLOW, (gr_funcptr) gr_series_set_shallow},
    {GR_METHOD_RANDTEST,    (gr_funcptr) gr_series_randtest},
    {GR_METHOD_WRITE,       (gr_funcptr) gr_series_write},
    {GR_METHOD_ZERO,        (gr_funcptr) gr_series_zero},
    {GR_METHOD_ONE,         (gr_funcptr) gr_series_one},
    {GR_METHOD_IS_ZERO,     (gr_funcptr) gr_series_is_zero},
    {GR_METHOD_IS_ONE,      (gr_funcptr) gr_series_is_one},
    {GR_METHOD_EQUAL,       (gr_funcptr) gr_series_equal},
    {GR_METHOD_GEN,         (gr_funcptr) gr_series_gen},
    {GR_METHOD_GENS,        (gr_funcptr) gr_generic_gens_single},
    {GR_METHOD_GENS_RECURSIVE,  (gr_funcptr) gr_series_gens_recursive},
    {GR_METHOD_SET,         (gr_funcptr) gr_series_set},
    {GR_METHOD_SET_UI,      (gr_funcptr) gr_series_set_ui},
    {GR_METHOD_SET_SI,      (gr_funcptr) gr_series_set_si},
    {GR_METHOD_SET_FMPZ,    (gr_funcptr) gr_series_set_fmpz},
    {GR_METHOD_SET_FMPQ,    (gr_funcptr) gr_series_set_fmpq},
    {GR_METHOD_SET_OTHER,   (gr_funcptr) gr_series_set_other},
    {GR_METHOD_SET_STR,     (gr_funcptr) gr_generic_set_str_balance_additions},
    {GR_METHOD_NEG,         (gr_funcptr) gr_series_neg},
    {GR_METHOD_ADD,         (gr_funcptr) gr_series_add},
    {GR_METHOD_SUB,         (gr_funcptr) gr_series_sub},
    {GR_METHOD_MUL,         (gr_funcptr) gr_series_mul},
    {GR_METHOD_INV,         (gr_funcptr) gr_series_inv},
    {GR_METHOD_DIV,         (gr_funcptr) gr_series_div},
    {GR_METHOD_DIVEXACT,    (gr_funcptr) gr_series_divexact},
    {GR_METHOD_SQRT,        (gr_funcptr) gr_series_sqrt},
    {GR_METHOD_RSQRT,       (gr_funcptr) gr_series_rsqrt},
    {GR_METHOD_EXP,         (gr_funcptr) gr_series_exp},
    {GR_METHOD_LOG,         (gr_funcptr) gr_series_log},
    {GR_METHOD_TAN,         (gr_funcptr) gr_series_tan},
    {GR_METHOD_ASIN,        (gr_funcptr) gr_series_asin},
    {GR_METHOD_ACOS,        (gr_funcptr) gr_series_acos},
    {GR_METHOD_ATAN,        (gr_funcptr) gr_series_atan},
    {GR_METHOD_ASINH,       (gr_funcptr) gr_series_asinh},
    {GR_METHOD_ACOSH,       (gr_funcptr) gr_series_acosh},
    {GR_METHOD_ATANH,       (gr_funcptr) gr_series_atanh},
    {GR_METHOD_GAMMA,       (gr_funcptr) gr_series_gamma},
    {GR_METHOD_RGAMMA,      (gr_funcptr) gr_series_rgamma},
    {GR_METHOD_LGAMMA,      (gr_funcptr) gr_series_lgamma},
    {GR_METHOD_DIGAMMA,     (gr_funcptr) gr_series_digamma},
    {GR_METHOD_ERF,         (gr_funcptr) gr_series_erf},
    {GR_METHOD_ERFC,        (gr_funcptr) gr_series_erfc},
    {GR_METHOD_ERFI,        (gr_funcptr) gr_series_erfi},
    {GR_METHOD_FRESNEL,     (gr_funcptr) gr_series_fresnel},
    {GR_METHOD_FRESNEL_S,   (gr_funcptr) gr_series_fresnel_s},
    {GR_METHOD_FRESNEL_C,   (gr_funcptr) gr_series_fresnel_c},
    {GR_METHOD_AIRY,        (gr_funcptr) gr_series_airy},
    {GR_METHOD_AIRY_AI,        (gr_funcptr) gr_series_airy_ai},
    {GR_METHOD_AIRY_AI_PRIME,  (gr_funcptr) gr_series_airy_ai_prime},
    {GR_METHOD_AIRY_BI,        (gr_funcptr) gr_series_airy_bi},
    {GR_METHOD_AIRY_BI_PRIME,  (gr_funcptr) gr_series_airy_bi_prime},
    {GR_METHOD_EXP_INTEGRAL_EI,       (gr_funcptr) gr_series_exp_integral_ei},
    {GR_METHOD_COS_INTEGRAL,          (gr_funcptr) gr_series_cos_integral},
    {GR_METHOD_COSH_INTEGRAL,         (gr_funcptr) gr_series_cosh_integral},
    {GR_METHOD_SIN_INTEGRAL,          (gr_funcptr) gr_series_sin_integral},
    {GR_METHOD_SINH_INTEGRAL,         (gr_funcptr) gr_series_sinh_integral},
    {GR_METHOD_LOG_INTEGRAL,          (gr_funcptr) gr_series_log_integral},
    {GR_METHOD_GAMMA_UPPER,           (gr_funcptr) gr_series_gamma_upper},
    {GR_METHOD_GAMMA_LOWER,           (gr_funcptr) gr_series_gamma_lower},
    {GR_METHOD_BETA_LOWER,            (gr_funcptr) gr_series_beta_lower},
    {GR_METHOD_HYPGEOM_PFQ,           (gr_funcptr) gr_series_hypgeom_pfq},
    {GR_METHOD_HURWITZ_ZETA,          (gr_funcptr) gr_series_hurwitz_zeta},
    {GR_METHOD_POLYLOG,               (gr_funcptr) gr_series_polylog},
    {GR_METHOD_DIRICHLET_L,           (gr_funcptr) gr_series_dirichlet_l},
    {GR_METHOD_DIRICHLET_HARDY_Z,     (gr_funcptr) gr_series_dirichlet_hardy_z},
    {GR_METHOD_DIRICHLET_HARDY_THETA, (gr_funcptr) gr_series_dirichlet_hardy_theta},
    {GR_METHOD_JACOBI_THETA,          (gr_funcptr) gr_series_jacobi_theta},
    {GR_METHOD_JACOBI_THETA_1,          (gr_funcptr) gr_series_jacobi_theta_1},
    {GR_METHOD_JACOBI_THETA_2,          (gr_funcptr) gr_series_jacobi_theta_2},
    {GR_METHOD_JACOBI_THETA_3,          (gr_funcptr) gr_series_jacobi_theta_3},
    {GR_METHOD_JACOBI_THETA_4,          (gr_funcptr) gr_series_jacobi_theta_4},
    {GR_METHOD_AGM1,                   (gr_funcptr) gr_series_agm1},
    {GR_METHOD_ELLIPTIC_K,             (gr_funcptr) gr_series_elliptic_k},
    {GR_METHOD_WEIERSTRASS_P,          (gr_funcptr) gr_series_weierstrass_p},

    {0,                     (gr_funcptr) NULL},
};

void
gr_series_ctx_init(gr_ctx_t ctx, gr_ctx_t base_ring, slong prec)
{
    ctx->which_ring = GR_CTX_GR_SERIES;
    ctx->sizeof_elem = sizeof(gr_series_struct);
    ctx->size_limit = WORD_MAX;

    prec = FLINT_MIN(prec, GR_SERIES_ERR_MAX);
    prec = FLINT_MAX(prec, 0);

    GR_SERIES_CTX(ctx)->base_ring = (gr_ctx_struct *) base_ring;
    GR_SERIES_CTX(ctx)->prec = prec;
    GR_SERIES_CTX(ctx)->var = (char *) default_var;

    ctx->methods = _gr_series_methods;

    if (!_gr_series_methods_initialized)
    {
        gr_method_tab_init(_gr_series_methods, _gr_series_methods_input);
        _gr_series_methods_initialized = 1;
    }
}

/* compatibility */
void
gr_ctx_init_gr_series(gr_ctx_t ctx, gr_ctx_t base_ring, slong n)
{
    gr_series_ctx_init(ctx, base_ring, n);
}

