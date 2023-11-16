/*
    Copyright (C) 2023 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "gr_special.h"
#include "gr_vec.h"
#include "gr_poly.h"

static int
_gr_poly_inv_trig_series(gr_ptr res, gr_srcptr f, slong flen, slong len, gr_ctx_t ctx, int function)
{
    gr_ptr c;
    int status = GR_SUCCESS;
    slong sz = ctx->sizeof_elem;

    flen = FLINT_MIN(flen, len);

    if (flen == 0)
        return GR_UNABLE;

    GR_TMP_INIT(c, ctx);

    status |= (((gr_method_unary_op *) ctx->methods)[function])(c, f, ctx);

    if (status == GR_SUCCESS)
    {
        if (flen == 1)
        {
            status |= _gr_vec_zero(GR_ENTRY(res, 1, sz), len - 1, ctx);
        }
        else
        {
            gr_ptr t, u;
            slong ulen;

            ulen = FLINT_MIN(len, 2 * flen - 1);

            GR_TMP_INIT_VEC(t, len + ulen, ctx);
            u = GR_ENTRY(t, len, sz);

            /* asin(f(x)) = integral(f'(x)/sqrt(1-f(x)^2)) */
            /* acos(f(x)) = integral(-f'(x)/sqrt(1-f(x)^2)) */
            /* asinh(f(x)) = integral(f'(x)/sqrt(f(x)^2+1)) */
            /* acosh(f(x)) = integral(f'(x)/sqrt(f(x)^2-1)) */
            /* atan(f(x)) = integral(f'(x)/(1+f(x)^2)) */
            /* atanh(f(x)) = integral(f'(x)/(1-f(x)^2)) */
            status |= _gr_poly_mullow(u, f, flen, f, flen, ulen, ctx);

            if (function == GR_METHOD_ASINH || function == GR_METHOD_ATAN)
                status |= gr_add_ui(u, u, 1, ctx);
            else
                status |= gr_sub_ui(u, u, 1, ctx);

            if (function == GR_METHOD_ASIN || function == GR_METHOD_ACOS || function == GR_METHOD_ATANH)
                status |= _gr_vec_neg(u, u, ulen, ctx);

            if (function == GR_METHOD_ATAN || function == GR_METHOD_ATANH)
            {
                status |= _gr_poly_derivative(t, f, flen, ctx);
                status |= _gr_poly_div_series(res, t, flen - 1, u, ulen, len, ctx);
            }
            else
            {
                status |= _gr_poly_rsqrt_series(t, u, ulen, len, ctx);
                status |= _gr_poly_derivative(u, f, flen, ctx);
                status |= _gr_poly_mullow(res, t, len, u, flen - 1, len, ctx);
            }

            status |= _gr_poly_integral(res, res, len, ctx);
            if (function == GR_METHOD_ACOS)
                status |= _gr_vec_neg(res, res, len, ctx);

            GR_TMP_CLEAR_VEC(t, len + ulen, ctx);
        }

        gr_swap(res, c, ctx);
    }

    GR_TMP_CLEAR(c, ctx);

    return status;
}

int
_gr_poly_asin_series(gr_ptr res, gr_srcptr f, slong flen, slong len, gr_ctx_t ctx)
{
    return _gr_poly_inv_trig_series(res, f, flen, len, ctx, GR_METHOD_ASIN);
}

int
_gr_poly_acos_series(gr_ptr res, gr_srcptr f, slong flen, slong len, gr_ctx_t ctx)
{
    return _gr_poly_inv_trig_series(res, f, flen, len, ctx, GR_METHOD_ACOS);
}

int
_gr_poly_asinh_series(gr_ptr res, gr_srcptr f, slong flen, slong len, gr_ctx_t ctx)
{
    return _gr_poly_inv_trig_series(res, f, flen, len, ctx, GR_METHOD_ASINH);
}

int
_gr_poly_acosh_series(gr_ptr res, gr_srcptr f, slong flen, slong len, gr_ctx_t ctx)
{
    return _gr_poly_inv_trig_series(res, f, flen, len, ctx, GR_METHOD_ACOSH);
}

int
_gr_poly_atan_series(gr_ptr res, gr_srcptr f, slong flen, slong len, gr_ctx_t ctx)
{
    return _gr_poly_inv_trig_series(res, f, flen, len, ctx, GR_METHOD_ATAN);
}

int
_gr_poly_atanh_series(gr_ptr res, gr_srcptr f, slong flen, slong len, gr_ctx_t ctx)
{
    return _gr_poly_inv_trig_series(res, f, flen, len, ctx, GR_METHOD_ATANH);
}

int
gr_poly_asin_series(gr_poly_t res, const gr_poly_t f, slong len, gr_ctx_t ctx)
{
    slong flen = f->length;
    int status = GR_SUCCESS;

    if (flen == 0 || len == 0)
        return gr_poly_zero(res, ctx);

    if (flen == 1)
        len = 1;

    gr_poly_fit_length(res, len, ctx);
    status |= _gr_poly_asin_series(res->coeffs, f->coeffs, flen, len, ctx);
    _gr_poly_set_length(res, len, ctx);
    _gr_poly_normalise(res, ctx);
    return status;
}

int
gr_poly_asinh_series(gr_poly_t res, const gr_poly_t f, slong len, gr_ctx_t ctx)
{
    slong flen = f->length;
    int status = GR_SUCCESS;

    if (flen == 0 || len == 0)
        return gr_poly_zero(res, ctx);

    if (flen == 1)
        len = 1;

    gr_poly_fit_length(res, len, ctx);
    status |= _gr_poly_asinh_series(res->coeffs, f->coeffs, flen, len, ctx);
    _gr_poly_set_length(res, len, ctx);
    _gr_poly_normalise(res, ctx);
    return status;
}

int
gr_poly_acos_series(gr_poly_t res, const gr_poly_t f, slong len, gr_ctx_t ctx)
{
    slong flen = f->length;
    int status = GR_SUCCESS;

    if (len == 0)
        return gr_poly_zero(res, ctx);

    if (flen == 1)
        len = 1;

    gr_poly_fit_length(res, len, ctx);

    if (flen == 0)
    {
        status |= gr_zero(res->coeffs, ctx);
        status |= _gr_poly_acos_series(res->coeffs, res->coeffs, 1, len, ctx);
    }
    else
    {
        status |= _gr_poly_acos_series(res->coeffs, f->coeffs, flen, len, ctx);
    }

    _gr_poly_set_length(res, len, ctx);
    _gr_poly_normalise(res, ctx);
    return status;
}

int
gr_poly_acosh_series(gr_poly_t res, const gr_poly_t f, slong len, gr_ctx_t ctx)
{
    slong flen = f->length;
    int status = GR_SUCCESS;

    if (len == 0)
        return gr_poly_zero(res, ctx);

    if (flen == 1)
        len = 1;

    gr_poly_fit_length(res, len, ctx);

    if (flen == 0)
    {
        status |= gr_zero(res->coeffs, ctx);
        status |= _gr_poly_acosh_series(res->coeffs, res->coeffs, 1, len, ctx);
    }
    else
    {
        status |= _gr_poly_acosh_series(res->coeffs, f->coeffs, flen, len, ctx);
    }

    _gr_poly_set_length(res, len, ctx);
    _gr_poly_normalise(res, ctx);
    return status;
}

int
gr_poly_atan_series(gr_poly_t res, const gr_poly_t f, slong len, gr_ctx_t ctx)
{
    slong flen = f->length;
    int status = GR_SUCCESS;

    if (flen == 0 || len == 0)
        return gr_poly_zero(res, ctx);

    if (flen == 1)
        len = 1;

    gr_poly_fit_length(res, len, ctx);
    status |= _gr_poly_atan_series(res->coeffs, f->coeffs, flen, len, ctx);
    _gr_poly_set_length(res, len, ctx);
    _gr_poly_normalise(res, ctx);
    return status;
}

int
gr_poly_atanh_series(gr_poly_t res, const gr_poly_t f, slong len, gr_ctx_t ctx)
{
    slong flen = f->length;
    int status = GR_SUCCESS;

    if (flen == 0 || len == 0)
        return gr_poly_zero(res, ctx);

    if (flen == 1)
        len = 1;

    gr_poly_fit_length(res, len, ctx);
    status |= _gr_poly_atanh_series(res->coeffs, f->coeffs, flen, len, ctx);
    _gr_poly_set_length(res, len, ctx);
    _gr_poly_normalise(res, ctx);
    return status;
}
