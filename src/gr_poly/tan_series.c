/*
    Copyright (C) 2023, 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "gr_special.h"
#include "gr_poly.h"

/*
with sign/scale changes to compute the following functions:

    func = 0 -> tan
    func = 1 -> tanh
    func = 2 -> cot
    func = 3 -> coth
    func = 4 -> tan_pi
    func = 5 -> tanh_pi (not yet implemented/used; scalar function missing)
    func = 6 -> cot_pi
    func = 7 -> coth_pi (not yet implemented/used; scalar function missing)
*/

static int
_gr_poly_tan_series_default(gr_ptr res, gr_srcptr f, slong flen, slong len, int func, gr_ctx_t ctx)
{
    slong cutoff;

    flen = FLINT_MIN(flen, len);

    if (flen == 1)
        return _gr_poly_tan_series_basecase(res, f, flen, len, func, ctx);

    if (gr_ctx_is_finite(ctx) == T_TRUE)
    {
        /* Reasonable tuning for nmod */
        cutoff = 500;
    }
    else if (gr_ctx_has_real_prec(ctx) == T_TRUE)
    {
        gr_ptr t, u;
        int cmp = 0, want_exponential = 0;

        /* For large imaginary arguments (or large real arguments for the
           hyperbolic functions), favor using a quotient of
           (1 + {tiny exponential}) to optimize relative accuracy, even if
           this is slower than a direct tangent method. TODO: make this
           precision-dependent; if the precision loss is small compared to
           the overall precision, use the tangent method. */

        GR_TMP_INIT2(t, u, ctx);
        GR_IGNORE(gr_one(u, ctx));

        /* Hyperbolic */
        if (func & 1)
            GR_IGNORE(gr_re(t, f, ctx));
        else
            GR_IGNORE(gr_im(t, f, ctx));

        GR_IGNORE(gr_cmpabs(&cmp, t, u, ctx));
        want_exponential = (cmp == 1);

        GR_TMP_CLEAR2(t, u, ctx);

        if (want_exponential)
            return _gr_poly_tan_series_exponential(res, f, flen, len, func, ctx);

        /* Reasonable tuning for arb */
        GR_MUST_SUCCEED(gr_ctx_get_real_prec(&cutoff, ctx));
        cutoff = FLINT_MAX(cutoff, 1);
        cutoff = FLINT_MIN(500, 20 + 100000 / cutoff);
    }
    else
    {
        /* Reasonable tuning checked for fmpq */
        cutoff = 25;
    }

    if (flen < cutoff)
        return _gr_poly_tan_series_basecase(res, f, flen, len, func, ctx);
    else
        return _gr_poly_tan_series_newton(res, f, flen, len, cutoff, func, ctx);
}

int
_gr_poly_tan_series(gr_ptr res, gr_srcptr f, slong flen, slong len, gr_ctx_t ctx)
{
    return _gr_poly_tan_series_default(res, f, flen, len, 0, ctx);
}

int
_gr_poly_tanh_series(gr_ptr res, gr_srcptr f, slong flen, slong len, gr_ctx_t ctx)
{
    return _gr_poly_tan_series_default(res, f, flen, len, 1, ctx);
}

int
_gr_poly_cot_series(gr_ptr res, gr_srcptr f, slong flen, slong len, gr_ctx_t ctx)
{
    return _gr_poly_tan_series_default(res, f, flen, len, 2, ctx);
}

int
_gr_poly_coth_series(gr_ptr res, gr_srcptr f, slong flen, slong len, gr_ctx_t ctx)
{
    return _gr_poly_tan_series_default(res, f, flen, len, 3, ctx);
}

int
_gr_poly_tan_pi_series(gr_ptr res, gr_srcptr f, slong flen, slong len, gr_ctx_t ctx)
{
    return _gr_poly_tan_series_default(res, f, flen, len, 4, ctx);
}

int
_gr_poly_cot_pi_series(gr_ptr res, gr_srcptr f, slong flen, slong len, gr_ctx_t ctx)
{
    return _gr_poly_tan_series_default(res, f, flen, len, 6, ctx);
}

int
gr_poly_tan_series(gr_poly_t res, const gr_poly_t f, slong len, gr_ctx_t ctx)
{
    slong flen = f->length;
    int status = GR_SUCCESS;

    if (flen == 0 || len == 0)
        return gr_poly_zero(res, ctx);

    if (flen == 1)
        len = 1;

    gr_poly_fit_length(res, len, ctx);
    status |= _gr_poly_tan_series(res->coeffs, f->coeffs, flen, len, ctx);
    _gr_poly_set_length(res, len, ctx);
    _gr_poly_normalise(res, ctx);
    return status;
}

int
gr_poly_tanh_series(gr_poly_t res, const gr_poly_t f, slong len, gr_ctx_t ctx)
{
    slong flen = f->length;
    int status = GR_SUCCESS;

    if (flen == 0 || len == 0)
        return gr_poly_zero(res, ctx);

    if (flen == 1)
        len = 1;

    gr_poly_fit_length(res, len, ctx);
    status |= _gr_poly_tanh_series(res->coeffs, f->coeffs, flen, len, ctx);
    _gr_poly_set_length(res, len, ctx);
    _gr_poly_normalise(res, ctx);
    return status;
}

int
gr_poly_cot_series(gr_poly_t res, const gr_poly_t f, slong len, gr_ctx_t ctx)
{
    slong flen = f->length;
    int status = GR_SUCCESS;

    if (f->length == 0)
        return GR_DOMAIN;

    if (len == 0)
        return gr_poly_zero(res, ctx);

    if (flen == 1)
        len = 1;

    gr_poly_fit_length(res, len, ctx);
    status |= _gr_poly_cot_series(res->coeffs, f->coeffs, flen, len, ctx);
    _gr_poly_set_length(res, len, ctx);
    _gr_poly_normalise(res, ctx);
    return status;
}

int
gr_poly_coth_series(gr_poly_t res, const gr_poly_t f, slong len, gr_ctx_t ctx)
{
    slong flen = f->length;
    int status = GR_SUCCESS;

    if (f->length == 0)
        return GR_DOMAIN;

    if (len == 0)
        return gr_poly_zero(res, ctx);

    if (flen == 1)
        len = 1;

    gr_poly_fit_length(res, len, ctx);
    status |= _gr_poly_coth_series(res->coeffs, f->coeffs, flen, len, ctx);
    _gr_poly_set_length(res, len, ctx);
    _gr_poly_normalise(res, ctx);
    return status;
}

int
gr_poly_tan_pi_series(gr_poly_t res, const gr_poly_t f, slong len, gr_ctx_t ctx)
{
    slong flen = f->length;
    int status = GR_SUCCESS;

    if (flen == 0 || len == 0)
        return gr_poly_zero(res, ctx);

    if (flen == 1)
        len = 1;

    gr_poly_fit_length(res, len, ctx);
    status |= _gr_poly_tan_pi_series(res->coeffs, f->coeffs, flen, len, ctx);
    _gr_poly_set_length(res, len, ctx);
    _gr_poly_normalise(res, ctx);
    return status;
}

int
gr_poly_cot_pi_series(gr_poly_t res, const gr_poly_t f, slong len, gr_ctx_t ctx)
{
    slong flen = f->length;
    int status = GR_SUCCESS;

    if (f->length == 0)
        return GR_DOMAIN;

    if (len == 0)
        return gr_poly_zero(res, ctx);

    if (flen == 1)
        len = 1;

    gr_poly_fit_length(res, len, ctx);
    status |= _gr_poly_cot_pi_series(res->coeffs, f->coeffs, flen, len, ctx);
    _gr_poly_set_length(res, len, ctx);
    _gr_poly_normalise(res, ctx);
    return status;
}

