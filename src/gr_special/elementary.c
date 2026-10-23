/*
    Copyright (C) 2023 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "gr_special.h"

int
gr_generic_exp(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)
{
    if (gr_is_zero(x, ctx) == T_TRUE)
        return gr_one(res, ctx);

    return GR_UNABLE;
}

int
gr_generic_expm1(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)
{
    int status = gr_exp(res, x, ctx);
    status |= gr_sub_ui(res, res, 1, ctx);
    return status;
}

int
gr_generic_exp2(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)
{
    gr_ptr t;
    int status = GR_SUCCESS;

    GR_TMP_INIT(t, ctx);
    status |= gr_set_ui(t, 2, ctx);
    status |= gr_pow(res, t, x, ctx);
    GR_TMP_CLEAR(t, ctx);

    return status;
}

int
gr_generic_exp10(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)
{
    gr_ptr t;
    int status = GR_SUCCESS;

    GR_TMP_INIT(t, ctx);
    status |= gr_set_ui(t, 10, ctx);
    status |= gr_pow(res, t, x, ctx);
    GR_TMP_CLEAR(t, ctx);

    return status;
}

int
gr_generic_log(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)
{
    if (gr_is_one(x, ctx) == T_TRUE)
        return gr_zero(res, ctx);

    return GR_UNABLE;
}

int
gr_generic_log1p(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)
{
    int status = gr_add_ui(res, x, 1, ctx);
    status |= gr_log(res, res, ctx);
    return status;
}

int
gr_generic_log2(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)
{
    gr_ptr t;
    int status = GR_SUCCESS;

    GR_TMP_INIT(t, ctx);
    status |= gr_set_ui(t, 2, ctx);
    status |= gr_log(t, t, ctx);
    status |= gr_log(res, x, ctx);
    status |= gr_div(res, res, t, ctx);
    GR_TMP_CLEAR(t, ctx);

    if (status != GR_SUCCESS)
        status = GR_UNABLE;

    return status;
}

int
gr_generic_log10(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)
{
    gr_ptr t;
    int status = GR_SUCCESS;

    GR_TMP_INIT(t, ctx);
    status |= gr_set_ui(t, 10, ctx);
    status |= gr_log(t, t, ctx);
    status |= gr_log(res, x, ctx);
    status |= gr_div(res, res, t, ctx);
    GR_TMP_CLEAR(t, ctx);

    if (status != GR_SUCCESS)
        status = GR_UNABLE;

    return status;
}

int
gr_generic_sin(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)
{
    if (gr_is_zero(x, ctx) == T_TRUE)
        return gr_zero(res, ctx);

    return GR_UNABLE;
}

int
gr_generic_cos(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)
{
    if (gr_is_zero(x, ctx) == T_TRUE)
        return gr_one(res, ctx);

    return GR_UNABLE;
}

int
gr_generic_sin_cos(gr_ptr res1, gr_ptr res2, gr_srcptr x, gr_ctx_t ctx)
{
    if (gr_is_zero(x, ctx) == T_TRUE)
        return gr_zero(res1, ctx) | gr_one(res2, ctx);
    else
        return GR_UNABLE;
}

int
gr_generic_tan(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)
{
    if (gr_is_zero(x, ctx) == T_TRUE)
        return gr_zero(res, ctx);

    return GR_UNABLE;
}

int
gr_generic_tanh(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)
{
    if (gr_is_zero(x, ctx) == T_TRUE)
        return gr_zero(res, ctx);

    return GR_UNABLE;
}

int
gr_generic_asin(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)
{
    if (gr_is_zero(x, ctx) == T_TRUE)
        return gr_zero(res, ctx);

    return GR_UNABLE;
}

int
gr_generic_asinh(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)
{
    if (gr_is_zero(x, ctx) == T_TRUE)
        return gr_zero(res, ctx);

    return GR_UNABLE;
}

int
gr_generic_atan(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)
{
    if (gr_is_zero(x, ctx) == T_TRUE)
        return gr_zero(res, ctx);

    return GR_UNABLE;
}

int
gr_generic_atanh(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)
{
    if (gr_is_zero(x, ctx) == T_TRUE)
        return gr_zero(res, ctx);

    return GR_UNABLE;
}

int
gr_generic_sec_pi(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)
{
    int status = GR_SUCCESS;

    status |= gr_cos_pi(res, x, ctx);
    status |= gr_inv(res, res, ctx);

    return status;
}

int
gr_generic_csc_pi(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)
{
    int status = GR_SUCCESS;

    status |= gr_sin_pi(res, x, ctx);
    status |= gr_inv(res, res, ctx);

    return status;
}

/* f(x) / pi for an inverse trigonometric function f */
static int
_gr_generic_div_pi(gr_ptr res, int status, gr_ctx_t ctx)
{
    gr_ptr t;

    if (status != GR_SUCCESS)
        return status;

    GR_TMP_INIT(t, ctx);
    status = gr_pi(t, ctx);
    status |= gr_div(res, res, t, ctx);
    GR_TMP_CLEAR(t, ctx);

    return status;
}

int gr_generic_asin_pi(gr_ptr res, gr_srcptr x, gr_ctx_t ctx) { return _gr_generic_div_pi(res, gr_asin(res, x, ctx), ctx); }
int gr_generic_acos_pi(gr_ptr res, gr_srcptr x, gr_ctx_t ctx) { return _gr_generic_div_pi(res, gr_acos(res, x, ctx), ctx); }
int gr_generic_atan_pi(gr_ptr res, gr_srcptr x, gr_ctx_t ctx) { return _gr_generic_div_pi(res, gr_atan(res, x, ctx), ctx); }
int gr_generic_acot_pi(gr_ptr res, gr_srcptr x, gr_ctx_t ctx) { return _gr_generic_div_pi(res, gr_acot(res, x, ctx), ctx); }
int gr_generic_asec_pi(gr_ptr res, gr_srcptr x, gr_ctx_t ctx) { return _gr_generic_div_pi(res, gr_asec(res, x, ctx), ctx); }
int gr_generic_acsc_pi(gr_ptr res, gr_srcptr x, gr_ctx_t ctx) { return _gr_generic_div_pi(res, gr_acsc(res, x, ctx), ctx); }

int
gr_generic_acot(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)
{
    int status = GR_SUCCESS;

    status |= gr_inv(res, x, ctx);
    status |= gr_atan(res, res, ctx);

    return status;
}

int
gr_generic_asec(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)
{
    int status = GR_SUCCESS;

    status |= gr_inv(res, x, ctx);
    status |= gr_acos(res, res, ctx);

    return status;
}

int
gr_generic_acsc(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)
{
    int status = GR_SUCCESS;

    status |= gr_inv(res, x, ctx);
    status |= gr_asin(res, res, ctx);

    return status;
}

int
gr_generic_acoth(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)
{
    int status = GR_SUCCESS;

    status |= gr_inv(res, x, ctx);
    status |= gr_atanh(res, res, ctx);

    return status;
}

int
gr_generic_asech(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)
{
    int status = GR_SUCCESS;

    status |= gr_inv(res, x, ctx);
    status |= gr_acosh(res, res, ctx);

    return status;
}

int
gr_generic_acsch(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)
{
    int status = GR_SUCCESS;

    status |= gr_inv(res, x, ctx);
    status |= gr_asinh(res, res, ctx);

    return status;
}
