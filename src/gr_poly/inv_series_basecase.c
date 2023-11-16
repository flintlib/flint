/*
    Copyright (C) 2023 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "gr_vec.h"
#include "gr_poly.h"

/* Ainv is the inverse constant term of A; it is allowed to be aliased
   with res[0]. */
int
_gr_poly_inv_series_basecase_preinv1(gr_ptr res, gr_srcptr A, slong Alen, gr_srcptr Ainv, slong len, gr_ctx_t ctx)
{
    slong sz = ctx->sizeof_elem;
    int status = GR_SUCCESS;

    if (len == 0)
        return GR_SUCCESS;

    if (Alen == 0)
        return GR_DOMAIN;

    Alen = FLINT_MIN(Alen, len);

    if (Ainv != res)
        status |= gr_set(res, Ainv, ctx);

    if (Alen == 1)
    {
        status |= _gr_vec_zero(GR_ENTRY(res, 1, sz), len - 1, ctx);
    }
    else if (len == 2)
    {
        status |= gr_mul(GR_ENTRY(res, 1, sz), res, res, ctx);
        status |= gr_mul(GR_ENTRY(res, 1, sz), GR_ENTRY(res, 1, sz), GR_ENTRY(A, 1, sz), ctx);
        status |= gr_neg(GR_ENTRY(res, 1, sz), GR_ENTRY(res, 1, sz), ctx);
    }
    else if (Alen == 2 || _gr_vec_is_zero(GR_ENTRY(A, 1, sz), Alen - 2, ctx) == T_TRUE)
    {
        /* Special-case binomials */
        /* todo: implement using vector functions (geometric series + inflate) */
        slong i, j, step;

        step = Alen - 1;

        if (gr_is_one(res, ctx) == T_TRUE)
        {
            status |= gr_neg(GR_ENTRY(res, step, sz), GR_ENTRY(A, step, sz), ctx);
            for (i = 2 * step; i < len; i += step)
                status |= gr_mul(GR_ENTRY(res, i, sz), GR_ENTRY(res, i - step, sz), GR_ENTRY(res, step, sz), ctx);
        }
        else if (gr_is_neg_one(res, ctx) == T_TRUE)
        {
            status |= gr_neg(GR_ENTRY(res, step, sz), GR_ENTRY(A, step, sz), ctx);
            for (i = 2 * step; i < len; i += step)
                status |= gr_mul(GR_ENTRY(res, i, sz), GR_ENTRY(res, i - step, sz), GR_ENTRY(A, step, sz), ctx);
        }
        else
        {
            gr_ptr t;
            GR_TMP_INIT(t, ctx);

            status |= gr_mul(t, res, GR_ENTRY(A, step, sz), ctx);
            status |= gr_neg(t, t, ctx);

            for (i = step; i < len; i += step)
                status |= gr_mul(GR_ENTRY(res, i, sz), GR_ENTRY(res, i - step, sz), t, ctx);

            GR_TMP_CLEAR(t, ctx);
        }

        for (i = 0; i < len; i += step)
            for (j = i + 1; j < FLINT_MIN(len, i + step); j++)
                status |= gr_zero(GR_ENTRY(res, j, sz), ctx);
    }
    else
    {
        int is_one;
        slong i, l;

        is_one = (gr_is_one(res, ctx) == T_TRUE);

        for (i = 1; i < len; i++)
        {
            l = FLINT_MIN(i, Alen - 1);

            status |= _gr_vec_dot_rev(GR_ENTRY(res, i, sz), NULL, 1, GR_ENTRY(A, 1, sz), GR_ENTRY(res, i - l, sz), l, ctx);

            if (!is_one)
                status |= gr_mul(GR_ENTRY(res, i, sz), GR_ENTRY(res, i, sz), res, ctx);
        }
    }

    return status;
}

int
_gr_poly_inv_series_basecase_generic(gr_ptr res, gr_srcptr A, slong Alen, slong len, gr_ctx_t ctx)
{
    int status = GR_SUCCESS;

    if (len == 0)
        return GR_SUCCESS;

    if (Alen == 0)
        return GR_DOMAIN;

    Alen = FLINT_MIN(Alen, len);

    status = gr_inv(res, A, ctx);

    if (status != GR_SUCCESS)
        return status;

    return _gr_poly_inv_series_basecase_preinv1(res, A, Alen, res, len, ctx);
}

int
gr_poly_inv_series_basecase(gr_poly_t res, const gr_poly_t A, slong len, gr_ctx_t ctx)
{
    int status = GR_SUCCESS;
    slong Alen;

    if (len == 0)
        return gr_poly_zero(res, ctx);

    Alen = A->length;

    if (Alen == 0)
        return GR_DOMAIN;

    if (Alen == 1)
        len = 1;

    if (res == A)
    {
        gr_poly_t t;
        gr_poly_init(t, ctx);
        status = gr_poly_inv_series_basecase(t, A, len, ctx);
        gr_poly_swap(res, t, ctx);
        gr_poly_clear(t, ctx);
        return status;
    }

    gr_poly_fit_length(res, len, ctx);
    status |= _gr_poly_inv_series_basecase(res->coeffs, A->coeffs, A->length, len, ctx);
    _gr_poly_set_length(res, len, ctx);
    _gr_poly_normalise(res, ctx);
    return status;
}
