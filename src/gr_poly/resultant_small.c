/*
    Copyright (C) 2023 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "gr_vec.h"
#include "gr_poly.h"

int
_gr_poly_resultant_small(gr_ptr res, gr_srcptr poly1, slong len1,
                    gr_srcptr poly2, slong len2, gr_ctx_t ctx)
{
    if (poly1 == poly2 && len1 == len2)
    {
        return gr_zero(res, ctx);
    }
    else if (len2 == 1)
    {
        if (len1 == 1)
            return gr_one(res, ctx);
        else if (len1 == 2)
            return gr_set(res, poly2, ctx);
        else
            return gr_pow_ui(res, poly2, len1 - 1, ctx);
    }
    else if (len1 == 2 && len2 == 2)
    {
        /* resultant(a0 + a1*x, b0 + b1*x) = a1*b0 - a0*b1 */
        /* todo: use gr_fmms when available */
        gr_ptr t, u;
        int status = GR_SUCCESS;
        slong sz = ctx->sizeof_elem;

        GR_TMP_INIT2(t, u, ctx);

        status |= gr_mul(t, GR_ENTRY(poly1, 1, sz), GR_ENTRY(poly2, 0, sz), ctx);
        status |= gr_mul(u, GR_ENTRY(poly1, 0, sz), GR_ENTRY(poly2, 1, sz), ctx);
        status |= gr_sub(res, t, u, ctx);

        GR_TMP_CLEAR2(t, u, ctx);

        return status;
    }
    else if (len1 == 3 && len2 == 2)
    {
        /* res(a0 + a1*x + a2*x**2, b0 + b1*x) = b1*(a0*b1 - a1*b0) + a2*b0^2 */
        gr_ptr t, u;
        int status = GR_SUCCESS;
        slong sz = ctx->sizeof_elem;

        GR_TMP_INIT2(t, u, ctx);

        status |= gr_mul(t, GR_ENTRY(poly1, 0, sz), GR_ENTRY(poly2, 1, sz), ctx);
        status |= gr_mul(u, GR_ENTRY(poly1, 1, sz), GR_ENTRY(poly2, 0, sz), ctx);
        status |= gr_sub(t, t, u, ctx);
        status |= gr_mul(t, t, GR_ENTRY(poly2, 1, sz), ctx);
        status |= gr_sqr(u, GR_ENTRY(poly2, 0, sz), ctx);
        status |= gr_mul(u, u, GR_ENTRY(poly1, 2, sz), ctx);
        status |= gr_add(res, t, u, ctx);

        GR_TMP_CLEAR2(t, u, ctx);

        return status;
    }
    else if (len2 == 2)
    {
        gr_ptr t, u, v;
        int status = GR_SUCCESS;
        slong sz = ctx->sizeof_elem;
        slong i;

        GR_TMP_INIT_VEC(t, len1 + 2, ctx);
        u = GR_ENTRY(t, len1, sz);
        v = GR_ENTRY(t, len1 + 1, sz);

        status |= _gr_vec_set_powers(t, poly2, len1, ctx);

        status |= gr_neg(u, GR_ENTRY(poly2, 1, sz), ctx);
        status |= gr_set(v, u, ctx);

        for (i = 1; i < len1; i++)
        {
            status |= gr_mul(GR_ENTRY(t, len1 - i - 1, sz), GR_ENTRY(t, len1 - i - 1, sz), u, ctx);
            if (i < len1 - 1)
                status |= gr_mul(u, u, v, ctx);
        }

        status |= _gr_vec_dot(res, NULL, 0, poly1, t, len1, ctx);

        GR_TMP_CLEAR_VEC(t, len1 + 2, ctx);

        return status;
    }
    else if (len1 == 3 && len2 == 3)
    {
        /* res(a0 + a1*x + a2*x^2, b0 + b1*x + b2*x^2, x) =

            t0=a0*a2; t1=b0*b2; t2=a0*b2; t3=a2*b0; t7=a1*b1;
            t0*(b1**2-2*t1) + t1*a1**2 + t2*(t2-t7) + t3*(t3-t7)
        */
        gr_ptr t;
        int status = GR_SUCCESS;
        slong sz = ctx->sizeof_elem;

        GR_TMP_INIT_VEC(t, 8, ctx);

        status |= gr_mul(GR_ENTRY(t, 0, sz), GR_ENTRY(poly1, 0, sz), GR_ENTRY(poly1, 2, sz), ctx);
        status |= gr_mul(GR_ENTRY(t, 1, sz), GR_ENTRY(poly2, 0, sz), GR_ENTRY(poly2, 2, sz), ctx);
        status |= gr_mul(GR_ENTRY(t, 2, sz), GR_ENTRY(poly1, 0, sz), GR_ENTRY(poly2, 2, sz), ctx);
        status |= gr_mul(GR_ENTRY(t, 3, sz), GR_ENTRY(poly1, 2, sz), GR_ENTRY(poly2, 0, sz), ctx);
        status |= gr_mul(GR_ENTRY(t, 7, sz), GR_ENTRY(poly1, 1, sz), GR_ENTRY(poly2, 1, sz), ctx);
        status |= gr_sqr(GR_ENTRY(t, 4, sz), GR_ENTRY(poly2, 1, sz), ctx);
        status |= gr_submul_ui(GR_ENTRY(t, 4, sz), GR_ENTRY(t, 1, sz), 2, ctx);
        status |= gr_sqr(GR_ENTRY(t, 5, sz), GR_ENTRY(poly1, 1, sz), ctx);
        status |= gr_sub(GR_ENTRY(t, 6, sz), GR_ENTRY(t, 2, sz), GR_ENTRY(t, 7, sz), ctx);
        status |= gr_sub(GR_ENTRY(t, 7, sz), GR_ENTRY(t, 3, sz), GR_ENTRY(t, 7, sz), ctx);
        status |= _gr_vec_dot(res, NULL, 0, t, GR_ENTRY(t, 4, sz), 4, ctx);

        GR_TMP_CLEAR_VEC(t, 8, ctx);

        return status;
    }
    else
    {
        return GR_UNABLE;
    }
}

int gr_poly_resultant_small(gr_ptr r, const gr_poly_t f,
                             const gr_poly_t g, gr_ctx_t ctx)
{
    slong len1 = f->length;
    slong len2 = g->length;
    int status = GR_SUCCESS;
    slong sz = ctx->sizeof_elem;

    if (len1 == 0 || len2 == 0)
    {
        return gr_zero(r, ctx);
    }

    if (gr_is_zero(GR_ENTRY(f->coeffs, len1 - 1, sz), ctx) != T_FALSE ||
        gr_is_zero(GR_ENTRY(g->coeffs, len2 - 1, sz), ctx) != T_FALSE)
    {
        return GR_UNABLE;
    }

    if (len1 >= len2)
    {
        status |= _gr_poly_resultant_small(r, f->coeffs, len1,  g->coeffs, len2, ctx);
    }
    else
    {
        status |= _gr_poly_resultant_small(r, g->coeffs, len2, f->coeffs, len1, ctx);

        if (((len1 | len2) & 1) == 0)
            status |= gr_neg(r, r, ctx);
    }

    return status;
}
