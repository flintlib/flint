/*
    Copyright (C) 2007, 2008 William Hart
    Copyright (C) 2011 Sebastian Pancratz
    Copyright (C) 2014 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "gr_vec.h"
#include "gr_poly.h"

#define GR_VEC_NORM(status, R, lenR, sz, ctx) \
    do { \
        (void) sz; \
        (status) |= _gr_vec_normalise(&(lenR), (R), (lenR), (ctx)); \
    } while (0)

int
_gr_poly_resultant_euclidean(gr_ptr res, gr_srcptr poly1, slong len1,
                    gr_srcptr poly2, slong len2, gr_ctx_t ctx)
{
    gr_ptr u, v, r, t, w, q;
    slong l0, l1, l2;
    gr_ptr lc;
    int status = GR_SUCCESS;
    slong sz = ctx->sizeof_elem;

    if (len2 == 1)
        return _gr_poly_resultant_small(res, poly1, len1, poly2, len2, ctx);

    /* todo: we can skip the q tmp allocation when we have a dedicated rem() */
    GR_TMP_INIT_VEC(w, 4 * len1 + 1, ctx);

    q = w;
    u = GR_ENTRY(q, len1, sz);
    v = GR_ENTRY(u, len1, sz);
    r = GR_ENTRY(v, len1, sz);
    lc = GR_ENTRY(r, len1, sz);

    status |= gr_one(res, ctx);

    status |= _gr_vec_set(u, poly1, len1, ctx);
    status |= _gr_vec_set(v, poly2, len2, ctx);
    l1 = len1;
    l2 = len2;

    do
    {
        l0 = l1;
        l1 = l2;
        status |= gr_set(lc, GR_ENTRY(v, l1 - 1, sz), ctx);
        /* todo: just rem */
        status |= _gr_poly_divrem(q, r, u, l0, v, l1, ctx);

        if (status != GR_SUCCESS)
            break;

        l2 = l1 - 1;

        GR_VEC_NORM(status, r, l2, sz, ctx);

        {
            t = u;
            u = v;
            v = r;
            r = t;
        }

        if (l2 >= 1)
        {
            status |= gr_pow_ui(lc, lc, l0 - l2, ctx);
            status |= gr_mul(res, res, lc, ctx);

            if (((l0 | l1) & 1) == 0)
                status |= gr_neg(res, res, ctx);
        }
        else
        {
            if (l1 == 1)
            {
                status |= gr_pow_ui(lc, lc, l0 - 1, ctx);
                status |= gr_mul(res, res, lc, ctx);
            }
            else
            {
                status |= gr_zero(res, ctx);
            }
        }
    }
    while (l2 > 0);

    GR_TMP_CLEAR_VEC(w, 4 * len1 + 1, ctx);

    return status;
}

int gr_poly_resultant_euclidean(gr_ptr r, const gr_poly_t f,
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
        status |= _gr_poly_resultant_euclidean(r, f->coeffs, len1,  g->coeffs, len2, ctx);
    }
    else
    {
        status |= _gr_poly_resultant_euclidean(r, g->coeffs, len2, f->coeffs, len1, ctx);

        if (((len1 | len2) & 1) == 0)
            status |= gr_neg(r, r, ctx);
    }

    return status;
}
