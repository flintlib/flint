/*
    Copyright (C) 2025 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <string.h>
#include "gr.h"
#include "gr_vec.h"
#include "gr_poly.h"

#define RE(xx) (xx)
#define IM(xx) (GR_ENTRY(xx, 1, real_sz))

int
_gr_poly_mullow_complex_reorder(gr_ptr res,
    gr_srcptr poly1, slong len1,
    gr_srcptr poly2, slong len2, slong n,
    int karatsuba, gr_ctx_t ctx, gr_ctx_t real_ctx)
{
    gr_ptr a, b, c, d, e, f, w;
    gr_ptr t, u, v;
    slong i, alloc;
    slong sz = ctx->sizeof_elem;
    slong real_sz = real_ctx->sizeof_elem;
    int status = GR_SUCCESS;
    int squaring;

    len1 = FLINT_MIN(len1, n);
    len2 = FLINT_MIN(len2, n);

    squaring = (poly1 == poly2) && (len1 == len2);

    /* Shallow reorder arrays */
    alloc = 2 * (len1 + (len2 * !squaring) + n);
    w = GR_TMP_ALLOC(alloc * real_sz);
    a = w;
    b = GR_ENTRY(a, len1, real_sz);
    c = squaring ? a : GR_ENTRY(b, len1, real_sz);
    d = squaring ? b : GR_ENTRY(c, len2, real_sz);
    e = GR_ENTRY(d, len2, real_sz);
    f = GR_ENTRY(e, n, real_sz);

    for (i = 0; i < len1; i++)
    {
        gr_set_shallow(GR_ENTRY(a, i, real_sz), RE(GR_ENTRY(poly1, i, sz)), real_ctx);
        gr_set_shallow(GR_ENTRY(b, i, real_sz), IM(GR_ENTRY(poly1, i, sz)), real_ctx);
    }

    for (i = 0; i < n; i++)
    {
        gr_set_shallow(GR_ENTRY(e, i, real_sz), RE(GR_ENTRY(res, i, sz)), real_ctx);
        gr_set_shallow(GR_ENTRY(f, i, real_sz), IM(GR_ENTRY(res, i, sz)), real_ctx);
    }

    if (squaring)
    {
        if (_gr_vec_is_zero(b, len1, real_ctx) == T_TRUE)
        {
            status |= _gr_poly_mullow(e, a, len1, a, len1, n, real_ctx);
            status |= _gr_vec_zero(f, n, real_ctx);
        }
        else if (_gr_vec_is_zero(a, len1, real_ctx) == T_TRUE)
        {
            status |= _gr_poly_mullow(e, b, len1, b, len1, n, real_ctx);
            status |= _gr_vec_neg(e, e, n, real_ctx);
            status |= _gr_vec_zero(f, n, real_ctx);
        }
        else if (karatsuba)
        {
            GR_TMP_INIT_VEC(t, 3 * n, real_ctx);
            u = GR_ENTRY(t, n, real_sz);
            v = GR_ENTRY(u, n, real_sz);

            status |= _gr_vec_add(t, a, b, len1, real_ctx);
            status |= _gr_poly_mullow(v, t, len1, t, len1, n, real_ctx);
            status |= _gr_poly_mullow(t, a, len1, a, len1, n, real_ctx);
            status |= _gr_poly_mullow(u, b, len1, b, len1, n, real_ctx);
            status |= _gr_vec_sub(e, t, u, n, real_ctx);
            status |= _gr_vec_sub(f, v, t, n, real_ctx);
            status |= _gr_vec_sub(f, f, u, n, real_ctx);

            GR_TMP_CLEAR_VEC(t, 3 * n, real_ctx);
        }
        else
        {
            status |= _gr_poly_mullow(e, a, len1, a, len1, n, real_ctx);
            status |= _gr_poly_mullow(f, b, len1, b, len1, n, real_ctx);
            status |= _gr_vec_sub(e, e, f, n, real_ctx);
            status |= _gr_poly_mullow(f, a, len1, b, len1, n, real_ctx);
            status |= _gr_vec_mul_scalar_2exp_si(f, f, n, 1, real_ctx);
        }
    }
    else
    {
        for (i = 0; i < len2; i++)
        {
            gr_set_shallow(GR_ENTRY(c, i, real_sz), RE(GR_ENTRY(poly2, i, sz)), real_ctx);
            gr_set_shallow(GR_ENTRY(d, i, real_sz), IM(GR_ENTRY(poly2, i, sz)), real_ctx);
        }

        if (_gr_vec_is_zero(b, len1, real_ctx) == T_TRUE)
        {
            if (_gr_vec_is_zero(d, len2, real_ctx) == T_TRUE)
                status |= _gr_vec_zero(f, n, real_ctx);
            else
                status |= _gr_poly_mullow(f, a, len1, d, len2, n, real_ctx);

            if (_gr_vec_is_zero(c, len2, real_ctx) == T_TRUE)
                status |= _gr_vec_zero(e, n, real_ctx);
            else
                status |= _gr_poly_mullow(e, a, len1, c, len2, n, real_ctx);
        }
        else if (_gr_vec_is_zero(a, len1, real_ctx) == T_TRUE)
        {
            if (_gr_vec_is_zero(d, len2, real_ctx) == T_TRUE)
                status |= _gr_vec_zero(e, n, real_ctx);
            else
            {
                status |= _gr_poly_mullow(e, b, len1, d, len2, n, real_ctx);
                status |= _gr_vec_neg(e, e, n, real_ctx);
            }

            if (_gr_vec_is_zero(c, len2, real_ctx) == T_TRUE)
                status |= _gr_vec_zero(f, n, real_ctx);
            else
                status |= _gr_poly_mullow(f, b, len1, c, len2, n, real_ctx);
        }
        else if (_gr_vec_is_zero(d, len2, real_ctx) == T_TRUE)
        {
            status |= _gr_poly_mullow(e, a, len1, c, len2, n, real_ctx);
            status |= _gr_poly_mullow(f, b, len1, c, len2, n, real_ctx);
        }
        else if (_gr_vec_is_zero(c, len2, real_ctx) == T_TRUE)
        {
            status |= _gr_poly_mullow(e, b, len1, d, len2, n, real_ctx);
            status |= _gr_vec_neg(e, e, n, real_ctx);
            status |= _gr_poly_mullow(f, a, len1, d, len2, n, real_ctx);
        }
        else if (karatsuba)
        {
            GR_TMP_INIT_VEC(t, 3 * n, real_ctx);
            u = GR_ENTRY(t, n, real_sz);
            v = GR_ENTRY(u, n, real_sz);

            status |= _gr_vec_add(t, a, b, len1, real_ctx);
            status |= _gr_vec_add(u, c, d, len2, real_ctx);
            status |= _gr_poly_mullow(v, t, len1, u, len2, n, real_ctx);
            status |= _gr_poly_mullow(t, a, len1, c, len2, n, real_ctx);
            status |= _gr_poly_mullow(u, b, len1, d, len2, n, real_ctx);
            status |= _gr_vec_sub(e, t, u, n, real_ctx);
            status |= _gr_vec_sub(f, v, t, n, real_ctx);
            status |= _gr_vec_sub(f, f, u, n, real_ctx);

            GR_TMP_CLEAR_VEC(t, 3 * n, real_ctx);
        }
        else
        {
            GR_TMP_INIT_VEC(t, n, real_ctx);

            status |= _gr_poly_mullow(e, a, len1, c, len2, n, real_ctx);
            status |= _gr_poly_mullow(t, b, len1, d, len2, n, real_ctx);
            status |= _gr_vec_sub(e, e, t, n, real_ctx);
            status |= _gr_poly_mullow(f, a, len1, d, len2, n, real_ctx);
            status |= _gr_poly_mullow(t, b, len1, c, len2, n, real_ctx);
            status |= _gr_vec_add(f, f, t, n, real_ctx);

            GR_TMP_CLEAR_VEC(t, n, real_ctx);
        }
    }

    for (i = 0; i < n; i++)
    {
        gr_set_shallow(RE(GR_ENTRY(res, i, sz)), GR_ENTRY(e, i, real_sz), real_ctx);
        gr_set_shallow(IM(GR_ENTRY(res, i, sz)), GR_ENTRY(f, i, real_sz), real_ctx);
    }

    GR_TMP_FREE(w, alloc * real_sz);

    return status;
}

int
gr_poly_mullow_complex_reorder(gr_poly_t res, const gr_poly_t poly1,
                                            const gr_poly_t poly2,
                                                slong n, int karatsuba, gr_ctx_t ctx, gr_ctx_t real_ctx)
{
    slong len_out;
    int status;

    if (poly1->length == 0 || poly2->length == 0 || n == 0)
        return gr_poly_zero(res, ctx);

    len_out = poly1->length + poly2->length - 1;
    n = FLINT_MIN(n, len_out);

    if (res == poly1 || res == poly2)
    {
        gr_poly_t t;
        gr_poly_init2(t, n, ctx);
        status = _gr_poly_mullow_complex_reorder(t->coeffs,
            poly1->coeffs, poly1->length, poly2->coeffs, poly2->length, n, karatsuba, ctx, real_ctx);
        gr_poly_swap(res, t, ctx);
        gr_poly_clear(t, ctx);
    }
    else
    {
        gr_poly_fit_length(res, n, ctx);
        status = _gr_poly_mullow_complex_reorder(res->coeffs,
            poly1->coeffs, poly1->length, poly2->coeffs, poly2->length, n, karatsuba, ctx, real_ctx);
    }

    _gr_poly_set_length(res, n, ctx);
    _gr_poly_normalise(res, ctx);
    return status;
}
