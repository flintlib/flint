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
_gr_poly_mulmid_complex_reorder(gr_ptr res,
    gr_srcptr poly1, slong len1,
    gr_srcptr poly2, slong len2, slong nlo, slong nhi,
    int karatsuba, gr_ctx_t ctx, gr_ctx_t real_ctx)
{
    gr_ptr a, b, c, d, e, f, w;
    gr_ptr t, u, v;
    slong i, alloc;
    slong sz = ctx->sizeof_elem;
    slong real_sz = real_ctx->sizeof_elem;
    slong n = nhi - nlo;
    int status = GR_SUCCESS;
    int squaring;

    len1 = FLINT_MIN(len1, nhi);
    len2 = FLINT_MIN(len2, nhi);

    if (nlo != 0)
    {
        slong nlo2 = (len1 + len2 - 1) - nlo;

        if (len1 > nlo2)
        {
            slong trunc = len1 - nlo2;
            poly1 = GR_ENTRY(poly1, trunc, sz);
            len1 -= trunc;
            nlo -= trunc;
            nhi -= trunc;
        }

        if (len2 > nlo2)
        {
            slong trunc = len2 - nlo2;
            poly2 = GR_ENTRY(poly2, trunc, sz);
            len2 -= trunc;
            nlo -= trunc;
            nhi -= trunc;
        }
    }

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
            status |= _gr_poly_mulmid(e, a, len1, a, len1, nlo, nhi, real_ctx);
            status |= _gr_vec_zero(f, n, real_ctx);
        }
        else if (_gr_vec_is_zero(a, len1, real_ctx) == T_TRUE)
        {
            status |= _gr_poly_mulmid(e, b, len1, b, len1, nlo, nhi, real_ctx);
            status |= _gr_vec_neg(e, e, n, real_ctx);
            status |= _gr_vec_zero(f, n, real_ctx);
        }
        else if (karatsuba)
        {
            slong tn = FLINT_MAX(len1, n);

            GR_TMP_INIT_VEC(t, 3 * tn, real_ctx);
            u = GR_ENTRY(t, tn, real_sz);
            v = GR_ENTRY(u, tn, real_sz);

            status |= _gr_vec_add(t, a, b, len1, real_ctx);
            status |= _gr_poly_mulmid(v, t, len1, t, len1, nlo, nhi, real_ctx);
            status |= _gr_poly_mulmid(t, a, len1, a, len1, nlo, nhi, real_ctx);
            status |= _gr_poly_mulmid(u, b, len1, b, len1, nlo, nhi, real_ctx);
            status |= _gr_vec_sub(e, t, u, n, real_ctx);
            status |= _gr_vec_sub(f, v, t, n, real_ctx);
            status |= _gr_vec_sub(f, f, u, n, real_ctx);

            GR_TMP_CLEAR_VEC(t, 3 * tn, real_ctx);
        }
        else
        {
            status |= _gr_poly_mulmid(e, a, len1, a, len1, nlo, nhi, real_ctx);
            status |= _gr_poly_mulmid(f, b, len1, b, len1, nlo, nhi, real_ctx);
            status |= _gr_vec_sub(e, e, f, n, real_ctx);
            status |= _gr_poly_mulmid(f, a, len1, b, len1, nlo, nhi, real_ctx);
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
                status |= _gr_poly_mulmid(f, a, len1, d, len2, nlo, nhi, real_ctx);

            if (_gr_vec_is_zero(c, len2, real_ctx) == T_TRUE)
                status |= _gr_vec_zero(e, n, real_ctx);
            else
                status |= _gr_poly_mulmid(e, a, len1, c, len2, nlo, nhi, real_ctx);
        }
        else if (_gr_vec_is_zero(a, len1, real_ctx) == T_TRUE)
        {
            if (_gr_vec_is_zero(d, len2, real_ctx) == T_TRUE)
                status |= _gr_vec_zero(e, n, real_ctx);
            else
            {
                status |= _gr_poly_mulmid(e, b, len1, d, len2, nlo, nhi, real_ctx);
                status |= _gr_vec_neg(e, e, n, real_ctx);
            }

            if (_gr_vec_is_zero(c, len2, real_ctx) == T_TRUE)
                status |= _gr_vec_zero(f, n, real_ctx);
            else
                status |= _gr_poly_mulmid(f, b, len1, c, len2, nlo, nhi, real_ctx);
        }
        else if (_gr_vec_is_zero(d, len2, real_ctx) == T_TRUE)
        {
            status |= _gr_poly_mulmid(e, a, len1, c, len2, nlo, nhi, real_ctx);
            status |= _gr_poly_mulmid(f, b, len1, c, len2, nlo, nhi, real_ctx);
        }
        else if (_gr_vec_is_zero(c, len2, real_ctx) == T_TRUE)
        {
            status |= _gr_poly_mulmid(e, b, len1, d, len2, nlo, nhi, real_ctx);
            status |= _gr_vec_neg(e, e, n, real_ctx);
            status |= _gr_poly_mulmid(f, a, len1, d, len2, nlo, nhi, real_ctx);
        }
        else if (karatsuba)
        {
            slong tn = FLINT_MAX(FLINT_MAX(len1, len2), n);

            GR_TMP_INIT_VEC(t, 3 * tn, real_ctx);
            u = GR_ENTRY(t, tn, real_sz);
            v = GR_ENTRY(u, tn, real_sz);

            status |= _gr_vec_add(t, a, b, len1, real_ctx);
            status |= _gr_vec_add(u, c, d, len2, real_ctx);
            status |= _gr_poly_mulmid(v, t, len1, u, len2, nlo, nhi, real_ctx);
            status |= _gr_poly_mulmid(t, a, len1, c, len2, nlo, nhi, real_ctx);
            status |= _gr_poly_mulmid(u, b, len1, d, len2, nlo, nhi, real_ctx);
            status |= _gr_vec_sub(e, t, u, n, real_ctx);
            status |= _gr_vec_sub(f, v, t, n, real_ctx);
            status |= _gr_vec_sub(f, f, u, n, real_ctx);

            GR_TMP_CLEAR_VEC(t, 3 * tn, real_ctx);
        }
        else
        {
            GR_TMP_INIT_VEC(t, n, real_ctx);

            status |= _gr_poly_mulmid(e, a, len1, c, len2, nlo, nhi, real_ctx);
            status |= _gr_poly_mulmid(t, b, len1, d, len2, nlo, nhi, real_ctx);
            status |= _gr_vec_sub(e, e, t, n, real_ctx);
            status |= _gr_poly_mulmid(f, a, len1, d, len2, nlo, nhi, real_ctx);
            status |= _gr_poly_mulmid(t, b, len1, c, len2, nlo, nhi, real_ctx);
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
_gr_poly_mullow_complex_reorder(gr_ptr res,
    gr_srcptr poly1, slong len1,
    gr_srcptr poly2, slong len2, slong n,
    int karatsuba, gr_ctx_t ctx, gr_ctx_t real_ctx)
{
    return _gr_poly_mulmid_complex_reorder(res, poly1, len1, poly2, len2, 0, n,
        karatsuba, ctx, real_ctx);
}

int
gr_poly_mulmid_complex_reorder(gr_poly_t res, const gr_poly_t poly1,
                                            const gr_poly_t poly2,
                                                slong nlo, slong nhi, int karatsuba, gr_ctx_t ctx, gr_ctx_t real_ctx)
{
    slong len1 = poly1->length;
    slong len2 = poly2->length;
    int status;
    slong len;

    FLINT_ASSERT(nlo >= 0);
    FLINT_ASSERT(nhi >= 0);

    if (len1 == 0 || len2 == 0 || nlo >= FLINT_MIN(nhi, len1 + len2 - 1))
        return gr_poly_zero(res, ctx);

    nhi = FLINT_MIN(nhi, len1 + len2 - 1);
    len = nhi - nlo;

    if (res == poly1 || res == poly2)
    {
        gr_poly_t t;
        gr_poly_init2(t, len, ctx);
        status = _gr_poly_mulmid_complex_reorder(t->coeffs, poly1->coeffs, len1, poly2->coeffs, len2, nlo, nhi, karatsuba, ctx, real_ctx);
        gr_poly_swap(res, t, ctx);
        gr_poly_clear(t, ctx);
    }
    else
    {
        gr_poly_fit_length(res, len, ctx);
        status = _gr_poly_mulmid_complex_reorder(res->coeffs, poly1->coeffs, len1, poly2->coeffs, len2, nlo, nhi, karatsuba, ctx, real_ctx);
    }

    _gr_poly_set_length(res, len, ctx);
    _gr_poly_normalise(res, ctx);
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
