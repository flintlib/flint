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
_gr_poly_mulmid_classical(gr_ptr res,
    gr_srcptr poly1, slong len1,
    gr_srcptr poly2, slong len2, slong nlo, slong nhi, gr_ctx_t ctx)
{
    int status = GR_SUCCESS;
    slong sz = ctx->sizeof_elem;

    FLINT_ASSERT(len1 != 0);
    FLINT_ASSERT(len2 != 0);
    FLINT_ASSERT(nhi != 0);
    FLINT_ASSERT(nlo < nhi);
    FLINT_ASSERT(nlo >= 0);
    FLINT_ASSERT(nhi <= len1 + len2 - 1);

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

    if (nhi == 1)
        return gr_mul(res, poly1, poly2, ctx);

    if (len1 == 1)
        return _gr_scalar_mul_vec(res, poly1, GR_ENTRY(poly2, nlo, sz), nhi - nlo, ctx);

    if (len2 == 1)
        return _gr_vec_mul_scalar(res, GR_ENTRY(poly1, nlo, sz), nhi - nlo, poly2, ctx);

    res = GR_ENTRY(res, -nlo, sz);

    /* Squaring */
    if (poly1 == poly2 && len1 == len2 &&
        (gr_ctx_is_commutative_ring(ctx) == T_TRUE
        || gr_ctx_is_approx_commutative_ring(ctx) == T_TRUE))
    {
        slong i, start, stop;

        /* todo: double, square, addmul */

        if (nlo == 0)
            status |= gr_sqr(res, poly1, ctx);

        for (i = FLINT_MAX(1, nlo); i < FLINT_MIN(nhi, 2 * len1 - 2); i++)
        {
            start = FLINT_MAX(0, i - len1 + 1);
            stop = FLINT_MIN(len1 - 1, (i + 1) / 2 - 1);

            status |= _gr_vec_dot_rev(GR_ENTRY(res, i, sz), NULL, 0, GR_ENTRY(poly1, start, sz), GR_ENTRY(poly1, i - stop, sz), stop - start + 1, ctx);
            status |= gr_mul_two(GR_ENTRY(res, i, sz), GR_ENTRY(res, i, sz), ctx);

            if (i % 2 == 0 && i / 2 < len1)   /* should be addsqr? */
                status |= gr_addmul(GR_ENTRY(res, i, sz), GR_ENTRY(poly1, i / 2, sz), GR_ENTRY(poly1, i / 2, sz), ctx);
        }

        if (nhi >= 2 * len1 - 1)
            status |= gr_sqr(GR_ENTRY(res, 2 * len1 - 2, sz), GR_ENTRY(poly1, len1 - 1, sz), ctx);

        return status;
    }
    else
    {
        slong i, top1, top2;

        if (nlo == 0)
            status = gr_mul(res, poly1, poly2, ctx);

        for (i = FLINT_MAX(1, nlo); i < FLINT_MIN(nhi, len1 + len2 - 2); i++)
        {
            top1 = FLINT_MIN(len1 - 1, i);
            top2 = FLINT_MIN(len2 - 1, i);

            status |= _gr_vec_dot_rev(GR_ENTRY(res, i, sz), NULL, 0, GR_ENTRY(poly1, i - top2, sz), GR_ENTRY(poly2, i - top1, sz), top1 + top2 - i + 1, ctx);
        }

        if (nhi >= len1 + len2 - 1)
            status |= gr_mul(GR_ENTRY(res, len1 + len2 - 2, sz), GR_ENTRY(poly1, len1 - 1, sz), GR_ENTRY(poly2, len2 - 1, sz), ctx);

        return status;
    }
}

int
gr_poly_mulmid_classical(gr_poly_t res, const gr_poly_t poly1,
                                            const gr_poly_t poly2,
                                                slong nlo, slong nhi, gr_ctx_t ctx)
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
        status = _gr_poly_mulmid_classical(t->coeffs, poly1->coeffs, len1, poly2->coeffs, len2, nlo, nhi, ctx);
        gr_poly_swap(res, t, ctx);
        gr_poly_clear(t, ctx);
    }
    else
    {
        gr_poly_fit_length(res, len, ctx);
        status = _gr_poly_mulmid_classical(res->coeffs, poly1->coeffs, len1, poly2->coeffs, len2, nlo, nhi, ctx);
    }

    _gr_poly_set_length(res, len, ctx);
    _gr_poly_normalise(res, ctx);
    return status;
}


int
gr_poly_mullow_classical(gr_poly_t res, const gr_poly_t poly1,
                                            const gr_poly_t poly2,
                                                slong n, gr_ctx_t ctx)
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
        status = _gr_poly_mullow_classical(t->coeffs, poly1->coeffs, poly1->length, poly2->coeffs, poly2->length, n, ctx);
        gr_poly_swap(res, t, ctx);
        gr_poly_clear(t, ctx);
    }
    else
    {
        gr_poly_fit_length(res, n, ctx);
        status = _gr_poly_mullow_classical(res->coeffs, poly1->coeffs, poly1->length, poly2->coeffs, poly2->length, n, ctx);
    }

    _gr_poly_set_length(res, n, ctx);
    _gr_poly_normalise(res, ctx);
    return status;
}
