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

int
_gr_poly_mullow_bivariate_KS(gr_ptr res,
    gr_srcptr poly1, slong len1,
    gr_srcptr poly2, slong len2, slong n, gr_ctx_t ctx)
{
    slong i, l, max_len1, max_len2, inner_len;
    gr_ctx_struct * cctx;
    gr_ctx_t tmp_ctx;
    const gr_poly_struct * ppoly1, * ppoly2;
    gr_poly_struct * pres;
    gr_ptr R, P1, P2, pp, tmp_zero;
    slong csz;
    int status = GR_SUCCESS;
    int squaring;

    /* Todo: other base types */
    if (ctx->which_ring == GR_CTX_GR_POLY)
    {
        cctx = POLYNOMIAL_ELEM_CTX(ctx);
    }
    else if (ctx->which_ring == GR_CTX_FMPZ_POLY)
    {
        gr_ctx_init_fmpz(tmp_ctx);
        cctx = tmp_ctx;
    }
    else
        return GR_UNABLE;

    csz = cctx->sizeof_elem;

    len1 = FLINT_MIN(len1, n);
    len2 = FLINT_MIN(len2, n);

    squaring = (poly1 == poly2) && (len1 == len2);
    ppoly1 = poly1;
    ppoly2 = poly2;
    pres = res;

    max_len1 = 0;
    max_len2 = 0;

    for (i = 0; i < len1; i++)
    {
        l = ppoly1[i].length;
        max_len1 = FLINT_MAX(l, max_len1);
    }

    if (squaring)
    {
        max_len2 = max_len1;
    }
    else
    {
        for (i = 0; i < len2; i++)
        {
            l = ppoly2[i].length;
            max_len2 = FLINT_MAX(l, max_len2);
        }
    }

    if (max_len1 == 0 || max_len2 == 0)
        return _gr_vec_zero(res, n, ctx);

    inner_len = max_len1 + max_len2 - 1;

    GR_TMP_INIT_VEC(R, n * inner_len, cctx);
    GR_TMP_INIT_VEC(tmp_zero, inner_len, cctx);
    P1 = GR_TMP_ALLOC(len1 * inner_len * cctx->sizeof_elem);

    if (squaring)
        P2 = P1;
    else
        P2 = GR_TMP_ALLOC(len2 * inner_len * cctx->sizeof_elem);

    for (i = 0; i < len1; i++)
    {
        l = ppoly1[i].length;
        memcpy(GR_ENTRY(P1, i * inner_len, csz), ppoly1[i].coeffs, l * csz);
        memcpy(GR_ENTRY(P1, i * inner_len + l, csz), tmp_zero, (inner_len - l) * csz);
    }

    if (!squaring)
    {
        for (i = 0; i < len2; i++)
        {
            l = ppoly2[i].length;
            memcpy(GR_ENTRY(P2, i * inner_len, csz), ppoly2[i].coeffs, l * csz);
            memcpy(GR_ENTRY(P2, i * inner_len + l, csz), tmp_zero, (inner_len - l) * csz);
        }
    }

    status = _gr_poly_mullow(R, P1, len1 * inner_len, P2, len2 * inner_len, n * inner_len, cctx);

    for (i = 0; i < n; i++)
    {
        l = inner_len;
        pp = GR_ENTRY(R, i * inner_len, csz);
        while (l > 0 && gr_is_zero(GR_ENTRY(pp, l - 1, csz), cctx) == T_TRUE)
            l--;
        gr_poly_fit_length(pres + i, l, cctx);
        _gr_poly_set_length(pres + i, l, cctx);
        _gr_vec_swap(pres[i].coeffs, pp, l, cctx);
    }

    GR_TMP_CLEAR_VEC(R, n * inner_len, cctx);
    GR_TMP_CLEAR_VEC(tmp_zero, inner_len, cctx);
    GR_TMP_FREE(P1, len1 * inner_len * cctx->sizeof_elem);
    if (!squaring)
        GR_TMP_FREE(P2, len2 * inner_len * cctx->sizeof_elem);

    return status;
}

int
gr_poly_mullow_bivariate_KS(gr_poly_t res, const gr_poly_t poly1,
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
        status = _gr_poly_mullow_bivariate_KS(t->coeffs, poly1->coeffs, poly1->length, poly2->coeffs, poly2->length, n, ctx);
        gr_poly_swap(res, t, ctx);
        gr_poly_clear(t, ctx);
    }
    else
    {
        gr_poly_fit_length(res, n, ctx);
        status = _gr_poly_mullow_bivariate_KS(res->coeffs, poly1->coeffs, poly1->length, poly2->coeffs, poly2->length, n, ctx);
    }

    _gr_poly_set_length(res, n, ctx);
    _gr_poly_normalise(res, ctx);
    return status;
}

