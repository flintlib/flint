/*
    Copyright (C) 2011, 2025 Fredrik Johansson
    Copyright (C) 2012 Lina Kulakova
    Copyright (C) 2013 Mike Hansen

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "gr_vec.h"
#include "gr_poly.h"

int
_gr_poly_compose_mod_horner_preinv(gr_ptr res,
    gr_srcptr f, slong lenf,
    gr_srcptr g,
    gr_srcptr h, slong lenh,
    gr_srcptr hinv, slong lenhinv,
    gr_ctx_t ctx)
{
    slong i, len;
    gr_ptr t;
    int status = GR_SUCCESS;
    slong sz = ctx->sizeof_elem;

    if (lenh == 1)
        return status;

    if (lenf == 1)
        return gr_set(res, f, ctx);

    if (lenh == 2)
        return _gr_poly_evaluate(res, f, lenf, g, ctx);

    len = lenh - 1;
    i = lenf - 1;
    GR_TMP_INIT_VEC(t, 2 * lenh - 3, ctx);

    status |= _gr_vec_mul_scalar(res, g, len, GR_ENTRY(f, i, sz), ctx);
    i--;
    if (i >= 0)
        status |= gr_add(res, res, GR_ENTRY(f, i, sz), ctx);

    while (i > 0)
    {
        i--;
        status |= _gr_poly_mulmod_preinv(t, res, len, g, len, h, lenh, hinv, lenhinv, ctx);
        status |= _gr_poly_add(res, t, len, GR_ENTRY(f, i, sz), 1, ctx);
    }

    GR_TMP_CLEAR_VEC(t, 2 * lenh - 3, ctx);

    return status;
}

int
gr_poly_compose_mod_horner_preinv(gr_poly_t res,
                                      const gr_poly_t poly1,
                                      const gr_poly_t poly2,
                                      const gr_poly_t poly3,
                                      const gr_poly_t poly3inv,
                                      gr_ctx_t ctx)
{
    slong len1 = poly1->length;
    slong len2 = poly2->length;
    slong len3 = poly3->length;
    slong len3inv = poly3inv->length;
    slong len = len3 - 1;
    slong vec_len = FLINT_MAX(len3 - 1, len2);
    gr_ptr ptr2;
    slong sz = ctx->sizeof_elem;
    int status = GR_SUCCESS;

    if (len3 == 0)
        return GR_DOMAIN;

    if (len1 == 0 || len3 == 1)
        return gr_poly_zero(res, ctx);

    if (len1 == 1)
        return gr_poly_set(res, poly1, ctx);

    if (res == poly3 || res == poly1)
    {
        gr_poly_t tmp;
        gr_poly_init(tmp, ctx);
        status = gr_poly_compose_mod_horner_preinv(tmp, poly1, poly2, poly3, poly3inv, ctx);
        gr_poly_swap(tmp, res, ctx);
        gr_poly_clear(tmp, ctx);
        return status;
    }

    GR_TMP_INIT_VEC(ptr2, vec_len, ctx);

    if (len2 <= len3 - 1)
    {
        status |= _gr_vec_set(ptr2, poly2->coeffs, len2, ctx);
        status |= _gr_vec_zero(GR_ENTRY(ptr2, len2, sz), vec_len - len2, ctx);
    }
    else
    {
        /* todo: take advantage of preinv */
        status |= _gr_poly_rem(ptr2, poly2->coeffs, len2, poly3->coeffs, len3, ctx);
    }

    gr_poly_fit_length(res, len, ctx);
    status |= _gr_poly_compose_mod_horner_preinv(res->coeffs, poly1->coeffs, len1, ptr2, poly3->coeffs, len3, poly3inv->coeffs, len3inv, ctx);
    _gr_poly_set_length(res, len, ctx);
    _gr_poly_normalise(res, ctx);

    GR_TMP_CLEAR_VEC(ptr2, vec_len, ctx);
    return status;
}
