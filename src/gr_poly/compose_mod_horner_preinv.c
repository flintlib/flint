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
    GR_TMP_INIT_VEC(t, len, ctx);

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

    GR_TMP_CLEAR_VEC(t, len, ctx);

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
    return gr_poly_compose_mod_preinv_wrapper(_gr_poly_compose_mod_horner_preinv, res, poly1, poly2, poly3, poly3inv, ctx);
}
