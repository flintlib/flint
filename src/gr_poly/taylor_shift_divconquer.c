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

/* todo: fmpz_poly has a better algorithm (assuming that generating
   binomial coefficients is fast) */
int
_gr_poly_taylor_shift_divconquer(gr_ptr res, gr_srcptr poly, slong len, gr_srcptr c, gr_ctx_t ctx)
{
    int status;
    slong sz = ctx->sizeof_elem;

    status = GR_SUCCESS;

    if (res != poly)
        status |= _gr_vec_set(res, poly, len, ctx);

    if (len <= 1 || gr_is_zero(c, ctx) == T_TRUE)
        return status;

    if (len == 2)
        return gr_addmul(res, GR_ENTRY(res, 1, sz), c, ctx);

    {
        gr_ptr t;

        GR_TMP_INIT_VEC(t, 2, ctx);

        status |= gr_set(t, c, ctx);
        status |= gr_one(GR_ENTRY(t, 1, sz), ctx);
        status |= _gr_poly_compose_divconquer(res, res, len, t, 2, ctx);

        GR_TMP_CLEAR_VEC(t, 2, ctx);

        return status;
    }
}

int
gr_poly_taylor_shift_divconquer(gr_poly_t res, const gr_poly_t f, gr_srcptr c, gr_ctx_t ctx)
{
    int status = GR_SUCCESS;

    if (res != f)
        status |= gr_poly_set(res, f, ctx);

    status |= _gr_poly_taylor_shift_divconquer(res->coeffs, res->coeffs, res->length, c, ctx);
    return status;
}
