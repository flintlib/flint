/*
    Copyright (C) 2026 Joel Dahne

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "gr_special.h"
#include "gr_vec.h"
#include "gr_poly.h"

int
_gr_poly_bessel_j_series(gr_ptr res, gr_srcptr nu, gr_srcptr z, slong zlen, slong len, gr_ctx_t ctx)
{
    gr_ptr t, u;
    slong sz = ctx->sizeof_elem;
    int status = GR_SUCCESS;

    if (zlen == 0)
    {
        /* z = 0, constant term is bessel_j(nu, 0) rest is zero */
        status |= _gr_vec_zero(res, len, ctx);
        if (len > 0)
            status |= gr_bessel_j(GR_ENTRY(res, 0, sz), nu, GR_ENTRY(res, 0, sz), ctx);
        return status;
    }

    t = gr_heap_init_vec(len, ctx);

    /* Compute jet for constant part */
    status |= gr_bessel_j_jet(t, nu, GR_ENTRY(z, 0, sz), len, ctx);

    if (len < 2 || zlen == 1)
    {
        _gr_vec_swap(res, t, len, ctx);
        gr_heap_clear_vec(t, len, ctx);
        return status;
    }

    u = gr_heap_init_vec(zlen, ctx);

    /* Compose with nonconstant part */
    status |= _gr_vec_set(GR_ENTRY(u, 1, sz), GR_ENTRY(z, 1, sz), zlen - 1, ctx);
    status |=  _gr_poly_compose_series(res, t, len, u, zlen, len, ctx);

    gr_heap_clear_vec(t, len, ctx);
    gr_heap_clear_vec(u, zlen, ctx);

    return status;
}

int
gr_poly_bessel_j_series(gr_poly_t res, gr_srcptr nu, const gr_poly_t z, slong len, gr_ctx_t ctx)
{
    int status = GR_SUCCESS;

    if (len == 0)
        return gr_poly_zero(res, ctx);

    if (z->length < 2)
        len = 1;

    gr_poly_fit_length(res, len, ctx);
    status |= _gr_poly_bessel_j_series(res->coeffs, nu, z->coeffs, z->length, len, ctx);
    _gr_poly_set_length(res, len, ctx);
    _gr_poly_normalise(res, ctx);
    return status;
}
