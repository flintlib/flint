/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "gr_vec.h"
#include "gr_poly.h"

truth_t
gr_poly_is_squarefree(const gr_poly_t f, gr_ctx_t ctx)
{
    gr_poly_t d, g;
    truth_t res;
    int status = GR_SUCCESS;

    if (f->length == 0)
        return T_FALSE;

    if (gr_is_zero(gr_poly_coeff_srcptr(f, f->length - 1, ctx), ctx) != T_FALSE)
        return T_UNKNOWN;

    if (f->length <= 2)
        return T_TRUE;

    /* Over a field, f is squarefree iff gcd(f, f') = 1, except in
       positive characteristic where f' = 0 is possible: then f is a p-th
       power over a perfect field (in particular over a finite field). */
    if (gr_ctx_is_field(ctx) != T_TRUE)
        return T_UNKNOWN;

    gr_poly_init(d, ctx);
    gr_poly_init(g, ctx);

    status |= gr_poly_derivative(d, f, ctx);

    if (status != GR_SUCCESS || (d->length > 0 &&
            gr_is_zero(gr_poly_coeff_srcptr(d, d->length - 1, ctx), ctx) != T_FALSE))
    {
        res = T_UNKNOWN;
    }
    else if (d->length == 0)
    {
        if (gr_ctx_is_finite_characteristic(ctx) == T_FALSE)
            res = T_UNKNOWN;  /* should not happen */
        else if (gr_ctx_is_finite(ctx) == T_TRUE)
            res = T_FALSE;
        else
            res = T_UNKNOWN;  /* possibly imperfect field */
    }
    else
    {
        status |= gr_poly_gcd(g, f, d, ctx);

        if (status != GR_SUCCESS)
            res = T_UNKNOWN;
        else if (g->length == 1)
            res = T_TRUE;
        else if (g->length > 1 && gr_is_zero(gr_poly_coeff_srcptr(g, g->length - 1, ctx), ctx) == T_FALSE)
            res = T_FALSE;
        else
            res = T_UNKNOWN;
    }

    gr_poly_clear(d, ctx);
    gr_poly_clear(g, ctx);

    return res;
}
