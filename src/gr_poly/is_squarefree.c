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

/* f of length len >= 3 with nonzero leading coefficient, over a field */
static truth_t
_gr_poly_is_squarefree_field(gr_srcptr f, slong len, gr_ctx_t ctx)
{
    gr_ptr d, g;
    slong dlen, glen, sz = ctx->sizeof_elem;
    truth_t res;
    int status = GR_SUCCESS;

    /* derivative and gcd share one temporary vector */
    GR_TMP_INIT_VEC(d, 2 * (len - 1), ctx);
    g = GR_ENTRY(d, len - 1, sz);

    status |= _gr_poly_derivative(d, f, len, ctx);
    dlen = len - 1;
    status |= _gr_vec_normalise(&dlen, d, dlen, ctx);

    if (status != GR_SUCCESS || (dlen > 0 &&
            gr_is_zero(GR_ENTRY(d, dlen - 1, sz), ctx) != T_FALSE))
    {
        res = T_UNKNOWN;
    }
    else if (dlen == 0)
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
        /* the gcd need not be made monic: only its length is used
           (both inputs have nonzero leading coefficients) */
        status |= _gr_poly_gcd(g, &glen, f, len, d, dlen, ctx);

        if (status != GR_SUCCESS)
            res = T_UNKNOWN;
        else if (glen == 1)
            res = T_TRUE;
        else if (glen > 1 && gr_is_zero(GR_ENTRY(g, glen - 1, sz), ctx) == T_FALSE)
            res = T_FALSE;
        else
            res = T_UNKNOWN;
    }

    GR_TMP_CLEAR_VEC(d, 2 * (len - 1), ctx);

    return res;
}

truth_t
gr_poly_is_squarefree(const gr_poly_t f, gr_ctx_t ctx)
{
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

    return _gr_poly_is_squarefree_field(f->coeffs, f->length, ctx);
}
