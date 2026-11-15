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
gr_poly_randtest(gr_poly_t poly, flint_rand_t state, slong len, gr_ctx_t ctx)
{
    int status = GR_SUCCESS;
    gr_poly_fit_length(poly, len, ctx);
    status |= _gr_vec_randtest(poly->coeffs, state, len, ctx);
    _gr_poly_set_length_normalise(poly, len, ctx);
    return status;
}

/* generated using Claude Opus 5.5 */
int
gr_poly_randtest_shiftful(gr_poly_t poly, flint_rand_t state, slong len, gr_ctx_t ctx)
{
    gr_poly_t u;
    gr_ptr a;
    int status = GR_SUCCESS;

    if (len < 3)
        return gr_poly_randtest(poly, state, len, ctx);

    gr_poly_init(u, ctx);
    a = gr_heap_init(ctx);

    status |= gr_poly_one(poly, ctx);

    for (slong iter = 0; iter < 10 && status == GR_SUCCESS; iter++)
    {
        slong room = len - poly->length + 1;

        if (poly->length == 0 || room < 2)
            break;

        status |= gr_poly_randtest(u, state,
                2 + n_randint(state, FLINT_MIN(room, 4) - 1), ctx);
        if (u->length < 2)
            continue;

        slong mult = 1 + n_randint(state, 5);
        for (slong j = 0; j < mult && status == GR_SUCCESS; j++)
        {
            if (poly->length + u->length - 1 > len)
                break;

            if (n_randint(state, 2))
            {
                slong c = n_randint(state, 4) ? n_randint(state, 8)
                                              : n_randint(state, 1 << 20);
                if (n_randint(state, 2))
                    c = -c;
                status |= gr_set_si(a, c, ctx);
                status |= gr_poly_taylor_shift(u, u, a, ctx);
            }

            if (n_randint(state, 4))
                status |= gr_poly_mul(poly, poly, u, ctx);
        }

        if (n_randint(state, 4) == 0)
            break;
    }

    gr_heap_clear(a, ctx);
    gr_poly_clear(u, ctx);

    return status;
}
