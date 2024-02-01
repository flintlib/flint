/*
    Copyright (C) 2022 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "gr_mat.h"

int
gr_mat_adjugate_charpoly(gr_mat_t adj, gr_ptr det, const gr_mat_t A, gr_ctx_t ctx)
{
    int status = GR_SUCCESS;
    gr_ptr pol;
    slong n;
    slong sz = ctx->sizeof_elem;

    n = gr_mat_nrows(A, ctx);

    if (n != gr_mat_ncols(A, ctx))
        return GR_DOMAIN;

    if (n == 0)
    {
        return gr_one(det, ctx);
    }
    else
    {
        GR_TMP_INIT_VEC(pol, n + 1, ctx);

        status |= _gr_mat_charpoly(pol, A, ctx);

        if (n % 2)
            status |= gr_neg(det, pol, ctx);
        else
            gr_swap(det, pol, ctx);

        /* todo: verify that _gr_mat_gr_poly_evaluate supports aliasing */
        status |= _gr_mat_gr_poly_evaluate(adj, GR_ENTRY(pol, 1, sz), n, A, ctx);

        if (n % 2 == 0)
            status |= gr_mat_neg(adj, adj, ctx);

        GR_TMP_CLEAR_VEC(pol, n + 1, ctx);
    }

    return status;
}
