/*
    Copyright (C) 2024 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "mpn_mod.h"

int
mpn_mod_mat_det(nn_ptr res, const gr_mat_t A, gr_ctx_t ctx)
{
    slong n = A->r;

    if (A->r != A->c)
        return GR_DOMAIN;

    if (n <= 4)
        return gr_mat_det_cofactor(res, A, ctx);

    if (n == 5)
        return gr_mat_det_berkowitz(res, A, ctx);

    if (gr_mat_det_lu(res, A, ctx) != GR_SUCCESS)
    {
        /* Fall back on division-free algorithm if we encountered an impossible inverse */
        /* Could try something else here: faddeev_bsgs (O(n^3.5)) or Howell form. */
        GR_MUST_SUCCEED(gr_mat_det_berkowitz(res, A, ctx));
    }

    return GR_SUCCESS;
}
