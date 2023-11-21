/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ca_mat.h"
#include "gr_mat.h"

void
_ca_mat_charpoly_berkowitz(ca_ptr cp, const ca_mat_t mat, ca_ctx_t ctx)
{
    const slong n = mat->r;

    if (n == 0)
    {
        ca_one(cp, ctx);
    }
    else if (n == 1)
    {
        ca_neg(cp + 0, ca_mat_entry(mat, 0, 0), ctx);
        ca_one(cp + 1, ctx);
    }
    else if (n == 2)
    {
        ca_mat_det_cofactor(cp, mat, ctx);
        ca_add(cp + 1, ca_mat_entry(mat, 0, 0), ca_mat_entry(mat, 1, 1), ctx);
        ca_neg(cp + 1, cp + 1, ctx);
        ca_one(cp + 2, ctx);
    }
    else
    {
        gr_ctx_t gr_ctx;

        _gr_ctx_init_ca_from_ref(gr_ctx, GR_CTX_CC_CA, ctx);
        GR_MUST_SUCCEED(_gr_mat_charpoly_berkowitz(cp, (const gr_mat_struct *) mat, gr_ctx));
    }
}

void ca_mat_charpoly_berkowitz(ca_poly_t cp, const ca_mat_t mat, ca_ctx_t ctx)
{
    ca_poly_fit_length(cp, mat->r + 1, ctx);
    _ca_poly_set_length(cp, mat->r + 1, ctx);
    _ca_mat_charpoly_berkowitz(cp->coeffs, mat, ctx);
}
