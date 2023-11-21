/*
    Copyright (C) 2015 William Hart
    Copyright (C) 2020 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ca_mat.h"
#include "gr_mat.h"

int
_ca_mat_charpoly_danilevsky_inplace(ca_ptr p, ca_mat_t A, ca_ctx_t ctx)
{
    int success;
    gr_ctx_t gr_ctx;

    _gr_ctx_init_ca_from_ref(gr_ctx, GR_CTX_CC_CA, ctx);
    success = (_gr_mat_charpoly_danilevsky_inplace(p, (gr_mat_struct *) A, gr_ctx) == GR_SUCCESS);

    return success;
}

int
_ca_mat_charpoly_danilevsky(ca_ptr cp, const ca_mat_t A, ca_ctx_t ctx)
{
    ca_mat_t T;
    int success;

    ca_mat_init(T, ca_mat_nrows(A), ca_mat_nrows(A), ctx);
    ca_mat_set(T, A, ctx);
    success = _ca_mat_charpoly_danilevsky_inplace(cp, T, ctx);
    ca_mat_clear(T, ctx);
    return success;
}

int
ca_mat_charpoly_danilevsky(ca_poly_t cp, const ca_mat_t mat, ca_ctx_t ctx)
{
    ca_poly_fit_length(cp, mat->r + 1, ctx);
    _ca_poly_set_length(cp, mat->r + 1, ctx);
    return _ca_mat_charpoly_danilevsky(cp->coeffs, mat, ctx);
}
