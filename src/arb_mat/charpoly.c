/*
    Copyright (C) 2012 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "gr.h"
#include "gr_mat.h"
#include "arb_poly.h"
#include "arb_mat.h"

/* todo: other algorithms */
void _arb_mat_charpoly(arb_ptr cp, const arb_mat_t mat, slong prec)
{
    if (!arb_mat_is_finite(mat))
    {
        _arb_vec_indeterminate(cp, mat->r + 1);
    }
    else
    {
        gr_ctx_t ctx;
        gr_ctx_init_real_arb(ctx, prec);
        if (_gr_mat_charpoly_berkowitz(cp, (const gr_mat_struct *) mat, ctx) != GR_SUCCESS)
            _arb_vec_indeterminate(cp, mat->r + 1);
    }
}

void arb_mat_charpoly(arb_poly_t cp, const arb_mat_t mat, slong prec)
{
    if (mat->r != mat->c)
    {
        flint_throw(FLINT_ERROR, "Exception (arb_mat_charpoly).  Non-square matrix.\n");
    }

    arb_poly_fit_length(cp, mat->r + 1);
    _arb_poly_set_length(cp, mat->r + 1);
    _arb_mat_charpoly(cp->coeffs, mat, prec);
}
