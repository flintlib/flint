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
#include "acb_poly.h"
#include "acb_mat.h"

void _acb_mat_charpoly(acb_ptr cp, const acb_mat_t mat, slong prec)
{
    if (!acb_mat_is_finite(mat))
    {
        _acb_vec_indeterminate(cp, mat->r + 1);
    }
    else
    {
        gr_ctx_t ctx;
        gr_ctx_init_complex_acb(ctx, prec);
        if (_gr_mat_charpoly_berkowitz(cp, (const gr_mat_struct *) mat, ctx) != GR_SUCCESS)
            _acb_vec_indeterminate(cp, mat->r + 1);
    }
}

void acb_mat_charpoly(acb_poly_t cp, const acb_mat_t mat, slong prec)
{
    if (mat->r != mat->c)
    {
        flint_throw(FLINT_ERROR, "Exception (acb_mat_charpoly).  Non-square matrix.\n");
    }

    acb_poly_fit_length(cp, mat->r + 1);
    _acb_poly_set_length(cp, mat->r + 1);
    _acb_mat_charpoly(cp->coeffs, mat, prec);
}
