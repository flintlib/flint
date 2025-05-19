/*
    Copyright (C) 2018 Fredrik Johansson

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

void
_acb_mat_companion(acb_mat_t A, acb_srcptr poly, slong prec)
{
    int status;
    gr_ctx_t ctx;
    gr_ctx_init_complex_acb(ctx, prec);
    status = _gr_mat_companion((gr_mat_struct *) A, (gr_srcptr) poly, ctx);
    gr_ctx_clear(ctx);
    if (status != GR_SUCCESS)
        acb_mat_indeterminate(A);
}

void
acb_mat_companion(acb_mat_t A, const acb_poly_t poly, slong prec)
{
    int status;
    gr_ctx_t ctx;
    gr_ctx_init_complex_acb(ctx, prec);
    status = gr_mat_companion((gr_mat_struct *) A, (const gr_poly_struct *) poly, ctx);
    gr_ctx_clear(ctx);
    if (status != GR_SUCCESS)
        acb_mat_indeterminate(A);
}
