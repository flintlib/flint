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
#include "arb_poly.h"
#include "arb_mat.h"

void
_arb_mat_companion(arb_mat_t A, arb_srcptr poly, slong prec)
{
    int status;
    gr_ctx_t ctx;
    gr_ctx_init_real_arb(ctx, prec);
    status = _gr_mat_companion((gr_mat_struct *) A, (gr_srcptr) poly, ctx);
    gr_ctx_clear(ctx);
    if (status != GR_SUCCESS)
        arb_mat_indeterminate(A);
}

void
arb_mat_companion(arb_mat_t A, const arb_poly_t poly, slong prec)
{
    int status;
    gr_ctx_t ctx;
    gr_ctx_init_real_arb(ctx, prec);
    status = gr_mat_companion((gr_mat_struct *) A, (const gr_poly_struct *) poly, ctx);
    gr_ctx_clear(ctx);
    if (status != GR_SUCCESS)
        arb_mat_indeterminate(A);
}
