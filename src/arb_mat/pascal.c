/*
    Copyright (C) 2018 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb_mat.h"
#include "gr.h"
#include "gr_mat.h"

void
arb_mat_pascal(arb_mat_t mat, int triangular, slong prec)
{
    gr_ctx_t ctx;
    gr_ctx_init_real_arb(ctx, prec);
    GR_MUST_SUCCEED(gr_mat_pascal((gr_mat_struct *) mat, triangular, ctx));
}
