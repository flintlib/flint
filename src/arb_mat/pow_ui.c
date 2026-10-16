/*
    Copyright (C) 2012 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb.h"
#include "arb_mat.h"
#include "gr.h"
#include "gr_mat.h"

void
arb_mat_pow_ui(arb_mat_t B, const arb_mat_t A, ulong exp, slong prec)
{
    gr_ctx_t ctx;
    gr_ctx_init_real_arb(ctx, prec);
    GR_MUST_SUCCEED(gr_mat_pow_ui((gr_mat_struct *) B, (const gr_mat_struct *) A, exp, ctx));
    gr_ctx_clear(ctx);
}
