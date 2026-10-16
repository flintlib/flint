/*
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_poly_mat.h"
#include "gr.h"
#include "gr_mat.h"

void
fmpz_poly_mat_trace(fmpz_poly_t trace, const fmpz_poly_mat_t mat)
{
    gr_ctx_t ctx;
    gr_ctx_init_fmpz_poly(ctx);
    GR_MUST_SUCCEED(gr_mat_trace(trace, (const gr_mat_struct *) mat, ctx));
}
