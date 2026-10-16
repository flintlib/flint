/*
    Copyright (C) 2010,2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mat.h"
#include "gr.h"
#include "gr_mat.h"

int
fmpz_mat_randpermdiag(fmpz_mat_t mat, flint_rand_t state,
                      const fmpz * diag, slong n)
{
    int parity;
    gr_ctx_t ctx;
    gr_ctx_init_fmpz(ctx);
    GR_MUST_SUCCEED(gr_mat_randpermdiag(&parity, (gr_mat_struct *) mat, state, (gr_ptr) diag, n, ctx));
    return parity;
}
