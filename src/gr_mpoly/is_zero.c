/*
    Copyright (C) 2020 Daniel Schultz
    Copyright (C) 2025 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "mpoly.h"
#include "gr_mpoly.h"

truth_t gr_mpoly_is_zero(const gr_mpoly_t A, gr_mpoly_ctx_t ctx)
{
    if (A->length == 0)
        return T_TRUE;

    if (gr_ctx_is_canonical(GR_MPOLY_CCTX(ctx)) == T_TRUE)
        return T_FALSE;
    else
        return _gr_vec_is_zero(A->coeffs, A->length, GR_MPOLY_CCTX(ctx));
}
