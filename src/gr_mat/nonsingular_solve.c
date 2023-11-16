/*
    Copyright (C) 2022 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "gr_mat.h"

/* todo: solve_adjugate for n <= 4 */
/* todo: algorithm selection */
int
gr_mat_nonsingular_solve(gr_mat_t X, const gr_mat_t A, const gr_mat_t B, gr_ctx_t ctx)
{
    if (gr_ctx_is_field(ctx) == T_TRUE || gr_ctx_is_exact(ctx) == T_FALSE)
        return gr_mat_nonsingular_solve_lu(X, A, B, ctx);
    else
        return gr_mat_nonsingular_solve_fflu(X, A, B, ctx);
}
