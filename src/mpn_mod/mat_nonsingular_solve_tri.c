/*
    Copyright (C) 2024 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "mpn_mod.h"

static const int mpn_mod_mat_solve_tri_cutoff[MPN_MOD_MAX_LIMBS + 1] =
{
    0, 0, 700, 200, 160, 80, 80, 40, 30, 30, 30, 30, 20, 20, 20, 20, 20,
};

int
mpn_mod_mat_nonsingular_solve_tril(gr_mat_t X, const gr_mat_t L, const gr_mat_t B, int unit, gr_ctx_t ctx)
{
    slong cutoff = mpn_mod_mat_solve_tri_cutoff[MPN_MOD_CTX_NLIMBS(ctx)];

    if (B->r < cutoff || B->c < cutoff)
        return gr_mat_nonsingular_solve_tril_classical(X, L, B, unit, ctx);
    else
        return gr_mat_nonsingular_solve_tril_recursive(X, L, B, unit, ctx);
}

int
mpn_mod_mat_nonsingular_solve_triu(gr_mat_t X, const gr_mat_t U, const gr_mat_t B, int unit, gr_ctx_t ctx)
{
    slong cutoff = mpn_mod_mat_solve_tri_cutoff[MPN_MOD_CTX_NLIMBS(ctx)];

    if (B->r < cutoff || B->c < cutoff)
        return gr_mat_nonsingular_solve_triu_classical(X, U, B, unit, ctx);
    else
        return gr_mat_nonsingular_solve_triu_recursive(X, U, B, unit, ctx);
}
