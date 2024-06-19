/*
    Copyright (C) 2024 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "mpn_mod.h"
#include "gr_mat.h"

/* todo: retune this when arithmetic is more optimized */

static const unsigned int mpn_mod_mat_lu_delayed_cutoff[MPN_MOD_MAX_LIMBS + 1] =
{
    0, 0, 23, 8, 7, 7, 7, 7, 6, 6, 6, 6, 5, 5, 5, 5, 5,
};

static const unsigned int mpn_mod_mat_lu_recursive_cutoff[MPN_MOD_MAX_LIMBS + 1] =
{
    0, 0, 200, 180, 100, 50, 40, 40, 30, 20, 20, 20, 20, 20, 20, 20, 20,
};

int
mpn_mod_mat_lu(slong * rank, slong * P, gr_mat_t LU, const gr_mat_t A, int rank_check, gr_ctx_t ctx)
{
    slong d, cutoff, nlimbs = MPN_MOD_CTX_NLIMBS(ctx);

    d = FLINT_MIN(A->r, A->c);

    cutoff = mpn_mod_mat_lu_delayed_cutoff[nlimbs];
    if (d < cutoff)
        return gr_mat_lu_classical(rank, P, LU, A, rank_check, ctx);

    cutoff = mpn_mod_mat_lu_recursive_cutoff[nlimbs];
    if (d < cutoff)
        return mpn_mod_mat_lu_classical_delayed(rank, P, LU, A, rank_check, ctx);

    return gr_mat_lu_recursive(rank, P, LU, A, rank_check, ctx);
}
