/*
    Copyright (C) 2024 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "mpn_mod.h"

/* todo: retune this when arithmetic is more optimized */
/* note: mul_strassen can be very slightly faster than classical and multi_mod,
   but the range is so narrow that we don't bother */

/* note: same cutoff as multi_mod for n = 2 means we don't use it */
static const slong mpn_mod_mat_mul_waksman_cutoff[MPN_MOD_MAX_LIMBS + 1] =
{
    0, 0, 300, 20, 10, 8, 8, 8, 6, 5, 5, 4, 4, 4, 4, 4, 4,
};

static const int mpn_mod_mat_mul_multi_mod_cutoff[MPN_MOD_MAX_LIMBS + 1] =
{
    0, 0, 300, 80, 80, 90, 100, 90, 90, 90, 80, 80, 70, 70, 70, 70, 70,
};

int
mpn_mod_mat_mul(gr_mat_t C, const gr_mat_t A, const gr_mat_t B, gr_ctx_t ctx)
{
    slong ar = A->r;

    if (ar <= 3)
        return gr_mat_mul_classical(C, A, B, ctx);

    slong ac = A->c;
    slong bc = B->c;
    slong cutoff;
    slong n = MPN_MOD_CTX_NLIMBS(ctx);

    cutoff = mpn_mod_mat_mul_waksman_cutoff[n];

    if (ar < cutoff || ac < cutoff || bc < cutoff)
        return gr_mat_mul_classical(C, A, B, ctx);

    cutoff = mpn_mod_mat_mul_multi_mod_cutoff[n];

    if (ar < cutoff || ac < cutoff || bc < cutoff)
        return mpn_mod_mat_mul_waksman(C, A, B, ctx);

    return mpn_mod_mat_mul_multi_mod(C, A, B, ctx);
}
