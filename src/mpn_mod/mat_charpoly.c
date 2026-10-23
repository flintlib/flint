/*
    Copyright (C) 2025 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "mpn_mod.h"
#include "gr.h"
#include "gr_mat.h"

int
_mpn_mod_mat_charpoly(nn_ptr res, const gr_mat_t mat, gr_ctx_t ctx)
{
    slong n = mat->r;

    /* todo: depends weakly on limb count; could be fine-tuned */
    /* Danilevsky needs inverses, so it can fail for a composite modulus;
       Berkowitz is division-free and always succeeds */
    if (n > 11 && _gr_mat_charpoly_danilevsky(res, (const gr_mat_struct *) mat, ctx) == GR_SUCCESS)
        return GR_SUCCESS;

    return _gr_mat_charpoly_berkowitz(res, (const gr_mat_struct *) mat, ctx);
}
