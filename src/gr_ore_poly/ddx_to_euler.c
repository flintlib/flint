/*
    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/* generated using Claude Opus 4.8 */

#include "fmpz.h"
#include "gr.h"
#include "gr_mat.h"
#include "gr_vec.h"
#include "gr_ore_poly.h"

/* x^(len-1) * sum_i op[i] D^i, expanded via x^i D^i = sum_k s(i,k) theta^k with
   s(i,k) the signed Stirling numbers of the first kind. Here x is the generator
   of index var. */
int
_gr_ore_poly_ddx_to_euler(gr_ptr res, gr_srcptr op, slong len, slong var, gr_ctx_t ctx)
{
    int status = GR_SUCCESS;
    slong sz = ctx->sizeof_elem;
    slong i, k;
    gr_ctx_t zz;
    gr_mat_t S;
    gr_vec_t gens;
    gr_ptr x, xp, c;

    if (len <= 0)
        return GR_SUCCESS;

    gr_ctx_init_fmpz(zz);
    gr_mat_init(S, len, len, zz);
    /* s(i,k) = 0 for k > i, so only the lower triangle of S is used below */
    status |= gr_mat_stirling(S, 1, zz);

    GR_TMP_INIT3(x, xp, c, ctx);

    status |= _gr_vec_zero(res, len, ctx);

    gr_vec_init(gens, 0, ctx);
    status |= gr_gens(gens, ctx);
    if (var >= 0 && var < gens->length)
        status |= gr_set(x, gr_vec_entry_srcptr(gens, var, ctx), ctx);
    else
        status |= GR_UNABLE;
    gr_vec_clear(gens, ctx);

    status |= gr_one(xp, ctx);

    for (i = len - 1; i >= 0; i--)
    {
        status |= gr_mul(c, xp, GR_ENTRY(op, i, sz), ctx);
        for (k = 0; k <= i; k++)
            status |= gr_addmul_fmpz(GR_ENTRY(res, k, sz), c,
                                     gr_mat_entry_srcptr(S, i, k, zz), ctx);
        if (i > 0)
            status |= gr_mul(xp, xp, x, ctx);
    }

    GR_TMP_CLEAR3(x, xp, c, ctx);
    gr_mat_clear(S, zz);
    gr_ctx_clear(zz);

    return status;
}
