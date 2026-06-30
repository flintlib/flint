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

/* sum_k op[k] theta^k expanded via theta^k = sum_i S(k,i) x^i D^i with S(k,i)
   the Stirling numbers of the second kind. Here x is the generator of index
   var. */
int
_gr_ore_poly_euler_to_ddx(gr_ptr res, gr_srcptr op, slong len, slong var, gr_ctx_t ctx)
{
    int status = GR_SUCCESS;
    slong sz = ctx->sizeof_elem;
    slong i, k;
    gr_ctx_t zz;
    gr_mat_t S;
    gr_vec_t gens;
    gr_ptr x, xp, t;

    if (len <= 0)
        return GR_SUCCESS;

    gr_ctx_init_fmpz(zz);
    gr_mat_init(S, len, len, zz);
    /* S(k,i) = 0 for i > k, so only the lower triangle of S is used below */
    status |= gr_mat_stirling(S, 2, zz);

    GR_TMP_INIT3(x, xp, t, ctx);

    gr_vec_init(gens, 0, ctx);
    status |= gr_gens(gens, ctx);
    if (var >= 0 && var < gens->length)
        status |= gr_set(x, gr_vec_entry_srcptr(gens, var, ctx), ctx);
    else
        status |= GR_UNABLE;
    gr_vec_clear(gens, ctx);

    status |= gr_one(xp, ctx);

    for (i = 0; i < len; i++)
    {
        status |= gr_zero(t, ctx);
        for (k = i; k < len; k++)
            status |= gr_addmul_fmpz(t, GR_ENTRY(op, k, sz),
                                     gr_mat_entry_srcptr(S, k, i, zz), ctx);
        status |= gr_mul(GR_ENTRY(res, i, sz), xp, t, ctx);
        status |= gr_mul(xp, xp, x, ctx);
    }

    GR_TMP_CLEAR3(x, xp, t, ctx);
    gr_mat_clear(S, zz);
    gr_ctx_clear(zz);

    return status;
}
