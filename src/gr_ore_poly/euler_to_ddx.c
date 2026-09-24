/*
    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/* generated using Claude Opus 4.8 */

#include "gr.h"
#include "gr_mat.h"
#include "gr_vec.h"
#include "gr_ore_poly.h"

int
_gr_ore_poly_euler_to_ddx(gr_ptr res, gr_srcptr op, slong len, slong var, gr_ctx_t ctx)
{
    int status = GR_SUCCESS;
    slong sz = ctx->sizeof_elem;
    slong i, k;
    gr_ctx_t ZZ;
    gr_mat_t stirling;
    gr_vec_t gens;
    gr_ptr x, xp, t;

    if (len == 0)
        return GR_SUCCESS;

    gr_ctx_init_fmpz(ZZ);
    gr_mat_init(stirling, len, len, ZZ);
    status |= gr_mat_stirling(stirling, 2, ZZ);

    GR_TMP_INIT3(x, xp, t, ctx);

    /* todo: gr_nth_gen? */
    gr_vec_init(gens, 0, ctx);
    status |= gr_gens(gens, ctx);
    if (var >= 0 && var < gens->length)
        status |= gr_set(x, gr_vec_entry_srcptr(gens, var, ctx), ctx);
    else
        status |= GR_UNABLE;
    gr_vec_clear(gens, ctx);

    status |= gr_one(xp, ctx);

    /* theta^k = sum_i S(k,i) x^i D^i, so res[i] = (sum_k S(k, i) op[k]) x^i */
    for (i = 0; i < len; i++)
    {
        status |= gr_zero(t, ctx);
        for (k = i; k < len; k++)
            status |= gr_addmul_fmpz(t, GR_ENTRY(op, k, sz),
                                     gr_mat_entry_srcptr(stirling, k, i, ZZ),
                                     ctx);
        status |= gr_mul(GR_ENTRY(res, i, sz), xp, t, ctx);
        status |= gr_mul(xp, xp, x, ctx);
    }

    GR_TMP_CLEAR3(x, xp, t, ctx);
    gr_mat_clear(stirling, ZZ);
    gr_ctx_clear(ZZ);

    return status;
}
