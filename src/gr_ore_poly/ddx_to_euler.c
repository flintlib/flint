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
_gr_ore_poly_ddx_to_euler(gr_ptr res, gr_srcptr op, slong len, slong var, gr_ctx_t ctx)
{
    int status = GR_SUCCESS;
    slong sz = ctx->sizeof_elem;
    slong i, k;
    gr_ctx_t ZZ;
    gr_mat_t stirling;
    gr_vec_t gens;
    gr_ptr x, xp, c;

    if (len <= 0)
        return GR_SUCCESS;

    gr_ctx_init_fmpz(ZZ);
    gr_mat_init(stirling, len, len, ZZ);
    status |= gr_mat_stirling(stirling, 1, ZZ);

    GR_TMP_INIT3(x, xp, c, ctx);

    status |= _gr_vec_zero(res, len, ctx);

    /* todo: gr_nth_gen? */
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
        for (k = 0; k <= i; k++)
        {
            status |= gr_mul_fmpz(c, GR_ENTRY(op, i, sz),
                                  gr_mat_entry_srcptr(stirling, i, k, ZZ), ctx);
            status |= gr_addmul(GR_ENTRY(res, k, sz), c, xp, ctx);
        }
        if (i > 0)
            status |= gr_mul(xp, xp, x, ctx);
    }

    GR_TMP_CLEAR3(x, xp, c, ctx);
    gr_mat_clear(stirling, ZZ);
    gr_ctx_clear(ZZ);

    return status;
}
