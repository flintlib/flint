/*
    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/* generated using Claude Opus 4.8 */

#include "gr.h"
#include "gr_ore_poly.h"

int
gr_ore_poly_shift_to_differential(gr_ore_poly_t res, slong * p,
                                  const gr_ore_poly_t op,
                                  gr_ore_poly_ctx_t res_ctx,
                                  gr_ore_poly_ctx_t op_ctx)
{
    gr_ore_poly_ctx_t rec_ctx, eul_ctx;
    gr_ore_poly_t rec, eul;

    int status = GR_SUCCESS;

    gr_ctx_struct * base = GR_ORE_POLY_ELEM_CTX(op_ctx);
    if (GR_ORE_POLY_ELEM_CTX(res_ctx) != base || base->which_ring != GR_CTX_GR_POLY)
        return GR_UNABLE;

    ore_algebra_t da = GR_ORE_POLY_CTX(res_ctx)->which_algebra;
    if (da != ORE_ALGEBRA_DERIVATIVE && da != ORE_ALGEBRA_EULER_DERIVATIVE)
        return GR_UNABLE;

    slong len = op->length;
    slong var = GR_ORE_POLY_ORE_DATA(res_ctx)->base_var;

    gr_ore_poly_ctx_init(rec_ctx, base, var, ORE_ALGEBRA_BACKWARD_SHIFT);
    gr_ore_poly_ctx_init(eul_ctx, base, var, ORE_ALGEBRA_EULER_DERIVATIVE);

    gr_ore_poly_init(rec, rec_ctx);
    gr_ore_poly_init(eul, eul_ctx);

    /* this is where the input algebra type is validated */
    status |= gr_ore_poly_convert(rec, p, op, rec_ctx, op_ctx);
    *p = -*p;

    if (len == 0)
    {
        status |= gr_ore_poly_zero(res, res_ctx);
        goto cleanup;
    }

    slong elen = 0;
    for (slong k = 0; k < len; k++)
    {
        const gr_poly_struct * q = GR_ENTRY(rec->coeffs, k, base->sizeof_elem);
        if (q->length > elen)
            elen = q->length;
    }
    if (elen == 0)
    {
        status |= gr_ore_poly_zero(res, res_ctx);
        goto cleanup;
    }

    gr_ore_poly_fit_length(res, elen, res_ctx);
    if (da == ORE_ALGEBRA_EULER_DERIVATIVE)
        status |= _gr_ore_poly_backshift_to_euler_univar(res->coeffs, elen, rec->coeffs, len, base);
    else
    {
        gr_ore_poly_fit_length(eul, elen, res_ctx);
        status |= _gr_ore_poly_backshift_to_euler_univar(eul->coeffs, elen, rec->coeffs, len, base);
        status |= _gr_ore_poly_euler_to_ddx(res->coeffs, eul->coeffs, elen, var, base);
    }
    _gr_ore_poly_set_length(res, elen, res_ctx);
    _gr_ore_poly_normalise(res, res_ctx);

cleanup:

    gr_ore_poly_clear(rec, rec_ctx);
    gr_ore_poly_clear(eul, eul_ctx);
    gr_ore_poly_ctx_clear(rec_ctx);
    gr_ore_poly_ctx_clear(eul_ctx);

    return status;
}
