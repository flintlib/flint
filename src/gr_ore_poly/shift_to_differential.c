/*
    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/* generated using Claude Opus 4.8 */

#include "gr.h"
#include "gr_poly.h"
#include "gr_vec.h"
#include "gr_ore_poly.h"

int
gr_ore_poly_shift_to_differential(gr_ore_poly_t res, slong * p, const gr_ore_poly_t op,
                                  gr_ore_poly_ctx_t res_ctx, gr_ore_poly_ctx_t op_ctx)
{
    gr_ctx_struct * base = GR_ORE_POLY_ELEM_CTX(op_ctx);
    ore_algebra_t sa = GR_ORE_POLY_CTX(op_ctx)->which_algebra;
    ore_algebra_t da = GR_ORE_POLY_CTX(res_ctx)->which_algebra;
    slong var = GR_ORE_POLY_ORE_DATA(res_ctx)->base_var;
    slong bsz = base->sizeof_elem;
    int status = GR_SUCCESS;
    slong len, elen = 0, sp = 0, k;
    gr_ptr rec = NULL, eul = NULL;

    *p = 0;

    if (GR_ORE_POLY_ELEM_CTX(res_ctx) != base || base->which_ring != GR_CTX_GR_POLY)
        return GR_UNABLE;

    if (da != ORE_ALGEBRA_DERIVATIVE && da != ORE_ALGEBRA_EULER_DERIVATIVE)
        return GR_DOMAIN;

    len = op->length;
    if (len == 0)
        return gr_ore_poly_zero(res, res_ctx);

    /* source shift type -> backward shift (the algebra is validated here) */
    GR_TMP_INIT_VEC(rec, len, base);
    status |= _gr_ore_poly_shift_convert(rec, &sp, op->coeffs, len,
                                         sa, ORE_ALGEBRA_BACKWARD_SHIFT, var, base);
    if (status != GR_SUCCESS)
        goto cleanup;

    for (k = 0; k < len; k++)
    {
        const gr_poly_struct * q = (const gr_poly_struct *) GR_ENTRY(rec, k, bsz);
        if (q->length > elen)
            elen = q->length;
    }
    if (elen == 0)
    {
        status |= gr_ore_poly_zero(res, res_ctx);
        goto cleanup;
    }

    /* backward shift -> Euler -> destination differential type. For the Euler
       destination, backshift_to_euler_univar writes the result in place. */
    gr_ore_poly_fit_length(res, elen, res_ctx);
    if (da == ORE_ALGEBRA_EULER_DERIVATIVE)
    {
        status |= _gr_ore_poly_backshift_to_euler_univar(res->coeffs, elen, rec, len, base);
    }
    else
    {
        GR_TMP_INIT_VEC(eul, elen, base);
        status |= _gr_ore_poly_backshift_to_euler_univar(eul, elen, rec, len, base);
        status |= _gr_ore_poly_euler_to_ddx(res->coeffs, eul, elen, var, base);
    }
    _gr_ore_poly_set_length(res, elen, res_ctx);
    _gr_ore_poly_normalise(res, res_ctx);
    *p = -sp;

cleanup:
    if (rec != NULL)
        GR_TMP_CLEAR_VEC(rec, len, base);
    if (eul != NULL)
        GR_TMP_CLEAR_VEC(eul, elen, base);

    return status;
}
