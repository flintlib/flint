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
gr_ore_poly_differential_to_shift(gr_ore_poly_t res, slong * p, const gr_ore_poly_t op,
                                  gr_ore_poly_ctx_t res_ctx, gr_ore_poly_ctx_t op_ctx)
{
    gr_ctx_struct * base = GR_ORE_POLY_ELEM_CTX(op_ctx);
    ore_algebra_t sa = GR_ORE_POLY_CTX(op_ctx)->which_algebra;
    ore_algebra_t da = GR_ORE_POLY_CTX(res_ctx)->which_algebra;
    slong var = GR_ORE_POLY_ORE_DATA(op_ctx)->base_var;
    slong bsz = base->sizeof_elem;
    int status = GR_SUCCESS;
    slong len, rlen = 0, a = 0, sp = 0, i;
    gr_srcptr eul;
    gr_ptr cur = NULL, rec = NULL;

    *p = 0;

    if (GR_ORE_POLY_ELEM_CTX(res_ctx) != base || base->which_ring != GR_CTX_GR_POLY)
        return GR_UNABLE;

    len = op->length;
    if (len == 0)
        return gr_ore_poly_zero(res, res_ctx);

    /* differential -> Euler (the Euler operator is op itself when sa is already
       the Euler derivative) */
    if (sa == ORE_ALGEBRA_DERIVATIVE)
    {
        GR_TMP_INIT_VEC(cur, len, base);
        status |= _gr_ore_poly_ddx_to_euler(cur, op->coeffs, len, var, base);
        eul = cur;
        a = len - 1;
    }
    else if (sa == ORE_ALGEBRA_EULER_DERIVATIVE)
    {
        eul = op->coeffs;
        a = 0;
    }
    else
    {
        status = GR_DOMAIN;
        goto cleanup;
    }

    for (i = 0; i < len; i++)
    {
        const gr_poly_struct * q = (const gr_poly_struct *) GR_ENTRY(eul, i, bsz);
        if (q->length > rlen)
            rlen = q->length;
    }
    if (rlen == 0)
    {
        status |= gr_ore_poly_zero(res, res_ctx);
        goto cleanup;
    }

    /* Euler -> backward shift, then backward shift -> destination shift type */
    GR_TMP_INIT_VEC(rec, rlen, base);
    status |= _gr_ore_poly_euler_to_backshift_univar(rec, rlen, eul, len, base);

    gr_ore_poly_fit_length(res, rlen, res_ctx);
    status |= _gr_ore_poly_shift_convert(res->coeffs, &sp, rec, rlen,
                                         ORE_ALGEBRA_BACKWARD_SHIFT, da, var, base);
    if (status == GR_SUCCESS)
    {
        _gr_ore_poly_set_length(res, rlen, res_ctx);
        _gr_ore_poly_normalise(res, res_ctx);
        *p = sp + a;
    }

cleanup:
    if (cur != NULL)
        GR_TMP_CLEAR_VEC(cur, len, base);
    if (rec != NULL)
        GR_TMP_CLEAR_VEC(rec, rlen, base);

    return status;
}
