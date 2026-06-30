/*
    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/* generated using Claude Opus 4.8 */

#include "gr.h"
#include "gr_vec.h"
#include "gr_ore_poly.h"

int
gr_ore_poly_convert(gr_ore_poly_t res, slong * p, const gr_ore_poly_t op,
                    gr_ore_poly_ctx_t res_ctx, gr_ore_poly_ctx_t op_ctx)
{
    gr_ctx_struct * base = GR_ORE_POLY_ELEM_CTX(op_ctx);
    ore_algebra_t sa = GR_ORE_POLY_CTX(op_ctx)->which_algebra;
    ore_algebra_t da = GR_ORE_POLY_CTX(res_ctx)->which_algebra;
    int status = GR_SUCCESS;
    slong len, sp = 0;

    *p = 0;

    if (GR_ORE_POLY_ELEM_CTX(res_ctx) != base)
        return GR_UNABLE;

    len = op->length;
    if (len == 0)
        return gr_ore_poly_zero(res, res_ctx);

    gr_ore_poly_fit_length(res, len, res_ctx);

    if ((sa == ORE_ALGEBRA_DERIVATIVE || sa == ORE_ALGEBRA_EULER_DERIVATIVE)
        && (da == ORE_ALGEBRA_DERIVATIVE || da == ORE_ALGEBRA_EULER_DERIVATIVE))
    {
        /* differential family {d/dx, theta} */
        slong var = GR_ORE_POLY_ORE_DATA(op_ctx)->base_var;

        if (sa == da)
            status |= _gr_vec_set(res->coeffs, op->coeffs, len, base);
        else if (sa == ORE_ALGEBRA_DERIVATIVE)      /* d/dx -> theta */
        {
            status |= _gr_ore_poly_ddx_to_euler(res->coeffs, op->coeffs, len, var, base);
            /* res = x^(len-1) * op, i.e. x^p * res = op with p = -(len-1) */
            *p = -(len - 1);
        }
        else                                        /* theta -> d/dx */
            status |= _gr_ore_poly_euler_to_ddx(res->coeffs, op->coeffs, len, var, base);
    }
    else
    {
        /* _gr_ore_poly_shift_convert handles the shift/difference family (the
           same-algebra case included) and returns GR_DOMAIN for anything that is
           not a shift/difference type, which is exactly the remaining unsupported
           case: a boundary-crossing pair, or an algebra such as Q_SHIFT, MAHLER
           or COMMUTATIVE (whose identity we deliberately do not claim, since its
           parameters could differ between the two contexts). */
        /* a custom algebra has no ore_data; shift_convert then returns GR_DOMAIN
           regardless of var, so the placeholder 0 is harmless */
        gr_ore_poly_ore_data_t * od = GR_ORE_POLY_ORE_DATA(op_ctx);
        slong var = (od != NULL) ? od->base_var : 0;

        status = _gr_ore_poly_shift_convert(res->coeffs, &sp, op->coeffs, len,
                                            sa, da, var, base);
        if (status == GR_DOMAIN)
            status = GR_UNABLE;
        else
            *p = sp;
    }

    if (status == GR_SUCCESS)
    {
        _gr_ore_poly_set_length(res, len, res_ctx);
        _gr_ore_poly_normalise(res, res_ctx);
    }

    return status;
}
