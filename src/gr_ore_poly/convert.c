/*
    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/* generated using Claude Opus 4.8 */

#include "gr_vec.h"
#include "gr_ore_poly.h"

static int
is_diff_case(ore_algebra_t a)
{
    return a == ORE_ALGEBRA_DERIVATIVE || a == ORE_ALGEBRA_EULER_DERIVATIVE;
}

static int
is_shift_case(ore_algebra_t a)
{
    return a == ORE_ALGEBRA_FORWARD_SHIFT || a == ORE_ALGEBRA_BACKWARD_SHIFT
        || a == ORE_ALGEBRA_FORWARD_DIFFERENCE || a == ORE_ALGEBRA_BACKWARD_DIFFERENCE;
}

int
gr_ore_poly_convert(gr_ore_poly_t res, slong * p, const gr_ore_poly_t op,
                    gr_ore_poly_ctx_t res_ctx, gr_ore_poly_ctx_t op_ctx)
{
    gr_ctx_struct * base = GR_ORE_POLY_ELEM_CTX(op_ctx);
    ore_algebra_t sa = GR_ORE_POLY_CTX(op_ctx)->which_algebra;
    ore_algebra_t da = GR_ORE_POLY_CTX(res_ctx)->which_algebra;
    int status = GR_SUCCESS;
    slong len;

    *p = 0;

    if (GR_ORE_POLY_ELEM_CTX(res_ctx) != base)
        return GR_UNABLE;

    len = op->length;
    if (len == 0)
        return gr_ore_poly_zero(res, res_ctx);

    gr_ore_poly_ore_data_t * od = GR_ORE_POLY_ORE_DATA(op_ctx);
    gr_ore_poly_ore_data_t * rd = GR_ORE_POLY_ORE_DATA(res_ctx);
    slong var = (od != NULL) ? od->base_var : -1;

    gr_ore_poly_fit_length(res, len, res_ctx);

    if (res_ctx == op_ctx)
    {
        status |= _gr_vec_set(res->coeffs, op->coeffs, len, base);
        return GR_SUCCESS;
    }

    if (!(is_diff_case(sa) || is_shift_case(sa))
        || !(is_diff_case(da) || is_shift_case(da)))
        return GR_UNABLE;

    if (od->base_var != rd->base_var)
        return GR_UNABLE;

    if (sa == da)
        status |= _gr_vec_set(res->coeffs, op->coeffs, len, base);
    else if (sa == ORE_ALGEBRA_DERIVATIVE && da == ORE_ALGEBRA_EULER_DERIVATIVE)
    {
        status |= _gr_ore_poly_ddx_to_euler(res->coeffs, op->coeffs, len, var, base);
        *p = -(len - 1);
    }
    else if (sa == ORE_ALGEBRA_EULER_DERIVATIVE && da == ORE_ALGEBRA_DERIVATIVE)
        status |= _gr_ore_poly_euler_to_ddx(res->coeffs, op->coeffs, len, var, base);
    else if (sa == ORE_ALGEBRA_FORWARD_DIFFERENCE
             && da == ORE_ALGEBRA_BACKWARD_DIFFERENCE)
        status |= _gr_ore_poly_shift_convert_difference(res->coeffs, p,
                                                        op->coeffs, len, 1, var,
                                                        base);
    else if (sa == ORE_ALGEBRA_BACKWARD_DIFFERENCE
             && da == ORE_ALGEBRA_FORWARD_DIFFERENCE)
        status |= _gr_ore_poly_shift_convert_difference(res->coeffs, p,
                                                        op->coeffs, len, 0, var,
                                                        base);
    else if (is_shift_case(sa) && is_shift_case(da))
        status = _gr_ore_poly_shift_convert(res->coeffs, p, op->coeffs, len,
                                            sa, da, var, base);
    else
        return GR_UNABLE;

    if (status == GR_SUCCESS)
    {
        _gr_ore_poly_set_length(res, len, res_ctx);
        _gr_ore_poly_normalise(res, res_ctx);
    }

    return status;
}
