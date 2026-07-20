/*
    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/* generated using Claude Opus 4.8 */

#include "flint.h"
#include "gr.h"
#include "gr_ore_poly.h"

int
gr_ore_poly_apply_custom(gr_ptr res, const gr_ore_poly_t P, gr_srcptr f, gr_srcptr d1, gr_ore_poly_ctx_t ctx)
{
    gr_ctx_struct * base = GR_ORE_POLY_ELEM_CTX(ctx);
    slong el = base->sizeof_elem;
    slong len = P->length;
    int status = GR_SUCCESS;
    truth_t d1_is_zero, d1_is_one;
    gr_ptr g, acc, sig, del, term;

    if (len == 0)
        return gr_zero(res, base);

    d1_is_zero = gr_is_zero(d1, base);
    d1_is_one = gr_is_one(d1, base);

    GR_TMP_INIT5(g, acc, sig, del, term, base);

    status |= gr_set(g, f, base);
    status |= gr_zero(acc, base);  /* todo: add underscore version? */

    for (slong i = 0; i < len; i++)
    {
        /* acc += p_i * (D^i f) */
        status |= gr_mul(term, GR_ENTRY(P->coeffs, i, el), g, base);
        status |= gr_add(acc, acc, term, base);

        /* g <- D(g) = sigma(g)*d1 + delta(g) for the next iteration */
        if (i + 1 < len)
        {
            if (d1_is_zero == T_TRUE)
                status |= gr_ore_poly_delta(g, g, ctx);
            else
            {
                status |= gr_ore_poly_sigma_delta(sig, del, g, ctx);
                if (d1_is_one != T_TRUE)
                    status |= gr_mul(sig, sig, d1, base);
                status |= gr_add(g, sig, del, base);
            }
        }
    }

    status |= gr_set(res, acc, base);

    GR_TMP_CLEAR5(g, acc, sig, del, term, base);

    return status;
}

int
gr_ore_poly_apply(gr_ptr res, const gr_ore_poly_t P, gr_srcptr f, gr_ore_poly_ctx_t ctx)
{
    gr_ctx_struct * base = GR_ORE_POLY_ELEM_CTX(ctx);
    int status = GR_SUCCESS;
    gr_ptr d1;

    GR_TMP_INIT(d1, base);

    switch (GR_ORE_POLY_CTX(ctx)->which_algebra)
    {
        case ORE_ALGEBRA_DERIVATIVE:
        case ORE_ALGEBRA_EULER_DERIVATIVE:
        case ORE_ALGEBRA_FORWARD_DIFFERENCE:
        case ORE_ALGEBRA_BACKWARD_DIFFERENCE:
            status |= gr_zero(d1, base);
            break;
        case ORE_ALGEBRA_FORWARD_SHIFT:
        case ORE_ALGEBRA_BACKWARD_SHIFT:
        case ORE_ALGEBRA_Q_SHIFT:
        case ORE_ALGEBRA_MAHLER:
        case ORE_ALGEBRA_FROBENIUS:
            status |= gr_one(d1, base);
            break;
        default:
            status = GR_UNABLE;
    }

    if (status == GR_SUCCESS)
        status |= gr_ore_poly_apply_custom(res, P, f, d1, ctx);

    GR_TMP_CLEAR(d1, base);

    return status;
}
