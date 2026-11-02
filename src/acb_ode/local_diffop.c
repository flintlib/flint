/*
    Copyright (C) 2026 Marc Mezzarobba

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb_ode.h"
#include "gr_ore_poly.h"
#include "gr_poly.h"

void
acb_ode_local_diffop_ctx_init(gr_ctx_t Pol, gr_ore_poly_ctx_t Euler,
                              gr_ctx_t Dop, gr_ctx_t Cst)
{
    int status;
    char * gen = NULL;

    status = gr_ore_poly_ctx_over_gr_poly_base_ptrs(NULL, NULL, Dop);
    GR_MUST_SUCCEED(status);

    FLINT_ASSERT(GR_ORE_POLY_ORE_DATA(Dop)->base_var == 0);

    gr_ctx_init_gr_poly(Pol, Cst);
    gr_ore_poly_ctx_init(Euler, Pol, 0, ORE_ALGEBRA_EULER_DERIVATIVE);

    status |= gr_ctx_gen_name(&gen, 0, GR_ORE_POLY_CTX(Dop)->base_ring);
    status |= gr_ctx_set_gen_name(Pol, gen);

    status |= gr_ctx_set_gen_name(Euler, GR_ORE_POLY_CTX(Dop)->var);
    if (GR_ORE_POLY_CTX(Euler)->var != NULL
        && GR_ORE_POLY_CTX(Euler)->var[0] == 'D')
        GR_ORE_POLY_CTX(Euler)->var[0] = 'T';

    GR_MUST_SUCCEED(status);

    flint_free(gen);
}

WARN_UNUSED_RESULT int
_acb_ode_local_diffop(gr_poly_struct * ldop,
                      const gr_poly_struct * dop, slong len,
                      ore_algebra_t src_alg, gr_ctx_t Pol)
{
    int status = GR_SUCCESS;

    FLINT_ASSERT(Pol->which_ring == GR_CTX_GR_POLY);
    gr_ctx_struct * Cst = POLYNOMIAL_ELEM_CTX(Pol);

    if (src_alg == ORE_ALGEBRA_EULER_DERIVATIVE)
        status |= _gr_vec_set(ldop, dop, len, Pol);
    else if (src_alg == ORE_ALGEBRA_DERIVATIVE)
        status |= _gr_ore_poly_ddx_to_euler(ldop, dop, len, 0, Pol);
    else
        return GR_UNABLE;

    slong xval = WORD_MAX;
    for (slong i = 0; xval > 0 && i < len; i++)
    {
        slong v = 0;
        while (v < xval && v < ldop[i].length
               && gr_is_zero(gr_poly_coeff_ptr(ldop + i, v, Cst), Cst) == T_TRUE)
            v++;
        if (v < ldop[i].length)
            xval = v;
    }

    for (slong i = 0; i < len; i++)
        status |= gr_poly_shift_right(ldop + i, ldop + i, xval, Cst);

    return status;
}

WARN_UNUSED_RESULT int
acb_ode_local_diffop(gr_ore_poly_t ldop, const gr_ore_poly_t dop,
                     gr_ctx_t Dop, gr_ore_poly_ctx_t Euler)
{
    gr_ctx_struct * Cst0, * Pol0, * Cst1, * Pol1;
    int status = GR_SUCCESS;

    ore_algebra_t src_alg = GR_ORE_POLY_CTX(Dop)->which_algebra;
    if (GR_ORE_POLY_CTX(Euler)->which_algebra != ORE_ALGEBRA_EULER_DERIVATIVE)
        return GR_UNABLE;

    status |= gr_ore_poly_ctx_over_gr_poly_base_ptrs(&Cst0, &Pol0, Dop);
    status |= gr_ore_poly_ctx_over_gr_poly_base_ptrs(&Cst1, &Pol1, Euler);

    if (dop->length == 0)
        status |= gr_zero(ldop, Euler);

    if (status != GR_SUCCESS)
        return status;

    /* todo: may eventually be subsumed by a more powerful gr_ore_poly
       implementation of gr_set_other */
    gr_ore_poly_fit_length(ldop, dop->length, Euler);
    _gr_ore_poly_set_length(ldop, dop->length, Euler);

    if (Cst0->which_ring == GR_CTX_CC_ACB && Cst1->which_ring == GR_CTX_CC_ACB)
        /* Avoid rounding to the destination precision.
           todo: something similar for RR_ARB. */
        status |= _acb_ode_local_diffop(ldop->coeffs, dop->coeffs, dop->length,
                                        src_alg, Pol0);
    else
    {
        gr_poly_struct * ldop0;
        GR_TMP_INIT_VEC(ldop0, dop->length, Pol0);

        status |= _acb_ode_local_diffop(ldop0, dop->coeffs, dop->length,
                                        src_alg, Pol0);
        status |= _gr_vec_set_other(ldop->coeffs, ldop0, Pol0, dop->length, Pol1);

        GR_TMP_CLEAR_VEC(ldop0, dop->length, Pol0);
    }

    return status;
}
