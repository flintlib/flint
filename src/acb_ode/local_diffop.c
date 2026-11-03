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

static slong
x_valuation(const gr_ore_poly_t dop, gr_ctx_t Cst)
{
    slong val = WORD_MAX;

    for (slong i = 0; val > 0 && i < dop->length; i++)
    {
        gr_poly_struct * pol = ((gr_poly_struct *) (dop->coeffs)) + i;
        slong v = 0;
        while (v < val && v < pol->length
               && gr_is_zero(gr_poly_coeff_ptr(pol, v, Cst), Cst) == T_TRUE)
            v++;
        if (v < pol->length)
            val = v;
    }

    return val;
}

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

/* assumes Dop and Euler have the correct shape and the same base ring */
WARN_UNUSED_RESULT static int
_acb_ode_local_diffop(gr_ore_poly_t ldop, const gr_ore_poly_t dop,
                      gr_ctx_t Dop, gr_ore_poly_ctx_t Euler)
{
    slong p;
    gr_ctx_struct * Cst;
    int status = GR_SUCCESS;

    status |= gr_ore_poly_convert(ldop, &p, dop, Euler, Dop);

    status |= gr_ore_poly_ctx_over_gr_poly_base_ptrs(&Cst, NULL, Euler);
    slong xval = x_valuation(ldop, Cst);
    for (slong i = 0; i < ldop->length; i++)
    {
        gr_ptr pol = gr_ore_poly_coeff_ptr(ldop, i, Dop);
        status |= gr_poly_shift_right(pol, pol, xval, Cst);
    }

    return status;
}

WARN_UNUSED_RESULT int
acb_ode_local_diffop(gr_ore_poly_t ldop, ulong * flags, const gr_ore_poly_t dop,
                     gr_ctx_t Dop, gr_ore_poly_ctx_t Euler)
{
    gr_ctx_struct * Cst0, * Pol0, * Cst1, * Pol1;
    gr_ore_poly_struct * ldop0;
    int status = GR_SUCCESS;

    if (GR_ORE_POLY_CTX(Euler)->which_algebra != ORE_ALGEBRA_EULER_DERIVATIVE)
        return GR_UNABLE;

    if (dop->length == 0)
    {
        status |= gr_zero(ldop, Euler);
        FLINT_ASSERT(!(*flags & ACB_ODE_NONZERO));
        return status;
    }

    status |= gr_ore_poly_ctx_over_gr_poly_base_ptrs(&Cst0, &Pol0, Dop);
    status |= gr_ore_poly_ctx_over_gr_poly_base_ptrs(&Cst1, &Pol1, Euler);

    if (status != GR_SUCCESS)
        return status;

    slong order = dop->length - 1;
    gr_poly_struct * lc0 = (gr_poly_struct *) (dop->coeffs) + order;

    if (!(*flags & ACB_ODE_ORDINARY)
        && GR_ORE_POLY_CTX(Dop)->which_algebra == ORE_ALGEBRA_DERIVATIVE
        && gr_is_zero(lc0->coeffs, Cst0) == T_FALSE)
        *flags |= ACB_ODE_ORDINARY | ACB_ODE_REGULAR | ACB_ODE_NONZERO;

    /* Check REGULAR and NONZERO before changing the constant ring?
       But current callers are not interested in cases where one of these
       predicates would be true without being checkable on ldop. */

    {
        gr_ctx_t Euler0;

        gr_ore_poly_ctx_init(Euler0, Pol0, 0, ORE_ALGEBRA_EULER_DERIVATIVE);
        GR_TMP_INIT(ldop0, Euler0);

        /* todo: may eventually be subsumed by a more powerful gr_ore_poly
           implementation of gr_set_other */
        status |= _acb_ode_local_diffop(ldop0, dop, Dop, Euler0);
        if (Cst0->which_ring == GR_CTX_CC_ACB
            && Cst1->which_ring == GR_CTX_CC_ACB)
            /* Avoid rounding to the destination precision.
               todo: something similar for RR_ARB. */
            status |= gr_poly_set((gr_ptr) ldop, (gr_ptr) ldop0, Pol0);
        else
            status |= gr_poly_set_gr_poly_other((gr_ptr) ldop, (gr_ptr) ldop0,
                                                Pol0, Pol1);

        GR_TMP_CLEAR(ldop0, Euler0);
        gr_ctx_clear(Euler0);
    }

    if (status != GR_SUCCESS)
        return status;

    gr_poly_struct * lc = (gr_poly_struct *) (ldop->coeffs) + order;

    if (!(*flags & ACB_ODE_REGULAR)
        && gr_is_zero(lc->coeffs, Cst1) == T_FALSE)
        *flags |= ACB_ODE_REGULAR | ACB_ODE_NONZERO;

    if (!(*flags & ACB_ODE_NONZERO)
        && gr_is_zero(lc, Pol1) == T_FALSE)
        *flags |= ACB_ODE_NONZERO;

    return status;
}
