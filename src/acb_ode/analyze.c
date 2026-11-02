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

ulong
acb_ode_analyze(const gr_ore_poly_t dop, ulong flags, gr_ctx_t Dop)
{
    gr_ctx_struct * Cst, * Pol;

    flags &= ACB_ODE_ORDINARY | ACB_ODE_REGULAR | ACB_ODE_NONZERO;

    if (gr_ore_poly_ctx_over_gr_poly_base_ptrs(&Cst, &Pol, Dop) != GR_SUCCESS)
        return flags;

    if (dop->length == 0)
    {
        FLINT_ASSERT(!(flags & ACB_ODE_NONZERO));
        return flags;
    }

    ore_algebra_t alg = GR_ORE_POLY_CTX(Dop)->which_algebra;
    slong order = dop->length - 1;
    gr_poly_struct * lc = (gr_poly_struct *) (dop->coeffs) + order;

    if (!(flags & ACB_ODE_ORDINARY)
        && alg == ORE_ALGEBRA_DERIVATIVE
        && gr_is_zero(lc->coeffs, Cst) == T_FALSE)
        flags |= ACB_ODE_ORDINARY | ACB_ODE_REGULAR | ACB_ODE_NONZERO;

    if (!(flags & ACB_ODE_REGULAR)
        && alg == ORE_ALGEBRA_EULER_DERIVATIVE
        && gr_is_zero(lc->coeffs, Cst) == T_FALSE)
        flags |= ACB_ODE_REGULAR | ACB_ODE_NONZERO;

    if (!(flags & ACB_ODE_NONZERO)
        && gr_is_zero(lc, Pol) == T_FALSE)
        flags |= ACB_ODE_NONZERO;

    return flags;
}
