/*
    Copyright (C) 2026 Marc Mezzarobba

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb.h"
#include "acb_ode.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "gr.h"
#include "gr_poly.h"
#include "gr_ore_poly.h"
#include "gr_vec.h"

WARN_UNUSED_RESULT int
acb_ode_lc_roots(acb_ptr lcrt,
                 const gr_ore_poly_t dop, gr_ctx_t Dop,
                 int include_zero,
                 gr_ctx_t CC)
{
    gr_ctx_struct * Cst;
    gr_vec_t sing;
    gr_poly_t pol;
    fmpz_vec_t sing_mult;
    slong val = 0;
    int status = GR_SUCCESS;

    status |= gr_ore_poly_ctx_over_gr_poly_base_ptrs(&Cst, NULL, Dop);
    if (status != GR_SUCCESS)
        return status;

    if (CC->which_ring != GR_CTX_CC_ACB)
        return GR_UNABLE;

    if (dop->length == 0)
        return GR_DOMAIN;

    gr_poly_struct * lc = (gr_poly_struct *) (dop->coeffs) + dop->length - 1;

    if (!include_zero)
    {
        for (; val < lc->length; val++)
        {
            int is_zero = gr_is_zero(gr_poly_coeff_ptr(lc, val, Cst), Cst);
            if (is_zero == T_FALSE)
                break;
            else if (is_zero == T_UNKNOWN)
                return GR_UNABLE;
        }
    }

    gr_poly_init(pol, Cst);
    gr_vec_init(sing, 0, CC);
    fmpz_vec_init(sing_mult, 0);

    status |= gr_poly_shift_right(pol, lc, val, Cst);

    /* todo: should use (multiple) root enclosures supporting polynomials
       with interval coefficients */
    int st1 = gr_poly_roots_other(sing, sing_mult, pol, Cst, 0, CC);
    for (slong i = 0; st1 == GR_UNABLE; i++)  /* should gr_poly_root_other for acb try this? */
    {
        gr_ctx_t ctx;
        gr_vec_t tmp;

        if (i == 0)
            gr_ctx_init_complex_algebraic_ca(ctx);
        else if (i == 1)
            gr_ctx_init_complex_qqbar(ctx);
        else
            break;

        gr_vec_init(tmp, 0, ctx);

        st1 = GR_SUCCESS;

        st1 |= gr_poly_roots_other(tmp, sing_mult, pol, Cst, 0, ctx);
        gr_vec_set_length(sing, tmp->length, CC);
        st1 |= _gr_vec_set_other(sing->entries, tmp->entries, ctx, tmp->length, CC);

        gr_vec_clear(tmp, ctx);
        gr_ctx_clear(ctx);
    }
    status |= st1;

    if (status == GR_SUCCESS)
    {
        for (slong j0 = 0, r = 0; j0 < sing->length; j0++)
        {
            slong m = fmpz_get_si((fmpz *) (sing_mult->entries) + j0);
            for (slong j1 = 0; j1 < m; j1++)
                acb_set(lcrt + r + j1, ((acb_ptr) sing->entries) + j0);
            r += m;
        }
    }

    gr_poly_clear(pol, Cst);
    gr_vec_clear(sing, CC);
    fmpz_vec_clear(sing_mult);

    return status;
}
