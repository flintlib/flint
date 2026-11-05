/*
    Copyright (C) 2026 Marc Mezzarobba

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb_poly.h"
#include "acb_ode.h"
#include "gr.h"
#include "gr_ore_poly.h"
#include "gr_poly.h"
#include "ulong_extras.h"

int
acb_ode_solution_series(gr_poly_vec_t series,
                        const gr_ore_poly_t dop, gr_ore_poly_ctx_t Dop,
                        const acb_ode_group_t group,
                        const gr_vec_t ini,
                        slong n,
                        ulong flags,
                        gr_ctx_t CC)
{
    gr_ctx_t wCC, wCCx, Euler;
    gr_ore_poly_struct * ldop;
    slong prec;

    int status = GR_SUCCESS;

    status |= gr_ore_poly_ctx_over_gr_poly_base_ptrs(NULL, NULL, Dop);
    if (status != GR_SUCCESS)
        return status;

    if (CC->which_ring != GR_CTX_CC_ACB)
        return GR_UNABLE;

    slong dop_clen = 0;
    for (slong i = 0; i < dop->length; i++)
    {
        const gr_poly_struct * c = gr_ore_poly_coeff_srcptr(dop, i, Dop);
        dop_clen = FLINT_MAX(dop_clen, c->length);
    }

    /* todo: better estimate, linear error analysis */
    status |= gr_ctx_get_real_prec(&prec, CC);
    slong wp = 8 + prec + n * n_clog(FLINT_MAX(1, dop_clen - 3), 2);

    gr_ctx_init_complex_acb(wCC, wp);
    acb_ode_local_diffop_ctx_init(wCCx, Euler, Dop, wCC);
    GR_TMP_INIT(ldop, Euler);

    status |= acb_ode_local_diffop(ldop, &flags, dop, Dop, Euler);

    if (!(flags & ACB_ODE_REGULAR))
        status |= GR_UNABLE;

    if (status != GR_SUCCESS)
        ;
    else if (ldop->length == 1)
    {
        gr_poly_vec_set_length(series, 0, CC);
    }
    else
    {
        acb_ode_sum_t sum;
        mag_t dummy;

        mag_init(dummy);
        acb_ode_sum_init(sum, ldop->length, 0, 1, 0);

        acb_ode_sum_set_diffop(sum, ldop->coeffs, ldop->length, dummy);
        acb_ode_sum_set_group(sum, group);
        acb_ode_sol_set_ini(sum->sol, ini->entries, group);
        if (flags & ACB_ODE_ORDINARY)
            sum->sol->future_logs = 1;
        sum->flags = ACB_ODE_WANT_SERIES;

        acb_ode_sum_divconquer(sum, n, NULL, NULL, wp);

        FLINT_ASSERT(sum->n0 == 0);
        FLINT_ASSERT(sum->n == n);

        gr_poly_vec_set_length(series, sum->sol->nlogs, CC);
        for (slong k = 0; k < sum->sol->nlogs; k++)
        {
            acb_poly_truncate(sum->sol->series + k, n);
            _acb_poly_normalise(sum->sol->series + k);
            acb_poly_swap((void *) gr_poly_vec_entry_ptr(series, k, CC),
                          sum->sol->series + k);
        }

        acb_ode_sum_clear(sum);
        mag_clear(dummy);
    }

    GR_TMP_CLEAR(ldop, Euler);
    gr_ctx_clear(wCCx);
    gr_ctx_clear(wCC);
    gr_ctx_clear(Euler);

    return status;
}
