/*
    Copyright (C) 2025, 2026 Marc Mezzarobba

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb_ode.h"

void
acb_ode_sum_divconquer(acb_ode_sum_t sum, slong nterms, acb_ode_bound_t bound,
                       acb_ode_group_bound_t gbound, slong prec)
{
    FLINT_ASSERT((bound == NULL) == (gbound == NULL));

    if (nterms < 0)
    {
        if (bound == NULL || sum->nder == 0)
            flint_throw(FLINT_ERROR, "Unlimited terms, no stopping criterion");
        nterms = WORD_MAX;
    }

    acb_ode_sum_precompute(sum);
    acb_ode_sum_zero(sum);
    /* todo: take into account nterms and the APPROX flag */
    /* todo: make it possible to customize the working precision at this level? */
    sum->sum_wp = acb_ode_choose_prec(&sum->wp, sum->dop, sum->dop_len,
                                      sum->mag, sum->cvrad, prec);
    for (slong m = 0; m < sum->nsols; m++)
    {
        sum->sol[m].cvest->accuracy = sum->wp;
        sum->sol[m].cvest->loss_rate = sum->dop_len * sum->dop_clen;
    }

    slong block_len = FLINT_MAX(1, sum->dop_clen - 1);
    slong stride = block_len;  /* todo: less when block_len is large */

    for (slong m = 0; m < sum->nsols; m++)
        acb_ode_sol_fit_length(sum->sol + m, 2 * block_len);

    _acb_ode_sum_forward_divconquer(sum, nterms, block_len, stride, bound,
                                    gbound, prec);
    _acb_ode_sum_fix(sum);
}
