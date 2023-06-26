/*
    Copyright (C) 2023 Jean Kieffer

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "acb_theta.h"

void
acb_theta_eld_round(slong * r, const arb_mat_t v)
{
    slong g = arb_mat_nrows(v);
    slong j;

    for (j = 0; j < g; j++)
    {
        if (!arb_is_finite(arb_mat_entry(v, j, 0))
            || arf_cmpabs_ui(arb_midref(arb_mat_entry(v, j, 0)), WORD_MAX) > 0)
        {
            flint_printf("acb_theta_eld_round: Error (impossible rounding)\n");
            fflush(stdout);
            flint_abort();
        }
        r[j] = arf_get_si(arb_midref(arb_mat_entry(v, j, 0)), ARF_RND_NEAR);
    }
}
