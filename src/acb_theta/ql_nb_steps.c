/*
    Copyright (C) 2023 Jean Kieffer

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb.h"
#include "arb_mat.h"
#include "acb_mat.h"
#include "acb_theta.h"

/* This is the all-important function to increase performance. */

int
acb_theta_ql_nb_steps(slong * pattern, const arb_mat_t cho, slong prec)
{
    slong g = arb_mat_nrows(cho);
    slong lp = ACB_THETA_LOW_PREC;
    arb_t x, t;
    slong s, nb;

    arb_init(x);
    arb_init(t);

    for (s = 0; s < g; s++)
    {
        arb_sqr(x, arb_mat_entry(cho, s, s), lp);
        arb_const_log2(t, lp);
        arb_div(x, x, t, lp);
        arb_div_si(x, x, prec, lp);
        arb_log(x, x, lp);
        arb_div(x, x, t, lp);

        if (!arb_is_finite(x) || arf_cmpabs_2exp_si(arb_midref(x), FLINT_BITS - 4) > 0)
        {
            arb_clear(x);
            arb_clear(t);
            return 0;
        }

        nb =  -arf_get_si(arb_midref(x), ARF_RND_NEAR);
        if (s == 0)
        {
            if (g == 1)
            {
                nb -= 7;
            }
            else if (g == 2)
            {
                nb -= 3;
            }
            else if (g <= 5)
            {
                nb -= 1;
            }
        }
        else
        {
            nb += 1;
        }
        pattern[s] = FLINT_MAX(0, nb);
    }

    arb_clear(x);
    arb_clear(t);
    return 1;
}
