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
acb_theta_ql_nb_steps(slong * pattern, const acb_mat_t tau, int cst, slong prec)
{
    slong g = acb_mat_nrows(tau);
    slong lp = ACB_THETA_LOW_PREC;
    arb_mat_t cho, yinv;
    arb_t x, t;
    slong s, nb;

    arb_init(x);
    arb_init(t);
    arb_mat_init(cho, g, g);
    arb_mat_init(yinv, g, g);

    acb_siegel_cho_yinv(cho, yinv, tau, lp);

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

        /* Adapt pattern in light of experimental performance */
        /* See /path/to/flint/build/acb_theta/profile/p-acb_theta_ql_exact */
        if (s == 0)
        {
            /* Rationale is: acb_modular_theta_sum is so fast that we don't
               need so many duplication steps. */
            if (g == 1 && cst)
            {
                nb -= 8;
            }
            else
            {
                nb -= 5;
            }
            /* One less step if at least 9 */
            /* Rationale is: more balance with summation phase */
            if (nb > 8)
            {
                nb -= 1;
            }
            /* Never make <= 2 duplication steps. */
            /* Rationale is: this avoids having the slight overhead of
               computing distances, etc. */
            if (nb <= 2)
            {
                nb = 0;
            }
        }
        if (s == 1)
        {
            /* Rationale is: summation in genus 2 is also quite efficient. */
            nb -= 1;
            if (nb < pattern[0])
            {
                /* Rationale is: if we use the dimension-lowering strategy, then
                   one more duplication step will nicely decrease the number of
                   auxiliary points. */
                nb = FLINT_MIN(pattern[0] - 1, nb + 2);
            }
            /* Never make <= 2 duplication steps. */
            /* Rationale is: this avoids having the slight overhead of
               computing distances, etc. */
            if (nb <= 2)
            {
                nb = 0;
            }
        }
        pattern[s] = FLINT_MAX(0, nb);
    }

    /* Clean up: make pattern a nonincreasing vector */
    for (s = g - 1; s >= 1; s--)
    {
        pattern[s - 1] = FLINT_MAX(pattern[s - 1], pattern[s]);
    }
    /* This seems to work experimentally. */
    if (g >= 2 && pattern[1] == 0)
    {
        pattern[0] = FLINT_MAX(0, pattern[0] - 1);
    }

    arb_clear(x);
    arb_clear(t);
    arb_mat_clear(cho);
    arb_mat_clear(yinv);
    return 1;
}
