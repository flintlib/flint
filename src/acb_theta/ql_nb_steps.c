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
    slong s, nb, j, nb_max;

    arb_init(x);
    arb_init(t);
    arb_mat_init(cho, g, g);
    arb_mat_init(yinv, g, g);

    acb_siegel_cho_yinv(cho, yinv, tau, lp);

    /* Compute rough pattern (could be negative) */
    /* Note: we could be more precise by scaling x by something else than
       powers of 2. */
    nb_max = 0;
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

        nb = -arf_get_si(arb_midref(x), ARF_RND_NEAR);
        if (s == 0)
        {
            /* acb_modular_theta_sum is so fast that we don't need so many
               duplication steps. */
            nb -= 5;
        }
        else if (s == 1)
        {
            /* summation in genus 2 is also quite efficient. */
            nb -= 3;
        }
        else if (s == 2)
        {
            nb -= 1;
        }
        else
        {
            nb += 1;
        }
        /* One less step if at least 9 for more balance with summation phase */
        if (nb > 8)
        {
            nb -= 1;
        }

        pattern[s] = nb;
        nb_max = FLINT_MAX(nb_max, nb);
    }

    /* Start adapting the pattern from s = g-1 downwards. This is because the
      choice of whether to trigger dimension-lowering formulas in low
      dimensions will depend on whether or not duplications/dimension-lowerings
      have already been applied. */
    /* See /path/to/flint/build/acb_theta/profile/p-acb_theta_ql_exact */
    for (s = g - 1; s >= 0; s--)
    {
        /* Find out how many duplication steps have been performed so far
           (could be negative) */
        nb = -10;
        for (j = s + 1; j < g; j++)
        {
            nb = FLINT_MAX(nb, pattern[j]);
        }

        /* Force trigger dimension-lowering at that point ? We only do this if
           nb is negative as the ellipsoid really contains very few points. */
        if (nb < 0 && (s == 0 || pattern[s] == 0))
        {
            pattern[s] = FLINT_MAX(1, pattern[s]);
        }

        /* Force more duplication steps ? We only do this if there will be a
           dimension-lowering later on. */
        if (pattern[s] < nb_max)
        {
            pattern[s] = FLINT_MIN(nb_max - 1, pattern[s] + 3);
        }

        /* Avoid making any duplication steps ? We only do this if the
           suggested number of steps for this s and the total number of steps
           are small. */
        if (nb == 0 && pattern[s] <= 2)
        {
            if (s >= 2 && nb_max <= 1)
            {
                for (j = 0; j <= s; j++)
                {
                    pattern[s] = 0;
                }
            }
            else if (s == 1 && nb_max + pattern[s] <= 5)
            {
                pattern[0] = FLINT_MAX(pattern[0] - 2, 0);
                pattern[s] == 0;
            }
            else if (s == 0)
            {
                pattern[s] = 0;
            }
        }
        /* In the case of genus 1 theta constants, we are more aggressive in
           avoiding duplication, because acb_modular_theta_sum is faster for constants */
        if (s == 0 && nb == 0 && cst)
        {
            pattern[0] = pattern[0] - 3;
            if (pattern[0] <= 2)
            {
                pattern[0] = 0;
            }
        }
        /* One less duplication step in for s = 0 if it doesn't mess with
           dimension lowering */
        if (s == 0 && pattern[0] > nb + 1)
        {
            pattern[0] -= 1;
        }

        /* Make pattern a nonincreasing vector */
        pattern[s] = FLINT_MAX(0, nb);
    }

    /* Clean up: make pattern a nonnegative vector */
    for (s = 0; s < g; s++)
    {
        pattern[s] = FLINT_MAX(0, pattern[s]);
    }

    arb_clear(x);
    arb_clear(t);
    arb_mat_clear(cho);
    arb_mat_clear(yinv);
    return 1;
}
