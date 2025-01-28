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
            nb -= 3;
            if (g == 1)
            {
                if (nb > 8)
                {
                    nb -= 1;
                }
            }
        }
        else if (s == 1)
        {
            /* summation in genus 2 is also quite efficient. */
            nb -= 2;
        }
        else if (s == 2)
        {
            nb -= 1;
        }
        else
        {
            nb += 1;
        }

        pattern[s] = nb;
        nb_max = FLINT_MAX(nb_max, nb);
    }

    flint_printf("(ql_nb_steps) rough pattern:");
    for (s = 0; s < g; s++)
    {
        flint_printf(" %wd", pattern[s]);
    }
    flint_printf("\n");

    /* Start adapting the pattern from s = g-1 downwards. This is because the
      choice of whether to trigger dimension-lowering formulas in low
      dimensions will depend on whether or not duplications/dimension-lowerings
      have already been applied. */
    /* See /path/to/flint/build/acb_theta/profile/p-acb_theta_ql_exact */
    for (s = g - 1; s >= 0; s--)
    {
        /* Find out how many duplication steps have been performed so far
           (could be negative if s < g - 1) */
        nb = -10;
        for (j = s + 1; j < g; j++)
        {
            nb = FLINT_MAX(nb, pattern[j]);
        }
        if (s == g - 1)
        {
            nb = 0;
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
            j = FLINT_MIN(nb_max - 1, pattern[s] + 3);
            if (j <= 1)
            {
                j = -1; /* this is to still force nb < 0 */
            }
            pattern[s] = j;
        }

        /* Remove further duplication steps in genus 1 if it doesn't mess with
           the dimension-lowering strategy (keep it triggered if nb < 0) */
        if (s == 0 && pattern[s] > nb + 1)
        {
            pattern[s] = FLINT_MAX(FLINT_MAX(0, nb) + 1, pattern[s] - 3);
            nb_max = pattern[s];
        }
        /* In the case of genus 1 theta constants, we are even more aggressive
           in avoiding duplication as acb_modular_theta_sum is even faster */
        if (s == 0 && nb == 0 && cst)
        {
            pattern[s] = FLINT_MAX(FLINT_MAX(0, nb) + 1, pattern[s] - 2);
            nb_max = pattern[s];
        }

        /* Avoid making any duplication steps ? We only do this if the
           suggested number of steps for this s and the total number of steps
           are small. (If nb < 0, we instead forced duplication to happen.) */
        if (nb == 0 && pattern[s] >= 1 && pattern[s] <= 2)
        {
            if ((s >= 2 && nb_max <= 1)
                || (s <= 1 && nb_max <= 2))
            {
                for (j = 0; j <= s; j++)
                {
                    pattern[s] = 0;
                }
            }
            if (s == 1 && pattern[s] <= 1)
            {
                pattern[s] = 0;
            }
        }

        /* Make pattern a nonincreasing vector */
        if (s < g - 1)
        {
            pattern[s] = FLINT_MAX(nb, pattern[s]);
        }
        flint_printf("s = %wd, set pattern to %wd\n", s, pattern[s]);
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
