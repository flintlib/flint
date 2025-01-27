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
acb_theta_ql_nb_step_new(slong * pattern, const acb_mat_t tau, int cst, slong prec)
{
    slong g = acb_mat_nrows(tau);
    slong lp = ACB_THETA_LOW_PREC;
    arb_mat_t cho, yinv;
    arb_t x, t;
    slong s, nb, j;

    arb_init(x);
    arb_init(t);
    arb_mat_init(cho, g, g);
    arb_mat_init(yinv, g, g);

    acb_siegel_cho_yinv(cho, yinv, tau, lp);

    /* Compute rough pattern (could be negative) */
    /* Note: we could be more precise by scaling x by something else than
       powers of 2 */
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

        pattern[s] = -arf_get_si(arb_midref(x), ARF_RND_NEAR);
    }

    flint_printf("(ql_nb_steps) rough pattern:");
    for (s = 0; s < g; s++)
    {
        flint_printf(" %wd", pattern[s]);
    }
    flint_printf("\n");

    /* Experimental data from p-acb_theta_ql_exact: rough -> desired pattern */
    /* Genus 1 theta constants */
    /* 0, ..., 9 -> 0
       10 -> 3
       11, 12 -> 3 or 4
       13, 14 -> 5
       15 -> 6
       16 -> 7 */
    /* Genus 1 general theta values */
    /* 0, ..., 6 -> 0
       7, 8 -> 2 or 3 or 0
       9 -> 3 or 4
       10, 11 -> 5
       12 -> 6
       13 -> 6 or 7
       14 -> 7
       15 -> 8 */
    /* Genus 2 */
    /* 4 4, 5 5 -> 0 0
       6 6 -> 4 4
       7 7 -> 5 5 etc. up to 14 14
       4 2, 5 3 -> 0 0
       6 4 -> 2 2 or 3 3
       7 5 -> 3 3 or 4 4
       8 6 -> 5 5
       10 8 -> 7 7 etc.
       4 0, 6 2 -> 0 0
       7 3 -> 4 3
       8 4 -> 5 4 or 6 5
       9 5 -> 6 5 or 7 6
       11 7 -> 8 7 or 9 8
       13 9, 14 10 -> 11 10
       4 -2, 5 -1, 6 0 -> 0 0 or 1 0
       7 1 -> 1 0
       8 2 -> 4 3
       10 4 -> 6 5
       11 5 -> 7 6
       4 -4, ..., 7 -1 -> 1 0
       8 0 -> 2 1
       9 1 -> 4 2
       10 2 -> 5 2
       11 3 -> 6 3 or 8 5
       12 4 -> 8 4 or 6 5
       13 5 -> 7 5 or 8 6
       14 6 -> 9 6 or 9 7 */
    /* Genus 3 */
    /* 4 4 4 -> 4 4 4
       6 6 6 -> 6 6 6
       7 7 7 -> 7 7 7
       10 10 9 -> 10 10 10
       12 12 12 -> 12 12 12
       13 13 12 -> 13 13 13
       4 4 2 -> 2 2 2 or 3 3 3
       5 5 3 -> 4 4 4
       7 7 5 -> 6 6 6
       8 8 6 -> 8 8 7
       9 9 7 -> 9 9 8
       10 10 8 -> 10 0 9
       11 11 9 -> 11 11 10
       13 13 10 -> 12 12 11
       14 14 12 -> 13 13 13
       4 4 0 -> 2 2 2
       5 5 1 -> 4 4 4 or 4 4 3
       6 6 1, 6 6 2 -> 4 4 4 or 5 5 4
       7 7 3 -> 6 6 5
       8 8 4, 8 9 4 -> 8 8 7
       9 9 5 -> 8 8 7
       10 10 6 -> 9 9 8
       13 13 9 -> 12 12 12 or 12 12 11
       4 4 -2 -> 1 0 0 or 1 1 0
       5 5 -1 -> 3 3 0 or 2 2 0
       6 6 0 -> 4 4 2 or 3 3 1
       7 7 1 -> 6 6 3 or 5 5 3
       8 8 2 -> 6 6 5
       8 9 2 -> 7 7 5 or 7 7 4
       9 9 3, 10 10 3 -> 7 7 5
       10 10 4 -> 8 8 5 (?)
       11 11 5 -> 9 9 8
       12 12 6 -> 10 10 9
       13 13 7 -> 11 11 10
       4 4 -4 -> 1 1 0
       5 5 -3, 6 6 -3 -> 3 3 0
       7 7 -1 -> 5 5 0
       8 8 0 -> 6 6 3
       9 9 1 -> 7 7 4
       11 11 3 -> 9 9 6
       13 13 5 -> 11 11 8
       4 2 2 -> 2 2 2
       5 3 3 -> 3 3 3
       6 4 4 -> 4 4 4
       7 5 5 -> 5 5 5 or 6 6 6
       8 6 6 -> 8 7 7 or 7 6 6
       9 7 7 -> 7 7 7 or 8 8 8 or 9 8 8
       10 8 7 -> 8 8 8 or 9 8 8
       11 9 9 -> 9 9 9 or 10 10 10
       12 10 10 -> 10 10 10 or 11 11 11
       13 11 10, 13 11 11 -> 12 12 12
       4 0 0 -> 2 2 2 or 3 2 2
       5 1 1 -> 4 3 3
       6 2 1 -> 4 3 3
       7 3 3 -> 6 5 5
       8 4 4 -> 7 6 6
       9 5 5 -> 8 7 7
       13 9 8 -> 11 10 10
       4 -4 -4, ..., 8 0 0 -> 1 0 0
       9 1 1 -> 4 3 3
       10 2 2 -> 5 4 4 or 5 3 3
       11 3 3 -> 6 5 5 or 6 4 4 */

    /* Adjustments for s = 0 */
    if (g == 0 && cst)
    {
        pattern[0] -= 8;
        if (pattern[0] >= 6)
        {
            pattern[0] -= 1;
        }
        if (pattern[0] <= 2)
        {
            pattern[0] = 0;
        }
    }
    else if (g == 0 && !cst)
    {
        pattern[0] -= 6;
        if (pattern[0] <= 2)
        {
            pattern[0] = 0;
        }
    }
    else /* higher genus */
    {
        nb = 0;
        for (j = 1; j < g; j++)
        {
            nb = FLINT_MAX(nb, pattern[j]);
        }
        if (pattern[0] <= nb + 2) /* no dimension-lowering */
        {
            pattern[0] = FLINT_MIN(nb, pattern[0]);
        }
        else /* keep dimension-lowering */
        {
            pattern[0] = FLINT_MAX(pattern[0] - 5, FLINT_MAX(1, nb + 3));
        }
    }

    /* Adjustments for other s */
    for (s = 1; s < g; s++)
    {
        
    }
    
    for (s = 0; s < g; s++)
    {
        nb = pattern[s];
        if (s == 0)
        {
            nb -= 2;
            if (g == 1)
            {
                nb -= 1;
                if (nb > 8)
                {
                    /* a bit more balance between summation and duplication */
                    nb -= 1;
                }
            }
        }
        else if (s == 1)
        {
            /* summation in genus 2 is also quite efficient. */
            nb -= 2;
        }
        else if (s >= 3)
        {
            nb += 1;
        }
        pattern[s] = nb;
    }

    flint_printf("(ql_nb_steps) modified pattern:");
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
        if (s == g - 1)
        {
            nb = 0;
        }
        else
        {
            nb = pattern[s + 1];
        }

        /* Force trigger dimension-lowering at that point ? We only do this if
           s = 0 and nb is negative as the ellipsoid really contains very few
           points. */
        if (s == 0 && nb < 0)
        {
            pattern[s] = FLINT_MAX(1, pattern[s]);
        }

        /* Force more duplication steps ? We only do this for s = 1 if there
           will be a dimension-lowering later on, or for s >= 2. */
        if (s == 1 && pattern[s] < pattern[0])
        {
            j = FLINT_MIN(pattern[0] - 2, pattern[s] + 3);
            pattern[s] = j;
        }
        else if (s >= 2 && pattern[s] <= pattern[0] + 1 && pattern[s] >= 0)
        {
            pattern[s] += 2;
            if (pattern[s] >= pattern[0] - 1)
            {
                pattern[s] += 1;
            }
            pattern[s] = FLINT_MIN(pattern[s], pattern[0] + 1);
        }

        /* Remove further duplication steps in genus 1 if it doesn't mess with
           the dimension-lowering strategy */
        if (s == 0 && pattern[s] > nb + 1)
        {
            pattern[s] = FLINT_MAX(FLINT_MAX(0, nb) + 1, pattern[s] - 3);
        }
        /* In the case of genus 1 theta constants, be even more aggressive in
           avoiding duplication as acb_modular_theta_sum is even faster */
        if (s == 0 && nb == 0 && cst)
        {
            pattern[s] = FLINT_MAX(FLINT_MAX(0, nb) + 1, pattern[s] - 2);
        }

        /* Avoid making any duplication steps ? We only do this if the
           suggested number of steps for this s and the total number of steps
           are small. If s = 0 and nb < 0, we instead forced duplication to
           happen. */
        if (pattern[s] >= 1)
        {
            if ((s >= 2 && pattern[s] <= 1)
                || (s == 1 && pattern[s] <= 2)
                || (s == 0 && pattern[s] <= 3))
            {
                for (j = 0; j < g; j++)
                {
                    pattern[s] = 0;
                }
            }
            else if (s == 1 && pattern[s] <= 1)
            {
                pattern[s] = 0;
            }
        }

        /* Make pattern a nonincreasing vector */
        if (s < g - 1)
        {
            pattern[s] = FLINT_MAX(nb, pattern[s]);
        }
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
