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
    slong s, nb, dupl;
    slong * rough;

    arb_init(x);
    arb_init(t);
    arb_mat_init(cho, g, g);
    arb_mat_init(yinv, g, g);
    rough = flint_malloc(g * sizeof(slong));

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
            /* Should not happen in tests */
            arb_clear(x);
            arb_clear(t);
            return 0;
        }

        rough[s] = -arf_get_si(arb_midref(x), ARF_RND_NEAR);
    }

    /* Experimental data from p-acb_theta_ql_exact: rough -> desired pattern */
    /* Genus 1 theta constants */
    /* 0, ..., 9 -> 0
       10 -> 3
       11 -> 4
       12, 13 -> 5
       14, 15 -> 6
       16 -> 7
       17 -> 8 */
    /* Genus 1 general theta values */
    /* 0, ..., 8 -> 0
       9 -> 0 or 4
       10 -> 5
       11 -> 6, ...,
       13, 14 -> 8
       15, 16 -> 9
       17 -> 10 */
    /* Genus 2 */
    /* 4 3, 5 3 -> 0 0
       6 4 -> 3 3 or 4 4
       6 5 -> 4 4 or 5 5
       8 6 -> 6 6
       11 10 -> 10 10
       12 10 -> 10 10 or 11 11
       16 14 -> 14 14 or 15 15
       8 2 -> 3 2 or 4 2
       9 3 -> 4 2
       10 4 -> 6 5 or 7 6
       11 5 -> 7 6
       12 6 -> 8 7
       15 9 -> 11 10
       8 0 -> 2 0 or 0 0
       9 1 -> 4 2
       10 2 -> 6 2
       12 4 -> 8 6
       14 6 -> 8 7 */

    /* Start adapting the pattern from s = g-1 downwards. This is because the
      choice of whether to trigger dimension-lowering formulas in low
      dimensions will depend on whether or not duplications/dimension-lowerings
      have already been applied. */
    /* See /path/to/flint/build/acb_theta/profile/p-acb_theta_ql_exact */
    /* Some of these branches will not show up in tests. */
    for (s = g - 1; s >= 0; s--)
    {
        /* Find out how many duplication steps have been performed so far
           (could be negative if s < g - 1) */
        if (s == g - 1)
        {
            dupl = 0;
        }
        else
        {
            dupl = pattern[s + 1];
        }
        pattern[s] = rough[s];

        /* Force trigger dimension-lowering at that point ? We only do this if
           s = 0 and dupl is negative as the ellipsoid really contains very few
           points. */
        if (s == 0 && dupl < 0)
        {
            pattern[s] = FLINT_MAX(1, pattern[s]);
        }

        /* Force more duplication steps ? We only do this for s = 1 if there
           will be a dimension-lowering later on, or for s >= 2. */
        if (s == 1 && pattern[s] < rough[0] && pattern[s] >= 1)
        {
            pattern[s] = FLINT_MIN(rough[0] - 1, pattern[s] + 2);
        }
        else if (s >= 2 && pattern[s] >= 1)
        {
            pattern[s] += 2;
        }

        /* Remove further duplication steps in genus 1 if it doesn't mess with
           the dimension-lowering strategy. */
        if (s == 0)
        {
            nb = pattern[s] - 5;
            if (nb >= 10)
            {
                nb -= 2;
            }
            else if (nb >= 8)
            {
                nb -= 1;
            }
            if (g == 1 && cst)
            {
                nb -= 2;
            }
            pattern[s] = FLINT_MAX(FLINT_MAX(0, dupl) + 1, nb);
        }

        /* Remove duplication steps to avoid dimension-lowering altogether ? We
           only do this if the number of steps is <= dupl + 2. In genus 1, we
           additionally demand that rough[s] be less than dupl + 1, unless dupl
           == 0 */
        if (pattern[s] <= dupl + 2
            && (s > 0 || dupl == 0))
        {
            pattern[s] = dupl;
        }
        else if (s == 0 && rough[s] <= dupl + 1)
        {
            pattern[s] = dupl;
        }

        /* In any case, make pattern a nonincreasing vector */
        if (s < g - 1)
        {
            pattern[s] = FLINT_MAX(dupl, pattern[s]);
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
    flint_free(rough);
    return 1;
}
