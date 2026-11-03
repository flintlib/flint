/*
    Copyright (C) 2026 Marc Mezzarobba

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <math.h>

#include "acb_ode.h"
#include "arf.h"
#include "longlong.h"
#include "mag.h"

int
_acb_ode_sum_done(acb_ode_sum_struct * sum, slong stride,
                  /* in the future sum_done may refine bound, gbound */
                  acb_ode_bound_t bound, acb_ode_group_bound_t gbound,
                  slong prec)
{
    mag_t eps, est, sum_mag, sum_rad;

    if (sum->nder == 0)
        return 0;

    mag_init(eps);
    mag_init(est);
    mag_init(sum_mag);
    mag_init(sum_rad);

    slong deg = sum->dop_clen - 1;
    FLINT_ASSERT(sum->n >= sum->n0 + deg);  /* todo: lift this limitation */

    /* detect solutions that have converged or will not improve */

    double progress = NAN;
    int have_tail_bounds = 0;

    int stop = 1;
    for (slong m = 0; m < sum->nsols; stop = sum->sol[m++].done && stop)
    {
        acb_ode_sol_struct * sol = sum->sol + m;

        if (sol->done)
            continue;

        /* todo: non-rigorous check with less than deg terms?
           (maybe using known terms from the previous residuals in addition to
           the available coefficients) */
        slong off = sum->n - deg - sum->n0;
        slong accuracy = acb_ode_sol_estimate_terms(est, sol, off, deg,
                                                    sum->magpow);
        slong sum_accuracy = acb_ode_sol_estimate_sums(sum_mag, sum_rad, sol);
        acb_ode_cvest_update(sol->cvest, sol->cvest,
                             est, accuracy, stride, prec, sum->wp);

        /* Attempt to detect cases where it is better to bail out even though,
           technically, we may still be able to improve the results (typically
           entire functions with a large hump to climb).
           XXX Not sure if this belongs here; maybe find a better criterion. */
        if (mag_cmp_2exp_si(sum_rad, 4 * prec) >= 0
                && (sum_accuracy < 0 || sum->n > prec * FLINT_BIT_COUNT(prec)))
            mag_inf(eps);
        /* sum_accuracy may be <= 0 and still have a chance to improve,
           typically when the terms are still increasing.

           If accuracy <= 0, the approximation of the sum is unlikely to
           improve, but the error bound still might (e.g., if the evaluation
           point contains zero), roughly until the radius of the partial sum
           exceeds the tail bound. */
        else if (accuracy >= 0 || sum_accuracy >= 0)
        {
            /* stopping condition: ~prec correct bits OR working precision
               insufficient to do better */
            mag_mul_2exp_si(eps, sum_mag, -prec);
            mag_max(eps, eps, sum_rad);
            mag_mul_2exp_si(eps, eps, -8);
        }
        else
            /* todo: The point about the bound having a chance to improve still
               holds. A better heuristic (maybe based on the rate at which the
               tail bounds are improving?) should improve evaluations on wide
               intervals (typically containing 0). */
            mag_inf(eps);

        slong max_shift = sum->group->shifts[sum->group->nshifts - 1].n;

        /* quick non-rigorous check */
        if (accuracy >= 0 && mag_cmp(est, eps) > 0)
        {
            sol->done = 0;
            continue;
        }
        else if (sum->flags & ACB_ODE_APPROX)
        {
            /* todo: test initial conditions at the remaining shifts */
            if (sum->n < max_shift + sum->dop_clen)
                sol->done = 0;
            else
                sol->done = 16;
            continue;
        }

        /* rigorous check */
        if (!have_tail_bounds)
        {
            /* XXX wasteful when only some sols have converged */
            progress = _acb_ode_sum_update_tail_bounds(sum, bound, gbound,
                                                       MAG_BITS);
            have_tail_bounds = 1;
        }

        int all_tb_ok = 1;
        for (slong i = 0; i < sum->nder; i++)
            if (mag_cmp(sol->tb + i, eps) > 0)
                all_tb_ok = 0;

        if (all_tb_ok)
            sol->done = 1;
        else if (sum_accuracy == -ARF_PREC_EXACT)
            sol->done = 2;
        /* the following cases are heuristics that should only be applied past
           the last singular index */
        else if (sum->n < max_shift + sum->dop_clen)
            sol->done = 0;
        else if (progress < 0.)  /* bounds are worsening */
            sol->done = 3;
        /* Crude heuristic to detect situations where the tail bound might as
           well be infinite. Note that the previous tail bound may actually be
           more than stride terms old.
           XXX Maybe should take into account indicial roots from other groups.
         */
        else if (sum->nder > 0
                 && mag_get_d_log2_approx(sol->tb + 0)
                        > mag_get_d_log2_approx(eps) + prec
                          + sum->n * progress / stride)
            sol->done = 4;
        else
            sol->done = 0;
    }

    mag_clear(eps);
    mag_clear(est);
    mag_clear(sum_mag);
    mag_clear(sum_rad);

    return stop;
}
