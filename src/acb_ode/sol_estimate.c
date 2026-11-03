/*
    Copyright (C) 2026 Marc Mezzarobba

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb.h"
#include "acb_ode.h"
#include "arf.h"
#include "mag.h"

slong
acb_ode_sol_estimate_terms(mag_t est,
                           const acb_ode_sol_t sol, slong off, slong len,
                           const mag_t radpow)
{
    mag_t m;
    int all_bad = 1;

    mag_init(m);

    slong accuracy = ARF_PREC_EXACT;
    mag_zero(est);

    for (slong k = 0; k < sol->nlogs; k++)
    {
        for (slong i = 0; i < len; i++)
        {
            acb_ptr c = (sol->series + k)->coeffs + off + i;

            acb_get_mag(m, c);
            mag_add(est, est, m);

            slong acc = acb_rel_accuracy_bits(c);
            /* This is a bit of a kludge to prevent zeros in the coefficient
               sequence from making the overall accuracy estimate negative.
               Maybe returning the absolute accuracy of c * radpow instead
               would actually be more useful (what exactly does cvest_update
               need?). */
            if (acc >= 0)
                all_bad = 0;
            else
                acc = acb_rel_one_accuracy_bits(c);
            accuracy = FLINT_MIN(accuracy, acc);
        }
    }

    if (all_bad)
        accuracy = - ARF_PREC_EXACT;

    mag_mul(est, est, radpow);

    mag_clear(m);
    return accuracy;
}

slong
acb_ode_sol_estimate_sums(mag_t mag, mag_t rad, const acb_ode_sol_t sol)
{
    mag_t tmp;

    mag_init(tmp);

    mag_zero(mag);
    mag_zero(rad);
    slong accuracy = ARF_PREC_EXACT;

    for (slong p = 0; p < sol->npts; p++)
        for (slong k = 0; k < sol->nlogs; k++)
            for (slong i = 0; i < sol->nder; i++)
            {
                acb_ptr c = acb_ode_sol_sum_ptr(sol, p, k, i);

                acb_get_mag(tmp, c);
                mag_max(mag, mag, tmp);

                mag_max(rad, rad, &c->real.rad);
                mag_max(rad, rad, &c->imag.rad);

                accuracy = FLINT_MIN(accuracy, acb_rel_accuracy_bits(c));
            }

    mag_clear(tmp);

    return accuracy;
}
