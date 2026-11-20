/*
    Copyright (C) 2026 Marc Mezzarobba

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/* generated using claude-fable-5, claude-opus-5-5 */

#include "test_helpers.h"
#include "acb.h"
#include "acb_poly.h"
#include "acb_ode.h"

TEST_FUNCTION_START(acb_ode_sum_forward, state)
{
    for (slong iter = 0; iter < 100 * flint_test_multiplier(); iter++)
    {
        acb_ode_exponents_t expos;
        acb_ode_sum_struct sum[2];
        mag_t cvrad, mag;
        ulong seed1, seed2;

        flint_rand_get_seed(&seed1, &seed2, state);

        slong dop_len = 2 + n_randint(state, 4);
        slong lcdeg = n_randint(state, 3);
        slong prec = 2 + n_randint(state, 200);

        acb_poly_struct * dop = _acb_poly_vec_init(dop_len);
        acb_ptr lcroots = _acb_vec_init(lcdeg);
        acb_ode_exponents_init(expos);
        mag_init(cvrad);
        mag_init(mag);

        acb_ode_randtest_acb(dop, lcroots, expos, state, n_randint(state, 11),
                             lcdeg, 1 + n_randint(state, 6), dop_len, prec);
        _acb_vec_get_mag_lower(cvrad, lcroots, lcdeg);
        slong dop_clen = _acb_poly_vec_length(dop, dop_len);
        slong d = dop_clen - 1;

        acb_ode_group_struct * group = expos->groups
                                     + n_randint(state, expos->ngroups);
        int highest = n_randint(state, 2);
        slong maxsols = highest ? group->nshifts : acb_ode_group_length(group);
        slong nsols = n_randint(state, maxsols + 1);
        slong npts = mag_is_zero(cvrad) ? 0 : n_randint(state, 4);
        slong nder = n_randint(state, 4);
        slong nterms = n_randint(state, 8 * dop_clen + 1);
        slong block_len = FLINT_MAX(1, d) + n_randint(state, 4);
        ulong flags = (n_randint(state, 2) ? ACB_ODE_WANT_SERIES : 0)
                    | (n_randint(state, 2) ? ACB_ODE_APPROX : 0);

        acb_ptr pts = _acb_vec_init(npts);
        for (slong p = 0; p < npts; p++)
        {
            acb_randtest(pts + p, state, prec, 2);
            acb_get_mag(mag, pts + p);
            mag_div(mag, mag, cvrad);
            if (mag_cmp_2exp_si(mag, 0) > 0)
                acb_mul_2exp_si(pts + p, pts + p,
                                -fmpz_get_si(MAG_EXPREF(mag)));
        }

        for (int i = 0; i < 2; i++)
        {
            acb_ode_sum_init(sum + i, dop_len, npts, nsols, nder);
            acb_ode_sum_set_diffop(sum + i, dop, dop_len, cvrad);
            acb_ode_sum_set_group(sum + i, group);
            if (highest)
                acb_ode_sum_set_ini_highest(sum + i);
            else
                acb_ode_sum_set_ini_echelon(sum + i);
            acb_ode_sum_set_points(sum + i, pts, npts);
            sum[i].flags = flags;

            acb_ode_sum_precompute(sum + i);
            sum[i].wp = sum[i].sum_wp = prec;
            acb_ode_sum_fit_length(sum + i, nterms + dop_clen);
        }

        for (slong n = 0; n < nterms; n++)
        {
            _acb_ode_sum_forward_1(sum + 0);
            _acb_ode_sum_update_residuals(sum + 0, n, n + 1, n + 1 + d);
        }

        /* random mix of forward_1 steps and divconquer blocks, followed in
           both cases by the completion of the residual coefficients n..n+d-1 */
        while (sum[1].n < nterms)
        {
            slong n = sum[1].n;

            if (nterms - n > d && n_randint(state, 2))
            {
                slong high = n_randint(state, 2) ? nterms
                           : n + d + 1 + n_randint(state, nterms - n - d);
                _acb_ode_sum_forward_divconquer(sum + 1, high, block_len,
                                                block_len, NULL, NULL, prec);
                _acb_ode_sum_update_residuals(sum + 1, high - d, high,
                                              high + d);
            }
            else
            {
                _acb_ode_sum_forward_1(sum + 1);
                _acb_ode_sum_update_residuals(sum + 1, n, n + 1, n + 1 + d);
            }
        }

        _acb_ode_sum_fix(sum + 0);
        _acb_ode_sum_fix(sum + 1);

        if (sum[0].n != nterms || sum[1].n != nterms)
        {
            flint_printf("FAIL (nterms)\n\n");
            goto fail;
        }

        for (slong m = 0; m < nsols && !(flags & ACB_ODE_APPROX); m++)
        {
            acb_ode_sol_struct * sol0 = sum[0].sol + m;
            acb_ode_sol_struct * sol1 = sum[1].sol + m;
            slong n0 = sum[1].n0;

            for (slong k = 0; k < dop_len - 1; k++)
            {
                if (!_acb_vec_overlaps(sol0->series[k].coeffs + n0,
                                       sol1->series[k].coeffs, nterms - n0))
                {
                    flint_printf("FAIL (series)\n\n");
                    flint_printf("m = %wd, k = %wd, n0 = %wd\n", m, k, n0);
                    goto fail;
                }
            }

            if (!_acb_vec_overlaps(sol0->sums, sol1->sums,
                                   npts * (dop_len - 1) * nder))
            {
                flint_printf("FAIL (sums)\n\n");
                flint_printf("m = %wd\n", m);
                goto fail;
            }
        }

        acb_ode_sum_clear(sum + 0);
        acb_ode_sum_clear(sum + 1);
        _acb_vec_clear(pts, npts);
        _acb_vec_clear(lcroots, lcdeg);
        _acb_poly_vec_clear(dop, dop_len);
        acb_ode_exponents_clear(expos);
        mag_clear(cvrad);
        mag_clear(mag);

        continue;

fail:
        flint_printf("flint_rand_set_seed(state, 0x%wx, 0x%wx);\n",
                     seed1, seed2);
        for (slong i = 0; i < dop_len; i++)
            flint_printf("dop[%wd] = %{acb_poly}\n", i, dop + i);
        flint_printf("expos = ");
        acb_ode_exponents_println(expos);
        flint_printf("pts = %{acb*}\n", pts, npts);
        flint_printf("highest = %d, nsols = %wd, nder = %wd\n",
                     highest, nsols, nder);
        flint_printf("nterms = %wd, n = %wd/%wd, block_len = %wd, "
                     "flags = %wu, prec = %wd\n", nterms, sum[0].n,
                     sum[1].n, block_len, flags, prec);
        flint_abort();
    }

    TEST_FUNCTION_END(state);
}
