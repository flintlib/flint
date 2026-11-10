/*
    Copyright (C) 2025, 2026 Marc Mezzarobba

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

TEST_FUNCTION_START(acb_ode_sum_divconquer, state)
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
        slong prec[2] = { 2 + n_randint(state, 200), 2 + n_randint(state, 200) };

        acb_poly_struct * dop = _acb_poly_vec_init(dop_len);
        acb_ptr lcroots = _acb_vec_init(lcdeg);
        acb_ode_exponents_init(expos);
        mag_init(cvrad);
        mag_init(mag);

        acb_ode_randtest_acb(dop, lcroots, expos, state, n_randint(state, 11),
                             lcdeg, 1 + n_randint(state, 6), dop_len, prec[0]);
        _acb_vec_get_mag_lower(cvrad, lcroots, lcdeg);
        slong dop_clen = _acb_poly_vec_length(dop, dop_len);

        acb_ode_group_struct * group = expos->groups
                                     + n_randint(state, expos->ngroups);
        int highest = n_randint(state, 2);
        slong maxsols = highest ? group->nshifts : acb_ode_group_length(group);
        slong nsols = n_randint(state, maxsols + 1);
        slong npts = mag_is_zero(cvrad) ? 0 : n_randint(state, 4);
        slong nder = n_randint(state, 4);
        slong nterms = n_randint(state, 8 * dop_clen + 1);

        acb_ptr pts = _acb_vec_init(npts);
        for (slong p = 0; p < npts; p++)
        {
            acb_randtest(pts + p, state, prec[0], 2);
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
        }

        sum[0].flags = ACB_ODE_WANT_SERIES;
        sum[1].flags = (n_randint(state, 2) ? ACB_ODE_WANT_SERIES : 0)
                     | (n_randint(state, 2) ? ACB_ODE_APPROX : 0);

        /* leftovers from a previous run must not matter */
        acb_ode_sum_divconquer(sum + 1, n_randint(state, 30), NULL, NULL,
                               prec[1]);

        for (int i = 0; i < 2; i++)
            acb_ode_sum_divconquer(sum + i, nterms, NULL, NULL, prec[i]);

        if (sum[0].n != nterms || sum[1].n != nterms)
        {
            flint_printf("FAIL (nterms)\n\n");
            goto fail;
        }

        acb_poly_struct * res = _acb_poly_vec_init(dop_len - 1);
        _acb_poly_vec_fit_length(res, dop_len - 1, nterms);

        for (slong m = 0; m < nsols; m++)
        {
            acb_ode_sol_struct * sol0 = sum[0].sol + m;
            acb_ode_sol_struct * sol1 = sum[1].sol + m;

            _acb_ode_apply_diffop_polmul(res, 0, dop, dop_len, group->leader, 0,
                                         sol0->series, 0, nterms, sol0->nlogs,
                                         0, nterms, prec[0]);
            for (slong k = 0; k < sol0->nlogs; k++)
            {
                for (slong j = 0; j < nterms; j++)
                {
                    if (!acb_contains_zero(res[k].coeffs + j))
                    {
                        flint_printf("FAIL (residual)\n\n");
                        flint_printf("m = %wd, k = %wd, j = %wd\n", m, k, j);
                        goto fail;
                    }
                }
            }

            if (sum[1].flags & ACB_ODE_APPROX)
                continue;

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

        _acb_poly_vec_clear(res, dop_len - 1);
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
        flint_printf("nterms = %wd, n = %wd/%wd, flags = %wu, "
                     "prec = %wd/%wd\n", nterms, sum[0].n, sum[1].n,
                     sum[1].flags, prec[0], prec[1]);
        flint_abort();
    }

    TEST_FUNCTION_END(state);
}
