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
#include "acb_mat.h"
#include "acb_poly.h"
#include "acb_ode.h"
#include "arb_poly.h"

TEST_FUNCTION_START(acb_ode_sum_update_tail_bounds, state)
{
    for (slong iter = 0; iter < 50 * flint_test_multiplier(); iter++)
    {
        acb_ode_exponents_t expos;
        acb_ode_bound_t bound;
        acb_ode_group_bound_t gbound;
        acb_ode_sum_struct sum[2];
        acb_poly_t val, refval;
        mag_t cvrad, rad, mag;
        ulong seed1, seed2;

        flint_rand_get_seed(&seed1, &seed2, state);

        slong dop_len = 2 + n_randint(state, 4);
        slong lcdeg = n_randint(state, 3);
        slong prec = 2 + n_randint(state, 100);

        acb_poly_struct * dop = _acb_poly_vec_init(dop_len);
        acb_ptr lcroots = _acb_vec_init(lcdeg);
        acb_ode_exponents_init(expos);
        acb_ode_bound_init(bound);
        acb_ode_group_bound_init(gbound);
        acb_poly_init(val);
        acb_poly_init(refval);
        mag_init(cvrad);
        mag_init(rad);
        mag_init(mag);

        /* bound_precompute does not support an irregular singular point at the
           origin */
        do
        {
            acb_ode_randtest_acb(dop, lcroots, expos, state,
                                 n_randint(state, 11), lcdeg,
                                 1 + n_randint(state, 6), dop_len, prec);
            _acb_vec_get_mag_lower(cvrad, lcroots, lcdeg);
        }
        while (mag_is_zero(cvrad));
        slong dop_clen = _acb_poly_vec_length(dop, dop_len);

        slong g = n_randint(state, expos->ngroups);
        acb_ode_group_struct * group = expos->groups + g;
        int highest = n_randint(state, 2);
        slong maxsols = highest ? group->nshifts : acb_ode_group_length(group);
        slong nsols = n_randint(state, maxsols + 1);
        slong npts = 1 + n_randint(state, 2);
        slong nder = 1 + n_randint(state, 3);
        slong nterms = n_randint(state, 40);
        slong pol_part_len = 1 + n_randint(state, dop_clen + 4);

        acb_ptr pts = _acb_vec_init(npts);
        mag_mul_2exp_si(rad, cvrad, -4);
        for (slong p = 0; p < npts; p++)
        {
            acb_randtest(pts + p, state, prec, 2);
            acb_get_mag(mag, pts + p);
            mag_div(mag, mag, rad);
            if (mag_cmp_2exp_si(mag, 0) > 0)
                acb_mul_2exp_si(pts + p, pts + p,
                                -fmpz_get_si(MAG_EXPREF(mag)));
        }

        acb_ode_bound_precompute(bound, dop, dop_len, lcroots, pol_part_len,
                                 MAG_BITS);
        acb_ode_group_bound_precompute(gbound, dop, dop_len, expos, g, bound,
                                       bound->prec);

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

        acb_ode_sum_divconquer(sum + 0, 1000, bound, gbound, prec);

        acb_ode_sum_precompute(sum + 1);
        sum[1].wp = sum[1].sum_wp = prec;
        for (slong p = 0; p < npts; p++)
            sum[1].shifted_sums[p] = 0;
        acb_ode_sum_fit_length(sum + 1, nterms + dop_clen);

        for (slong n = 0; n < nterms; n++)
        {
            _acb_ode_sum_forward_1(sum + 1);
            _acb_ode_sum_update_residuals(sum + 1, n, n + 1, n + dop_clen);
            _acb_ode_sum_update_tail_bounds(sum + 1, bound, gbound, MAG_BITS);

            for (slong m = 0; m < nsols; m++)
            {
                for (slong p = 0; p < npts; p++)
                {
                    acb_ode_sol_jet(val, group->leader, sum[1].sol + m, p,
                                    pts + p, nder, 1, prec);
                    acb_ode_sol_jet(refval, group->leader, sum[0].sol + m, p,
                                    pts + p, nder, 1, prec);

                    if (!acb_poly_overlaps(val, refval))
                    {
                        flint_printf("FAIL (jet)\n\n");
                        flint_printf("n = %wd, m = %wd, p = %wd\n", n, m, p);
                        flint_printf("tb = %{mag*}\n", sum[1].sol[m].tb, nder);
                        flint_printf("val = %{acb_poly}\n", val);
                        flint_printf("ref = %{acb_poly}\n", refval);
                        goto fail;
                    }
                }
            }
        }

        acb_ode_sum_clear(sum + 0);
        acb_ode_sum_clear(sum + 1);
        _acb_vec_clear(pts, npts);
        _acb_vec_clear(lcroots, lcdeg);
        _acb_poly_vec_clear(dop, dop_len);
        acb_ode_group_bound_clear(gbound);
        acb_ode_bound_clear(bound);
        acb_ode_exponents_clear(expos);
        acb_poly_clear(val);
        acb_poly_clear(refval);
        mag_clear(cvrad);
        mag_clear(rad);
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
        flint_printf("pol_part_len = %wd, prec = %wd\n", pol_part_len, prec);
        flint_abort();
    }

    /* At zero, tail bounds are bounds on the series coefficients themselves.
       For dop = Q(θ)·lc, the numerators vanish, so that the bounds are close to
       tight, especially when the roots of lc form a small angle.

       Test generated by Claude Opus 5.5, not reviewed in detail. */
    for (slong iter = 0; iter < 20 * flint_test_multiplier(); iter++)
    {
        acb_ode_bound_t bound;
        acb_ode_group_bound_t gbound;
        acb_ode_sum_struct sum[3];
        acb_poly_t ind, lc, shifted;
        arb_poly_t nres_maj, jet_tb;
        acb_t c;
        arb_t a, rad;
        mag_t cvrad, mag, tb;
        ulong seed1, seed2;

        flint_rand_get_seed(&seed1, &seed2, state);

        slong prec = 128;
        slong N = 14 + n_randint(state, 20);
        slong s = N + 2 + n_randint(state, N / 2);
        slong nder = s + N;
        slong n0 = 1 + n_randint(state, s);
        slong pol_part_len = 1 + n_randint(state, 7);

        acb_ode_shift_struct shifts[2] = {
            { .n = 0, .mult = 1 },
            { .n = s, .mult = 1 + n_randint(state, 2) },
        };
        acb_ode_group_struct group[1] = {{ .shifts = shifts, .nshifts = 2 }};
        acb_ode_exponents_struct expos[1] = {{ .groups = group, .ngroups = 1 }};
        slong dop_len = 2 + shifts[1].mult;

        acb_poly_struct * dop = _acb_poly_vec_init(dop_len);
        acb_ptr lcroots = _acb_vec_init(2);
        acb_init(group->leader);
        acb_ode_bound_init(bound);
        acb_ode_group_bound_init(gbound);
        acb_poly_init(ind);
        acb_poly_init(lc);
        acb_poly_init(shifted);
        arb_poly_init(nres_maj);
        arb_poly_init(jet_tb);
        acb_init(c);
        arb_init(a);
        arb_init(rad);
        mag_init(cvrad);
        mag_init(mag);
        mag_init(tb);

        acb_one(lcroots + 0);
        acb_unit_root(lcroots + 1, N, prec);
        acb_poly_product_roots(lc, lcroots, 2, prec);
        slong dop_clen = lc->length;

        acb_ode_indicial_polynomial_from_exponents(ind, expos, prec);
        for (slong j = 0; j < dop_clen; j++)
        {
            acb_set_si(c, j);
            acb_poly_taylor_shift(shifted, ind, c, prec);
            acb_poly_scalar_mul(shifted, shifted, lc->coeffs + j, prec);
            for (slong i = 0; i < dop_len; i++)
                acb_poly_set_coeff_acb(dop + i, j, shifted->coeffs + i);
        }

        acb_ode_bound_precompute(bound, dop, dop_len, lcroots, pol_part_len,
                                 MAG_BITS);
        acb_ode_group_bound_precompute(gbound, dop, dop_len, expos, 0, bound,
                                       bound->prec);

        for (int i = 0; i < 3; i++)
        {
            acb_ode_sum_init(sum + i, dop_len, 0, 1, nder);
            acb_ode_sum_set_diffop(sum + i, dop, dop_len, cvrad);
            acb_ode_sum_set_group(sum + i, group);
            acb_ode_sum_set_ini_highest(sum + i);
            acb_ode_sum_precompute(sum + i);
            sum[i].flags = ACB_ODE_WANT_SERIES;
            sum[i].wp = sum[i].sum_wp = prec;
            acb_ode_sum_fit_length(sum + i, nder + dop_clen);
        }

        acb_ode_sum_divconquer(sum + 0, nder, NULL, NULL, prec);

        for (slong n = 0; n < nder; n++)
        {
            _acb_ode_sum_forward_1(sum + 1);
            _acb_ode_sum_update_residuals(sum + 1, n, n + 1, n + dop_clen);
            _acb_ode_sum_update_tail_bounds(sum + 1, bound, gbound, MAG_BITS);

            for (slong k = 0; k < sum[0].sol->nlogs; k++)
            {
                for (slong j = n + 1; j < nder; j++)
                {
                    acb_srcptr cj = sum[0].sol->series[k].coeffs + j;
                    acb_get_mag_lower(mag, cj);
                    if (mag_cmp(mag, sum[1].sol->tb + j) > 0)
                    {
                        flint_printf("FAIL (coefficient vs tail bound)\n\n");
                        flint_printf("n = %wd, k = %wd, j = %wd\n",
                                     n + 1, k, j);
                        flint_printf("c = %{acb}\n", cj);
                        flint_printf("tb = %{mag}\n", sum[1].sol->tb + j);
                        goto fail2;
                    }
                }
            }
        }

        /* tail_bound_jet, for a solution such that lc·y_{n0:} has no term at
           the positions of the initial values */

        for (slong n = 0; n < nder; n++)
        {
            acb_poly_struct * series = sum[2].sol->series;

            if (n == n0)
            {
                for (slong i = 0; i < dop_clen - 1; i++)
                {
                    if (acb_ode_group_multiplicity(group, n + i))
                        continue;
                    acb_set_si(c, n + i);
                    acb_poly_evaluate(c, sum[2].ind, c, prec);
                    acb_div(c, series->coeffs + n + i, c, prec);
                    acb_abs(a, c, prec);
                    arb_poly_set_coeff_arb(nres_maj, i, a);
                }
                acb_ode_tail_bound_jet(jet_tb, group, bound, gbound, n, 1,
                                       nres_maj, rad, nder, MAG_BITS);
            }

            if (n == s && n > n0)
            {
                for (slong k = 0; k < shifts[1].mult; k++)
                {
                    acb_zero(c);
                    for (slong i = 1; i < dop_clen && n - i >= n0; i++)
                        acb_addmul(c, lc->coeffs + i,
                                   series[k].coeffs + n - i, prec);
                    acb_div(c, c, lc->coeffs, prec);
                    acb_neg(acb_mat_entry(sum[2].sol->extini, 1, k), c);
                }
            }

            _acb_ode_sum_forward_1(sum + 2);
            _acb_ode_sum_update_residuals(sum + 2, n, n + 1, n + dop_clen);
        }

        for (slong k = 0; k < sum[2].sol->nlogs; k++)
        {
            for (slong j = n0; j < nder; j++)
            {
                acb_srcptr cj = sum[2].sol->series[k].coeffs + j;
                acb_get_mag_lower(mag, cj);
                arb_poly_get_coeff_arb(a, jet_tb, j);
                arb_get_mag(tb, a);
                if (mag_cmp(mag, tb) > 0)
                {
                    flint_printf("FAIL (coefficient vs tail_bound_jet)\n\n");
                    flint_printf("n0 = %wd, k = %wd, j = %wd\n", n0, k, j);
                    flint_printf("c = %{acb}\n", cj);
                    flint_printf("tb = %{mag}\n", tb);
                    goto fail2;
                }
            }
        }

        for (int i = 0; i < 3; i++)
            acb_ode_sum_clear(sum + i);
        _acb_vec_clear(lcroots, 2);
        _acb_poly_vec_clear(dop, dop_len);
        acb_clear(group->leader);
        acb_ode_group_bound_clear(gbound);
        acb_ode_bound_clear(bound);
        acb_poly_clear(ind);
        acb_poly_clear(lc);
        acb_poly_clear(shifted);
        arb_poly_clear(nres_maj);
        arb_poly_clear(jet_tb);
        acb_clear(c);
        arb_clear(a);
        arb_clear(rad);
        mag_clear(cvrad);
        mag_clear(mag);
        mag_clear(tb);

        continue;

fail2:
        flint_printf("flint_rand_set_seed(state, 0x%wx, 0x%wx);\n",
                     seed1, seed2);
        for (slong i = 0; i < dop_len; i++)
            flint_printf("dop[%wd] = %{acb_poly}\n", i, dop + i);
        flint_printf("expos = ");
        acb_ode_exponents_println(expos);
        flint_printf("N = %wd, pol_part_len = %wd\n", N, pol_part_len);
        flint_abort();
    }

    TEST_FUNCTION_END(state);
}
