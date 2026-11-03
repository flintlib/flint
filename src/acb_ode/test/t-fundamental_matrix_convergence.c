/*
    Copyright (C) 2026 Marc Mezzarobba

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/* generated using Claude Code (claude-opus-5-5) */

#include "test_helpers.h"
#include "acb.h"
#include "acb_mat.h"
#include "acb_ode.h"
#include "fmpz.h"
#include "gr_poly.h"
#include "gr_ore_poly.h"

/* Test that for exact inputs, the radii of acb_ode_fundamental_matrix results
   tend to zero as prec tends to infinity. */

#define MAX_RETRIES 4

TEST_FUNCTION_START(acb_ode_fundamental_matrix_convergence, state)
{
    gr_ctx_t ZZ;
    int flags = 0 * GR_TEST_VERBOSE;
    slong nretries[2][MAX_RETRIES + 2] = {{0}}, nskip[2] = {0, 0};
    slong nfail = 0;

    gr_ctx_init_fmpz(ZZ);

    slong niter = 100 * flint_test_multiplier();
    for (slong iter = 0; iter < niter; iter++)
    {
        gr_ctx_t CC, Pol;
        gr_ore_poly_ctx_t Dop;
        gr_ore_poly_t dop;
        acb_ode_exponents_t expos;
        acb_ptr lcroots;
        slong lcdeg;
        mag_t rho, m;
        int status = GR_SUCCESS;

        ulong seed1, seed2;
        flint_rand_get_seed(&seed1, &seed2, state);

        slong prec = 32 + n_randint(state, 300);
        slong target = 3 * prec / 4;
        slong hp = (prec << MAX_RETRIES) + 8;

        int exact = n_randint(state, 2);
        ore_algebra_t alg = (!exact || n_randint(state, 2))
                          ? ORE_ALGEBRA_EULER_DERIVATIVE
                          : ORE_ALGEBRA_DERIVATIVE;
        slong dop_len = 2 + n_randint(state, 4);
        slong order = dop_len - 1;

        gr_ctx_init_complex_acb(CC, hp);
        gr_ctx_init_gr_poly(Pol, exact ? ZZ : CC);
        gr_ore_poly_ctx_init(Dop, Pol, 0, alg);
        gr_ore_poly_init(dop, Dop);
        acb_ode_exponents_init(expos);
        mag_init(rho);
        mag_init(m);

        gr_ore_poly_fit_length(dop, dop_len, Dop);

        if (exact)
        {
            /* small coefficients to avoid solutions of huge magnitude */
            slong clen = 1 + n_randint(state, 5);
            for (slong i = 0; i < dop_len; i++)
            {
                gr_poly_struct * c = gr_ore_poly_coeff_ptr(dop, i, Dop);
                gr_poly_fit_length(c, clen, ZZ);
                for (slong j = 0; j < clen; j++)
                    fmpz_randtest(gr_poly_coeff_ptr(c, j, ZZ), state, 4);
                if (i == dop_len - 1 && fmpz_is_zero(c->coeffs))
                    fmpz_one(c->coeffs);
                _gr_poly_set_length(c, clen, ZZ);
                _gr_poly_normalise(c, ZZ);
            }
            _gr_ore_poly_set_length(dop, dop_len, Dop);

            lcdeg = gr_poly_length(gr_ore_poly_coeff_ptr(dop, order, Dop), ZZ) - 1;
            lcroots = _acb_vec_init(lcdeg);
            status |= acb_ode_lc_roots(lcroots, dop, Dop, 0, CC);
            status |= acb_ode_exponents(expos, dop, Dop, CC);
        }
        else
        {
            /* todo: switch to a (yet to be written) variant of
               acb_ode_randtest_acb that generates exact operators with an
               interesting exponent structure, making the midpoint and
               widening hacks below unnecessary */
            slong clen = 1 + n_randint(state, 6);
            lcdeg = n_randint(state, clen);
            lcroots = _acb_vec_init(lcdeg);
            acb_ode_randtest_acb(dop->coeffs, lcroots, expos, state, 10,
                                 lcdeg, clen, dop_len, hp);
            _gr_ore_poly_set_length(dop, dop_len, Dop);
            for (slong i = 0; i < dop_len; i++)
            {
                acb_poly_struct * c = gr_ore_poly_coeff_ptr(dop, i, Dop);
                for (slong j = 0; j < c->length; j++)
                    acb_get_mid(c->coeffs + j, c->coeffs + j);
            }
        }

        /* large exponents make the computation slow; imprecise singularities
           or exponents limit the accuracy of the result */
        for (slong g = 0; g < expos->ngroups; g++)
        {
            acb_ode_group_struct * grp = expos->groups + g;
            acb_get_mag(m, grp->leader);
            if (mag_cmp_2exp_si(m, 5) > 0
                || grp->shifts[grp->nshifts - 1].n > 5)
                status = GR_UNABLE;
            mag_max(m, arb_radref(acb_realref(grp->leader)),
                       arb_radref(acb_imagref(grp->leader)));
            if (mag_cmp_2exp_si(m, -hp / 2) > 0)
                status = GR_UNABLE;
        }
        for (slong i = 0; i < lcdeg; i++)
            if (acb_rel_accuracy_bits(lcroots + i) < 3 * hp / 4)
                status = GR_UNABLE;

        /* the roots of the leading coefficient of the midpoint operator
           differ from the ones returned by acb_ode_randtest_acb by up to about
           eps^(1/mult)

           WARNING: this is only a heuristic; the widened lcroots are not
           guaranteed to be rigorous enclosures, so neither are the output
           balls (this test only checks their accuracy) */
        if (!exact)
            for (slong i = 0; i < lcdeg; i++)
            {
                acb_get_mag(m, lcroots + i);
                mag_mul_2exp_si(m, m, -hp / lcdeg);
                acb_add_error_mag(lcroots + i, m);
            }

        /* exact evaluation points at most 2^-mult of the way to the nearest
           singularity, where mult is the largest multiplicity of a root of
           the leading coefficient (convergence is slow otherwise) */

        slong mult = 1;
        for (slong i = 0; i < lcdeg; i++)
        {
            slong m = 1;
            for (slong j = i + 1; j < lcdeg; j++)
                m += acb_overlaps(lcroots + i, lcroots + j);
            mult = FLINT_MAX(mult, m);
        }

        _acb_vec_get_mag_lower(rho, lcroots, lcdeg);
        if (mag_is_inf(rho))
            mag_one(rho);

        slong npts = 1 + n_randint(state, 2);
        acb_ptr pts = _acb_vec_init(npts);
        for (slong p = 0; p < npts; p++)
        {
            /* avoid evaluating at a singular point */
            do
            {
                acb_randtest_precise(pts + p, state, prec, 2);
                acb_get_mid(pts + p, pts + p);
            }
            while (alg == ORE_ALGEBRA_EULER_DERIVATIVE && acb_is_zero(pts + p));

            acb_get_mag(m, pts + p);
            mag_div(m, m, rho);
            while (mag_cmp_2exp_si(m, -mult) > 0)
            {
                acb_mul_2exp_si(pts + p, pts + p, -1);
                mag_mul_2exp_si(m, m, -1);
            }
        }

        slong nrows = 1 + n_randint(state, order + 1);
        acb_mat_struct * mat = flint_malloc(sizeof(acb_mat_struct) * npts);
        for (slong p = 0; p < npts; p++)
            acb_mat_init(mat + p, nrows, order);

        if (status != GR_SUCCESS)
        {
            nskip[exact]++;
            goto cleanup;
        }

        slong retries, acc = -WORD_MAX;
        for (retries = 0; retries <= MAX_RETRIES; retries++)
        {
            status = acb_ode_fundamental_matrix_vec(mat, dop, Dop, expos,
                                                    lcroots, pts, npts,
                                                    ACB_ODE_BASIS_ECHELON, 0,
                                                    prec << retries);

            acc = WORD_MAX;
            for (slong p = 0; p < npts; p++)
                for (slong i = 0; i < nrows; i++)
                    for (slong j = 0; j < order; j++)
                    {
                        acb_srcptr c = acb_mat_entry(mat + p, i, j);
                        slong a = acb_rel_accuracy_bits(c);
                        mag_max(m, arb_radref(acb_realref(c)),
                                   arb_radref(acb_imagref(c)));
                        if (mag_is_zero(m))
                            a = WORD_MAX;
                        else if (mag_is_finite(m))
                            a = FLINT_MAX(a, -fmpz_get_si(MAG_EXPREF(m)));
                        acc = FLINT_MIN(acc, a);
                    }

            if (status == GR_SUCCESS && acc >= target)
                break;
        }

        nretries[exact][retries]++;

        if (retries > MAX_RETRIES)
        {
            nfail++;

            /*
            if (flags & GR_TEST_VERBOSE)
            {
                flint_printf("flint_rand_set_seed(state, 0x%wx, 0x%wx);\n",
                             seed1, seed2);
                flint_printf("prec = %wd, accuracy = %wd, status = %d, "
                             "alg = %d, nrows = %wd\n",
                             prec, acc, status, alg, nrows);
                flint_printf("dop = %{gr}\n", dop, Dop);
                flint_printf("expos = ");
                acb_ode_exponents_println(expos);
                flint_printf("lcroots = %{acb*}\n", lcroots, lcdeg);
                flint_printf("pts = %{acb*}\n", pts, npts);
                flint_printf("\n");
            }
            */

            if (nfail * 20 > niter)
            {
                flint_printf("FAIL: insufficient accuracy after %d retries\n\n",
                             MAX_RETRIES);
                flint_printf("%wd failures after %wd of %wd iterations\n",
                             nfail, iter + 1, niter);
                flint_printf("flint_rand_set_seed(state, 0x%wx, 0x%wx);\n",
                             seed1, seed2);
                flint_abort();
            }
        }

cleanup:
        for (slong p = 0; p < npts; p++)
            acb_mat_clear(mat + p);
        flint_free(mat);
        _acb_vec_clear(pts, npts);
        _acb_vec_clear(lcroots, lcdeg);
        mag_clear(rho);
        mag_clear(m);
        acb_ode_exponents_clear(expos);
        gr_ore_poly_clear(dop, Dop);
        gr_ctx_clear(Dop);
        gr_ctx_clear(Pol);
        gr_ctx_clear(CC);
    }

    if (flags & GR_TEST_VERBOSE)
    {
        for (int exact = 0; exact < 2; exact++)
        {
            flint_printf("%s operators: skipped %wd, retries:",
                         exact ? "exact" : "inexact", nskip[exact]);
            for (slong r = 0; r <= MAX_RETRIES; r++)
                flint_printf(" %wd", nretries[exact][r]);
            flint_printf(", failed %wd\n", nretries[exact][MAX_RETRIES + 1]);
        }
    }

    gr_ctx_clear(ZZ);

    TEST_FUNCTION_END(state);
}
