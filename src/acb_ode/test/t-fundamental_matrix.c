/*
    Copyright (C) 2026 Marc Mezzarobba

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "acb_types.h"
#include "acb_mat.h"
#include "acb_ode.h"
#include "acb_poly.h"
#include "acb.h"
#include "gr_poly.h"

/* originally human-written, revised using Claude Code (claude-opus-5-5) */

FLINT_DLL extern gr_static_method_table _ca_methods,  _qqbar_methods;

TEST_GR_FUNCTION_START(acb_ode_fundamental_matrix, state, count_success, count_domain, count_unable)
{
    gr_ctx_t ZZ, QQ, ZZi, QQbar, CaCC;

    gr_ctx_init_fmpz(ZZ);
    gr_ctx_init_fmpq(QQ);
    gr_ctx_init_fmpzi(ZZi);
    gr_ctx_init_complex_qqbar(QQbar);
    gr_ctx_init_complex_ca(CaCC);

    for (slong iter = 0; iter < 100 * flint_test_multiplier(); iter++)
    {
        ulong seed1, seed2;
        flint_rand_get_seed(&seed1, &seed2, state);

        int status = GR_SUCCESS;
        /* status of the optional parts of the input (used for the second
           comparison only) */
        int status_opt = GR_SUCCESS;
        int expect_success = 0;

        gr_ctx_struct * Cst;
        gr_ctx_t CC, Pol;
        gr_ore_poly_ctx_t Dop;
        gr_ptr pert;
        gr_ore_poly_t dop, dop1;
        acb_t r;
        acb_ode_exponents_t expos;
        acb_ode_exponents_struct * expos1;
        acb_ptr lcroots = NULL;
        slong lcdeg1;
        mag_t refmag, mag;

        acb_init(r);
        acb_ode_exponents_init(expos);
        mag_init(refmag);
        mag_init(mag);

        /* target precision */

        slong prec = 2 + n_randint(state, 100);
        slong prec1 = 2 + n_randint(state, prec);  /* prec1 < prec */
        gr_ctx_init_complex_acb(CC, prec + 10);

        /* Ore polynomial ring (often but not always differential case) */

        if (!n_randint(state, 16))
            gr_ore_poly_ctx_init_randtest2(Pol, Dop, state);
        else
        {
            slong rnd = n_randint(state, 12);
            if (rnd == 0)
                gr_ctx_init_random_poly(Pol, state);
            else if (rnd == 1)
                gr_ctx_init_gr_poly(Pol, QQ);
            else if (rnd == 2)
                gr_ctx_init_gr_poly(Pol, ZZi);
            else if (rnd == 3)
                gr_ctx_init_gr_poly(Pol, QQbar);
            else if (rnd == 4)
                gr_ctx_init_gr_poly(Pol, CaCC);
            else if (rnd <= 8)
            {
                gr_ctx_init_gr_poly(Pol, ZZ);
                expect_success = 1;
            }
            else
                gr_ctx_init_gr_poly(Pol, CC);

            ore_algebra_t alg = n_randint(state, 4)
                ? ORE_ALGEBRA_EULER_DERIVATIVE
                : ORE_ALGEBRA_DERIVATIVE;
            gr_ore_poly_ctx_init(Dop, Pol, 0, alg);
        }

        /* operator */

        GR_TMP_INIT(pert, Pol);
        gr_ore_poly_init(dop,  Dop);
        gr_ore_poly_init(dop1, Dop);

        slong dop_len_bound = 8;
        if (Pol->which_ring == GR_CTX_GR_POLY
            && (POLYNOMIAL_ELEM_CTX(Pol)->methods == _ca_methods
                || POLYNOMIAL_ELEM_CTX(Pol)->methods ==  _qqbar_methods))
            dop_len_bound = 4;
        if (n_randint(state, 4))
            dop_len_bound = FLINT_MIN(dop_len_bound, 5);
        slong dop_len = n_randint(state, 8)
                      ? 2 + n_randint(state, dop_len_bound - 2)
                      : n_randint(state, 2);

        slong disp = n_randint(state, 16) ? 10 : 1000;

        if (GR_ORE_POLY_CTX(Dop)->which_algebra == ORE_ALGEBRA_EULER_DERIVATIVE
            && Pol->which_ring == GR_CTX_GR_POLY
            /* todo: acb_ode_randtest for exact coefficients */
            && POLYNOMIAL_ELEM_CTX(Pol)->which_ring == GR_CTX_CC_ACB
            && n_randint(state, 8))
        {
            /* todo: increase once sum_forward_divconquer is able to stop in the
               middle of a block */
            slong clen = 1 + n_randint(state, 10);
            slong lcdeg = n_randint(state, clen);
            lcdeg1 = lcdeg + n_randint(state, 4);
            lcroots = _acb_vec_init(lcdeg1);

            gr_ore_poly_fit_length(dop, dop_len, Dop);
            acb_ode_randtest_acb(dop->coeffs, lcroots, expos, state, disp,
                                 lcdeg, clen, dop_len, prec * 4);
            _gr_ore_poly_set_length(dop, dop_len, Dop);

            for (slong i = lcdeg; i < lcdeg1; i++)  /* extra lcroots for dop1 */
            {
                if (n_randint(state, 8))
                    acb_randtest_precise(lcroots + i, state, prec * 4, 3);
                else
                    acb_randtest(lcroots + i, state, prec * 4, 3);
                if (acb_is_zero(lcroots + i))
                    acb_one(lcroots + i);
            }
            acb_poly_product_roots(pert, lcroots + lcdeg, lcdeg1 - lcdeg,
                                   prec * 4);

            /* depends on the fact that lcroots is provided */
            expect_success = 1;
        }
        else
        {
            int ordinary = 0;

            status |= gr_ore_poly_randtest(dop, state, dop_len, Dop);
            /* most of the time, force the order to be dop_len - 1 and make
               zero a regular point */
            if (Pol->which_ring == GR_CTX_GR_POLY && dop_len > 0
                && n_randint(state, 8))
            {
                gr_ore_poly_fit_length(dop, dop_len, Dop);
                gr_ptr lc = gr_ore_poly_coeff_ptr(dop, dop_len - 1, Dop);
                status |= gr_randtest_not_zero(lc, state, Pol);
                status |= gr_randtest_not_zero(
                        gr_poly_coeff_ptr(lc, 0, POLYNOMIAL_ELEM_CTX(Pol)),
                        state,
                        POLYNOMIAL_ELEM_CTX(Pol));
                _gr_ore_poly_set_length(dop, dop_len, Dop);
                _gr_ore_poly_normalise(dop, Dop);
                ordinary = (status == GR_SUCCESS
                            && GR_ORE_POLY_CTX(Dop)->which_algebra
                               == ORE_ALGEBRA_DERIVATIVE);
            }

            /* precompute the exponents and pass them to the main function so
               that exponent groups are ordered consistently across calls */
            if (ordinary && n_randint(state, 2))
                acb_ode_exponents_ordinary(expos, dop_len - 1);
            else
                status |= acb_ode_exponents(expos, dop, Dop, CC);

            /* skip large exponents (slow)
               todo: relax this once we have better support for this case */
            if (status == GR_SUCCESS)
                for (slong g = 0; g < expos->ngroups; g++)
                {
                    acb_ode_group_struct * grp = expos->groups + g;
                    acb_get_mag(mag, grp->leader);
                    if (mag_cmp_2exp_si(mag, 10) > 0)
                        status = GR_UNABLE;
                    for (slong s = 0; s < grp->nshifts; s++)
                        if (grp->shifts[s].n > disp)
                            status = GR_UNABLE;
                }

            status_opt |= gr_randtest_not_zero(pert, state, Pol);
            /* ensure however that dop1 has no extra singular point at zero in
               interesting cases */
            if (Pol->which_ring == GR_CTX_GR_POLY)
                status_opt |= gr_randtest_not_zero(
                        gr_poly_coeff_ptr(pert, 0, POLYNOMIAL_ELEM_CTX(Pol)),
                        state,
                        POLYNOMIAL_ELEM_CTX(Pol));
        }

        slong order = dop->length - 1;

        if (expos->ngroups == 1 && n_randint(state, 4))
            expos1 = NULL;
        else
            expos1 = expos;

        status_opt |= gr_ore_poly_other_mul(dop1, pert, Pol, dop, Dop);

        /* evaluation points */

        mag_inf(refmag);
        if (lcroots)
            _acb_vec_get_mag_lower(refmag, lcroots, lcdeg1);
        else if (status_opt == GR_SUCCESS && Pol->which_ring == GR_CTX_GR_POLY
                 && dop1->length > 0)
        {
            slong lcdeg = gr_poly_length(gr_ore_poly_coeff_ptr(
                    dop1, dop1->length - 1, Dop), POLYNOMIAL_ELEM_CTX(Pol)) - 1;
            acb_ptr sing = _acb_vec_init(lcdeg);
            if (acb_ode_lc_roots(sing, dop1, Dop, 0, CC) == GR_SUCCESS)
                _acb_vec_get_mag_lower(refmag, sing, lcdeg);
            _acb_vec_clear(sing, lcdeg);
        }
        if (mag_is_inf(refmag))
            mag_set_ui(refmag, 3);

        slong npts = n_randint(state, 8) ? 1 + n_randint(state, 2)
                                         : n_randint(state, 20);
        acb_ptr pts = _acb_vec_init(npts * 2);
        acb_ptr pts1 = pts + npts;
        for (slong p = 0; p < npts; p++)
        {
            acb_randtest_precise(pts + p, state, prec + 30, 3);
            acb_get_mag(mag, pts + p);
            mag_div(mag, mag, refmag);
            slong s = FLINT_MAX(0., mag_get_d_log2_approx(mag) + 1.5);
            acb_mul_2exp_si(pts + p, pts + p, -s);

            if (n_randint(state, 16))
                acb_set(pts1 + p, pts + p);
            else
            {
                acb_randtest_special(r, state, prec, 3);
                acb_add(pts1 + p, pts  + p, r, prec);
                acb_sub(pts1 + p, pts1 + p, r, prec);
            }
        }

        acb_ode_basis_t basis = n_randint(state, ACB_ODE_NUM_BASES);

        acb_mat_struct * mat = flint_malloc(sizeof(acb_mat_struct) * npts * 3);
        acb_mat_struct * mat1 = mat + npts;
        acb_mat_struct * mat2 = mat1 + npts;
        for (slong p = 0; p < npts; p++)
        {
            slong nrows = (order >= 0 && n_randint(state, 4))
                        ? 1 + n_randint(state, order + 1)
                        : n_randint(state, order + 2);
            acb_mat_init(mat + p,  nrows, order);
            acb_mat_init(mat1 + p, nrows, order);
            acb_mat_init(mat2 + p, nrows, order);
        }

        /* test */

        if (status != GR_SUCCESS)
        {
            status = GR_SUCCESS;  /* don't count test setup failures */
            goto end;
        }

        if (expect_success
            && (gr_ore_poly_ctx_over_gr_poly_base_ptrs(&Cst, NULL, Dop) != GR_SUCCESS
                || dop->length == 0
                || gr_is_zero(
                        gr_poly_coeff_ptr(
                                gr_ore_poly_coeff_ptr(dop, dop->length - 1, Dop),
                                0, Cst),
                        Cst) != T_FALSE))
        {
            expect_success = 0;
        }

        status |= acb_ode_fundamental_matrix_vec(mat, dop, Dop, expos, lcroots,
                                                 pts, npts, basis, 0, prec + 30);

        if (expect_success && status == GR_UNABLE)
        {
            flint_printf("FAIL: unexpected failure\n\n");
            flint_printf("flint_rand_set_seed(state, 0x%wx, 0x%wx);\n", seed1, seed2);
            flint_printf("Dop = %{gr_ctx}\n", Dop);
            flint_printf("dop = %{gr}\n", dop, Dop);
            flint_printf("expos = "); acb_ode_exponents_println(expos);
            flint_printf("pts = %{acb*}\n", pts, npts);
            flint_printf("prec = %wd\n", prec);
            flint_abort();
        }

        if (status != GR_SUCCESS)
            goto end;

        status |= acb_ode_fundamental_matrix_vec(mat1, dop, Dop, expos1, lcroots,
                                                 pts1, npts, basis, 0, prec1);

        if (status == GR_SUCCESS)
        {
            for (slong p = 0; p < npts; p++)
            {
                if (!acb_mat_overlaps(mat + p, mat1 + p))
                {
                    flint_printf("FAIL: overlap 1\n\n");
                    flint_printf("flint_rand_set_seed(state, 0x%wx, 0x%wx);\n", seed1, seed2);
                    flint_printf("Dop = %{gr_ctx}\n", Dop);
                    flint_printf("dop = %{gr}\n", dop, Dop);
                    flint_printf("expos = "); acb_ode_exponents_println(expos);
                    flint_printf("prec = %wd, prec1 = %wd\n", prec, prec1);
                    flint_printf("pts[%wd] = %{acb}\n", p, pts + p);
                    flint_printf("pts1[%wd] = %{acb}\n", p, pts1 + p);
                    flint_printf("mat = %{acb_mat}\n", mat + p);
                    flint_printf("mat1 = %{acb_mat}\n", mat1 + p);
                    flint_abort();
                }
            }
        }

        if (status_opt != GR_SUCCESS
            || (Pol->which_ring == GR_CTX_GR_POLY
                && POLYNOMIAL_ELEM_CTX(Pol)->methods == _ca_methods))
            goto end;

        status |= acb_ode_fundamental_matrix_vec(mat2, dop1, Dop, expos, lcroots,
                                                 pts, npts, basis, 0, prec);

        if (status == GR_SUCCESS)
        {
            for (slong p = 0; p < npts; p++)
            {
                if (!acb_mat_overlaps(mat + p, mat2 + p))
                {
                    flint_printf("FAIL: overlap 2\n\n");
                    flint_printf("flint_rand_set_seed(state, 0x%wx, 0x%wx);\n", seed1, seed2);
                    flint_printf("Dop = %{gr_ctx}\n", Dop);
                    flint_printf("dop = %{gr}\n", dop, Dop);
                    flint_printf("dop1 = %{gr}\n", dop1, Dop);
                    flint_printf("expos = "); acb_ode_exponents_println(expos);
                    flint_printf("prec = %wd\n", prec);
                    flint_printf("pts[%wd] = %{acb}\n", p, pts + p);
                    flint_printf("mat = %{acb_mat}\n", mat + p);
                    flint_printf("mat2 = %{acb_mat}\n", mat2 + p);
                    flint_abort();
                }
            }
        }

end:
        count_success += (status == GR_SUCCESS);
        count_domain += ((status & GR_DOMAIN) != 0);
        count_unable += ((status & GR_UNABLE) != 0);

        for (slong p = 0; p < npts; p++)
        {
            acb_mat_clear(mat2 + p);
            acb_mat_clear(mat1 + p);
            acb_mat_clear(mat + p);
        }
        flint_free(mat);
        _acb_vec_clear(pts, npts * 2);
        if (lcroots != NULL)
            _acb_vec_clear(lcroots, lcdeg1);
        acb_ode_exponents_clear(expos);
        gr_ore_poly_clear(dop1, Dop);
        gr_ore_poly_clear(dop, Dop);
        GR_TMP_CLEAR(pert, Pol);
        gr_ctx_clear(Dop);
        gr_ctx_clear(Pol);
        gr_ctx_clear(CC);
        acb_clear(r);
        mag_clear(refmag);
        mag_clear(mag);
    }

    gr_ctx_clear(CaCC);
    gr_ctx_clear(QQbar);
    gr_ctx_clear(ZZi);
    gr_ctx_clear(QQ);
    gr_ctx_clear(ZZ);

    TEST_GR_FUNCTION_END(state, count_success, count_domain, count_unable);
}
