#include "test_helpers.h"
#include "acb_types.h"
#include "acb_mat.h"
#include "acb_ode.h"
#include "acb_poly.h"
#include "acb.h"
#include "arf.h"
#include "gr_poly.h"

FLINT_DLL extern gr_static_method_table _ca_methods,  _qqbar_methods;

TEST_FUNCTION_START(acb_ode_fundamental_matrix, state)
{
    gr_ctx_t ZZ, QQ, ZZi, QQbar, CaCC;

    gr_ctx_init_fmpz(ZZ);
    gr_ctx_init_fmpq(QQ);
    gr_ctx_init_fmpzi(ZZi);
    gr_ctx_init_complex_qqbar(QQbar);
    gr_ctx_init_complex_ca(CaCC);

    for (slong iter = 0; iter < 1000 * flint_test_multiplier(); iter++)
    {
        int status = GR_SUCCESS;
        int expect_success = 0;

        gr_ctx_t CC, Pol;
        gr_ore_poly_ctx_t Dop;
        gr_ptr pert;
        gr_ore_poly_t dop, dop1;
        acb_t r;
        arf_t frac;
        acb_ode_exponents_t expos;
        acb_ode_exponents_struct * expos1;
        acb_ptr lcroots = NULL;
        slong lcdeg1;
        mag_t refmag, mag;

        acb_init(r);
        arf_init(frac);
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
            slong rnd = n_randint(state, 8);
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
            {
                gr_ctx_init_gr_poly(Pol, CC);
            }

            ore_algebra_t alg = ORE_ALGEBRA_EULER_DERIVATIVE;
            /* todo: support standard derivative */
            /* ore_algebra_t alg = n_randint(state, 4)
                ? ORE_ALGEBRA_EULER_DERIVATIVE
                : ORE_ALGEBRA_DERIVATIVE; */
            gr_ore_poly_ctx_init(Dop, Pol, 0, alg);
        }

        /* operator */

        GR_TMP_INIT(pert, Pol);
        gr_ore_poly_init(dop,  Dop);
        gr_ore_poly_init(dop1, Dop);

        slong dop_len_bound = 20;
        if (Pol->which_ring == GR_CTX_GR_POLY
            && (POLYNOMIAL_ELEM_CTX(Pol)->methods == _ca_methods
                || POLYNOMIAL_ELEM_CTX(Pol)->methods ==  _qqbar_methods))
            dop_len_bound = 4;
        slong dop_len = n_randint(state, dop_len_bound);

        slong disp = n_randint(state, 16) ? 10 : 1000;

        if (GR_ORE_POLY_CTX(Dop)->which_algebra == ORE_ALGEBRA_EULER_DERIVATIVE
            && Pol->which_ring == GR_CTX_GR_POLY
            /* todo: acb_ode_randtest for exact coefficients */
            && POLYNOMIAL_ELEM_CTX(Pol)->which_ring == GR_CTX_CC_ACB
            && n_randint(state, 8))
        {
            /* todo: increase once sum_forward_divconquer is able to stop in the
               middle of a block */
            slong clen = 1 + (n_randint(state, 8) ? n_randint(state, 10)
                                                  : n_randint(state, 100));
            slong lcdeg = n_randint(state, clen);
            lcdeg1 = lcdeg + n_randint(state, 4);
            lcroots = _acb_vec_init(lcdeg1);

            gr_ore_poly_fit_length(dop, dop_len, Dop);
            acb_ode_randtest_acb(dop->coeffs, lcroots, expos, state, disp,
                                 lcdeg, clen, dop_len, prec);
            _gr_ore_poly_set_length(dop, dop_len, Dop);

            for (slong i = lcdeg; i < lcdeg1; i++)  /* extra lcroots for dop1 */
            {
                acb_randtest(lcroots + i, state, prec, 3);
                if (acb_is_zero(lcroots + i))
                    acb_one(lcroots + i);
            }
            acb_poly_product_roots(pert, lcroots + lcdeg, lcdeg1 - lcdeg, prec);

            expect_success = 1;
        }
        else
        {
            status |= gr_ore_poly_randtest(dop, state, dop_len, Dop);
            /* precompute the exponents and pass them to the main function so
               that exponent groups are ordered consistently across calls */
            status |= acb_ode_exponents(expos, dop, Dop, CC);

            for (slong g = 0; g < expos->ngroups; g++)
            {
                acb_ode_group_struct * grp = expos->groups + g;
                for (slong s = 0; s < grp->nshifts; s++)
                    if (grp->shifts[s].n > disp)
                        status = GR_UNABLE;
            }

            status |= gr_randtest_not_zero(pert, state, Pol);
            /* ensure however that dop1 has no extra singular point at zero in
               interesting cases */
            if (Pol->which_ring == GR_CTX_GR_POLY)
                status |= gr_randtest_not_zero(
                        gr_poly_coeff_ptr(pert, 0, POLYNOMIAL_ELEM_CTX(Pol)),
                        state,
                        POLYNOMIAL_ELEM_CTX(Pol));
        }

        if (expos->ngroups == 1 && n_randint(state, 4))
            expos1 = NULL;
        else
            expos1 = expos;

        status |= gr_ore_poly_other_mul(dop1, pert, Pol, dop, Dop);
        slong order = dop->length - 1;

        /* evaluation points */

        mag_inf(refmag);
        if (lcroots)
            _acb_vec_get_mag_lower(refmag, lcroots, lcdeg1);
        if (mag_is_inf(refmag))
            mag_set_ui(refmag, 3);

        // XXX large npts hits hard entire function case too often
        // slong npts = n_randint(state, n_randint(state, 8) ? 3 : 100);
        slong npts = n_randint(state, 3);
        acb_ptr pts = _acb_vec_init(npts * 2);
        acb_ptr pts1 = pts + npts;
        for (slong p = 0; p < npts; p++)
        {
            acb_randtest_precise(pts + p, state, prec + 30, 3);
            acb_get_mag(mag, pts + p);
            mag_div(mag, mag, refmag);
            slong s = FLINT_MAX(0., mag_get_d_log2_approx(mag) + 1.5);
            acb_mul_2exp_si(pts + p, pts + p, -s);

            if (n_randint(state, 2))
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
            slong nrows = n_randint(state, order + 2);
            acb_mat_init(mat + p,  nrows, order);
            acb_mat_init(mat1 + p, nrows, order);
            acb_mat_init(mat2 + p, nrows, order);
        }

        if (status != GR_SUCCESS)
            goto cleanup;

        flint_printf("******** testing: ********\n", Dop);
        flint_printf("Dop = %{gr_ctx}\n", Dop);
        flint_printf("dop = %{gr}\n", dop, Dop);
        flint_printf("dop1 = %{gr}\n", dop1, Dop);
        flint_printf("pts = %{acb*}\n", pts, npts);
        flint_printf("pts1 = %{acb*}\n", pts1, npts);
        flint_printf("expos = "); acb_ode_exponents_println(expos);
        flint_printf("expos1 = %p\n", expos1);
        flint_printf("prec = %wd, prec1 = %wd\n", prec, prec1);

        flint_printf("--- run 0 ---\n");

        status |= acb_ode_fundamental_matrix_vec(mat, dop, Dop, expos, lcroots,
                                                 pts, npts, basis, prec + 30);
        flint_printf("status0 = %wd\n", status);

        if (expect_success && status == GR_UNABLE)
        {
            flint_printf("FAIL: unexpected failure\n");
            flint_abort();
        }

        if (status != GR_SUCCESS)
            goto cleanup;

        // for (slong p = 0; p < npts; p++)
        //     for (slong i = 0; i < acb_mat_nrows(mat + p); i++)
        //         for (slong j = 0; j < order; j++)
        //         {
        //             slong acc = acb_rel_accuracy_bits(
        //                     acb_mat_entry_ptr(mat + p, i, j));
        //             if (acc < prec)
        //                 flint_printf("poor accuracy @pt%wd (%wd,%wd) = %wd\n",
        //                              p, i, j, acc);
        //         }

        // for (slong p = 0; p < npts; p++)
        // {
        //     flint_printf("res@%{acb} =\n", pts + p);
        //     acb_mat_printd(mat + p, prec+4);
        // }

        flint_printf("--- run 1 ---\n");

        status |= acb_ode_fundamental_matrix_vec(mat1, dop, Dop, expos, lcroots,
                                                 pts1, npts, basis, prec1);
        flint_printf("status1 = %wd\n", status);

        // if (status == GR_SUCCESS)
        //     for (slong p = 0; p < npts; p++)
        //     {
        //         flint_printf("res1@%{acb} =\n", pts1 + p);
        //         acb_mat_printd(mat1 + p, prec+4);
        //     }

        if (status == GR_SUCCESS)
        {
            for (slong p = 0; p < npts; p++)
            {
                slong failed = 0;
                for (slong i = 0; i < mat[p].r; i++)
                    for (slong j = 0; j < mat[p].c; j++)
                    {
                        acb_ptr c = acb_mat_entry_ptr(mat + p, i, j);
                        acb_ptr c1 = acb_mat_entry_ptr(mat1 + p, i, j);
                        if (!acb_overlaps(c, c1))
                        {
                            flint_printf("\nres[%ld,%ld] ", i, j);
                            acb_printn(c, prec/3, 0);
                            flint_printf(" ≠ ");
                            acb_printn(c1, prec/3, 0);
                            flint_printf("\n");
                            // acb_print(c);
                            // flint_printf(" ≠ ");
                            // acb_print(c1);
                            // flint_printf("\n");
                            failed = 1;
                        }
                    }
                if (failed)
                {
                    flint_printf("FAIL: overlap 1\n\n");
                    flint_printf("pts[%ld] = %{acb}\n", p, pts + p);
                    flint_printf("pts1[%ld] = %{acb}\n", p, pts1 + p);
                    flint_abort();
                }
            }
        }

        if (Pol->which_ring == GR_CTX_GR_POLY
            && POLYNOMIAL_ELEM_CTX(Pol)->methods == _ca_methods)
            goto cleanup;

        flint_printf("--- run 2 ---\n");

        status |= acb_ode_fundamental_matrix_vec(mat2, dop1, Dop, expos, lcroots,
                                                 pts, npts, basis, prec);
        flint_printf("status2 = %wd\n", status);

        if (status == GR_SUCCESS)
        {
            for (slong p = 0; p < npts; p++)
            {
                slong failed = 0;
                for (slong i = 0; i < mat[p].r; i++)
                    for (slong j = 0; j < mat[p].c; j++)
                    {
                        acb_ptr c = acb_mat_entry_ptr(mat + p, i, j);
                        acb_ptr c2 = acb_mat_entry_ptr(mat2 + p, i, j);
                        if (!acb_overlaps(c, c2))
                        {
                            flint_printf("\nres[%ld,%ld] ", i, j);
                            acb_printn(c, prec/3, 0);
                            flint_printf(" ≠ ");
                            acb_printn(c2, prec/3, 0);
                            flint_printf("\n");
                            // acb_print(c);
                            // flint_printf(" ≠ ");
                            // acb_print(c1);
                            // flint_printf("\n");
                            failed = 1;
                        }
                    }
                if (failed)
                {
                    flint_printf("FAIL: overlap 2\n\n");
                    flint_abort();
                }
            }
        }

cleanup:
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
        arf_clear(frac);
        acb_clear(r);
        mag_clear(refmag);
        mag_clear(mag);
    }

    gr_ctx_clear(CaCC);
    gr_ctx_clear(QQbar);
    gr_ctx_clear(ZZi);
    gr_ctx_clear(QQ);
    gr_ctx_clear(ZZ);

    TEST_FUNCTION_END(state);
}

