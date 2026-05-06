#include "test_helpers.h"
#include "acb_types.h"
#include "acb_mat.h"
#include "acb_ode.h"
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

        gr_ctx_t CC, Pol;
        gr_ore_poly_ctx_t Dop;
        gr_ptr pert;
        gr_ore_poly_t dop, dop1;
        acb_t r;
        arf_t frac;

        acb_init(r);
        arf_init(frac);

        slong prec = 2 + n_randint(state, 100);
        slong prec1 = 2 + n_randint(state, prec);
        gr_ctx_init_complex_acb(CC, prec);

        int expect_success = 0;
        if (!n_randint(state, 8))
            gr_ore_poly_ctx_init_randtest2(Pol, Dop, state);
        else
        {
            switch (n_randint(state, 8))
            {
                case 0:
                    gr_ctx_init_random_poly(Pol, state);
                    break;
                case 1:
                    gr_ctx_init_gr_poly(Pol, CC);
                    break;
                case 2:
                    gr_ctx_init_gr_poly(Pol, QQ);
                    break;
                case 3:
                    gr_ctx_init_gr_poly(Pol, ZZi);
                    break;
                case 4:
                    gr_ctx_init_gr_poly(Pol, QQbar);
                    break;
                case 5:
                    gr_ctx_init_gr_poly(Pol, CaCC);
                    break;
                default:
                    gr_ctx_init_gr_poly(Pol, ZZ);
                    expect_success = 1;
                    break;
            }
            ore_algebra_t alg = n_randint(state, 4)
                ? ORE_ALGEBRA_EULER_DERIVATIVE
                : ORE_ALGEBRA_DERIVATIVE;
            if (alg != ORE_ALGEBRA_EULER_DERIVATIVE)  // XXX
                expect_success = 0;
            gr_ore_poly_ctx_init(Dop, Pol, 0, alg);
        }

        GR_TMP_INIT(pert, Pol);
        gr_ore_poly_init(dop,  Dop);
        gr_ore_poly_init(dop1, Dop);

        /* todo: Test expos != NULL and lcrt != NULL (by generating the leading
           coefficient and/or indicial polynomial from its roots, cf.
           t-shiftless_decomposition). Update expect_success as needed.

           todo: With reasonable probability, generate evaluation points in the
           disk of convergence. */

        /* also generate operators with coefficients of higher degree? */
        slong dop_len_bound = 20;
        if (Pol->which_ring == GR_CTX_GR_POLY
            && (POLYNOMIAL_ELEM_CTX(Pol)->methods == _ca_methods
                || POLYNOMIAL_ELEM_CTX(Pol)->methods ==  _qqbar_methods))
            dop_len_bound = 4;
        status |= gr_ore_poly_randtest(dop, state,
                                       n_randint(state, dop_len_bound), Dop);
        status |= gr_randtest_not_zero(pert, state, Pol);
        /* dop1 may have extra singular points; ensure however that there is no
           extra singular point at zero in interesting cases */
        if (Pol->which_ring == GR_CTX_GR_POLY)
            status |= gr_randtest_not_zero(
                    gr_poly_coeff_ptr(pert, 0, POLYNOMIAL_ELEM_CTX(Pol)),
                    state,
                    POLYNOMIAL_ELEM_CTX(Pol));
        status |= gr_ore_poly_other_mul(dop1, pert, Pol, dop, Dop);
        slong order = dop->length - 1;

        // XXX large npts hits hard entire function case too often
        // slong npts = n_randint(state, n_randint(state, 8) ? 3 : 100);
        slong npts = n_randint(state, 3);
        acb_ptr pts = _acb_vec_init(npts * 2);
        acb_ptr pts1 = pts + npts;
        for (slong p = 0; p < npts; p++)
        {
            acb_randtest_precise(pts + p, state, prec + 30, 3);
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
        flint_printf("prec = %wd, prec1 = %wd\n", prec, prec1);

        status |= acb_ode_fundamental_matrix_vec(mat, dop, Dop, NULL, NULL,
                                                 pts, npts, basis, prec + 30);
        flint_printf("status0 = %wd\n", status);

        if (expect_success && status == GR_UNABLE)
        {
            flint_printf("FAIL: unable\n");
            flint_abort();
        }

        if (status != GR_SUCCESS)
        {
            goto cleanup;
        }

        int test0_is_accurate = 1;
        for (slong p = 0; p < npts; p++)
            for (slong i = 0; i < acb_mat_nrows(mat + p); i++)
                for (slong j = 0; j < order; j++)
                {
                    slong acc = acb_rel_accuracy_bits(
                            acb_mat_entry_ptr(mat + p, i, j));
                    if (acc < prec)
                    {
                        flint_printf("poor accuracy @pt%wd (%wd,%wd) = %wd\n",
                                     p, i, j, acc);
                        test0_is_accurate = 0;
                    }
                }

        for (slong p = 0; p < npts; p++)  // tmp
        {
            flint_printf("res@%{acb} =\n", pts + p);
            acb_mat_printd(mat + p, prec+4);
        }

        status |= acb_ode_fundamental_matrix_vec(mat1, dop, Dop, NULL, NULL,
                                                 pts1, npts, basis, prec1);
        flint_printf("status1 = %wd\n", status);

        if (status == GR_SUCCESS)  // tmp
            for (slong p = 0; p < npts; p++)
            {
                flint_printf("res1@%{acb} =\n", pts1 + p);
                acb_mat_printd(mat1 + p, prec+4);
            }

        if (status == GR_SUCCESS)
        {
            for (slong p = 0; p < npts; p++)
                if (test0_is_accurate)
                {
                    if (!acb_mat_contains(mat1 + p, mat + p))
                    {
                        flint_printf("FAIL: containment (@pt%ld)\n\n", p);
                        flint_abort();
                    }
                }
                else if (!acb_mat_overlaps(mat1 + p, mat + p))
                {
                    flint_printf("FAIL: overlap 1 (@pt%ld)\n\n", p);
                    flint_abort();
                }
        }

        if (Pol->which_ring == GR_CTX_GR_POLY
            && POLYNOMIAL_ELEM_CTX(Pol)->methods == _ca_methods)
            goto cleanup;

        status |= acb_ode_fundamental_matrix_vec(mat2, dop1, Dop, NULL, NULL,
                                                 pts, npts, basis, prec);
        flint_printf("status2 = %wd\n", status);

        if (status == GR_SUCCESS)
        {
            for (slong p = 0; p < npts; p++)
                if (!acb_mat_overlaps(mat2 + p, mat + p))
                {
                    flint_printf("FAIL: overlap 2\n\n");
                    flint_abort();
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
        gr_ore_poly_clear(dop1, Dop);
        gr_ore_poly_clear(dop, Dop);
        GR_TMP_CLEAR(pert, Pol);
        gr_ctx_clear(Dop);
        gr_ctx_clear(Pol);
        gr_ctx_clear(CC);
        arf_clear(frac);
        acb_clear(r);
    }

    gr_ctx_clear(CaCC);
    gr_ctx_clear(QQbar);
    gr_ctx_clear(ZZi);
    gr_ctx_clear(QQ);
    gr_ctx_clear(ZZ);

    TEST_FUNCTION_END(state);
}

