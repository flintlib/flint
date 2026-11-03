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
#include "acb_types.h"
#include "acb.h"
#include "acb_poly.h"
#include "acb_ode.h"
#include "gr.h"
#include "gr_vec.h"
#include "gr_poly.h"
#include "gr_ore_poly.h"

TEST_GR_FUNCTION_START(acb_ode_solution_series, state, count_success, count_domain, count_unable)
{
    gr_ctx_t ZZ;
    gr_ctx_init_fmpz(ZZ);

    for (slong iter = 0; iter < 200 * flint_test_multiplier(); iter++)
    {
        ulong seed1, seed2;
        flint_rand_get_seed(&seed1, &seed2, state);

        int status = GR_SUCCESS;
        int expect_success = 0, ordinary = 0;
        const char * fail = NULL;
        slong n = 0, n2 = 0;
        ulong flags = 0, lflags = 0;

        slong prec = 2 + n_randint(state, 200);

        gr_ctx_t CC, PolCC, Pol, PolE;
        gr_ore_poly_ctx_t Dop, Euler;
        gr_ctx_init_complex_acb(CC, prec);
        gr_ctx_init_gr_poly(PolCC, CC);
        if (n_randint(state, 2))
            gr_ctx_init_gr_poly(Pol, ZZ);
        else
            gr_ctx_init_gr_poly(Pol, CC);
        ore_algebra_t alg = n_randint(state, 4) ? ORE_ALGEBRA_EULER_DERIVATIVE
                                                : ORE_ALGEBRA_DERIVATIVE;
        gr_ore_poly_ctx_init(Dop, Pol, 0, alg);
        acb_ode_local_diffop_ctx_init(PolE, Euler, Dop, CC);

        gr_ore_poly_t dop, ldop;
        acb_ode_exponents_t expos;
        gr_vec_t ini;
        gr_poly_vec_t series, series2, res;
        acb_poly_t t;
        gr_ore_poly_init(dop, Dop);
        gr_ore_poly_init(ldop, Euler);
        acb_ode_exponents_init(expos);
        gr_vec_init(ini, 0, CC);
        gr_poly_vec_init(series, 0, CC);
        gr_poly_vec_init(series2, 0, CC);
        gr_poly_vec_init(res, 0, CC);
        acb_poly_init(t);

        /* operator */

        slong dop_len = n_randint(state, 8) ? 2 + n_randint(state, 4)
                                            : n_randint(state, 2);
        gr_ctx_struct * Cst = POLYNOMIAL_ELEM_CTX(Pol);

        if (alg == ORE_ALGEBRA_EULER_DERIVATIVE
            && Cst->which_ring == GR_CTX_CC_ACB && n_randint(state, 8))
        {
            slong clen = 1 + n_randint(state, 10);
            slong lcdeg = n_randint(state, clen);
            acb_ptr lcroots = _acb_vec_init(lcdeg);
            gr_ore_poly_fit_length(dop, dop_len, Dop);
            acb_ode_randtest_acb(dop->coeffs, lcroots, expos, state,
                                 n_randint(state, 4) ? 10 : 100,
                                 lcdeg, clen, dop_len, prec * 4);
            _gr_ore_poly_set_length(dop, dop_len, Dop);
            _acb_vec_clear(lcroots, lcdeg);
            expect_success = 1;
        }
        else
        {
            status |= gr_ore_poly_randtest(dop, state, dop_len, Dop);
            if (dop_len > 0 && n_randint(state, 8))
            {
                gr_ore_poly_fit_length(dop, dop_len, Dop);
                gr_ptr lc = gr_ore_poly_coeff_ptr(dop, dop_len - 1, Dop);
                status |= gr_randtest_not_zero(lc, state, Pol);
                status |= gr_randtest_not_zero(gr_poly_coeff_ptr(lc, 0, Cst),
                                               state, Cst);
                _gr_ore_poly_set_length(dop, dop_len, Dop);
                ordinary = (alg == ORE_ALGEBRA_DERIVATIVE);
                expect_success = (Cst->which_ring == GR_CTX_FMPZ);
            }
            status |= acb_ode_exponents(expos, dop, Dop, CC);
        }

        status |= acb_ode_local_diffop(ldop, &lflags, dop, Dop, Euler);

        if (status != GR_SUCCESS || expos->ngroups == 0)
        {
            status = GR_SUCCESS;  /* don't count test setup failures */
            goto end;
        }

        /* parameters */

        acb_ode_group_struct * group =
            expos->groups + n_randint(state, expos->ngroups);
        slong glen = acb_ode_group_length(group);

        gr_vec_set_length(ini, glen, CC);
        if (n_randint(state, 4))
            status |= _gr_vec_randtest(ini->entries, state, glen, CC);
        else if (glen > 0)
            acb_one((acb_ptr) ini->entries + n_randint(state, glen));

        n = n_randint(state, 4) ? n_randint(state, 20) : n_randint(state, 60);
        slong nlogs = acb_ode_group_nlogs(group, n - 1);
        flags = (ordinary && n_randint(state, 2)) ? ACB_ODE_ORDINARY : 0;

        /* test */

        /* wide enclosures of the roots of lc can make lc(0) contain zero */
        if (expect_success
            && gr_is_zero(gr_poly_coeff_ptr(
                              gr_ore_poly_coeff_ptr(dop, dop->length - 1, Dop),
                              0, Cst),
                          Cst) != T_FALSE)
            expect_success = 0;

        status |= acb_ode_solution_series(series, dop, Dop, group, ini, n,
                                          flags, CC);

        if (expect_success && status != GR_SUCCESS)
        {
            fail = "unexpected failure";
            goto end;
        }

        if (status != GR_SUCCESS)
            goto end;

        if (series->length > nlogs)
        {
            fail = "too many logs";
            goto end;
        }

        gr_poly_vec_set_length(series, nlogs, CC);
        acb_poly_struct * sers = (void *) series->entries;

        for (slong k = 0; k < nlogs; k++)
        {
            if (sers[k].length > n)
            {
                fail = "series too long";
                goto end;
            }
        }

        for (slong s = 0, j = 0; s < group->nshifts && group->shifts[s].n < n; s++)
        {
            slong mult = group->shifts[s].mult;
            j += mult;
            for (slong k = 0; k < mult; k++)
            {
                acb_srcptr c = acb_poly_get_coeff_ptr(sers + k, group->shifts[s].n);
                acb_srcptr a = (acb_srcptr) ini->entries + j - 1 - k;
                if (c == NULL ? !acb_contains_zero(a) : !acb_overlaps(c, a))
                {
                    fail = "initial values";
                    goto end;
                }
            }
        }

        gr_poly_vec_set_length(res, nlogs, CC);
        acb_poly_struct * r = (void *) res->entries;
        acb_ode_apply_diffop_polmul(r, ldop->coeffs, ldop->length,
                                    group->leader, 0, sers, nlogs, 0, n, prec);
        for (slong k = 0; k < nlogs; k++)
            for (slong i = 0; i < r[k].length; i++)
                if (!acb_contains_zero(r[k].coeffs + i))
                    fail = "residual";
        if (fail != NULL)
            goto end;

        /* shorter expansion, starting from a nonempty output vector */

        n2 = n_randint(state, n + 1);
        slong nlogs2 = acb_ode_group_nlogs(group, n2 - 1);

        gr_poly_vec_set_length(series2, n_randint(state, 4), CC);
        status |= _gr_vec_randtest(series2->entries, state, series2->length,
                                   PolCC);
        status |= acb_ode_solution_series(series2, dop, Dop, group, ini, n2,
                (ordinary && n_randint(state, 2)) ? ACB_ODE_ORDINARY : 0, CC);

        if (status != GR_SUCCESS)
            goto end;

        if (series2->length > nlogs2)
        {
            fail = "too many logs (2)";
            goto end;
        }

        gr_poly_vec_set_length(series2, nlogs2, CC);
        for (slong k = 0; k < nlogs2; k++)
        {
            acb_poly_set_trunc(t, sers + k, n2);
            if (!acb_poly_overlaps(t, (acb_poly_struct *) series2->entries + k))
            {
                fail = "consistency";
                goto end;
            }
        }

        /* unsupported coefficient domain */

        if (n_randint(state, 4) == 0)
        {
            gr_ctx_t RR, PolRR;
            gr_vec_t ini3;
            gr_poly_vec_t series3;
            gr_ctx_init_real_arb(RR, prec);
            gr_ctx_init_gr_poly(PolRR, RR);
            gr_vec_init(ini3, glen, RR);
            gr_poly_vec_init(series3, 0, RR);

            if (acb_ode_solution_series(series3, dop, Dop, group, ini3, n,
                                        flags, RR) == GR_SUCCESS)
                fail = "unsupported CC";

            gr_poly_vec_clear(series3, RR);
            gr_vec_clear(ini3, RR);
            gr_ctx_clear(PolRR);
            gr_ctx_clear(RR);
        }

end:
        if (fail != NULL)
        {
            flint_printf("FAIL: %s\n\n", fail);
            flint_printf("flint_rand_set_seed(state, 0x%wx, 0x%wx);\n", seed1, seed2);
            flint_printf("dop = %{gr} in %{gr_ctx}\n", dop, Dop, Dop);
            flint_printf("expos = "); acb_ode_exponents_println(expos);
            flint_printf("ini = %{acb*}\n", ini->entries, ini->length);
            flint_printf("n = %wd, n2 = %wd, flags = %wu, prec = %wd\n",
                         n, n2, flags, prec);
            flint_printf("series = %{gr*}\n", series->entries, series->length, PolCC);
            flint_printf("series2 = %{gr*}\n", series2->entries, series2->length, PolCC);
            flint_abort();
        }

        count_success += (status == GR_SUCCESS);
        count_domain += ((status & GR_DOMAIN) != 0);
        count_unable += ((status & GR_UNABLE) != 0);

        acb_poly_clear(t);
        gr_poly_vec_clear(res, CC);
        gr_poly_vec_clear(series2, CC);
        gr_poly_vec_clear(series, CC);
        gr_vec_clear(ini, CC);
        acb_ode_exponents_clear(expos);
        gr_ore_poly_clear(ldop, Euler);
        gr_ore_poly_clear(dop, Dop);
        gr_ctx_clear(Euler);
        gr_ctx_clear(PolE);
        gr_ctx_clear(Dop);
        gr_ctx_clear(Pol);
        gr_ctx_clear(PolCC);
        gr_ctx_clear(CC);
    }

    gr_ctx_clear(ZZ);

    TEST_GR_FUNCTION_END(state, count_success, count_domain, count_unable);
}
