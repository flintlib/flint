/*
    Copyright (C) 2023 Jean Kieffer

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "acb_theta.h"

int main(void)
{
    slong iter;
    flint_rand_t state;

    flint_printf("g2_psi4....");
    fflush(stdout);

    flint_randinit(state);

    /* Test: agrees with polynomial expression in terms of chi6m2 */
    for (iter = 0; iter < 10 * flint_test_multiplier(); iter++)
    {
        slong g = 2;
        slong n = 1 << (2 * g);
        slong nb = acb_theta_jet_nb(1, g + 1);
        slong prec = 100 + n_randint(state, 100);
        slong bits = n_randint(state, 4);
        char str[] = "1620*a6^2*a0^2 + ((300*a5^2 - 504*a4*a6)*a2 + ((-540*a1*a6 - 180*a4*a3)*a5 + (324*a3^2*a6 + 48*a4^3)))*a0 + (48*a6*a2^3 + (-12*a3*a5 + 4*a4^2)*a2^2 + (4*a4*a1*a5 - 180*a3*a1*a6)*a2 + (-80*a1^2*a5^2 + 36*a3^2*a1*a5 + (300*a4*a1^2*a6 - 12*a4^2*a3*a1)))";
        char* vars[] = {"a0", "a1", "a2", "a3", "a4", "a5", "a6"};
        fmpz_mpoly_ctx_t ctx;
        fmpz_mpoly_t pol;
        acb_mat_t tau;
        acb_ptr z, th2, dth, val;
        acb_poly_t chi;
        acb_t psi4, test;
        slong k;

        fmpz_mpoly_ctx_init(ctx, 7, ORD_LEX);
        fmpz_mpoly_init(pol, ctx);
        acb_mat_init(tau, g, g);
        z = _acb_vec_init(g);
        th2 = _acb_vec_init(n);
        dth = _acb_vec_init(n * nb);
        val = _acb_vec_init(7);
        acb_poly_init(chi);
        acb_init(psi4);
        acb_init(test);

        acb_siegel_randtest_reduced(tau, state, prec, bits);

        acb_theta_jet_all(dth, z, tau, 1, prec);
        acb_theta_all(th2, z, tau, 1, prec);
        acb_theta_g2_psi4(psi4, th2, prec);

        acb_theta_g2_chi6m2(chi, dth, prec);
        for (k = 0; k <= 6; k++)
        {
            acb_poly_get_coeff_acb(&val[k], chi, 6 - k);
        }
        fmpz_mpoly_set_str_pretty(pol, str, (const char**) vars, ctx);
        acb_eval_fmpz_mpoly(test, pol, val, ctx, prec);
        acb_mul_2exp_si(test, test, -2);

        if (!acb_overlaps(psi4, test))
        {
            flint_printf("FAIL\n");
            flint_printf("chi6m2:\n");
            _acb_vec_printd(val, 7, 5);
            flint_printf("pol:\n");
            fmpz_mpoly_print_pretty(pol, (const char**) vars, ctx);
            flint_printf("\npsi4, test:\n");
            acb_printd(psi4, 10);
            flint_printf("\n");
            acb_printd(test, 10);
            flint_printf("\n");
            flint_abort();
        }

        fmpz_mpoly_clear(pol, ctx);
        fmpz_mpoly_ctx_clear(ctx);
        acb_mat_clear(tau);
        _acb_vec_clear(z, g);
        _acb_vec_clear(th2, n);
        _acb_vec_clear(dth, n * nb);
        _acb_vec_clear(val, 7);
        acb_poly_clear(chi);
        acb_clear(psi4);
        acb_clear(test);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return 0;
}
