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
    for (iter = 0; iter < 20 * flint_test_multiplier(); iter++)
    {
        slong g = 2;
        slong n = 1 << (2 * g);
        slong nb = acb_theta_jet_nb(1, g + 1);
        slong prec = 100 + n_randint(state, 100);
        slong bits = n_randint(state, 4);
        /* char str[] = "1620*g^2*a^2 + (-540*g*f*b + ((-504*g*e + 300*f^2)*c + (324*g*d^2 - 180*f*e*d + 48*e^3)))*a + ((300*g*e - 80*f^2)*b^2 + ((-180*g*d + 4*f*e)*c + (36*f*d^2 - 12*e^2*d))*b + (48*g*c^3 + (-12*f*d + 4*e^2)*c^2))"; */
        /* char vars[][1] = {"a", "b", "c", "d", "e", "f", "g"};*/
        char str[] = "1620*a6^2*a0^2 + ((300*a5^2 - 504*a4*a6)*a2 + ((-540*a1*a6 - 180*a4*a3)*a5 + (324*a3^2*a6 + 48*a4^3)))*a0 + (48*a6*a2^3 + (-12*a3*a5 + 4*a4^2)*a2^2 + (4*a4*a1*a5 - 180*a3*a1*a6)*a2 + (-80*a1^2*a5^2 + 36*a3^2*a1*a5 + (300*a4*a1^2*a6 - 12*a4^2*a3*a1)))";
        /* char vars[][2] = {"a0", "a1", "a2", "a3", "a4", "a5", "a6"};*/
        char** vars;
        fmpz_mpoly_ctx_t ctx;
        fmpz_mpoly_t pol;
        acb_mat_t tau;
        acb_ptr z, th2, dth, val;
        acb_t psi4, test;
        slong k;

        fmpz_mpoly_ctx_init(ctx, 7, ORD_LEX);
        vars = flint_malloc(7 * sizeof(char*));
        for (k = 0; k < 7; k++)
        {
            vars[k] = flint_malloc(2 * sizeof(char*));
            flint_sprintf(vars[k], "a%wd", k);
        }
        fmpz_mpoly_init(pol, ctx);
        acb_mat_init(tau, g, g);
        z = _acb_vec_init(g);
        th2 = _acb_vec_init(n);
        dth = _acb_vec_init(n * nb);
        val = _acb_vec_init(7);
        acb_init(psi4);
        acb_init(test);

        acb_siegel_randtest_reduced(tau, state, prec, bits);

        acb_theta_jet_all(dth, z, tau, 1, prec);
        acb_theta_all(th2, z, tau, 1, prec);
        acb_theta_g2_psi4(psi4, th2, prec);

        acb_mat_printd(tau, 5);
        flint_printf("dth, th2:\n");
        _acb_vec_printd(dth, n * nb, 5);
        _acb_vec_printd(th2, n, 5);

        acb_theta_g2_chi6m2(val, dth, prec);
        fmpz_mpoly_set_str_pretty(pol, str, (const char**) vars, ctx);
        acb_eval_fmpz_mpoly(test, pol, val, ctx, prec);
        acb_mul_2exp_si(test, test, -2);

            flint_printf("chi6m2:\n");
            _acb_vec_printd(val, 7, 5);
            flint_printf("pol:\n");
            fmpz_mpoly_print_pretty(pol, (const char**) vars, ctx);
            flint_printf("\npsi4, test:\n");
            acb_printd(psi4, 10);
            flint_printf("\n");
            acb_printd(test, 10);
            flint_printf("\n");
            
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
        for (k = 0; k < 7; k++)
        {
            flint_free(vars[k]);
        }
        flint_free(vars);
        acb_mat_clear(tau);
        _acb_vec_clear(z, g);
        _acb_vec_clear(th2, n);
        _acb_vec_clear(dth, n * nb);
        _acb_vec_clear(val, 7);
        acb_clear(psi4);
        acb_clear(test);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return 0;
}
