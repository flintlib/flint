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

    flint_printf("g2_covariant....");
    fflush(stdout);

    flint_randinit(state);

    /* Test: agrees with g2_psi4 using psi4 = -(Co20 - 3*Co40)/20 */
    for (iter = 0; iter < 5 * flint_test_multiplier(); iter++)
    {
        slong prec = 100 + n_randint(state, 100);
        slong g = 2;
        slong n = 1 << (2 * g);
        acb_mat_t tau;
        acb_ptr z, th2;
        acb_poly_struct* r;
        acb_poly_t u;
        fmpz_mpoly_ctx_t ctx;
        fmpz_mpoly_t pol, v;
        acb_t psi4, test;
        slong k;

        acb_mat_init(tau, g, g);
        z = _acb_vec_init(g);
        th2 = _acb_vec_init(n);
        r = flint_malloc(26 * sizeof(acb_poly_struct));
        for (k = 0; k < 26; k++)
        {
            acb_poly_init(&r[k]);
        }
        acb_poly_init(u);
        fmpz_mpoly_ctx_init(ctx, 26, ORD_LEX);
        fmpz_mpoly_init(pol, ctx);
        fmpz_mpoly_init(v, ctx);
        acb_init(psi4);
        acb_init(test);

        acb_siegel_randtest_nice(tau, state, prec);
        acb_theta_all(th2, z, tau, 1, prec);
        acb_theta_g2_psi4(psi4, th2, prec);
        acb_mul_si(psi4, psi4, -20, prec);

        acb_theta_g2_basic_covariants(r, tau, prec);
        fmpz_mpoly_gen(pol, 1, ctx);
        fmpz_mpoly_mul(pol, pol, pol, ctx);
        fmpz_mpoly_gen(v, 8, ctx);
        fmpz_mpoly_scalar_mul_si(v, v, -3, ctx);
        fmpz_mpoly_add(pol, pol, v, ctx);

        acb_theta_g2_covariant(u, pol, r, ctx, prec);
        acb_poly_get_coeff_acb(test, u, 0);

        if (!acb_overlaps(psi4, test))
        {
            flint_printf("FAIL\n");
            acb_mat_printd(tau, 5);
            flint_printf("psi4, test:\n");
            acb_printd(psi4, 10);
            flint_printf("\n");
            acb_printd(test, 10);
            flint_printf("\nu:\n");
            acb_poly_printd(u, 5);
            flint_printf("\ncovariants:\n");
            for (k = 0; k < 26; k++)
            {
                acb_poly_printd(&r[k], 5);
                flint_printf("\n");
            }
            flint_abort();
        }

        acb_mat_clear(tau);
        _acb_vec_clear(z, g);
        _acb_vec_clear(th2, n);
        for (k = 0; k < 26; k++)
        {
            acb_poly_clear(&r[k]);
        }
        flint_free(r);
        acb_poly_clear(u);
        fmpz_mpoly_clear(pol, ctx);
        fmpz_mpoly_clear(v, ctx);
        fmpz_mpoly_ctx_clear(ctx);
        acb_clear(psi4);
        acb_clear(test);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return 0;
}
