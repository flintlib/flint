/*
    Copyright (C) 2023 Jean Kieffer

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "acb_theta.h"

int
main()
{
    slong iter;
    flint_rand_t state;

    flint_printf("agm_ctx_set....");
    fflush(stdout);

    flint_randinit(state);

    /* Test: agm context must be valid for g=1, 2 and tau in fundamental domain */
    for (iter = 0; iter < 5 * arb_test_multiplier(); iter++)
    {
        slong g = 1 + n_randint(state, 2);
        slong prec = ACB_THETA_AGM_BASEPREC + n_randint(state, 1000);
        acb_mat_t tau;
        acb_ptr z;
        arf_t rad;
        acb_theta_agm_ctx_t ctx;
        int res;
        slong k;

        acb_mat_init(tau, g, g);
        z = _acb_vec_init(g);
        arf_init(rad);

        acb_siegel_randtest_fund(tau, state, prec);

        _acb_vec_zero(z, g);
        arf_one(rad);
        arf_mul_2exp_si(rad, rad, -4);
        for (k = 0; k < g; k++)
            acb_randtest_disk(&z[k], &z[k], rad, state, prec);

        if (iter % 2 == 0)
            acb_theta_agm_ctx_init(ctx, tau);
        else
            acb_theta_agm_ctx_init_ext(ctx, z, tau);

        res = acb_theta_agm_ctx_set(ctx, prec);

        if (!res)
        {
            flint_printf("FAIL\n");
            fflush(stdout);
            flint_abort();
        }

        acb_mat_clear(tau);
        _acb_vec_clear(z, g);
        arf_clear(rad);
        acb_theta_agm_ctx_clear(ctx);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}
