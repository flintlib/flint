/*
    Copyright (C) 2024 Jean Kieffer

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "acb.h"
#include "arb_mat.h"
#include "acb_theta.h"

TEST_FUNCTION_START(acb_theta_ctx_z_set, state)
{
    slong iter;

    /* Test: c corresponds to term for n = -a (when n fits in slongs) */
    for (iter = 0; iter < 200 * flint_test_multiplier(); iter++)
    {
        slong g = 1 + n_randint(state, 4);
        slong prec = 100 + n_randint(state, 200);
        slong mag_bits = n_randint(state, 10);
        acb_mat_t tau;
        acb_ptr z;
        acb_theta_ctx_tau_t ctx_tau;
        acb_theta_ctx_z_t ctx;
        acb_t t;
        slong * n = flint_malloc(g * sizeof(slong));
        fmpz_t m;
        slong k;
        int small = 1;

        acb_mat_init(tau, g, g);
        z = _acb_vec_init(g);
        acb_theta_ctx_tau_init(ctx_tau, g);
        acb_theta_ctx_z_init(ctx, g);
        acb_init(t);
        fmpz_init(m);

        acb_siegel_randtest(tau, state, prec, mag_bits);
        acb_siegel_randtest_vec(z, state, g, prec);
        acb_theta_ctx_tau_set(ctx_tau, tau, prec);
        acb_theta_ctx_z_set(ctx, z, ctx_tau, prec);

        for (k = 0; k < g; k++)
        {
            arb_get_unique_fmpz(m, &acb_theta_ctx_r(ctx)[k]);
            small = small && fmpz_fits_si(m);
            n[k] = (small ? -fmpz_get_si(m) : 0);
        }
        acb_theta_naive_term(t, z, tau, NULL, n, prec);

        if (small && !acb_overlaps(t, acb_theta_ctx_c(ctx)))
        {
            flint_printf("FAIL\n");
            acb_mat_printd(tau, 5);
            _acb_vec_printd(z, g, 5);
            flint_printf("\nn =");
            for (k = 0; k < g; k++)
            {
                flint_printf(" %wd", n[k]);
            }
            flint_printf("\nu = ");
            arb_printd(acb_theta_ctx_u(ctx), 5);
            if (g > 1)
            {
                flint_printf("\nv = ");
                _arb_vec_printd(acb_theta_ctx_v(ctx), g, 5);
            }
            flint_printf("t = ");
            acb_printd(t, 5);
            flint_printf("\nc = ");
            acb_printd(acb_theta_ctx_c(ctx), 5);
            flint_printf("\n");
            flint_abort();
        }

        flint_free(n);
        acb_mat_clear(tau);
        _acb_vec_clear(z, g);
        acb_theta_ctx_tau_clear(ctx_tau);
        acb_theta_ctx_z_clear(ctx);
        acb_clear(t);
        fmpz_clear(m);
    }

    TEST_FUNCTION_END(state);
}
