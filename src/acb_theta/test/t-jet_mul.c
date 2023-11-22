/*
    Copyright (C) 2023 Jean Kieffer

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz_mpoly.h"
#include "acb.h"
#include "acb_theta.h"

TEST_FUNCTION_START(acb_theta_jet_mul, state)
{
    slong iter;

    /* Test: matches multiplication of fmpz_mpoly_t */
    for (iter = 0; iter < 100 * flint_test_multiplier(); iter++)
    {
        slong g = 1 + n_randint(state, 4);
        slong ord = n_randint(state, 10);
        slong nb = acb_theta_jet_nb(ord, g);
        slong prec = 100;
        slong * tups;
        fmpz_mpoly_ctx_t ctx;
        fmpz_mpoly_t p1, p2, p3;
        fmpz_t c;
        acb_ptr v1, v2, v3;
        acb_t x;
        slong k, t;

        tups = flint_malloc(g * nb * sizeof(slong));
        fmpz_mpoly_ctx_init(ctx, g, ORD_LEX);
        fmpz_mpoly_init(p1, ctx);
        fmpz_mpoly_init(p2, ctx);
        fmpz_mpoly_init(p3, ctx);
        fmpz_init(c);
        v1 = _acb_vec_init(nb);
        v2 = _acb_vec_init(nb);
        v3 = _acb_vec_init(nb);
        acb_init(x);

        acb_theta_jet_tuples(tups, ord, g);
        for (k = 0; k < nb; k++)
        {
            t = n_randint(state, 100);
            acb_set_si(&v1[k], t);
            fmpz_mpoly_set_coeff_si_ui(p1, t, (ulong *) tups + k * g, ctx);

            t = n_randint(state, 100);
            acb_set_si(&v2[k], t);
            fmpz_mpoly_set_coeff_si_ui(p2, t, (ulong *) tups + k * g, ctx);
        }

        acb_theta_jet_mul(v3, v1, v2, ord, g, prec);
        fmpz_mpoly_mul(p3, p1, p2, ctx);

        for (k = 0; k < nb; k++)
        {
            fmpz_mpoly_get_coeff_fmpz_ui(c, p3, (ulong *) tups + k * g, ctx);
            acb_set_fmpz(x, c);
            if (!acb_eq(x, &v3[k]))
            {
                flint_printf("FAIL\n");
                flint_printf("g = %wd, ord = %wd, k = %wd, vectors:\n", g, ord, k);
                _acb_vec_printd(v1, nb, 5);
                _acb_vec_printd(v2, nb, 5);
                _acb_vec_printd(v3, nb, 5);
                flint_abort();
            }
        }

        flint_free(tups);
        fmpz_mpoly_clear(p1, ctx);
        fmpz_mpoly_clear(p2, ctx);
        fmpz_mpoly_clear(p3, ctx);
        fmpz_mpoly_ctx_clear(ctx);
        fmpz_clear(c);
        _acb_vec_clear(v1, nb);
        _acb_vec_clear(v2, nb);
        _acb_vec_clear(v3, nb);
        acb_clear(x);
    }

    TEST_FUNCTION_END(state);
}
