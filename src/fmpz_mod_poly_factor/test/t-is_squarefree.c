/*
    Copyright (C) 2011 Fredrik Johansson
    Copyright (C) 2012 Lina Kulakova

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "ulong_extras.h"
#include "fmpz.h"
#include "fmpz_mod.h"
#include "fmpz_mod_poly.h"
#include "fmpz_mod_poly_factor.h"

TEST_FUNCTION_START(fmpz_mod_poly_factor_is_squarefree, state)
{
    int iter;
    fmpz_mod_ctx_t ctx;

    fmpz_mod_ctx_init_ui(ctx, 2);

    for (iter = 0; iter < 200 * flint_test_multiplier(); iter++)
    {
        fmpz_mod_poly_t poly, Q, R, t;
        fmpz_t modulus;
        mp_limb_t mod;
        slong i, num_factors, exp, max_exp;
        int v, result;

        mod = n_randtest_prime(state, 0);
        fmpz_init_set_ui(modulus, mod);
        fmpz_mod_ctx_set_modulus(ctx, modulus);

        fmpz_mod_poly_init(poly, ctx);
        fmpz_mod_poly_init(t, ctx);
        fmpz_mod_poly_init(Q, ctx);
        fmpz_mod_poly_init(R, ctx);

        fmpz_mod_poly_set_coeff_ui(poly, 0, n_randint(state, mod), ctx);
        num_factors = n_randint(state, 5);

        max_exp = 0;
        for (i = 0; i < num_factors; i++)
        {
            do {
                fmpz_mod_poly_randtest(t, state, n_randint(state, 10), ctx);
            } while (!fmpz_mod_poly_is_irreducible(t, ctx) ||
                    (fmpz_mod_poly_length(t, ctx) < 2));

            exp = n_randint(state, 4) + 1;
            if (n_randint(state, 2) == 0)
                exp = 1;

            fmpz_mod_poly_divrem(Q, R, poly, t, ctx);
            if (!fmpz_mod_poly_is_zero(R, ctx))
            {
                fmpz_mod_poly_pow(t, t, exp, ctx);
                fmpz_mod_poly_mul(poly, poly, t, ctx);
                max_exp = FLINT_MAX(exp, max_exp);
            }
        }

        v = fmpz_mod_poly_is_squarefree(poly, ctx);

        if (v == 1)
            result = (max_exp <= 1 && !fmpz_mod_poly_is_zero(poly, ctx));
        else
            result = (max_exp > 1 || fmpz_mod_poly_is_zero(poly, ctx));

        if (!result)
        {
            flint_printf("FAIL: ");
            fmpz_print(modulus);
            flint_printf(" %wd, %d\n", max_exp, v);
            fmpz_mod_poly_print(poly, ctx); flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_clear(modulus);
        fmpz_mod_poly_clear(poly, ctx);
        fmpz_mod_poly_clear(t, ctx);
        fmpz_mod_poly_clear(Q, ctx);
        fmpz_mod_poly_clear(R, ctx);
    }

    fmpz_mod_ctx_clear(ctx);

    TEST_FUNCTION_END(state);
}
