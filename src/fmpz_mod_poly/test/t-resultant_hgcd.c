/*
    Copyright (C) 2011 Sebastian Pancratz
    Copyright (C) 2014 William Hart

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

TEST_FUNCTION_START(fmpz_mod_poly_resultant_hgcd, state)
{
    int i, result;
    fmpz_mod_ctx_t ctx;

    fmpz_mod_ctx_init_ui(ctx, 2);

    /* Check res(f, g) == (-1)^(deg f deg g) res(g, f) */
    for (i = 0; i < 60 * flint_test_multiplier(); i++)
    {
        fmpz_t x, y, n;
        fmpz_mod_poly_t f, g;

        fmpz_init(n);
        fmpz_init(x);
        fmpz_init(y);

        fmpz_set_ui(n, n_randtest_prime(state, 0));
        fmpz_mod_ctx_set_modulus(ctx, n);

        fmpz_mod_poly_init(f, ctx);
        fmpz_mod_poly_init(g, ctx);

        fmpz_mod_poly_randtest(f, state, n_randint(state, 300), ctx);
        fmpz_mod_poly_randtest(g, state, n_randint(state, 300), ctx);

        fmpz_mod_poly_resultant_hgcd(x, f, g, ctx);
        fmpz_mod_poly_resultant_hgcd(y, g, f, ctx);

        if ((fmpz_mod_poly_degree(f, ctx) * fmpz_mod_poly_degree(g, ctx)) % 2)
           fmpz_mod_neg(y, y, ctx);

        result = (fmpz_equal(x, y));
        if (!result)
        {
            flint_printf("FAIL (res(f, g) == (-1)^(deg f deg g) res(g, f)):\n");
            fmpz_mod_poly_print(f, ctx), flint_printf("\n\n");
            fmpz_mod_poly_print(g, ctx), flint_printf("\n\n");
            printf("x = "); fmpz_print(x); printf("\n");
            printf("y = "); fmpz_print(y); printf("\n");
            printf("n = "); fmpz_print(n); printf("\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_mod_poly_clear(f, ctx);
        fmpz_mod_poly_clear(g, ctx);
        fmpz_clear(n);
        fmpz_clear(x);
        fmpz_clear(y);
    }

    /* Check res(f h, g) == res(f, g) res(h, g) */
    for (i = 0; i < 60 * flint_test_multiplier(); i++)
    {
        fmpz_t x, y, z, n;
        fmpz_mod_poly_t f, g, h;

        fmpz_init(n);
        fmpz_init(x);
        fmpz_init(y);
        fmpz_init(z);

        fmpz_set_ui(n, n_randtest_prime(state, 0));
        fmpz_mod_ctx_set_modulus(ctx, n);

        fmpz_mod_poly_init(f, ctx);
        fmpz_mod_poly_init(g, ctx);
        fmpz_mod_poly_init(h, ctx);

        fmpz_mod_poly_randtest(f, state, n_randint(state, 200), ctx);
        fmpz_mod_poly_randtest(g, state, n_randint(state, 200), ctx);
        fmpz_mod_poly_randtest(h, state, n_randint(state, 100), ctx);

        fmpz_mod_poly_resultant_hgcd(y, f, g, ctx);
        fmpz_mod_poly_resultant_hgcd(z, h, g, ctx);
        fmpz_mod_mul(y, y, z, ctx);
        fmpz_mod_poly_mul(f, f, h, ctx);
        fmpz_mod_poly_resultant_hgcd(x, f, g, ctx);

        result = (fmpz_equal(x, y));
        if (!result)
        {
            flint_printf("FAIL (res(f h, g) == res(f, g) res(h, g)):\n");
            fmpz_mod_poly_print(f, ctx), flint_printf("\n\n");
            fmpz_mod_poly_print(g, ctx), flint_printf("\n\n");
            fmpz_mod_poly_print(h, ctx), flint_printf("\n\n");
            printf("x = "); fmpz_print(x); printf("\n");
            printf("y = "); fmpz_print(y); printf("\n");
            printf("z = "); fmpz_print(z); printf("\n");
            printf("n = "); fmpz_print(n); printf("\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_mod_poly_clear(f, ctx);
        fmpz_mod_poly_clear(g, ctx);
        fmpz_mod_poly_clear(h, ctx);
        fmpz_clear(n);
        fmpz_clear(x);
        fmpz_clear(y);
        fmpz_clear(z);
    }

    fmpz_mod_ctx_clear(ctx);

    TEST_FUNCTION_END(state);
}
