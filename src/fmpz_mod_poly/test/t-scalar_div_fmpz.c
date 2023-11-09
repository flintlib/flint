/*
    Copyright (C) 2009 William Hart
    Copyright (C) 2011 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz.h"
#include "fmpz_mod.h"
#include "fmpz_mod_poly.h"

TEST_FUNCTION_START(fmpz_mod_poly_scalar_div_fmpz, state)
{
    int i, result;
    fmpz_mod_ctx_t ctx;

    fmpz_mod_ctx_init_ui(ctx, 2);

    /* Check aliasing of a and b */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpz_t n, p, g;
        fmpz_mod_poly_t a, b;

        fmpz_init(p);
        fmpz_init(g);
        fmpz_randtest_unsigned(p, state, 2 * FLINT_BITS);
        fmpz_add_ui(p, p, 2);
        fmpz_mod_ctx_set_modulus(ctx, p);

        fmpz_init(n);
        do {
           fmpz_randtest(n, state, 200);
           fmpz_gcd(g, n, p);
        } while (!fmpz_is_one(g));

        fmpz_mod_poly_init(a, ctx);
        fmpz_mod_poly_init(b, ctx);
        fmpz_mod_poly_randtest(a, state, n_randint(state, 100), ctx);

        fmpz_mod_poly_scalar_div_fmpz(b, a, n, ctx);
        fmpz_mod_poly_scalar_div_fmpz(a, a, n, ctx);

        result = (fmpz_mod_poly_equal(a, b, ctx));
        if (!result)
        {
            flint_printf("FAIL:\n");
            fmpz_mod_poly_print(a, ctx), flint_printf("\n\n");
            fmpz_mod_poly_print(b, ctx), flint_printf("\n\n");
            fmpz_print(n), flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_clear(n);
        fmpz_mod_poly_clear(a, ctx);
        fmpz_mod_poly_clear(b, ctx);
        fmpz_clear(g);
        fmpz_clear(p);
    }

    /* Check (a*n)/n = a */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpz_t n, p, g;
        fmpz_mod_poly_t a, b, c;

        fmpz_init(p);
        fmpz_init(g);
        fmpz_randtest_unsigned(p, state, 2 * FLINT_BITS);
        fmpz_add_ui(p, p, 2);
        fmpz_mod_ctx_set_modulus(ctx, p);

        fmpz_init(n);
        do {
           fmpz_randtest(n, state, 200);
           fmpz_gcd(g, n, p);
        } while (!fmpz_is_one(g));

        fmpz_mod_poly_init(a, ctx);
        fmpz_mod_poly_init(b, ctx);
        fmpz_mod_poly_init(c, ctx);
        fmpz_mod_poly_randtest(a, state, n_randint(state, 100), ctx);

        fmpz_mod_poly_scalar_mul_fmpz(b, a, n, ctx);
        fmpz_mod_poly_scalar_div_fmpz(c, b, n, ctx);

        result = (fmpz_mod_poly_equal(a, c, ctx));
        if (!result)
        {
            flint_printf("FAIL:\n");
            fmpz_mod_poly_print(a, ctx), flint_printf("\n\n");
            fmpz_mod_poly_print(b, ctx), flint_printf("\n\n");
            fmpz_mod_poly_print(c, ctx), flint_printf("\n\n");
            fmpz_print(n), flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_clear(n);
        fmpz_mod_poly_clear(a, ctx);
        fmpz_mod_poly_clear(b, ctx);
        fmpz_mod_poly_clear(c, ctx);
        fmpz_clear(g);
        fmpz_clear(p);
    }

    fmpz_mod_ctx_clear(ctx);

    TEST_FUNCTION_END(state);
}
