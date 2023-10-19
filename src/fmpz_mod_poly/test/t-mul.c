/*
    Copyright (C) 2009 William Hart
    Copyright (C) 2010, 2011 Sebastian Pancratz

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

TEST_FUNCTION_START(fmpz_mod_poly_mul, state)
{
    int i, result;
    fmpz_mod_ctx_t ctx;

    fmpz_mod_ctx_init_ui(ctx, 2);

    /* Check aliasing of a and b */
    for (i = 0; i < 200 * flint_test_multiplier(); i++)
    {
        fmpz_t p;
        fmpz_mod_poly_t a, b, c;

        fmpz_init(p);
        fmpz_randtest_unsigned(p, state, 2 * FLINT_BITS);
        fmpz_add_ui(p, p, 2);
        fmpz_mod_ctx_set_modulus(ctx, p);

        fmpz_mod_poly_init(a, ctx);
        fmpz_mod_poly_init(b, ctx);
        fmpz_mod_poly_init(c, ctx);
        fmpz_mod_poly_randtest(b, state, n_randint(state, 50), ctx);
        fmpz_mod_poly_randtest(c, state, n_randint(state, 50), ctx);

        fmpz_mod_poly_mul(a, b, c, ctx);
        fmpz_mod_poly_mul(b, b, c, ctx);

        result = (fmpz_mod_poly_equal(a, b, ctx));
        if (!result)
        {
            flint_printf("FAIL:\n");
            fmpz_mod_poly_print(a, ctx), flint_printf("\n\n");
            fmpz_mod_poly_print(b, ctx), flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_mod_poly_clear(a, ctx);
        fmpz_mod_poly_clear(b, ctx);
        fmpz_mod_poly_clear(c, ctx);
        fmpz_clear(p);
    }

    /* Check aliasing of a and c */
    for (i = 0; i < 200 * flint_test_multiplier(); i++)
    {
        fmpz_t p;
        fmpz_mod_poly_t a, b, c;

        fmpz_init(p);
        fmpz_randtest_unsigned(p, state, 2 * FLINT_BITS);
        fmpz_add_ui(p, p, 2);
        fmpz_mod_ctx_set_modulus(ctx, p);

        fmpz_mod_poly_init(a, ctx);
        fmpz_mod_poly_init(b, ctx);
        fmpz_mod_poly_init(c, ctx);
        fmpz_mod_poly_randtest(b, state, n_randint(state, 50), ctx);
        fmpz_mod_poly_randtest(c, state, n_randint(state, 50), ctx);

        fmpz_mod_poly_mul(a, b, c, ctx);
        fmpz_mod_poly_mul(c, b, c, ctx);

        result = (fmpz_mod_poly_equal(a, c, ctx));
        if (!result)
        {
            flint_printf("FAIL:\n");
            fmpz_mod_poly_print(a, ctx), flint_printf("\n\n");
            fmpz_mod_poly_print(c, ctx), flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_mod_poly_clear(a, ctx);
        fmpz_mod_poly_clear(b, ctx);
        fmpz_mod_poly_clear(c, ctx);
        fmpz_clear(p);
    }

    /* Check (b*c)+(b*d) = b*(c+d) */
    for (i = 0; i < 200 * flint_test_multiplier(); i++)
    {
        fmpz_mod_poly_t a1, a2, b, c, d;
        fmpz_t p;

        fmpz_init(p);

        if (n_randint(state, 100) == 0)
        {
            flint_set_num_threads(2 + n_randint(state, 3));
            fmpz_randtest_unsigned(p, state, 40000);
        }
        else
        {
            flint_set_num_threads(1);
            fmpz_randtest_unsigned(p, state, 2 * FLINT_BITS);
        }

        fmpz_add_ui(p, p, 2);
        fmpz_mod_ctx_set_modulus(ctx, p);

        fmpz_mod_poly_init(a1, ctx);
        fmpz_mod_poly_init(a2, ctx);
        fmpz_mod_poly_init(b, ctx);
        fmpz_mod_poly_init(c, ctx);
        fmpz_mod_poly_init(d, ctx);
        fmpz_mod_poly_randtest(b, state, n_randint(state, 100), ctx);
        fmpz_mod_poly_randtest(c, state, n_randint(state, 100), ctx);
        fmpz_mod_poly_randtest(d, state, n_randint(state, 100), ctx);

        fmpz_mod_poly_mul(a1, b, c, ctx);
        fmpz_mod_poly_mul(a2, b, d, ctx);
        fmpz_mod_poly_add(a1, a1, a2, ctx);

        fmpz_mod_poly_add(c, c, d, ctx);
        fmpz_mod_poly_mul(a2, b, c, ctx);

        result = (fmpz_mod_poly_equal(a1, a2, ctx));
        if (!result)
        {
            flint_printf("FAIL:\n");
            fmpz_mod_poly_print(a1, ctx), flint_printf("\n\n");
            fmpz_mod_poly_print(a2, ctx), flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_mod_poly_clear(a1, ctx);
        fmpz_mod_poly_clear(a2, ctx);
        fmpz_mod_poly_clear(b, ctx);
        fmpz_mod_poly_clear(c, ctx);
        fmpz_mod_poly_clear(d, ctx);
        fmpz_clear(p);
    }

    fmpz_mod_ctx_clear(ctx);

    TEST_FUNCTION_END(state);
}
