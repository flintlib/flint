/*
    Copyright (C) 2012 Sebastian Pancratz

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

TEST_FUNCTION_START(fmpz_mod_poly_gcdinv_euclidean, state)
{
    int i, result;
    fmpz_mod_ctx_t ctx;

    fmpz_mod_ctx_init_ui(ctx, 2);

    /* Generic case, most likely co-prime arguments ******************************/

    /* Check s*a == g mod b */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fmpz_t p;
        fmpz_mod_poly_t a, b, g, s, u;

        fmpz_init(p);
        fmpz_set_ui(p, n_randtest_prime(state, 0));
        fmpz_mod_ctx_set_modulus(ctx, p);

        fmpz_mod_poly_init(a, ctx);
        fmpz_mod_poly_init(b, ctx);
        fmpz_mod_poly_init(g, ctx);
        fmpz_mod_poly_init(s, ctx);
        fmpz_mod_poly_init(u, ctx);

        fmpz_mod_poly_randtest(a, state, n_randint(state, 100), ctx);
        do
            fmpz_mod_poly_randtest(b, state, n_randint(state, 100), ctx);
        while (b->length < 2);

        fmpz_mod_poly_gcdinv_euclidean(g, s, a, b, ctx);
        fmpz_mod_poly_mul(u, s, a, ctx);
        fmpz_mod_poly_rem(u, u, b, ctx);

        result = (fmpz_mod_poly_equal(g, u, ctx) || fmpz_mod_poly_is_zero(g, ctx));
        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("a = "), fmpz_mod_poly_print(a, ctx), flint_printf("\n\n");
            flint_printf("b = "), fmpz_mod_poly_print(b, ctx), flint_printf("\n\n");
            flint_printf("g = "), fmpz_mod_poly_print(g, ctx), flint_printf("\n\n");
            flint_printf("s = "), fmpz_mod_poly_print(s, ctx), flint_printf("\n\n");
            flint_printf("u = "), fmpz_mod_poly_print(u, ctx), flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_mod_poly_clear(a, ctx);
        fmpz_mod_poly_clear(b, ctx);
        fmpz_mod_poly_clear(g, ctx);
        fmpz_mod_poly_clear(s, ctx);
        fmpz_mod_poly_clear(u, ctx);
        fmpz_clear(p);
    }

    /* Special case, arguments share a factor ********************************/

    /* Check s*a == g mod b */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fmpz_t p;
        fmpz_mod_poly_t a, b, f, g, s, u;

        fmpz_init(p);
        fmpz_set_ui(p, n_randtest_prime(state, 0));
        fmpz_mod_ctx_set_modulus(ctx, p);

        fmpz_mod_poly_init(a, ctx);
        fmpz_mod_poly_init(b, ctx);
        fmpz_mod_poly_init(f, ctx);
        fmpz_mod_poly_init(g, ctx);
        fmpz_mod_poly_init(s, ctx);
        fmpz_mod_poly_init(u, ctx);

        fmpz_mod_poly_randtest(a, state, n_randint(state, 100), ctx);
        do
            fmpz_mod_poly_randtest(b, state, n_randint(state, 100), ctx);
        while (b->length < 2);
        fmpz_mod_poly_randtest_not_zero(f, state, n_randint(state, 20) + 1, ctx);
        fmpz_mod_poly_mul(a, f, a, ctx);
        fmpz_mod_poly_mul(b, f, b, ctx);

        fmpz_mod_poly_gcdinv_euclidean(g, s, a, b, ctx);
        fmpz_mod_poly_mul(u, s, a, ctx);
        fmpz_mod_poly_rem(u, u, b, ctx);

        result = (fmpz_mod_poly_equal(u, g, ctx) || fmpz_mod_poly_is_zero(g, ctx));
        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("a = "), fmpz_mod_poly_print(a, ctx), flint_printf("\n\n");
            flint_printf("b = "), fmpz_mod_poly_print(b, ctx), flint_printf("\n\n");
            flint_printf("f = "), fmpz_mod_poly_print(f, ctx), flint_printf("\n\n");
            flint_printf("g = "), fmpz_mod_poly_print(g, ctx), flint_printf("\n\n");
            flint_printf("s = "), fmpz_mod_poly_print(s, ctx), flint_printf("\n\n");
            flint_printf("u = "), fmpz_mod_poly_print(u, ctx), flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_mod_poly_clear(a, ctx);
        fmpz_mod_poly_clear(b, ctx);
        fmpz_mod_poly_clear(f, ctx);
        fmpz_mod_poly_clear(g, ctx);
        fmpz_mod_poly_clear(s, ctx);
        fmpz_mod_poly_clear(u, ctx);
        fmpz_clear(p);
    }

    fmpz_mod_ctx_clear(ctx);

    TEST_FUNCTION_END(state);
}
