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

TEST_FUNCTION_START(fmpz_mod_poly_compose_mod, state)
{
    int i;
    fmpz_mod_ctx_t ctx;

    fmpz_mod_ctx_init_ui(ctx, 2);

    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fmpz_mod_poly_t a, b, c, d, e;
        fmpz_t p;

        fmpz_init(p);
        fmpz_set_ui(p, n_randtest_prime(state, 0));
        fmpz_mod_ctx_set_modulus(ctx, p);

        fmpz_mod_poly_init(a, ctx);
        fmpz_mod_poly_init(b, ctx);
        fmpz_mod_poly_init(c, ctx);
        fmpz_mod_poly_init(d, ctx);
        fmpz_mod_poly_init(e, ctx);

        fmpz_mod_poly_randtest(a, state, n_randint(state, 20) + 1, ctx);
        fmpz_mod_poly_randtest(b, state, n_randint(state, 20) + 1, ctx);
        fmpz_mod_poly_randtest_not_zero(c, state, n_randint(state, 20) + 1, ctx);

        fmpz_mod_poly_compose_mod(d, a, b, c, ctx);
        fmpz_mod_poly_compose(e, a, b, ctx);
        fmpz_mod_poly_rem(e, e, c, ctx);

        if (!fmpz_mod_poly_equal(d, e, ctx))
        {
            flint_printf("FAIL (composition):\n");
            flint_printf("a:\n"); fmpz_mod_poly_print(a, ctx); flint_printf("\n");
            flint_printf("b:\n"); fmpz_mod_poly_print(b, ctx); flint_printf("\n");
            flint_printf("c:\n"); fmpz_mod_poly_print(c, ctx); flint_printf("\n");
            flint_printf("d:\n"); fmpz_mod_poly_print(d, ctx); flint_printf("\n");
            flint_printf("e:\n"); fmpz_mod_poly_print(e, ctx); flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_clear(p);
        fmpz_mod_poly_clear(a, ctx);
        fmpz_mod_poly_clear(b, ctx);
        fmpz_mod_poly_clear(c, ctx);
        fmpz_mod_poly_clear(d, ctx);
        fmpz_mod_poly_clear(e, ctx);
    }

    /* Test aliasing of res and a */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fmpz_mod_poly_t a, b, c, d;
        fmpz_t p;

        fmpz_init(p);
        fmpz_set_ui(p, n_randtest_prime(state, 0));
        fmpz_mod_ctx_set_modulus(ctx, p);

        fmpz_mod_poly_init(a, ctx);
        fmpz_mod_poly_init(b, ctx);
        fmpz_mod_poly_init(c, ctx);
        fmpz_mod_poly_init(d, ctx);

        fmpz_mod_poly_randtest(a, state, n_randint(state, 20) + 1, ctx);
        fmpz_mod_poly_randtest(b, state, n_randint(state, 20) + 1, ctx);
        fmpz_mod_poly_randtest_not_zero(c, state, n_randint(state, 20) + 1, ctx);

        fmpz_mod_poly_compose_mod(d, a, b, c, ctx);
        fmpz_mod_poly_compose_mod(a, a, b, c, ctx);

        if (!fmpz_mod_poly_equal(d, a, ctx))
        {
            flint_printf("FAIL (aliasing a):\n");
            flint_printf("a:\n"); fmpz_mod_poly_print(a, ctx); flint_printf("\n");
            flint_printf("b:\n"); fmpz_mod_poly_print(b, ctx); flint_printf("\n");
            flint_printf("c:\n"); fmpz_mod_poly_print(c, ctx); flint_printf("\n");
            flint_printf("d:\n"); fmpz_mod_poly_print(d, ctx); flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_clear(p);
        fmpz_mod_poly_clear(a, ctx);
        fmpz_mod_poly_clear(b, ctx);
        fmpz_mod_poly_clear(c, ctx);
        fmpz_mod_poly_clear(d, ctx);
    }

    /* Test aliasing of res and b */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fmpz_mod_poly_t a, b, c, d;
        fmpz_t p;

        fmpz_init(p);
        fmpz_set_ui(p, n_randtest_prime(state, 0));
        fmpz_mod_ctx_set_modulus(ctx, p);

        fmpz_mod_poly_init(a, ctx);
        fmpz_mod_poly_init(b, ctx);
        fmpz_mod_poly_init(c, ctx);
        fmpz_mod_poly_init(d, ctx);

        fmpz_mod_poly_randtest(a, state, n_randint(state, 20) + 1, ctx);
        fmpz_mod_poly_randtest(b, state, n_randint(state, 20) + 1, ctx);
        fmpz_mod_poly_randtest_not_zero(c, state, n_randint(state, 20) + 1, ctx);

        fmpz_mod_poly_compose_mod(d, a, b, c, ctx);
        fmpz_mod_poly_compose_mod(b, a, b, c, ctx);

        if (!fmpz_mod_poly_equal(d, b, ctx))
        {
            flint_printf("FAIL (aliasing b)\n");
            flint_printf("a:\n"); fmpz_mod_poly_print(a, ctx); flint_printf("\n");
            flint_printf("b:\n"); fmpz_mod_poly_print(b, ctx); flint_printf("\n");
            flint_printf("c:\n"); fmpz_mod_poly_print(c, ctx); flint_printf("\n");
            flint_printf("d:\n"); fmpz_mod_poly_print(d, ctx); flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_clear(p);
        fmpz_mod_poly_clear(a, ctx);
        fmpz_mod_poly_clear(b, ctx);
        fmpz_mod_poly_clear(c, ctx);
        fmpz_mod_poly_clear(d, ctx);
    }

    /* Test aliasing of res and c */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fmpz_mod_poly_t a, b, c, d;
        fmpz_t p;

        fmpz_init(p);
        fmpz_set_ui(p, n_randtest_prime(state, 0));
        fmpz_mod_ctx_set_modulus(ctx, p);

        fmpz_mod_poly_init(a, ctx);
        fmpz_mod_poly_init(b, ctx);
        fmpz_mod_poly_init(c, ctx);
        fmpz_mod_poly_init(d, ctx);

        fmpz_mod_poly_randtest(a, state, n_randint(state, 20) + 1, ctx);
        fmpz_mod_poly_randtest(b, state, n_randint(state, 20) + 1, ctx);
        fmpz_mod_poly_randtest_not_zero(c, state, n_randint(state, 20) + 1, ctx);

        fmpz_mod_poly_compose_mod(d, a, b, c, ctx);
        fmpz_mod_poly_compose_mod(c, a, b, c, ctx);

        if (!fmpz_mod_poly_equal(d, c, ctx))
        {
            flint_printf("FAIL (aliasing c)\n");
            flint_printf("a:\n"); fmpz_mod_poly_print(a, ctx); flint_printf("\n");
            flint_printf("b:\n"); fmpz_mod_poly_print(b, ctx); flint_printf("\n");
            flint_printf("c:\n"); fmpz_mod_poly_print(c, ctx); flint_printf("\n");
            flint_printf("d:\n"); fmpz_mod_poly_print(d, ctx); flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_clear(p);
        fmpz_mod_poly_clear(a, ctx);
        fmpz_mod_poly_clear(b, ctx);
        fmpz_mod_poly_clear(c, ctx);
        fmpz_mod_poly_clear(d, ctx);
    }

    fmpz_mod_ctx_clear(ctx);

    TEST_FUNCTION_END(state);
}
