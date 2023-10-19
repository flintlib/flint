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
#include "fmpz_poly.h"
#include "fmpz_mod.h"
#include "fmpz_mod_poly.h"

TEST_FUNCTION_START(fmpz_mod_poly_evaluate_fmpz, state)
{
    int i, result;
    fmpz_mod_ctx_t ctx;

    fmpz_mod_ctx_init_ui(ctx, 2);

    /* Check aliasing */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpz_t a, b, p;
        fmpz_mod_poly_t f;

        fmpz_init(p);
        fmpz_randtest_unsigned(p, state, 2 * FLINT_BITS);
        fmpz_add_ui(p, p, 2);
        fmpz_mod_ctx_set_modulus(ctx, p);

        fmpz_init(a);
        fmpz_init(b);
        fmpz_randm(a, state, p);
        fmpz_mod_poly_init(f, ctx);
        fmpz_mod_poly_randtest(f, state, n_randint(state, 100), ctx);

        fmpz_mod_poly_evaluate_fmpz(b, f, a, ctx);
        fmpz_mod_poly_evaluate_fmpz(a, f, a, ctx);

        result = (fmpz_equal(a, b));
        if (!result)
        {
            flint_printf("FAIL:\n");
            fmpz_mod_poly_print(f, ctx), flint_printf("\n\n");
            fmpz_print(a), flint_printf("\n\n");
            fmpz_print(b), flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_mod_poly_clear(f, ctx);
        fmpz_clear(a);
        fmpz_clear(b);
        fmpz_clear(p);
    }

    /* Check that the result agrees with Z[X] */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpz_t a, b, c, p;
        fmpz_mod_poly_t f;
        fmpz_poly_t g;

        fmpz_init(p);
        fmpz_randtest_unsigned(p, state, 2 * FLINT_BITS);
        fmpz_add_ui(p, p, 2);
        fmpz_mod_ctx_set_modulus(ctx, p);

        fmpz_init(a);
        fmpz_init(b);
        fmpz_init(c);
        fmpz_mod_poly_init(f, ctx);
        fmpz_poly_init(g);
        fmpz_randm(a, state, p);
        fmpz_mod_poly_randtest(f, state, n_randint(state, 100), ctx);
        fmpz_mod_poly_get_fmpz_poly(g, f, ctx);

        fmpz_mod_poly_evaluate_fmpz(b, f, a, ctx);
        fmpz_poly_evaluate_fmpz(c, g, a);
        fmpz_mod(c, c, p);

        result = (fmpz_equal(b, c));
        if (!result)
        {
            flint_printf("FAIL (cmp with fmpz_poly):\n");
            fmpz_mod_poly_print(f, ctx), flint_printf("\n\n");
            fmpz_poly_print(g), flint_printf("\n\n");
            fmpz_print(a), flint_printf("\n\n");
            fmpz_print(b), flint_printf("\n\n");
            fmpz_print(c), flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_clear(a);
        fmpz_clear(b);
        fmpz_clear(c);
        fmpz_clear(p);
        fmpz_mod_poly_clear(f, ctx);
        fmpz_poly_clear(g);
    }

    fmpz_mod_ctx_clear(ctx);

    TEST_FUNCTION_END(state);
}
