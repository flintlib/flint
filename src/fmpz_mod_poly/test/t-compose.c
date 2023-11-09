/*
    Copyright (C) 2011, 2010 Sebastian Pancratz
    Copyright (C) 2009 William Hart

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

TEST_FUNCTION_START(fmpz_mod_poly_compose, state)
{
    int i, result;
    fmpz_mod_ctx_t ctx;

    fmpz_mod_ctx_init_ui(ctx, 2);

    /* Check aliasing of a and c */
    for (i = 0; i < 10 * flint_test_multiplier(); i++)
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
        fmpz_mod_poly_randtest(a, state, n_randint(state, 80), ctx);
        fmpz_mod_poly_randtest(b, state, n_randint(state, 30), ctx);

        fmpz_mod_poly_compose(c, a, b, ctx);
        fmpz_mod_poly_compose(a, a, b, ctx);

        result = (fmpz_mod_poly_equal(a, c, ctx));
        if (!result)
        {
            flint_printf("FAIL:\n");
            fmpz_mod_poly_print(a, ctx), flint_printf("\n\n");
            fmpz_mod_poly_print(b, ctx), flint_printf("\n\n");
            fmpz_mod_poly_print(c, ctx), flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_mod_poly_clear(a, ctx);
        fmpz_mod_poly_clear(b, ctx);
        fmpz_mod_poly_clear(c, ctx);
        fmpz_clear(p);
    }

    /* Check aliasing of b and c */
    for (i = 0; i < 10 * flint_test_multiplier(); i++)
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
        fmpz_mod_poly_randtest(a, state, n_randint(state, 80), ctx);
        fmpz_mod_poly_randtest(b, state, n_randint(state, 30), ctx);

        fmpz_mod_poly_compose(c, a, b, ctx);
        fmpz_mod_poly_compose(b, a, b, ctx);

        result = (fmpz_mod_poly_equal(b, c, ctx));
        if (!result)
        {
            flint_printf("FAIL:\n");
            fmpz_mod_poly_print(a, ctx), flint_printf("\n\n");
            fmpz_mod_poly_print(b, ctx), flint_printf("\n\n");
            fmpz_mod_poly_print(c, ctx), flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_mod_poly_clear(a, ctx);
        fmpz_mod_poly_clear(b, ctx);
        fmpz_mod_poly_clear(c, ctx);
        fmpz_clear(p);
    }

    /* Compare with composition over Z[X] */
    for (i = 0; i < 10 * flint_test_multiplier(); i++)
    {
        fmpz_t p;
        fmpz_mod_poly_t a, b, c, d;
        fmpz_poly_t A, B, C;

        fmpz_init(p);
        fmpz_randtest_unsigned(p, state, 2 * FLINT_BITS);
        fmpz_add_ui(p, p, 2);
        fmpz_mod_ctx_set_modulus(ctx, p);

        fmpz_mod_poly_init(a, ctx);
        fmpz_mod_poly_init(b, ctx);
        fmpz_mod_poly_init(c, ctx);
        fmpz_mod_poly_init(d, ctx);
        fmpz_mod_poly_randtest(a, state, n_randint(state, 80), ctx);
        fmpz_mod_poly_randtest(b, state, n_randint(state, 30), ctx);

        fmpz_poly_init(A);
        fmpz_poly_init(B);
        fmpz_poly_init(C);
        fmpz_mod_poly_get_fmpz_poly(A, a, ctx);
        fmpz_mod_poly_get_fmpz_poly(B, b, ctx);

        fmpz_mod_poly_compose(c, a, b, ctx);
        fmpz_poly_compose(C, A, B);
        fmpz_mod_poly_set_fmpz_poly(d, C, ctx);

        result = (fmpz_mod_poly_equal(c, d, ctx));
        if (!result)
        {
            flint_printf("FAIL:\n");
            fmpz_mod_poly_print(a, ctx), flint_printf("\n\n");
            fmpz_mod_poly_print(b, ctx), flint_printf("\n\n");
            fmpz_mod_poly_print(c, ctx), flint_printf("\n\n");
            fmpz_mod_poly_print(d, ctx), flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_mod_poly_clear(a, ctx);
        fmpz_mod_poly_clear(b, ctx);
        fmpz_mod_poly_clear(c, ctx);
        fmpz_mod_poly_clear(d, ctx);
        fmpz_poly_clear(A);
        fmpz_poly_clear(B);
        fmpz_poly_clear(C);
        fmpz_clear(p);
    }

    fmpz_mod_ctx_clear(ctx);

    TEST_FUNCTION_END(state);
}
