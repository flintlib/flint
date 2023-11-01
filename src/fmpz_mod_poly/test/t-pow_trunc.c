/*
    Copyright (C) 2010 William Hart
    Copyright (C) 2012 Lina Kulakova

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

TEST_FUNCTION_START(fmpz_mod_poly_pow_trunc, state)
{
    int i, result;
    fmpz_mod_ctx_t ctx;

    fmpz_mod_ctx_init_ui(ctx, 2);

    /* Check aliasing */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpz_mod_poly_t a, b, c;
        fmpz_t p;
        slong e, trunc;

        fmpz_init(p);
        fmpz_randtest_unsigned(p, state, 2 * FLINT_BITS);
        fmpz_add_ui(p, p, 2);
        fmpz_mod_ctx_set_modulus(ctx, p);

        fmpz_mod_poly_init(a, ctx);
        fmpz_mod_poly_init(b, ctx);
        fmpz_mod_poly_init(c, ctx);

        fmpz_mod_poly_randtest(a, state, n_randint(state, 30), ctx);
        e = n_randint(state, 20);
        trunc = n_randint(state, 30);

        fmpz_mod_poly_set(c, a, ctx);

        fmpz_mod_poly_pow_trunc(b, a, e, trunc, ctx);
        fmpz_mod_poly_pow_trunc(c, c, e, trunc, ctx);

        result = (fmpz_mod_poly_equal(b, c, ctx));
        if (!result)
        {
            flint_printf("FAIL aliasing:\n");
            flint_printf("a->length = %wd, exp = %wd, trunc = %wd\n",
                a->length, e, trunc);
            flint_printf("a:\n"); fmpz_mod_poly_print(a, ctx), flint_printf("\n\n");
            flint_printf("b:\n"); fmpz_mod_poly_print(b, ctx), flint_printf("\n\n");
            flint_printf("c:\n"); fmpz_mod_poly_print(c, ctx), flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_clear(p);
        fmpz_mod_poly_clear(a, ctx);
        fmpz_mod_poly_clear(b, ctx);
        fmpz_mod_poly_clear(c, ctx);
    }

    /* Check powering against naive method */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpz_mod_poly_t a, b, c;
        fmpz_t p;
        slong e, trunc;

        fmpz_init(p);
        fmpz_randtest_unsigned(p, state, 2 * FLINT_BITS);
        fmpz_add_ui(p, p, 2);
        fmpz_mod_ctx_set_modulus(ctx, p);

        fmpz_mod_poly_init(a, ctx);
        fmpz_mod_poly_init(b, ctx);
        fmpz_mod_poly_init(c, ctx);
        fmpz_mod_poly_randtest(a, state, n_randint(state, 30), ctx);
        e = n_randint(state, 20);
        trunc = n_randint(state, 30);

        fmpz_mod_poly_pow_trunc(b, a, e, trunc, ctx);
        fmpz_mod_poly_pow(c, a, e, ctx);
        fmpz_mod_poly_truncate(c, trunc, ctx);

        result = (fmpz_mod_poly_equal(b, c, ctx)
            || (a->length == 0 && e == 0 && c->length == 1 && c->coeffs[0] == 1));
        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("a->length = %wd, exp = %wd, trunc = %wd\n",
                a->length, e, trunc);
            flint_printf("a:\n"); fmpz_mod_poly_print(a, ctx), flint_printf("\n\n");
            flint_printf("b:\n"); fmpz_mod_poly_print(b, ctx), flint_printf("\n\n");
            flint_printf("c:\n"); fmpz_mod_poly_print(c, ctx), flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_clear(p);
        fmpz_mod_poly_clear(a, ctx);
        fmpz_mod_poly_clear(b, ctx);
        fmpz_mod_poly_clear(c, ctx);
    }

    fmpz_mod_ctx_clear(ctx);

    TEST_FUNCTION_END(state);
}
