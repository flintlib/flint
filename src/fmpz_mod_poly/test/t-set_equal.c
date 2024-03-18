/*
    Copyright (C) 2009 William Hart
    Copyright (C) 2011 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz.h"
#include "fmpz_mod.h"
#include "fmpz_mod_poly.h"

TEST_FUNCTION_START(fmpz_mod_poly_set_equal, state)
{
    int i, result;
    fmpz_mod_ctx_t ctx;

    fmpz_mod_ctx_init_ui(ctx, 2);

    /* Equal polynomials */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpz_t p;
        fmpz_mod_poly_t a, b;

        fmpz_init(p);
        fmpz_randtest_unsigned(p, state, 2 * FLINT_BITS);
        fmpz_add_ui(p, p, 2);
        fmpz_mod_ctx_set_modulus(ctx, p);

        fmpz_mod_poly_init(a, ctx);
        fmpz_mod_poly_init(b, ctx);
        fmpz_mod_poly_randtest(a, state, n_randint(state, 100), ctx);

        fmpz_mod_poly_set(b, a, ctx);

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

        fmpz_clear(p);
    }

    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpz_t p;
        fmpz_mod_poly_t a, b;
        slong coeff = n_randint(state, 100);
        fmpz_t x;

        fmpz_init(p);
        fmpz_randtest_unsigned(p, state, 2 * FLINT_BITS);
        fmpz_add_ui(p, p, 2);
        fmpz_mod_ctx_set_modulus(ctx, p);

        fmpz_init(x);
        fmpz_mod_poly_init(a, ctx);
        fmpz_mod_poly_init(b, ctx);
        fmpz_mod_poly_randtest(a, state, n_randint(state, 100), ctx);

        fmpz_mod_poly_set(b, a, ctx);

        fmpz_mod_poly_get_coeff_fmpz(x, b, coeff, ctx);
        fmpz_sub_ui(x, x, 1);
        fmpz_mod_poly_set_coeff_fmpz(b, coeff, x, ctx);

        result = (!fmpz_mod_poly_equal(a, b, ctx));
        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("p = "), fmpz_print(p), flint_printf("\n\n");
            fmpz_mod_poly_print(a, ctx), flint_printf("\n\n");
            fmpz_mod_poly_print(b, ctx), flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_clear(x);
        fmpz_mod_poly_clear(a, ctx);
        fmpz_mod_poly_clear(b, ctx);

        fmpz_clear(p);
    }

    fmpz_mod_ctx_clear(ctx);

    TEST_FUNCTION_END(state);
}
