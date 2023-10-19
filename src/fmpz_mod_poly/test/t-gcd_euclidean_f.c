/*
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

TEST_FUNCTION_START(fmpz_mod_poly_gcd_euclidean_f, state)
{
    int i, result;
    fmpz_mod_ctx_t ctx;

    fmpz_mod_ctx_init_ui(ctx, 2);

    /*
        Compare with the usual GCD function.

        N.B.  I checked by hand that this test shows both outcomes,
        i.e. trivial and non-trivial factors, sufficiently frequently.
     */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fmpz_t p, f;
        fmpz_mod_poly_t a, b, c, d;

        fmpz_init(f);
        fmpz_init(p);
        fmpz_randtest_unsigned(p, state, 100);
        fmpz_add_ui(p, p, 2);
        fmpz_mod_ctx_set_modulus(ctx, p);

        fmpz_mod_poly_init(a, ctx);
        fmpz_mod_poly_init(b, ctx);
        fmpz_mod_poly_init(c, ctx);
        fmpz_mod_poly_init(d, ctx);
        fmpz_mod_poly_randtest(a, state, n_randint(state, 60), ctx);
        fmpz_mod_poly_randtest(b, state, n_randint(state, 60), ctx);

        fmpz_mod_poly_gcd_euclidean_f(f, c, a, b, ctx);
        if (!fmpz_is_one(f))
        {
            result = 1;
        }
        else
        {
            fmpz_mod_poly_gcd(d, a, b, ctx);
            result = fmpz_mod_poly_equal(c, d, ctx);
        }

        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("p = "), fmpz_print(p), flint_printf("\n\n");
            flint_printf("a = "), fmpz_mod_poly_print(a, ctx), flint_printf("\n\n");
            flint_printf("b = "), fmpz_mod_poly_print(b, ctx), flint_printf("\n\n");
            flint_printf("c = "), fmpz_mod_poly_print(c, ctx), flint_printf("\n\n");
            flint_printf("d = "), fmpz_mod_poly_print(d, ctx), flint_printf("\n\n");
            flint_printf("f = "), fmpz_print(f), flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_mod_poly_clear(a, ctx);
        fmpz_mod_poly_clear(b, ctx);
        fmpz_mod_poly_clear(c, ctx);
        fmpz_mod_poly_clear(d, ctx);
        fmpz_clear(f);
        fmpz_clear(p);
    }

    fmpz_mod_ctx_clear(ctx);

    TEST_FUNCTION_END(state);
}
