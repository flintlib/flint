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

TEST_FUNCTION_START(fmpz_mod_poly_inv_series, state)
{
    int i, result;

    /* Check Q^{-1} * Q is congruent 1 mod t^n */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fmpz_t p;
        fmpz_mod_ctx_t ctx;
        fmpz_mod_poly_t a, b, c, one;
        slong n = n_randint(state, 80) + 1;

        fmpz_init(p);
        fmpz_randtest_unsigned(p, state, 2 * FLINT_BITS);
        fmpz_add_ui(p, p, 2);
        fmpz_mod_ctx_init(ctx, p);

        fmpz_mod_poly_init(a, ctx);
        fmpz_mod_poly_init(b, ctx);
        fmpz_mod_poly_init(c, ctx);
        fmpz_mod_poly_init(one, ctx);

        fmpz_mod_poly_randtest_not_zero(a, state, n_randint(state, 80) + 1, ctx);
        {
            fmpz_t d;

            fmpz_init(d);
            fmpz_gcd(d, a->coeffs, p);
            while (!fmpz_is_one(d))
            {
                fmpz_randm(a->coeffs, state, p);
                fmpz_gcd(d, a->coeffs, p);
            }
            fmpz_clear(d);
        }

        fmpz_mod_poly_set_ui(one, 1, ctx);

        fmpz_mod_poly_inv_series(b, a, n, ctx);
        fmpz_mod_poly_mullow(c, a, b, n, ctx);

        result = (fmpz_mod_poly_equal(c, one, ctx));
        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("a = "), fmpz_mod_poly_print(a, ctx), flint_printf("\n\n");
            flint_printf("b = "), fmpz_mod_poly_print(b, ctx), flint_printf("\n\n");
            flint_printf("c = "), fmpz_mod_poly_print(c, ctx), flint_printf("\n\n");
            flint_printf("p = "), fmpz_print(p), flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_mod_poly_clear(a, ctx);
        fmpz_mod_poly_clear(b, ctx);
        fmpz_mod_poly_clear(c, ctx);
        fmpz_mod_poly_clear(one, ctx);
        fmpz_clear(p);
        fmpz_mod_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
