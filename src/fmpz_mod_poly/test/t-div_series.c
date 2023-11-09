/*
    Copyright (C) 2009 William Hart
    Copyright (C) 2011 Sebastian Pancratz
    Copyright (C) 2014 William Hart

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

TEST_FUNCTION_START(fmpz_mod_poly_div_series, state)
{
    slong i, j;
    fmpz_mod_ctx_t ctx;

    fmpz_mod_ctx_init_ui(ctx, 2);

    /* Check A*B^{-1} * B is congruent A mod t^n */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpz_t p;
        fmpz_mod_poly_t a, b, c, d;
        slong n = n_randint(state, 80) + 1;

        fmpz_init(p);
        fmpz_init_set_ui(p, n_randtest_prime(state, 0));
        fmpz_mod_ctx_set_modulus(ctx, p);

        fmpz_mod_poly_init(a, ctx);
        fmpz_mod_poly_init(b, ctx);
        fmpz_mod_poly_init(c, ctx);
        fmpz_mod_poly_init(d, ctx);

        for (j = 0; j < 10; j++)
        {
            fmpz_mod_poly_randtest(a, state, n_randint(state, 80), ctx);
            fmpz_mod_poly_randtest_not_zero(b, state, n_randint(state, 80) + 1, ctx);
            if (fmpz_is_zero(b->coeffs + 0))
                fmpz_add_ui(b->coeffs + 0, b->coeffs + 0, 1);

            fmpz_mod_poly_div_series(c, a, b, n, ctx);
            fmpz_mod_poly_mullow(d, c, b, n, ctx);

            if (!fmpz_mod_poly_equal_trunc(d, a, n, ctx))
            {
                flint_printf("FAIL:\n");
                flint_printf("a = "), fmpz_mod_poly_print(a, ctx), flint_printf("\n\n");
                flint_printf("b = "), fmpz_mod_poly_print(b, ctx), flint_printf("\n\n");
                flint_printf("c = "), fmpz_mod_poly_print(c, ctx), flint_printf("\n\n");
                flint_printf("d = "), fmpz_mod_poly_print(d, ctx), flint_printf("\n\n");
                flint_printf("n = %wd\n", n);
                flint_printf("p = "), fmpz_print(p), flint_printf("\n\n");
                fflush(stdout);
                flint_abort();
            }
        }

        fmpz_mod_poly_clear(a, ctx);
        fmpz_mod_poly_clear(b, ctx);
        fmpz_mod_poly_clear(c, ctx);
        fmpz_mod_poly_clear(d, ctx);
        fmpz_clear(p);
    }

    fmpz_mod_ctx_clear(ctx);

    TEST_FUNCTION_END(state);
}
