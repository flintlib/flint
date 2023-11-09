/*
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

TEST_FUNCTION_START(fmpz_mod_poly_discriminant, state)
{
    int i, result;
    fmpz_mod_ctx_t ctx;

    fmpz_mod_ctx_init_ui(ctx, 2);

    /* Check disc(fg) == disc(f) * disc(g) * res(f, g)^2 */
    for (i = 0; i < 50 * flint_test_multiplier(); i++)
    {
        fmpz_mod_poly_t f, g, h;
        fmpz_t x, y, z, r;
        fmpz_t n;

        fmpz_init(n);
        fmpz_init(x);
        fmpz_init(y);
        fmpz_init(z);
        fmpz_init(r);

        fmpz_set_ui(n, n_randtest_prime(state, 0));
        fmpz_mod_ctx_set_modulus(ctx, n);

        fmpz_mod_poly_init(f, ctx);
        fmpz_mod_poly_init(g, ctx);
        fmpz_mod_poly_init(h, ctx);

        do {
           fmpz_mod_poly_randtest(f, state, n_randint(state, 200), ctx);
        } while (f->length < 2);
        do {
           fmpz_mod_poly_randtest(g, state, n_randint(state, 200), ctx);
        } while (g->length < 2);

        fmpz_mod_poly_discriminant(y, f, ctx);
        fmpz_mod_poly_discriminant(z, g, ctx);
        fmpz_mod_mul(y, y, z, ctx);
        fmpz_mod_poly_resultant(r, f, g, ctx);
        fmpz_mod_mul(r, r, r, ctx);
        fmpz_mod_mul(y, y, r, ctx);
        fmpz_mod_poly_mul(h, f, g, ctx);
        fmpz_mod_poly_discriminant(x, h, ctx);

        result = (fmpz_equal(x, y));
        if (!result)
        {
            flint_printf("FAIL (disc(fg) == res(f, g)^2 * disc(f) * disc(g):\n");
            fmpz_mod_poly_print(f, ctx), flint_printf("\n\n");
            fmpz_mod_poly_print(g, ctx), flint_printf("\n\n");
            fmpz_mod_poly_print(h, ctx), flint_printf("\n\n");
            flint_printf("x = "); fmpz_print(x); printf("\n");
            flint_printf("y = "); fmpz_print(y); printf("\n");
            flint_printf("z = "); fmpz_print(z); printf("\n");
            flint_printf("n = "); fmpz_print(n); printf("\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_mod_poly_clear(f, ctx);
        fmpz_mod_poly_clear(g, ctx);
        fmpz_mod_poly_clear(h, ctx);

        fmpz_clear(n);
        fmpz_clear(x);
        fmpz_clear(y);
        fmpz_clear(z);
        fmpz_clear(r);
    }

    /* Check disc(f) == 0 for length < 2 */
    for (i = 0; i < 50 * flint_test_multiplier(); i++)
    {
        fmpz_mod_poly_t f;
        fmpz_t y;
        fmpz_t n;

        fmpz_init(y);
        fmpz_init(n);

        fmpz_set_ui(n, n_randtest_prime(state, 0));
        fmpz_mod_ctx_set_modulus(ctx, n);

        fmpz_mod_poly_init(f, ctx);

        fmpz_mod_poly_randtest(f, state, 1, ctx);

        fmpz_mod_poly_discriminant(y, f, ctx);

        result = fmpz_is_zero(y);
        if (!result)
        {
            flint_printf("FAIL disc(f) == 0 for len f < 2:\n");
            fmpz_mod_poly_print(f, ctx), flint_printf("\n\n");
            flint_printf("y = "); fmpz_print(y); printf("\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_clear(n);
        fmpz_clear(y);

        fmpz_mod_poly_clear(f, ctx);
    }

    TEST_FUNCTION_END(state);
}
