/*
    Copyright (C) 2021 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mod_poly.h"
#include "test_helpers.h"

int main(void)
{
    slong i;
    FLINT_TEST_INIT(state);

    flint_printf("deflate_deflation_inflate....");
    fflush(stdout);

    /* Compare with left truncated product of a and b */
    for (i = 0; i < 20 * flint_test_multiplier(); i++)
    {
        fmpz_mod_ctx_t ctx;
        fmpz_mod_poly_t a, b, c;
        ulong k;

        fmpz_mod_ctx_init_rand_bits(ctx, state, 200);

        fmpz_mod_poly_init(a, ctx);
        fmpz_mod_poly_init(b, ctx);
        fmpz_mod_poly_init(c, ctx);

        fmpz_mod_poly_randtest(a, state, n_randint(state, 50), ctx);
        fmpz_mod_poly_randtest(b, state, n_randint(state, 50), ctx);
        fmpz_mod_poly_randtest(c, state, n_randint(state, 50), ctx);

        fmpz_mod_poly_inflate(b, b, n_randint(state, 5), ctx);

        k = fmpz_mod_poly_deflation(b, ctx);
        if (k > 0)
        {
            fmpz_mod_poly_deflate(a, b, k, ctx);
            fmpz_mod_poly_inflate(c, a, k, ctx);
            FLINT_TEST(fmpz_mod_poly_equal(b, c, ctx));

            fmpz_mod_poly_inflate(a, a, k, ctx);
            FLINT_TEST(fmpz_mod_poly_equal(a, c, ctx));
        }
        else
        {
            FLINT_TEST(fmpz_mod_poly_is_zero(b, ctx));
        }

        fmpz_mod_poly_clear(a, ctx);
        fmpz_mod_poly_clear(b, ctx);
        fmpz_mod_poly_clear(c, ctx);

        fmpz_mod_ctx_clear(ctx);
    }

    FLINT_TEST_CLEANUP(state);

    flint_printf("PASS\n");
    return 0;
}

