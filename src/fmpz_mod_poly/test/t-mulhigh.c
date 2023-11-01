/*
    Copyright (C) 2021 Daniel Schultz

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

TEST_FUNCTION_START(fmpz_mod_poly_mulhigh, state)
{
    slong i;

    /* Compare with left truncated product of a and b */
    for (i = 0; i < 20 * flint_test_multiplier(); i++)
    {
        fmpz_mod_ctx_t ctx;
        fmpz_mod_poly_t a, b, c;
        slong j, n;

        fmpz_mod_ctx_init_rand_bits(ctx, state, 200);

        fmpz_mod_poly_init(a, ctx);
        fmpz_mod_poly_init(b, ctx);
        fmpz_mod_poly_init(c, ctx);
        n = n_randint(state, 50);
        fmpz_mod_poly_randtest(b, state, n, ctx);
        fmpz_mod_poly_randtest(c, state, n, ctx);

        fmpz_mod_poly_mulhigh(a, b, c, n, ctx);
        fmpz_mod_poly_mul(b, b, c, ctx);
        for (j = 0; j < n; j++)
        {
            if (j < a->length)
                fmpz_zero(a->coeffs + j);
            if (j < b->length)
                fmpz_zero(b->coeffs + j);
        }
        _fmpz_mod_poly_normalise(a);
        _fmpz_mod_poly_normalise(b);

        FLINT_TEST(fmpz_mod_poly_equal(a, b, ctx));

        fmpz_mod_poly_clear(a, ctx);
        fmpz_mod_poly_clear(b, ctx);
        fmpz_mod_poly_clear(c, ctx);

        fmpz_mod_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
