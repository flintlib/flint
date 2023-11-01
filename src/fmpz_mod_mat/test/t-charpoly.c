/*
    Copyright (C) 2010 Fredrik Johansson
    Copyright (C) 2013 Mike Hansen
    Copyright (C) 2015 William Hart
    Copyright (C) 2021 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz_mod.h"
#include "fmpz_mod_mat.h"
#include "fmpz_mod_poly.h"

TEST_FUNCTION_START(fmpz_mod_mat_charpoly, state)
{
    fmpz_mod_ctx_t ctx;
    fmpz_mod_mat_t A, B, C;
    fmpz_mod_poly_t p1, p2;
    slong i, n;

    /* charpoly(AB) == charpoly(BA) */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        n = n_randint(state, 10);
        fmpz_mod_ctx_init_rand_bits(ctx, state, 200);

        fmpz_mod_poly_init(p1, ctx);
        fmpz_mod_poly_init(p2, ctx);
        fmpz_mod_mat_init(A, n, n, fmpz_mod_ctx_modulus(ctx));
        fmpz_mod_mat_init(B, n, n, fmpz_mod_ctx_modulus(ctx));
        fmpz_mod_mat_init(C, n, n, fmpz_mod_ctx_modulus(ctx));
        fmpz_mod_mat_randtest(A, state);
        fmpz_mod_mat_randtest(B, state);

        fmpz_mod_mat_mul(C, A, B);
        fmpz_mod_mat_charpoly(p1, C, ctx);
        fmpz_mod_mat_mul(C, B, A);
        fmpz_mod_mat_charpoly(p2, C, ctx);
        FLINT_TEST(fmpz_mod_poly_degree(p1, ctx) == n);
        FLINT_TEST(fmpz_mod_poly_equal(p1, p2, ctx));

        fmpz_mod_mat_clear(A);
        fmpz_mod_mat_clear(B);
        fmpz_mod_mat_clear(C);
        fmpz_mod_poly_clear(p1, ctx);
        fmpz_mod_poly_clear(p2, ctx);

        fmpz_mod_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
