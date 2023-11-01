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
#include "fmpz.h"
#include "fmpz_mod.h"
#include "fmpz_mod_mat.h"
#include "fmpz_mod_poly.h"

TEST_FUNCTION_START(fmpz_mod_mat_minpoly, state)
{
    fmpz_mod_ctx_t ctx;
    fmpz_t t;
    fmpz_mod_mat_t A, B;
    fmpz_mod_poly_t p1, p2, q, r;
    slong i, j, k, m, n;

    /* minpoly(A) divides charpoly(A) */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        m = n_randint(state, 10);
        n = m;
        fmpz_mod_ctx_init_rand_bits_prime(ctx, state, 200);

        fmpz_mod_poly_init(p1, ctx);
        fmpz_mod_poly_init(p2, ctx);
        fmpz_mod_poly_init(q, ctx);
        fmpz_mod_poly_init(r, ctx);
        fmpz_mod_mat_init(A, m, n, fmpz_mod_ctx_modulus(ctx));
        fmpz_mod_mat_randtest(A, state);

        fmpz_mod_mat_charpoly(p1, A, ctx);
        fmpz_mod_mat_minpoly(p2, A, ctx);

        fmpz_mod_poly_divrem(q, r, p1, p2, ctx);
        FLINT_TEST(fmpz_mod_poly_is_zero(r, ctx));

        fmpz_mod_mat_clear(A);
        fmpz_mod_poly_clear(p1, ctx);
        fmpz_mod_poly_clear(p2, ctx);
        fmpz_mod_poly_clear(q, ctx);
        fmpz_mod_poly_clear(r, ctx);

        fmpz_mod_ctx_clear(ctx);
    }

    /* minpoly(P^{-1}AP) == minpoly(A) */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        m = n_randint(state, 10);
        n = m;
        fmpz_mod_ctx_init_rand_bits_prime(ctx, state, 200);

        fmpz_init(t);
        fmpz_mod_poly_init(p1, ctx);
        fmpz_mod_poly_init(p2, ctx);
        fmpz_mod_mat_init(A, m, n, fmpz_mod_ctx_modulus(ctx));
        fmpz_mod_mat_init(B, m, n, fmpz_mod_ctx_modulus(ctx));
        fmpz_mod_mat_randtest(A, state);

        for (j = 0; j < n/2; j++)
        {
           for (k = 0; k < n/2; k++)
           {
              fmpz_zero(fmpz_mod_mat_entry(A, j + n/2, k));
              fmpz_zero(fmpz_mod_mat_entry(A, j, k + n/2));
              fmpz_set(fmpz_mod_mat_entry(A, j + n/2, k + n/2),
                       fmpz_mod_mat_entry(A, j, k));
           }
        }

        fmpz_mod_mat_set(B, A);
        fmpz_mod_mat_minpoly(p1, A, ctx);

        for (j = 0; j < n; j++)
        {
           fmpz_mod_set_ui(t, n_randint(state, 6) - 3, ctx);
           fmpz_mod_mat_similarity(B, n_randint(state, n), t);
        }

        fmpz_mod_mat_minpoly(p2, B, ctx);
        FLINT_TEST(fmpz_mod_poly_equal(p1, p2, ctx));

        fmpz_clear(t);
        fmpz_mod_mat_clear(A);
        fmpz_mod_mat_clear(B);
        fmpz_mod_poly_clear(p1, ctx);
        fmpz_mod_poly_clear(p2, ctx);

        fmpz_mod_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
