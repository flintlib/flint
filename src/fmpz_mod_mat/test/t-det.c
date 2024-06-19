/*
    Copyright (C) 2024 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz.h"
#include "fmpz_mod.h"
#include "fmpz_mod_mat.h"
#include "fmpz_mod_poly.h"
#include "fmpz_mat.h"

TEST_FUNCTION_START(fmpz_mod_mat_det, state)
{
    fmpz_mod_ctx_t ctx;
    fmpz_mod_mat_t A;
    fmpz_mat_t B;
    fmpz_t d1, d2;
    slong i, n;

    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        n = n_randint(state, 7);
        fmpz_mod_ctx_init_rand_bits(ctx, state, 200);

        fmpz_init(d1);
        fmpz_init(d2);

        fmpz_mod_mat_init(A, n, n, ctx);
        fmpz_mat_init(B, n, n);

        fmpz_mod_mat_randtest(A, state, ctx);
        fmpz_mod_mat_get_fmpz_mat(B, A, ctx);

        fmpz_mod_mat_det(d1, A, ctx);
        fmpz_mat_det(d2, B);
        fmpz_mod_set_fmpz(d2, d2, ctx);

        FLINT_TEST(fmpz_equal(d1, d2));

        fmpz_mod_mat_clear(A, ctx);
        fmpz_mat_clear(B);
        fmpz_clear(d1);
        fmpz_clear(d2);

        fmpz_mod_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
