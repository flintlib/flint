/*
    Copyright (C) 2011 Fredrik Johansson
    Copyright (C) 2013 Mike Hansen

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <string.h>
#include "test_helpers.h"
#include "fmpz.h"
#include "fmpz_mod.h"
#include "fmpz_mod_mat.h"

void perm(fmpz_mod_mat_t A, slong * P)
{
    slong i;
    fmpz * tmp;

    if (A->c == 0 || A->r == 0)
        return;

    tmp = flint_malloc(sizeof(fmpz) * A->r * A->c);

    for (i = 0; i < A->r; i++) memcpy(tmp + P[i] * A->c, fmpz_mat_entry(A, i, 0), A->c * sizeof(fmpz));
    for (i = 0; i < A->r; i++) memcpy(fmpz_mat_entry(A, i, 0), tmp + i * A->c, A->c * sizeof(fmpz));

    flint_free(tmp);
}

void check(slong * P, fmpz_mod_mat_t LU, const fmpz_mod_mat_t A, slong rank,
           const fmpz_mod_ctx_t ctx)
{
    fmpz_mod_mat_t B, L, U;
    slong m, n, i, j;

    m = A->r;
    n = A->c;

    fmpz_mod_mat_init(B, m, n, ctx);
    fmpz_mod_mat_init(L, m, m, ctx);
    fmpz_mod_mat_init(U, m, n, ctx);

    rank = FLINT_ABS(rank);

    for (i = rank; i < FLINT_MIN(m, n); i++)
    {
        for (j = i; j < n; j++)
        {
            FLINT_TEST(fmpz_is_zero(fmpz_mod_mat_entry(LU, i, j)));
        }
    }

    for (i = 0; i < m; i++)
    {
        for (j = 0; j < FLINT_MIN(i, n); j++)
            fmpz_mod_mat_set_entry(L, i, j, fmpz_mod_mat_entry(LU, i, j), ctx);
        if (i < rank)
            fmpz_one(fmpz_mod_mat_entry(L, i, i));
        for (j = i; j < n; j++)
            fmpz_mod_mat_set_entry(U, i, j, fmpz_mod_mat_entry(LU, i, j), ctx);
    }

    fmpz_mod_mat_mul(B, L, U, ctx);
    perm(B, P);

    FLINT_TEST(fmpz_mod_mat_equal(A, B, ctx));

    fmpz_mod_mat_clear(B, ctx);
    fmpz_mod_mat_clear(L, ctx);
    fmpz_mod_mat_clear(U, ctx);
}

TEST_FUNCTION_START(fmpz_mod_mat_lu, state)
{
    slong i;

    for (i = 0; i < 10 * flint_test_multiplier(); i++)
    {
        fmpz_mod_ctx_t ctx;
        fmpz_mod_mat_t A, LU;
        slong m, n, r, d, rank;
        slong *P;

        fmpz_mod_ctx_init_rand_bits_prime(ctx, state, 200);

        m = n_randint(state, 20);
        n = n_randint(state, 20);

        for (r = 0; r <= FLINT_MIN(m, n); r++)
        {
            fmpz_mod_mat_init(A, m, n, ctx);
            fmpz_mod_mat_randrank(A, state, r, ctx);

            if (n_randint(state, 2))
            {
                d = n_randint(state, 2 * m * n + 1);
                fmpz_mod_mat_randops(A, state, d, ctx);
            }

            fmpz_mod_mat_init_set(LU, A, ctx);
            P = flint_malloc(sizeof(slong) * m);

            rank = fmpz_mod_mat_lu(P, LU, 0, ctx);

            FLINT_TEST(r == rank);

            check(P, LU, A, rank, ctx);

            fmpz_mod_mat_clear(A, ctx);
            fmpz_mod_mat_clear(LU, ctx);
            flint_free(P);
        }

        fmpz_mod_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
