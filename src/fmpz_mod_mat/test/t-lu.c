/*
    Copyright (C) 2011 Fredrik Johansson
    Copyright (C) 2013 Mike Hansen

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

void perm(fmpz_mod_mat_t A, slong * P)
{
    slong i;
    fmpz ** tmp;

    if (A->mat->c == 0 || A->mat->r == 0)
        return;

    tmp = flint_malloc(sizeof(fmpz *) * A->mat->r);

    for (i = 0; i < A->mat->r; i++)
        tmp[P[i]] = A->mat->rows[i];
    for (i = 0; i < A->mat->r; i++)
        A->mat->rows[i] = tmp[i];

    flint_free(tmp);
}

void check(slong * P, fmpz_mod_mat_t LU, const fmpz_mod_mat_t A, slong rank,
           const fmpz_mod_ctx_t ctx)
{
    fmpz_mod_mat_t B, L, U;
    slong m, n, i, j;

    m = A->mat->r;
    n = A->mat->c;

    fmpz_mod_mat_init(B, m, n, fmpz_mod_ctx_modulus(ctx));
    fmpz_mod_mat_init(L, m, m, fmpz_mod_ctx_modulus(ctx));
    fmpz_mod_mat_init(U, m, n, fmpz_mod_ctx_modulus(ctx));

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
            fmpz_mod_mat_set_entry(L, i, j, fmpz_mod_mat_entry(LU, i, j));
        if (i < rank)
            fmpz_one(fmpz_mod_mat_entry(L, i, i));
        for (j = i; j < n; j++)
            fmpz_mod_mat_set_entry(U, i, j, fmpz_mod_mat_entry(LU, i, j));
    }

    fmpz_mod_mat_mul(B, L, U);
    perm(B, P);

    FLINT_TEST(fmpz_mod_mat_equal(A, B));

    fmpz_mod_mat_clear(B);
    fmpz_mod_mat_clear(L);
    fmpz_mod_mat_clear(U);
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
            fmpz_mod_mat_init(A, m, n, fmpz_mod_ctx_modulus(ctx));
            fmpz_mod_mat_randrank(A, state, r);

            if (n_randint(state, 2))
            {
                d = n_randint(state, 2 * m * n + 1);
                fmpz_mod_mat_randops(A, d, state);
            }

            fmpz_mod_mat_init_set(LU, A);
            P = flint_malloc(sizeof(slong) * m);

            rank = fmpz_mod_mat_lu(P, LU, 0);

            FLINT_TEST(r == rank);

            check(P, LU, A, rank, ctx);

            fmpz_mod_mat_clear(A);
            fmpz_mod_mat_clear(LU);
            flint_free(P);
        }

        fmpz_mod_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
