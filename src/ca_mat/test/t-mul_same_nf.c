/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "ca_mat.h"

void
ca_mat_randtest_same_nf(ca_mat_t res, flint_rand_t state, const ca_t x, slong bits, slong den_bits, ca_ctx_t ctx)
{
    slong i, j;
    fmpz_t t;

    for (i = 0; i < ca_mat_nrows(res); i++)
        for (j = 0; j < ca_mat_ncols(res); j++)
            ca_randtest_same_nf(ca_mat_entry(res, i, j), state, x, bits, 1, ctx);

    fmpz_init(t);
    fmpz_randtest_not_zero(t, state, den_bits);
    ca_mat_div_fmpz(res, res, t, ctx);
    fmpz_clear(t);
}

TEST_FUNCTION_START(ca_mat_mul_same_nf, state)
{
    slong iter;

    for (iter = 0; iter < 1000 * 0.1 * flint_test_multiplier(); iter++)
    {
        ca_ctx_t ctx;
        ca_mat_t A, B, C, D;
        qqbar_t t;
        ca_t x;
        slong m, n, k;
        ca_field_ptr K;

        /* Test (A*B)*C = A*(B*C) */
        m = n_randint(state, 5);
        n = n_randint(state, 5);
        k = n_randint(state, 5);

        ca_ctx_init(ctx);

        qqbar_init(t);
        ca_init(x, ctx);

        do {
            qqbar_randtest(t, state, 8, 10);
        } while (qqbar_is_rational(t));
        ca_set_qqbar(x, t, ctx);

        ca_mat_init(A, m, n, ctx);
        ca_mat_init(B, n, k, ctx);
        ca_mat_init(C, m, k, ctx);
        ca_mat_init(D, m, k, ctx);

        ca_mat_randtest_same_nf(A, state, x, 10, 10, ctx);
        ca_mat_randtest_same_nf(B, state, x, 10, 10, ctx);
        ca_mat_randtest(C, state, 1, 5, ctx);
        ca_mat_randtest(D, state, 1, 5, ctx);

        K = _ca_mat_same_field(A, ctx);

        if (K != NULL && CA_FIELD_IS_NF(K))
        {
            if (n_randint(state, 2) && (m == n && n == k))
            {   /* test aliasing */
                ca_mat_set(C, A, ctx);
                ca_mat_mul_same_nf(C, C, B, K, ctx);
            }
            else
            {
                ca_mat_mul_same_nf(C, A, B, K, ctx);
            }

            ca_mat_mul_classical(D, A, B, ctx);

            if (ca_mat_check_equal(C, D, ctx) == T_FALSE)
            {
                flint_printf("FAIL\n\n");
                flint_printf("A = "); ca_mat_print(A, ctx); flint_printf("\n");
                flint_printf("B = "); ca_mat_print(B, ctx); flint_printf("\n");
                flint_printf("C = "); ca_mat_print(C, ctx); flint_printf("\n");
                flint_printf("D = "); ca_mat_print(D, ctx); flint_printf("\n");
                flint_abort();
            }
        }

        ca_mat_clear(A, ctx);
        ca_mat_clear(B, ctx);
        ca_mat_clear(C, ctx);
        ca_mat_clear(D, ctx);

        ca_clear(x, ctx);
        qqbar_clear(t);
        ca_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
