/*
    Copyright (C) 2024 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "ulong_extras.h"
#include "fmpz.h"
#include "gr.h"
#include "gr_mat.h"

TEST_FUNCTION_START(mpn_mod, state)
{
    gr_ctx_t ZZn, VecZZn, MatZZn, ZZnx;
    fmpz_t n;
    slong iter;
    /* int flags = GR_TEST_ALWAYS_ABLE; */
    int flags = 0;

    fmpz_init(n);

    /* test prime close to the supported size limit */
    {
        if (FLINT_BITS == 64)
        {
            fmpz_ui_pow_ui(n, 2, 1024);
            fmpz_sub_ui(n, n, 105);
        }
        else
        {
            fmpz_ui_pow_ui(n, 2, 512);
            fmpz_sub_ui(n, n, 569);
        }

        GR_MUST_SUCCEED(gr_ctx_init_mpn_mod(ZZn, n));
        GR_MUST_SUCCEED(gr_ctx_set_is_field(ZZn, T_TRUE));
        gr_test_ring(ZZn, 100 * flint_test_multiplier(), flags);
        gr_ctx_clear(ZZn);
    }

    /* something close to the two-word limit */
    {
        if (FLINT_BITS == 64)
        {
            fmpz_ui_pow_ui(n, 2, 128);
            fmpz_sub_ui(n, n, 159);
        }
        else
        {
            fmpz_ui_pow_ui(n, 2, 64);
            fmpz_sub_ui(n, n, 59);
        }

        GR_MUST_SUCCEED(gr_ctx_init_mpn_mod(ZZn, n));
        GR_MUST_SUCCEED(gr_ctx_set_is_field(ZZn, T_TRUE));
        gr_test_ring(ZZn, 100 * flint_test_multiplier(), flags);
        gr_ctx_clear(ZZn);
    }

    for (iter = 0; iter < 100 * flint_test_multiplier(); iter++)
    {
        for (;;)
        {
            fmpz_randtest_not_zero(n, state, 600);
            fmpz_abs(n, n);
            if (gr_ctx_init_mpn_mod(ZZn, n) == GR_SUCCESS)
                break;
        }

        if (n_randint(state, 2))
            GR_MUST_SUCCEED(gr_ctx_set_is_field(ZZn, fmpz_is_probabprime(n) ? T_TRUE : T_FALSE));

        /* test Waksman mul */
        {
            gr_mat_t A, B, C, D;
            slong a, b, c;
            int status = GR_SUCCESS;

            a = n_randint(state, 8);
            b = n_randint(state, 8);
            c = n_randint(state, 8);

            gr_mat_init(A, a, b, ZZn);
            gr_mat_init(B, b, c, ZZn);
            gr_mat_init(C, a, c, ZZn);
            gr_mat_init(D, a, c, ZZn);

            status |= gr_mat_randtest(A, state, ZZn);
            status |= gr_mat_randtest(B, state, ZZn);
            status |= gr_mat_randtest(C, state, ZZn);
            status |= gr_mat_randtest(D, state, ZZn);

            if (b == c && n_randint(state, 2))
            {
                status |= gr_mat_set(C, A, ZZn);
                status |= gr_mat_mul(C, C, B, ZZn);
            }
            else if (a == b && n_randint(state, 2))
            {
                status |= gr_mat_set(C, B, ZZn);
                status |= gr_mat_mul(C, A, C, ZZn);
            }
            else
            {
                status |= gr_mat_mul(C, A, B, ZZn);
            }

            status |= gr_mat_mul_classical(D, A, B, ZZn);

            if (status != GR_SUCCESS && gr_mat_equal(C, D, ZZn) == T_FALSE)
            {
                flint_printf("FAIL:\n");
                gr_ctx_println(ZZn);
                flint_printf("A:\n"); gr_mat_print(A, ZZn); flint_printf("\n\n");
                flint_printf("B:\n"); gr_mat_print(B, ZZn); flint_printf("\n\n");
                flint_printf("C:\n"); gr_mat_print(C, ZZn); flint_printf("\n\n");
                flint_printf("D:\n"); gr_mat_print(D, ZZn); flint_printf("\n\n");
                flint_abort();
            }

            gr_mat_clear(A, ZZn);
            gr_mat_clear(B, ZZn);
            gr_mat_clear(C, ZZn);
            gr_mat_clear(D, ZZn);
        }

        gr_test_ring(ZZn, 5, flags);

        /* test matrices */
        if (n_randint(state, 10) == 0)
        {
            gr_ctx_init_matrix_ring(MatZZn, ZZn, 1 + n_randint(state, 5));
            gr_test_ring(MatZZn, 3, flags);
            gr_ctx_clear(MatZZn);
        }

        /* test vectors */
        if (n_randint(state, 10) == 0)
        {
            gr_ctx_init_vector_space_gr_vec(VecZZn, ZZn, 1 + n_randint(state, 5));
            gr_test_ring(VecZZn, 3, flags);
            gr_ctx_clear(VecZZn);
        }

        /* test polynomials */
        if (n_randint(state, 10) == 0)
        {
            gr_ctx_init_gr_poly(ZZnx, ZZn);
            gr_test_ring(ZZnx, 3, flags);
            gr_ctx_clear(ZZnx);
        }


        gr_ctx_clear(ZZn);
    }

    fmpz_clear(n);

    TEST_FUNCTION_END(state);
}
