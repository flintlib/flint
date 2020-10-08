/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of Calcium.

    Calcium is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "ca_poly.h"

int main()
{
    slong iter;
    flint_rand_t state;

    flint_printf("roots....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 1000 * calcium_test_multiplier(); iter++)
    {
        ca_ctx_t ctx;
        ca_vec_t R;
        ca_poly_t A, B;

        ca_ctx_init(ctx);

        ca_vec_init(R, 0, ctx);
        ca_poly_init(A, ctx);
        ca_poly_init(B, ctx);

        ca_poly_randtest(A, state, 5, 2, 10, ctx);

        if (ca_poly_roots(R, A, ctx))
        {
            ca_poly_set_roots(B, R, ctx);

            if (A->length)
                ca_poly_mul_ca(B, B, A->coeffs + A->length - 1, ctx);

            if (ca_poly_check_equal(A, B, ctx) == T_FALSE)
            {
                flint_printf("FAIL\n\n");
                flint_printf("A = "); ca_poly_print(A, ctx); flint_printf("\n");
                flint_printf("R = "); ca_vec_print(R, ctx); flint_printf("\n");
                flint_printf("B = "); ca_poly_print(B, ctx); flint_printf("\n");
                flint_abort();
            }

            if (0)
            {
                printf("=================================================================\n\n");
                printf("EQUAL: "); truth_print(ca_poly_check_equal(A, B, ctx)); printf("\n\n");
                flint_printf("A = "); ca_poly_print(A, ctx); flint_printf("\n\n");
                flint_printf("R = "); ca_vec_print(R, ctx); flint_printf("\n\n");
                ca_poly_sub(B, A, B, ctx);
                flint_printf("B = "); ca_poly_print(B, ctx); flint_printf("\n\n");
            }

        }

        ca_vec_clear(R, ctx);
        ca_poly_clear(A, ctx);
        ca_poly_clear(B, ctx);

        ca_ctx_clear(ctx);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}

