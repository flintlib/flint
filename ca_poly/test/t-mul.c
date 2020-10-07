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

    flint_printf("mul....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 300 * calcium_test_multiplier(); iter++)
    {
        ca_ctx_t ctx;
        ca_poly_t A, B, C, ABC, AB, AC, ABAC;

        /* Test A*(B+C) = A*B + A*C */
        ca_ctx_init(ctx);

        ca_poly_init(A, ctx);
        ca_poly_init(B, ctx);
        ca_poly_init(C, ctx);
        ca_poly_init(ABC, ctx);
        ca_poly_init(AB, ctx);
        ca_poly_init(AC, ctx);
        ca_poly_init(ABAC, ctx);

        ca_poly_randtest(A, state, 4, 2, 10, ctx);
        ca_poly_randtest(B, state, 4, 2, 10, ctx);
        ca_poly_randtest(C, state, 4, 2, 10, ctx);
        ca_poly_randtest(ABC, state, 4, 2, 10, ctx);
        ca_poly_randtest(AB, state, 4, 2, 10, ctx);
        ca_poly_randtest(AC, state, 4, 2, 10, ctx);
        ca_poly_randtest(ABAC, state, 4, 2, 10, ctx);

        ca_poly_add(ABC, B, C, ctx);
        if (n_randint(state, 2))
            ca_poly_mul(ABC, A, ABC, ctx);
        else
            ca_poly_mul(ABC, ABC, A, ctx);

        ca_poly_mul(AB, A, B, ctx);
        ca_poly_mul(AC, A, C, ctx);
        ca_poly_add(ABAC, AB, AC, ctx);

        if (ca_poly_check_equal(ABC, ABAC, ctx) == T_FALSE)
        {
            flint_printf("FAIL\n\n");
            flint_printf("A = "); ca_poly_print(A, ctx); flint_printf("\n");
            flint_printf("B = "); ca_poly_print(B, ctx); flint_printf("\n");
            flint_printf("C = "); ca_poly_print(C, ctx); flint_printf("\n");
            flint_printf("ABC = "); ca_poly_print(ABC, ctx); flint_printf("\n");
            flint_printf("ABAC = "); ca_poly_print(ABAC, ctx); flint_printf("\n");
            flint_abort();
        }

        /* Test A^2 = A * A */
        ca_poly_set(B, A, ctx);
        ca_poly_mul(AB, A, B, ctx);
        ca_poly_mul(AC, A, A, ctx);

        if (ca_poly_check_equal(AB, AC, ctx) == T_FALSE)
        {
            flint_printf("FAIL (squaring)\n\n");
            flint_printf("A = "); ca_poly_print(A, ctx); flint_printf("\n");
            flint_printf("AB = "); ca_poly_print(AB, ctx); flint_printf("\n");
            flint_printf("AC = "); ca_poly_print(AC, ctx); flint_printf("\n");
            flint_abort();
        }

        ca_poly_clear(A, ctx);
        ca_poly_clear(B, ctx);
        ca_poly_clear(C, ctx);
        ca_poly_clear(ABC, ctx);
        ca_poly_clear(AB, ctx);
        ca_poly_clear(AC, ctx);
        ca_poly_clear(ABAC, ctx);

        ca_ctx_clear(ctx);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}

