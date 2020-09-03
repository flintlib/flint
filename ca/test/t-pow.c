/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of Calcium.

    Calcium is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "ca.h"

int main()
{
    slong iter;
    flint_rand_t state;

    flint_printf("pow....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 1000 * calcium_test_multiplier(); iter++)
    {
        ca_ctx_t ctx;
        ca_t x, a, b, xa, xb, xaxb, ab, xab;
        truth_t equal;

        ca_ctx_init(ctx);
        ca_init(x, ctx);
        ca_init(a, ctx);
        ca_init(b, ctx);
        ca_init(xa, ctx);
        ca_init(xb, ctx);
        ca_init(xaxb, ctx);
        ca_init(ab, ctx);
        ca_init(xab, ctx);

        /* x^a * x^b = x^(a+b) */

        do {
            ca_randtest(x, state, 5, 5, ctx);
        } while (ca_check_is_zero(x, ctx) != T_FALSE);

        ca_randtest(a, state, 5, 5, ctx);
        ca_randtest(b, state, 5, 5, ctx);

        ca_pow(xa, x, a, ctx);
        ca_pow(xb, x, b, ctx);
        ca_mul(xaxb, xa, xb, ctx);

        ca_add(ab, a, b, ctx);
        ca_pow(xab, x, ab, ctx);

        equal = ca_check_equal(xab, xaxb, ctx);

        if (equal == T_FALSE)
        {
            flint_printf("FAIL x^a * x^b = x^(a+b)\n\n");
            flint_printf("x = "); ca_print(x, ctx); flint_printf(" ~= "); ca_printn(x, 10, ARB_STR_NO_RADIUS, ctx); flint_printf("\n\n");
            flint_printf("a = "); ca_print(a, ctx); flint_printf(" ~= "); ca_printn(a, 10, ARB_STR_NO_RADIUS, ctx); flint_printf("\n\n");
            flint_printf("b = "); ca_print(b, ctx); flint_printf(" ~= "); ca_printn(b, 10, ARB_STR_NO_RADIUS, ctx); flint_printf("\n\n");
            flint_printf("xa = "); ca_print(xa, ctx); flint_printf(" ~= "); ca_printn(xa, 10, ARB_STR_NO_RADIUS, ctx); flint_printf("\n\n");
            flint_printf("xb = "); ca_print(xb, ctx); flint_printf(" ~= "); ca_printn(xb, 10, ARB_STR_NO_RADIUS, ctx); flint_printf("\n\n");
            flint_printf("xaxb = "); ca_print(xaxb, ctx); flint_printf(" ~= "); ca_printn(xaxb, 10, ARB_STR_NO_RADIUS, ctx); flint_printf("\n\n");
            flint_printf("ab = "); ca_print(ab, ctx); flint_printf(" ~= "); ca_printn(ab, 10, ARB_STR_NO_RADIUS, ctx); flint_printf("\n\n");
            flint_printf("xab = "); ca_print(xab, ctx); flint_printf(" ~= "); ca_printn(xab, 10, ARB_STR_NO_RADIUS, ctx); flint_printf("\n\n");
            flint_abort();
        }

        ca_clear(x, ctx);
        ca_clear(a, ctx);
        ca_clear(b, ctx);
        ca_clear(xa, ctx);
        ca_clear(xb, ctx);
        ca_clear(xaxb, ctx);
        ca_clear(ab, ctx);
        ca_clear(xab, ctx);
        ca_ctx_clear(ctx);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}

