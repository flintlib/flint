/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "ca_poly.h"

TEST_FUNCTION_START(ca_poly_factor_squarefree, state)
{
    slong iter;

    for (iter = 0; iter < 1000 * 0.1 * flint_test_multiplier(); iter++)
    {
        ca_ctx_t ctx;
        ca_poly_t A, B, C, P, Q, G1, G2, G3;
        ca_t c;
        ca_poly_vec_t F;
        int success;
        ulong ea, eb, ec, maxexp;
        ulong * exp;
        slong i, j;

        ca_ctx_init(ctx);

        ca_poly_init(A, ctx);
        ca_poly_init(B, ctx);
        ca_poly_init(C, ctx);
        ca_poly_init(G1, ctx);
        ca_poly_init(G2, ctx);
        ca_poly_init(G3, ctx);
        ca_poly_init(P, ctx);
        ca_poly_init(Q, ctx);
        ca_poly_vec_init(F, 0, ctx);
        ca_init(c, ctx);

        ea = 1 + n_randint(state, 2);
        eb = 1 + n_randint(state, 2);
        ec = 1 + n_randint(state, 4);

        do {
            if (n_randint(state, 2) == 0)
            {
                ca_poly_randtest(A, state, 2, 0, 3, ctx);
                ca_poly_randtest_rational(B, state, 2, 3, ctx);
                ca_poly_randtest_rational(C, state, 2, 3, ctx);
            }
            else
            {
                ca_poly_randtest_rational(A, state, 4, 3, ctx);
                ca_poly_randtest_rational(B, state, 4, 3, ctx);
                ca_poly_randtest_rational(C, state, 4, 3, ctx);
            }
        }
        while (A->length < 2 || B->length < 2 || C->length < 2 ||
               !ca_poly_gcd(G1, A, B, ctx) || (ca_poly_check_is_one(G1, ctx) != T_TRUE) ||
               !ca_poly_gcd(G2, A, C, ctx) || (ca_poly_check_is_one(G2, ctx) != T_TRUE) ||
               !ca_poly_gcd(G3, B, C, ctx) || (ca_poly_check_is_one(G3, ctx) != T_TRUE));

        ca_poly_one(P, ctx);
        for (i = 0; i < ea; i++)
            ca_poly_mul(P, P, A, ctx);
        for (i = 0; i < eb; i++)
            ca_poly_mul(P, P, B, ctx);
        for (i = 0; i < ec; i++)
            ca_poly_mul(P, P, C, ctx);

        exp = flint_malloc(sizeof(ulong) * P->length);

        success = ca_poly_factor_squarefree(c, F, exp, P, ctx);

        if (success)
        {
            ca_poly_one(Q, ctx);
            for (i = 0; i < F->length; i++)
                for (j = 0; j < exp[i]; j++)
                    ca_poly_mul(Q, Q, F->entries + i, ctx);

            ca_poly_mul_ca(Q, Q, c, ctx);

            maxexp = 0;
            for (i = 0; i < F->length; i++)
                maxexp = FLINT_MAX(maxexp, exp[i]);

            if (ca_poly_check_equal(P, Q, ctx) == T_FALSE || maxexp < FLINT_MAX(FLINT_MAX(ea, eb), ec))
            {
                flint_printf("FAIL (product)\n\n");
                flint_printf("A = "); ca_poly_print(A, ctx); flint_printf("\n");
                flint_printf("ea = %wu\n\n", ea);
                flint_printf("B = "); ca_poly_print(B, ctx); flint_printf("\n");
                flint_printf("eb = %wu\n\n", eb);
                flint_printf("C = "); ca_poly_print(C, ctx); flint_printf("\n");
                flint_printf("ec = %wu\n\n", ec);
                flint_printf("P = "); ca_poly_print(P, ctx); flint_printf("\n");
                flint_printf("Q = "); ca_poly_print(Q, ctx); flint_printf("\n");

                for (i = 0; i < F->length; i++)
                {
                    flint_printf("Multiplicity %wu: ", exp[i]);
                    ca_poly_print(F->entries + i, ctx);
                }

                flint_abort();
            }
        }

        flint_free(exp);

        ca_poly_clear(A, ctx);
        ca_poly_clear(B, ctx);
        ca_poly_clear(C, ctx);
        ca_poly_clear(G1, ctx);
        ca_poly_clear(G2, ctx);
        ca_poly_clear(G3, ctx);
        ca_poly_clear(P, ctx);
        ca_poly_clear(Q, ctx);
        ca_poly_vec_clear(F, ctx);
        ca_clear(c, ctx);

        ca_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
