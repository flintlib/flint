/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz_mod_mpoly_q.h"

TEST_FUNCTION_START(fmpz_mod_mpoly_q_add, state)
{
    slong iter;
    
    fmpz_t m;
    fmpq_t xq;
    fmpz_mod_mpoly_q_t t, x, y, z;
    fmpz_mod_mpoly_ctx_t ctx;
    fmpz_init_set_si(m, 16552213);
    fmpz_mod_mpoly_ctx_init(ctx, 2, ORD_LEX, m);
    fmpz_init_set_si(m, 12225573);
    fmpz_mod_mpoly_q_init(t, ctx);
    fmpz_mod_mpoly_q_init(x, ctx);
    fmpz_mod_mpoly_q_init(y, ctx);
    fmpz_mod_mpoly_q_init(z, ctx);
    fmpz_mod_mpoly_q_set_str_pretty(x, "14793323*x1^2", NULL, ctx);
    fmpz_mod_mpoly_q_set_str_pretty(y, "15598365", NULL, ctx);
    fmpz_mod_mpoly_q_set_str_pretty(z, "12182006", NULL, ctx);
    fmpz_init_set(fmpq_numref(xq),fmpz_mod_mpoly_q_numref(x)->coeffs);fmpz_init_set(fmpq_denref(xq),fmpz_mod_mpoly_q_denref(x)->coeffs);
    fmpz_print(fmpq_numref(xq));flint_printf("/");fmpz_print(fmpq_denref(xq));flint_printf("\n");
    fmpz_mod_mpoly_q_sub(t,x,y,ctx);fmpz_mod_mpoly_q_print_pretty(t, NULL, ctx);flint_printf("\n");
    fmpz_mod_mpoly_q_add_fmpq(z,y,xq,ctx);fmpz_mod_mpoly_q_print_pretty(z, NULL, ctx);flint_printf("\n");    
    fmpz_mod_mpoly_q_print_pretty(x, NULL, ctx);flint_printf("\n");
    fmpz_mod_mpoly_q_print_pretty(y, NULL, ctx);flint_printf("\n");

    fmpz_clear(m);
    fmpz_mod_mpoly_q_clear(x, ctx);
    fmpz_mod_mpoly_q_clear(y, ctx);
    fmpz_mod_mpoly_q_clear(z, ctx);
    fmpz_mod_mpoly_q_clear(t, ctx);
    fmpz_mod_mpoly_ctx_clear(ctx);



    for (iter = 0; iter < 100000 * 0.1 * flint_test_multiplier(); iter++)
    {
        fmpz_mod_mpoly_ctx_t ctx;
        fmpz_mod_mpoly_q_t A, B, C, D;
        fmpz_mod_mpoly_t t, u;
        fmpz_t m;

        fmpz_init(m);
        fmpz_randtest_unsigned(m, state, 200);
        fmpz_add_ui(m, m, 20);
        fmpz_nextprime(m, m, 0);
        fmpz_mod_mpoly_ctx_init(ctx, 1 + n_randint(state, 4), ORD_LEX, m);

        fmpz_mod_mpoly_q_init(A, ctx);
        fmpz_mod_mpoly_q_init(B, ctx);
        fmpz_mod_mpoly_q_init(C, ctx);
        fmpz_mod_mpoly_q_init(D, ctx);
        fmpz_mod_mpoly_init(t, ctx);
        fmpz_mod_mpoly_init(u, ctx);

        fmpz_mod_mpoly_q_randtest(A, state, 5, 5, ctx);
        fmpz_mod_mpoly_q_randtest(B, state, 5, 5, ctx);

        fmpz_mod_mpoly_q_add(C, A, B, ctx);

        fmpz_mod_mpoly_mul(t, fmpz_mod_mpoly_q_numref(A), fmpz_mod_mpoly_q_denref(B), ctx);
        fmpz_mod_mpoly_mul(u, fmpz_mod_mpoly_q_numref(B), fmpz_mod_mpoly_q_denref(A), ctx);
        fmpz_mod_mpoly_add(fmpz_mod_mpoly_q_numref(D), t, u, ctx);
        fmpz_mod_mpoly_mul(fmpz_mod_mpoly_q_denref(D), fmpz_mod_mpoly_q_denref(A), fmpz_mod_mpoly_q_denref(B), ctx);

        fmpz_mod_mpoly_q_canonicalise(D, ctx);

        if (!fmpz_mod_mpoly_q_equal(C, D, ctx))
        {
            flint_printf("FAIL\n");
            flint_printf("A = "); fmpz_mod_mpoly_q_print_pretty(A, NULL, ctx); flint_printf("\n\n");
            flint_printf("B = "); fmpz_mod_mpoly_q_print_pretty(B, NULL, ctx); flint_printf("\n\n");
            flint_printf("C = "); fmpz_mod_mpoly_q_print_pretty(C, NULL, ctx); flint_printf("\n\n");
            flint_printf("D = "); fmpz_mod_mpoly_q_print_pretty(D, NULL, ctx); flint_printf("\n\n");
            flint_printf("mod = ");fmpz_print(ctx->ffinfo->n);flint_printf("\n\n");
            flint_abort();
        }

        if (n_randint(state, 2))
        {
            fmpz_mod_mpoly_q_set(C, A, ctx);
            fmpz_mod_mpoly_q_add(C, C, B, ctx);
        }
        else
        {
            fmpz_mod_mpoly_q_set(C, B, ctx);
            fmpz_mod_mpoly_q_add(C, A, C, ctx);
        }

        if (!fmpz_mod_mpoly_q_equal(C, D, ctx))
        {
            flint_printf("FAIL (aliasing)\n");
            flint_printf("A = "); fmpz_mod_mpoly_q_print_pretty(A, NULL, ctx); flint_printf("\n\n");
            flint_printf("B = "); fmpz_mod_mpoly_q_print_pretty(B, NULL, ctx); flint_printf("\n\n");
            flint_printf("C = "); fmpz_mod_mpoly_q_print_pretty(C, NULL, ctx); flint_printf("\n\n");
            flint_printf("D = "); fmpz_mod_mpoly_q_print_pretty(D, NULL, ctx); flint_printf("\n\n");
            flint_printf("mod = ");fmpz_print(ctx->ffinfo->n);flint_printf("\n\n");
            flint_abort();
        }

        fmpz_mod_mpoly_q_clear(A, ctx);
        fmpz_mod_mpoly_q_clear(B, ctx);
        fmpz_mod_mpoly_q_clear(C, ctx);
        fmpz_mod_mpoly_q_clear(D, ctx);
        fmpz_mod_mpoly_clear(t, ctx);
        fmpz_mod_mpoly_clear(u, ctx);
        fmpz_clear(m);

        fmpz_mod_mpoly_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
