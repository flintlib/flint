/*
    Copyright (C) 2023 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "ulong_extras.h"
#include "gr_vec.h"
#include "gr_poly.h"

FLINT_DLL extern gr_static_method_table _ca_methods;

TEST_FUNCTION_START(gr_poly_factor_squarefree, state)
{
    slong iter;

    for (iter = 0; iter < 1000; iter++)
    {
        gr_ctx_t ctx;
        gr_poly_t A, B, C, P, Q, G1, G2, G3;
        gr_ptr c;
        gr_vec_t F;
        gr_ctx_t poly_ctx, fmpz_ctx;
        ulong ea, eb, ec, maxexp;
        fmpz * expc;
        gr_vec_t exp;
        slong i, j;
        int status = GR_SUCCESS;
        int attempts = 0;

        gr_ctx_init_random(ctx, state);

        while (gr_ctx_is_field(ctx) != T_TRUE || ctx->methods == _ca_methods)
        {
            gr_ctx_clear(ctx);
            gr_ctx_init_random(ctx, state);
        }

        gr_ctx_init_gr_poly(poly_ctx, ctx);
        gr_ctx_init_fmpz(fmpz_ctx);

        gr_poly_init(A, ctx);
        gr_poly_init(B, ctx);
        gr_poly_init(C, ctx);
        gr_poly_init(G1, ctx);
        gr_poly_init(G2, ctx);
        gr_poly_init(G3, ctx);
        gr_poly_init(P, ctx);
        gr_poly_init(Q, ctx);

        gr_vec_init(F, 0, poly_ctx);
        gr_vec_init(exp, 0, fmpz_ctx);

        c = gr_heap_init(ctx);

        attempts = 0;

        do
        {
            attempts++;

            if (ctx->methods == _ca_methods)
            {
                ea = 1 + n_randint(state, 2);
                eb = 1 + n_randint(state, 2);
                ec = 1 + n_randint(state, 2);

                status |= gr_poly_randtest(A, state, 2, ctx);
                status |= gr_poly_randtest(B, state, 2, ctx);
                status |= gr_poly_randtest(C, state, 2, ctx);
            }
            else
            {
                ea = 1 + n_randint(state, 2);
                eb = 1 + n_randint(state, 2);
                ec = 1 + n_randint(state, 4);

                status |= gr_poly_randtest(A, state, 3, ctx);
                status |= gr_poly_randtest(B, state, 3, ctx);
                status |= gr_poly_randtest(C, state, 3, ctx);
            }

            if (ctx->which_ring != GR_CTX_FMPQ && attempts > 4)
            {
                status = GR_UNABLE;
                break;
            }
        }
        while (A->length < 2 || B->length < 2 || C->length < 2 ||
               gr_poly_gcd(G1, A, B, ctx) != GR_SUCCESS || (gr_poly_is_one(G1, ctx) != T_TRUE) ||
               gr_poly_gcd(G2, A, C, ctx) != GR_SUCCESS || (gr_poly_is_one(G2, ctx) != T_TRUE) ||
               gr_poly_gcd(G3, B, C, ctx) != GR_SUCCESS || (gr_poly_is_one(G3, ctx) != T_TRUE));

        status |= gr_poly_one(P, ctx);
        for (i = 0; i < ea; i++)
            status |= gr_poly_mul(P, P, A, ctx);
        for (i = 0; i < eb; i++)
            status |= gr_poly_mul(P, P, B, ctx);
        for (i = 0; i < ec; i++)
            status |= gr_poly_mul(P, P, C, ctx);

        status |= gr_poly_factor_squarefree(c, F, exp, P, ctx);

        if (ctx->which_ring == GR_CTX_FMPQ && status != GR_SUCCESS)
        {
            flint_printf("FAIL (unexpected failure)\n\n");
            flint_printf("P = "); gr_poly_print(P, ctx); flint_printf("\n");
            flint_abort();
        }

        if (status == GR_SUCCESS)
        {
            expc = exp->entries;

/*
            printf("LENGTHS %ld %ld\n", F->length, exp->length);
            gr_vec_print(F, poly_ctx); printf("\n\n");
            gr_vec_print(exp, fmpz_ctx); printf("\n\n");
*/

            status |= gr_poly_one(Q, ctx);
            for (i = 0; i < F->length; i++)
                for (j = 0; j < expc[i]; j++)
                    status |= gr_poly_mul(Q, Q, gr_vec_entry_ptr(F, i, poly_ctx), ctx);

            status |= gr_poly_mul_scalar(Q, Q, c, ctx);

            maxexp = 0;
            for (i = 0; i < F->length; i++)
                maxexp = FLINT_MAX(maxexp, expc[i]);

            if (status == GR_SUCCESS && (gr_poly_equal(P, Q, ctx) == T_FALSE || maxexp < FLINT_MAX(FLINT_MAX(ea, eb), ec)))
            {
                flint_printf("FAIL (product)\n\n");
                flint_printf("A = "); gr_poly_print(A, ctx); flint_printf("\n");
                flint_printf("ea = %wu\n\n", ea);
                flint_printf("B = "); gr_poly_print(B, ctx); flint_printf("\n");
                flint_printf("eb = %wu\n\n", eb);
                flint_printf("C = "); gr_poly_print(C, ctx); flint_printf("\n");
                flint_printf("ec = %wu\n\n", ec);
                flint_printf("P = "); gr_poly_print(P, ctx); flint_printf("\n");
                flint_printf("Q = "); gr_poly_print(Q, ctx); flint_printf("\n");

                for (i = 0; i < F->length; i++)
                {
                    flint_printf("Multiplicity %wu: ", exp[i]);
                    gr_poly_print(gr_vec_entry_ptr(F, i, poly_ctx), ctx);
                }

                flint_abort();
            }
        }

        gr_poly_clear(A, ctx);
        gr_poly_clear(B, ctx);
        gr_poly_clear(C, ctx);
        gr_poly_clear(G1, ctx);
        gr_poly_clear(G2, ctx);
        gr_poly_clear(G3, ctx);
        gr_poly_clear(P, ctx);
        gr_poly_clear(Q, ctx);

        gr_vec_clear(F, poly_ctx);
        gr_vec_clear(exp, fmpz_ctx);

        gr_heap_clear(c, ctx);

        gr_ctx_clear(poly_ctx);
        gr_ctx_clear(fmpz_ctx);
        gr_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
