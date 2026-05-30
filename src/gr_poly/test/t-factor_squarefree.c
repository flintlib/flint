/*
    Copyright (C) 2023, 2025, 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
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

    for (iter = 0; iter < 100 * flint_test_multiplier(); iter++)
    {
        gr_ctx_t ctx;
        gr_poly_t A, B, P, Q, G1, G2, G3;
        gr_ptr c;
        gr_vec_t fac, exp;
        gr_ctx_t poly_ctx, fmpz_ctx;
        ulong e, emax, maxexp;
        fmpz *expc;
        slong i, j, exp_bound, factor_len_bound, factor_bound;
        int status = GR_SUCCESS;
        int must_succeed;
        slong ngens;

        if (n_randint(state, 2))
        {
            gr_ctx_init_random_finite_field(ctx, state);
        }
        else
        {
            gr_ctx_init_random_commutative_ring(ctx, state);

            while (gr_ctx_is_integral_domain(ctx) != T_TRUE || ctx->methods == _ca_methods)
            {
                gr_ctx_clear(ctx);
                gr_ctx_init_random(ctx, state);
            }
        }

        must_succeed = (ctx->which_ring == GR_CTX_FMPQ  ||
                        ctx->which_ring == GR_CTX_FMPZ  ||
                        ctx->which_ring == GR_CTX_FMPZI ||
                        ctx->which_ring == GR_CTX_FMPZ_POLY ||
                        (ctx->which_ring == GR_CTX_NMOD && gr_ctx_is_field(ctx) == T_TRUE) ||
                        (ctx->which_ring == GR_CTX_NMOD8 && gr_ctx_is_field(ctx) == T_TRUE) ||
                        (ctx->which_ring == GR_CTX_MPN_MOD && gr_ctx_is_field(ctx) == T_TRUE) ||
                        (ctx->which_ring == GR_CTX_FMPZ_MOD && gr_ctx_is_field(ctx) == T_TRUE) ||
                        ctx->which_ring == GR_CTX_FQ_NMOD ||
                        ctx->which_ring == GR_CTX_FQ_ZECH);

        gr_ctx_init_gr_poly(poly_ctx, ctx);
        gr_ctx_init_fmpz(fmpz_ctx);

        gr_poly_init(A, ctx);
        gr_poly_init(B, ctx);
        gr_poly_init(G1, ctx);
        gr_poly_init(G2, ctx);
        gr_poly_init(G3, ctx);
        gr_poly_init(P, ctx);
        gr_poly_init(Q, ctx);

        gr_vec_init(fac, 0, poly_ctx);
        gr_vec_init(exp, 0, fmpz_ctx);

        c = gr_heap_init(ctx);

        if (gr_ctx_is_finite(ctx) == T_TRUE)
        {
            factor_bound = 5;
            factor_len_bound = 5;
            exp_bound = 5;
        }
        else
        {
            factor_bound = 3;
            factor_len_bound = 3;
            exp_bound = 3;

            /* Keep tests fast with multivariates */
            if (gr_ctx_ngens(&ngens, ctx) == GR_SUCCESS)
            {
                if (ngens > 1)
                {
                    factor_bound = 2;
                    factor_len_bound = 2;
                    exp_bound = 2;
                }
            }
        }

        status |= gr_poly_one(P, ctx);
        emax = 0;

        for (i = 0; i < factor_bound; i++)
        {
            do {
                status |= gr_poly_randtest(A, state, factor_len_bound, ctx);
            } while (A->length < 2);

            e = 1 + n_randint(state, exp_bound);

            if (gr_poly_gcd(G1, P, A, ctx) == GR_SUCCESS && gr_poly_is_scalar(G1, ctx) == T_TRUE)
            {
                status |= gr_poly_pow_ui(B, A, e, ctx);
                status |= gr_poly_mul(P, P, B, ctx);
                emax = FLINT_MAX(emax, e);
            }
        }

        if (status != GR_SUCCESS && !must_succeed)
            goto cleanup;

        status = gr_poly_factor_squarefree(c, fac, exp, P, ctx);

        if (must_succeed && status != GR_SUCCESS)
        {
            flint_printf("FAIL (unexpected GR_UNABLE)\n\n");
            gr_ctx_println(ctx);
            flint_printf("P = "); gr_poly_print(P, ctx); flint_printf("\n");
            flint_abort();
        }

        if (status == GR_SUCCESS)
        {
            int status1 = GR_SUCCESS;

            for (i = 0; i < fac->length; i++)
            {
                status1 |= gr_poly_derivative(Q, gr_vec_entry_ptr(fac, i, poly_ctx), ctx);
                status1 |= gr_poly_gcd(Q, Q, gr_vec_entry_ptr(fac, i, poly_ctx), ctx);

                if (status1 == GR_SUCCESS && gr_poly_is_scalar(Q, ctx) == T_FALSE)
                {
                    flint_printf("FAIL (factor not squarefree)\n\n");
                    gr_ctx_println(ctx);
                    flint_printf("emax = %wd\n\n", emax);
                    flint_printf("fac = %{gr*}\n\n", fac->entries, fac->length, poly_ctx);
                    flint_printf("exp = %{gr*}\n\n", exp->entries, exp->length, fmpz_ctx);
                    flint_printf("c = %{gr}\n\n", c, ctx);
                    flint_printf("P = "); gr_poly_print(P, ctx); flint_printf("\n");
                    flint_printf("Q = "); gr_poly_print(Q, ctx); flint_printf("\n");
                    flint_abort();
                }
            }
        }

        if (status == GR_SUCCESS)
        {
            expc = exp->entries;

            status |= gr_poly_one(Q, ctx);
            for (i = 0; i < fac->length; i++)
                for (j = 0; j < expc[i]; j++)
                    status |= gr_poly_mul(Q, Q,
                        gr_vec_entry_ptr(fac, i, poly_ctx), ctx);
            status |= gr_poly_mul_scalar(Q, Q, c, ctx);

            maxexp = 0;
            for (i = 0; i < fac->length; i++)
                maxexp = FLINT_MAX(maxexp, expc[i]);

            if (status == GR_SUCCESS &&
                (gr_poly_equal(P, Q, ctx) == T_FALSE ||
                 maxexp < emax))
            {
                flint_printf("FAIL (product)\n\n");
                gr_ctx_println(ctx);
                flint_printf("emax = %wd\n\n", emax);
                flint_printf("fac = %{gr*}\n\n", fac->entries, fac->length, poly_ctx);
                flint_printf("exp = %{gr*}\n\n", exp->entries, exp->length, fmpz_ctx);
                flint_printf("c = %{gr}\n\n", c, ctx);
                flint_printf("P = "); gr_poly_print(P, ctx); flint_printf("\n");
                flint_printf("Q = "); gr_poly_print(Q, ctx); flint_printf("\n");
                flint_abort();
            }
        }

        gr_poly_clear(A, ctx);
        gr_poly_clear(B, ctx);
        gr_poly_clear(G1, ctx);
        gr_poly_clear(G2, ctx);
        gr_poly_clear(G3, ctx);
        gr_poly_clear(P, ctx);
        gr_poly_clear(Q, ctx);

cleanup:
        gr_vec_clear(fac, poly_ctx);
        gr_vec_clear(exp, fmpz_ctx);
        gr_heap_clear(c, ctx);
        gr_ctx_clear(poly_ctx);
        gr_ctx_clear(fmpz_ctx);
        gr_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
