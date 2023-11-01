/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fq_nmod_mpoly.h"

/* fq_nmod_poly_resultant is missing */
#if 0
void test_resultant(
    const fq_nmod_mpoly_t fx,
    const fq_nmod_mpoly_t gx,
    const fq_nmod_mpoly_ctx_t ctx)
{
    fq_nmod_mpoly_univar_t F, G;
    fq_nmod_poly_t f, g;
    fq_nmod_mpoly_t R;
    fq_nmod_t r;

    fq_nmod_mpoly_univar_init(F, ctx);
    fq_nmod_mpoly_univar_init(G, ctx);
    fq_nmod_poly_init(f, ctx->fqctx);
    fq_nmod_poly_init(g, ctx->fqctx);
    fq_nmod_mpoly_init(R, ctx);
    fq_nmod_init(r, ctx->fqctx);

    fq_nmod_mpoly_get_fq_nmod_poly(f, fx, 0, ctx);
    fq_nmod_mpoly_get_fq_nmod_poly(g, gx, 0, ctx);
    fq_nmod_poly_resultant(r, f, g, ctx->fqctx);

    fq_nmod_mpoly_to_univar(F, fx, 0, ctx);
    fq_nmod_mpoly_to_univar(G, gx, 0, ctx);
    fq_nmod_mpoly_univar_resultant(R, F, G, ctx);

    if (!fq_nmod_mpoly_equal_fq_nmod(R, r, ctx))
    {
        flint_printf("FAIL: Check resultant against univariate\n");
        flint_printf("fx: ");
        fq_nmod_mpoly_print_pretty(fx, NULL, ctx);
        flint_printf("\n");
        flint_printf("gx: ");
        fq_nmod_mpoly_print_pretty(gx, NULL, ctx);
        flint_printf("\n");
        flint_printf("R: ");
        fq_nmod_mpoly_print_pretty(R, NULL, ctx);
        flint_printf("\n");
        flint_printf("r: %wu\n", r);
        fflush(stdout);
        flint_abort();
    }

    fq_nmod_clear(r, ctx->fqctx);
    fq_nmod_mpoly_univar_clear(F, ctx);
    fq_nmod_mpoly_univar_clear(G, ctx);
    fq_nmod_poly_clear(f, ctx->fqctx);
    fq_nmod_poly_clear(g, ctx->fqctx);
    fq_nmod_mpoly_clear(R, ctx);
}
#endif

TEST_FUNCTION_START(fq_nmod_mpoly_univar_resultant, state)
{
    {
        fq_nmod_mpoly_ctx_t ctx;
        fq_nmod_mpoly_t t, s;
        fq_nmod_mpoly_univar_t f;
        const char * vars[] = {"a", "b", "c", "d"};

        fq_nmod_mpoly_ctx_init_deg(ctx, 4, ORD_DEGLEX, 3, 4);
        fq_nmod_mpoly_init(t, ctx);
        fq_nmod_mpoly_init(s, ctx);
        fq_nmod_mpoly_univar_init(f, ctx);

        fq_nmod_mpoly_univar_zero(f, ctx);
        fq_nmod_mpoly_set_str_pretty(t, "a", vars, ctx);
        fq_nmod_mpoly_univar_set_coeff_ui(f, 1, t, ctx);
        fq_nmod_mpoly_set_str_pretty(t, "b", vars, ctx);
        fq_nmod_mpoly_univar_set_coeff_ui(f, 0, t, ctx);
        fq_nmod_mpoly_set_str_pretty(t, "1", vars, ctx);
        fq_nmod_mpoly_univar_discriminant(s, f, ctx);
        if (!fq_nmod_mpoly_equal(s, t, ctx))
        {
            flint_printf("FAIL: check linear discriminant\n");
            fflush(stdout);
            flint_abort();
        }

        fq_nmod_mpoly_univar_zero(f, ctx);
        fq_nmod_mpoly_set_str_pretty(t, "a", vars, ctx);
        fq_nmod_mpoly_univar_set_coeff_ui(f, 2, t, ctx);
        fq_nmod_mpoly_set_str_pretty(t, "b", vars, ctx);
        fq_nmod_mpoly_univar_set_coeff_ui(f, 1, t, ctx);
        fq_nmod_mpoly_set_str_pretty(t, "c", vars, ctx);
        fq_nmod_mpoly_univar_set_coeff_ui(f, 0, t, ctx);
        fq_nmod_mpoly_set_str_pretty(t, "b^2-4*a*c", vars, ctx);
        fq_nmod_mpoly_univar_discriminant(s, f, ctx);
        if (!fq_nmod_mpoly_equal(s, t, ctx))
        {
            flint_printf("FAIL: check quadratic discriminant\n");
            fflush(stdout);
            flint_abort();
        }

        fq_nmod_mpoly_univar_zero(f, ctx);
        fq_nmod_mpoly_set_str_pretty(t, "a", vars, ctx);
        fq_nmod_mpoly_univar_set_coeff_ui(f, 3, t, ctx);
        fq_nmod_mpoly_set_str_pretty(t, "b", vars, ctx);
        fq_nmod_mpoly_univar_set_coeff_ui(f, 2, t, ctx);
        fq_nmod_mpoly_set_str_pretty(t, "c", vars, ctx);
        fq_nmod_mpoly_univar_set_coeff_ui(f, 1, t, ctx);
        fq_nmod_mpoly_set_str_pretty(t, "d", vars, ctx);
        fq_nmod_mpoly_univar_set_coeff_ui(f, 0, t, ctx);
        fq_nmod_mpoly_set_str_pretty(t, "b^2*c^2-4*a*c^3-4*b^3*d-27*a^2*d^2+18*a*b*c*d", vars, ctx);
        fq_nmod_mpoly_univar_discriminant(s, f, ctx);
        if (!fq_nmod_mpoly_equal(s, t, ctx))
        {
            flint_printf("FAIL: check cubic discriminant\n");
            fflush(stdout);
            flint_abort();
        }

        fq_nmod_mpoly_clear(t, ctx);
        fq_nmod_mpoly_clear(s, ctx);
        fq_nmod_mpoly_univar_clear(f, ctx);
        fq_nmod_mpoly_ctx_clear(ctx);
    }

#if 0
    for (i = 0; i < 10 * flint_test_multiplier(); i++)
    {
        fq_nmod_mpoly_ctx_t ctx;
        fq_nmod_mpoly_t f, g, t;
        mp_limb_t modulus;

        modulus = n_randint(state, FLINT_BITS - 1) + 1;
        modulus = n_randbits(state, modulus);
        modulus = n_nextprime(modulus, 1);
        fq_nmod_mpoly_ctx_init_rand(ctx, state, 1, modulus);
        fq_nmod_mpoly_init(f, ctx);
        fq_nmod_mpoly_init(g, ctx);
        fq_nmod_mpoly_init(t, ctx);

        for (j = 0; j < 10; j++)
        {
            fq_nmod_mpoly_randtest_bound(f, state, 3, 3, ctx);
            if (fq_nmod_mpoly_is_zero(f, ctx))
                fq_nmod_mpoly_one(f, ctx);

            fq_nmod_mpoly_zero(g, ctx);

            while (fq_nmod_mpoly_degree_si(f, 0, ctx) < 100)
            {
                fq_nmod_mpoly_randtest_bound(t, state, 5, 10, ctx);
                fq_nmod_mpoly_mul(t, t, f, ctx);
                fq_nmod_mpoly_add(g, g, t, ctx);
                fq_nmod_mpoly_swap(f, g, ctx);
            }

            fq_nmod_mpoly_randtest_bound(f, state, 20, 50, ctx);
            fq_nmod_mpoly_randtest_bound(g, state, 20, 50, ctx);
            test_resultant(f, g, ctx);
        }

        fq_nmod_mpoly_clear(f, ctx);
        fq_nmod_mpoly_clear(g, ctx);
        fq_nmod_mpoly_clear(t, ctx);
        fq_nmod_mpoly_ctx_clear(ctx);
    }
#endif

    TEST_FUNCTION_END(state);
}
