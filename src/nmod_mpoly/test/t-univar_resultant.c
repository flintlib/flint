/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "nmod_mpoly.h"

void test_resultant(
    const nmod_mpoly_t fx,
    const nmod_mpoly_t gx,
    const nmod_mpoly_ctx_t ctx)
{
    nmod_mpoly_univar_t F, G;
    nmod_poly_t f, g;
    nmod_mpoly_t R;
    mp_limb_t r;

    nmod_mpoly_univar_init(F, ctx);
    nmod_mpoly_univar_init(G, ctx);
    nmod_poly_init(f, nmod_mpoly_ctx_modulus(ctx));
    nmod_poly_init(g, nmod_mpoly_ctx_modulus(ctx));
    nmod_mpoly_init(R, ctx);

    nmod_mpoly_get_nmod_poly(f, fx, 0, ctx);
    nmod_mpoly_get_nmod_poly(g, gx, 0, ctx);
    r = nmod_poly_resultant(f, g);

    nmod_mpoly_to_univar(F, fx, 0, ctx);
    nmod_mpoly_to_univar(G, gx, 0, ctx);
    nmod_mpoly_univar_resultant(R, F, G, ctx);

    if (!nmod_mpoly_equal_ui(R, r, ctx))
    {
        flint_printf("FAIL: Check resultant against univariate\n");
        flint_printf("fx: ");
        nmod_mpoly_print_pretty(fx, NULL, ctx);
        flint_printf("\n");
        flint_printf("gx: ");
        nmod_mpoly_print_pretty(gx, NULL, ctx);
        flint_printf("\n");
        flint_printf("R: ");
        nmod_mpoly_print_pretty(R, NULL, ctx);
        flint_printf("\n");
        flint_printf("r: %wu\n", r);
        fflush(stdout);
        flint_abort();
    }

    nmod_mpoly_univar_clear(F, ctx);
    nmod_mpoly_univar_clear(G, ctx);
    nmod_poly_clear(f);
    nmod_poly_clear(g);
    nmod_mpoly_clear(R, ctx);
}

TEST_FUNCTION_START(nmod_mpoly_univar_resultant, state)
{
    slong i, j;

    {
        nmod_mpoly_ctx_t ctx;
        nmod_mpoly_t t, s;
        nmod_mpoly_univar_t f;
        const char * vars[] = {"a", "b", "c", "d"};

        nmod_mpoly_ctx_init(ctx, 4, ORD_DEGLEX, 3);
        nmod_mpoly_init(t, ctx);
        nmod_mpoly_init(s, ctx);
        nmod_mpoly_univar_init(f, ctx);

        nmod_mpoly_univar_zero(f, ctx);
        nmod_mpoly_set_str_pretty(t, "a", vars, ctx);
        nmod_mpoly_univar_set_coeff_ui(f, 1, t, ctx);
        nmod_mpoly_set_str_pretty(t, "b", vars, ctx);
        nmod_mpoly_univar_set_coeff_ui(f, 0, t, ctx);
        nmod_mpoly_set_str_pretty(t, "1", vars, ctx);
        nmod_mpoly_univar_discriminant(s, f, ctx);
        if (!nmod_mpoly_equal(s, t, ctx))
        {
            flint_printf("FAIL: check linear discriminant\n");
            fflush(stdout);
            flint_abort();
        }

        nmod_mpoly_univar_zero(f, ctx);
        nmod_mpoly_set_str_pretty(t, "a", vars, ctx);
        nmod_mpoly_univar_set_coeff_ui(f, 2, t, ctx);
        nmod_mpoly_set_str_pretty(t, "b", vars, ctx);
        nmod_mpoly_univar_set_coeff_ui(f, 1, t, ctx);
        nmod_mpoly_set_str_pretty(t, "c", vars, ctx);
        nmod_mpoly_univar_set_coeff_ui(f, 0, t, ctx);
        nmod_mpoly_set_str_pretty(t, "b^2-4*a*c", vars, ctx);
        nmod_mpoly_univar_discriminant(s, f, ctx);
        if (!nmod_mpoly_equal(s, t, ctx))
        {
            flint_printf("FAIL: check quadratic discriminant\n");
            fflush(stdout);
            flint_abort();
        }

        nmod_mpoly_univar_zero(f, ctx);
        nmod_mpoly_set_str_pretty(t, "a", vars, ctx);
        nmod_mpoly_univar_set_coeff_ui(f, 3, t, ctx);
        nmod_mpoly_set_str_pretty(t, "b", vars, ctx);
        nmod_mpoly_univar_set_coeff_ui(f, 2, t, ctx);
        nmod_mpoly_set_str_pretty(t, "c", vars, ctx);
        nmod_mpoly_univar_set_coeff_ui(f, 1, t, ctx);
        nmod_mpoly_set_str_pretty(t, "d", vars, ctx);
        nmod_mpoly_univar_set_coeff_ui(f, 0, t, ctx);
        nmod_mpoly_set_str_pretty(t, "b^2*c^2-4*a*c^3-4*b^3*d-27*a^2*d^2+18*a*b*c*d", vars, ctx);
        nmod_mpoly_univar_discriminant(s, f, ctx);
        if (!nmod_mpoly_equal(s, t, ctx))
        {
            flint_printf("FAIL: check cubic discriminant\n");
            fflush(stdout);
            flint_abort();
        }

        nmod_mpoly_clear(t, ctx);
        nmod_mpoly_clear(s, ctx);
        nmod_mpoly_univar_clear(f, ctx);
        nmod_mpoly_ctx_clear(ctx);
    }

    for (i = 0; i < 10 * flint_test_multiplier(); i++)
    {
        nmod_mpoly_ctx_t ctx;
        nmod_mpoly_t f, g, t;
        mp_limb_t modulus;

        modulus = n_randint(state, FLINT_BITS - 1) + 1;
        modulus = n_randbits(state, modulus);
        modulus = n_nextprime(modulus, 1);
        nmod_mpoly_ctx_init(ctx, 1, mpoly_ordering_randtest(state), modulus);
        nmod_mpoly_init(f, ctx);
        nmod_mpoly_init(g, ctx);
        nmod_mpoly_init(t, ctx);

        for (j = 0; j < 10; j++)
        {
            nmod_mpoly_randtest_bound(f, state, 3, 3, ctx);
            if (nmod_mpoly_is_zero(f, ctx))
                nmod_mpoly_one(f, ctx);

            nmod_mpoly_zero(g, ctx);

            while (nmod_mpoly_degree_si(f, 0, ctx) < 100)
            {
                nmod_mpoly_randtest_bound(t, state, 5, 10, ctx);
                nmod_mpoly_mul(t, t, f, ctx);
                nmod_mpoly_add(g, g, t, ctx);
                nmod_mpoly_swap(f, g, ctx);
            }

            nmod_mpoly_randtest_bound(f, state, 20, 50, ctx);
            nmod_mpoly_randtest_bound(g, state, 20, 50, ctx);
            test_resultant(f, g, ctx);
        }

        nmod_mpoly_clear(f, ctx);
        nmod_mpoly_clear(g, ctx);
        nmod_mpoly_clear(t, ctx);
        nmod_mpoly_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
