/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include <stdlib.h>
#include "fmpz_mod_mpoly.h"

void test_resultant(
    const fmpz_mod_mpoly_t fx,
    const fmpz_mod_mpoly_t gx,
    const fmpz_mod_mpoly_ctx_t ctx)
{
    fmpz_mod_mpoly_univar_t F, G;
    fmpz_mod_poly_t f, g;
    fmpz_mod_mpoly_t R;
    fmpz_t r;

    fmpz_mod_mpoly_univar_init(F, ctx);
    fmpz_mod_mpoly_univar_init(G, ctx);
    fmpz_mod_poly_init(f, ctx->ffinfo);
    fmpz_mod_poly_init(g, ctx->ffinfo);
    fmpz_mod_mpoly_init(R, ctx);
    fmpz_init(r);

    fmpz_mod_mpoly_get_fmpz_mod_poly(f, fx, 0, ctx);
    fmpz_mod_mpoly_get_fmpz_mod_poly(g, gx, 0, ctx);
    fmpz_mod_poly_resultant(r, f, g, ctx->ffinfo);

    fmpz_mod_mpoly_to_univar(F, fx, 0, ctx);
    fmpz_mod_mpoly_to_univar(G, gx, 0, ctx);
    fmpz_mod_mpoly_univar_resultant(R, F, G, ctx);

    if (!fmpz_mod_mpoly_equal_fmpz(R, r, ctx))
    {
        flint_printf("FAIL: Check resultant against univariate\n");
        flint_printf("fx: ");
        fmpz_mod_mpoly_print_pretty(fx, NULL, ctx);
        flint_printf("\n");
        flint_printf("gx: ");
        fmpz_mod_mpoly_print_pretty(gx, NULL, ctx);
        flint_printf("\n");
        flint_printf("R: ");
        fmpz_mod_mpoly_print_pretty(R, NULL, ctx);
        flint_printf("\n");
        flint_printf("r: ");
        fmpz_print(r);
        flint_printf("\n");
        fflush(stdout);
        flint_abort();
    }

    fmpz_mod_mpoly_univar_clear(F, ctx);
    fmpz_mod_mpoly_univar_clear(G, ctx);
    fmpz_mod_poly_clear(f, ctx->ffinfo);
    fmpz_mod_poly_clear(g, ctx->ffinfo);
    fmpz_mod_mpoly_clear(R, ctx);
    fmpz_clear(r);
}

int
main(void)
{
    slong i, j;
    FLINT_TEST_INIT(state);

    flint_printf("univar_resultant....");
    fflush(stdout);

    {
        fmpz p = 1000003;
        fmpz_mod_mpoly_ctx_t ctx;
        fmpz_mod_mpoly_t t, s;
        fmpz_mod_mpoly_univar_t f;
        const char * vars[] = {"a", "b", "c", "d"};

        fmpz_mod_mpoly_ctx_init(ctx, 4, ORD_DEGLEX, &p);
        fmpz_mod_mpoly_init(t, ctx);
        fmpz_mod_mpoly_init(s, ctx);
        fmpz_mod_mpoly_univar_init(f, ctx);

        fmpz_mod_mpoly_univar_zero(f, ctx);
        fmpz_mod_mpoly_set_str_pretty(t, "a", vars, ctx);
        fmpz_mod_mpoly_univar_set_coeff_ui(f, 1, t, ctx);
        fmpz_mod_mpoly_set_str_pretty(t, "b", vars, ctx);
        fmpz_mod_mpoly_univar_set_coeff_ui(f, 0, t, ctx);
        fmpz_mod_mpoly_set_str_pretty(t, "1", vars, ctx);
        fmpz_mod_mpoly_univar_discriminant(s, f, ctx);
        if (!fmpz_mod_mpoly_equal(s, t, ctx))
        {
            flint_printf("FAIL: check linear discriminant\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_mod_mpoly_univar_zero(f, ctx);
        fmpz_mod_mpoly_set_str_pretty(t, "a", vars, ctx);
        fmpz_mod_mpoly_univar_set_coeff_ui(f, 2, t, ctx);
        fmpz_mod_mpoly_set_str_pretty(t, "b", vars, ctx);
        fmpz_mod_mpoly_univar_set_coeff_ui(f, 1, t, ctx);
        fmpz_mod_mpoly_set_str_pretty(t, "c", vars, ctx);
        fmpz_mod_mpoly_univar_set_coeff_ui(f, 0, t, ctx);
        fmpz_mod_mpoly_set_str_pretty(t, "b^2-4*a*c", vars, ctx);
        fmpz_mod_mpoly_univar_discriminant(s, f, ctx);
        if (!fmpz_mod_mpoly_equal(s, t, ctx))
        {
            flint_printf("FAIL: check quadratic discriminant\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_mod_mpoly_univar_zero(f, ctx);
        fmpz_mod_mpoly_set_str_pretty(t, "a", vars, ctx);
        fmpz_mod_mpoly_univar_set_coeff_ui(f, 3, t, ctx);
        fmpz_mod_mpoly_set_str_pretty(t, "b", vars, ctx);
        fmpz_mod_mpoly_univar_set_coeff_ui(f, 2, t, ctx);
        fmpz_mod_mpoly_set_str_pretty(t, "c", vars, ctx);
        fmpz_mod_mpoly_univar_set_coeff_ui(f, 1, t, ctx);
        fmpz_mod_mpoly_set_str_pretty(t, "d", vars, ctx);
        fmpz_mod_mpoly_univar_set_coeff_ui(f, 0, t, ctx);
        fmpz_mod_mpoly_set_str_pretty(t, "b^2*c^2-4*a*c^3-4*b^3*d-27*a^2*d^2+18*a*b*c*d", vars, ctx);
        fmpz_mod_mpoly_univar_discriminant(s, f, ctx);
        if (!fmpz_mod_mpoly_equal(s, t, ctx))
        {
            flint_printf("FAIL: check cubic discriminant\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_mod_mpoly_clear(t, ctx);
        fmpz_mod_mpoly_clear(s, ctx);
        fmpz_mod_mpoly_univar_clear(f, ctx);
        fmpz_mod_mpoly_ctx_clear(ctx);
    }

    for (i = 0; i < 10 * flint_test_multiplier(); i++)
    {
        fmpz_mod_mpoly_ctx_t ctx;
        fmpz_mod_mpoly_t f, g, t;

        fmpz_mod_mpoly_ctx_init_rand_bits_prime(ctx, state, 1, 200);
        if (fmpz_mod_mpoly_ctx_nvars(ctx) < 1)
        {
            fmpz_mod_mpoly_ctx_clear(ctx);
            continue;
        }

        fmpz_mod_mpoly_init(f, ctx);
        fmpz_mod_mpoly_init(g, ctx);
        fmpz_mod_mpoly_init(t, ctx);

        for (j = 0; j < 10; j++)
        {
            fmpz_mod_mpoly_randtest_bound(f, state, 3, 3, ctx);
            if (fmpz_mod_mpoly_is_zero(f, ctx))
                fmpz_mod_mpoly_one(f, ctx);

            fmpz_mod_mpoly_zero(g, ctx);

            while (fmpz_mod_mpoly_degree_si(f, 0, ctx) < 100)
            {
                fmpz_mod_mpoly_randtest_bound(t, state, 5, 10, ctx);
                fmpz_mod_mpoly_mul(t, t, f, ctx);
                fmpz_mod_mpoly_add(g, g, t, ctx);
                fmpz_mod_mpoly_swap(f, g, ctx);
            }

            fmpz_mod_mpoly_randtest_bound(f, state, 20, 50, ctx);
            fmpz_mod_mpoly_randtest_bound(g, state, 20, 50, ctx);
            test_resultant(f, g, ctx);
        }

        fmpz_mod_mpoly_clear(f, ctx);
        fmpz_mod_mpoly_clear(g, ctx);
        fmpz_mod_mpoly_clear(t, ctx);
        fmpz_mod_mpoly_ctx_clear(ctx);  
    }

    FLINT_TEST_CLEANUP(state);

    flint_printf("PASS\n");
    return 0;
}

