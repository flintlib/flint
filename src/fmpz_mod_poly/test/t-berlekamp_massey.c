/*
    Copyright (C) 2019 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "ulong_extras.h"
#include "fmpz.h"
#include "fmpz_mod.h"
#include "fmpz_mod_poly.h"

TEST_FUNCTION_START(fmpz_mod_poly_berlekamp_massey, state)
{
    slong i, j, k, l;
    fmpz_mod_ctx_t ctx;

    fmpz_mod_ctx_init_ui(ctx, 101);

    {
        fmpz_mod_berlekamp_massey_t B;
        fmpz_mod_berlekamp_massey_init(B, ctx);
        fmpz_mod_berlekamp_massey_add_point_ui(B, 1, ctx);
        fmpz_mod_berlekamp_massey_add_point_ui(B, 1, ctx);
        fmpz_mod_berlekamp_massey_add_point_ui(B, 2, ctx);
        fmpz_mod_berlekamp_massey_add_point_ui(B, 3, ctx);
        fmpz_mod_berlekamp_massey_add_point_ui(B, 5, ctx);
        fmpz_mod_berlekamp_massey_reduce(B, ctx);
        if (2 != fmpz_mod_poly_degree(fmpz_mod_berlekamp_massey_V_poly(B), ctx))
        {
            printf("FAIL\n");
            flint_printf("check fibonacci\n");
            fflush(stdout);
            flint_abort();
        }
        fmpz_mod_berlekamp_massey_clear(B, ctx);
    }

    for (i = 0; i < 15 * flint_test_multiplier(); i++)
    {
        fmpz_mod_berlekamp_massey_t B1, B2;
        fmpz_mod_berlekamp_massey_init(B1, ctx);
        fmpz_mod_berlekamp_massey_init(B2, ctx);
        for (j = 0; j < 10; j++)
        {
            fmpz_t p, apoint;

            fmpz_init(apoint);
            fmpz_init_set_ui(p, n_randtest_prime(state, 1));
            fmpz_mod_ctx_set_modulus(ctx, p);

            fmpz_mod_berlekamp_massey_start_over(B1, ctx);
            fmpz_mod_berlekamp_massey_start_over(B2, ctx);

            /* check intermediate reductions match */
            for (k = 0; k < 10; k++)
            {
                fmpz_randm(apoint, state, p);
                fmpz_mod_berlekamp_massey_add_point(B1, apoint, ctx);
                fmpz_mod_berlekamp_massey_add_zeros(B1, n_randint(state, 5), ctx);
                fmpz_randm(apoint, state, p);
                fmpz_mod_berlekamp_massey_add_point(B1, apoint, ctx);
                if (n_randlimb(state) & 1)
                    fmpz_mod_berlekamp_massey_reduce(B1, ctx);
            }
            fmpz_mod_berlekamp_massey_add_points(B2,
                                    fmpz_mod_berlekamp_massey_points(B1),
                                    fmpz_mod_berlekamp_massey_point_count(B1), ctx);
            fmpz_mod_berlekamp_massey_reduce(B2, ctx);
            fmpz_mod_berlekamp_massey_reduce(B1, ctx);
            if (!fmpz_mod_poly_equal(fmpz_mod_berlekamp_massey_V_poly(B1),
                                     fmpz_mod_berlekamp_massey_V_poly(B2), ctx))
            {
                printf("FAIL\n");
                flint_printf("check intermediate reductions match\n"
                                                   "i = %wd, j = %wd\n", i, j);
                fflush(stdout);
                flint_abort();
            }

            /*
                Check berlekamp-massey does its job - 2k coefficients of

                    u     a1    a2          a(2k)
                   --- = --- + --- + ... + ------ + ...
                    v     x    x^2         x^(2k)

                should be sufficient to reconstruct a divisor of v
            */
            for (k = 0; k < 15; k++)
            {
                fmpz_mod_poly_t u, v, s, q, r;
                fmpz_mod_poly_init(u, ctx);
                fmpz_mod_poly_init(v, ctx);
                fmpz_mod_poly_init(s, ctx);
                fmpz_mod_poly_init(q, ctx);
                fmpz_mod_poly_init(r, ctx);

                /* deg(u) < deg(v), deg(v) = k */
                fmpz_mod_poly_randtest(u, state, k, ctx);
                fmpz_mod_poly_randtest_monic(v, state, k + 1, ctx);

                /* q has enough coefficients of expansion of u/v at infty */
                fmpz_mod_poly_shift_left(s, u, 2*k, ctx);
                fmpz_mod_poly_divrem(q, r, s, v, ctx);

                fmpz_mod_berlekamp_massey_start_over(B1, ctx);
                for (l = 2*k - 1; l >= 0 ; l--)
                {
                    fmpz_mod_poly_get_coeff_fmpz(apoint, q, l, ctx);
                    fmpz_mod_berlekamp_massey_add_point(B1, apoint, ctx);
                }
                fmpz_mod_berlekamp_massey_reduce(B1, ctx);
                fmpz_mod_poly_divrem(q, r, v, fmpz_mod_berlekamp_massey_V_poly(B1), ctx);
                if (!fmpz_mod_poly_is_zero(r, ctx))
                {
                    flint_printf("check berlekamp_massey does its job\n"
                                       "i = %wd, j = %wd, k = %wd\n", i, j, k);
                    printf("v: "); fmpz_mod_poly_print_pretty(v, "#", ctx); printf("\n");
                    printf("B: "); fmpz_mod_berlekamp_massey_print(B1, ctx); printf("\n");
                    fflush(stdout);
                    flint_abort();
                }
                if (fmpz_mod_poly_degree(fmpz_mod_berlekamp_massey_R_poly(B1), ctx) >=
                    fmpz_mod_poly_degree(fmpz_mod_berlekamp_massey_V_poly(B1), ctx))
                {
                    flint_printf("check discrepancy\n"
                                       "i = %wd, j = %wd, k = %wd\n", i, j, k);
                    printf("v: "); fmpz_mod_poly_print_pretty(v, "#", ctx); printf("\n");
                    printf("B: "); fmpz_mod_berlekamp_massey_print(B1, ctx); printf("\n");
                    fflush(stdout);
                    flint_abort();
                }
                fmpz_mod_poly_clear(u, ctx);
                fmpz_mod_poly_clear(v, ctx);
                fmpz_mod_poly_clear(s, ctx);
                fmpz_mod_poly_clear(q, ctx);
                fmpz_mod_poly_clear(r, ctx);
            }

            fmpz_clear(p);
            fmpz_clear(apoint);
        }
        fmpz_mod_berlekamp_massey_clear(B1, ctx);
        fmpz_mod_berlekamp_massey_clear(B2, ctx);
    }

    fmpz_mod_ctx_clear(ctx);

    TEST_FUNCTION_END(state);
}
