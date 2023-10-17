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
#include "nmod_poly.h"

TEST_FUNCTION_START(nmod_poly_berlekamp_massey, state)
{
    slong i, j, k, l;

    {
        int changed;
        nmod_berlekamp_massey_t B;
        nmod_berlekamp_massey_init(B, 101);
        nmod_berlekamp_massey_add_point(B, 1);
        nmod_berlekamp_massey_add_point(B, 1);
        nmod_berlekamp_massey_add_point(B, 2);
        nmod_berlekamp_massey_add_point(B, 3);
        changed = nmod_berlekamp_massey_reduce(B);
        if (   changed != 1
            || 2 != nmod_poly_degree(nmod_berlekamp_massey_V_poly(B))
            || nmod_poly_degree(nmod_berlekamp_massey_R_poly(B)) > 1)
        {
            printf("FAIL\n");
            flint_printf("check fibonacci 1\n");
            fflush(stdout);
            flint_abort();
        }
        nmod_berlekamp_massey_add_point(B, 5);
        nmod_berlekamp_massey_add_point(B, 8);
        changed = nmod_berlekamp_massey_reduce(B);
        if (   changed != 0
            || 2 != nmod_poly_degree(nmod_berlekamp_massey_V_poly(B))
            || nmod_poly_degree(nmod_berlekamp_massey_R_poly(B)) > 1)
        {
            printf("FAIL\n");
            flint_printf("check fibonacci 2\n");
            fflush(stdout);
            flint_abort();
        }
        nmod_berlekamp_massey_clear(B);
    }

    for (i = 0; i < 20 * flint_test_multiplier(); i++)
    {
        nmod_berlekamp_massey_t B1, B2;
        nmod_berlekamp_massey_init(B1, 2);
        nmod_berlekamp_massey_init(B2, 2);
        for (j = 0; j < 10; j++)
        {
            mp_limb_t p;

            p = n_randtest_prime(state, 1);

            nmod_berlekamp_massey_set_prime(B1, p);
            nmod_berlekamp_massey_set_prime(B2, p);

            /* check intermediate reductions match */
            for (k = 0; k < 10; k++)
            {
                nmod_berlekamp_massey_add_point(B1, n_randint(state, p));
                nmod_berlekamp_massey_add_zeros(B1, n_randint(state, 5));
                nmod_berlekamp_massey_add_point(B1, n_randint(state, p));
                if (n_randlimb(state) & 1)
                {
                    nmod_berlekamp_massey_reduce(B1);
                }
            }
            nmod_berlekamp_massey_add_points(B2, nmod_berlekamp_massey_points(B1),
                                    nmod_berlekamp_massey_point_count(B1));
            nmod_berlekamp_massey_reduce(B2);
            nmod_berlekamp_massey_reduce(B1);
            if (!nmod_poly_equal(nmod_berlekamp_massey_V_poly(B1),
                                 nmod_berlekamp_massey_V_poly(B2)))
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
                nmod_poly_t u, v, s, q, r;
                nmod_poly_init(u, p);
                nmod_poly_init(v, p);
                nmod_poly_init(s, p);
                nmod_poly_init(q, p);
                nmod_poly_init(r, p);

                /* deg(u) < deg(v), deg(v) = k */
                nmod_poly_randtest(u, state, k);
                nmod_poly_randtest_monic(v, state, k + 1);

                /* q has enough coefficients of expansion of u/v at infty */
                nmod_poly_shift_left(s, u, 2*k);
                nmod_poly_divrem(q, r, s, v);

                nmod_berlekamp_massey_start_over(B1);
                for (l = 2*k - 1; l >= 0 ; l--)
                {
                    nmod_berlekamp_massey_add_point(B1, nmod_poly_get_coeff_ui(q, l));
                }
                nmod_berlekamp_massey_reduce(B1);
                nmod_poly_divrem(q, r, v, nmod_berlekamp_massey_V_poly(B1));
                if (!nmod_poly_is_zero(r))
                {
                    flint_printf("check berlekamp_massey does its job\n"
                                       "i = %wd, j = %wd, k = %wd\n", i, j, k);
                    printf("v: "); nmod_poly_print_pretty(v, "#"); printf("\n");
                    printf("B: "); nmod_berlekamp_massey_print(B1); printf("\n");
                    fflush(stdout);
                    flint_abort();
                }
                if (   nmod_poly_degree(nmod_berlekamp_massey_R_poly(B1))
                    >= nmod_poly_degree(nmod_berlekamp_massey_V_poly(B1)))
                {
                    flint_printf("check discrepancy\n"
                                       "i = %wd, j = %wd, k = %wd\n", i, j, k);
                    printf("v: "); nmod_poly_print_pretty(v, "#"); printf("\n");
                    printf("B: "); nmod_berlekamp_massey_print(B1); printf("\n");
                    fflush(stdout);
                    flint_abort();
                }

                nmod_poly_clear(u);
                nmod_poly_clear(v);
                nmod_poly_clear(s);
                nmod_poly_clear(q);
                nmod_poly_clear(r);
            }
        }
        nmod_berlekamp_massey_clear(B1);
        nmod_berlekamp_massey_clear(B2);
    }

    TEST_FUNCTION_END(state);
}
