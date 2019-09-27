/*
    Copyright (C) 2019 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include <stdlib.h>
#include "fmpz_mod_poly.h"

int
main(void)
{
    slong i, j, k, l;
    FLINT_TEST_INIT(state);

    flint_printf("berlekamp_massey....");
    fflush(stdout);

    {
        fmpz_mod_berlekamp_massey_t B;
        fmpz_mod_berlekamp_massey_init_ui(B, 101);
        fmpz_mod_berlekamp_massey_add_point_ui(B, 1);
        fmpz_mod_berlekamp_massey_add_point_ui(B, 1);
        fmpz_mod_berlekamp_massey_add_point_ui(B, 2);
        fmpz_mod_berlekamp_massey_add_point_ui(B, 3);
        fmpz_mod_berlekamp_massey_add_point_ui(B, 5);
        fmpz_mod_berlekamp_massey_reduce(B);
        if (2 != fmpz_mod_poly_degree(fmpz_mod_berlekamp_massey_V_poly(B)))
        {
            printf("FAIL\n");
            flint_printf("check fibonacci\n");
            flint_abort();
        }
        fmpz_mod_berlekamp_massey_clear(B);
    }

    for (i = 0; i < 15 * flint_test_multiplier(); i++)
    {
        fmpz_mod_berlekamp_massey_t B1, B2;
        fmpz_mod_berlekamp_massey_init_ui(B1, 2);
        fmpz_mod_berlekamp_massey_init_ui(B2, 2);
        for (j = 0; j < 10; j++)
        {
            fmpz_t p, apoint;

            fmpz_init(apoint);
            fmpz_init_set_ui(p, n_randtest_prime(state, 1));

            fmpz_mod_berlekamp_massey_set_prime(B1, p);
            fmpz_mod_berlekamp_massey_set_prime(B2, p);

            /* check intermediate reductions match */
            for (k = 0; k < 10; k++)
            {
                fmpz_randm(apoint, state, p);
                fmpz_mod_berlekamp_massey_add_point(B1, apoint);
                fmpz_mod_berlekamp_massey_add_zeros(B1, n_randint(state, 5));
                fmpz_randm(apoint, state, p);
                fmpz_mod_berlekamp_massey_add_point(B1, apoint);
                if (n_randlimb(state) & 1)
                {
                    fmpz_mod_berlekamp_massey_reduce(B1);
                }
            }
            fmpz_mod_berlekamp_massey_add_points(B2,
                                    fmpz_mod_berlekamp_massey_points(B1),
                                    fmpz_mod_berlekamp_massey_point_count(B1));
            fmpz_mod_berlekamp_massey_reduce(B2);
            fmpz_mod_berlekamp_massey_reduce(B1);
            if (!fmpz_mod_poly_equal(fmpz_mod_berlekamp_massey_V_poly(B1),
                                     fmpz_mod_berlekamp_massey_V_poly(B2)))
            {
                printf("FAIL\n");
                flint_printf("check intermediate reductions match\n"
                                                   "i = %wd, j = %wd\n", i, j);
                flint_abort();
            }

            /*
                Check berlekamp-massey does its job - 2k coeffcients of

                    u     a1    a2          a(2k)
                   --- = --- + --- + ... + ------ + ...
                    v     x    x^2         x^(2k)

                should be sufficient to reconstruct a divisor of v
            */
            for (k = 0; k < 15; k++)
            {
                fmpz_mod_poly_t u, v, s, q, r;
                fmpz_mod_poly_init(u, p);
                fmpz_mod_poly_init(v, p);
                fmpz_mod_poly_init(s, p);
                fmpz_mod_poly_init(q, p);
                fmpz_mod_poly_init(r, p);

                /* deg(u) < deg(v), deg(v) = k */
                fmpz_mod_poly_randtest(u, state, k);
                fmpz_mod_poly_randtest_monic(v, state, k + 1);

                /* q has enough coefficients of expansion of u/v at infty */
                fmpz_mod_poly_shift_left(s, u, 2*k);
                fmpz_mod_poly_divrem(q, r, s, v);
                
                fmpz_mod_berlekamp_massey_start_over(B1);
                for (l = 2*k - 1; l >= 0 ; l--)
                {
                    fmpz_mod_poly_get_coeff_fmpz(apoint, q, l);
                    fmpz_mod_berlekamp_massey_add_point(B1, apoint);
                }
                fmpz_mod_berlekamp_massey_reduce(B1);
                fmpz_mod_poly_divrem(q, r, v, fmpz_mod_berlekamp_massey_V_poly(B1));
                if (!fmpz_mod_poly_is_zero(r))
                {
                    flint_printf("check berlekamp_massey does its job\n"
                                       "i = %wd, j = %wd, k = %wd\n", i, j, k);
                    printf("v: "); fmpz_mod_poly_print_pretty(v, "#"); printf("\n");
                    printf("B: "); fmpz_mod_berlekamp_massey_print(B1); printf("\n");
                    flint_abort();
                }
                if (   fmpz_mod_poly_degree(fmpz_mod_berlekamp_massey_R_poly(B1))
                    >= fmpz_mod_poly_degree(fmpz_mod_berlekamp_massey_V_poly(B1)))
                {
                    flint_printf("check discrepancy\n"
                                       "i = %wd, j = %wd, k = %wd\n", i, j, k);
                    printf("v: "); fmpz_mod_poly_print_pretty(v, "#"); printf("\n");
                    printf("B: "); fmpz_mod_berlekamp_massey_print(B1); printf("\n");
                    flint_abort();
                }
                fmpz_mod_poly_clear(u);
                fmpz_mod_poly_clear(v);
                fmpz_mod_poly_clear(s);
                fmpz_mod_poly_clear(q);
                fmpz_mod_poly_clear(r);
            }

            fmpz_clear(p);
            fmpz_clear(apoint);
        }
        fmpz_mod_berlekamp_massey_clear(B1);
        fmpz_mod_berlekamp_massey_clear(B2);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}
