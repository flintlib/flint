/*
    Copyright (C) 2016 Vincent Delecroix

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "nmod.h"
#include "nmod_poly.h"

TEST_FUNCTION_START(nmod_poly_power_sums_naive, state)
{
    int l, result;
    mp_limb_t i, j, k, tot;

    /* Check that it is valid in degree 3 with integer roots, ie */
    /* for polynomials of the form (x-i)(x-j)(x-k)               */
    for (i = 0; i < 4; i++)
        for (j = 0; j < 4; j++)
            for (k = 0; k < 4; k++)
            {
                mp_limb_t n;
                nmod_t mod;
                nmod_poly_t a, b, c, d;

                do {
                    n = n_randtest_prime(state, 1);
                } while(n < 4);
                nmod_init(&mod, n);

                nmod_poly_init(a, n);
                nmod_poly_init(b, n);
                nmod_poly_init(c, n);
                nmod_poly_init(d, n);

                nmod_poly_randtest(b, state, 40);
                nmod_poly_randtest(c, state, 40);
                nmod_poly_randtest(d, state, 40);

                nmod_poly_set_coeff_ui(a, 0, nmod_neg((i*j*k) % n, mod));
                nmod_poly_set_coeff_ui(a, 1, i * j + i * k + j * k);
                nmod_poly_set_coeff_ui(a, 2, nmod_neg((i + j + k) % n, mod));
                nmod_poly_set_coeff_ui(a, 3, 1);

                nmod_poly_power_sums_naive(b, a, 20);
                nmod_poly_set(d, a);
                nmod_poly_power_sums_naive(d, d, 20);

                for (l = 0; l < FLINT_MIN(20, nmod_poly_length(b)); l++)
                {
                    tot = nmod_add(nmod_pow_ui(i, l, mod),
                                   nmod_pow_ui(j, l, mod), mod);
                    tot = nmod_add(tot, nmod_pow_ui(k, l, mod), mod);

                    result = nmod_poly_get_coeff_ui(b, l) == tot &&
                        nmod_poly_get_coeff_ui(d, l) == tot;
                    if (!result)
                    {
                        flint_printf("FAIL: power sums integral root\n");
                        flint_printf("%d %d %d %d\n", i, j, k, l);
                        flint_printf("a = "), nmod_poly_print(a),
                            flint_printf("\n");
                        flint_printf("b = "), nmod_poly_print(b),
                            flint_printf("\n");
                        flint_printf("d = "), nmod_poly_print(d),
                            flint_printf("\n");
                        fflush(stdout);
                        flint_abort();
                    }
                }

                nmod_poly_power_sums_to_poly_naive(c, b);
                nmod_poly_set(d, b);
                nmod_poly_power_sums_to_poly_naive(d, d);
                result = nmod_poly_equal(a, c) && nmod_poly_equal(a, d);
                if (!result)
                {
                    flint_printf("FAIL: power sums to poly naive\n");
                    flint_printf("a = "), nmod_poly_print(a),
                        flint_printf("\n\n");
                    flint_printf("b = "), nmod_poly_print(b),
                        flint_printf("\n\n");
                    flint_printf("c = "), nmod_poly_print(c),
                        flint_printf("\n\n");
                    flint_printf("d = "), nmod_poly_print(c),
                        flint_printf("\n\n");
                    fflush(stdout);
                    flint_abort();
                }

                nmod_poly_clear(a);
                nmod_poly_clear(b);
                nmod_poly_clear(c);
                nmod_poly_clear(d);
            }

    /* Check that going back and forth between the power sums representation gives the identity */
    for (i = 0; i < 50 * flint_test_multiplier(); i++)
    {
        nmod_poly_t a, b, c, d;
        mp_limb_t n;

        do{
            n = n_randtest_prime(state, 1);
        }while(n < 50);

        nmod_poly_init(a, n);
        nmod_poly_init(b, n);
        nmod_poly_init(c, n);
        nmod_poly_init(d, n);

        nmod_poly_randtest(c, state, 50);
        nmod_poly_randtest(d, state, 50);

        nmod_poly_randtest_not_zero(a, state, 1 + n_randint(state, 20));
        nmod_poly_make_monic(a, a);

        nmod_poly_power_sums_naive(b, a, 30);
        nmod_poly_power_sums_to_poly_naive(c, b);

        nmod_poly_set(d, a);
        nmod_poly_power_sums_naive(d, d, 30);
        nmod_poly_power_sums_to_poly_naive(d, d);

        result = nmod_poly_equal(a, c) && nmod_poly_equal(a, d);
        if (!result)
        {
            flint_printf("FAIL: power sums - power sums to poly\n");
            flint_printf("a = "), nmod_poly_print(a), flint_printf("\n\n");
            flint_printf("b = "), nmod_poly_print(b), flint_printf("\n\n");
            flint_printf("c = "), nmod_poly_print(c), flint_printf("\n\n");
            flint_printf("d = "), nmod_poly_print(c), flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        nmod_poly_clear(a);
        nmod_poly_clear(b);
        nmod_poly_clear(c);
        nmod_poly_clear(d);
    }

    /* Check that the product of polynomials correspond to the sum of Power sums series */
    for (i = 0; i < 20 * flint_test_multiplier(); i++)
    {
        nmod_poly_t a, b, c, d;
        mp_limb_t n;

        do{
            n = n_randtest_prime(state, 1);
        }while(n < 20);

        nmod_poly_init(a, n);
        nmod_poly_init(b, n);
        nmod_poly_init(c, n);
        nmod_poly_init(d, n);

        nmod_poly_randtest_not_zero(a, state, 1 + n_randint(state, 10));
        nmod_poly_randtest_not_zero(b, state, 1 + n_randint(state, 10));
        nmod_poly_randtest(c, state, 30);
        nmod_poly_randtest(d, state, 30);

        nmod_poly_mul(c, a, b);
        nmod_poly_power_sums_naive(c, c, 20);

        /* NOTE: the code path is not the same if the polynomial is monic. We let only a be monic */
        nmod_poly_make_monic(a, a);
        nmod_poly_power_sums_naive(a, a, 20);
        nmod_poly_power_sums_naive(b, b, 20);
        nmod_poly_add(d, a, b);

        result = nmod_poly_equal(c, d);
        if (!result)
        {
            flint_printf
                ("FAIL: PowerSums(p1 p2) = PowerSums(p1) + PowerSums(p2)\n");
            flint_printf("a = ");
            nmod_poly_print(a), flint_printf("\n");
            flint_printf("b = ");
            nmod_poly_print(b), flint_printf("\n");
            flint_printf("c = ");
            nmod_poly_print(c), flint_printf("\n");
            flint_printf("d = ");
            nmod_poly_print(d), flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }

        nmod_poly_clear(a);
        nmod_poly_clear(b);
        nmod_poly_clear(c);
        nmod_poly_clear(d);
    }

    TEST_FUNCTION_END(state);
}
