/*
    Copyright (C) 2016 Vincent Delecroix

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "ulong_extras.h"
#include "nmod_poly.h"

TEST_FUNCTION_START(nmod_poly_power_sums, state)
{
    int i, result;

    /* Check that the different version coincide and aliasing in nmod_poly_power_sums */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        mp_limb_t n;
        nmod_poly_t a, b, c, d, e;

        do{
            n = n_randtest_prime(state, 1);
        }while(n < 10);

        nmod_poly_init(a, n);
        nmod_poly_init(b, n);
        nmod_poly_init(c, n);
        nmod_poly_init(d, n);
        nmod_poly_init(e, n);

        nmod_poly_randtest_not_zero(a, state, 1 + n_randint(state, 10));

        nmod_poly_power_sums_naive(b, a, 20);
        nmod_poly_power_sums_schoenhage(c, a, 20);
        nmod_poly_power_sums(d, a, 20);

        nmod_poly_set(e, a);
        nmod_poly_power_sums(e, e, 20);

        result = nmod_poly_equal(b, c) && nmod_poly_equal(b, d) &&
            nmod_poly_equal(b, e);
        if (!result)
        {
            flint_printf
                ("FAIL: PowerSums(p1 p2) = PowerSums(p1) + PowerSums(p2)\n");
            flint_printf("a = "), nmod_poly_print(a), flint_printf("\n");
            flint_printf("b = "), nmod_poly_print(b), flint_printf("\n");
            flint_printf("c = "), nmod_poly_print(c), flint_printf("\n");
            flint_printf("d = "), nmod_poly_print(d), flint_printf("\n");
            flint_printf("e = "), nmod_poly_print(e), flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }

        nmod_poly_make_monic(a, a);

        nmod_poly_power_sums_to_poly_schoenhage(b, b);
        nmod_poly_power_sums_to_poly(d, c);
        nmod_poly_power_sums_to_poly(c, c);
        result = nmod_poly_equal(a, b) && nmod_poly_equal(a, c) &&
            nmod_poly_equal(a, d);
        if (!result)
        {
            flint_printf("FAIL: power_sums_to_poly\n");
            flint_printf("e = "), nmod_poly_print(e), flint_printf("\n");
            flint_printf("a = "), nmod_poly_print(a), flint_printf("\n");
            flint_printf("b = "), nmod_poly_print(b), flint_printf("\n");
            flint_printf("c = "), nmod_poly_print(c), flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }

        nmod_poly_clear(a);
        nmod_poly_clear(b);
        nmod_poly_clear(c);
        nmod_poly_clear(d);
        nmod_poly_clear(e);
    }

    TEST_FUNCTION_END(state);
}
