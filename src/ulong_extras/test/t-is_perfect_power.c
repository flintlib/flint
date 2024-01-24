/*
    Copyright (C) 2009, 2017 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "gmpcompat.h"
#include "ulong_extras.h"

TEST_FUNCTION_START(n_is_perfect_power, state)
{
    int i, result;
    ulong bits, root, hi, n;
    mp_limb_t d;

    for (i = 0; i < 1000 * flint_test_multiplier(); i++) /* Test that squares pass the test */
    {
        bits = n_randint(state, FLINT_BITS/2) + 1;
        d = n_randtest_bits(state, bits);

        result = n_is_perfect_power(&root, n_pow(d, 2));
        if (result == 0)
            TEST_FUNCTION_FAIL("d^2 = %wu is declared not a perfect power\n", d*d);

        if (n_pow(root, result) != n_pow(d, 2))
            TEST_FUNCTION_FAIL("%wu^%wu != %wu\n", root, result, d*d);
    }

    for (i = 0; i < 1000 * flint_test_multiplier(); i++) /* Test that cubes pass the test */
    {
        bits = n_randint(state, FLINT_BITS/3) + 1;
        d = n_randtest_bits(state, bits);

        result = n_is_perfect_power(&root, n_pow(d, 3));
        if (result == 0)
            TEST_FUNCTION_FAIL("d^3 = %wu is declared not a perfect power\n", d*d*d);

        if (n_pow(root, result) != n_pow(d, 3))
            TEST_FUNCTION_FAIL("%wu^%wu != %wu\n", root, result, d*d*d);
    }

    for (i = 0; i < 1000 * flint_test_multiplier(); i++) /* Test that fifth powers pass the test */
    {
        bits = n_randint(state, FLINT_BITS/5) + 1;
        d = n_randtest_bits(state, bits);

        result = n_is_perfect_power(&root, n_pow(d, 5));
        if (result == 0)
            TEST_FUNCTION_FAIL("d^5 = %wu is declared not a perfect power\n", d*d*d*d*d);

        if (n_pow(root, result) != n_pow(d, 5))
            TEST_FUNCTION_FAIL("%wu^%wu != %wu\n", root, result, d*d*d*d*d);
    }

    /* exhaustively test all other powers */
    for (d = 2; d < (UWORD(1) << (FLINT_BITS/5)); d++)
    {
        hi = 0;
        n = d*d;

        while (hi == 0)
        {
            result = n_is_perfect_power(&root, n);
            if (n_pow(root, result) != n)
                TEST_FUNCTION_FAIL("%wu^%wu != %wu\n", root, result, n);

            umul_ppmm(hi, n, n, d);
        }
    }

    for (i = 0; i < 10000 * flint_test_multiplier(); i++) /* Test that non perfect powers fail */
    {
        mpz_t d_m;
        mpz_init(d_m);

        do
        {
            d = n_randtest(state);
            flint_mpz_set_ui(d_m, d);
        } while (mpz_perfect_power_p(d_m));

        result = !n_is_perfect_power(&root, d);
        if (!result)
            TEST_FUNCTION_FAIL("d = %wu is declared a perfect power\n", d);

        mpz_clear(d_m);
    }

    TEST_FUNCTION_END(state);
}
