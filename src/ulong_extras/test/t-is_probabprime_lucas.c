/*
    Copyright (C) 2009 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "gmpcompat.h"
#include "ulong_extras.h"

TEST_FUNCTION_START(n_is_probabprime_lucas, state)
{
    int i, result;
    ulong count = UWORD(0);
    ulong d;
    mpz_t d_m;
    slong test_multiplier;

    test_multiplier = FLINT_MAX(1, flint_test_multiplier());

    {
        /* A selection of Lucas pseudoprimes. These are known composites for
           which the Lucas test *should* claim PRIME, but which a buggy
           implementation may claim COMPOSITE for spurious reasons. */
        ulong pseudo[] = {
            /* Initial entries from Dana Jacobsen's table (http://ntheory.org/pseudoprimes.html) */
            323, 377, 1159, 1829, 3827, 5459, 5777, 9071, 9179, 10877, 11419, 11663,
            13919, 14839, 16109, 16211, 18407, 18971, 19043, 22499, 23407, 24569,
            25199, 25877, 26069, 27323, 32759, 34943, 35207, 39059, 39203, 39689,
            /* Some bigger entries (24 bits) */
            8711699, 10402139, 12124559, 13695947, 13826231,
            /* Even bigger entries (32 bits) */
            UWORD(2877330059), UWORD(3816635327), UWORD(3964016279),
#if FLINT_BITS == 64
            /* Check that #2168 is fixed. */
            UWORD(9508976851322519),
#endif
            0
        };

        for (i = 0; (d = pseudo[i]) != 0; i++)
        {
            if (!n_is_probabprime_lucas(d))
                TEST_FUNCTION_FAIL("Lucas pseudoprime %wu claimed composite by Lucas test\n", d);
        }
    }

    for (i = 0; i < 10000 * test_multiplier; i++) /* Test that primes pass the test */
    {
        mpz_init(d_m);

        do
        {
            d = n_randtest_not_zero(state);
            flint_mpz_set_ui(d_m, d);
            mpz_nextprime(d_m, d_m);
            d = flint_mpz_get_ui(d_m);
        } while (mpz_size(d_m) > 1);

        result = n_is_probabprime_lucas(d);
        if (!result)
            TEST_FUNCTION_FAIL("d = %wu is declared composite\n", d);

        mpz_clear(d_m);
    }

    for (i = 0; i < 10000 * test_multiplier; i++) /* Test that not too many composites pass */
    {
        mpz_init(d_m);

        do
        {
            d = n_randtest(state);
            flint_mpz_set_ui(d_m, d);
        } while (mpz_probab_prime_p(d_m, 12));

        if (n_is_probabprime_lucas(d) == 1) count++;

        mpz_clear(d_m);
    }

    result = (count < 20 * test_multiplier);
    if (!result)
        TEST_FUNCTION_FAIL("%wu composites declared prime\n", count);

    TEST_FUNCTION_END(state);
}
