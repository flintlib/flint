/*
    Authored 2015 by Daniel S. Roche; US Government work in the public domain.
    Copyright (C) 2024 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz.h"

TEST_FUNCTION_START(fmpz_randprime, state)
{
    slong ix;
    int result;

    for (ix = 0; ix < 500 * flint_test_multiplier(); ix++)
    {
        /* Test that output is prime and that it has said number of bits */
        fmpz_t p;
        flint_bitcnt_t bits;
        int proved;

        fmpz_init(p);

        bits = 2 + n_randint(state, 100);
        proved = n_randint(state, 2);

        fmpz_randprime(p, state, bits, proved);

        result = (fmpz_is_prime(p) && (fmpz_bits(p) == bits));

        if (!result)
            TEST_FUNCTION_FAIL(
                    "bits = %{ulong}\n"
                    "p = %{fmpz}\n",
                    bits, p);

        fmpz_clear(p);
    }

    TEST_FUNCTION_END(state);
}
