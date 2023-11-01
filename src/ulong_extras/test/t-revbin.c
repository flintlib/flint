/*
    Copyright (C) 2015 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "ulong_extras.h"

ulong
n_revbin_naive(ulong n, ulong b)
{
    ulong r = 0, i;

    for (i = 0; i < b; i++)
    {
        r <<= 1;
        r += (n & 1);
        n >>= 1;
    }

    return r;
}

TEST_FUNCTION_START(n_revbin, state)
{
    int i, result;

    /* 0 bits */
    for (i = 0; i < 10 * flint_test_multiplier(); i++)
    {
        ulong n, b, r;

        n = n_randlimb(state);
        b = 0;

        r = n_revbin(n, b);

        result = (r == 0);
        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("b = %wu\n", b);
            flint_printf("n = %wx\n", n);
            flint_printf("r = %wx\n", r);
            fflush(stdout);
            flint_abort();
        }
    }

    /* at most 8 bits */
    for (i = 0; i < 10 * flint_test_multiplier(); i++)
    {
        ulong n, b, d, r1, r2, mask;

        for (b = 1; b <= 8; b++)
        {
            mask = ~((UWORD(1) << b) - UWORD(1));   /* 1111..11100..0 */

            for (d = 0; d < UWORD(1) << b; d++)
            {
                /* garbage in top FLINT_BITS - b bits */
                n = (n_randlimb(state) & mask);
                n |= d;

                r1 = n_revbin(n, b);
                r2 = n_revbin_naive(d, b);

                result = (r1 == r2);
                if (!result)
                {
                    flint_printf("FAIL:\n");
                    flint_printf("b = %wu\n", b);
                    flint_printf("n = %wx\n", n);
                    flint_printf("d = %wx\n", d);
                    flint_printf("r1 = %wx\n", r1);
                    flint_printf("r2 = %wx\n", r2);
                    fflush(stdout);
                    flint_abort();
                }
            }
        }
    }

    /* random number of bits */
    for (i = 0; i < 100; i++)
    {
        ulong n, b, d, j, r1, r2, mask;

        b = n_randint(state, FLINT_BITS) + 1;

        /* 1111..11100..0 */
        mask = b == FLINT_BITS ? UWORD(0) : ~((UWORD(1) << b) - UWORD(1));

        for (j = 0; j < 1000 * flint_test_multiplier(); j++)
        {
            /* garbage in top FLINT_BITS - b bits */
            n = (n_randlimb(state) & mask);
            d = n_randbits(state, n_randint(state, b) + 1);
            n |= d;

            r1 = n_revbin(n, b);
            r2 = n_revbin_naive(d, b);

            result = (r1 == r2);
            if (!result)
            {
                flint_printf("FAIL:\n");
                flint_printf("b = %wu\n", b);
                flint_printf("n = %wx\n", n);
                flint_printf("d = %wx\n", d);
                flint_printf("r1 = %wx\n", r1);
                flint_printf("r2 = %wx\n", r2);
                fflush(stdout);
                flint_abort();
            }
        }
    }

    TEST_FUNCTION_END(state);
}
