/*
    Authored 2015 by Daniel S. Roche; US Government work in the public domain.

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz.h"

const char * const manual_tests[] = {
    /* anything <= 1 should produce 2 */
    "1", "2",
    "0", "2",
    "-13842090335966649306", "2",
    "-819856806963901485525117166411138935256434837055304839519773", "2",
    /* some small checks */
    "2", "3",
    "3", "5",
    "4", "5",
    /* around the 32-bit boundary */
    "2147483630", "2147483647",
    "2147483648", "2147483659",
    "4294967279", "4294967291",
    "4294967296", "4294967311",
    /* around the 64-bit boundary */
    "9223372036854775643", "9223372036854775783",
    "9223372036854775808", "9223372036854775837",
    "18446744073709551533", "18446744073709551557",
    "18446744073709551616", "18446744073709551629",
    "\0"
};

TEST_FUNCTION_START(fmpz_nextprime, state)
{
    int i;

    for (i=0; manual_tests[i][0]; i += 2)
    {
        fmpz_t start, expected, actual;

        fmpz_init(start);
        fmpz_init(expected);
        fmpz_init(actual);

        fmpz_set_str(start, manual_tests[i], 10);
        fmpz_nextprime(actual, start, 0);
        fmpz_set_str(expected, manual_tests[i+1], 10);

        if (!fmpz_equal(actual, expected) || !_fmpz_is_canonical(actual))
        {
            flint_printf("FAIL:\n");
            fmpz_print(start); flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_clear(start);
        fmpz_clear(expected);
        fmpz_clear(actual);
    }

    for (i=0; i < 100 * flint_test_multiplier(); ++i)
    {
        fmpz_t start, res, iter;

        fmpz_init(start);
        fmpz_init(res);
        fmpz_init(iter);

        fmpz_randtest_unsigned(start, state, 160);
        fmpz_nextprime(res, start, 0);

        fmpz_set(iter, start);
        do fmpz_add_ui(iter, iter, UWORD(1));
        while (!fmpz_is_probabprime(iter));

        if (!fmpz_equal(res, iter) || !_fmpz_is_canonical(res))
        {
            flint_printf("FAIL:\n");
            fmpz_print(start); flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_clear(start);
        fmpz_clear(res);
        fmpz_clear(iter);
    }

    TEST_FUNCTION_END(state);
}
