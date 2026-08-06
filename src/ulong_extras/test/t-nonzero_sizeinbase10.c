/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <string.h>
#include "test_helpers.h"
#include "gmpcompat.h"
#include "ulong_extras.h"

TEST_FUNCTION_START(n_nonzero_sizeinbase10, state)
{
    ulong n;
    slong rep, size1, size2;
    mpz_t t;
    char * str;

    mpz_init(t);
    str = flint_malloc((FLINT_BITS + 1) * sizeof(char));

    for (rep = 0; rep < 10000 * flint_test_multiplier(); rep++)
    {
        n = n_randtest_not_zero(state);

        size1 = n_nonzero_sizeinbase10(n);

        flint_mpz_set_ui(t, n);
        mpz_get_str(str, 10, t);
        size2 = strlen(str);

        if (size1 != size2)
            TEST_FUNCTION_FAIL(
                    "n = %wu\n"
                    "n_nonzero_sizeinbase10: %wd, strlen: %wd\n",
                    n, size1, size2);
    }

    flint_free(str);
    mpz_clear(t);

    TEST_FUNCTION_END(state);
}
