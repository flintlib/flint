/*
    Copyright (C) 2014 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz_mat.h"

int should_have_hadamard(int n)
{
    if (n <= 2)
        return 1;

    if (n % 4 != 0)
        return 0;

    if (n <= 300 && n != 92 && n != 116 && n != 156 && n != 172 && n != 184 &&
        n != 188 && n != 232 && n != 236 && n != 260 && n != 268 && n != 292)
        return 1;

    return 0;
}

TEST_FUNCTION_START(fmpz_mat_hadamard, state)
{
    int n;

    for (n = 0; n <= 300; n++)
    {
        fmpz_mat_t h;
        int success;

        fmpz_mat_init(h, n, n);
        success = fmpz_mat_hadamard(h);

        if (success)
        {
            if (!fmpz_mat_is_hadamard(h))
            {
                printf("FAIL: output is not a Hadamard matrix\n");
                printf("n = %d\n\n", n);
                fmpz_mat_print_pretty(h); printf("\n\n");
                fflush(stdout);
                flint_abort();
            }
        }
        else if (should_have_hadamard(n))
        {
            printf("FAIL: expected Hadamard matrix of size %d to work\n\n", n);
            fflush(stdout);
            flint_abort();
        }

        fmpz_mat_clear(h);
    }

    TEST_FUNCTION_END(state);
}
