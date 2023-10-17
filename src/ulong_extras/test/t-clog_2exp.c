/*
    Copyright (C) 2021 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "ulong_extras.h"

TEST_FUNCTION_START(n_clog_2exp, state)
{
    slong i;
    ulong t[][3] = {{1, 2, 1},
                    {2, 2, 2},
                    {698, 127, 100},
                    {699, 127, 101},
                    {700, 127, 101},
                    {701, 127, 101},
                    {699, 128, 100},
                    {700, 128, 100},
                    {701, 128, 101},
                    {694, 129,  99},
                    {695, 129, 100},
                    {696, 129, 100},
                    {697, 129, 100},
                    {698, 129, 100},
                    {699, 129, 100},
                    {700, 129, 100},
                    {701, 129, 100},
                    {702, 129, 101},
                    {0, UWORD_MAX, 0},
                    {1, UWORD_MAX, 1},
                    {2, UWORD_MAX, 1},
                    {FLINT_BITS - 1, UWORD_MAX, 1},
                    {FLINT_BITS, UWORD_MAX, 2}};

    for (i = 0; i < sizeof(t)/sizeof(t[0]); i++)
    {
        ulong r = n_clog_2exp(t[i][0], t[i][1]);

        if (r != t[i][2])
        {
            flint_printf("FAIL:\n");
            flint_printf("clog_2exp(%wu, %wu) = %wu\n", t[i][0], t[i][1], t[i][2]);
            flint_printf("but computed %wu\n", r);
        }
    }

    TEST_FUNCTION_END(state);
}
