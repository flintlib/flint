/*
    Copyright (C) 2019 Daniel Schultz
    This file is part of FLINT.
    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "ulong_extras.h"

int main(void)
{
    int i, j, result;
    FLINT_TEST_INIT(state);

    flint_printf("add_ssssaaaaaaaa....");
    fflush(stdout);

    for (i = 0; i < 1000000; i++)
    {
        mp_limb_t s[4], t[4], a[4], b[4];

        for (j = 0; j < 4; j++)
        {
            s[j] = n_randtest(state);
            t[j] = n_randtest(state);
            a[j] = n_randtest(state);
            b[j] = n_randtest(state);
        }

        add_ssssaaaaaaaa(s[3], s[2], s[1], s[0], a[3], a[2], a[1], a[0],
                                                 b[3], b[2], b[1], b[0]);

        mpn_add_n(t, a, b, 4);

        result = ((s[3] == t[3]) && (s[2] == t[2]) && (s[1] == t[1]) && (s[0] == t[0]));
        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("a[3] = %wu, a[2] = %wu, a[1] = %wu, a[0] = %wu\n", a[3], a[2], a[1], a[0]);
            flint_printf("b[3] = %wu, b[2] = %wu, b[1] = %wu, b[0] = %wu\n", b[3], b[2], b[1], b[0]);
            flint_printf("s[3] = %wu, s[2] = %wu, s[1] = %wu, s[0] = %wu\n", s[3], s[2], s[1], s[0]);
            flint_printf("t[3] = %wu, t[2] = %wu, t[1] = %wu, t[0] = %wu\n", t[3], t[2], t[1], t[0]);
            flint_abort();
        }

        for (j = 0; j < 4; j++)
        {
            s[j] = n_randtest(state);
            t[j] = n_randtest(state);
            a[j] = n_randtest(state);
            b[j] = n_randtest(state);
        }

        mpn_add_n(t, s, b, 4);

        add_ssssaaaaaaaa(s[3], s[2], s[1], s[0], s[3], s[2], s[1], s[0],
                                                 b[3], b[2], b[1], b[0]);

        result = ((s[3] == t[3]) && (s[2] == t[2]) && (s[1] == t[1]) && (s[0] == t[0]));
        if (!result)
        {
            flint_printf("FAIL:\naliasing\n");
            flint_printf("a[3] = %wu, a[2] = %wu, a[1] = %wu, a[0] = %wu\n", a[3], a[2], a[1], a[0]);
            flint_printf("b[3] = %wu, b[2] = %wu, b[1] = %wu, b[0] = %wu\n", b[3], b[2], b[1], b[0]);
            flint_printf("s[3] = %wu, s[2] = %wu, s[1] = %wu, s[0] = %wu\n", s[3], s[2], s[1], s[0]);
            flint_printf("t[3] = %wu, t[2] = %wu, t[1] = %wu, t[0] = %wu\n", t[3], t[2], t[1], t[0]);
            flint_abort();
        }
    }

    FLINT_TEST_CLEANUP(state);

    flint_printf("PASS\n");
    return 0;
}
