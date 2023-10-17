/*
    Copyright (C) 2010 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "ulong_extras.h"

#ifndef check
#define check check
void check(mp_limb_t n, int s1, int s2)
{
    if (s1 != s2)
    {
        flint_printf("FAIL:\n");
        flint_printf("%wu: got %d instead of %d\n", n, s1, s2);
        fflush(stdout);
        flint_abort();
    }
}
#endif

TEST_FUNCTION_START(n_is_squarefree, state)
{
    int s, k;

    check(0, n_is_squarefree(0), 0);
    check(1, n_is_squarefree(1), 1);
    check(2, n_is_squarefree(2), 1);
    check(3, n_is_squarefree(3), 1);
    check(4, n_is_squarefree(4), 0);
    check(5, n_is_squarefree(5), 1);

    check(16, n_is_squarefree(16), 0);
    check(25, n_is_squarefree(25), 0);
    check(49, n_is_squarefree(49), 0);
    check(16*3, n_is_squarefree(16*3), 0);
    check(25*3, n_is_squarefree(25*3), 0);
    check(49*3, n_is_squarefree(49*3), 0);

    check(101*103, n_is_squarefree(101*103), 1);
    check(101*101, n_is_squarefree(101*101), 0);
    check(101*103*4, n_is_squarefree(101*103*4), 0);
    check(101*103*5, n_is_squarefree(101*103*5), 1);
    check(101*103*103*5, n_is_squarefree(101*103*103*5), 0);
    check(101*103*25, n_is_squarefree(101*103*25), 0);

    s = 0;
    for (k = 0; k <= 10000; k++)
        s += n_is_squarefree(k);

    if (s != 6083)
    {
        flint_printf("FAIL:\n");
        flint_printf("expected %d squarefree numbers <= 10000 (got %d)\n", 6083, s);
        fflush(stdout);
        flint_abort();
    }

    TEST_FUNCTION_END(state);
}
