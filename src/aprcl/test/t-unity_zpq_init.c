/*
    Copyright (C) 2015 Vladimir Glazachev

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz.h"
#include "aprcl.h"

TEST_FUNCTION_START(aprcl_unity_zpq_init, state)
{
    int i;

    for (i = 0; i < 10 * flint_test_multiplier(); i++)
    {
        ulong p, q;
        fmpz_t n;
        unity_zpq f;

        fmpz_init(n);

        p = n_randprime(state, 2 + n_randint(state, 16), 0);
        q = n_randprime(state, 2 + n_randint(state, 16), 0);

        fmpz_randtest_unsigned(n, state, 200);
        fmpz_add_ui(n, n, 1);

        unity_zpq_init(f, q, p, n);
        unity_zpq_clear(f);

        fmpz_clear(n);
    }

    TEST_FUNCTION_END(state);
}
