/*
    Copyright (C) 2015 Kushagra Singh

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "gmpcompat.h"
#include "ulong_extras.h"

TEST_FUNCTION_START(n_cbrtrem, state)
{
    int i, result;

    for (i = 0; i < 10000 * flint_test_multiplier(); i++)
    {
        mp_limb_t a, b, c, i, j;
        mpz_t e, f, g;

        mpz_init(e);
        mpz_init(f);
        mpz_init(g);

        c = n_randint(state, 0);    /*number */
        flint_mpz_set_ui(g, c);

        a = n_cbrtrem(&b, c);
        mpz_rootrem(e, f, g, 3);

        i = flint_mpz_get_ui(e);
        j = flint_mpz_get_ui(f);

        result = ((a == i) && (b == j));

        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("Passed Parameters : n = %wu", c);
            flint_printf("Answer generated : base = %wu remainder = %wu", a, b);
            flint_printf("Expected answer : base = %wu remainder = %wu", i, j);
            fflush(stdout);
            flint_abort();
        }
        mpz_clear(e);
        mpz_clear(f);
        mpz_clear(g);
    }

    TEST_FUNCTION_END(state);
}
