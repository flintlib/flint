/*
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz.h"
#include "arith.h"

TEST_FUNCTION_START(arith_bell_number_multi_mod, state)
{
    slong i;

    for (i = 0; i < 100; i++)
    {
        slong n;
        fmpz_t b1, b2;

        fmpz_init(b1);
        fmpz_init(b2);

        n = n_randint(state, 500);

        arith_bell_number_dobinski(b1, n);
        arith_bell_number_multi_mod(b2, n);

        if (!fmpz_equal(b1, b2))
        {
            flint_printf("FAIL:\n");
            flint_printf("n = %wd\n", n);
            fflush(stdout);
            flint_abort();
        }

        fmpz_clear(b1);
        fmpz_clear(b2);
    }

    TEST_FUNCTION_END(state);
}
