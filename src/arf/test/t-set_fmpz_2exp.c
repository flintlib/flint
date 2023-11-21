/*
    Copyright (C) 2012 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "arf.h"

TEST_FUNCTION_START(arf_set_fmpz_2exp, state)
{
    slong iter;

    /* test exact roundtrip R -> Q -> R */
    for (iter = 0; iter < 100000 * 0.1 * flint_test_multiplier(); iter++)
    {
        slong bits;
        arf_t x, z;
        fmpz_t y, e;

        bits = 2 + n_randint(state, 200);

        arf_init(x);
        arf_init(z);
        fmpz_init(y);
        fmpz_init(e);

        arf_randtest(x, state, bits, 1 + n_randint(state, 100));
        arf_randtest(z, state, bits, 1 + n_randint(state, 100));

        arf_get_fmpz_2exp(y, e, x);
        arf_set_fmpz_2exp(z, y, e);

        if (!arf_equal(x, z))
        {
            flint_printf("FAIL\n\n");
            flint_printf("bits: %wd\n", bits);
            flint_printf("x = "); arf_print(x); flint_printf("\n\n");
            flint_printf("y = "); fmpz_print(y); flint_printf("\n\n");
            flint_printf("e = "); fmpz_print(e); flint_printf("\n\n");
            flint_printf("z = "); arf_print(z); flint_printf("\n\n");
            flint_abort();
        }

        arf_clear(x);
        arf_clear(z);
        fmpz_clear(y);
        fmpz_clear(e);
    }

    TEST_FUNCTION_END(state);
}
