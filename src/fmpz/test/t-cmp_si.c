/*
    Copyright (C) 2010 Sebastian Pancratz
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "long_extras.h"
#include "fmpz.h"

TEST_FUNCTION_START(fmpz_cmp_si, state)
{
    int i, result;

    /* Compare with fmpz_cmp */
    for (i = 0; i < 10000 * flint_test_multiplier(); i++)
    {
        fmpz_t a, b;
        slong n;
        int lhs, rhs;

        fmpz_init(a);
        fmpz_init(b);

        fmpz_randtest(a, state, 200);

        n = z_randtest(state);
        fmpz_set_si(b, n);

        lhs = fmpz_cmp(a, b);
        rhs = fmpz_cmp_si(a, n);

        result = (lhs < 0) ? (rhs < 0) : ((lhs > 0) ? (rhs > 0) : (rhs == 0));
        if (result == 0)
        {
            flint_printf("FAIL:\n");
            flint_printf("a = "), fmpz_print(a), flint_printf("\n");
            flint_printf("b = "), fmpz_print(b), flint_printf("\n");
            flint_printf("n = %wd\n", n);
            flint_printf("cmp(a, b) = %d\n", fmpz_cmp(a, b));
            flint_printf("cmp_si(a, n) = %d\n", fmpz_cmp_si(a, n));
            fflush(stdout);
            flint_abort();
        }

        fmpz_clear(a);
        fmpz_clear(b);
    }

    TEST_FUNCTION_END(state);
}
