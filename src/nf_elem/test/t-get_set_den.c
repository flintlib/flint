/*
    Copyright (C) 2013 William Hart
                  2020 Julian Rüth

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "nf_elem.h"

TEST_FUNCTION_START(nf_elem_get_set_den, state)
{
    int i, result;

    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        nf_t nf;
        nf_elem_t a;
        fmpz_t d, d2;

        nf_init_randtest(nf, state, 40, 200);

        fmpz_init(d);
        fmpz_init(d2);

        nf_elem_init(a, nf);
        nf_elem_randtest(a, state, 200, nf);

        fmpz_randtest_not_zero(d, state, 200);

        nf_elem_set_den(a, d, nf);
        nf_elem_get_den(d2, a, nf);

        result = fmpz_equal(d, d2);
        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("d = "); fmpz_print(d); printf("\n");
            flint_printf("d2 = "); fmpz_print(d2); printf("\n");
            flint_printf("a = "); nf_elem_print_pretty(a, nf, "x");
            flint_abort();
        }

        nf_elem_clear(a, nf);

        nf_clear(nf);

        fmpz_clear(d);
        fmpz_clear(d2);
    }

    TEST_FUNCTION_END(state);
}
