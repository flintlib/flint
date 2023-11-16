/*
    Copyright (C) 2019 D.H.J Polymath

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "acb_dirichlet.h"

TEST_FUNCTION_START(acb_dirichlet_turing_method_bound, state)
{
    slong iter;

    for (iter = 0; iter < 500 * 0.1 * flint_test_multiplier(); iter++)
    {
        ulong u, v;
        fmpz_t a, b;

        fmpz_init(a);
        fmpz_init(b);

        fmpz_randtest_unsigned(a, state, 500);
        fmpz_randtest_unsigned(b, state, 500);
        fmpz_sub_ui(a, a, 1);
        fmpz_sub_ui(b, b, 1);

        u = acb_dirichlet_turing_method_bound(a);
        v = acb_dirichlet_turing_method_bound(b);

        if ((fmpz_cmp(a, b) < 0 && u > v) ||
            (fmpz_cmp(a, b) > 0 && u < v))
        {
            flint_printf("FAIL: increasing on p >= -1\n\n");
            flint_printf("a = "); fmpz_print(a); flint_printf("\n\n");
            flint_printf("b = "); fmpz_print(b); flint_printf("\n\n");
            flint_printf("u = %wu\n", u);
            flint_printf("v = %wu\n", v);
            flint_abort();
        }

        fmpz_clear(a);
        fmpz_clear(b);
    }

    TEST_FUNCTION_END(state);
}
