/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpq.h"
#include "qqbar.h"

TEST_FUNCTION_START(qqbar_conjugates, state)
{
    slong iter;

    for (iter = 0; iter < 1000 * 0.1 * flint_test_multiplier(); iter++)
    {
        qqbar_t x, y, z;
        qqbar_ptr r;
        fmpq_t s;
        slong i, d;

        qqbar_init(x);
        qqbar_init(y);
        qqbar_init(z);
        fmpq_init(s);

        qqbar_randtest(x, state, 4, 10);
        d = qqbar_degree(x);

        r = _qqbar_vec_init(d);

        qqbar_conjugates(r, x);

        for (i = 0; i < d; i++)
            qqbar_add(y, y, r + i);

        fmpq_set_fmpz_frac(s, QQBAR_COEFFS(x) + d - 1, QQBAR_COEFFS(x) + d);
        fmpq_neg(s, s);

        qqbar_set_fmpq(z, s);

        if (!qqbar_equal(y, z))
        {
            flint_printf("FAIL! %d\n");
            flint_printf("x = "); qqbar_print(x); flint_printf("\n\n");
            flint_printf("y = "); qqbar_print(y); flint_printf("\n\n");
            flint_printf("z = "); qqbar_print(z); flint_printf("\n\n");
            for (i = 0; i < d; i++)
            {
                flint_printf("r%wd = ", i); qqbar_print(r + i); flint_printf("\n\n");
            }
            flint_abort();
        }

        qqbar_clear(x);
        qqbar_clear(y);
        qqbar_clear(z);
        _qqbar_vec_clear(r, d);
        fmpq_clear(s);
    }

    TEST_FUNCTION_END(state);
}
