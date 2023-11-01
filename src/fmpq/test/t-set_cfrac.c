/*
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz_vec.h"
#include "fmpq.h"

TEST_FUNCTION_START(fmpq_set_cfrac, state)
{
    int i;

    for (i = 0; i < 10000; i++)
    {
        fmpq_t x, y, r;
        fmpz * c;
        slong n, bound;

        fmpq_init(x);
        fmpq_init(y);
        fmpq_init(r);

        fmpq_randtest(x, state, 1 + n_randint(state, 1000));
        bound = fmpq_cfrac_bound(x);

        c = _fmpz_vec_init(bound);

        n = fmpq_get_cfrac(c, r, x, bound);
        fmpq_set_cfrac(y, c, n);

        if (!fmpq_equal(x, y))
        {
            flint_printf("FAIL: x != y\n");
            flint_printf("x = "); fmpq_print(x); flint_printf("\n");
            flint_printf("y = "); fmpq_print(y); flint_printf("\n");
            flint_printf("c = "); _fmpz_vec_print(c, n); flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        _fmpz_vec_clear(c, bound);
        fmpq_clear(x);
        fmpq_clear(y);
        fmpq_clear(r);
    }

    TEST_FUNCTION_END(state);
}
