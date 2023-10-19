/*
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "long_extras.h"
#include "fmpq.h"

TEST_FUNCTION_START(fmpq_set_si, state)
{
    int i;

    for (i = 0; i < 10000; i++)
    {
        fmpq_t x, y;
        fmpz_t p, q;
        slong P, Q;

        fmpq_init(x);
        fmpq_init(y);
        fmpz_init(p);
        fmpz_init(q);

        P = z_randtest(state);
        Q = n_randtest_not_zero(state);

        fmpz_set_si(p, P);
        fmpz_set_ui(q, Q);

        fmpq_set_fmpz_frac(x, p, q);
        fmpq_set_si(y, P, Q);

        if (!fmpq_is_canonical(y) || !fmpq_equal(x, y))
        {
            flint_printf("FAIL");
            flint_printf("p: "); fmpz_print(p); flint_printf("\n");
            flint_printf("q: "); fmpz_print(q); flint_printf("\n");
            flint_printf("x: "); fmpq_print(x); flint_printf("\n");
            flint_printf("y: "); fmpq_print(y); flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }

        fmpq_clear(x);
        fmpq_clear(y);
        fmpz_clear(p);
        fmpz_clear(q);
    }

    TEST_FUNCTION_END(state);
}
