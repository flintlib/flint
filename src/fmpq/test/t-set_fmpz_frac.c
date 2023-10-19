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

TEST_FUNCTION_START(fmpq_set_fmpz_frac, state)
{
    int i;

    for (i = 0; i < 10000; i++)
    {
        fmpq_t x, y;
        fmpz_t p, q;

        fmpq_init(x);
        fmpq_init(y);
        fmpz_init(p);
        fmpz_init(q);

        fmpz_randtest(p, state, 100);
        fmpz_randtest_not_zero(q, state, 100);

        fmpq_set_fmpz_frac(x, p, q);

        fmpz_set(fmpq_numref(y), p);
        fmpz_set(fmpq_denref(y), q);
        fmpq_canonicalise(y);

        if (!fmpq_is_canonical(x) || !fmpq_equal(x, y))
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
