/*
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpq.h"

TEST_FUNCTION_START(fmpq_add, state)
{
    int i, result;

    for (i = 0; i < 10000 * flint_test_multiplier(); i++)
    {
        fmpq_t a, b, c;
        mpq_t d, e, f, g;
        int aliasing;

        fmpq_init(a);
        fmpq_init(b);
        fmpq_init(c);

        mpq_init(d);
        mpq_init(e);
        mpq_init(f);
        mpq_init(g);

        fmpq_randtest(a, state, 200);
        fmpq_randtest(b, state, 200);

        fmpq_get_mpq(d, a);
        fmpq_get_mpq(e, b);

        aliasing = n_randint(state, 4);

        if (aliasing == 0)
        {
            fmpq_add(c, a, b);
        }
        else if (aliasing == 1)
        {
            fmpq_set(a, b);
            mpq_set(d, e);
            fmpq_add(c, a, a);
        }
        else if (aliasing == 2)
        {
            fmpq_set(c, a);
            fmpq_add(c, c, b);
        }
        else
        {
            fmpq_set(c, b);
            fmpq_add(c, a, c);
        }

        mpq_add(f, d, e);

        fmpq_get_mpq(g, c);

        result = (mpq_cmp(f, g) == 0) && fmpq_is_canonical(c);
        if (!result)
        {
            flint_printf("FAIL:\n");
            gmp_printf("d = %Qd, e = %Qd, f = %Qd, g = %Qd\n", d, e, f, g);
            fflush(stdout);
            flint_abort();
        }

        fmpq_clear(a);
        fmpq_clear(b);
        fmpq_clear(c);

        mpq_clear(d);
        mpq_clear(e);
        mpq_clear(f);
        mpq_clear(g);
    }

    TEST_FUNCTION_END(state);
}
