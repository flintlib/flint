/*
    Copyright (C) 2009 William Hart
    Copyright (C) 2013 Tom Bachmann (adapt for fmpq)

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpq.h"

TEST_FUNCTION_START(fmpq_abs, state)
{
    int i, result;

    for (i = 0; i < 10000 * flint_test_multiplier(); i++)
    {
        fmpq_t a, b;
        mpq_t c, d;

        fmpq_init(a);
        fmpq_init(b);
        mpq_init(c);
        mpq_init(d);

        fmpq_randtest(a, state, 200);

        fmpq_get_mpq(c, a);

        if (n_randint(state, 2)) /* test aliasing */
        {
            fmpq_set(b, a);
            fmpq_abs(b, b);
        }
        else
        {
            fmpq_abs(b, a);
        }

        mpq_abs(c, c);

        fmpq_get_mpq(d, b);

        result = (mpq_cmp(c, d) == 0) && fmpq_is_canonical(b);
        if (!result)
        {
            flint_printf("FAIL:\n");
            gmp_printf("c = %Qd, d = %Qd\n", c, d);
            fflush(stdout);
            flint_abort();
        }

        fmpq_clear(a);
        fmpq_clear(b);
        mpq_clear(c);
        mpq_clear(d);
    }

    TEST_FUNCTION_END(state);
}
