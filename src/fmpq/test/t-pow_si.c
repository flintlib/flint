/*
    Copyright (C) 2009 William Hart
    Copyright (C) 2011 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "long_extras.h"
#include "fmpq.h"

TEST_FUNCTION_START(fmpq_pow_si, state)
{
    int i, result;

    /* Check aliasing of a and b */
    for (i = 0; i < 100000; i++)
    {
        fmpq_t a, b, c;
        slong e;

        fmpq_init(a);
        fmpq_init(b);
        fmpq_init(c);

        fmpq_randtest(a, state, 100);
        fmpq_set(b, a);

        e = z_randint(state, 20);

        if (fmpq_is_zero(b) && e < 0)
            e = -e;

        fmpq_pow_si(c, b, e);
        fmpq_pow_si(b, b, e);

        result = fmpq_equal(b, c);
        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("a = "), fmpq_print(a), flint_printf("\n\n");
            flint_printf("b = "), fmpq_print(b), flint_printf("\n\n");
            flint_printf("c = "), fmpq_print(c), flint_printf("\n\n");
            flint_printf("e = %wd\n", e);
            fflush(stdout);
            flint_abort();
        }

        fmpq_clear(a);
        fmpq_clear(b);
        fmpq_clear(c);
    }

    /* Compare with repeated multiplication */
    for (i = 0; i < 10000; i++)
    {
        fmpq_t a, b, c;
        slong j, e;

        fmpq_init(a);
        fmpq_init(b);
        fmpq_init(c);

        fmpq_randtest(a, state, 50);

        e = z_randint(state, 20);

        if (fmpq_is_zero(a) &&  e < 0)
            e = -e;

        fmpq_pow_si(b, a, e);

        fmpq_one(c);
        for (j = 0; j < FLINT_ABS(e); j++)
            fmpq_mul(c, c, a);
        if (e < 0)
            fmpq_inv(c, c);

        result = fmpq_equal(b, c);
        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("a = "), fmpq_print(a), flint_printf("\n\n");
            flint_printf("b = "), fmpq_print(b), flint_printf("\n\n");
            flint_printf("c = "), fmpq_print(c), flint_printf("\n\n");
            flint_printf("e = %wd\n", e);
            fflush(stdout);
            flint_abort();
        }

        fmpq_clear(a);
        fmpq_clear(b);
        fmpq_clear(c);
    }

    TEST_FUNCTION_END(state);
}
