/*
    Copyright (C) 2009 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpq_poly.h"

TEST_FUNCTION_START(fmpq_poly_mullow, state)
{
    int i, result;
    ulong cflags = UWORD(0);

    /* Compare with truncated product of a and b */
    for (i = 0; i < 200 * flint_test_multiplier(); i++)
    {
        fmpq_poly_t a, b, c;
        slong trunc;

        fmpq_poly_init(a);
        fmpq_poly_init(b);
        fmpq_poly_init(c);
        trunc = n_randint(state, 50);
        fmpq_poly_randtest(b, state, trunc, 200);
        fmpq_poly_randtest(c, state, trunc, 200);

        fmpq_poly_mullow(a, b, c, trunc);
        fmpq_poly_mul(b, b, c);
        fmpq_poly_truncate(b, trunc);

        cflags |= fmpq_poly_is_canonical(a) ? 0 : 1;
        cflags |= fmpq_poly_is_canonical(b) ? 0 : 2;
        result = (fmpq_poly_equal(a, b) && !cflags);
        if (!result)
        {
            flint_printf("FAIL:\n");
            fmpq_poly_debug(a), flint_printf("\n\n");
            fmpq_poly_debug(b), flint_printf("\n\n");
            flint_printf("cflags = %wu\n\n", cflags);
            fflush(stdout);
            flint_abort();
        }

        fmpq_poly_clear(a);
        fmpq_poly_clear(b);
        fmpq_poly_clear(c);
    }

    TEST_FUNCTION_END(state);
}
