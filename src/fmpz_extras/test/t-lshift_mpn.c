/*
    Copyright (C) 2014 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz_extras.h"

TEST_FUNCTION_START(fmpz_lshift_mpn, state)
{
    slong iter;

    for (iter = 0; iter < 100000 * 0.1 * flint_test_multiplier(); iter++)
    {
        fmpz_t a, b, c;
        ulong e;
        mp_limb_t atmp;
        mp_ptr aptr;
        mp_size_t an;
        int asgnbit;

        fmpz_init(a);
        fmpz_init(b);
        fmpz_init(c);

        fmpz_randtest_not_zero(a, state, 1 + n_randint(state, 2000));
        fmpz_randtest(b, state, 1 + n_randint(state, 2000));
        fmpz_randtest(c, state, 1 + n_randint(state, 2000));

        e = n_randint(state, 1000);

        FMPZ_GET_MPN_READONLY(asgnbit, an, aptr, atmp, *a)
        fmpz_lshift_mpn(b, aptr, an, asgnbit, e);

        fmpz_mul_2exp(c, a, e);

        if (!fmpz_equal(b, c))
        {
            flint_printf("FAIL\n");
            fmpz_print(a); flint_printf("\n\n");
            fmpz_print(b); flint_printf("\n\n");
            fmpz_print(c); flint_printf("\n\n");
            flint_abort();
        }

        fmpz_clear(a);
        fmpz_clear(b);
        fmpz_clear(c);
    }

    TEST_FUNCTION_END(state);
}
