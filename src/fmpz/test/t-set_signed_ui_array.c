/*
    Copyright (C) 2021 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "ulong_extras.h"
#include "fmpz.h"

TEST_FUNCTION_START(fmpz_set_signed_ui_array, state)
{
    int i;
    slong max_limbs = 20;
    ulong * limbs;
    fmpz_t a, b, c;

    fmpz_init(a);
    fmpz_init(b);
    fmpz_init(c);
    limbs = (ulong *) flint_malloc(max_limbs*sizeof(ulong));

    for (i = 0; i < 100000 * flint_test_multiplier(); i++)
    {
        slong j, n;

        n = n_randint(state, max_limbs) + 1;

        for (j = 0; j < n; j++)
        {
            limbs[j] = n_randlimb(state);

            if (n_randint(state, 10) == 0)
                limbs[j] = 0;

            if (n_randint(state, 10) == 0)
                limbs[j] = -UWORD(1);
        }

        fmpz_set_ui_array(a, limbs, n);
        fmpz_set_signed_ui_array(b, limbs, n);

        fmpz_sub(a, a, b);

        fmpz_one(c);
        fmpz_mul_2exp(c, c, n*FLINT_BITS);

        if (!fmpz_divisible(a, c) || !_fmpz_is_canonical(b))
        {
            flint_printf("FAIL: check answer mod 2^(n*FLINT_BITS)\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_one(c);
        fmpz_mul_2exp(c, c, n*FLINT_BITS - 1);

        if (fmpz_cmp(b, c) >= 0)
        {
            flint_printf("FAIL: check answer < 2^(n*FLINT_BITS - 1)\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_neg(c, c);
        if (fmpz_cmp(b, c) < 0)
        {
            flint_printf("FAIL: check answer >= -2^(n*FLINT_BITS - 1)\n");
            fflush(stdout);
            flint_abort();
        }
    }

    fmpz_clear(a);
    fmpz_clear(b);
    fmpz_clear(c);
    flint_free(limbs);

    TEST_FUNCTION_END(state);
}
