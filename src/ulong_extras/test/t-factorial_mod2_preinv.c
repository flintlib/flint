/*
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "gmpcompat.h"
#include "ulong_extras.h"

#ifndef n_factorial_mod2_foolproof
#define n_factorial_mod2_foolproof n_factorial_mod2_foolproof
static mp_limb_t
n_factorial_mod2_foolproof(ulong n, mp_limb_t p, mp_limb_t pinv)
{
    mp_limb_t prod = UWORD(1) % p;

    while (n)
    {
        prod = n_mulmod2_preinv(prod, n, p, pinv);
        n--;
    }

    return prod;
}
#endif

TEST_FUNCTION_START(n_factorial_mod2_preinv, state)
{
    int ix;

    /* n is small (n < 1000) */
    for (ix = 0; ix < 100 * flint_test_multiplier(); ix++)
    {
        mp_limb_t n, p, pinv, x, y;

        n = n_randint(state, 1000);
        p = n_randtest_not_zero(state);
        pinv = n_preinvert_limb(p);

        x = n_factorial_mod2_preinv(n, p, pinv);
        y = n_factorial_mod2_foolproof(n, p, pinv);

        if (x != y)
        {
            flint_printf("FAIL:\n"
                         "n = %wu\n"
                         "p = %wu\n"
                         "x = %wu\n"
                         "y = %wu\n",
                         n, p, x, y);
            flint_abort();
        }
    }

    /* FIXME: Fix some fast test that tests this for random big n */
    /* n is big */
#if FLINT64
    {
        mp_limb_t n, p, pinv, x, y;

        n = UWORD(1000003);
        p = UWORD(187263871632876172);
        pinv = n_preinvert_limb(p);

        x = n_factorial_mod2_preinv(n, p, pinv);
        y = UWORD(31386460113503852);

        if (x != y)
        {
            flint_printf("FAIL:\n"
                         "n = %wu\n"
                         "p = %wu\n"
                         "x = %wu\n"
                         "y = %wu\n",
                         n, p, x, y);
            flint_abort();
        }
    }
#endif

    TEST_FUNCTION_END(state);
}
