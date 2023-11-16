/*
    Copyright (C) 2014 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "arf.h"

TEST_FUNCTION_START(arf_set_round_uiui, state)
{
    slong iter;

    for (iter = 0; iter < 1000000 * 0.1 * flint_test_multiplier(); iter++)
    {
        arf_t x, y;
        slong prec, fix1, fix2;
        int ret1, ret2, sgnbit;
        mp_limb_t t[2];
        arf_rnd_t rnd;

        prec = 2 + n_randint(state, 1000);

        arf_init(x);
        arf_init(y);

        arf_randtest_special(x, state, 1 + n_randint(state, 200), 1 + n_randint(state, 100));
        arf_randtest_special(y, state, 1 + n_randint(state, 200), 1 + n_randint(state, 100));

        do {
            t[0] = n_randtest(state);
            t[1] = n_randtest(state);
        } while (t[0] == 0 && t[1] == 0);

        sgnbit = n_randint(state, 2);

        switch (n_randint(state, 10))
        {
            case 0: rnd = ARF_RND_DOWN; break;
            case 1: rnd = ARF_RND_UP; break;
            case 2: rnd = ARF_RND_FLOOR; break;
            case 3: rnd = ARF_RND_CEIL; break;
            default: rnd = ARF_RND_NEAR; break;
        }

        if (t[1] != 0)
        {
            ret1 = _arf_set_round_mpn(x, &fix1, t, 2, sgnbit, prec, rnd);
            fmpz_set_si(ARF_EXPREF(x), 2 * FLINT_BITS + fix1);
        }
        else
        {
            ret1 = _arf_set_round_mpn(x, &fix1, t, 1, sgnbit, prec, rnd);
            fmpz_set_si(ARF_EXPREF(x), FLINT_BITS + fix1);
        }

        ret2 = _arf_set_round_uiui(y, &fix2, t[1], t[0], sgnbit, prec, rnd);
        fmpz_set_si(ARF_EXPREF(y), 2 * FLINT_BITS + fix2);

        if (!arf_equal(x, y) || (ret1 != ret2))
        {
            flint_printf("FAIL\n\n");
            flint_printf("prec = %wd", prec); flint_printf("\n\n");
            flint_printf("hi = %wu, lo = %wu\n\n", t[1], t[0]);
            flint_printf("x = "); arf_print(x); flint_printf("\n\n");
            flint_printf("y = "); arf_print(y); flint_printf("\n\n");
            flint_printf("ret1 = %d, ret2 = %d\n\n", ret1, ret2);
            flint_abort();
        }

        arf_clear(x);
        arf_clear(y);
    }

    TEST_FUNCTION_END(state);
}
