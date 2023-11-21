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

TEST_FUNCTION_START(arf_set_round_ui, state)
{
    slong iter;

    for (iter = 0; iter < 1000000 * 0.1 * flint_test_multiplier(); iter++)
    {
        arf_t x, y;
        slong prec, fix1;
        int ret1, ret2;
        mp_limb_t t;
        arf_rnd_t rnd;

        prec = 2 + n_randint(state, 1000);

        arf_init(x);
        arf_init(y);

        arf_randtest_special(x, state, 1 + n_randint(state, 200), 1 + n_randint(state, 100));
        arf_randtest_special(y, state, 1 + n_randint(state, 200), 1 + n_randint(state, 100));

        t = n_randtest(state);

        switch (n_randint(state, 5))
        {
            case 0: rnd = ARF_RND_DOWN; break;
            case 1: rnd = ARF_RND_UP; break;
            case 2: rnd = ARF_RND_FLOOR; break;
            case 3: rnd = ARF_RND_CEIL; break;
            default: rnd = ARF_RND_NEAR; break;
        }

        if (t == 0)
        {
            arf_zero(x);
            ret1 = 0;
        }
        else
        {
            ret1 = _arf_set_round_mpn(x, &fix1, &t, 1, 0, prec, rnd);
            fmpz_set_si(ARF_EXPREF(x), FLINT_BITS + fix1);
        }

        ret2 = arf_set_round_ui(y, t, prec, rnd);

        if (!arf_equal(x, y) || (ret1 != ret2) || (ret2 == 0 && !arf_equal_ui(y, t)))
        {
            flint_printf("FAIL\n\n");
            flint_printf("prec = %wd", prec); flint_printf("\n\n");
            flint_printf("t = %wu\n\n", t);
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
