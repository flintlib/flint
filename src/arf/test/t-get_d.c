/*
    Copyright (C) 2012 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "mpfr.h"
#include "arf.h"

TEST_FUNCTION_START(arf_get_d, state)
{
    slong iter;

    /* test exact roundtrip */
    for (iter = 0; iter < 100000 * 0.1 * flint_test_multiplier(); iter++)
    {
        arf_t x, z;
        double y;
        arf_rnd_t rnd;

        arf_init(x);
        arf_init(z);

        switch (n_randint(state, 4))
        {
            case 0:  rnd = ARF_RND_DOWN; break;
            case 1:  rnd = ARF_RND_UP; break;
            case 2:  rnd = ARF_RND_FLOOR; break;
            case 3:  rnd = ARF_RND_CEIL; break;
            default: rnd = ARF_RND_NEAR; break;
        }

        arf_randtest_special(x, state, 53, 8);
        y = arf_get_d(x, rnd);
        arf_set_d(z, y);

        if (!arf_equal(x, z))
        {
            flint_printf("FAIL:\n\n");
            flint_printf("x = "); arf_print(x); flint_printf("\n\n");
            flint_printf("y = %.17g\n\n", y);
            flint_printf("z = "); arf_print(z); flint_printf("\n\n");
            flint_abort();
        }

        arf_clear(x);
        arf_clear(z);
    }

    /* test rounding */
    for (iter = 0; iter < 100000 * 0.1 * flint_test_multiplier(); iter++)
    {
        arf_t x, z, w;
        arf_rnd_t rnd;
        double y;

        arf_init(x);
        arf_init(z);
        arf_init(w);

        arf_randtest_special(x, state, 300, 8);

        switch (n_randint(state, 4))
        {
            case 0:  rnd = ARF_RND_DOWN; break;
            case 1:  rnd = ARF_RND_UP; break;
            case 2:  rnd = ARF_RND_FLOOR; break;
            default: rnd = ARF_RND_CEIL; break;
        }

        y = arf_get_d(x, rnd);
        arf_set_d(w, y);

        arf_set_round(z, x, 53, rnd);

        if (!arf_equal(w, z))
        {
            flint_printf("FAIL:\n\n");
            flint_printf("x = "); arf_print(x); flint_printf("\n\n");
            flint_printf("y = %.17g\n\n", y);
            flint_printf("z = "); arf_print(z); flint_printf("\n\n");
            flint_printf("w = "); arf_print(w); flint_printf("\n\n");
            flint_abort();
        }

        arf_clear(x);
        arf_clear(z);
        arf_clear(w);
    }

    /* compare with mpfr */
    for (iter = 0; iter < 100000 * 0.1 * flint_test_multiplier(); iter++)
    {
        arf_t x, r1, r2;
        arf_rnd_t rnd;
        mpfr_t t;
        double d1, d2;

        arf_init(x);
        arf_init(r1);
        arf_init(r2);
        mpfr_init2(t, 300);

        arf_randtest_special(x, state, 300, 20);
        arf_get_mpfr(t, x, MPFR_RNDD);

        switch (n_randint(state, 4))
        {
            case 0:  rnd = ARF_RND_DOWN; break;
            case 1:  rnd = ARF_RND_UP; break;
            case 2:  rnd = ARF_RND_FLOOR; break;
            case 3:  rnd = ARF_RND_CEIL; break;
            default: rnd = ARF_RND_NEAR; break;
        }

        d1 = arf_get_d(x, rnd);
        d2 = mpfr_get_d(t, arf_rnd_to_mpfr(rnd));

        arf_set_d(r1, d1);
        arf_set_d(r2, d2);

        if (!arf_equal(r1, r2))
        {
            flint_printf("FAIL:\n\n");
            flint_printf("rnd = %i\n\n", rnd);
            flint_printf("x = "); arf_print(x); flint_printf("\n\n");
            flint_printf("d1 = %.17g\n\n", d1);
            flint_printf("d2 = %.17g\n\n", d2);
            flint_printf("r1 = "); arf_print(r1); flint_printf("\n\n");
            flint_printf("r2 = "); arf_print(r2); flint_printf("\n\n");
            flint_abort();
        }

        arf_clear(x);
        arf_clear(r1);
        arf_clear(r2);
        mpfr_clear(t);
    }

    TEST_FUNCTION_END(state);
}
