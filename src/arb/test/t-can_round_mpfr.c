/*
    Copyright (C) 2016 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <mpfr.h>
#include "test_helpers.h"
#include "arb.h"

TEST_FUNCTION_START(arb_can_round_mpfr, state)
{
    slong iter;

    for (iter = 0; iter < 1000000 * 0.1 * flint_test_multiplier(); iter++)
    {
        mpfr_t x, y1, y2;
        int r1, r2;
        arb_t t;
        slong prec;
        mpfr_rnd_t rnd;

        prec = 2 + n_randint(state, 300);

        mpfr_init2(x, 2 + n_randint(state, 300));
        mpfr_init2(y1, prec);
        mpfr_init2(y2, prec);

        arb_init(t);

        switch (n_randint(state, 5))
        {
            case 0: rnd = MPFR_RNDN; break;
            case 1: rnd = MPFR_RNDZ; break;
            case 2: rnd = MPFR_RNDU; break;
            case 3: rnd = MPFR_RNDD; break;
            default: rnd = MPFR_RNDA;
        }

        arf_randtest(arb_midref(t), state, mpfr_get_prec(x), 1 + n_randint(state, 10));
        arf_abs(arb_midref(t), arb_midref(t));
        arf_get_mpfr(x, arb_midref(t), MPFR_RNDN);
        arb_root_ui(t, t, 4, 2 + n_randint(state, 300));

        if (arb_can_round_mpfr(t, prec, rnd))
        {
            r1 = mpfr_rootn_ui(y1, x, 4, rnd);
            r2 = arf_get_mpfr(y2, arb_midref(t), rnd);

            if (r1 != r2 || !mpfr_equal_p(y1, y2))
            {
                flint_printf("FAIL!\n");
                flint_printf("r1 = %d, r2 = %d, prec = %wd\n", r1, r2, prec);
                flint_printf("x = "); mpfr_dump(x); flint_printf("\n");
                flint_printf("y1 = "); mpfr_dump(y1); flint_printf("\n");
                flint_printf("y2 = "); mpfr_dump(y2); flint_printf("\n");
                flint_abort();
            }
        }

        arb_clear(t);
        mpfr_clear(x);
        mpfr_clear(y1);
        mpfr_clear(y2);
    }

    TEST_FUNCTION_END(state);
}
