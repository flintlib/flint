/*
    Copyright (C) 2023 Jean Kieffer

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "acb.h"

TEST_FUNCTION_START(acb_sqrts, state)
{
    slong iter;

    /* Test: - acb_sqrts on y = x^2 gives both x and -x
       - acb_sqrts on a precise number does not lose precision */
    for (iter = 0; iter < 1000 * flint_test_multiplier(); iter++)
    {
        acb_t x, y1, y2, t;
        arf_t e;
        slong prec = 20 + n_randint(state, 1000);
        slong mag_bits = n_randint(state, 10);

        acb_init(x);
        acb_init(y1);
        acb_init(y2);
        acb_init(t);
        arf_init(e);

        acb_randtest(x, state, prec, mag_bits);
        acb_sqr(y1, x, prec);
        acb_sqrts(y1, y2, y1, prec);
        acb_neg(t, x);

        if (!(acb_contains(y1, x) && acb_contains(y2, t))
            && !(acb_contains(y1, t) && acb_contains(y2, x)))
        {
            flint_printf("FAIL (containment)\n");
            acb_printd(x, 10);
            flint_printf("\n");
            acb_printd(y1, 10);
            flint_printf("\n");
            acb_printd(y2, 10);
            flint_printf("\n");
            flint_abort();
        }

        acb_urandom(x, state, prec);
        acb_sqrts(y1, y2, x, prec);

        arf_one(e);
        arf_mul_2exp_si(e, e, -prec / 2 + 10);

        acb_get_mid(t, y1);
        acb_add_error_arf(t, e);
        if (!acb_contains(t, y1))
        {
            flint_printf("FAIL (precision)\n");
            flint_printf("prec = %wd, y1, x:\n", prec);
            acb_printd(y1, 10);
            flint_printf("\n");
            acb_printd(x, 10);
            flint_printf("\n");
            flint_abort();
        }

        acb_clear(x);
        acb_clear(y1);
        acb_clear(y2);
        acb_clear(t);
        arf_clear(e);
    }

    TEST_FUNCTION_END(state);
}
