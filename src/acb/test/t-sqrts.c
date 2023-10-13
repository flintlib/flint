/*
    Copyright (C) 2023 Jean Kieffer

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "acb.h"

int main(void)
{
    slong iter;
    flint_rand_t state;

    flint_printf("sqrts....");
    fflush(stdout);

    flint_randinit(state);

    /* check sqrt overlaps one of the results, y1 = -y2, no precision loss */
    for (iter = 0; iter < 1000 * flint_test_multiplier(); iter++)
    {
        acb_t x, y1, y2, t;
        arf_t e;
        slong prec = 20 + n_randint(state, 1000);

        acb_init(x);
        acb_init(y1);
        acb_init(y2);
        acb_init(t);
        arf_init(e);

        acb_urandom(x, state, prec);
        acb_sqrt(t, x, prec);
        acb_sqrts(y1, y2, x, prec);

        if (!acb_overlaps(y1, t) && !acb_overlaps(y2, t))
        {
            flint_printf("FAIL (overlap)\n");
            acb_printd(x, 10);
            flint_printf("\n");
            acb_printd(y1, 10);
            flint_printf("\n");
            acb_printd(y2, 10);
            flint_printf("\n");
            flint_abort();
        }

        acb_neg(y2, y2);
        if (!acb_equal(y1, y2))
        {
            flint_printf("FAIL (negative)\n");
            acb_printd(y1, 10);
            flint_printf("\n");
            acb_printd(y2, 10);
            flint_printf("\n");
            flint_abort();
        }

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

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return 0;
}

