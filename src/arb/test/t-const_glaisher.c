/*
    Copyright (C) 2013 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "arb.h"

TEST_FUNCTION_START(arb_const_glaisher, state)
{
    slong iter;

    for (iter = 0; iter < 250 * 0.1 * flint_test_multiplier(); iter++)
    {
        arb_t r, s, t;
        fmpz_t v;
        slong accuracy, prec;

        prec = 2 + n_randint(state, 2000);

        arb_init(r);
        arb_init(s);
        arb_init(t);
        fmpz_init(v);

        arb_const_glaisher(r, prec);
        arb_const_glaisher(s, prec + 100);

        if (!arb_overlaps(r, s))
        {
            flint_printf("FAIL: containment\n\n");
            flint_printf("prec = %wd\n", prec);
            flint_printf("r = "); arb_printd(r, prec / 3.33); flint_printf("\n\n");
            flint_abort();
        }

        accuracy = arb_rel_accuracy_bits(r);

        if (accuracy < prec - 4)
        {
            flint_printf("FAIL: poor accuracy\n\n");
            flint_printf("prec = %wd\n", prec);
            flint_printf("r = "); arb_printd(r, prec / 3.33); flint_printf("\n\n");
            flint_abort();
        }

        if (n_randint(state, 30) == 0)
        {
            flint_cleanup();
        }

        fmpz_set_str(v, "128242712910062263687534256886979172776768892732500", 10);
        arb_set_fmpz(t, v);
        mag_one(arb_radref(t));
        fmpz_ui_pow_ui(v, 10, 50);
        arb_div_fmpz(t, t, v, 170);

        if (!arb_overlaps(r, t))
        {
            flint_printf("FAIL: reference value\n\n");
            flint_printf("prec = %wd\n", prec);
            flint_printf("r = "); arb_printd(r, prec / 3.33); flint_printf("\n\n");
            flint_abort();
        }

        arb_clear(r);
        arb_clear(s);
        arb_clear(t);
        fmpz_clear(v);
    }

    TEST_FUNCTION_END(state);
}
