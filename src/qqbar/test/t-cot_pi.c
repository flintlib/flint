/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpq.h"
#include "qqbar.h"

TEST_FUNCTION_START(qqbar_cot_pi, state)
{
    slong iter;

    for (iter = 0; iter < 200 * 0.1 * flint_test_multiplier(); iter++)
    {
        qqbar_t x;
        arb_t z, w;
        fmpq_t t;
        slong p;
        ulong q;
        slong prec;

        qqbar_init(x);
        arb_init(z);
        arb_init(w);
        fmpq_init(t);

        q = 1 + n_randint(state, 20);
        p = n_randint(state, 1000);
        p -= 500;
        prec = 2 + n_randint(state, 1000);

        fmpq_set_si(t, p, q);
        arb_sin_pi_fmpq(w, t, prec);
        arb_cos_pi_fmpq(z, t, prec);
        arb_div(w, z, w, prec);

        if (arb_is_finite(w))
        {
            qqbar_cot_pi(x, p, q);
            qqbar_get_arb(z, x, prec);

            if (!arb_overlaps(z, w))
            {
                flint_printf("FAIL!\n");
                flint_printf("x = "); qqbar_print(x); flint_printf("\n\n");
                flint_printf("z = "); arb_printn(z, 200, 0); flint_printf("\n\n");
                flint_printf("w = "); arb_printn(w, 200, 0); flint_printf("\n\n");
                flint_printf("p, q = %wd %wu\n\n", p, q);
                flint_abort();
            }
        }

        qqbar_clear(x);
        arb_clear(z);
        arb_clear(w);
        fmpq_clear(t);
    }

    TEST_FUNCTION_END(state);
}
