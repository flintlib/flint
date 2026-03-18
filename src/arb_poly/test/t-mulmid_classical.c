/*
    Copyright (C) 2012 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpq_poly.h"
#include "arb_poly.h"

TEST_FUNCTION_START(arb_poly_mulmid_classical, state)
{
    slong iter;

    for (iter = 0; iter < 1000 * 0.1 * flint_test_multiplier(); iter++)
    {
        slong rbits1, rbits2, rbits3, nlo, nhi;
        arb_poly_t a, b, c, d;
        int aliasing = n_randint(state, 5);

        rbits1 = 2 + n_randint(state, 300);
        rbits2 = 2 + n_randint(state, 300);
        rbits3 = 2 + n_randint(state, 300);

        nlo = n_randint(state, 10);
        nhi = n_randint(state, 10);

        arb_poly_init(a);
        arb_poly_init(b);
        arb_poly_init(c);
        arb_poly_init(d);

        arb_poly_randtest(a, state, 1 + n_randint(state, 10), rbits1, 1 + n_randint(state, 100));
        arb_poly_randtest(b, state, 1 + n_randint(state, 10), rbits2, 1 + n_randint(state, 100));
        arb_poly_randtest(c, state, 1 + n_randint(state, 10), rbits2, 1 + n_randint(state, 100));

        if (aliasing == 0)
        {
            arb_poly_mulmid_classical(c, a, b, nlo, nhi, rbits3);
        }
        else if (aliasing == 1)
        {
            arb_poly_set(c, a);
            arb_poly_mulmid_classical(c, c, b, nlo, nhi, rbits3);
        }
        else if (aliasing == 2)
        {
            arb_poly_set(c, b);
            arb_poly_mulmid_classical(c, a, c, nlo, nhi, rbits3);
        }
        else if (aliasing == 3)
        {
            arb_poly_set(b, a);
            arb_poly_mulmid_classical(c, a, a, nlo, nhi, rbits3);
        }
        else
        {
            arb_poly_set(b, a);
            arb_poly_set(c, a);
            arb_poly_mulmid_classical(c, c, c, nlo, nhi, rbits3);
        }

        if (nlo == 0)
            arb_poly_mullow_classical(d, a, b, nhi, rbits3);
        else
            arb_poly_mullow_block(d, a, b, nhi, rbits3);
        arb_poly_shift_right(d, d, nlo);

        if (!arb_poly_overlaps(c, d))
        {
            flint_printf("FAIL: arb_poly_mulmid_classical\n\n");
            flint_printf("bits3 = %wd\n", rbits3);
            flint_printf("nlo = %wd, nhi = %wd\n", nlo, nhi);

            flint_printf("a = "); arb_poly_printd(a, 15); flint_printf("\n\n");
            flint_printf("b = "); arb_poly_printd(b, 15); flint_printf("\n\n");
            flint_printf("c = "); arb_poly_printd(c, 15); flint_printf("\n\n");
            flint_printf("d = "); arb_poly_printd(c, 15); flint_printf("\n\n");

            flint_abort();
        }

        arb_poly_clear(a);
        arb_poly_clear(b);
        arb_poly_clear(c);
        arb_poly_clear(d);
    }

    TEST_FUNCTION_END(state);
}
