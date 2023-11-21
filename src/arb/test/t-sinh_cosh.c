/*
    Copyright (C) 2012 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include <mpfr.h>
#include "fmpq.h"
#include "arb.h"

TEST_FUNCTION_START(arb_sinh_cosh, state)
{
    slong iter;

    for (iter = 0; iter < 100000 * 0.1 * flint_test_multiplier(); iter++)
    {
        arb_t a, b, c;
        fmpq_t q;
        mpfr_t t, u;
        slong prec = 2 + n_randint(state, 200);

        arb_init(a);
        arb_init(b);
        arb_init(c);
        fmpq_init(q);
        mpfr_init2(t, prec + 100);
        mpfr_init2(u, prec + 100);

        arb_randtest(a, state, 1 + n_randint(state, 200), 3);
        arb_randtest(b, state, 1 + n_randint(state, 200), 3);
        arb_randtest(c, state, 1 + n_randint(state, 200), 3);
        arb_get_rand_fmpq(q, state, a, 1 + n_randint(state, 200));

        fmpq_get_mpfr(t, q, MPFR_RNDN);
        mpfr_sinh_cosh(t, u, t, MPFR_RNDN);

        arb_sinh_cosh(b, c, a, prec);

        if (!arb_contains_mpfr(b, t))
        {
            flint_printf("FAIL: containment (sin)\n\n");
            flint_printf("a = "); arb_print(a); flint_printf("\n\n");
            flint_printf("b = "); arb_print(b); flint_printf("\n\n");
            flint_abort();
        }

        if (!arb_contains_mpfr(c, u))
        {
            flint_printf("FAIL: containment (cos)\n\n");
            flint_printf("a = "); arb_print(a); flint_printf("\n\n");
            flint_printf("c = "); arb_print(c); flint_printf("\n\n");
            flint_abort();
        }

        arb_clear(a);
        arb_clear(b);
        arb_clear(c);
        fmpq_clear(q);
        mpfr_clear(t);
        mpfr_clear(u);
    }

    TEST_FUNCTION_END(state);
}
