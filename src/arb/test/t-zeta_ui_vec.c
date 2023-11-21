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
#include "arb.h"

TEST_FUNCTION_START(arb_zeta_ui_vec, state)
{
    slong iter;

    for (iter = 0; iter < 100 * 0.1 * flint_test_multiplier(); iter++)
    {
        arb_ptr r;
        ulong n;
        slong i, num;
        mpfr_t s;
        slong prec, accuracy;

        prec = 2 + n_randint(state, 1 << n_randint(state, 13));
        num = 1 + n_randint(state, 20);

        r = _arb_vec_init(num);
        mpfr_init2(s, prec + 100);

        do { n = n_randint(state, 1 << n_randint(state, 10)); } while (n < 2);

        arb_zeta_ui_vec(r, n, num, prec);

        for (i = 0; i < num; i++)
        {
            mpfr_zeta_ui(s, n + i, MPFR_RNDN);

            if (!arb_contains_mpfr(r + i, s))
            {
                flint_printf("FAIL: containment\n\n");
                flint_printf("n = %wu\n\n", n + i);
                flint_printf("r = "); arb_printd(r + i, prec / 3.33); flint_printf("\n\n");
                flint_printf("s = "); mpfr_printf("%.275Rf\n", s); flint_printf("\n\n");
                flint_abort();
            }

            accuracy = arb_rel_accuracy_bits(r + i);

            if (accuracy < prec - 4)
            {
                flint_printf("FAIL: accuracy = %wd, prec = %wd\n\n", accuracy, prec);
                flint_printf("n = %wu\n\n", n + i);
                flint_printf("r = "); arb_printd(r + i, prec / 3.33); flint_printf("\n\n");
                flint_abort();
            }
        }

        _arb_vec_clear(r, num);
        mpfr_clear(s);
    }

    TEST_FUNCTION_END(state);
}
