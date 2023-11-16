/*
    Copyright (C) 2022 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifdef __GNUC__
# define fabs __builtin_fabs
#else
# include <math.h>
#endif

#include "test_helpers.h"
#include "arb.h"

TEST_FUNCTION_START(arb_atan_frac_bsplit, state)
{
    slong iter;

    for (iter = 0; iter < 1000 * 0.1 * flint_test_multiplier(); iter++)
    {
        arb_t s, t;
        fmpz_t p, q;
        slong prec;
        int hyperbolic;

        arb_init(s);
        arb_init(t);
        fmpz_init(p);
        fmpz_init(q);

        prec = 2 + n_randint(state, 800);
        hyperbolic = n_randint(state, 2);

        fmpz_randtest(p, state, 100);
        fmpz_randtest_not_zero(q, state, 100);

        arb_atan_frac_bsplit(s, p, q, hyperbolic, prec);

        if (arb_is_finite(s))
        {
            arb_set_fmpz(t, p);
            arb_div_fmpz(t, t, q, prec);

            if (hyperbolic)
                arb_atanh(t, t, prec);
            else
                arb_atan(t, t, prec);

            if (!arb_overlaps(s, t))
            {
                flint_printf("FAIL: overlap\n\n");
                flint_printf("s = "); arb_printd(s, 100); flint_printf("\n\n");
                flint_printf("t = "); arb_printd(t, 100); flint_printf("\n\n");
                flint_abort();
            }
        }

        if (fabs(fmpz_get_d(p) / fmpz_get_d(q)) < 0.75 && arb_rel_accuracy_bits(s) < prec - 2)
        {
            flint_printf("FAIL: accuracy\n\n");
            flint_printf("s = "); arb_printd(s, 100); flint_printf("\n\n");
            flint_printf("%wd, %wd\n", prec, arb_rel_accuracy_bits(s));
            flint_abort();
        }

        arb_clear(s);
        arb_clear(t);
        fmpz_clear(p);
        fmpz_clear(q);
    }

    TEST_FUNCTION_END(state);
}
