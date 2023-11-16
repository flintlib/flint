/*
    Copyright (C) 2012 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "arb.h"

TEST_FUNCTION_START(arb_get_fmpz_mid_rad_10exp, state)
{
    slong iter;

    for (iter = 0; iter < 5000 * 0.1 * flint_test_multiplier(); iter++)
    {
        arb_t x, y, t;
        fmpz_t mid, rad, exp;
        slong n, prec;

        arb_init(x);
        arb_init(y);
        arb_init(t);
        fmpz_init(mid);
        fmpz_init(rad);
        fmpz_init(exp);

        arb_randtest_special(x, state, 1 + n_randint(state, 500), 1 + n_randint(state, 500));
        n = 1 + n_randint(state, 500);
        prec = 2 + n_randint(state, 1500);

        arb_get_fmpz_mid_rad_10exp(mid, rad, exp, x, n);

        arf_set_fmpz(arb_midref(y), mid);
        mag_set_fmpz(arb_radref(y), rad);
        arb_set_ui(t, 10);
        arb_pow_fmpz(t, t, exp, prec);
        arb_mul(y, y, t, prec);

        if (arb_is_finite(x) && !arb_is_zero(x) &&
            fmpz_sizeinbase(mid, 10) < n && fmpz_sizeinbase(rad, 10) < n)
        {
            flint_printf("FAIL (too few digits):\n\n");
            flint_printf("x = "); arb_printd(x, 50); flint_printf("\n\n");
            flint_printf("y = "); arb_printd(y, 50); flint_printf("\n\n");
            flint_printf("mid = "); fmpz_print(mid); flint_printf("\n\n");
            flint_printf("rad = "); fmpz_print(rad); flint_printf("\n\n");
            flint_printf("exp = "); fmpz_print(exp); flint_printf("\n\n");
            flint_abort();
        }

        if (arb_is_finite(x) && !arb_contains(y, x))
        {
            flint_printf("FAIL (containment):\n\n");
            flint_printf("x = "); arb_printd(x, 50); flint_printf("\n\n");
            flint_printf("y = "); arb_printd(y, 50); flint_printf("\n\n");
            flint_printf("mid = "); fmpz_print(mid); flint_printf("\n\n");
            flint_printf("rad = "); fmpz_print(rad); flint_printf("\n\n");
            flint_printf("exp = "); fmpz_print(exp); flint_printf("\n\n");
            flint_abort();
        }

        arb_clear(x);
        arb_clear(y);
        arb_clear(t);
        fmpz_clear(mid);
        fmpz_clear(rad);
        fmpz_clear(exp);
    }

    TEST_FUNCTION_END(state);
}
