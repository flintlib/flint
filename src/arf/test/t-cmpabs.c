/*
    Copyright (C) 2012, 2014 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include <mpfr.h>
#include "arf.h"

TEST_FUNCTION_START(arf_cmpabs, state)
{
    slong iter;

    /* compare with fmpz */
    {
        arf_t x, y;
        fmpz_t X, Y;

        arf_init(x);
        arf_init(y);
        fmpz_init(X);
        fmpz_init(Y);

        for (iter = 0; iter < 100000 * 0.1 * flint_test_multiplier(); iter++)
        {
            int cmp1, cmp2;

            fmpz_randtest(X, state, 1 + n_randint(state, 1000));

            switch (n_randint(state, 8))
            {
                case 0:
                    fmpz_neg(Y, X);
                    break;
                case 1:
                    fmpz_set(Y, X);
                    break;
                default:
                    fmpz_randtest(Y, state, 1 + n_randint(state, 1000));
            }

            arf_set_fmpz(x, X);
            arf_set_fmpz(y, Y);

            cmp1 = arf_cmpabs(x, y);

            cmp2 = fmpz_cmpabs(X, Y);
            cmp2 = (cmp2 > 0) - (cmp2 < 0);

            if (cmp1 != cmp2)
            {
                flint_printf("FAIL\n\n");
                flint_printf("x = "); arf_debug(x); flint_printf("\n\n");
                flint_printf("y = "); arf_debug(y); flint_printf("\n\n");
                flint_printf("X = "); fmpz_print(X); flint_printf("\n\n");
                flint_printf("Y = "); fmpz_print(Y); flint_printf("\n\n");
                flint_printf("cmp1 = %d, cmp2 = %d\n\n", cmp1, cmp2);
                flint_abort();
            }
        }

        arf_clear(x);
        arf_clear(y);
        fmpz_clear(X);
        fmpz_clear(Y);
    }

    /* compare with mpfr */
    for (iter = 0; iter < 100000 * 0.1 * flint_test_multiplier(); iter++)
    {
        slong bits;
        arf_t x, y;
        mpfr_t X, Y;
        int cmp1, cmp2;

        bits = 2 + n_randint(state, 200);

        arf_init(x);
        arf_init(y);

        mpfr_init2(X, bits);
        mpfr_init2(Y, bits);

        arf_randtest_special(x, state, bits, 10);
        arf_randtest_special(y, state, bits, 10);

        arf_get_mpfr(X, x, MPFR_RNDN);
        arf_get_mpfr(Y, y, MPFR_RNDN);

        mpfr_abs(X, X, MPFR_RNDN);
        mpfr_abs(Y, Y, MPFR_RNDN);

        cmp1 = arf_cmpabs(x, y);
        cmp2 = mpfr_cmp(X, Y);

        if (cmp1 != cmp2)
        {
            flint_printf("FAIL\n\n");
            flint_printf("x = "); arf_print(x); flint_printf("\n\n");
            flint_printf("y = "); arf_print(y); flint_printf("\n\n");
            flint_printf("cmp1 = %d, cmp2 = %d\n\n", cmp1, cmp2);
            flint_abort();
        }

        arf_clear(x);
        arf_clear(y);

        mpfr_clear(X);
        mpfr_clear(Y);
    }

    TEST_FUNCTION_END(state);
}
