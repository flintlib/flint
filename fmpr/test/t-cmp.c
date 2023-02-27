/*
    Copyright (C) 2012 Fredrik Johansson

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fmpr.h"


int main()
{
    slong iter;
    flint_rand_t state;

    flint_printf("cmp....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 10000 * arb_test_multiplier(); iter++)
    {
        slong bits;
        fmpr_t x, y;
        mpfr_t X, Y;
        int cmp1, cmp2;

        bits = 2 + n_randint(state, 200);

        fmpr_init(x);
        fmpr_init(y);

        mpfr_init2(X, bits);
        mpfr_init2(Y, bits);

        fmpr_randtest_special(x, state, bits, 10);
        fmpr_randtest_special(y, state, bits, 10);

        fmpr_get_mpfr(X, x, MPFR_RNDN);
        fmpr_get_mpfr(Y, y, MPFR_RNDN);

        cmp1 = fmpr_cmp(x, y);
        cmp2 = mpfr_cmp(X, Y);

        if (cmp1 != cmp2)
        {
            flint_printf("FAIL\n\n");
            flint_printf("x = "); fmpr_print(x); flint_printf("\n\n");
            flint_printf("y = "); fmpr_print(y); flint_printf("\n\n");
            flint_printf("cmp1 = %d, cmp2 = %d\n\n", cmp1, cmp2);
            flint_abort();
        }

        fmpr_clear(x);
        fmpr_clear(y);

        mpfr_clear(X);
        mpfr_clear(Y);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}
