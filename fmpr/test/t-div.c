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

    flint_printf("div....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 10000 * arb_test_multiplier(); iter++)
    {
        slong bits, r1, r2;
        fmpr_t x, y, z, w;
        mpfr_t X, Y, Z;

        bits = 2 + n_randint(state, 200);

        fmpr_init(x);
        fmpr_init(y);
        fmpr_init(z);
        fmpr_init(w);

        mpfr_init2(X, bits + 100);
        mpfr_init2(Y, bits + 100);
        mpfr_init2(Z, bits);

        fmpr_randtest_special(x, state, bits + n_randint(state, 100), 10);
        fmpr_randtest_special(y, state, bits + n_randint(state, 100), 10);
        fmpr_randtest_special(z, state, bits + n_randint(state, 100), 10);

        fmpr_get_mpfr(X, x, MPFR_RNDN);
        fmpr_get_mpfr(Y, y, MPFR_RNDN);

        switch (n_randint(state, 4))
        {
            case 0:
                r1 = mpfr_div(Z, X, Y, MPFR_RNDZ);
                r2 = fmpr_div(z, x, y, bits, FMPR_RND_DOWN);
                break;
            case 1:
                r1 = mpfr_div(Z, X, Y, MPFR_RNDA);
                r2 = fmpr_div(z, x, y, bits, FMPR_RND_UP);
                break;
            case 2:
                r1 = mpfr_div(Z, X, Y, MPFR_RNDD);
                r2 = fmpr_div(z, x, y, bits, FMPR_RND_FLOOR);
                break;
            default:
                r1 = mpfr_div(Z, X, Y, MPFR_RNDU);
                r2 = fmpr_div(z, x, y, bits, FMPR_RND_CEIL);
                break;
        }

        /* we use slightly different semantics for special values (?) */
        if (mpfr_zero_p(Y))
            mpfr_set_nan(Z);
        if (mpfr_inf_p(X) && mpfr_zero_p(Y))
            mpfr_set_nan(Z);

        fmpr_set_mpfr(w, Z);

        if (!fmpr_equal(z, w) || ((r1 == 0) != (r2 == FMPR_RESULT_EXACT))
            || !fmpr_check_ulp(z, r2, bits))
        {
            flint_printf("FAIL\n\n");
            flint_printf("bits = %wd\n", bits);
            flint_printf("r1 = %wd, r2 = %wd\n", r1, r2);
            flint_printf("x = "); fmpr_print(x); flint_printf("\n\n");
            flint_printf("y = "); fmpr_print(y); flint_printf("\n\n");
            flint_printf("z = "); fmpr_print(z); flint_printf("\n\n");
            flint_printf("w = "); fmpr_print(w); flint_printf("\n\n");
            flint_abort();
        }

        fmpr_clear(x);
        fmpr_clear(y);
        fmpr_clear(z);
        fmpr_clear(w);

        mpfr_clear(X);
        mpfr_clear(Y);
        mpfr_clear(Z);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}
