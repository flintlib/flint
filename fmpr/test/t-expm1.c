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

    flint_printf("expm1....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 10000 * arb_test_multiplier(); iter++)
    {
        slong bits;
        fmpr_t x, z, w;
        mpfr_t X, Z;

        bits = 2 + n_randint(state, 200);

        fmpr_init(x);
        fmpr_init(z);
        fmpr_init(w);

        mpfr_init2(X, bits + 100);
        mpfr_init2(Z, bits);

        fmpr_randtest_special(x, state, bits + n_randint(state, 100), 3);
        fmpr_randtest_special(z, state, bits + n_randint(state, 100), 3);

        fmpr_get_mpfr(X, x, MPFR_RNDN);

        switch (n_randint(state, 4))
        {
            case 0:
                mpfr_expm1(Z, X, MPFR_RNDZ);
                fmpr_expm1(z, x, bits, FMPR_RND_DOWN);
                break;
            case 1:
                mpfr_expm1(Z, X, MPFR_RNDA);
                fmpr_expm1(z, x, bits, FMPR_RND_UP);
                break;
            case 2:
                mpfr_expm1(Z, X, MPFR_RNDD);
                fmpr_expm1(z, x, bits, FMPR_RND_FLOOR);
                break;
            case 3:
                mpfr_expm1(Z, X, MPFR_RNDU);
                fmpr_expm1(z, x, bits, FMPR_RND_CEIL);
                break;
        }

        fmpr_set_mpfr(w, Z);

        if (!fmpr_equal(z, w))
        {
            flint_printf("FAIL\n\n");
            flint_printf("bits = %wd\n", bits);
            flint_printf("x = "); fmpr_print(x); flint_printf("\n\n");
            flint_printf("z = "); fmpr_print(z); flint_printf("\n\n");
            flint_printf("w = "); fmpr_print(w); flint_printf("\n\n");
            flint_abort();
        }

        fmpr_clear(x);
        fmpr_clear(z);
        fmpr_clear(w);

        mpfr_clear(X);
        mpfr_clear(Z);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}
