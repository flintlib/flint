/*
    Copyright (C) 2012, 2013 Fredrik Johansson

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fmpr.h"

#if MPFR_VERSION_MAJOR >= 4
#define func_mpfr_root mpfr_rootn_ui
#else
#define func_mpfr_root mpfr_root
#endif

int main()
{
    slong iter;
    flint_rand_t state;

    flint_printf("root....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 10000 * arb_test_multiplier(); iter++)
    {
        slong bits;
        ulong k;
        fmpr_t x, z, w;
        mpfr_t X, Z;

        bits = 2 + n_randint(state, 200);
        k = 1 + n_randint(state, 20);

        fmpr_init(x);
        fmpr_init(z);
        fmpr_init(w);

        mpfr_init2(X, k * (bits + 100));
        mpfr_init2(Z, bits);

        fmpr_randtest_special(x, state, bits + n_randint(state, 100), 10);
        fmpr_abs(x, x);
        fmpr_randtest_special(z, state, bits + n_randint(state, 100), 10);

        /* occasionally produce perfect powers */
        if (n_randint(state, 4) == 0)
            fmpr_mul(x, x, x, k * bits, FMPR_RND_DOWN);

        fmpr_get_mpfr(X, x, MPFR_RNDN);

        switch (n_randint(state, 4))
        {
            case 0:
                func_mpfr_root(Z, X, k, MPFR_RNDZ);
                fmpr_root(z, x, k, bits, FMPR_RND_DOWN);
                break;
            case 1:
                func_mpfr_root(Z, X, k, MPFR_RNDA);
                fmpr_root(z, x, k, bits, FMPR_RND_UP);
                break;
            case 2:
                func_mpfr_root(Z, X, k, MPFR_RNDD);
                fmpr_root(z, x, k, bits, FMPR_RND_FLOOR);
                break;
            case 3:
                func_mpfr_root(Z, X, k, MPFR_RNDU);
                fmpr_root(z, x, k, bits, FMPR_RND_CEIL);
                break;
        }

        fmpr_set_mpfr(w, Z);

        if (!fmpr_equal(z, w))
        {
            flint_printf("FAIL\n\n");
            flint_printf("bits = %wd, k = %wu\n", bits, k);
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
