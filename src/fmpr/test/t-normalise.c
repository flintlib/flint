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

    flint_printf("normalise....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 100000 * arb_test_multiplier(); iter++)
    {
        fmpr_t x, xcopy, y, err_bound, err;
        slong prec, ret1, ret2, bits, ebits;
        fmpr_rnd_t rnd;

        fmpr_init(x);
        fmpr_init(xcopy);
        fmpr_init(y);
        fmpr_init(err_bound);
        fmpr_init(err);

        bits = 2 + n_randint(state, 2000);
        ebits = 2 + n_randint(state, 200);

        fmpz_randtest_not_zero(fmpr_manref(x), state, bits);
        fmpz_randtest(fmpr_expref(x), state, ebits);

        fmpz_set(fmpr_manref(y), fmpr_manref(x));
        fmpz_set(fmpr_expref(y), fmpr_expref(x));

        fmpr_set(xcopy, x);

        switch (n_randint(state, 4))
        {
            case 0: rnd = FMPR_RND_DOWN; break;
            case 1: rnd = FMPR_RND_UP; break;
            case 2: rnd = FMPR_RND_FLOOR; break;
            default: rnd = FMPR_RND_CEIL; break;
        }

        prec = 2 + n_randint(state, 2000);

        ret1 = _fmpr_normalise(fmpr_manref(x), fmpr_expref(x), prec, rnd);
        ret2 = _fmpr_normalise_naive(fmpr_manref(y), fmpr_expref(y), prec, rnd);

        if (!fmpr_equal(x, y) || ret1 != ret2)
        {
            flint_printf("FAIL!\n");
            flint_printf("x = "); fmpr_print(x); flint_printf("\n\n");
            flint_printf("y = "); fmpr_print(y); flint_printf("\n\n");
            flint_printf("ret1 = %wd, ret2 = %wd\n", ret1, ret2);
            flint_abort();
        }

        fmpr_sub(err, x, xcopy, FMPR_PREC_EXACT, FMPR_RND_DOWN);
        fmpr_abs(err, err);
        fmpr_set_error_result(err_bound, x, ret1);

        if (fmpr_cmp(err, err_bound) > 0)
        {
            flint_printf("FAIL (error bound)!\n");
            flint_printf("x (original) = "); fmpr_print(xcopy); flint_printf("\n\n");
            flint_printf("x = "); fmpr_print(x); flint_printf("\n\n");
            flint_printf("error: "); fmpr_print(err); flint_printf("\n\n");
            flint_printf("error bound: "); fmpr_print(err_bound); flint_printf("\n\n");
            flint_printf("ret = %wd\n", ret1);
            flint_abort();
        }

        fmpr_clear(x);
        fmpr_clear(xcopy);
        fmpr_clear(y);
        fmpr_clear(err_bound);
        fmpr_clear(err);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}

