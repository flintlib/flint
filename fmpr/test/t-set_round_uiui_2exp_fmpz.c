/*
    Copyright (C) 2013 Fredrik Johansson

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

    flint_printf("set_round_uiui_2exp_fmpz....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 10000 * arb_test_multiplier(); iter++)
    {
        slong prec, ret1, ret2;
        fmpz_t man, exp;
        fmpr_t x, y;
        mp_limb_t hi, lo;
        fmpr_rnd_t rnd;
        int negative;

        fmpz_init(man);
        fmpz_init(exp);
        fmpr_init(x);
        fmpr_init(y);

        prec = 2 + n_randint(state, 1000);
        if (n_randint(state, 10) == 0)
            prec = FMPR_PREC_EXACT;

        negative = n_randint(state, 2);
        hi = n_randtest(state);
        lo = n_randtest(state);
        fmpz_randtest(exp, state, 1 + n_randint(state, 100));
        if (negative)
            fmpz_neg_uiui(man, hi, lo);
        else
            fmpz_set_uiui(man, hi, lo);

        switch (n_randint(state, 4))
        {
            case 0: rnd = FMPR_RND_DOWN; break;
            case 1: rnd = FMPR_RND_UP; break;
            case 2: rnd = FMPR_RND_FLOOR; break;
            default: rnd = FMPR_RND_CEIL; break;
        }

        ret1 = fmpr_set_round_uiui_2exp_fmpz(x, hi, lo, exp, negative, prec, rnd);

        fmpr_set_fmpz_2exp(y, man, exp);
        ret2 = fmpr_set_round(y, y, prec, rnd);

        if (!fmpr_equal(x, y) || ret1 != ret2 ||
            !fmpr_check_ulp(x, ret1, prec) || !fmpr_check_ulp(y, ret2, prec))
        {
            flint_printf("FAIL\n\n");
            flint_printf("prec: %wd\n", prec);
            flint_printf("hi = %wu, lo = %wu\n", hi, lo);
            flint_printf("man = "); fmpz_print(man); flint_printf("\n\n");
            flint_printf("exp = "); fmpz_print(exp); flint_printf("\n\n");
            flint_printf("x = "); fmpr_print(x); flint_printf("\n\n");
            flint_printf("y = "); fmpr_print(y); flint_printf("\n\n");
            flint_printf("ret1 = %wd, ret2 = %wd\n\n", ret1, ret2);
            flint_abort();
        }

        fmpz_clear(man);
        fmpz_clear(exp);
        fmpr_clear(x);
        fmpr_clear(y);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}

