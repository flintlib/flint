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

    flint_printf("set_round_mpn....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 100000 * arb_test_multiplier(); iter++)
    {
        slong prec, bits, shift, ret1, ret2;
        fmpz_t man, exp;
        fmpr_t y;
        mpz_t x;
        fmpz_t fx;
        fmpr_rnd_t rnd;

        fmpz_init(man);
        fmpz_init(exp);
        mpz_init(x);
        fmpz_init(fx);
        fmpr_init(y);

        bits = 1 + n_randint(state, 1000);
        prec = 2 + n_randint(state, 1000);
        if (n_randint(state, 10) == 0)
            prec = FMPR_PREC_EXACT;

        fmpz_randtest_not_zero(fx, state, bits);
        fmpz_get_mpz(x, fx);

        switch (n_randint(state, 4))
        {
            case 0: rnd = FMPR_RND_DOWN; break;
            case 1: rnd = FMPR_RND_UP; break;
            case 2: rnd = FMPR_RND_FLOOR; break;
            default: rnd = FMPR_RND_CEIL; break;
        }

        ret1 = _fmpr_set_round_mpn(&shift, man, x->_mp_d,
            FLINT_ABS(x->_mp_size), (x->_mp_size < 0) ? 1 : 0, prec, rnd);
        fmpz_set_si(exp, shift);

        ret2 = fmpr_set_round_fmpz(y, fx, prec, rnd);

        if (!fmpz_equal(fmpr_manref(y), man) ||
            !fmpz_equal(fmpr_expref(y), exp) || ret1 != ret2 ||
            !fmpr_check_ulp(y, ret2, prec))
        {
            flint_printf("FAIL\n\n");
            flint_printf("bits: %wd\n", bits);
            flint_printf("prec: %wd\n", prec);
            flint_printf("x = "); fmpz_print(fx); flint_printf("\n\n");
            flint_printf("man = "); fmpz_print(man); flint_printf("\n\n");
            flint_printf("exp = "); fmpz_print(exp); flint_printf("\n\n");
            flint_printf("y = "); fmpr_print(y); flint_printf("\n\n");
            flint_printf("ret1 = %wd, ret2 = %wd\n\n", ret1, ret2);
            flint_abort();
        }

        fmpz_clear(man);
        fmpz_clear(exp);
        mpz_clear(x);
        fmpz_clear(fx);
        fmpr_clear(y);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}
