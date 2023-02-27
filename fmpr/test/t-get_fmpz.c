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

    flint_printf("get_fmpz....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 10000 * arb_test_multiplier(); iter++)
    {
        slong bits;
        fmpr_t x;
        mpfr_t y;
        fmpz_t z, z2;
        mpz_t w;

        bits = 2 + n_randint(state, 1000);

        fmpr_init(x);
        mpfr_init2(y, bits);
        fmpz_init(z);
        fmpz_init(z2);
        mpz_init(w);

        fmpr_randtest(x, state, bits, 10);

        fmpr_get_mpfr(y, x, MPFR_RNDN);

        switch (n_randint(state, 5))
        {
            case 0:
                fmpr_get_fmpz(z, x, FMPR_RND_FLOOR);
                mpfr_get_z(w, y, MPFR_RNDD);
                break;
            case 1:
                fmpr_get_fmpz(z, x, FMPR_RND_CEIL);
                mpfr_get_z(w, y, MPFR_RNDU);
                break;
            case 2:
                fmpr_get_fmpz(z, x, FMPR_RND_DOWN);
                mpfr_get_z(w, y, MPFR_RNDZ);
                break;
            case 3:
                fmpr_get_fmpz(z, x, FMPR_RND_UP);
                mpfr_get_z(w, y, MPFR_RNDA);
                break;
            default:
                fmpr_get_fmpz(z, x, FMPR_RND_NEAR);
                mpfr_get_z(w, y, MPFR_RNDN);
                break;
        }

        fmpz_set_mpz(z2, w);

        if (!fmpz_equal(z, z2))
        {
            flint_printf("FAIL\n\n");
            flint_printf("x = "); fmpr_print(x); flint_printf("\n\n");
            flint_printf("z = "); fmpz_print(z); flint_printf("\n\n");
            flint_printf("z2 = "); fmpz_print(z2); flint_printf("\n\n");
            flint_abort();
        }

        fmpr_clear(x);
        mpfr_clear(y);
        fmpz_clear(z);
        fmpz_clear(z2);
        mpz_clear(w);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}
