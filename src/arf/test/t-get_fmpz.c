/*
    Copyright (C) 2012 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include <mpfr.h>
#include "arf.h"

TEST_FUNCTION_START(arf_get_fmpz, state)
{
    slong iter;

    for (iter = 0; iter < 10000 * 0.1 * flint_test_multiplier(); iter++)
    {
        arf_t x, x2;
        fmpz_t z, z2, e;
        int ret1, ret2;

        arf_init(x);
        arf_init(x2);
        fmpz_init(z);
        fmpz_init(z2);
        fmpz_init(e);

        arf_randtest(x, state, 2 + n_randint(state, 1000), 10);
        fmpz_randtest(z, state, 1 + n_randint(state, 1000));
        fmpz_randtest(z2, state, 1 + n_randint(state, 1000));
        fmpz_randtest(e, state, 1 + n_randint(state, 200));
        arf_mul_2exp_fmpz(x2, x, e);

        ret1 = arf_get_fmpz(z, x, ARF_RND_DOWN);
        ret2 = arf_get_fmpz_fixed_fmpz(z2, x2, e);

        if (!fmpz_equal(z, z2) || (ret1 != ret2))
        {
            flint_printf("FAIL (fixed_fmpz)\n\n");
            flint_printf("x = "); arf_print(x); flint_printf("\n\n");
            flint_printf("x2 = "); arf_print(x2); flint_printf("\n\n");
            flint_printf("z = "); fmpz_print(z); flint_printf("\n\n");
            flint_printf("z2 = "); fmpz_print(z2); flint_printf("\n\n");
            flint_printf("ret1 = %d, ret2 = %d\n\n", ret1, ret2);
            flint_abort();
        }

        arf_clear(x);
        arf_clear(x2);
        fmpz_clear(z);
        fmpz_clear(z2);
        fmpz_clear(e);
    }

    for (iter = 0; iter < 10000 * 0.1 * flint_test_multiplier(); iter++)
    {
        arf_t x, x2;
        fmpz_t z, z2;
        slong e;
        int ret1, ret2;

        arf_init(x);
        arf_init(x2);
        fmpz_init(z);
        fmpz_init(z2);

        arf_randtest(x, state, 2 + n_randint(state, 1000), 10);
        fmpz_randtest(z, state, 1 + n_randint(state, 1000));
        fmpz_randtest(z2, state, 1 + n_randint(state, 1000));
        e = n_randtest(state);
        arf_mul_2exp_si(x2, x, e);

        ret1 = arf_get_fmpz(z, x, ARF_RND_DOWN);
        ret2 = arf_get_fmpz_fixed_si(z2, x2, e);

        if (!fmpz_equal(z, z2) || (ret1 != ret2))
        {
            flint_printf("FAIL (fixed_si)\n\n");
            flint_printf("x = "); arf_print(x); flint_printf("\n\n");
            flint_printf("x2 = "); arf_print(x2); flint_printf("\n\n");
            flint_printf("z = "); fmpz_print(z); flint_printf("\n\n");
            flint_printf("z2 = "); fmpz_print(z2); flint_printf("\n\n");
            flint_printf("ret1 = %d, ret2 = %d\n\n", ret1, ret2);
            flint_abort();
        }

        arf_clear(x);
        arf_clear(x2);
        fmpz_clear(z);
        fmpz_clear(z2);
    }

    for (iter = 0; iter < 1000000 * 0.1 * flint_test_multiplier(); iter++)
    {
        slong bits;
        arf_t x;
        mpfr_t y;
        fmpz_t z, z2;
        mpz_t w;
        int ret1, ret2;

        bits = 2 + n_randint(state, 1000);

        arf_init(x);
        mpfr_init2(y, bits);
        fmpz_init(z);
        fmpz_init(z2);
        mpz_init(w);

        arf_randtest(x, state, bits, 10);
        fmpz_randtest(z, state, 1 + n_randint(state, 1000));

        arf_get_mpfr(y, x, MPFR_RNDN);

        switch (n_randint(state, 5))
        {
            case 0:
                ret1 = arf_get_fmpz(z, x, ARF_RND_FLOOR);
                ret2 = mpfr_get_z(w, y, MPFR_RNDD);
                break;
            case 1:
                ret1 = arf_get_fmpz(z, x, ARF_RND_CEIL);
                ret2 = mpfr_get_z(w, y, MPFR_RNDU);
                break;
            case 2:
                ret1 = arf_get_fmpz(z, x, ARF_RND_DOWN);
                ret2 = mpfr_get_z(w, y, MPFR_RNDZ);
                break;
            case 3:
                ret1 = arf_get_fmpz(z, x, ARF_RND_UP);
                ret2 = mpfr_get_z(w, y, MPFR_RNDA);
                break;
            default:
                ret1 = arf_get_fmpz(z, x, ARF_RND_NEAR);
                ret2 = mpfr_get_z(w, y, MPFR_RNDN);
                break;
        }

        fmpz_set_mpz(z2, w);

        if (!fmpz_equal(z, z2) || (ret1 != (ret2 != 0)))
        {
            flint_printf("FAIL\n\n");
            flint_printf("x = "); arf_print(x); flint_printf("\n\n");
            flint_printf("z = "); fmpz_print(z); flint_printf("\n\n");
            flint_printf("z2 = "); fmpz_print(z2); flint_printf("\n\n");
            flint_printf("ret1 = %d, ret2 = %d\n\n", ret1, ret2);
            flint_abort();
        }

        arf_clear(x);
        mpfr_clear(y);
        fmpz_clear(z);
        fmpz_clear(z2);
        mpz_clear(w);
    }

    TEST_FUNCTION_END(state);
}
