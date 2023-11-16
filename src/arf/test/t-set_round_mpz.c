/*
    Copyright (C) 2014 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include <mpfr.h>
#include "arf.h"

TEST_FUNCTION_START(arf_set_round_mpz, state)
{
    slong iter;

    {
        arf_t x, y;

        arf_init(x);
        arf_init(y);

        for (iter = 0; iter < 100000 * 0.1 * flint_test_multiplier(); iter++)
        {
            slong bits1, bits2;
            int ret1, ret2;
            fmpz_t a;
            mpz_t b;
            mpfr_t g;
            arf_rnd_t rnd;

            bits1 = 1 + n_randint(state, 1000);
            bits2 = 2 + n_randint(state, 1000);

            if (n_randint(state, 100) == 0)
                bits2 = ARF_PREC_EXACT;

            fmpz_init(a);
            mpz_init(b);
            mpfr_init2(g, FLINT_MIN(bits2, 10000));

            if (n_randint(state, 100) == 0)
            {
                arf_clear(x);
                arf_clear(y);
                arf_init(x);
                arf_init(y);
            }

            /* dirty output variables */
            if (n_randint(state, 2))
            {
                arf_randtest_special(x, state, 1 + n_randint(state, 1000),
                    1 + n_randint(state, 100));
                arf_randtest_special(y, state, 1 + n_randint(state, 1000),
                    1 + n_randint(state, 100));
            }

            fmpz_randtest(a, state, bits1);
            fmpz_get_mpz(b, a);

            switch (n_randint(state, 5))
            {
                case 0: rnd = ARF_RND_DOWN; break;
                case 1: rnd = ARF_RND_UP; break;
                case 2: rnd = ARF_RND_FLOOR; break;
                case 3: rnd = ARF_RND_CEIL; break;
                default: rnd = ARF_RND_NEAR; break;
            }

            ret1 = arf_set_round_mpz(x, b, bits2, rnd);
            ret2 = mpfr_set_z(g, b, arf_rnd_to_mpfr(rnd));
            arf_set_mpfr(y, g);
            arf_equal(x, y);

            if (!arf_equal(x, y) || ((ret1 == ARF_RESULT_EXACT) != (ret2 == 0)))
            {
                flint_printf("FAIL\n\n");
                flint_printf("bits1: %wd\n", bits1);
                flint_printf("bits2: %wd\n", bits2);
                flint_printf("a = "); fmpz_print(a); flint_printf("\n\n");
                flint_printf("x = "); arf_print(x); flint_printf("\n\n");
                flint_printf("y = "); arf_print(y); flint_printf("\n\n");
                flint_printf("ret1 = %d, ret2 = %d\n\n", ret1, ret2);
                flint_abort();
            }

            fmpz_clear(a);
            mpz_clear(b);
            mpfr_clear(g);
        }

        arf_clear(x);
        arf_clear(y);
    }

    TEST_FUNCTION_END(state);
}
