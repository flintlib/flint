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

TEST_FUNCTION_START(arf_neg_round, state)
{
    slong iter;

    {
        arf_t x, y, z;

        arf_init(x);
        arf_init(y);
        arf_init(z);

        for (iter = 0; iter < 100000 * 0.1 * flint_test_multiplier(); iter++)
        {
            slong bits1, bits2;
            int ret1, ret2;
            mpfr_t g1, g2;
            fmpz_t e;
            arf_rnd_t rnd;

            bits1 = 1 + n_randint(state, 1000);
            bits2 = 2 + n_randint(state, 1000);

            if (n_randint(state, 100) == 0)
                bits2 = ARF_PREC_EXACT;

            switch (n_randint(state, 5))
            {
                case 0: rnd = ARF_RND_DOWN; break;
                case 1: rnd = ARF_RND_UP; break;
                case 2: rnd = ARF_RND_FLOOR; break;
                case 3: rnd = ARF_RND_CEIL; break;
                default: rnd = ARF_RND_NEAR; break;
            }

            fmpz_init(e);
            mpfr_init2(g1, FLINT_MAX(2, bits1));
            mpfr_init2(g2, FLINT_MIN(bits2, 10000));

            if (n_randint(state, 100) == 0)
            {
                arf_clear(x); arf_clear(y); arf_clear(z);
                arf_init(x); arf_init(y); arf_init(z);
            }

            /* dirty output variables */
            if (n_randint(state, 2))
            {
                arf_randtest_special(y, state, 1 + n_randint(state, 1000),
                    1 + n_randint(state, 100));
                arf_randtest_special(z, state, 1 + n_randint(state, 1000),
                    1 + n_randint(state, 100));
            }

            arf_randtest_special(x, state, bits1, 1 + n_randint(state, 10));
            arf_get_mpfr(g1, x, MPFR_RNDD); /* exact */

            /* test large exponents */
            if (n_randint(state, 4) == 0)
                fmpz_randtest(e, state, 1 + n_randint(state, 100));

            if (!arf_is_special(x))
                fmpz_add(ARF_EXPREF(x), ARF_EXPREF(x), e);

            ret1 = arf_neg_round(y, x, bits2, rnd);
            ret2 = mpfr_neg(g2, g1, arf_rnd_to_mpfr(rnd));
            arf_set_mpfr(z, g2);

            if (!arf_is_special(y))
                fmpz_sub(ARF_EXPREF(y), ARF_EXPREF(y), e);

            if (!arf_equal(y, z) || ((ret1 == ARF_RESULT_EXACT) != (ret2 == 0)))
            {
                flint_printf("FAIL\n\n");
                flint_printf("bits1: %wd\n", bits1);
                flint_printf("bits2: %wd\n", bits2);
                flint_printf("x = "); arf_print(x); flint_printf("\n\n");
                flint_printf("y = "); arf_print(y); flint_printf("\n\n");
                flint_printf("z = "); arf_print(z); flint_printf("\n\n");
                flint_printf("ret1 = %d, ret2 = %d\n\n", ret1, ret2);
                flint_abort();
            }

            if (!arf_is_special(x))
                fmpz_add(ARF_EXPREF(x), ARF_EXPREF(x), e);

            ret1 = arf_neg_round(y, x, bits2, rnd);
            arf_set(z, x);
            ret2 = arf_neg_round(z, z, bits2, rnd);

            if (!arf_equal(y, z) || ret1 != ret2)
            {
                flint_printf("FAIL (aliasing)\n\n");
                flint_printf("bits1: %wd\n", bits1);
                flint_printf("bits2: %wd\n", bits2);
                flint_printf("x = "); arf_print(x); flint_printf("\n\n");
                flint_printf("y = "); arf_print(y); flint_printf("\n\n");
                flint_printf("z = "); arf_print(z); flint_printf("\n\n");
                flint_printf("ret1 = %d, ret2 = %d\n\n", ret1, ret2);
                flint_abort();
            }

            mpfr_clear(g1);
            mpfr_clear(g2);
            fmpz_clear(e);
        }

        arf_clear(x);
        arf_clear(y);
        arf_clear(z);
    }

    TEST_FUNCTION_END(state);
}
