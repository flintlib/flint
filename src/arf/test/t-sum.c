/*
    Copyright (C) 2013 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "arf.h"

TEST_FUNCTION_START(arf_sum, state)
{
    slong iter;

    for (iter = 0; iter < 1000000 * 0.1 * flint_test_multiplier(); iter++)
    {
        slong i, len, prec, bits, expbits;
        int res1, res2;
        arf_t s1, s2, s3, err;
        mag_t err_bound;
        arf_struct terms[20];
        arf_rnd_t rnd;

        len = n_randint(state, 20);
        bits = 2 + n_randint(state, 1000);
        prec = 2 + n_randint(state, 1000);
        expbits = n_randint(state, 14);

        arf_init(s1);
        arf_init(s2);
        arf_init(s3);
        arf_init(err);
        mag_init(err_bound);

        for (i = 0; i < len; i++)
        {
            arf_init(terms + i);
            arf_randtest_special(terms + i, state, bits, expbits);
        }

        switch (n_randint(state, 5))
        {
            case 0:  rnd = ARF_RND_DOWN; break;
            case 1:  rnd = ARF_RND_UP; break;
            case 2:  rnd = ARF_RND_FLOOR; break;
            case 3:  rnd = ARF_RND_CEIL; break;
            default: rnd = ARF_RND_NEAR; break;
        }

        res1 = arf_sum(s1, terms, len, prec, rnd);

        arf_zero(s2);
        for (i = 0; i < len; i++)
            arf_add(s2, s2, terms + i, ARF_PREC_EXACT, ARF_RND_DOWN);
        res2 = arf_set_round(s3, s2, prec, rnd);

        if (!arf_equal(s1, s3) || res1 != res2)
        {
            flint_printf("FAIL (%wd)\n\n", iter);
            flint_printf("prec = %wd\n\n", prec);
            for (i = 0; i < len; i++)
            {
                flint_printf("terms[%wd] = ", i); arf_print(terms + i); flint_printf("\n\n");
            }
            flint_printf("s1 = "); arf_print(s1); flint_printf("\n\n");
            flint_printf("s2 = "); arf_print(s2); flint_printf("\n\n");
            flint_printf("s3 = "); arf_print(s3); flint_printf("\n\n");
            flint_printf("res1 = %d, res2 = %d\n\n", res1, res2);
            flint_abort();
        }

        arf_sub(err, s1, s2, ARF_PREC_EXACT, ARF_RND_DOWN);
        arf_abs(err, err);

        if (res1)
            arf_mag_set_ulp(err_bound, s1, prec);
        else
            mag_zero(err_bound);

        if (arf_cmpabs_mag(err, err_bound) > 0)
        {
            flint_printf("FAIL (error bound)!\n");
            flint_printf("prec = %wd\n\n", prec);
            for (i = 0; i < len; i++)
            {
                flint_printf("terms[%wd] = ", i); arf_print(terms + i); flint_printf("\n\n");
            }
            flint_printf("s1 = "); arf_print(s1); flint_printf("\n\n");
            flint_printf("s2 = "); arf_print(s2); flint_printf("\n\n");
            flint_printf("s3 = "); arf_print(s3); flint_printf("\n\n");
            flint_printf("error: "); arf_print(err); flint_printf("\n\n");
            flint_printf("error bound: "); mag_print(err_bound); flint_printf("\n\n");
            flint_printf("res1 = %d, res2 = %d\n\n", res1, res2);
            flint_abort();
        }

        arf_clear(s1);
        arf_clear(s2);
        arf_clear(s3);
        arf_clear(err);
        mag_clear(err_bound);

        for (i = 0; i < len; i++)
            arf_clear(terms + i);
    }

    TEST_FUNCTION_END(state);
}
