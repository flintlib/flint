/*
    Copyright (C) 2014 Abhinav Baid

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "mpf_vec.h"

TEST_FUNCTION_START(mpf_vec_norm2, state)
{
    int i, result;

    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        mpf *a;
        mpf_t res1, res2, res3;
        slong len = n_randint(state, 100);
        if (!len)
            continue;

        a = _mpf_vec_init(len, 200);
        _mpf_vec_randtest(a, state, len, 200);

        mpf_inits(res1, res2, res3, NULL);
        _mpf_vec_norm2(res1, a, len - 1, 200);
        _mpf_vec_norm2(res2, a + len - 1, 1, 200);
        _mpf_vec_norm2(res3, a, len, 200);

        mpf_add(res1, res1, res2);
        result = mpf_cmp(res1, res3);
        if (result)
        {
            flint_printf("FAIL:\n");
            flint_printf("%d\n", len);
            mpf_out_str(stdout, 10, 0, res1);
            flint_printf("\n");
            mpf_out_str(stdout, 10, 0, res3);
            flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }

        _mpf_vec_clear(a, len);
        mpf_clears(res1, res2, res3, NULL);
    }

    TEST_FUNCTION_END(state);
}
