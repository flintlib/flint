/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "mpn_extras.h"
#include "fmpz.h"

TEST_FUNCTION_START(flint_mpn_fmms1, state)
{
    slong i, j;
    fmpz_t X1, X2, Y, T;
    mp_limb_t * y, a1, * x1, a2, * x2;
    mp_size_t yn, n;

    fmpz_init(X1);
    fmpz_init(X2);
    fmpz_init(Y);
    fmpz_init(T);

    for (i = 0; i < 10000*flint_test_multiplier(); i++)
    {
        n = n_randint(state, 100) + 1;

        x1 = (mp_limb_t *) flint_malloc(n*sizeof(mp_limb_t));
        x2 = (mp_limb_t *) flint_malloc(n*sizeof(mp_limb_t));
        y = (mp_limb_t *) flint_malloc(n*sizeof(mp_limb_t));

        a1 = n_randtest(state);
        a2 = n_randtest(state);
        for (j = 0; j < n; j++)
        {
            x1[j] = n_randtest(state);
            x2[j] = n_randtest(state);
        }

        fmpz_set_ui_array(X1, x1, n);
        fmpz_set_ui_array(X2, x2, n);
        fmpz_mul_ui(Y, X1, a1);
        fmpz_submul_ui(Y, X2, a2);

        yn = flint_mpn_fmms1(y, a1, x1, a2, x2, n);

        if (yn > 0)
        {
            fmpz_set_ui_array(T, y, yn);
            if (y[yn - 1] == 0 || !fmpz_equal(T, Y))
            {
                flint_printf("FAIL\n");
                flint_printf("check positive answer, i = %wd\n", i);
                fflush(stdout);
                flint_abort();
            }
        }
        else if (yn == 0)
        {
            if (!fmpz_is_zero(Y))
            {
                flint_printf("FAIL\n");
                flint_printf("check zero answer, i = %wd\n", i);
                fflush(stdout);
                flint_abort();
            }
        }
        else
        {
            if (fmpz_sgn(Y) >= 0 && fmpz_size(Y) <= n)
            {
                flint_printf("FAIL\n");
                flint_printf("check failed answer, i = %wd\n", i);
                fflush(stdout);
                flint_abort();
            }
        }

        flint_free(x1);
        flint_free(x2);
        flint_free(y);
    }

    fmpz_clear(X1);
    fmpz_clear(X2);
    fmpz_clear(Y);
    fmpz_clear(T);

    TEST_FUNCTION_END(state);
}
