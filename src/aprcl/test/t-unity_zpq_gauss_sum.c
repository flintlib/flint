/*
    Copyright (C) 2015 Vladimir Glazachev

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz.h"
#include "aprcl.h"

TEST_FUNCTION_START(aprcl_unity_zpq_gauss_sum, state)
{
    int i, j;

    for (i = 0; i < 10 * flint_test_multiplier(); i++)
    {
        int result;
        ulong p, q, pnum, ppow;
        fmpz_t n;
        unity_zpq gausssigma, gauss, gausspower;
        n_factor_t factors;

        n_factor_init(&factors);

        q = n_randprime(state, 6, 0);
        if (q == 2)
            q = 7;

        n_factor(&factors, q - 1, 0);

        pnum = n_randint(state, factors.num);
        p = factors.p[pnum];
        ppow = n_randint(state, factors.exp[pnum]);
        if (ppow == 0)
            ppow = 1;

        p = n_pow(p, ppow);

        fmpz_init_set_ui(n, n_randprime(state, 16, 0));

        unity_zpq_init(gausssigma, q, p, n);
        unity_zpq_init(gauss, q, p, n);
        unity_zpq_init(gausspower, q, p, n);

        unity_zpq_gauss_sum(gauss, q, p);
        unity_zpq_gauss_sum_sigma_pow(gausssigma, q, p);

        unity_zpq_pow(gausspower, gauss, n);

        result = 0;
        for (j = 0; j < p; j++)
        {
            unity_zpq_mul_unity_p_pow(gauss, gausspower, j);
            if (unity_zpq_equal(gauss, gausssigma))
            {
                result = 1;
                break;
            }
        }

        if (result == 0)
        {
            flint_printf("FAIL\n");
            fflush(stdout);
            flint_abort();
        }

        unity_zpq_clear(gausssigma);
        unity_zpq_clear(gauss);
        unity_zpq_clear(gausspower);
        fmpz_clear(n);

    }

    TEST_FUNCTION_END(state);
}
