/*
    Copyright (C) 2024 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "ulong_extras.h"

TEST_FUNCTION_START(n_factor_evaluate, state)
{
    int result;
    slong iter;

    for (iter = 0; iter < 10000 * flint_test_multiplier(); iter++)
    {
        ulong eval, fac_eval;
        n_factor_t factors;
        int type;

        n_factor_init(&factors);

        type = n_randint(state, 2);

        if (type == 0)
        {
            /* Test random limb */
            eval = n_randtest_not_zero(state);
            n_factor(&factors, eval, 0);
        }
        else
        {
            /* Test factorization whose evaluation cannot fit in a limb */
            ulong top;

            eval = 1;

            do
            {
                ulong prime = n_randtest_prime(state, 0);
                umul_ppmm(top, eval, eval, prime);
                n_factor_insert(&factors, prime, 1);
            } while (top == UWORD(0));

            /* And maybe add extra exponents */
            if (n_randint(state, 2))
                factors.exp[n_randint(state, factors.num)] += 1 + n_randint(state, 10);

            eval = 0; /* n_factor_evaluate should return 0 */
        }

        fac_eval = n_factor_evaluate(&factors);

        result = (eval == fac_eval);
        if (!result)
            TEST_FUNCTION_FAIL(
                    "eval     = %wu\n"
                    "fac_eval = %wu\n",
                    eval, fac_eval);
    }

    TEST_FUNCTION_END(state);
}
