/*
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "nmod_poly.h"
#include "nmod_poly_factor.h"

TEST_FUNCTION_START(nmod_poly_factor_is_squarefree, state)
{
    int iter;

    for (iter = 0; iter < 200 * flint_test_multiplier(); iter++)
    {
        nmod_poly_t poly, Q, R, t;
        mp_limb_t modulus;
        slong i, num_factors, exp, max_exp;
        int v, result;

        modulus = n_randtest_prime(state, 0);

        nmod_poly_init(poly, modulus);
        nmod_poly_init(t, modulus);
        nmod_poly_init(Q, modulus);
        nmod_poly_init(R, modulus);

        nmod_poly_set_coeff_ui(poly, 0, n_randint(state, modulus));
        num_factors = n_randint(state, 5);

        max_exp = 0;
        for (i = 0; i < num_factors; i++)
        {
            do {
                nmod_poly_randtest(t, state, n_randint(state, 10));
            } while (!nmod_poly_is_irreducible(t) ||
                    (nmod_poly_length(t) < 2));

            exp = n_randint(state, 4) + 1;
            if (n_randint(state, 2) == 0)
                exp = 1;

            nmod_poly_divrem(Q, R, poly, t);
            if (!nmod_poly_is_zero(R))
            {
                nmod_poly_pow(t, t, exp);
                nmod_poly_mul(poly, poly, t);
                max_exp = FLINT_MAX(exp, max_exp);
            }
        }

        v = nmod_poly_is_squarefree(poly);

        if (v == 1)
            result = (max_exp <= 1 && !nmod_poly_is_zero(poly));
        else
            result = (max_exp > 1 || nmod_poly_is_zero(poly));

        if (!result)
        {
            flint_printf("FAIL: %wu, %wd, %d\n", modulus, max_exp, v);
            nmod_poly_print(poly); flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }

        nmod_poly_clear(poly);
        nmod_poly_clear(t);
        nmod_poly_clear(Q);
        nmod_poly_clear(R);
    }

    TEST_FUNCTION_END(state);
}
