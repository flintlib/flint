/*
    Copyright (C) 2007 David Howden
    Copyright (C) 2007, 2008, 2009, 2010 William Hart
    Copyright (C) 2008 Richard Howell-Peak
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

TEST_FUNCTION_START(nmod_poly_factor_squarefree, state)
{
    int iter;

    for (iter = 0; iter < 30 * flint_test_multiplier(); iter++)
    {
        int result = 1;
        nmod_poly_t pol1, poly, quot, rem;
        nmod_poly_factor_t res;
        mp_limb_t modulus;
        slong exp[5], prod1;
        slong length, i, j, num;

        modulus = n_randtest_prime(state, 0);

        nmod_poly_init(pol1, modulus);
        nmod_poly_init(poly, modulus);
        nmod_poly_init(quot, modulus);
        nmod_poly_init(rem, modulus);

        nmod_poly_zero(pol1);
        nmod_poly_set_coeff_ui(pol1, 0, 1);

        length = n_randint(state, 7) + 2;

        do
        {
            nmod_poly_randtest(poly, state, length);
            if(!nmod_poly_is_zero(poly))
                nmod_poly_make_monic(poly, poly);
        }
        while ((!nmod_poly_is_irreducible(poly)) || (poly->length < 2));
        exp[0] = n_randprime(state, 5, 0);

        prod1 = exp[0];
        for (i = 0; i < exp[0]; i++)
            nmod_poly_mul(pol1, pol1, poly);

        num = n_randint(state, 5) + 1;
        for (i = 1; i < num; i++)
        {
            do
            {
                length = n_randint(state, 7) + 2;
                nmod_poly_randtest(poly, state, length);
                if (poly->length)
                {
                    nmod_poly_make_monic(poly, poly);
                    nmod_poly_divrem(quot, rem, pol1, poly);
                }
            }
            while ((!nmod_poly_is_irreducible(poly)) ||
                (poly->length < 2) || (rem->length == 0));

            do exp[i] = n_randprime(state, 5, 0);
            while (prod1 % exp[i] == 0);

            prod1 *= exp[i];
            for (j = 0; j < exp[i]; j++)
                nmod_poly_mul(pol1, pol1, poly);
        }

        nmod_poly_factor_init(res);
        nmod_poly_factor_squarefree(res, pol1);

        result &= (res->num == num);
        if (result)
        {
            ulong prod2 = 1;
            for (i = 0; i < num; i++)
                prod2 *= res->exp[i];
            result &= (prod1 == prod2);
        }

        if (!result)
        {
            flint_printf("Error: exp don't match. Modulus = %wu\n", modulus);
            for (i = 0; i < res->num; i++)
                flint_printf("%wd ", res->exp[i]);
            flint_printf("\n");
            for (i = 0; i < num; i++)
                flint_printf("%wd ", exp[i]);
            flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }

        nmod_poly_clear(quot);
        nmod_poly_clear(rem);
        nmod_poly_clear(pol1);
        nmod_poly_clear(poly);
        nmod_poly_factor_clear(res);
    }

    TEST_FUNCTION_END(state);
}
