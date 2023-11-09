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

TEST_FUNCTION_START(nmod_poly_factor_berlekamp, state)
{
    int iter;

    for (iter = 0; iter < 20 * flint_test_multiplier(); iter++)
    {
        int result = 1;
        nmod_poly_t pol1, poly, quot, rem;
        nmod_poly_factor_t res;
        mp_limb_t modulus;
        slong i, length, num;

        modulus = n_randtest_prime(state, 0);

        nmod_poly_init(pol1, modulus);
        nmod_poly_init(poly, modulus);
        nmod_poly_init(quot, modulus);
        nmod_poly_init(rem, modulus);

        length = n_randint(state, 10) + 2;
        do
        {
            nmod_poly_randtest(pol1, state, length);
            if (pol1->length)
                nmod_poly_make_monic(pol1, pol1);
        }
        while ((!nmod_poly_is_irreducible(pol1)) || (pol1->length < 2));

        num = n_randint(state, 5) + 1;
        for (i = 1; i < num; i++)
        {
            do
            {
                length = n_randint(state, 10) + 2;
                nmod_poly_randtest(poly, state, length);
                if (poly->length)
                {
                    nmod_poly_make_monic(poly, poly);
                    nmod_poly_divrem(quot, rem, pol1, poly);
                }
            }
            while ((!nmod_poly_is_irreducible(poly)) || (poly->length < 2)
                || (rem->length == 0));
            nmod_poly_mul(pol1, pol1, poly);
        }

        nmod_poly_factor_init(res);
        nmod_poly_factor_berlekamp(res, pol1);

        result = (res->num == num);
        if (!result)
        {
            flint_printf("FAIL: %wu, %wd, %wd\n", modulus, num, res->num);
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
