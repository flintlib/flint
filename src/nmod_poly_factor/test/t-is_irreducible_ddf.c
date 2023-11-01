/*
    Copyright (C) 2011 Fredrik Johansson
    Copyright (C) 2013 Martin Lee

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "nmod_poly.h"
#include "nmod_poly_factor.h"

TEST_FUNCTION_START(nmod_poly_factor_is_irreducible_ddf, state)
{
    int iter;

    for (iter = 0; iter < 200 * flint_test_multiplier(); iter++)
    {
        nmod_poly_t poly, poly2, poly3;
        mp_limb_t modulus;
        slong length, length2;
        int result = 1;

        modulus = n_randtest_prime(state, 0);

        nmod_poly_init(poly, modulus);
        nmod_poly_init(poly2, modulus);
        nmod_poly_init(poly3, modulus);

        length = n_randint(state, 10) + 2;

        do
        {
            nmod_poly_randtest(poly, state, length);
            if(!nmod_poly_is_zero(poly))
                nmod_poly_make_monic(poly, poly);
        }
        while ((poly->length < 2));

        result &= (nmod_poly_is_irreducible_rabin (poly) == nmod_poly_is_irreducible_ddf (poly));
        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("result of is_irreducible and is_irreducible_ddf does not coincide\n");
            flint_printf("poly = "), nmod_poly_print(poly), flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }

        length2 = n_randint(state, 10) + 2;

        do
        {
            nmod_poly_randtest(poly2, state, length2);
            if(!nmod_poly_is_zero(poly2))
                nmod_poly_make_monic(poly2, poly2);
        }
        while ((!nmod_poly_is_irreducible_rabin(poly2)) || (poly2->length < 2));

        nmod_poly_mul(poly3, poly, poly2);

        result &= !nmod_poly_is_irreducible_ddf(poly3);
        if (!result)
        {
            flint_printf("Error: reducible polynomial declared irreducible!\n");
            nmod_poly_print(poly3); flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }

        nmod_poly_clear(poly);
        nmod_poly_clear(poly2);
        nmod_poly_clear(poly3);
    }

    TEST_FUNCTION_END(state);
}
