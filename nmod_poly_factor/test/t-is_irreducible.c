/*
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "nmod_vec.h"
#include "nmod_poly.h"
#include "ulong_extras.h"

int
main(void)
{
    int iter;
    FLINT_TEST_INIT(state);
    

    flint_printf("is_irreducible....");
    fflush(stdout);

    for (iter = 0; iter < 200 * flint_test_multiplier(); iter++)
    {
        nmod_poly_t poly, poly2, poly3;
        nmod_poly_factor_t factors;
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
        while ((!nmod_poly_is_irreducible(poly)) || (poly->length < 2));

        nmod_poly_factor_init(factors);
        nmod_poly_factor_berlekamp(factors, poly);
        result &= (factors->num == 1);
        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("Irreducible polynomial should not have non-trivial factors!\n");
            flint_printf("poly = "), nmod_poly_print(poly), flint_printf("\n");
            abort();
        }
        nmod_poly_factor_clear(factors);

        length2 = n_randint(state, 10) + 2;

        do 
        {
            nmod_poly_randtest(poly2, state, length2); 
            if(!nmod_poly_is_zero(poly2))
                nmod_poly_make_monic(poly2, poly2);
        }
        while ((!nmod_poly_is_irreducible(poly2)) || (poly2->length < 2));

        nmod_poly_mul(poly3, poly, poly2);

        result &= !nmod_poly_is_irreducible(poly3);
        if (!result)
        {
            flint_printf("Error: reducible polynomial declared irreducible!\n");
            nmod_poly_print(poly3); flint_printf("\n");
            abort();
        }

        nmod_poly_clear(poly);
        nmod_poly_clear(poly2);
        nmod_poly_clear(poly3);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}
