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

#include <stdlib.h>
#include "nmod_poly.h"
#include "ulong_extras.h"

int
main(void)
{
    int iter;
    FLINT_TEST_INIT(state);
    

    flint_printf("factor_cantor_zassenhaus....");
    fflush(stdout);

    for (iter = 0; iter < 20 * flint_test_multiplier(); iter++)
    {
        int result = 1;
        nmod_poly_t pol1, poly, quot, rem;
        nmod_poly_t product;
        nmod_poly_factor_t res;
        mp_limb_t modulus, lead;
        slong i, j, length, num;
        slong prod1, exp[5];

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
        while ((poly->length < 2) || (!nmod_poly_is_irreducible(poly)));

        exp[0] = n_randint(state, 30) + 1;
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
                if(!nmod_poly_is_zero(poly))
                {
                    nmod_poly_make_monic(poly, poly);
                    nmod_poly_divrem(quot, rem, pol1, poly);
                }
            }
            while ((!nmod_poly_is_irreducible(poly)) || (poly->length < 2) ||
                (rem->length == 0));

            exp[i] = n_randint(state, 30) + 1;
            prod1 *= exp[i];
            for (j = 0; j < exp[i]; j++)
                nmod_poly_mul(pol1, pol1, poly);
        }

        nmod_poly_factor_init(res);
        nmod_poly_factor_cantor_zassenhaus(res, pol1);
        result &= (res->num == num);

        if (!result)
        {
            flint_printf("Error: number of factors incorrect, %wd, %wd\n",
                res->num, num);
        }

        nmod_poly_init(product, pol1->mod.n);
        nmod_poly_set_coeff_ui(product, 0, 1);
        for (i = 0; i < res->num; i++)
            for (j = 0; j < res->exp[i]; j++)
                nmod_poly_mul(product, product, res->p + i);

        lead = pol1->coeffs[pol1->length - 1];
        nmod_poly_scalar_mul_nmod(product, product, lead);
        result &= nmod_poly_equal(pol1, product);

        if (!result)
        {
            flint_printf("Error: product of factors does not equal original polynomial\n");
            nmod_poly_print(pol1); flint_printf("\n");
            nmod_poly_print(product); flint_printf("\n");
        }

        if (!result)
            abort();


        nmod_poly_clear(product);
        nmod_poly_clear(quot);
        nmod_poly_clear(rem);
        nmod_poly_clear(pol1);
        nmod_poly_clear(poly);
        nmod_poly_factor_clear(res);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}
