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

TEST_FUNCTION_START(nmod_poly_factor, state)
{
    int iter;

    /* Default algorithm */
    for (iter = 0; iter < 10 * flint_test_multiplier(); iter++)
    {
        int result = 1;
        nmod_poly_t pol1, poly, quot, rem, product;
        nmod_poly_factor_t res;
        mp_limb_t modulus, lead = 1;
        slong length, num, i, j;
        ulong exp[5];

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
            if (poly->length)
                nmod_poly_make_monic(poly, poly);
        }
        while ((!nmod_poly_is_irreducible(poly)) || (poly->length < 2));

        exp[0] = n_randint(state, 30) + 1;
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
            exp[i] = n_randint(state, 30) + 1;

            for (j = 0; j < exp[i]; j++)
                nmod_poly_mul(pol1, pol1, poly);
        }

        nmod_poly_factor_init(res);

        switch (n_randint(state, 3))
        {
            case 0:
                lead = nmod_poly_factor(res, pol1);
                break;
            case 1:
                lead = nmod_poly_factor_with_berlekamp(res, pol1);
                break;
            case 2:
                if (modulus == 2)
                    lead = nmod_poly_factor(res, pol1);
                else
                    lead = nmod_poly_factor_with_cantor_zassenhaus(res, pol1);
                break;
        }

        result &= (res->num == num);
        if (!result)
        {
            flint_printf("Error: number of factors incorrect, %wd, %wd\n",
                res->num, num);
            fflush(stdout);
            flint_abort();
        }

        nmod_poly_init(product, pol1->mod.n);
        nmod_poly_set_coeff_ui(product, 0, 1);
        for (i = 0; i < res->num; i++)
            for (j = 0; j < res->exp[i]; j++)
                nmod_poly_mul(product, product, res->p + i);
        nmod_poly_scalar_mul_nmod(product, product, lead);
        result &= nmod_poly_equal(pol1, product);
        if (!result)
        {
            flint_printf("Error: product of factors does not equal original polynomial\n");
            nmod_poly_print(pol1); flint_printf("\n");
            nmod_poly_print(product); flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }
        nmod_poly_clear(product);

        nmod_poly_clear(quot);
        nmod_poly_clear(rem);
        nmod_poly_clear(pol1);
        nmod_poly_clear(poly);
        nmod_poly_factor_clear(res);
    }

    /* Test deflation trick */
    for (iter = 0; iter < 10 * flint_test_multiplier(); iter++)
    {
        nmod_poly_t pol1, poly, quot, rem;
        nmod_poly_factor_t res, res2;
        mp_limb_t modulus;
        slong length, num, i, j;
        slong exp[5];
        ulong inflation;
        int found;

        do {
            modulus = n_randtest_prime(state, 0);
        } while (modulus == 2); /* To compare with CZ */

        nmod_poly_init(pol1, modulus);
        nmod_poly_init(poly, modulus);
        nmod_poly_init(quot, modulus);
        nmod_poly_init(rem, modulus);

        nmod_poly_zero(pol1);
        nmod_poly_set_coeff_ui(pol1, 0, 1);

        inflation = n_randint(state, 7) + 1;

        length = n_randint(state, 7) + 2;
        do
        {
            nmod_poly_randtest(poly, state, length);
            if (poly->length)
                nmod_poly_make_monic(poly, poly);
        }
        while ((!nmod_poly_is_irreducible(poly)) || (poly->length < 2));
        nmod_poly_inflate(poly, poly, inflation);

        exp[0] = n_randint(state, 6) + 1;
        for (i = 0; i < exp[0]; i++)
            nmod_poly_mul(pol1, pol1, poly);

        num = n_randint(state, 5) + 1;
        for (i = 1; i < num; i++)
        {
            do
            {
                length = n_randint(state, 6) + 2;
                nmod_poly_randtest(poly, state, length);
                if (poly->length)
                {
                    nmod_poly_make_monic(poly, poly);
                    nmod_poly_divrem(quot, rem, pol1, poly);
                }
            }
            while ((!nmod_poly_is_irreducible(poly)) ||
                (poly->length < 2) || (rem->length == 0));
            exp[i] = n_randint(state, 6) + 1;
            nmod_poly_inflate(poly, poly, inflation);

            for (j = 0; j < exp[i]; j++)
                nmod_poly_mul(pol1, pol1, poly);
        }

        nmod_poly_factor_init(res);
        nmod_poly_factor_init(res2);

        switch (n_randint(state, 3))
        {
            case 0:
                nmod_poly_factor(res, pol1);
                break;
            case 1:
                nmod_poly_factor_with_berlekamp(res, pol1);
                break;
            case 2:
                nmod_poly_factor_with_cantor_zassenhaus(res, pol1);
                break;
        }

        nmod_poly_factor_cantor_zassenhaus(res2, pol1);

        if (res->num != res2->num)
        {
            flint_printf("FAIL: different number of factors found\n");
            fflush(stdout);
            flint_abort();
        }

        for (i = 0; i < res->num; i++)
        {
            found = 0;
            for (j = 0; j < res2->num; j++)
            {
                if (nmod_poly_equal(res->p + i, res2->p + j) &&
                        res->exp[i] == res2->exp[j])
                {
                    found = 1;
                    break;
                }
            }

            if (!found)
            {
                flint_printf("FAIL: factor not found\n");
                fflush(stdout);
                flint_abort();
            }
        }

        nmod_poly_clear(quot);
        nmod_poly_clear(rem);
        nmod_poly_clear(pol1);
        nmod_poly_clear(poly);
        nmod_poly_factor_clear(res);
        nmod_poly_factor_clear(res2);
    }

    TEST_FUNCTION_END(state);
}
