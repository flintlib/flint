/*
    Copyright (C) 2007 David Howden
    Copyright (C) 2007, 2008, 2009, 2010 William Hart
    Copyright (C) 2008 Richard Howell-Peak
    Copyright (C) 2011 Fredrik Johansson
    Copyright (C) 2012 Lina Kulakova

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdlib.h>
#include "fmpz_mod_poly.h"
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
        fmpz_mod_poly_t poly1, poly, q, r, product;
        fmpz_mod_poly_factor_t res;
        fmpz_t modulus;
        slong i, j, length, num;
        slong exp[5];

        fmpz_init_set_ui(modulus, n_randtest_prime(state, 0));

        fmpz_mod_poly_init(poly1, modulus);
        fmpz_mod_poly_init(poly, modulus);
        fmpz_mod_poly_init(q, modulus);
        fmpz_mod_poly_init(r, modulus);

        fmpz_mod_poly_zero(poly1);
        fmpz_mod_poly_set_coeff_ui(poly1, 0, 1);

        length = n_randint(state, 7) + 2;
        do
        {
            fmpz_mod_poly_randtest(poly, state, length);
            if (poly->length)
                fmpz_mod_poly_make_monic(poly, poly);
        }
        while ((poly->length < 2) || (!fmpz_mod_poly_is_irreducible(poly)));

        exp[0] = n_randint(state, 30) + 1;
        for (i = 0; i < exp[0]; i++)
            fmpz_mod_poly_mul(poly1, poly1, poly);

        num = n_randint(state, 5) + 1;

        for (i = 1; i < num; i++)
        {
            do
            {
                length = n_randint(state, 7) + 2;
                fmpz_mod_poly_randtest(poly, state, length);
                if (poly->length)
                {
                    fmpz_mod_poly_make_monic(poly, poly);
                    fmpz_mod_poly_divrem(q, r, poly1, poly);
                }
            }
            while ((poly->length < 2) || (!fmpz_mod_poly_is_irreducible(poly))
                   || (r->length == 0));

            exp[i] = n_randint(state, 30) + 1;
            for (j = 0; j < exp[i]; j++)
                fmpz_mod_poly_mul(poly1, poly1, poly);
        }

        fmpz_mod_poly_factor_init(res);
        fmpz_mod_poly_factor_cantor_zassenhaus(res, poly1);

        if (res->num != num)
        {
            flint_printf("Error: number of factors incorrect: %wd != %wd\n",
                   res->num, num);
            abort();
        }

        fmpz_mod_poly_init(product, &poly1->p);
        fmpz_mod_poly_set_coeff_ui(product, 0, 1);
        for (i = 0; i < res->num; i++)
            for (j = 0; j < res->exp[i]; j++)
                fmpz_mod_poly_mul(product, product, res->poly + i);

        fmpz_mod_poly_scalar_mul_fmpz(product, product,
                                      &(poly1->coeffs[poly1->length - 1]));

        if (!fmpz_mod_poly_equal(poly1, product))
        {
            flint_printf
                ("Error: product of factors does not equal to the original polynomial\n");
            flint_printf("poly:\n");
            fmpz_mod_poly_print(poly1);
            flint_printf("\n");
            flint_printf("product:\n");
            fmpz_mod_poly_print(product);
            flint_printf("\n");
            abort();
        }

        fmpz_clear(modulus);
        fmpz_mod_poly_clear(product);
        fmpz_mod_poly_clear(q);
        fmpz_mod_poly_clear(r);
        fmpz_mod_poly_clear(poly1);
        fmpz_mod_poly_clear(poly);
        fmpz_mod_poly_factor_clear(res);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}
