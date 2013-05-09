/*=============================================================================

    This file is part of FLINT.

    FLINT is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    FLINT is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with FLINT; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301 USA

=============================================================================*/
/******************************************************************************

    Copyright (C) 2007 David Howden
    Copyright (C) 2007, 2008, 2009, 2010 William Hart
    Copyright (C) 2008 Richard Howell-Peak
    Copyright (C) 2011 Fredrik Johansson
    Copyright (C) 2012 Lina Kulakova

******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "fmpz_vec.h"
#include "fmpz_mod_poly_factor.h"
#include "ulong_extras.h"

int
main(void)
{
    int iter;
    flint_rand_t state;
    flint_randinit(state);

    printf("factor_squarefree....");
    fflush(stdout);

    for (iter = 0; iter < 300; iter++)
    {
        int result = 1;
        fmpz_mod_poly_t pol1, poly, quot, rem;
        fmpz_mod_poly_factor_t res;
        fmpz_t modulus;
        len_t exp[5], prod1;
        len_t length, i, j, num;

        fmpz_init_set_ui(modulus, n_randtest_prime(state, 0));

        fmpz_mod_poly_init(pol1, modulus);
        fmpz_mod_poly_init(poly, modulus);
        fmpz_mod_poly_init(quot, modulus);
        fmpz_mod_poly_init(rem, modulus);

        fmpz_mod_poly_zero(pol1);
        fmpz_mod_poly_set_coeff_ui(pol1, 0, 1);

        length = n_randint(state, 7) + 2;

        do
        {
            fmpz_mod_poly_randtest(poly, state, length);
            fmpz_mod_poly_make_monic(poly, poly);
        }
        while ((!fmpz_mod_poly_is_irreducible(poly)) || (poly->length < 2));
        exp[0] = n_randprime(state, 5, 0);

        prod1 = exp[0];
        for (i = 0; i < exp[0]; i++)
            fmpz_mod_poly_mul(pol1, pol1, poly);

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
                    fmpz_mod_poly_divrem(quot, rem, pol1, poly);
                }
            }
            while ((!fmpz_mod_poly_is_irreducible(poly)) ||
                   (poly->length < 2) || (rem->length == 0));

            do
                exp[i] = n_randprime(state, 5, 0);
            while (prod1 % exp[i] == 0);

            prod1 *= exp[i];
            for (j = 0; j < exp[i]; j++)
                fmpz_mod_poly_mul(pol1, pol1, poly);
        }

        fmpz_mod_poly_factor_init(res);
        fmpz_mod_poly_factor_squarefree(res, pol1);

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
            printf("Error: exp don't match. Modulus = ");
            fmpz_print(modulus);
            printf("\n");
            for (i = 0; i < res->num; i++)
                printf("%ld ", res->exp[i]);
            printf("\n");
            for (i = 0; i < num; i++)
                printf("%ld ", exp[i]);
            printf("\n");
            abort();
        }

        fmpz_clear(modulus);
        fmpz_mod_poly_clear(quot);
        fmpz_mod_poly_clear(rem);
        fmpz_mod_poly_clear(pol1);
        fmpz_mod_poly_clear(poly);
        fmpz_mod_poly_factor_clear(res);
    }

    flint_randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return 0;
}
