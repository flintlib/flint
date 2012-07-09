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

    Copyright (C) 2011 Fredrik Johansson
    Copyright (C) 2012 Lina Kulakova

******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <mpir.h>
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

    printf("is_irreducible....");
    fflush(stdout);

    for (iter = 0; iter < 100; iter++)
    {
        fmpz_mod_poly_t poly1, poly2;
        fmpz_t modulus;
        long length, test_modulus;
        int result = 1;
        int i, num;

        do
        {
            test_modulus = n_randtest_prime(state, 0);
        } while (test_modulus <= 2);
        fmpz_init_set_ui(modulus, test_modulus);

        fmpz_mod_poly_init(poly1, modulus);
        fmpz_mod_poly_init(poly2, modulus);

        length = n_randint(state, 10) + 2;
        do
        {
            fmpz_mod_poly_randtest(poly1, state, length);
            fmpz_mod_poly_make_monic(poly1, poly1);
        }
        while ((!fmpz_mod_poly_is_irreducible(poly1)) || (poly1->length < 2));

        num = n_randint(state, 5) + 1;

        for (i = 0; i < num; i++)
        {
            do
            {
                fmpz_mod_poly_randtest(poly2, state, length);
                fmpz_mod_poly_make_monic(poly2, poly2);
            }
            while ((!fmpz_mod_poly_is_irreducible(poly2))
                   || (poly2->length < 2));

            fmpz_mod_poly_mul(poly1, poly1, poly2);
        }

        result &= !fmpz_mod_poly_is_irreducible(poly1);
        if (!result)
        {
            printf("Error: reducible polynomial declared irreducible!\n");
            printf("poly:\n");
            fmpz_mod_poly_print(poly1);
            printf("\n");
            abort();
        }

        fmpz_clear(modulus);
        fmpz_mod_poly_clear(poly1);
        fmpz_mod_poly_clear(poly2);
    }

    flint_randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return 0;
}
