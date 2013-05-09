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

    Copyright (C) 2009 William Hart

******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_poly.h"
#include "ulong_extras.h"

int
main(void)
{
    int i, j, result;
    flint_rand_t state;

    printf("get/set_coeff_mpz....");
    fflush(stdout);

    flint_randinit(state);

    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fmpz_poly_t a;
        fmpz_t x1, x2;
        mpz_t y1, y2;
        len_t coeff, len;

        fmpz_poly_init(a);
        fmpz_init(x1);
        fmpz_init(x2);
        mpz_init(y1);
        mpz_init(y2);
        len = n_randint(state, 100) + 1;

        for (j = 0; j < 1000; j++)
        {
            fmpz_randtest(x1, state, 200);
            fmpz_get_mpz(y1, x1);
            coeff = n_randint(state, len);
            fmpz_poly_set_coeff_mpz(a, coeff, y1);
            fmpz_poly_get_coeff_mpz(y2, a, coeff);
            fmpz_set_mpz(x2, y2);

            result = (fmpz_equal(x1, x2));
            if (!result)
            {
                printf("FAIL:\n");
                printf("x1 = "), fmpz_print(x1), printf("\n");
                printf("x2 = "), fmpz_print(x2), printf("\n");
                printf("coeff = %ld, length = %ld\n", coeff, len);
                abort();
            }
        }

        fmpz_clear(x1);
        fmpz_clear(x2);
        mpz_clear(y1);
        mpz_clear(y2);
        fmpz_poly_clear(a);
    }

    flint_randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return 0;
}
