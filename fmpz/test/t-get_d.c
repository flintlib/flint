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

******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gmp.h>
#include "flint.h"
#include "ulong_extras.h"
#include "fmpz.h"

int
main(void)
{
    int i, result;
    flint_rand_t state;

    printf("get_d....");
    fflush(stdout);

    flint_randinit(state);

    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpz_t x, y;
        mpz_t z;
        double a, b;

        fmpz_init(x);
        fmpz_init(y);
        mpz_init(z);

        fmpz_randtest(x, state, 200);
        fmpz_get_mpz(z, x);

        a = fmpz_get_d(x);
        b = mpz_get_d(z);

        result = (a == b);
        if (!result)
        {
            printf("FAIL:\n");
            printf("x = "), fmpz_print(x), printf("\n");
            printf("a = %f\n", a);
            printf("b = %f\n", b);
            abort();
        }

        a = a * (n_randtest(state) / (double) n_randtest_not_zero(state));

        fmpz_set_d(x, a);
        mpz_set_d(z, a);

        fmpz_set_mpz(y, z);
        result = fmpz_equal(x, y);

        if (!result)
        {
            printf("FAIL:\n");
            printf("x = "), fmpz_print(x), printf("\n");
            printf("y = "), fmpz_print(y), printf("\n");
            printf("a = %f\n", a);
            abort();
        }

        fmpz_clear(x);
        fmpz_clear(y);
        mpz_clear(z);
    }

    flint_randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return 0;
}
