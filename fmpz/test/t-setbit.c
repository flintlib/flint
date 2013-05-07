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
    Copyright (C) 2012 Sebastian Pancratz

******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "ulong_extras.h"
#include "fmpz.h"

int
main(void)
{
    int i, result;
    flint_rand_t state;

    printf("setbit....");
    fflush(stdout);

    flint_randinit(state);

    for (i = 0; i < 100000 * flint_test_multiplier(); i++)
    {
        ulong j;
        fmpz_t a, b, c;
        mpz_t z;

        fmpz_init(a);
        fmpz_init(b);
        fmpz_init(c);
        mpz_init(z);

        fmpz_randtest(a, state, 2 * FLINT_BITS);
        fmpz_set(b, a);
        fmpz_get_mpz(z, b);
        j = n_randint(state, 3 * FLINT_BITS);

        fmpz_setbit(b, j);
        mpz_setbit(z, j);
        fmpz_set_mpz(c, z);

        result = (fmpz_equal(b, c));

        if (!result)
        {
            printf("FAIL:\n");
            printf("a = "), fmpz_print(a), printf("\n");
            printf("b = "), fmpz_print(b), printf("\n");
            printf("c = "), fmpz_print(c), printf("\n");
            gmp_printf("z = %Zd\n", z);
            printf("j = %ld\n", j);
            abort();
        }

        fmpz_clear(a);
        fmpz_clear(b);
        fmpz_clear(c);
        mpz_clear(z);
    }

    flint_randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return 0;
}
