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

    Copyright (C) 2012 William Hart

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

    printf("is_prime_pseudosquare....");
    fflush(stdout);

    flint_randinit(state);

    for (i = 0; i < 100000 * flint_test_multiplier(); i++)
    {
        fmpz_t p;
        int r1, r2;

        fmpz_init(p);

        fmpz_randtest_unsigned(p, state, n_randint(state, 94) + 1);

        r1 = fmpz_is_probabprime(p);
        r2 = fmpz_is_prime_pseudosquare(p);

        result = (r1 == r2);
        if (!result)
        {
            printf("FAIL:\n");
            fmpz_print(p);
            printf("r1 = %d, r2 = %d\n", r1, r2);
            abort();
        }

        fmpz_clear(p);
    }

    flint_randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return 0;
}
