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
    mp_bitcnt_t r1, r2;
    flint_rand_t state;

    printf("popcnt....");
    fflush(stdout);

    flint_randinit(state);

    for (i = 0; i < 10000 * flint_test_multiplier(); i++)
    {
        fmpz_t a;
        mpz_t b;

        fmpz_init(a);
        mpz_init(b);

        fmpz_randtest(a, state, 2 * FLINT_BITS);
        fmpz_get_mpz(b, a);

        r1 = fmpz_popcnt(a);
        r2 = mpz_popcount(b);

        result = r1 == r2;

        if (!result && fmpz_cmp_si(a,0) >= 0)
        {
            printf("FAIL:\n");
            printf("a = "), fmpz_print(a), printf("\n");
            gmp_printf("b = %Zd\n", b);
            printf("r1 = %lu  r2 = %lu\n", r1, r2);
            abort();
        }

        fmpz_clear(a);
        mpz_clear(b);
    }

    flint_randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return 0;
}
