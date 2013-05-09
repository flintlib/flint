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
#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "ulong_extras.h"

int
main(void)
{
    len_t i;
    ulong n, k;
    fmpz_t x, y;
    mpz_t z;
    flint_rand_t state;

    printf("bin_uiui....");
    fflush(stdout);

    flint_randinit(state);

    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpz_init(x);
        fmpz_init(y);
        mpz_init(z);

        n = n_randint(state, 1000);
        k = n_randint(state, 1000);

        fmpz_bin_uiui(x, n, k);
        mpz_bin_uiui(z, n, k);
        fmpz_set_mpz(y, z);

        if (!fmpz_equal(x, y))
        {
            printf("FAIL: n,k = %lu,%lu\n", n, k);
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
