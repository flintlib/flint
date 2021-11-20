/*
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "ulong_extras.h"

int
main(void)
{
    slong i;
    ulong n, k;
    fmpz_t x, y;
    mpz_t z;
    FLINT_TEST_INIT(state);

    flint_printf("bin_uiui....");
    fflush(stdout);

    

    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpz_init(x);
        fmpz_init(y);
        mpz_init(z);

        n = n_randint(state, 1000);
        k = n_randint(state, 1000);

        fmpz_bin_uiui(x, n, k);
        flint_mpz_bin_uiui(z, n, k);
        fmpz_set_mpz(y, z);

        if (!fmpz_equal(x, y))
        {
            flint_printf("FAIL: n,k = %wu,%wu\n", n, k);
            fflush(stdout);
            flint_abort();
        }

        fmpz_clear(x);
        fmpz_clear(y);
        mpz_clear(z);
    }

    

    FLINT_TEST_CLEANUP(state);
    flint_printf("PASS\n");
    return 0;
}
