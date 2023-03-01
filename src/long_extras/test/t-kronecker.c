/*
    Copyright (C) 2021 Daniel Schultz

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
#include "long_extras.h"
#include "ulong_extras.h"
#include "gmpcompat.h"

int main(void)
{
    slong i;
    mpz_t aa, nn;
    FLINT_TEST_INIT(state);

    flint_printf("kronecker....");
    fflush(stdout);

    mpz_init(aa);
    mpz_init(nn);

    for (i = 0; i < 300000 * flint_test_multiplier(); i++)
    {
        slong a = z_randtest(state);
        slong n = z_randtest(state);

        if (n_randlimb(state) & 1)
            a = -a;

        if (n_randlimb(state) & 1)
            n = -n;

        flint_mpz_set_si(aa, a);
        flint_mpz_set_si(nn, n);

        if (mpz_kronecker(aa, nn) != z_kronecker(a, n))
        {
            flint_printf("FAIL:\n");
            flint_printf("a = %wd, n = %wd\n", a, n);
            fflush(stdout);
            flint_abort();
        }
    }

    mpz_clear(aa);
    mpz_clear(nn);

    FLINT_TEST_CLEANUP(state);

    flint_printf("PASS\n");
    return 0;
}
