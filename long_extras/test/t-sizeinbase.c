/*
    Copyright (C) 2009 William Hart
    Copyright (C) 2011 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "long_extras.h"
#include "ulong_extras.h"

int
main(void)
{
    int i, result;
    FLINT_TEST_INIT(state);

    flint_printf("sizeinbase....");
    fflush(stdout);

    

    for (i = 0; i < 100000; i++)
    {
        slong a;
        mpz_t b;
        int base;
        size_t r1, r2;

        a = z_randtest(state);
        flint_mpz_init_set_si(b, a);
        base = (int) (n_randint(state, 61) + 2);

        r1 = z_sizeinbase(a, base);
        r2 = mpz_sizeinbase(b, base);
        result = (r1 == r2 || r1 + 1 == r2);

        if (!result)
        {
            flint_printf("FAIL:\n");
            gmp_printf("b = %Zd\n", b);
            flint_printf("base = %d\n", base);
            flint_printf("r1 = %wu\n, r2 = %wu\n", (ulong) r1, (ulong) r2);
            abort();
        }

        mpz_clear(b);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}
