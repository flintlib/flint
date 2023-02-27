/*
    Copyright (C) 2009 William Hart
    Copyright (C) 2010 Fredrik Johansson

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
#include "ulong_extras.h"
#include "fmpz.h"

int
main(void)
{
    slong i, n;
    fmpz_t x, y, z, w;

    FLINT_TEST_INIT(state);
    
    flint_printf("fib_ui....");
    fflush(stdout);

    fmpz_init(x);
    fmpz_init(y);
    fmpz_init(z);
    fmpz_init(w);

    /* Twice to check demotion */
    for (n = 0; n < 2; n++)
    {
        for (i = 0; i < 200; i++)
        {
            fmpz_fib_ui(x, i);
            fmpz_fib_ui(y, i+1);
            fmpz_fib_ui(z, i+2);
            fmpz_add(w, x, y);

            if (!fmpz_equal(w, z))
            {
                flint_printf("FAIL: %wd\n", i);
                fmpz_print(x);
                flint_printf("\n");
                fmpz_print(y);
                flint_printf("\n");
                fmpz_print(z);
                flint_printf("\n");
                fflush(stdout);
                flint_abort();
            }
        }
    }

    fmpz_clear(x);
    fmpz_clear(y);
    fmpz_clear(z);
    fmpz_clear(w);

    FLINT_TEST_CLEANUP(state);
    flint_printf("PASS\n");
    return 0;
}
