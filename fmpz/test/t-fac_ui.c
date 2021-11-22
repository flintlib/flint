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
    fmpz_t x;
    fmpz_t y;

    FLINT_TEST_INIT(state);
    
    flint_printf("fac_ui....");
    fflush(stdout);

    fmpz_init(x);
    fmpz_init(y);

    /* Twice to check demotion */
    for (n = 0; n < 2; n++)
    {
        fmpz_set_ui(y, UWORD(1));

        for (i = 0; i < 100; i++)
        {
            fmpz_fac_ui(x, i);
            fmpz_mul_ui(y, y, FLINT_MAX(1, i));
            if (!fmpz_equal(x, y))
            {
                flint_printf("FAIL: %wd\n", i);
                fmpz_print(x);
                flint_printf("\n");
                fmpz_print(y);
                flint_printf("\n");
                fflush(stdout);
                flint_abort();
            }
        }
    }

    fmpz_clear(x);
    fmpz_clear(y);

    FLINT_TEST_CLEANUP(state);
    flint_printf("PASS\n");
    return 0;
}
