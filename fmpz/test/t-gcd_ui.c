/*
    Copyright (C) 2021 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <gmp.h>
#include "flint.h"
#include "fmpz.h"

int
main(void)
{
    int i, result;
    FLINT_TEST_INIT(state);

    flint_printf("gcd_ui....");
    fflush(stdout);
    

    for (i = 0; i < 10000 * flint_test_multiplier(); i++)
    {
        fmpz_t r1, r2, a, fb;
        ulong b;

        fmpz_init(r1);
        fmpz_init(r2);
        fmpz_init(a);
        fmpz_init(fb);

        fmpz_randtest(a, state, 200);
        b = n_randtest(state);
        fmpz_set_ui(fb, b);

        fmpz_gcd(r1, a, fb);
        fmpz_gcd_ui(r2, a, b);

        result = (fmpz_cmp(r1, r2) == 0);

        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("a = ");
            fmpz_print(a);
            flint_printf("\nb = %wu", b);
            flint_printf("\nr1 = ");
            fmpz_print(r1);
            flint_printf("\nr2 = ");
            fmpz_print(r2);
            flint_printf("\n");
            fflush(stdout);
            abort();
        }

        fmpz_clear(r1);
        fmpz_clear(r2);
        fmpz_clear(a);
        fmpz_clear(fb);
    }

    /* Check aliasing of res and a */
    for (i = 0; i < 10000 * flint_test_multiplier(); i++)
    {
        fmpz_t r, a, fb;
        ulong b;

        fmpz_init(r);
        fmpz_init(a);
        fmpz_init(fb);

        fmpz_randtest(a, state, 200);
        b = n_randtest(state);
        fmpz_set_ui(fb, b);

        fmpz_gcd(r, a, fb);
        fmpz_gcd_ui(a, a, b);

        result = (fmpz_cmp(r, a) == 0);
        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("Aliasing with res and a.\n");
            flint_printf("a = ");
            fmpz_print(a);
            flint_printf("\nb = %wu", b);
            flint_printf("\nr = ");
            fmpz_print(r);
            flint_printf("\n");
            fflush(stdout);
            abort();
        }

        fmpz_clear(r);
        fmpz_clear(a);
        fmpz_clear(fb);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}
