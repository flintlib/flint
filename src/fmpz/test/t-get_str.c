/*
    Copyright (C) 2009 William Hart
    Copyright (C) 2010 Sebastian Pancratz

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
#include "ulong_extras.h"
#include "fmpz.h"

int
main(void)
{
    int i, result;
    FLINT_TEST_INIT(state);

    flint_printf("get_str....");
    fflush(stdout);

    

    for (i = 0; i < 10000 * flint_test_multiplier(); i++)
    {
        fmpz_t a;
        mpz_t b;
        int base, j;
        char *str1, *str2;

        fmpz_init(a);
        mpz_init(b);
        fmpz_randtest(a, state, 200);
        base = (int) (n_randint(state, 61) + 2);

        fmpz_get_mpz(b, a);

        str1 = fmpz_get_str(NULL, base, a);
        str2 = mpz_get_str(NULL, base, b);
        result = strlen(str1) == strlen(str2);
        if (result)
        {
            for (j = 0; result && j < strlen(str1); j++)
                if (str1[j] != str2[j])
                    result = 0;
        }

        if (!result)
        {
            flint_printf("FAIL:\n");
            gmp_printf("b = %Zd\n", b);
            flint_printf("base = %d\n", base);
            flint_printf("str1 = %s\n, str2 = %s\n", str1, str2);
            fflush(stdout);
            flint_abort();
        }

        flint_free(str1);
        flint_free(str2);

        fmpz_clear(a);
        mpz_clear(b);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}
