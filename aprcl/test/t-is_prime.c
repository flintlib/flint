/*
    Copyright (C) 2015 Vladimir Glazachev

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
#include "aprcl.h"

int main(void)
{
    int i;
    FLINT_TEST_INIT(state);
   
    flint_printf("is_prime....");
    fflush(stdout);

    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        int pbprime, cycloprime;
        fmpz_t n;
        fmpz_init(n);

        fmpz_randtest_unsigned(n, state, 1000);
        while (fmpz_cmp_ui(n, 100) <= 0)
            fmpz_randtest_unsigned(n, state, 1000);

        pbprime = fmpz_is_probabprime(n);
        cycloprime = aprcl_is_prime_jacobi(n);
        
        if (pbprime != cycloprime)
        {
            flint_printf("FAIL\n");
            flint_printf("Testing number = ");
            fmpz_print(n);
            flint_printf("\nis_probabprime = %i, aprcl_is_prime_jacobi = %i\n", pbprime, cycloprime);
            abort();
        }

        fmpz_clear(n);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}

