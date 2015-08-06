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

    Copyright (C) 2015 Vladimir Glazachev

******************************************************************************/

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

    for (i = 0; i < 10000; i++)
    {
        int pbprime, cycloprime;
        fmpz_t n;
        fmpz_init(n);

        fmpz_randtest_unsigned(n, state, 1000);
        while (fmpz_cmp_ui(n, 100) <= 0)
            fmpz_randtest_unsigned(n, state, 1000);

        pbprime = fmpz_is_probabprime(n);
        cycloprime = is_prime_jacobi(n);
        
        if (pbprime != cycloprime)
        {
            flint_printf("FAIL\n");
            flint_printf("Testing number = ");
            fmpz_print(n);
            flint_printf("\nis_probabprime = %i, is_prime_gauss = %i\n", pbprime, cycloprime);
            abort();
        }

        fmpz_clear(n);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}

