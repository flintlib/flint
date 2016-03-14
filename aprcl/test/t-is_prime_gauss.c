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
   
    flint_printf("is_prime_gauss....");
    fflush(stdout);

    for (i = 0; i < 100; i++)
    {
        int pbprime, cycloprime;
        fmpz_t n;
        fmpz_init(n);

        fmpz_randtest_unsigned(n, state, 50);
        while (fmpz_cmp_ui(n, 100) <= 0)
            fmpz_randtest_unsigned(n, state, 50);

        pbprime = fmpz_is_probabprime(n);
        cycloprime = is_prime_gauss(n);
        
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

    {
        int result;
        fmpz_t n;
        fmpz_init(n);
        result = 1;

        /* Test big primes. */
        fmpz_set_str(n, "40206835204840513073", 10);
        if (is_prime_gauss(n) == 0)
            result = 0;

        /* 521419622856657689423872613771 % 4 == 3 */
        fmpz_set_str(n, "521419622856657689423872613771", 10);
        if (is_prime_gauss(n) == 0)
            result = 0;

        /* 5991810554633396517767024967580894321153 % 4 == 1 */
        fmpz_set_str(n, "5991810554633396517767024967580894321153", 10);
        if (is_prime_gauss(n) == 0)
            result = 0;

        /* Test big composite. */
        /* 1500450271 * 5915587277 */
        fmpz_set_str(n, "8876044532898802067", 10);
        if (is_prime_gauss(n) == 1)
            result = 0;

        /* 5915587277 * 54673257461630679457 */
        fmpz_set_str(n, "323424426232167763068694468589", 10);
        if (is_prime_gauss(n) == 1)
            result = 0;

        /* 48112959837082048697 * 66405897020462343733 */
        fmpz_set_str(n, "3194984256290911228520362769161858765901", 10);
        if (is_prime_gauss(n) == 1)
            result = 0;

        if (result == 0)
        {
            flint_printf("FAIL\n");
            abort();
        }

        fmpz_clear(n);
    }

    FLINT_TEST_CLEANUP(state);

    flint_printf("PASS\n");
    return 0;
}

