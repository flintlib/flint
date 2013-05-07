/*============================================================================

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

===============================================================================*/
/****************************************************************************

   Copyright (C) 2011 Fredrik Johansson

*****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "ulong_extras.h"


int main()
{
    long i, j;
    int sign;

    fmpz_t input;
    fmpz_t result;
    fmpz_t r1;
    fmpz_t m1;
    fmpz_t mprod;
    ulong r2, m2;

    flint_rand_t state;

    printf("CRT_ui....");
    fflush(stdout);

    fmpz_init(input);
    fmpz_init(result);
    fmpz_init(r1);
    fmpz_init(m1);
    fmpz_init(mprod);
    flint_randinit(state);

    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        long nprimes;

        m2 = n_randtest_prime(state, 0);
        nprimes = 1 + n_randint(state, 4);

        fmpz_set_ui(m1, 1UL);
        for (j = 0; j < nprimes; )
        {
            ulong t = n_randtest_prime(state, 0);
            if (t != m2)
            {
                fmpz_mul_ui(m1, m1, t);
                j++;
            }
        }

        fmpz_mul_ui(mprod, m1, m2);

        sign = n_randint(state, 2);

        if (sign)
            fmpz_randtest_mod_signed(input, state, mprod);
        else
            fmpz_randtest_mod(input, state, mprod);

        fmpz_mod(r1, input, m1);
        r2 = fmpz_fdiv_ui(input, m2);

        fmpz_CRT_ui(result, r1, m1, r2, m2, sign);

        if (!fmpz_equal(result, input))
        {
            printf("FAIL:\n");
            printf("m1: "); fmpz_print(m1); printf("\n");
            printf("m2: %lu\n", m2);
            printf("m1*m2: "); fmpz_print(mprod); printf("\n");
            printf("input: "); fmpz_print(input); printf("\n");
            printf("r1: "); fmpz_print(r1); printf("\n");
            printf("r2: %lu\n", r2);
            printf("result: "); fmpz_print(result); printf("\n");
            printf("%ld Equalness: %d\n", i, fmpz_equal(result, input));
            printf("\n");
            abort();
        }
    }

    fmpz_clear(input);
    fmpz_clear(result);
    fmpz_clear(r1);
    fmpz_clear(m1);
    fmpz_clear(mprod);
    flint_randclear(state);

    _fmpz_cleanup();

    printf("PASS\n");
    return 0;
}
