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

    Copyright (C) 2011 Sebastian Pancratz

******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "ulong_extras.h"

int main(void)
{
    int i, result;
    flint_rand_t state;

    printf("sqrtmodn....");
    fflush(stdout);

    flint_randinit(state);
    
    for (i = 0; i < 1000 * flint_test_multiplier(); i++) /* Test random squares mod n */
    {
        mp_limb_t a, b, n, ninv;
        len_t num, i;
        mp_bitcnt_t bits;
        mp_limb_t * sqrt;
        int btest;
        n_factor_t fac;

        bits = n_randint(state, 18) + 2;
        n = n_randtest_bits(state, bits);
        if (n == 0) n = 1;
        b = n_randtest(state) % n;
        
        n_factor_init(&fac);
        n_factor(&fac, n, 0);

        ninv = n_preinvert_limb(n);
        a = n_mulmod2_preinv(b, b, n, ninv);

        num = n_sqrtmodn(&sqrt, a, &fac);
        
        btest = 0;
        for (i = 0; i < num; i++)
        {
            if (a != n_mulmod2_preinv(sqrt[i], sqrt[i], n, ninv))
                break;
            if (sqrt[i] == b)
                btest = 1;
        }

        result = btest & (i == num);
        if (!result)
        {
            printf("FAIL:\n");
            printf("n = %lu\n", n);
            printf("a = %lu\n", a);
            printf("b = %lu\n", b);
            printf("num = %ld\n", num);
            
            if (!btest)
                printf("Square root not found.\n");
            if (i != num)
                printf("%lu not a square root of %lu mod %lu\n", sqrt[i], a, n);

            abort();
        }

        flint_free(sqrt);
    }
    
    for (i = 0; i < 500 * flint_test_multiplier(); i++) /* test random nonsquares */
    {
        mp_limb_t a, b, n, ninv;
        mp_bitcnt_t bits;
        mp_limb_t * sqrt;
        n_factor_t fac;

        bits = n_randint(state, 18) + 2;
        n = n_randtest_bits(state, bits);
        if (n == 2) n++;
        n_factor_init(&fac);
        n_factor(&fac, n, 0);

        ninv = n_preinvert_limb(n);
        
        a = n_randtest(state) % n;
        while (n_sqrtmodn(&sqrt, a, &fac))
        {
            if (n_mulmod2_preinv(sqrt[0], sqrt[0], n, ninv) != a)
            {
                printf("FAIL:\n");
                printf("%lu^2 is not %lu mod %lu\n", sqrt[0], a, n);
                abort();
            }
            
            flint_free(sqrt);
            a = n_randtest(state) % n;
        }
        
        for (b = 0; b < n; b++)
        {
            if (n_mulmod2_preinv(b, b, n, ninv) == a)
                break;
        }

        result = (b == n);
        if (!result)
        {
            printf("FAIL:\n");
            printf("n = %lu\n", n);
            printf("a = %lu\n", a);
            printf("b = %lu\n", b);

            abort();
        }

        flint_free(sqrt);
    }

    flint_randclear(state);

    printf("PASS\n");
    return 0;
}
