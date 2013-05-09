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

    printf("sqrtmod_primepow....");
    fflush(stdout);

    flint_randinit(state);

    for (i = 0; i < 1000 * flint_test_multiplier(); i++) /* Test random squares mod a power of 2 */
    {
        mp_limb_t a, b, p, pow, pow2, pinv;
        len_t exp, num, i;
        mp_limb_t * sqrt;
        int btest;

        p = 2;
        exp = n_randint(state, FLINT_BITS - 1) + 1;
        pow = n_pow(p, exp);
        b = n_randtest(state) % pow;

        pow2 = p;
        while (FLINT_BIT_COUNT(p*pow2) <= 12)
            pow2 *= p;

        if ((b % (p*pow2)) == 0) 
        {
            b += pow2;
            b %= pow;
        }
        
        pinv = n_preinvert_limb(pow);
        a = n_mulmod2_preinv(b, b, pow, pinv);

        num = n_sqrtmod_primepow(&sqrt, a, p, exp);

        btest = 0;
        for (i = 0; i < num; i++)
        {
            if (a != n_mulmod2_preinv(sqrt[i], sqrt[i], pow, pinv))
                break;
            if (sqrt[i] == b)
                btest = 1;
        }

        result = btest & (i == num);
        if (!result)
        {
            printf("FAIL:\n");
            printf("p = %lu\n", p);
            printf("exp = %ld\n", exp);
            printf("a = %lu\n", a);
            printf("b = %lu\n", b);
            printf("num = %ld\n", num);
            
            if (!btest)
                printf("Square root not found.\n");
            if (i != num)
                printf("%lu not a square root of %lu mod %lu\n", sqrt[i], a, pow);

            abort();
        }

        flint_free(sqrt);
    }
    
    for (i = 0; i < 1000 * flint_test_multiplier(); i++) /* Test random squares mod other prime powers */
    {
        mp_limb_t a, b, p, pow, pow2, pinv;
        len_t exp, maxexp, num, i;
        mp_bitcnt_t bits;
        mp_limb_t * sqrt;
        int btest;

        bits = n_randint(state, 18) + 2;
        p = n_randprime(state, bits, 0);
        maxexp = FLINT_BITS/bits;
        exp = n_randint(state, maxexp) + 1;
        pow = n_pow(p, exp);
        b = n_randtest(state) % pow;
        
        if (bits <= FLINT_BITS/2)
        { 
            pow2 = p;
            while (FLINT_BIT_COUNT(p*pow2) <= 12)
                pow2 *= p;

            if ((b % (p*pow2)) == 0) 
                b += pow2;

            b %= pow;
        }

        pinv = n_preinvert_limb(pow);
        a = n_mulmod2_preinv(b, b, pow, pinv);

        num = n_sqrtmod_primepow(&sqrt, a, p, exp);

        btest = 0;
        for (i = 0; i < num; i++)
        {
            if (a != n_mulmod2_preinv(sqrt[i], sqrt[i], pow, pinv))
                break;
            if (sqrt[i] == b)
                btest = 1;
        }

        result = btest & (i == num);
        if (!result)
        {
            printf("FAIL:\n");
            printf("p = %lu\n", p);
            printf("exp = %ld\n", exp);
            printf("a = %lu\n", a);
            printf("b = %lu\n", b);
            printf("num = %ld\n", num);
            
            if (!btest)
                printf("Square root not found.\n");
            if (i != num)
                printf("%lu not a square root of %lu mod %lu\n", sqrt[i], a, pow);

            abort();
        }

        flint_free(sqrt);
    }

    for (i = 0; i < 500 * flint_test_multiplier(); i++) /* Test random nonsquares */
    {
        mp_limb_t a, b, p, pow, pinv;
        len_t exp, maxexp;
        mp_bitcnt_t bits;
        mp_limb_t * sqrt;
        
        bits = n_randint(state, 18) + 2;
        p = n_randprime(state, bits, 0);
        maxexp = 20/bits;
        exp = n_randint(state, maxexp) + 1 + (p == 2);
        pow = n_pow(p, exp);
        
        pinv = n_preinvert_limb(pow);
        
        a = n_randtest(state) % pow;
        while (n_sqrtmod_primepow(&sqrt, a, p, exp))
        {
            if (n_mulmod2_preinv(sqrt[0], sqrt[0], pow, pinv) != a)
            {
                printf("FAIL:\n");
                printf("%lu^2 is not %lu mod %lu\n", sqrt[0], a, pow);
                abort();
            }
            
            flint_free(sqrt);
            a = n_randtest(state) % pow;
        }
        
        for (b = 0; b < pow; b++)
        {
            if (n_mulmod2_preinv(b, b, pow, pinv) == a)
                break;
        }

        result = (b == pow);
        if (!result)
        {
            printf("FAIL:\n");
            printf("p = %lu\n", p);
            printf("exp = %ld\n", exp);
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
