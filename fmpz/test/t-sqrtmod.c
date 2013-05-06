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
#include "fmpz.h"
#include "ulong_extras.h"

int main(void)
{
    int i, result;
    flint_rand_t state;

    printf("sqrtmod....");
    fflush(stdout);

    flint_randinit(state);

    for (i = 0; i < 100 * flint_test_multiplier(); i++) /* Test random integers */
    {
        int ans;
        fmpz_t a, b, c, p;
        mp_limb_t prime;

        prime = n_randint(state, 1UL << (FLINT_BITS - 1));
        prime = n_nextprime(prime, 1);

        fmpz_init(a);
        fmpz_init(b);
        fmpz_init(c);
        fmpz_init(p);

        fmpz_set_ui(p, prime);
        fmpz_randm(a, state, p);

        ans = fmpz_sqrtmod(b, a, p);

        fmpz_mul(c, b, b);
        fmpz_mod(c, c, p);

        result = (ans == 0 || fmpz_equal(a, c));
        if (!result)
        {
            printf("FAIL (random):\n");
            printf("p = "), fmpz_print(p), printf("\n");
            printf("a = "), fmpz_print(a), printf("\n");
            printf("b = "), fmpz_print(b), printf("\n");
            printf("c = "), fmpz_print(c), printf("\n");
            abort();
        }

        fmpz_clear(a);
        fmpz_clear(b);
        fmpz_clear(c);
        fmpz_clear(p);
    }

    for (i = 0; i < 100 * flint_test_multiplier(); i++) /* Test random squares */
    {
        int ans;
        fmpz_t a, b, c, d, p;
        mp_limb_t prime;

        prime = n_randint(state, 1UL << (FLINT_BITS - 1));
        prime = n_nextprime(prime, 1);

        fmpz_init(a);
        fmpz_init(b);
        fmpz_init(c);
        fmpz_init(d);
        fmpz_init(p);

        fmpz_set_ui(p, prime);
        do 
            fmpz_randm(b, state, p);
        while (fmpz_is_zero(b));

        fmpz_mul(a, b, b);
        fmpz_mod(a, a, p);

        /* check a special case */
        if (i == 0)
        {
            fmpz_set_str(p, "15951355998396157", 10);
            fmpz_set_str(a, "7009303413761286", 10);
        }

        ans = fmpz_sqrtmod(c, a, p);

        fmpz_mul(d, c, c);
        fmpz_mod(d, d, p);

        result = (ans && fmpz_equal(a, d));
        if (!result)
        {
            printf("FAIL (squares):\n");
            printf("p            = "), fmpz_print(p), printf("\n");
            printf("a (= b^2)    = "), fmpz_print(a), printf("\n");
            printf("b            = "), fmpz_print(b), printf("\n");
            printf("c (= sqrt(a) = "), fmpz_print(c), printf("\n");
            printf("d (= c^2)    = "), fmpz_print(d), printf("\n");
            printf("ans          = %d\n", ans);
            abort();
        }

        fmpz_clear(a);
        fmpz_clear(b);
        fmpz_clear(c);
        fmpz_clear(d);
        fmpz_clear(p);
    }

    flint_randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return 0;
}
