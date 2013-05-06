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

    printf("sqrtmod....");
    fflush(stdout);

    flint_randinit(state);

    for (i = 0; i < 1000 * flint_test_multiplier(); i++) /* Test random integers */
    {
        mp_limb_t a, b, p, pinv;

        p = n_randtest_prime(state, 0);
        a = n_randtest(state) % p;

        b = n_sqrtmod(a, p);
        pinv = n_preinvert_limb(p);

        result = (b == 0 || n_mulmod2_preinv(b, b, p, pinv) == a);
        if (!result)
        {
            printf("FAIL:\n");
            printf("p = %lu\n", p);
            printf("a = %lu\n", a);
            printf("b = %lu\n", b);
            abort();
        }
    }

    for (i = 0; i < 1000 * flint_test_multiplier(); i++) /* Test random squares */
    {
        mp_limb_t a, b, p, pinv;

        p = n_randtest_prime(state, 0);

        do 
            b = n_randtest(state) % p;
        while (b == 0);

        pinv = n_preinvert_limb(p);
        a = n_mulmod2_preinv(b, b, p, pinv);

        b = n_sqrtmod(a, p);

        result = (n_mulmod2_preinv(b, b, p, pinv) == a);
        if (!result)
        {
            printf("FAIL:\n");
            printf("p = %lu\n", p);
            printf("a = %lu\n", a);
            printf("b = %lu\n", b);
            abort();
        }
    }

    flint_randclear(state);

    printf("PASS\n");
    return 0;
}
