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

    Copyright (C) 2011 Fredrik Johansson

******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include "flint.h"
#include "ulong_extras.h"


static mp_limb_t
n_factorial_mod2_foolproof(ulong n, mp_limb_t p, mp_limb_t pinv)
{
    mp_limb_t prod = 1UL % p;

    while (n)
    {
        prod = n_mulmod2_preinv(prod, n, p, pinv);
        n--;
    }

    return prod;
}

int main(void)
{
    flint_rand_t state;
    mp_limb_t n;
    int j;
    flint_randinit(state);

    printf("factorial_mod2_preinv....");
    fflush(stdout);

    for (n = 0; n < 100 * flint_test_multiplier(); n++)
    {
        mp_limb_t p, pinv, x, y;

        for (j = 0; j < 10; j++)
        {
            p = n_randtest_not_zero(state);
            pinv = n_preinvert_limb(p);
            x = n_factorial_mod2_preinv(n, p, pinv);
            y = n_factorial_mod2_foolproof(n, p, pinv);

            if (x != y)
            {
                printf("FAIL:\n");
                printf("n = %lu\np = %lu\nx = %lu\ny = %lu\n", n, p, x, y);
                abort();
            }
        }
    }

    flint_randclear(state);
    printf("PASS\n");
    return 0;
}
