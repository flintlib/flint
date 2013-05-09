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

    Copyright (C) 2012 Fredrik Johansson

******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "ulong_extras.h"

int main(void)
{
    flint_rand_t state;
    len_t j, k, l;

    printf("primes_jump_after....");
    fflush(stdout);

    flint_randinit(state);

    for (j = 0; j < 10; j++)
    {
        n_primes_t iter;

        n_primes_init(iter);

        for (k = 0; k < 100; k++)
        {
            mp_limb_t p, q;

            q = n_randtest(state) % 1000000000UL;

            n_primes_jump_after(iter, q);

            for (l = 0; l < 100; l++)
            {
                p = n_primes_next(iter);
                q = n_nextprime(q, 0);

                if (p != q)
                {
                    printf("FAIL\n");
                    printf("p = %lu, q = %lu\n", p, q);
                    abort();
                }
            }
        }

        n_primes_clear(iter);
    }

    flint_randclear(state);
    printf("PASS\n");
    return 0;
}
