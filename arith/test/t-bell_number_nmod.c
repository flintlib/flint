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
#include <gmp.h>
#include <mpfr.h>
#include "flint.h"
#include "arith.h"
#include "fmpz_vec.h"
#include "ulong_extras.h"

int main(void)
{
    flint_rand_t state;
    long i, j;

    printf("bell_number_nmod....");
    fflush(stdout);

    flint_randinit(state);

    for (i = 0; i < 10; i++)
    {
        mp_ptr b;
        long n;
        nmod_t mod;
        mp_limb_t p;

        n = n_randint(state, 1000);
        p = n_randtest_prime(state, 0);

        nmod_init(&mod, p);

        b = _nmod_vec_init(n + 1);
        arith_bell_number_nmod_vec(b, n + 1, mod);

        for (j = 0; j <= n; j++)
        {
            mp_limb_t u = arith_bell_number_nmod(j, mod);

            if (u != b[j])
            {
                printf("FAIL: p = %lu, i = %ld\n", p, j);
                abort();
            }
        }

        _nmod_vec_clear(b);
    }

    flint_randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return 0;
}
