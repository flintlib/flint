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

    Copyright (C) 2008, 2009, William Hart
    Copyright (C) 2010 Fredrik Johansson

******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "ulong_extras.h"
#include "fmpz.h"


int main()
{
    len_t i, j;
    mp_limb_t n;
    len_t num_primes;
    mp_limb_t * primes;
    mp_limb_t p;
    fmpz_comb_t comb;
    flint_rand_t state;
    flint_randinit(state);

    printf("comb_init/clear....");
    fflush(stdout);

    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        n = n_randint(state, 10);
        num_primes = (1L << n);
        primes = (mp_limb_t *) flint_malloc(num_primes * sizeof(mp_limb_t));
        p = n_nextprime((1UL << (FLINT_BITS-1)) - 10000000L, 0);

        for (j = 0; j < num_primes; j++)
        {
            primes[j] = p;
            p = n_nextprime(p, 0);
        }

        fmpz_comb_init(comb, primes, num_primes);
        fmpz_comb_clear(comb);
        flint_free(primes);
    }

    flint_randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return 0;
}
