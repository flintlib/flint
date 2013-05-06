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
    fmpz * b1;
    fmpz * b2;
    long n, k;

    const long maxn = 400;

    printf("bell_number....");
    fflush(stdout);

    flint_randinit(state);

    b1 = _fmpz_vec_init(maxn);

    /* Consistency test */
    for (n = 0; n < maxn; n++)
        arith_bell_number(b1 + n, n);

    for (n = 0; n < maxn; n++)
    {
        b2 = _fmpz_vec_init(n);
        arith_bell_number_vec(b2, n);

        if (!_fmpz_vec_equal(b1, b2, n))
        {
            printf("FAIL:\n");
            printf("n = %ld\n", n);
            abort();
        }

        _fmpz_vec_clear(b2, n);
    }

    /* Compare with B_n = sum of Stirling numbers of 2nd kind */
    for (n = 0; n < 1000; n += (n < 50) ? + 1 : n/4)
    {
        b2 = _fmpz_vec_init(n+1);

        arith_stirling_number_2_vec(b2, n, n+1);

        for (k = 1; k <= n; k++)
            fmpz_add(b2, b2, b2 + k);

        arith_bell_number(b1, n);

        if (!fmpz_equal(b1, b2))
        {
            printf("FAIL:\n");
            printf("n = %ld\n", n);
            fmpz_print(b1);
            printf("\n");
            fmpz_print(b2);
            printf("\n");
            abort();
        }

        /* Also check nmod value */
        {
            nmod_t mod;
            mp_limb_t bb;

            nmod_init(&mod, n_randtest_prime(state, 0));
            bb = arith_bell_number_nmod(n, mod);

            if (fmpz_fdiv_ui(b1, mod.n) != bb)
            {
                printf("FAIL:\n");
                printf("n = %ld\n", n);
                fmpz_print(b1);
                printf("\n");
                printf("should be %lu mod %lu\n", bb, mod.n);
                abort();
            }
        }

        _fmpz_vec_clear(b2, n+1);
    }

    _fmpz_vec_clear(b1, maxn);

    flint_randclear(state);

    mpfr_free_cache();
    _fmpz_cleanup();
    printf("PASS\n");
    return 0;
}
