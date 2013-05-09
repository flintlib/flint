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
#include <limits.h>
#include <gmp.h>
#include "flint.h"
#include "arith.h"
#include "fmpz_vec.h"
#include "ulong_extras.h"
#include "profiler.h"
#include "nmod_vec.h"

int main(void)
{
    flint_rand_t state;
    fmpz * p;
    mp_ptr pmod;
    len_t k, n;

    const len_t maxn = 1000;

    printf("number_of_partitions_vec....");
    fflush(stdout);

    flint_randinit(state);
    p = _fmpz_vec_init(maxn);
    pmod = _nmod_vec_init(maxn);

    for (n = 0; n < maxn; n += (n < 50) ? + 1 : n/4)
    {
        fmpz_t s, t;
        nmod_t mod;
        nmod_init(&mod, n_randtest_prime(state, 0));

        arith_number_of_partitions_vec(p, n);
        arith_number_of_partitions_nmod_vec(pmod, n, mod);

        for (k = 0; k < n; k++)
        {
            if (fmpz_fdiv_ui(p + k, mod.n) != pmod[k])
            {
                printf("FAIL:\n");
                printf("n = %ld, k = %ld\n", n, k);
                abort();
            }
        }

        if (n > 1)
        {
            fmpz_init(s);
            fmpz_init(t);

            for (k = 1; k < n; k++)
            {
                len_t j;

                j = n - 1 - k*(3*k - 1)/2;
                if (j >= 0)
                    fmpz_set(t, p + j);
                else
                    fmpz_zero(t);

                j = n - 1 - k*(3*k + 1)/2;
                if (j >= 0)
                    fmpz_add(t, t, p + j);

                if (k % 2)
                    fmpz_add(s, s, t);
                else
                    fmpz_sub(s, s, t);
            }

            if (!fmpz_equal(s, p + n - 1))
            {
                printf("FAIL:\n");
                printf("n = %ld\n", n);
                fmpz_print(s);
                printf("\n");
                fmpz_print(p + n - 1);
                printf("\n");
                abort();
            }

            fmpz_clear(s);
            fmpz_clear(t);
        }
    }

    _fmpz_vec_clear(p, maxn);
    _nmod_vec_clear(pmod);

    flint_randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return 0;
}
