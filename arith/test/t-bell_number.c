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
    fmpz * b1;
    fmpz * b2;
    slong n, k;

    const slong maxn = 400;

    FLINT_TEST_INIT(state);

    flint_printf("bell_number....");
    fflush(stdout);    

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
            flint_printf("FAIL:\n");
            flint_printf("n = %wd\n", n);
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
            flint_printf("FAIL:\n");
            flint_printf("n = %wd\n", n);
            fmpz_print(b1);
            flint_printf("\n");
            fmpz_print(b2);
            flint_printf("\n");
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
                flint_printf("FAIL:\n");
                flint_printf("n = %wd\n", n);
                fmpz_print(b1);
                flint_printf("\n");
                flint_printf("should be %wu mod %wu\n", bb, mod.n);
                abort();
            }
        }

        _fmpz_vec_clear(b2, n+1);
    }

    _fmpz_vec_clear(b1, maxn);

    FLINT_TEST_CLEANUP(state);

    flint_printf("PASS\n");
    return 0;
}
