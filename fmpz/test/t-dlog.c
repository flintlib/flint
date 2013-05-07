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
#include "ulong_extras.h"
#include "fmpz.h"

int
main(void)
{
    int i, result;
    flint_rand_t state;

    printf("dlog....");
    fflush(stdout);

    flint_randinit(state);

    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpz_t x;
        mpz_t z;
        mpfr_t r;
        double y, w;

        fmpz_init(x);
        mpz_init(z);
        mpfr_init2(r, 53);

        fmpz_randtest_not_zero(x, state, 10000);
        fmpz_abs(x, x);

        y = fmpz_dlog(x);

        fmpz_get_mpz(z, x);
        mpfr_set_z(r, z, MPFR_RNDN);

        mpfr_log(r, r, MPFR_RNDN);
        w = mpfr_get_d(r, MPFR_RNDN);

        result = (FLINT_ABS(y - w) <= w * 1e-13);

        if (!result)
        {
            printf("FAIL:\n");
            printf("x = "), fmpz_print(x), printf("\n");
            printf("y = %.20g\n", y);
            printf("w = %.20g\n", w);
            abort();
        }

        fmpz_clear(x);
        mpz_clear(z);
        mpfr_clear(r);
    }

    mpfr_free_cache();
    flint_randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}

