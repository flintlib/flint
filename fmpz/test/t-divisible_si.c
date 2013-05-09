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
#include <limits.h>
#include "flint.h"
#include "ulong_extras.h"
#include "long_extras.h"
#include "fmpz.h"

int
main(void)
{
    int i, result;
    flint_rand_t state;

    printf("divisible_si....");
    fflush(stdout);

    flint_randinit(state);

    /* Compare with MPIR:  random */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        len_t a;
        fmpz_t b;
        mpz_t d;
        int e, f;

        fmpz_init(b);
        mpz_init(d);

        a = z_randtest(state);
        if (a == LONG_MIN)
            a = 1;
        a = FLINT_ABS(a) + 1;
        fmpz_randtest(b, state, 200);

        fmpz_get_mpz(d, b);

        e = fmpz_divisible_si(b, a);
        f = mpz_divisible_ui_p(d, a);

        result = (e == f);
        if (!result)
        {
            printf("FAIL:\n");
            printf("a = %ld, b = ", a), fmpz_print(b), printf("\n");
            abort();
        }

        fmpz_clear(b);
        mpz_clear(d);
    }

    /* Compare with MPIR:  b a multiple of a */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        len_t a;
        fmpz_t b;
        mpz_t d;
        int e, f;

        fmpz_init(b);
        mpz_init(d);

        a = z_randtest(state);
        if (a == LONG_MIN)
            a = 1;
        a = FLINT_ABS(a) + 1;
        fmpz_randtest(b, state, 200);
        fmpz_mul_ui(b, b, a);

        fmpz_get_mpz(d, b);

        e = fmpz_divisible_si(b, a);
        f = mpz_divisible_ui_p(d, a);

        result = (e == f && e == 1);
        if (!result)
        {
            printf("FAIL:\n");
            printf("a = %ld, b = ", a), fmpz_print(b), printf("\n");
            abort();
        }

        fmpz_clear(b);
        mpz_clear(d);
    }

    flint_randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}

