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

    Copyright (C) 2009 William Hart

******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "ulong_extras.h"
#include "fmpz.h"

int
main(void)
{
    int i, result;
    flint_rand_t state;

    printf("abs....");
    fflush(stdout);

    flint_randinit(state);

    for (i = 0; i < 10000 * flint_test_multiplier(); i++)
    {
        fmpz_t a, b;
        mpz_t c, d;

        fmpz_init(a);
        fmpz_init(b);
        mpz_init(c);
        mpz_init(d);

        fmpz_randtest(a, state, 200);

        fmpz_get_mpz(c, a);

        fmpz_abs(b, a);
        mpz_abs(c, c);

        fmpz_get_mpz(d, b);

        result = (mpz_cmp(c, d) == 0);
        if (!result)
        {
            printf("FAIL:\n");
            gmp_printf("c = %Zd, d = %Zd\n", c, d);
            abort();
        }

        fmpz_clear(a);
        fmpz_clear(b);
        mpz_clear(c);
        mpz_clear(d);
    }

    /* Check aliasing */
    for (i = 0; i < 10000 * flint_test_multiplier(); i++)
    {
        fmpz_t a;
        mpz_t c, d;

        fmpz_init(a);
        mpz_init(c);
        mpz_init(d);

        fmpz_randtest(a, state, 200);

        fmpz_get_mpz(c, a);

        fmpz_abs(a, a);
        mpz_abs(c, c);

        fmpz_get_mpz(d, a);

        result = (mpz_cmp(c, d) == 0);
        if (!result)
        {
            printf("FAIL:\n");
            gmp_printf("c = %Zd, d = %Zd\n", c, d);
            abort();
        }

        fmpz_clear(a);
        mpz_clear(c);
        mpz_clear(d);
    }

    flint_randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return 0;
}
