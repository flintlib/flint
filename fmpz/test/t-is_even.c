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
#include "flint.h"
#include "ulong_extras.h"
#include "fmpz.h"

int
main(void)
{
    int i, result;
    flint_rand_t state;

    printf("is_even/odd....");
    fflush(stdout);

    flint_randinit(state);

    for (i = 0; i < 10000 * flint_test_multiplier(); i++)
    {
        fmpz_t f;
        mpz_t g;

        fmpz_init(f);
        mpz_init(g);

        fmpz_randtest(f, state, 100);
        fmpz_get_mpz(g, f);

        result = (fmpz_is_even(f) == mpz_even_p(g));
        if (!result)
        {
            printf("FAIL:\n");
            printf("f = "), fmpz_print(f), printf("\n");
            gmp_printf("g = %Zd\n", g);
            abort();
        }

        fmpz_clear(f);
        mpz_clear(g);
    }

    for (i = 0; i < 10000 * flint_test_multiplier(); i++)
    {
        fmpz_t f;
        mpz_t g;

        fmpz_init(f);
        mpz_init(g);

        fmpz_randtest(f, state, 100);
        fmpz_get_mpz(g, f);

        result = (fmpz_is_odd(f) == mpz_odd_p(g));
        if (!result)
        {
            printf("FAIL:\n");
            printf("f = "), fmpz_print(f), printf("\n");
            gmp_printf("g = %Zd\n", g);
            abort();
        }

        fmpz_clear(f);
        mpz_clear(g);
    }

    flint_randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return 0;
}
