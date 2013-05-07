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

    printf("divisible....");
    fflush(stdout);

    flint_randinit(state);

    /* Compare with MPIR:  random */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpz_t a, b;
        mpz_t c, d;
        int e, f;

        fmpz_init(a);
        fmpz_init(b);
        mpz_init(c);
        mpz_init(d);

        fmpz_randtest_unsigned(a, state, 200);
        fmpz_add_ui(a, a, 1);
        fmpz_randtest(b, state, 200);

        fmpz_get_mpz(c, a);
        fmpz_get_mpz(d, b);

        e = fmpz_divisible(b, a);
        f = mpz_divisible_p(d, c);

        result = (e == f);
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

    /* Compare with MPIR:  b a multiple of a */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpz_t a, b;
        mpz_t c, d;
        int e, f;

        fmpz_init(a);
        fmpz_init(b);
        mpz_init(c);
        mpz_init(d);

        fmpz_randtest_unsigned(a, state, 200);
        fmpz_add_ui(a, a, 1);
        fmpz_randtest(b, state, 200);
        fmpz_mul(b, a, b);

        fmpz_get_mpz(c, a);
        fmpz_get_mpz(d, b);

        e = fmpz_divisible(b, a);
        f = mpz_divisible_p(d, c);

        result = (e == f && e == 1);
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
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpz_t a;
        int b;

        fmpz_init(a);

        fmpz_randtest_unsigned(a, state, 200);
        fmpz_add_ui(a, a, 1);

        b = fmpz_divisible(a, a);

        result = (b == 1);
        if (!result)
        {
            printf("FAIL:\n");
            printf("a = "), fmpz_print(a);
            abort();
        }

        fmpz_clear(a);
    }

    flint_randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}

