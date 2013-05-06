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

    Copyright (C) 2010 Sebastian Pancratz
    Copyright (C) 2011 Fredrik Johansson

******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "long_extras.h"
#include "fmpz.h"

int
main(void)
{
    int i, result;
    flint_rand_t state;

    printf("equal_si....");
    fflush(stdout);

    flint_randinit(state);

    /* Compare with fmpz_equal, random values */
    for (i = 0; i < 10000 * flint_test_multiplier(); i++)
    {
        fmpz_t a, b;
        long n;
        int lhs, rhs;

        fmpz_init(a);
        fmpz_init(b);

        fmpz_randtest(a, state, 200);

        n = z_randtest(state);
        fmpz_set_si(b, n);

        lhs = fmpz_equal(a, b);
        rhs = fmpz_equal_si(a, n);

        result = (lhs == rhs);
        if (result == 0)
        {
            printf("FAIL:\n");
            printf("a = "), fmpz_print(a), printf("\n");
            printf("b = "), fmpz_print(b), printf("\n");
            printf("n = %ld\n", n);
            printf("equal(a, b) = %d\n", fmpz_equal(a, b));
            printf("equal_si(a, n) = %d\n", fmpz_equal_si(a, n));
            abort();
        }

        fmpz_clear(a);
        fmpz_clear(b);
    }

    /* Compare with fmpz_equal, equal values */
    for (i = 0; i < 10000 * flint_test_multiplier(); i++)
    {
        fmpz_t a, b;
        long n;
        int lhs, rhs;

        fmpz_init(a);
        fmpz_init(b);

        n = z_randtest(state);
        fmpz_set_si(a, n);
        fmpz_set_si(b, n);

        lhs = fmpz_equal(a, b);
        rhs = fmpz_equal_si(a, n);

        result = (lhs == rhs) && (lhs == 1);
        if (result == 0)
        {
            printf("FAIL:\n");
            printf("a = "), fmpz_print(a), printf("\n");
            printf("b = "), fmpz_print(b), printf("\n");
            printf("n = %ld\n", n);
            printf("equal(a, b) = %d\n", fmpz_equal(a, b));
            printf("equal_si(a, n) = %d\n", fmpz_equal_si(a, n));
            abort();
        }

        fmpz_clear(a);
        fmpz_clear(b);
    }

    flint_randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}
