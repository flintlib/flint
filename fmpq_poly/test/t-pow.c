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
    Copyright (C) 2010 Sebastian Pancratz

******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpq_poly.h"
#include "ulong_extras.h"

int
main(void)
{
    int i, result;
    flint_rand_t state;
    ulong cflags = 0UL;

    printf("pow....");
    fflush(stdout);

    flint_randinit(state);

    /* Check aliasing of a and b */
    for (i = 0; i < 200 * flint_test_multiplier(); i++)
    {
        fmpq_poly_t a, b;
        ulong exp;

        fmpq_poly_init(a);
        fmpq_poly_init(b);
        fmpq_poly_randtest(b, state, n_randint(state, 10), 100);

        exp = (ulong) n_randtest(state) % 20UL;

        fmpq_poly_pow(a, b, exp);
        fmpq_poly_pow(b, b, exp);

        cflags |= fmpq_poly_is_canonical(a) ? 0 : 1;
        cflags |= fmpq_poly_is_canonical(b) ? 0 : 2;
        result = (fmpq_poly_equal(a, b) && !cflags);
        if (!result)
        {
            printf("FAIL:\n");
            printf("exp = %lu\n", exp);
            printf("a = "), fmpq_poly_debug(a), printf("\n\n");
            printf("b = "), fmpq_poly_debug(b), printf("\n\n");
            printf("cflags = %lu\n\n", cflags);
            abort();
        }

        fmpq_poly_clear(a);
        fmpq_poly_clear(b);
    }

    /* Compare with repeated multiplications by the base */
    for (i = 0; i < 200 * flint_test_multiplier(); i++)
    {
        fmpq_poly_t a, b, c;
        ulong exp;

        fmpq_poly_init(a);
        fmpq_poly_init(b);
        fmpq_poly_init(c);
        fmpq_poly_randtest(b, state, n_randint(state, 10), 100);

        exp = (ulong) n_randtest(state) % 20UL;

        fmpq_poly_pow(a, b, exp);

        if (exp == 0)
        {
            fmpq_poly_set_ui(c, 1);
        }
        else
        {
            ulong j;
            fmpq_poly_set(c, b);

            for (j = 1; j < exp; j++)
                fmpq_poly_mul(c, c, b);
        }

        cflags |= fmpq_poly_is_canonical(a) ? 0 : 1;
        result = (fmpq_poly_equal(a, c) && !cflags);
        if (!result)
        {
            printf("FAIL:\n");
            printf("exp = %lu\n", exp);
            printf("a = "), fmpq_poly_debug(a), printf("\n\n");
            printf("c = "), fmpq_poly_debug(c), printf("\n\n");
            printf("cflags = %lu\n\n", cflags);
            abort();
        }

        fmpq_poly_clear(a);
        fmpq_poly_clear(b);
        fmpq_poly_clear(c);
    }

    flint_randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return 0;
}
