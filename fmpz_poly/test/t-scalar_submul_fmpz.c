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
#include "fmpz_poly.h"
#include "ulong_extras.h"

int
main(void)
{
    int i, result;
    flint_rand_t state;

    printf("scalar_submul_fmpz....");
    fflush(stdout);

    flint_randinit(state);

    /* Check aliasing of a and b */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpz_poly_t a, b, c;
        fmpz_t x;

        fmpz_init(x);
        fmpz_poly_init(a);
        fmpz_poly_init(b);
        fmpz_poly_init(c);
        fmpz_randtest(x, state, n_randint(state, 100));
        fmpz_poly_randtest(a, state, n_randint(state, 100), 200);
        fmpz_poly_set(b, a);
        fmpz_poly_set(c, a);

        fmpz_poly_scalar_submul_fmpz(b, a, x);
        fmpz_poly_scalar_submul_fmpz(a, a, x);

        result = (fmpz_poly_equal(a, b));
        if (!result)
        {
            printf("FAIL (1):\n");
            printf("a = "), fmpz_poly_print(a), printf("\n\n");
            printf("b = "), fmpz_poly_print(b), printf("\n\n");
            printf("c = "), fmpz_poly_print(c), printf("\n\n");
            printf("x = "), fmpz_print(x), printf("\n\n");
            abort();
        }

        fmpz_clear(x);
        fmpz_poly_clear(a);
        fmpz_poly_clear(b);
        fmpz_poly_clear(c);
    }

    /* Check that b += x*a is the same as c = b + x*a */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpz_poly_t a, b, c;
        fmpz_t x;

        fmpz_init(x);
        fmpz_poly_init(a);
        fmpz_poly_init(b);
        fmpz_poly_init(c);
        fmpz_randtest(x, state, n_randint(state, 100));
        fmpz_poly_randtest(a, state, n_randint(state, 100), 200);
        fmpz_poly_randtest(b, state, n_randint(state, 100), 200);

        fmpz_poly_scalar_mul_fmpz(c, a, x);
        fmpz_poly_sub(c, b, c);

        fmpz_poly_scalar_submul_fmpz(b, a, x);

        result = (fmpz_poly_equal(b, c));
        if (!result)
        {
            printf("FAIL (2):\n");
            printf("a = "), fmpz_poly_print(a), printf("\n\n");
            printf("b = "), fmpz_poly_print(b), printf("\n\n");
            printf("c = "), fmpz_poly_print(c), printf("\n\n");
            printf("x = "), fmpz_print(x), printf("\n\n");
            abort();
        }

        fmpz_clear(x);
        fmpz_poly_clear(a);
        fmpz_poly_clear(b);
        fmpz_poly_clear(c);
    }

    flint_randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return 0;
}
