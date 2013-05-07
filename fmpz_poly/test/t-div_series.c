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
#include "fmpz.h"
#include "fmpz_poly.h"
#include "ulong_extras.h"

int
main(void)
{
    int i, result;
    flint_rand_t state;

    printf("div_series....");
    fflush(stdout);

    flint_randinit(state);

    /* Check aliasing q and a */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fmpz_poly_t a, b, q;
        long n = n_randint(state, 50) + 1;

        fmpz_poly_init(a);
        fmpz_poly_init(b);
        fmpz_poly_init(q);

        fmpz_poly_randtest(a, state, n_randint(state, 50) + 1, 100);
        fmpz_poly_randtest_not_zero(b, state, n_randint(state, 50) + 1, 100);
        fmpz_poly_set_coeff_ui(b, 0, 1);

        fmpz_poly_div_series(q, a, b, n);
        fmpz_poly_div_series(a, a, b, n);

        result = (fmpz_poly_equal(q, a));
        if (!result)
        {
            printf("FAIL (alias q and a):\n");
            printf("a = "), fmpz_poly_print(a), printf("\n\n");
            printf("b = "), fmpz_poly_print(b), printf("\n\n");
            printf("q = "), fmpz_poly_print(q), printf("\n\n");
            abort();
        }

        fmpz_poly_clear(a);
        fmpz_poly_clear(b);
        fmpz_poly_clear(q);
    }

    /* Check aliasing q and b */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fmpz_poly_t a, b, q;
        long n = n_randint(state, 50) + 1;

        fmpz_poly_init(a);
        fmpz_poly_init(b);
        fmpz_poly_init(q);

        fmpz_poly_randtest(a, state, n_randint(state, 50) + 1, 100);
        fmpz_poly_randtest_not_zero(b, state, n_randint(state, 50) + 1, 100);
        fmpz_poly_set_coeff_ui(b, 0, 1);

        fmpz_poly_div_series(q, a, b, n);
        fmpz_poly_div_series(b, a, b, n);

        result = (fmpz_poly_equal(q, b));
        if (!result)
        {
            printf("FAIL (alias q and b):\n");
            printf("a = "), fmpz_poly_print(a), printf("\n\n");
            printf("b = "), fmpz_poly_print(b), printf("\n\n");
            printf("q = "), fmpz_poly_print(q), printf("\n\n");
            abort();
        }

        fmpz_poly_clear(a);
        fmpz_poly_clear(b);
        fmpz_poly_clear(q);
    }

    /* Check that Q * B == A */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fmpz_poly_t a, b, p, q;
        long n = n_randint(state, 50) + 1;

        fmpz_poly_init(a);
        fmpz_poly_init(b);
        fmpz_poly_init(p);
        fmpz_poly_init(q);

        fmpz_poly_randtest(a, state, n_randint(state, 50) + 1, 100);
        fmpz_poly_randtest_not_zero(b, state, n_randint(state, 50) + 1, 100);
        fmpz_poly_set_coeff_ui(b, 0, 1);

        fmpz_poly_div_series(q, a, b, n);
        fmpz_poly_mullow(p, q, b, n);

        fmpz_poly_truncate(a, n);

        result = (fmpz_poly_equal(p, a));
        if (!result)
        {
            printf("FAIL (check Q * B = A):\n");
            printf("a = "), fmpz_poly_print(a), printf("\n\n");
            printf("b = "), fmpz_poly_print(b), printf("\n\n");
            printf("p = "), fmpz_poly_print(p), printf("\n\n");
            printf("q = "), fmpz_poly_print(q), printf("\n\n");
            abort();
        }

        fmpz_poly_clear(a);
        fmpz_poly_clear(b);
        fmpz_poly_clear(p);
        fmpz_poly_clear(q);
    }

    flint_randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return 0;
}
