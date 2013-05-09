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

    Copyright (C) 2012 Fredrik Johansson

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
    int i;
    flint_rand_t state;

    printf("sqrt... ");
    fflush(stdout);

    flint_randinit(state);

    /* Test aliasing */
    for (i = 0; i < 200 * flint_test_multiplier(); i++)
    {
        fmpz_poly_t a, b;
        int square1, square2;

        fmpz_poly_init(a);
        fmpz_poly_init(b);

        fmpz_poly_randtest(a, state, 1 + n_randint(state, 20),
            1 + n_randint(state, 200));

        if (n_randint(state, 2))
            fmpz_poly_sqr(a, a);

        square1 = fmpz_poly_sqrt(b, a);
        square2 = fmpz_poly_sqrt(a, a);

        if ((square1 != square2) || (square1 && !fmpz_poly_equal(a, b)))
        {
            printf("FAIL: aliasing:\n");
            printf("square1 = %d, square2 = %d\n\n", square1, square2);
            printf("a: "); fmpz_poly_print(a); printf("\n\n");
            printf("b: "); fmpz_poly_print(b); printf("\n\n");
            abort();
        }

        fmpz_poly_clear(a);
        fmpz_poly_clear(b);
    }

    /* Test random squares */
    for (i = 0; i < 200 * flint_test_multiplier(); i++)
    {
        fmpz_poly_t a, b, c;
        int square;

        fmpz_poly_init(a);
        fmpz_poly_init(b);
        fmpz_poly_init(c);

        fmpz_poly_randtest(a, state, 1 + n_randint(state, 20),
            1 + n_randint(state, 200));
        fmpz_poly_sqr(b, a);
        square = fmpz_poly_sqrt(c, b);

        if (!square)
        {
            printf("FAIL: square reported nonsquare:\n");
            printf("a: "); fmpz_poly_print(a); printf("\n\n");
            printf("b: "); fmpz_poly_print(b); printf("\n\n");
            printf("c: "); fmpz_poly_print(c); printf("\n\n");
            abort();
        }

        if (!fmpz_poly_is_zero(c) &&
            fmpz_sgn(fmpz_poly_get_coeff_ptr(c, fmpz_poly_degree(c))) < 0)
        {
            printf("FAIL: leading coefficient not positive:\n");
            printf("a: "); fmpz_poly_print(a); printf("\n\n");
            printf("b: "); fmpz_poly_print(b); printf("\n\n");
            printf("c: "); fmpz_poly_print(c); printf("\n\n");
            abort();
        }

        fmpz_poly_sqr(c, c);
        if (!fmpz_poly_equal(c, b))
        {
            printf("FAIL: sqrt(b)^2 != b:\n");
            printf("a: "); fmpz_poly_print(a); printf("\n\n");
            printf("b: "); fmpz_poly_print(b); printf("\n\n");
            printf("c: "); fmpz_poly_print(c); printf("\n\n");
            abort();
        }

        fmpz_poly_clear(a);
        fmpz_poly_clear(b);
        fmpz_poly_clear(c);
    }

    /* Test "almost" squares */
    for (i = 0; i < 200 * flint_test_multiplier(); i++)
    {
        fmpz_poly_t a, b, c;
        fmpz_t t;
        len_t j;
        int square;

        fmpz_poly_init(a);
        fmpz_poly_init(b);
        fmpz_poly_init(c);
        fmpz_init(t);

        fmpz_poly_randtest_not_zero(a, state, 1 + n_randint(state, 20),
            1 + n_randint(state, 200));
        fmpz_poly_sqr(b, a);

        j = n_randint(state, fmpz_poly_length(b));
        fmpz_randtest_not_zero(t, state, 1 + n_randint(state, 100));
        fmpz_add(b->coeffs + j, b->coeffs + j, t);
        _fmpz_poly_normalise(b);

        square = fmpz_poly_sqrt(c, b);

        if (square)
        {
            fmpz_poly_sqr(c, c);
            if (!fmpz_poly_equal(c, b))
            {
                printf("FAIL: sqrt(b)^2 != b:\n");
                printf("a: "); fmpz_poly_print(a); printf("\n\n");
                printf("b: "); fmpz_poly_print(b); printf("\n\n");
                printf("c: "); fmpz_poly_print(c); printf("\n\n");
                abort();
            }
        }

        fmpz_clear(t);
        fmpz_poly_clear(a);
        fmpz_poly_clear(b);
        fmpz_poly_clear(c);
    }

    flint_randclear(state);
    printf("PASS\n");
    _fmpz_cleanup();
    return 0;
}
