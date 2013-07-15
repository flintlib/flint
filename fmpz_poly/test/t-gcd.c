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

    printf("gcd....");
    fflush(stdout);

    flint_randinit(state);

    /* Check aliasing of a and b */
    for (i = 0; i < 50 * flint_test_multiplier(); i++)
    {
        fmpz_poly_t a, b, c;

        fmpz_poly_init(a);
        fmpz_poly_init(b);
        fmpz_poly_init(c);
        fmpz_poly_randtest(b, state, n_randint(state, 40), 80);
        fmpz_poly_randtest(c, state, n_randint(state, 40), 80);

        fmpz_poly_gcd(a, b, c);
        fmpz_poly_gcd(b, b, c);

        result = (fmpz_poly_equal(a, b));
        if (!result)
        {
            printf("FAIL (aliasing a and b):\n");
            fmpz_poly_print(a), printf("\n\n");
            fmpz_poly_print(b), printf("\n\n");
            abort();
        }

        fmpz_poly_clear(a);
        fmpz_poly_clear(b);
        fmpz_poly_clear(c);
    }

    /* Check aliasing of a and c */
    for (i = 0; i < 50 * flint_test_multiplier(); i++)
    {
        fmpz_poly_t a, b, c;

        fmpz_poly_init(a);
        fmpz_poly_init(b);
        fmpz_poly_init(c);
        fmpz_poly_randtest(b, state, n_randint(state, 40), 80);
        fmpz_poly_randtest(c, state, n_randint(state, 40), 80);

        fmpz_poly_gcd(a, b, c);
        fmpz_poly_gcd(c, b, c);

        result = (fmpz_poly_equal(a, c));
        if (!result)
        {
            printf("FAIL (aliasing a and c):\n");
            fmpz_poly_print(a), printf("\n\n");
            fmpz_poly_print(c), printf("\n\n");
            abort();
        }

        fmpz_poly_clear(a);
        fmpz_poly_clear(b);
        fmpz_poly_clear(c);
    }

    /* Check that a divides GCD(af, ag) */
    for (i = 0; i < 50 * flint_test_multiplier(); i++)
    {
        fmpz_poly_t a, d, f, g, q, r;

        fmpz_poly_init(a);
        fmpz_poly_init(d);
        fmpz_poly_init(f);
        fmpz_poly_init(g);
        fmpz_poly_init(q);
        fmpz_poly_init(r);
        fmpz_poly_randtest_not_zero(a, state, n_randint(state, 24) + 1, 24);
        fmpz_poly_randtest(f, state, n_randint(state, 40), 80);
        fmpz_poly_randtest(g, state, n_randint(state, 40), 80);

        fmpz_poly_mul(f, a, f);
        fmpz_poly_mul(g, a, g);
        fmpz_poly_gcd(d, f, g);

        fmpz_poly_divrem_divconquer(q, r, d, a);

        result = (r->length == 0L);
        if (!result)
        {
            printf("FAIL (check a | gcd(af, ag)):\n");
            fmpz_poly_print(f), printf("\n");
            fmpz_poly_print(g), printf("\n");
            fmpz_poly_print(a), printf("\n");
            fmpz_poly_print(d), printf("\n");
            abort();
        }

        fmpz_poly_clear(a);
        fmpz_poly_clear(d);
        fmpz_poly_clear(f);
        fmpz_poly_clear(g);
        fmpz_poly_clear(q);
        fmpz_poly_clear(r);
    }

    flint_randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return 0;
}
