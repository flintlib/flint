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
    Copyright (C) 2010, 2011 Sebastian Pancratz

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

    printf("lcm....");
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

        fmpz_poly_lcm(a, b, c);
        fmpz_poly_lcm(b, b, c);

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

        fmpz_poly_lcm(a, b, c);
        fmpz_poly_lcm(c, b, c);

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

    /* Check that GCD(f, g) LCM(f, g) == f g */
    for (i = 0; i < 50 * flint_test_multiplier(); i++)
    {
        fmpz_poly_t f, g, gcd, lcm, lhs, rhs;

        fmpz_poly_init(f);
        fmpz_poly_init(g);
        fmpz_poly_init(gcd);
        fmpz_poly_init(lcm);
        fmpz_poly_init(lhs);
        fmpz_poly_init(rhs);

        fmpz_poly_randtest(f, state, n_randint(state, 40), 80);
        fmpz_poly_randtest(g, state, n_randint(state, 40), 80);

        fmpz_poly_gcd(gcd, f, g);
        fmpz_poly_lcm(lcm, f, g);
        fmpz_poly_mul(lhs, gcd, lcm);
        fmpz_poly_mul(rhs, f, g);
        if (!fmpz_poly_is_zero(rhs) && fmpz_sgn(fmpz_poly_lead(rhs)) < 0)
            fmpz_poly_neg(rhs, rhs);

        result = (fmpz_poly_equal(lhs, rhs));
        if (!result)
        {
            printf("FAIL (GCD(f, g) * LCM(f, g) == f * g):\n");
            fmpz_poly_print(f), printf("\n");
            fmpz_poly_print(g), printf("\n");
            fmpz_poly_print(gcd), printf("\n");
            fmpz_poly_print(lcm), printf("\n");
            fmpz_poly_print(lhs), printf("\n");
            fmpz_poly_print(rhs), printf("\n");
            abort();
        }

        fmpz_poly_clear(f);
        fmpz_poly_clear(g);
        fmpz_poly_clear(gcd);
        fmpz_poly_clear(lcm);
        fmpz_poly_clear(lhs);
        fmpz_poly_clear(rhs);
    }

    flint_randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}

