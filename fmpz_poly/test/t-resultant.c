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

    printf("resultant....");
    fflush(stdout);

    flint_randinit(state);

    /* Just one specific test */
    {
        fmpz_poly_t f, g;
        fmpz_t a, b;
        fmpz_poly_init(f);
        fmpz_poly_init(g);
        fmpz_init(a);
        fmpz_init(b);
        fmpz_poly_set_str(f, "11  -15 -2 -2 17 0 0 6 0 -5 1 -1");
        fmpz_poly_set_str(g, "9  2 1 1 1 1 1 0 -1 -2");
        fmpz_poly_resultant(a, f, g);
        fmpz_set_str(b, "-44081924855067", 10);

        result = (fmpz_equal(a, b));
        if (!result)
        {
            printf("FAIL:\n");
            printf("f = "), fmpz_poly_print(f), printf("\n\n");
            printf("g = "), fmpz_poly_print(g), printf("\n\n");
            printf("res(f, h)  = "), fmpz_print(a), printf("\n\n");
            abort();
        }

        fmpz_poly_clear(f);
        fmpz_poly_clear(g);
        fmpz_clear(a);
        fmpz_clear(b);
    }

    /* Check that R(fg, h) = R(f, h) R(g, h) */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fmpz_t a, b, c, d;
        fmpz_poly_t f, g, h, p;

        fmpz_init(a);
        fmpz_init(b);
        fmpz_init(c);
        fmpz_init(d);
        fmpz_poly_init(f);
        fmpz_poly_init(g);
        fmpz_poly_init(h);
        fmpz_poly_init(p);
        fmpz_poly_randtest(f, state, n_randint(state, 50), 100);
        fmpz_poly_randtest(g, state, n_randint(state, 50), 100);
        fmpz_poly_randtest(h, state, n_randint(state, 10), 100);

        fmpz_poly_resultant(a, f, h);
        fmpz_poly_resultant(b, g, h);
        fmpz_mul(c, a, b);
        fmpz_poly_mul(p, f, g);
        fmpz_poly_resultant(d, p, h);

        result = (fmpz_equal(c, d));
        if (!result)
        {
            printf("FAIL:\n");
            printf("f = "), fmpz_poly_print(f), printf("\n\n");
            printf("g = "), fmpz_poly_print(g), printf("\n\n");
            printf("h = "), fmpz_poly_print(h), printf("\n\n");
            printf("res(f, h)  = "), fmpz_print(a), printf("\n\n");
            printf("res(g, h)  = "), fmpz_print(b), printf("\n\n");
            printf("res(fg, h) = "), fmpz_print(d), printf("\n\n");
            abort();
        }

        fmpz_clear(a);
        fmpz_clear(b);
        fmpz_clear(c);
        fmpz_clear(d);
        fmpz_poly_clear(f);
        fmpz_poly_clear(g);
        fmpz_poly_clear(h);
        fmpz_poly_clear(p);
    }

    flint_randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return 0;
}
