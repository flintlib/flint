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
    Copyright (C) 2011 Fredrik Johansson

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

    printf("compose_series....");
    fflush(stdout);

    flint_randinit(state);

    /* Check aliasing of the first argument */
    for (i = 0; i < 10 * flint_test_multiplier(); i++)
    {
        fmpz_poly_t f, g, h;
        len_t n;

        fmpz_poly_init(f);
        fmpz_poly_init(g);
        fmpz_poly_init(h);
        fmpz_poly_randtest(g, state, n_randint(state, 40), 80);
        fmpz_poly_randtest(h, state, n_randint(state, 20), 50);
        fmpz_poly_set_coeff_ui(h, 0, 0);
        n = n_randint(state, 20);

        fmpz_poly_compose_series(f, g, h, n);
        fmpz_poly_compose_series(g, g, h, n);

        result = (fmpz_poly_equal(f, g));
        if (!result)
        {
            printf("FAIL (aliasing 1):\n");
            fmpz_poly_print(f), printf("\n\n");
            fmpz_poly_print(g), printf("\n\n");
            abort();
        }

        fmpz_poly_clear(f);
        fmpz_poly_clear(g);
        fmpz_poly_clear(h);
    }

    /* Check aliasing of the second argument */
    for (i = 0; i < 10 * flint_test_multiplier(); i++)
    {
        fmpz_poly_t f, g, h;
        len_t n;

        fmpz_poly_init(f);
        fmpz_poly_init(g);
        fmpz_poly_init(h);
        fmpz_poly_randtest(g, state, n_randint(state, 40), 80);
        fmpz_poly_randtest(h, state, n_randint(state, 20), 50);
        fmpz_poly_set_coeff_ui(h, 0, 0);
        n = n_randint(state, 20);

        fmpz_poly_compose_series(f, g, h, n);
        fmpz_poly_compose_series(h, g, h, n);

        result = (fmpz_poly_equal(f, h));
        if (!result)
        {
            printf("FAIL (aliasing 2):\n");
            fmpz_poly_print(f), printf("\n\n");
            fmpz_poly_print(h), printf("\n\n");
            abort();
        }

        fmpz_poly_clear(f);
        fmpz_poly_clear(g);
        fmpz_poly_clear(h);
    }

    /* Compare with compose */
    for (i = 0; i < 10 * flint_test_multiplier(); i++)
    {
        fmpz_poly_t f, g, h, s, t;
        len_t n;

        fmpz_poly_init(f);
        fmpz_poly_init(g);
        fmpz_poly_init(h);
        fmpz_poly_init(s);
        fmpz_poly_init(t);
        fmpz_poly_randtest(g, state, n_randint(state, 40), 80);
        fmpz_poly_randtest(h, state, n_randint(state, 20), 50);
        fmpz_poly_set_coeff_ui(h, 0, 0);
        n = n_randint(state, 10);

        fmpz_poly_compose(s, g, h);
        fmpz_poly_truncate(s, n);
        fmpz_poly_compose_series(f, g, h, n);

        result = (fmpz_poly_equal(f, s));
        if (!result)
        {
            printf("FAIL (comparison):\n");
            printf("n = %ld\n", n);
            printf("g = "), fmpz_poly_print(g), printf("\n\n");
            printf("h = "), fmpz_poly_print(h), printf("\n\n");
            printf("f = "), fmpz_poly_print(f), printf("\n\n");
            printf("s = "), fmpz_poly_print(s), printf("\n\n");
            abort();
        }

        fmpz_poly_clear(f);
        fmpz_poly_clear(g);
        fmpz_poly_clear(h);
        fmpz_poly_clear(s);
        fmpz_poly_clear(t);
    }

    flint_randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return 0;
}
