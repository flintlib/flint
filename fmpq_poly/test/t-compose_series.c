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
#include "fmpq_poly.h"
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
        fmpq_poly_t f, g, h;
        long n;

        fmpq_poly_init(f);
        fmpq_poly_init(g);
        fmpq_poly_init(h);
        fmpq_poly_randtest(g, state, n_randint(state, 40), 80);
        fmpq_poly_randtest(h, state, n_randint(state, 20), 50);
        fmpq_poly_set_coeff_ui(h, 0, 0);
        n = n_randint(state, 20);

        fmpq_poly_compose_series(f, g, h, n);
        fmpq_poly_compose_series(g, g, h, n);

        result = (fmpq_poly_equal(f, g));
        if (!result)
        {
            printf("FAIL (aliasing 1):\n");
            fmpq_poly_print(f), printf("\n\n");
            fmpq_poly_print(g), printf("\n\n");
            abort();
        }

        fmpq_poly_clear(f);
        fmpq_poly_clear(g);
        fmpq_poly_clear(h);
    }

    /* Check aliasing of the second argument */
    for (i = 0; i < 10 * flint_test_multiplier(); i++)
    {
        fmpq_poly_t f, g, h;
        long n;

        fmpq_poly_init(f);
        fmpq_poly_init(g);
        fmpq_poly_init(h);
        fmpq_poly_randtest(g, state, n_randint(state, 40), 80);
        fmpq_poly_randtest(h, state, n_randint(state, 20), 50);
        fmpq_poly_set_coeff_ui(h, 0, 0);
        n = n_randint(state, 20);

        fmpq_poly_compose_series(f, g, h, n);
        fmpq_poly_compose_series(h, g, h, n);

        result = (fmpq_poly_equal(f, h));
        if (!result)
        {
            printf("FAIL (aliasing 2):\n");
            fmpq_poly_print(f), printf("\n\n");
            fmpq_poly_print(h), printf("\n\n");
            abort();
        }

        fmpq_poly_clear(f);
        fmpq_poly_clear(g);
        fmpq_poly_clear(h);
    }

    /* Compare with compose */
    for (i = 0; i < 20 * flint_test_multiplier(); i++)
    {
        fmpq_poly_t f, g, h, s, t;
        long n;

        fmpq_poly_init(f);
        fmpq_poly_init(g);
        fmpq_poly_init(h);
        fmpq_poly_init(s);
        fmpq_poly_init(t);
        fmpq_poly_randtest(g, state, n_randint(state, 40), 80);
        fmpq_poly_randtest(h, state, n_randint(state, 20), 50);
        fmpq_poly_set_coeff_ui(h, 0, 0);
        n = n_randint(state, 20);

        fmpq_poly_compose(s, g, h);
        fmpq_poly_truncate(s, n);
        fmpq_poly_compose_series(f, g, h, n);

        result = (fmpq_poly_equal(f, s));
        if (!result)
        {
            printf("FAIL (comparison):\n");
            printf("n = %ld\n", n);
            printf("g = "), fmpq_poly_print(g), printf("\n\n");
            printf("h = "), fmpq_poly_print(h), printf("\n\n");
            printf("f = "), fmpq_poly_print(f), printf("\n\n");
            printf("s = "), fmpq_poly_print(s), printf("\n\n");
            abort();
        }

        fmpq_poly_clear(f);
        fmpq_poly_clear(g);
        fmpq_poly_clear(h);
        fmpq_poly_clear(s);
        fmpq_poly_clear(t);
    }

    flint_randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return 0;
}
