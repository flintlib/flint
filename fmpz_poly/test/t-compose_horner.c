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

    printf("compose_horner....");
    fflush(stdout);

    flint_randinit(state);

    /* Check aliasing of the first argument */
    for (i = 0; i < 8 * flint_test_multiplier(); i++)
    {
        fmpz_poly_t f, g, h;

        fmpz_poly_init(f);
        fmpz_poly_init(g);
        fmpz_poly_init(h);
        fmpz_poly_randtest(g, state, n_randint(state, 40), 80);
        fmpz_poly_randtest(h, state, n_randint(state, 20), 50);

        fmpz_poly_compose_horner(f, g, h);
        fmpz_poly_compose_horner(g, g, h);

        result = (fmpz_poly_equal(f, g));
        if (!result)
        {
            printf("FAIL:\n");
            fmpz_poly_print(f), printf("\n\n");
            fmpz_poly_print(g), printf("\n\n");
            abort();
        }

        fmpz_poly_clear(f);
        fmpz_poly_clear(g);
        fmpz_poly_clear(h);
    }

    /* Check aliasing of the second argument */
    for (i = 0; i < 8 * flint_test_multiplier(); i++)
    {
        fmpz_poly_t f, g, h;

        fmpz_poly_init(f);
        fmpz_poly_init(g);
        fmpz_poly_init(h);
        fmpz_poly_randtest(g, state, n_randint(state, 40), 80);
        fmpz_poly_randtest(h, state, n_randint(state, 20), 50);

        fmpz_poly_compose_horner(f, g, h);
        fmpz_poly_compose_horner(h, g, h);

        result = (fmpz_poly_equal(f, h));
        if (!result)
        {
            printf("FAIL:\n");
            fmpz_poly_print(f), printf("\n\n");
            fmpz_poly_print(h), printf("\n\n");
            abort();
        }

        fmpz_poly_clear(f);
        fmpz_poly_clear(g);
        fmpz_poly_clear(h);
    }

    /* Compare with the default method */
    for (i = 0; i < 8 * flint_test_multiplier(); i++)
    {
        fmpz_poly_t f1, f2, g, h;

        fmpz_poly_init(f1);
        fmpz_poly_init(f2);
        fmpz_poly_init(g);
        fmpz_poly_init(h);
        fmpz_poly_randtest(g, state, n_randint(state, 40), 80);
        fmpz_poly_randtest(h, state, n_randint(state, 20), 50);
        
        fmpz_poly_compose_horner(f1, g, h);
        fmpz_poly_compose(f2, g, h);

        result = (fmpz_poly_equal(f1, f2));
        if (!result)
        {
            printf("FAIL:\n");
            fmpz_poly_print(f1), printf("\n\n");
            fmpz_poly_print(f2), printf("\n\n");
            abort();
        }

        fmpz_poly_clear(f1);
        fmpz_poly_clear(f2);
        fmpz_poly_clear(g);
        fmpz_poly_clear(h);
    }

    flint_randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return 0;
}
