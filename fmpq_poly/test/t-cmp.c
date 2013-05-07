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

    printf("cmp....");
    fflush(stdout);

    flint_randinit(state);

    /* Check aliasing */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpq_poly_t f;

        fmpq_poly_init(f);
        fmpq_poly_randtest(f, state, n_randint(state, 100), 200);

        result = (fmpq_poly_cmp(f, f) == 0);
        if (!result)
        {
            printf("FAIL:\n");
            fmpq_poly_debug(f), printf("\n");
            abort();
        }

        fmpq_poly_clear(f);
    }

    /*
       Check transitivity, i.e. f <= g <= h implies f <= h, that is 
       NOT (f <= g <= h) OR f <= h
     */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpq_poly_t f, g, h;

        fmpq_poly_init(f);
        fmpq_poly_init(g);
        fmpq_poly_init(h);
        fmpq_poly_randtest(f, state, n_randint(state, 100), 200);
        fmpq_poly_randtest(g, state, n_randint(state, 100), 200);
        fmpq_poly_randtest(h, state, n_randint(state, 100), 200);

        result = !(fmpq_poly_cmp(f, g) <= 0) || !(fmpq_poly_cmp(g, h) <= 0)
            || (fmpq_poly_cmp(f, h) <= 0);
        if (!result)
        {
            printf("FAIL:\n");
            fmpq_poly_debug(f), printf("\n");
            fmpq_poly_debug(g), printf("\n");
            fmpq_poly_debug(h), printf("\n");
            abort();
        }

        fmpq_poly_clear(f);
        fmpq_poly_clear(g);
        fmpq_poly_clear(h);
    }

    /* Check that <, ==, or > */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpq_poly_t f, g;

        fmpq_poly_init(f);
        fmpq_poly_init(g);
        fmpq_poly_randtest(f, state, n_randint(state, 100), 200);
        fmpq_poly_randtest(g, state, n_randint(state, 100), 200);

        result = (fmpq_poly_cmp(f, g) < 0) || (fmpq_poly_equal(f, g))
            || (fmpq_poly_cmp(f, g) > 0);

        if (!result)
        {
            printf("FAIL:\n");
            fmpq_poly_debug(f), printf("\n");
            fmpq_poly_debug(g), printf("\n");
            printf("cmp(f,g) = %d\n", fmpq_poly_cmp(f, g));
            abort();
        }

        fmpq_poly_clear(f);
        fmpq_poly_clear(g);
    }

    flint_randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return 0;
}
