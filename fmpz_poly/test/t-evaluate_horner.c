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
#include <mpir.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_poly.h"
#include "ulong_extras.h"

int
main(void)
{
    int i, result;
    fmpz_randstate_t state;

    printf("evaluate_horner....");
    fflush(stdout);

    fmpz_poly_randinit(state);

    /* Check aliasing */
    for (i = 0; i < 10000; i++)
    {
        fmpz_t a, b;
        fmpz_poly_t f;

        fmpz_init(a);
        fmpz_init(b);
        fmpz_poly_init(f);
        fmpz_poly_randtest(f, state, n_randint(100), n_randint(200));
        fmpz_randtest(a, state, n_randint(100));

        fmpz_poly_evaluate_horner(b, f, a);
        fmpz_poly_evaluate_horner(a, f, a);

        result = (fmpz_equal(a, b));
        if (!result)
        {
            printf("FAIL:\n");
            fmpz_print(a), printf("\n\n");
            fmpz_print(b), printf("\n\n");
            abort();
        }

        fmpz_clear(a);
        fmpz_clear(b);
        fmpz_poly_clear(f);
    }

    /* Check that (f+g)(a) = f(a) + g(a) */
    for (i = 0; i < 10000; i++)
    {
        fmpz_t a, b, c;
        fmpz_poly_t f, g;

        fmpz_init(a);
        fmpz_init(b);
        fmpz_init(c);
        fmpz_poly_init(f);
        fmpz_poly_init(g);
        fmpz_poly_randtest(f, state, n_randint(100), n_randint(200));
        fmpz_poly_randtest(g, state, n_randint(100), n_randint(200));
        fmpz_randtest(a, state, n_randint(100));

        fmpz_poly_evaluate_horner(b, f, a);
        fmpz_poly_evaluate_horner(c, g, a);
        fmpz_add(b, b, c);
        fmpz_poly_add(f, f, g);
        fmpz_poly_evaluate_horner(c, f, a);

        result = (fmpz_equal(b, c));
        if (!result)
        {
            printf("FAIL:\n");
            fmpz_print(b), printf("\n\n");
            fmpz_print(c), printf("\n\n");
            abort();
        }

        fmpz_clear(a);
        fmpz_clear(b);
        fmpz_clear(c);
        fmpz_poly_clear(f);
        fmpz_poly_clear(g);
    }

    fmpz_poly_randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return 0;
}
