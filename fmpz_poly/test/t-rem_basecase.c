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

    printf("rem_basecase....");
    fflush(stdout);

    flint_randinit(state);

    /* Compare with full division, no aliasing */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fmpz_poly_t a, b, q, r, r2;

        fmpz_poly_init(a);
        fmpz_poly_init(b);
        fmpz_poly_init(q);
        fmpz_poly_init(r);
        fmpz_poly_init(r2);

        fmpz_poly_randtest(a, state, n_randint(state, 50), 100);
        fmpz_poly_randtest_not_zero(b, state, n_randint(state, a->length + 1) + 1, 100);

        fmpz_poly_divrem_basecase(q, r, a, b);
        fmpz_poly_rem_basecase(r2, a, b);

        result = (fmpz_poly_equal(r, r2));
        if (!result)
        {
            printf("FAIL:\n");
            fmpz_poly_print(a), printf("\n\n");
            fmpz_poly_print(b), printf("\n\n");
            fmpz_poly_print(q), printf("\n\n");
            fmpz_poly_print(r), printf("\n\n");
            fmpz_poly_print(r2), printf("\n\n");
            abort();
        }

        fmpz_poly_clear(a);
        fmpz_poly_clear(b);
        fmpz_poly_clear(q);
        fmpz_poly_clear(r);
        fmpz_poly_clear(r2);
    }

    /* Check r and a alias */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fmpz_poly_t a, b, r;

        fmpz_poly_init(a);
        fmpz_poly_init(b);
        fmpz_poly_init(r);
        fmpz_poly_randtest(a, state, n_randint(state, 50), 100);
        fmpz_poly_randtest_not_zero(b, state, n_randint(state, 50) + 1, 100);

        fmpz_poly_rem_basecase(r, a, b);
        fmpz_poly_rem_basecase(a, a, b);

        result = (fmpz_poly_equal(a, r));
        if (!result)
        {
            printf("FAIL:\n");
            fmpz_poly_print(a), printf("\n\n");
            fmpz_poly_print(r), printf("\n\n");
            abort();
        }

        fmpz_poly_clear(a);
        fmpz_poly_clear(b);
        fmpz_poly_clear(r);
    }

    /* Check r and b alias */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fmpz_poly_t a, b, r;

        fmpz_poly_init(a);
        fmpz_poly_init(b);
        fmpz_poly_init(r);
        fmpz_poly_randtest(a, state, n_randint(state, 50), 100);
        fmpz_poly_randtest_not_zero(b, state, n_randint(state, 50) + 1, 100);

        fmpz_poly_rem_basecase(r, a, b);
        fmpz_poly_rem_basecase(b, a, b);

        result = (fmpz_poly_equal(b, r));
        if (!result)
        {
            printf("FAIL:\n");
            fmpz_poly_print(a), printf("\n\n");
            fmpz_poly_print(r), printf("\n\n");
            abort();
        }

        fmpz_poly_clear(a);
        fmpz_poly_clear(b);
        fmpz_poly_clear(r);
    }


    flint_randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return 0;
}
