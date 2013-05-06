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

    printf("inv_series_newton....");
    fflush(stdout);

    flint_randinit(state);

    /* Check Q^{-1} * Q is congruent 1 mod t^n */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fmpz_poly_t a, b, c, one;
        long n = n_randint(state, 80) + 1;

        fmpz_poly_init(a);
        fmpz_poly_init(b);
        fmpz_poly_init(c);
        fmpz_poly_init(one);

        fmpz_poly_randtest_not_zero(a, state, n_randint(state, 80) + 1, 100);
        fmpz_poly_set_coeff_si(a, 0, n_randint(state, 2) ? 1 : -1);

        fmpz_poly_set_ui(one, 1);

        fmpz_poly_inv_series_newton(b, a, n);
        fmpz_poly_mullow(c, a, b, n);

        result = (fmpz_poly_equal(c, one));
        if (!result)
        {
            printf("FAIL:\n");
            printf("a = "), fmpz_poly_print(a), printf("\n\n");
            printf("b = "), fmpz_poly_print(b), printf("\n\n");
            printf("c = "), fmpz_poly_print(c), printf("\n\n");
            abort();
        }

        fmpz_poly_clear(a);
        fmpz_poly_clear(b);
        fmpz_poly_clear(c);
        fmpz_poly_clear(one);
    }

    /* Verify bug fix for the case Q = -1 mod (x) */
    {
        fmpz_poly_t a, b, c, one;
        long n = 1;

        fmpz_poly_init(a);
        fmpz_poly_init(b);
        fmpz_poly_init(c);
        fmpz_poly_init(one);

        fmpz_poly_set_si(a, -1);
        fmpz_poly_set_ui(one, 1);

        fmpz_poly_inv_series_newton(b, a, n);
        fmpz_poly_mullow(c, a, b, n);

        result = (fmpz_poly_equal(c, one));
        if (!result)
        {
            printf("FAIL:\n");
            printf("a = "), fmpz_poly_print(a), printf("\n\n");
            printf("b = "), fmpz_poly_print(b), printf("\n\n");
            printf("c = "), fmpz_poly_print(c), printf("\n\n");
            abort();
        }

        fmpz_poly_clear(a);
        fmpz_poly_clear(b);
        fmpz_poly_clear(c);
        fmpz_poly_clear(one);
    }

    flint_randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return 0;
}
