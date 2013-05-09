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
    Copyright (C) 2011 Sebastian Pancratz

******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_mod_poly.h"
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
    for (i = 0; i < 1000; i++)
    {
        fmpz_t p;
        fmpz_mod_poly_t a, b, c, one;
        len_t n = n_randint(state, 80) + 1;

        fmpz_init(p);
        fmpz_randtest_unsigned(p, state, 2 * FLINT_BITS);
        fmpz_add_ui(p, p, 2);

        fmpz_mod_poly_init(a, p);
        fmpz_mod_poly_init(b, p);
        fmpz_mod_poly_init(c, p);
        fmpz_mod_poly_init(one, p);

        fmpz_mod_poly_randtest_not_zero(a, state, n_randint(state, 80) + 1);
        {
            fmpz_t d;

            fmpz_init(d);
            fmpz_gcd(d, a->coeffs, p);
            while (!fmpz_is_one(d))
            {
                fmpz_randm(a->coeffs, state, p);
                fmpz_gcd(d, a->coeffs, p);
            }
            fmpz_clear(d);
        }

        fmpz_mod_poly_set_ui(one, 1);

        fmpz_mod_poly_inv_series_newton(b, a, n);
        fmpz_mod_poly_mullow(c, a, b, n);

        result = (fmpz_mod_poly_equal(c, one));
        if (!result)
        {
            printf("FAIL:\n");
            printf("a = "), fmpz_mod_poly_print(a), printf("\n\n");
            printf("b = "), fmpz_mod_poly_print(b), printf("\n\n");
            printf("c = "), fmpz_mod_poly_print(c), printf("\n\n");
            printf("p = "), fmpz_print(p), printf("\n\n");
            abort();
        }

        fmpz_mod_poly_clear(a);
        fmpz_mod_poly_clear(b);
        fmpz_mod_poly_clear(c);
        fmpz_mod_poly_clear(one);
        fmpz_clear(p);
    }

    flint_randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return 0;
}

