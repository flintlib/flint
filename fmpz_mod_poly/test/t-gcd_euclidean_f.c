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

    Copyright (C) 2011 Sebastian Pancratz

******************************************************************************/

#include "fmpz.h"
#include "fmpz_mod_poly.h"
#include "ulong_extras.h"

int
main(void)
{
    int i, result;
    flint_rand_t state;

    printf("gcd_euclidean_f....");
    fflush(stdout);

    flint_randinit(state);

    /*
        Compare with the usual GCD function.

        N.B.  I checked by hand that this test shows both outcomes, 
        i.e. trivial and non-trivial factors, sufficiently frequently.
     */
    for (i = 0; i < 1000; i++)
    {
        fmpz_t p, f;
        fmpz_mod_poly_t a, b, c, d;

        fmpz_init(f);
        fmpz_init(p);
        fmpz_randtest_unsigned(p, state, 100);
        fmpz_add_ui(p, p, 2);

        fmpz_mod_poly_init(a, p);
        fmpz_mod_poly_init(b, p);
        fmpz_mod_poly_init(c, p);
        fmpz_mod_poly_init(d, p);
        fmpz_mod_poly_randtest(a, state, n_randint(state, 60));
        fmpz_mod_poly_randtest(b, state, n_randint(state, 60));

        fmpz_mod_poly_gcd_euclidean_f(f, c, a, b);
        if (!fmpz_is_one(f))
        {
            result = 1;
        }
        else
        {
            fmpz_mod_poly_gcd(d, a, b);
            result = fmpz_mod_poly_equal(c, d);
        }

        if (!result)
        {
            printf("FAIL:\n");
            printf("p = "), fmpz_print(p), printf("\n\n");
            printf("a = "), fmpz_mod_poly_print(a), printf("\n\n");
            printf("b = "), fmpz_mod_poly_print(b), printf("\n\n");
            printf("c = "), fmpz_mod_poly_print(c), printf("\n\n");
            printf("d = "), fmpz_mod_poly_print(d), printf("\n\n");
            printf("f = "), fmpz_print(f), printf("\n\n");
            abort();
        }

        fmpz_mod_poly_clear(a);
        fmpz_mod_poly_clear(b);
        fmpz_mod_poly_clear(c);
        fmpz_mod_poly_clear(d);
        fmpz_clear(f);
        fmpz_clear(p);
    }

    flint_randclear(state);
    flint_cleanup();
    printf("PASS\n");
    return 0;
}

