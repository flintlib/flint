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

    Copyright (C) 2012 Sebastian Pancratz

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
    FLINT_TEST_INIT(state);

    flint_printf("gcdinv....");
    fflush(stdout);

    

    /* Generic case, most likely co-prime arguments ******************************/

    /* Compare with result from XGCD */
    for (i = 0; i < 1000; i++)
    {
        fmpz_t p;
        fmpz_mod_poly_t a, b, d, g, s, t, u;

        fmpz_init(p);
        fmpz_set_ui(p, n_randtest_prime(state, 0));

        fmpz_mod_poly_init(a, p);
        fmpz_mod_poly_init(b, p);
        fmpz_mod_poly_init(d, p);
        fmpz_mod_poly_init(g, p);
        fmpz_mod_poly_init(s, p);
        fmpz_mod_poly_init(t, p);
        fmpz_mod_poly_init(u, p);

        fmpz_mod_poly_randtest(a, state, n_randint(state, 100));
        do 
            fmpz_mod_poly_randtest(b, state, n_randint(state, 100));
        while (b->length < 2);

        fmpz_mod_poly_gcdinv(d, u, a, b);
        fmpz_mod_poly_xgcd(g, s, t, a, b);

        result = ((fmpz_mod_poly_equal(d, g) && fmpz_mod_poly_equal(u, s))
                  || (fmpz_mod_poly_is_zero(d)));
        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("a = "), fmpz_mod_poly_print(a), flint_printf("\n\n");
            flint_printf("b = "), fmpz_mod_poly_print(b), flint_printf("\n\n");
            flint_printf("d = "), fmpz_mod_poly_print(d), flint_printf("\n\n");
            flint_printf("g = "), fmpz_mod_poly_print(g), flint_printf("\n\n");
            flint_printf("s = "), fmpz_mod_poly_print(s), flint_printf("\n\n");
            flint_printf("t = "), fmpz_mod_poly_print(t), flint_printf("\n\n");
            flint_printf("u = "), fmpz_mod_poly_print(u), flint_printf("\n\n");
            abort();
        }

        fmpz_mod_poly_clear(a);
        fmpz_mod_poly_clear(b);
        fmpz_mod_poly_clear(d);
        fmpz_mod_poly_clear(g);
        fmpz_mod_poly_clear(s);
        fmpz_mod_poly_clear(t);
        fmpz_mod_poly_clear(u);
        fmpz_clear(p);
    }

    /* Special case, arguments share a factor ********************************/

    /* Compare with result from XGCD */
    for (i = 0; i < 1000; i++)
    {
        fmpz_t p;
        fmpz_mod_poly_t a, b, d, f, g, s, t, u;

        fmpz_init(p);
        fmpz_set_ui(p, n_randtest_prime(state, 0));

        fmpz_mod_poly_init(a, p);
        fmpz_mod_poly_init(b, p);
        fmpz_mod_poly_init(d, p);
        fmpz_mod_poly_init(f, p);
        fmpz_mod_poly_init(g, p);
        fmpz_mod_poly_init(s, p);
        fmpz_mod_poly_init(t, p);
        fmpz_mod_poly_init(u, p);

        fmpz_mod_poly_randtest(a, state, n_randint(state, 100));
        do 
            fmpz_mod_poly_randtest(b, state, n_randint(state, 100));
        while (b->length < 2);
        fmpz_mod_poly_randtest_not_zero(f, state, n_randint(state, 20) + 1);
        fmpz_mod_poly_mul(a, f, a);
        fmpz_mod_poly_mul(b, f, b);

        fmpz_mod_poly_gcdinv(d, u, a, b);
        fmpz_mod_poly_xgcd(g, s, t, a, b);

        result = ((fmpz_mod_poly_equal(d, g) && fmpz_mod_poly_equal(u, s))
                  || (fmpz_mod_poly_is_zero(d)));
        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("a = "), fmpz_mod_poly_print(a), flint_printf("\n\n");
            flint_printf("b = "), fmpz_mod_poly_print(b), flint_printf("\n\n");
            flint_printf("d = "), fmpz_mod_poly_print(d), flint_printf("\n\n");
            flint_printf("f = "), fmpz_mod_poly_print(f), flint_printf("\n\n");
            flint_printf("g = "), fmpz_mod_poly_print(g), flint_printf("\n\n");
            flint_printf("s = "), fmpz_mod_poly_print(s), flint_printf("\n\n");
            flint_printf("t = "), fmpz_mod_poly_print(t), flint_printf("\n\n");
            flint_printf("u = "), fmpz_mod_poly_print(u), flint_printf("\n\n");
            abort();
        }

        fmpz_mod_poly_clear(a);
        fmpz_mod_poly_clear(b);
        fmpz_mod_poly_clear(d);
        fmpz_mod_poly_clear(f);
        fmpz_mod_poly_clear(g);
        fmpz_mod_poly_clear(s);
        fmpz_mod_poly_clear(t);
        fmpz_mod_poly_clear(u);
        fmpz_clear(p);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}

