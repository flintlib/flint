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
    flint_rand_t state;

    printf("xgcd....");
    fflush(stdout);

    flint_randinit(state);

    /* Generic case, most likely co-prime arguments ******************************/

    /* Compare with result from GCD and check correctness */
    for (i = 0; i < 1000; i++)
    {
        fmpz_t p;
        fmpz_mod_poly_t a, b, d, g, s, t, v, w;

        fmpz_init(p);
        fmpz_set_ui(p, n_randtest_prime(state, 0));

        fmpz_mod_poly_init(a, p);
        fmpz_mod_poly_init(b, p);
        fmpz_mod_poly_init(d, p);
        fmpz_mod_poly_init(g, p);
        fmpz_mod_poly_init(s, p);
        fmpz_mod_poly_init(t, p);
        fmpz_mod_poly_init(v, p);
        fmpz_mod_poly_init(w, p);
        fmpz_mod_poly_randtest(a, state, n_randint(state, 100));
        fmpz_mod_poly_randtest(b, state, n_randint(state, 100));

        fmpz_mod_poly_gcd(d, a, b);
        fmpz_mod_poly_xgcd(g, s, t, a, b);

        fmpz_mod_poly_mul(v, s, a);
        fmpz_mod_poly_mul(w, t, b);
        fmpz_mod_poly_add(w, v, w);

        result = (fmpz_mod_poly_equal(d, g) && fmpz_mod_poly_equal(g, w));
        if (!result)
        {
            printf("FAIL:\n");
            fmpz_mod_poly_print(a), printf("\n\n");
            fmpz_mod_poly_print(b), printf("\n\n");
            fmpz_mod_poly_print(d), printf("\n\n");
            fmpz_mod_poly_print(g), printf("\n\n");
            fmpz_mod_poly_print(s), printf("\n\n");
            fmpz_mod_poly_print(t), printf("\n\n");
            fmpz_mod_poly_print(v), printf("\n\n");
            fmpz_mod_poly_print(w), printf("\n\n");
            abort();
        }

        fmpz_mod_poly_clear(a);
        fmpz_mod_poly_clear(b);
        fmpz_mod_poly_clear(d);
        fmpz_mod_poly_clear(g);
        fmpz_mod_poly_clear(s);
        fmpz_mod_poly_clear(t);
        fmpz_mod_poly_clear(v);
        fmpz_mod_poly_clear(w);
        fmpz_clear(p);
    }

    /* Special case, arguments share a factor ********************************/

    /* Compare with result from GCD and check correctness */
    for (i = 0; i < 1000; i++)
    {
        fmpz_t p;
        fmpz_mod_poly_t a, b, d, f, g, s, t, v, w;

        fmpz_init(p);
        fmpz_set_ui(p, n_randtest_prime(state, 0));

        fmpz_mod_poly_init(a, p);
        fmpz_mod_poly_init(b, p);
        fmpz_mod_poly_init(d, p);
        fmpz_mod_poly_init(f, p);
        fmpz_mod_poly_init(g, p);
        fmpz_mod_poly_init(s, p);
        fmpz_mod_poly_init(t, p);
        fmpz_mod_poly_init(v, p);
        fmpz_mod_poly_init(w, p);
        fmpz_mod_poly_randtest(a, state, n_randint(state, 100));
        fmpz_mod_poly_randtest(b, state, n_randint(state, 100));
        fmpz_mod_poly_randtest(f, state, n_randint(state, 20));
        fmpz_mod_poly_mul(a, a, f);
        fmpz_mod_poly_mul(b, b, f);

        fmpz_mod_poly_gcd(d, a, b);
        fmpz_mod_poly_xgcd(g, s, t, a, b);

        fmpz_mod_poly_mul(v, s, a);
        fmpz_mod_poly_mul(w, t, b);
        fmpz_mod_poly_add(w, v, w);

        result = (fmpz_mod_poly_equal(d, g) && fmpz_mod_poly_equal(g, w));
        if (!result)
        {
            printf("FAIL:\n");
            fmpz_mod_poly_print(a), printf("\n\n");
            fmpz_mod_poly_print(b), printf("\n\n");
            fmpz_mod_poly_print(d), printf("\n\n");
            fmpz_mod_poly_print(f), printf("\n\n");
            fmpz_mod_poly_print(g), printf("\n\n");
            fmpz_mod_poly_print(s), printf("\n\n");
            fmpz_mod_poly_print(t), printf("\n\n");
            fmpz_mod_poly_print(v), printf("\n\n");
            fmpz_mod_poly_print(w), printf("\n\n");
            abort();
        }

        fmpz_mod_poly_clear(a);
        fmpz_mod_poly_clear(b);
        fmpz_mod_poly_clear(d);
        fmpz_mod_poly_clear(f);
        fmpz_mod_poly_clear(g);
        fmpz_mod_poly_clear(s);
        fmpz_mod_poly_clear(t);
        fmpz_mod_poly_clear(v);
        fmpz_mod_poly_clear(w);
        fmpz_clear(p);
    }

    flint_randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return 0;
}

