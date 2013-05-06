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

    printf("gcd_euclidean....");
    fflush(stdout);

    flint_randinit(state);

    /* Generic case, most likely co-prime arguments ******************************/

    /* Check aliasing of a and c */
    for (i = 0; i < 1000; i++)
    {
        fmpz_t p;
        fmpz_mod_poly_t a, b, c;

        fmpz_init(p);
        fmpz_set_ui(p, n_randtest_prime(state, 0));

        fmpz_mod_poly_init(a, p);
        fmpz_mod_poly_init(b, p);
        fmpz_mod_poly_init(c, p);
        fmpz_mod_poly_randtest(a, state, n_randint(state, 100));
        fmpz_mod_poly_randtest(b, state, n_randint(state, 100));

        fmpz_mod_poly_gcd_euclidean(c, a, b);
        fmpz_mod_poly_gcd_euclidean(a, a, b);

        result = (fmpz_mod_poly_equal(a, c));
        if (!result)
        {
            printf("FAIL:\n");
            fmpz_mod_poly_print(a), printf("\n\n");
            fmpz_mod_poly_print(b), printf("\n\n");
            fmpz_mod_poly_print(c), printf("\n\n");
            abort();
        }

        fmpz_mod_poly_clear(a);
        fmpz_mod_poly_clear(b);
        fmpz_mod_poly_clear(c);
        fmpz_clear(p);
    }

    /* Check aliasing of b and c */
    for (i = 0; i < 1000; i++)
    {
        fmpz_t p;
        fmpz_mod_poly_t a, b, c;

        fmpz_init(p);
        fmpz_set_ui(p, n_randtest_prime(state, 0));

        fmpz_mod_poly_init(a, p);
        fmpz_mod_poly_init(b, p);
        fmpz_mod_poly_init(c, p);
        fmpz_mod_poly_randtest(a, state, n_randint(state, 100));
        fmpz_mod_poly_randtest(b, state, n_randint(state, 100));

        fmpz_mod_poly_gcd_euclidean(c, a, b);
        fmpz_mod_poly_gcd_euclidean(b, a, b);

        result = (fmpz_mod_poly_equal(b, c));
        if (!result)
        {
            printf("FAIL:\n");
            fmpz_mod_poly_print(a), printf("\n\n");
            fmpz_mod_poly_print(b), printf("\n\n");
            fmpz_mod_poly_print(c), printf("\n\n");
            abort();
        }

        fmpz_mod_poly_clear(a);
        fmpz_mod_poly_clear(b);
        fmpz_mod_poly_clear(c);
        fmpz_clear(p);
    }

    /* 
        Check that g = GCD(a,b) divides a and b, 
        and that 1 == GCD(a/g, b/g)
     */
    for (i = 0; i < 500; i++)
    {
        fmpz_t p;
        fmpz_mod_poly_t a, b, c, d, g, h, s, t;

        fmpz_init(p);
        fmpz_set_ui(p, n_randtest_prime(state, 0));

        fmpz_mod_poly_init(a, p);
        fmpz_mod_poly_init(b, p);
        fmpz_mod_poly_init(c, p);
        fmpz_mod_poly_init(d, p);
        fmpz_mod_poly_init(g, p);
        fmpz_mod_poly_init(h, p);
        fmpz_mod_poly_init(s, p);
        fmpz_mod_poly_init(t, p);
        fmpz_mod_poly_randtest(a, state, n_randint(state, 100));
        fmpz_mod_poly_randtest(b, state, n_randint(state, 100));

        fmpz_mod_poly_gcd_euclidean(g, a, b);

        if (fmpz_mod_poly_is_zero(g))
        {
            result = 1;
        }
        else
        {
            fmpz_mod_poly_divrem_basecase(c, s, a, g);
            fmpz_mod_poly_divrem_basecase(d, t, b, g);
            fmpz_mod_poly_gcd_euclidean(h, c, d);

            result = (fmpz_mod_poly_is_zero(s) && fmpz_mod_poly_is_zero(t) 
                      && (h->length == 1) && (fmpz_is_one(h->coeffs)));
        }

        if (!result)
        {
            printf("FAIL:\n");
            printf("p = "), fmpz_print(p), printf("\n\n");
            printf("a = "), fmpz_mod_poly_print(a), printf("\n\n");
            printf("b = "), fmpz_mod_poly_print(b), printf("\n\n");
            printf("c = "), fmpz_mod_poly_print(c), printf("\n\n");
            printf("d = "), fmpz_mod_poly_print(d), printf("\n\n");
            printf("g = "), fmpz_mod_poly_print(g), printf("\n\n");
            printf("h = "), fmpz_mod_poly_print(h), printf("\n\n");
            printf("s = "), fmpz_mod_poly_print(s), printf("\n\n");
            printf("t = "), fmpz_mod_poly_print(t), printf("\n\n");
            abort();
        }

        fmpz_mod_poly_clear(a);
        fmpz_mod_poly_clear(b);
        fmpz_mod_poly_clear(c);
        fmpz_mod_poly_clear(d);
        fmpz_mod_poly_clear(g);
        fmpz_mod_poly_clear(h);
        fmpz_mod_poly_clear(s);
        fmpz_mod_poly_clear(t);
        fmpz_clear(p);
    }

    /* Special case, arguments share a factor ********************************/

    /* Check aliasing of a and c */
    for (i = 0; i < 1000; i++)
    {
        fmpz_t p;
        fmpz_mod_poly_t a, b, c, f;

        fmpz_init(p);
        fmpz_set_ui(p, n_randtest_prime(state, 0));

        fmpz_mod_poly_init(a, p);
        fmpz_mod_poly_init(b, p);
        fmpz_mod_poly_init(c, p);
        fmpz_mod_poly_init(f, p);
        fmpz_mod_poly_randtest(a, state, n_randint(state, 100));
        fmpz_mod_poly_randtest(b, state, n_randint(state, 100));
        fmpz_mod_poly_randtest(f, state, n_randint(state, 20));
        fmpz_mod_poly_mul(a, a, f);
        fmpz_mod_poly_mul(b, b, f);

        fmpz_mod_poly_gcd_euclidean(c, a, b);
        fmpz_mod_poly_gcd_euclidean(a, a, b);

        result = (fmpz_mod_poly_equal(a, c));
        if (!result)
        {
            printf("FAIL:\n");
            fmpz_mod_poly_print(a), printf("\n\n");
            fmpz_mod_poly_print(b), printf("\n\n");
            fmpz_mod_poly_print(c), printf("\n\n");
            abort();
        }

        fmpz_mod_poly_clear(a);
        fmpz_mod_poly_clear(b);
        fmpz_mod_poly_clear(c);
        fmpz_mod_poly_clear(f);
        fmpz_clear(p);
    }

    /* Check aliasing of b and c */
    for (i = 0; i < 1000; i++)
    {
        fmpz_t p;
        fmpz_mod_poly_t a, b, c, f;

        fmpz_init(p);
        fmpz_set_ui(p, n_randtest_prime(state, 0));

        fmpz_mod_poly_init(a, p);
        fmpz_mod_poly_init(b, p);
        fmpz_mod_poly_init(c, p);
        fmpz_mod_poly_init(f, p);
        fmpz_mod_poly_randtest(a, state, n_randint(state, 100));
        fmpz_mod_poly_randtest(b, state, n_randint(state, 100));
        fmpz_mod_poly_randtest(f, state, n_randint(state, 20));
        fmpz_mod_poly_mul(a, a, f);
        fmpz_mod_poly_mul(b, b, f);

        fmpz_mod_poly_gcd_euclidean(c, a, b);
        fmpz_mod_poly_gcd_euclidean(b, a, b);

        result = (fmpz_mod_poly_equal(b, c));
        if (!result)
        {
            printf("FAIL:\n");
            fmpz_mod_poly_print(a), printf("\n\n");
            fmpz_mod_poly_print(b), printf("\n\n");
            fmpz_mod_poly_print(c), printf("\n\n");
            abort();
        }

        fmpz_mod_poly_clear(a);
        fmpz_mod_poly_clear(b);
        fmpz_mod_poly_clear(c);
        fmpz_mod_poly_clear(f);
        fmpz_clear(p);
    }

    /* 
        Check that g = GCD(a,b) divides a and b, 
        and that 1 == GCD(a/g, b/g)
     */
    for (i = 0; i < 500; i++)
    {
        fmpz_t p;
        fmpz_mod_poly_t a, b, c, d, f, g, h, s, t;

        fmpz_init(p);
        fmpz_set_ui(p, n_randtest_prime(state, 0));

        fmpz_mod_poly_init(a, p);
        fmpz_mod_poly_init(b, p);
        fmpz_mod_poly_init(c, p);
        fmpz_mod_poly_init(d, p);
        fmpz_mod_poly_init(f, p);
        fmpz_mod_poly_init(g, p);
        fmpz_mod_poly_init(h, p);
        fmpz_mod_poly_init(s, p);
        fmpz_mod_poly_init(t, p);
        fmpz_mod_poly_randtest(a, state, n_randint(state, 100));
        fmpz_mod_poly_randtest(b, state, n_randint(state, 100));
        fmpz_mod_poly_randtest(f, state, n_randint(state, 20));
        fmpz_mod_poly_mul(a, a, f);
        fmpz_mod_poly_mul(b, b, f);

        fmpz_mod_poly_gcd_euclidean(g, a, b);

        if (fmpz_mod_poly_is_zero(g))
        {
            result = 1;
        }
        else
        {
            fmpz_mod_poly_divrem_basecase(c, s, a, g);
            fmpz_mod_poly_divrem_basecase(d, t, b, g);
            fmpz_mod_poly_gcd_euclidean(h, c, d);

            result = (fmpz_mod_poly_is_zero(s) && fmpz_mod_poly_is_zero(t) 
                      && (h->length == 1) && (fmpz_is_one(h->coeffs)));
        }

        if (!result)
        {
            printf("FAIL:\n");
            printf("p = "), fmpz_print(p), printf("\n\n");
            printf("a = "), fmpz_mod_poly_print(a), printf("\n\n");
            printf("b = "), fmpz_mod_poly_print(b), printf("\n\n");
            printf("c = "), fmpz_mod_poly_print(c), printf("\n\n");
            printf("d = "), fmpz_mod_poly_print(d), printf("\n\n");
            printf("g = "), fmpz_mod_poly_print(g), printf("\n\n");
            printf("h = "), fmpz_mod_poly_print(h), printf("\n\n");
            printf("s = "), fmpz_mod_poly_print(s), printf("\n\n");
            printf("t = "), fmpz_mod_poly_print(t), printf("\n\n");
            abort();
        }

        fmpz_mod_poly_clear(a);
        fmpz_mod_poly_clear(b);
        fmpz_mod_poly_clear(c);
        fmpz_mod_poly_clear(d);
        fmpz_mod_poly_clear(f);
        fmpz_mod_poly_clear(g);
        fmpz_mod_poly_clear(h);
        fmpz_mod_poly_clear(s);
        fmpz_mod_poly_clear(t);
        fmpz_clear(p);
    }

    flint_randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return 0;
}

