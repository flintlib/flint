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

    printf("invmod....");
    fflush(stdout);

    flint_randinit(state);

    /* Test aliasing *************************************************************/

    /* Aliasing c and a */
    for (i = 0; i < 500; i++)
    {
        fmpz_t p;
        fmpz_mod_poly_t a, b, c;
        int ans1, ans2;

        fmpz_init(p);
        fmpz_set_ui(p, n_randtest_prime(state, 0));

        fmpz_mod_poly_init(a, p);
        fmpz_mod_poly_init(b, p);
        fmpz_mod_poly_init(c, p);

        do 
            fmpz_mod_poly_randtest(b, state, n_randint(state, 100));
        while (b->length < 3);
        fmpz_mod_poly_randtest(a, state, n_randint(state, 100));

        ans1 = fmpz_mod_poly_invmod(c, a, b);
        ans2 = fmpz_mod_poly_invmod(a, a, b);

        result = (ans1 == ans2 && fmpz_mod_poly_equal(a, c));
        if (!result)
        {
            printf("FAIL (alias a and c):\n");
            fmpz_mod_poly_print(a), printf("\n\n");
            fmpz_mod_poly_print(b), printf("\n\n");
            fmpz_mod_poly_print(c), printf("\n\n");
            printf("ans1 = %d\n\n", ans1);
            printf("ans2 = %d\n\n", ans2);
            abort();
        }

        fmpz_mod_poly_clear(a);
        fmpz_mod_poly_clear(b);
        fmpz_mod_poly_clear(c);
        fmpz_clear(p);
    }

    /* Aliasing c and b */
    for (i = 0; i < 500; i++)
    {
        fmpz_t p;
        fmpz_mod_poly_t a, b, c;
        int ans1, ans2;

        fmpz_init(p);
        fmpz_set_ui(p, n_randtest_prime(state, 0));

        fmpz_mod_poly_init(a, p);
        fmpz_mod_poly_init(b, p);
        fmpz_mod_poly_init(c, p);

        do 
            fmpz_mod_poly_randtest(b, state, n_randint(state, 100));
        while (b->length < 3);
        fmpz_mod_poly_randtest(a, state, n_randint(state, 100));

        ans1 = fmpz_mod_poly_invmod(c, a, b);
        ans2 = fmpz_mod_poly_invmod(b, a, b);

        result = ((ans1 == ans2) && fmpz_mod_poly_equal(b, c));
        if (!result)
        {
            printf("FAIL (alias b and c):\n");
            fmpz_mod_poly_print(a), printf("\n\n");
            fmpz_mod_poly_print(b), printf("\n\n");
            fmpz_mod_poly_print(c), printf("\n\n");
            printf("ans1 = %d\n\n", ans1);
            printf("ans2 = %d\n\n", ans2);
            abort();
        }

        fmpz_mod_poly_clear(a);
        fmpz_mod_poly_clear(b);
        fmpz_mod_poly_clear(c);
        fmpz_clear(p);
    }

    /* Compare with result from XGCD */
    for (i = 0; i < 1000; i++)
    {
        fmpz_t p;
        fmpz_mod_poly_t a, b, g, s, t, u;
        int ans;

        fmpz_init(p);
        fmpz_set_ui(p, n_randtest_prime(state, 0));

        fmpz_mod_poly_init(a, p);
        fmpz_mod_poly_init(b, p);
        fmpz_mod_poly_init(g, p);
        fmpz_mod_poly_init(s, p);
        fmpz_mod_poly_init(t, p);
        fmpz_mod_poly_init(u, p);

        do 
            fmpz_mod_poly_randtest(b, state, n_randint(state, 100));
        while (b->length < 3);
        fmpz_mod_poly_randtest(a, state, n_randint(state, 100));

        ans = fmpz_mod_poly_invmod(u, a, b);
        fmpz_mod_poly_xgcd(g, s, t, a, b);

        result = (((ans) && g->length == 1 
                        && fmpz_is_one(g->coeffs) && fmpz_mod_poly_equal(s, u)) 
                 || (!(ans) && g->length > 1));

        if (!result)
        {
            printf("FAIL:\n");
            fmpz_mod_poly_print(a), printf("\n\n");
            fmpz_mod_poly_print(b), printf("\n\n");
            fmpz_mod_poly_print(g), printf("\n\n");
            fmpz_mod_poly_print(s), printf("\n\n");
            fmpz_mod_poly_print(t), printf("\n\n");
            fmpz_mod_poly_print(u), printf("\n\n");
            printf("ans = %d\n\n", ans);
            abort();
        }

        fmpz_mod_poly_clear(a);
        fmpz_mod_poly_clear(b);
        fmpz_mod_poly_clear(g);
        fmpz_mod_poly_clear(s);
        fmpz_mod_poly_clear(t);
        fmpz_mod_poly_clear(u);
        fmpz_clear(p);
    }

    /* Special case, arguments share a factor ********************************/

    /* Check correctness */
    for (i = 0; i < 1000; i++)
    {
        fmpz_t p;
        fmpz_mod_poly_t a, b, f, u;
        int ans;

        fmpz_init(p);
        fmpz_set_ui(p, n_randtest_prime(state, 0));

        fmpz_mod_poly_init(a, p);
        fmpz_mod_poly_init(b, p);
        fmpz_mod_poly_init(f, p);
        fmpz_mod_poly_init(u, p);

        do 
            fmpz_mod_poly_randtest(b, state, n_randint(state, 100));
        while (b->length < 2);
        fmpz_mod_poly_randtest(a, state, n_randint(state, 100));
        do 
            fmpz_mod_poly_randtest_not_zero(f, state, n_randint(state, 20) + 1);
        while (f->length < 2);
        fmpz_mod_poly_mul(a, f, a);
        fmpz_mod_poly_mul(b, f, b);

        ans = fmpz_mod_poly_invmod(u, a, b);

        result = (!ans);
        if (!result)
        {
            printf("FAIL:\n");
            fmpz_mod_poly_print(a), printf("\n\n");
            fmpz_mod_poly_print(b), printf("\n\n");
            fmpz_mod_poly_print(f), printf("\n\n");
            fmpz_mod_poly_print(u), printf("\n\n");
            printf("ans = %d\n\n", ans);
            abort();
        }

        fmpz_mod_poly_clear(a);
        fmpz_mod_poly_clear(b);
        fmpz_mod_poly_clear(f);
        fmpz_mod_poly_clear(u);
        fmpz_clear(p);
    }

    flint_randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return 0;
}

