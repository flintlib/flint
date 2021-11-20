/*
    Copyright (C) 2013 Mike Hansen

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "nmod_poly.h"
#include "ulong_extras.h"

int
main(void)
{
    int i, result;
    
    FLINT_TEST_INIT(state);

    printf("invmod....");
    fflush(stdout);

    /* Test aliasing *************************************************************/

    /* Aliasing c and a */
    for (i = 0; i < 500; i++)
    {
        mp_limb_t p;
        nmod_poly_t a, b, c;
        int ans1, ans2;

        p = n_randtest_prime(state, 0);

        nmod_poly_init(a, p);
        nmod_poly_init(b, p);
        nmod_poly_init(c, p);

        do 
            nmod_poly_randtest(b, state, n_randint(state, 100));
        while (b->length < 3);
        nmod_poly_randtest(a, state, n_randint(state, 100));

        ans1 = nmod_poly_invmod(c, a, b);
        ans2 = nmod_poly_invmod(a, a, b);

        result = (ans1 == ans2 && nmod_poly_equal(a, c));
        if (!result)
        {
            printf("FAIL (alias a and c):\n");
            nmod_poly_print(a), printf("\n\n");
            nmod_poly_print(b), printf("\n\n");
            nmod_poly_print(c), printf("\n\n");
            printf("ans1 = %d\n\n", ans1);
            printf("ans2 = %d\n\n", ans2);
            fflush(stdout);
            flint_abort();
        }

        nmod_poly_clear(a);
        nmod_poly_clear(b);
        nmod_poly_clear(c);
    }

    /* Aliasing c and b */
    for (i = 0; i < 500; i++)
    {
        mp_limb_t p;
        nmod_poly_t a, b, c;
        int ans1, ans2;

        p = n_randtest_prime(state, 0);

        nmod_poly_init(a, p);
        nmod_poly_init(b, p);
        nmod_poly_init(c, p);

        do 
            nmod_poly_randtest(b, state, n_randint(state, 100));
        while (b->length < 3);
        nmod_poly_randtest(a, state, n_randint(state, 100));

        ans1 = nmod_poly_invmod(c, a, b);
        ans2 = nmod_poly_invmod(b, a, b);

        result = ((ans1 == ans2) && nmod_poly_equal(b, c));
        if (!result)
        {
            printf("FAIL (alias b and c):\n");
            nmod_poly_print(a), printf("\n\n");
            nmod_poly_print(b), printf("\n\n");
            nmod_poly_print(c), printf("\n\n");
            printf("ans1 = %d\n\n", ans1);
            printf("ans2 = %d\n\n", ans2);
            fflush(stdout);
            flint_abort();
        }

        nmod_poly_clear(a);
        nmod_poly_clear(b);
        nmod_poly_clear(c);
    }

    /* Compare with result from XGCD */
    for (i = 0; i < 1000; i++)
    {
        mp_limb_t p;
        nmod_poly_t a, b, g, s, t, u;
        int ans;

        p = n_randtest_prime(state, 0);

        nmod_poly_init(a, p);
        nmod_poly_init(b, p);
        nmod_poly_init(g, p);
        nmod_poly_init(s, p);
        nmod_poly_init(t, p);
        nmod_poly_init(u, p);

        do 
            nmod_poly_randtest(b, state, n_randint(state, 100));
        while (b->length < 3);
        nmod_poly_randtest(a, state, n_randint(state, 100));

        ans = nmod_poly_invmod(u, a, b);
        nmod_poly_xgcd(g, s, t, a, b);

        result = (((ans) && g->length == 1 
                   && g->coeffs[0] == WORD(1) && nmod_poly_equal(s, u)) 
                 || (!(ans) && g->length > 1));

        if (!result)
        {
            printf("FAIL:\n");
            nmod_poly_print(a), printf("\n\n");
            nmod_poly_print(b), printf("\n\n");
            nmod_poly_print(g), printf("\n\n");
            nmod_poly_print(s), printf("\n\n");
            nmod_poly_print(t), printf("\n\n");
            nmod_poly_print(u), printf("\n\n");
            printf("ans = %d\n\n", ans);
            fflush(stdout);
            flint_abort();
        }

        nmod_poly_clear(a);
        nmod_poly_clear(b);
        nmod_poly_clear(g);
        nmod_poly_clear(s);
        nmod_poly_clear(t);
        nmod_poly_clear(u);
    }

    /* Special case, arguments share a factor ********************************/

    /* Check correctness */
    for (i = 0; i < 1000; i++)
    {
        mp_limb_t p;
        nmod_poly_t a, b, f, u;
        int ans;

        p = n_randtest_prime(state, 0);

        nmod_poly_init(a, p);
        nmod_poly_init(b, p);
        nmod_poly_init(f, p);
        nmod_poly_init(u, p);

        do 
            nmod_poly_randtest(b, state, n_randint(state, 100));
        while (b->length < 2);
        nmod_poly_randtest(a, state, n_randint(state, 100));
        do 
            nmod_poly_randtest_not_zero(f, state, n_randint(state, 20) + 1);
        while (f->length < 2);
        nmod_poly_mul(a, f, a);
        nmod_poly_mul(b, f, b);

        ans = nmod_poly_invmod(u, a, b);

        result = (!ans);
        if (!result)
        {
            printf("FAIL:\n");
            nmod_poly_print(a), printf("\n\n");
            nmod_poly_print(b), printf("\n\n");
            nmod_poly_print(f), printf("\n\n");
            nmod_poly_print(u), printf("\n\n");
            printf("ans = %d\n\n", ans);
            fflush(stdout);
            flint_abort();
        }

        nmod_poly_clear(a);
        nmod_poly_clear(b);
        nmod_poly_clear(f);
        nmod_poly_clear(u);
    }

    FLINT_TEST_CLEANUP(state);
    printf("PASS\n");
    return 0;
}

