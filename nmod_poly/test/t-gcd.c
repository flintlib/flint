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

    Copyright (C) 2011 William Hart
    Copyright (C) 2012 Sebastian Pancratz

******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "nmod_poly.h"
#include "ulong_extras.h"

int
main(void)
{
    int i, result;
    flint_rand_t state;
    flint_randinit(state);

    printf("gcd....");
    fflush(stdout);

    /* 
       Find coprime polys, multiply by another poly 
       and check the GCD is that poly 
    */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        nmod_poly_t a, b, c, g;

        mp_limb_t n;
        do n = n_randtest_not_zero(state);
        while (!n_is_probabprime(n));

        nmod_poly_init(a, n);
        nmod_poly_init(b, n);
        nmod_poly_init(c, n);
        nmod_poly_init(g, n);
        
        do {
            nmod_poly_randtest(a, state, n_randint(state, 1000));
            nmod_poly_randtest(b, state, n_randint(state, 1000));
            nmod_poly_gcd_euclidean(g, a, b);
        } while (g->length != 1);

        do {
            nmod_poly_randtest(c, state, n_randint(state, 1000));
        } while (c->length < 2);
        nmod_poly_make_monic(c, c);
        
        nmod_poly_mul(a, a, c);
        nmod_poly_mul(b, b, c);

        nmod_poly_gcd_euclidean(g, a, b);

        result = (nmod_poly_equal(g, c));
        if (!result)
        {
            printf("FAIL:\n");
            nmod_poly_print(a), printf("\n\n");
            nmod_poly_print(b), printf("\n\n");
            nmod_poly_print(c), printf("\n\n");
            nmod_poly_print(g), printf("\n\n");
            printf("n = %ld\n", n);
            abort();
        }
        
        nmod_poly_clear(a);
        nmod_poly_clear(b);
        nmod_poly_clear(c);
        nmod_poly_clear(g);
    }

    /* Check aliasing of a and g */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        nmod_poly_t a, b, g;

        mp_limb_t n;
        do n = n_randtest(state);
        while (!n_is_probabprime(n));

        nmod_poly_init(a, n);
        nmod_poly_init(b, n);
        nmod_poly_init(g, n);
        nmod_poly_randtest(a, state, n_randint(state, 1000));
        nmod_poly_randtest(b, state, n_randint(state, 1000));
        
        nmod_poly_gcd_euclidean(g, a, b);
        nmod_poly_gcd_euclidean(a, a, b);

        result = (nmod_poly_equal(a, g));
        if (!result)
        {
            printf("FAIL:\n");
            nmod_poly_print(a), printf("\n\n");
            nmod_poly_print(b), printf("\n\n");
            nmod_poly_print(g), printf("\n\n");
            printf("n = %ld\n", n);
            abort();
        }

        nmod_poly_clear(a);
        nmod_poly_clear(b);
        nmod_poly_clear(g);
    }

    /* Check aliasing of b and g */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        nmod_poly_t a, b, g;

        mp_limb_t n;
        do n = n_randtest(state);
        while (!n_is_probabprime(n));

        nmod_poly_init(a, n);
        nmod_poly_init(b, n);
        nmod_poly_init(g, n);
        nmod_poly_randtest(a, state, n_randint(state, 1000));
        nmod_poly_randtest(b, state, n_randint(state, 1000));
       
        nmod_poly_gcd_euclidean(g, a, b);
        nmod_poly_gcd_euclidean(b, a, b);

        result = (nmod_poly_equal(b, g));
        if (!result)
        {
            printf("FAIL:\n");
            nmod_poly_print(a), printf("\n\n");
            nmod_poly_print(b), printf("\n\n");
            nmod_poly_print(g), printf("\n\n");
            printf("n = %ld\n", n);
            abort();
        }

        nmod_poly_clear(a);
        nmod_poly_clear(b);
        nmod_poly_clear(g);
    }

    flint_randclear(state);

    printf("PASS\n");
    return 0;
}
