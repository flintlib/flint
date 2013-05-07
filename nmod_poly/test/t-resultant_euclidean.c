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
#include "nmod_poly.h"
#include "ulong_extras.h"

int
main(void)
{
    int i, result;
    flint_rand_t state;
    flint_randinit(state);

    printf("resultant_euclidean....");
    fflush(stdout);

    /* Check res(f, g) == (-1)^(deg f deg g) res(g, f) */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        nmod_poly_t f, g;
        mp_limb_t x, y;
        mp_limb_t n;

        do n = n_randtest_not_zero(state);
        while (!n_is_probabprime(n));

        nmod_poly_init(f, n);
        nmod_poly_init(g, n);
        
        nmod_poly_randtest(f, state, n_randint(state, 200));
        nmod_poly_randtest(g, state, n_randint(state, 200));

        x = nmod_poly_resultant_euclidean(f, g);
        y = nmod_poly_resultant_euclidean(g, f);

        if ((nmod_poly_degree(f) * nmod_poly_degree(g)) % 2)
            y = nmod_neg(y, f->mod);

        result = (x == y);
        if (!result)
        {
            printf("FAIL (res(f, g) == (-1)^(deg f deg g) res(g, f)):\n");
            nmod_poly_print(f), printf("\n\n");
            nmod_poly_print(g), printf("\n\n");
            printf("x = %lu\n", x);
            printf("y = %lu\n", y);
            printf("n = %lu\n", n);
            abort();
        }
        
        nmod_poly_clear(f);
        nmod_poly_clear(g);
    }

    /* Check res(f h, g) == res(f, g) res(h, g) */
    for (i = 0; i < 50 * flint_test_multiplier(); i++)
    {
        nmod_poly_t f, g, h;
        mp_limb_t x, y, z;
        mp_limb_t n;

        do n = n_randtest_not_zero(state);
        while (!n_is_probabprime(n));

        nmod_poly_init(f, n);
        nmod_poly_init(g, n);
        nmod_poly_init(h, n);
        
        nmod_poly_randtest(f, state, n_randint(state, 200));
        nmod_poly_randtest(g, state, n_randint(state, 200));
        nmod_poly_randtest(h, state, n_randint(state, 200));

        y = nmod_poly_resultant_euclidean(f, g);
        z = nmod_poly_resultant_euclidean(h, g);
        y = nmod_mul(y, z, f->mod);
        nmod_poly_mul(f, f, h);
        x = nmod_poly_resultant_euclidean(f, g);

        result = (x == y);
        if (!result)
        {
            printf("FAIL (res(f h, g) == res(f, g) res(h, g)):\n");
            nmod_poly_print(f), printf("\n\n");
            nmod_poly_print(g), printf("\n\n");
            nmod_poly_print(h), printf("\n\n");
            printf("x = %lu\n", x);
            printf("y = %lu\n", y);
            printf("z = %ld\n", z);
            printf("n = %lu\n", n);
            abort();
        }
        
        nmod_poly_clear(f);
        nmod_poly_clear(g);
        nmod_poly_clear(h);
    }

    flint_randclear(state);

    printf("PASS\n");
    return 0;
}
