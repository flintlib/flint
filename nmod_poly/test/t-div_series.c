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

    printf("div_series....");
    fflush(stdout);

    /* Check A/B * B = A */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        nmod_poly_t q, a, b, prod;
        len_t m;

        mp_limb_t n;
        do n = n_randtest_not_zero(state);
        while (!n_is_probabprime(n));

        nmod_poly_init(prod, n);
        nmod_poly_init(a, n);
        nmod_poly_init(b, n);
        nmod_poly_init(q, n);
        
        nmod_poly_randtest(a, state, n_randint(state, 2000));
        do nmod_poly_randtest(b, state, n_randint(state, 2000));
        while (b->length == 0 || b->coeffs[0] == 0);

        m = n_randint(state, 2000) + 1;

        nmod_poly_div_series(q, a, b, m);
        nmod_poly_mullow(prod, q, b, m);
        nmod_poly_truncate(a, m);

        result = (nmod_poly_equal(a, prod));
        if (!result)
        {
            printf("FAIL:\n");
            nmod_poly_print(q), printf("\n\n");
            nmod_poly_print(b), printf("\n\n");
            nmod_poly_print(a), printf("\n\n");
            nmod_poly_print(prod), printf("\n\n");
            printf("n = %ld\n", n);
            abort();
        }
        
        nmod_poly_clear(q);
        nmod_poly_clear(a);
        nmod_poly_clear(b);
        nmod_poly_clear(prod);
    }

    /* Check aliasing of q and a */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        nmod_poly_t q, a, b;
        len_t m;

        mp_limb_t n;
        do n = n_randtest(state);
        while (!n_is_probabprime(n));

        nmod_poly_init(q, n);
        nmod_poly_init(a, n);
        nmod_poly_init(b, n);

        nmod_poly_randtest(a, state, n_randint(state, 1000));
        do nmod_poly_randtest(b, state, n_randint(state, 1000));
        while (b->length == 0 || b->coeffs[0] == 0);

        m = n_randint(state, 1000) + 1;

        nmod_poly_div_series(q, a, b, m);
        nmod_poly_div_series(a, a, b, m);
        
        result = (nmod_poly_equal(q, a));
        if (!result)
        {
            printf("FAIL:\n");
            nmod_poly_print(b), printf("\n\n");
            nmod_poly_print(q), printf("\n\n");
            nmod_poly_print(a), printf("\n\n");
            printf("n = %ld, m = %ld\n", n, m);
            abort();
        }

        nmod_poly_clear(q);
        nmod_poly_clear(a);
        nmod_poly_clear(b);
    }

    /* Check aliasing of q and b */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        nmod_poly_t q, a, b;
        len_t m;

        mp_limb_t n;
        do n = n_randtest(state);
        while (!n_is_probabprime(n));

        nmod_poly_init(q, n);
        nmod_poly_init(a, n);
        nmod_poly_init(b, n);

        nmod_poly_randtest(a, state, n_randint(state, 1000));
        do nmod_poly_randtest(b, state, n_randint(state, 1000));
        while (b->length == 0 || b->coeffs[0] == 0);

        m = n_randint(state, 1000) + 1;

        nmod_poly_div_series(q, a, b, m);
        nmod_poly_div_series(b, a, b, m);
        
        result = (nmod_poly_equal(q, b));
        if (!result)
        {
            printf("FAIL:\n");
            nmod_poly_print(a), printf("\n\n");
            nmod_poly_print(q), printf("\n\n");
            nmod_poly_print(b), printf("\n\n");
            printf("n = %ld, m = %ld\n", n, m);
            abort();
        }

        nmod_poly_clear(q);
        nmod_poly_clear(a);
        nmod_poly_clear(b);
    }

    flint_randclear(state);

    printf("PASS\n");
    return 0;
}
