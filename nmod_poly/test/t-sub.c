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

    Copyright (C) 2010 William Hart

******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "nmod_vec.h"
#include "nmod_poly.h"
#include "ulong_extras.h"

int
main(void)
{
    int i, result;
    flint_rand_t state;
    flint_randinit(state);

    printf("sub....");
    fflush(stdout);

    /* Check a - b = a + neg(b) */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        nmod_poly_t a, b, c, d;
        mp_limb_t n = n_randtest_not_zero(state);

        nmod_poly_init(a, n);
        nmod_poly_init(b, n);
        nmod_poly_init(c, n);
        nmod_poly_init(d, n);
        nmod_poly_randtest(a, state, n_randint(state, 100));
        nmod_poly_randtest(b, state, n_randint(state, 100));

        nmod_poly_sub(c, a, b);
        nmod_poly_neg(b, b);
        nmod_poly_add(d, a, b);

        result = (nmod_poly_equal(d, d));
        if (!result)
        {
            printf("FAIL:\n");
            nmod_poly_print(a), printf("\n\n");
            nmod_poly_print(b), printf("\n\n");
            nmod_poly_print(c), printf("\n\n");
            nmod_poly_print(d), printf("\n\n");
            abort();
        }

        nmod_poly_clear(a);
        nmod_poly_clear(b);
        nmod_poly_clear(c);
        nmod_poly_clear(d);
    }

    /* Check aliasing of a and c */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        nmod_poly_t a, b, c;
        mp_limb_t n = n_randtest_not_zero(state);

        nmod_poly_init(a, n);
        nmod_poly_init(b, n);
        nmod_poly_init(c, n);
        nmod_poly_randtest(a, state, n_randint(state, 100));
        nmod_poly_randtest(b, state, n_randint(state, 100));

        nmod_poly_sub(c, a, b);
        nmod_poly_sub(a, a, b);

        result = (nmod_poly_equal(a, c));
        if (!result)
        {
            printf("FAIL:\n");
            nmod_poly_print(a), printf("\n\n");
            nmod_poly_print(b), printf("\n\n");
            nmod_poly_print(c), printf("\n\n");
            abort();
        }

        nmod_poly_clear(a);
        nmod_poly_clear(b);
        nmod_poly_clear(c);
    }

    /* Check aliasing of b and c */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        nmod_poly_t a, b, c;
        mp_limb_t n = n_randtest_not_zero(state);

        nmod_poly_init(a, n);
        nmod_poly_init(b, n);
        nmod_poly_init(c, n);
        nmod_poly_randtest(a, state, n_randint(state, 100));
        nmod_poly_randtest(b, state, n_randint(state, 100));

        nmod_poly_sub(c, a, b);
        nmod_poly_sub(b, a, b);

        result = (nmod_poly_equal(b, c));
        if (!result)
        {
            printf("FAIL:\n");
            nmod_poly_print(a), printf("\n\n");
            nmod_poly_print(b), printf("\n\n");
            nmod_poly_print(c), printf("\n\n");
            abort();
        }

        nmod_poly_clear(a);
        nmod_poly_clear(b);
        nmod_poly_clear(c);
    }

    flint_randclear(state);

    printf("PASS\n");
    return 0;
}
