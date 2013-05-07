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
#include "nmod_poly.h"
#include "ulong_extras.h"

int
main(void)
{
    int i, result;
    flint_rand_t state;
    flint_randinit(state);

    printf("scalar_mul_nmod....");
    fflush(stdout);

    /* Check aliasing of a and b */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        nmod_poly_t a, b;
        mp_limb_t n = n_randtest_not_zero(state);
        mp_limb_t c = n_randint(state, n);

        nmod_poly_init(a, n);
        nmod_poly_init(b, n);
        nmod_poly_randtest(a, state, n_randint(state, 100));

        nmod_poly_scalar_mul_nmod(b, a, c);
        nmod_poly_scalar_mul_nmod(a, a, c);

        result = (nmod_poly_equal(a, b));
        if (!result)
        {
            printf("FAIL:\n");
            nmod_poly_print(a), printf("\n\n");
            nmod_poly_print(b), printf("\n\n");
            abort();
        }

        nmod_poly_clear(a);
        nmod_poly_clear(b);
    }

    /* Check (a + b)*c = a*c + b*c */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        nmod_poly_t a, b, d1, d2;
        mp_limb_t n = n_randtest_not_zero(state);
        mp_limb_t c = n_randint(state, n);

        nmod_poly_init(a, n);
        nmod_poly_init(b, n);
        nmod_poly_init(d1, n);
        nmod_poly_init(d2, n);
        nmod_poly_randtest(a, state, n_randint(state, 100));
        nmod_poly_randtest(b, state, n_randint(state, 100));

        nmod_poly_add(d1, a, b);
        nmod_poly_scalar_mul_nmod(d1, d1, c);

        nmod_poly_scalar_mul_nmod(d2, a, c);
        nmod_poly_scalar_mul_nmod(b, b, c);
        nmod_poly_add(d2, d2, b);

        result = (nmod_poly_equal(d1, d2));
        if (!result)
        {
            printf("FAIL:\n");
            nmod_poly_print(d1), printf("\n\n");
            nmod_poly_print(d2), printf("\n\n");
            abort();
        }

        nmod_poly_clear(a);
        nmod_poly_clear(b);
        nmod_poly_clear(d1);
        nmod_poly_clear(d2);
    }

    flint_randclear(state);

    printf("PASS\n");
    return 0;
}
