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

    Copyright (C) 2011 Fredrik Johansson

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
    int iter;
    FLINT_TEST_INIT(state);
    

    flint_printf("inflate....");
    fflush(stdout);

    for (iter = 0; iter < 100 * flint_test_multiplier(); iter++)
    {
        nmod_poly_t poly1, poly2, poly3, xp;
        mp_limb_t modulus;
        ulong inflation;

        modulus = n_randtest_prime(state, 0);

        nmod_poly_init(poly1, modulus);
        nmod_poly_init(poly2, modulus);
        nmod_poly_init(poly3, modulus);
        nmod_poly_init(xp, modulus);

        nmod_poly_randtest(poly1, state, n_randint(state, 20));
        inflation = n_randint(state, 10);

        nmod_poly_inflate(poly2, poly1, inflation);

        nmod_poly_set_coeff_ui(xp, inflation, 1);
        nmod_poly_compose(poly3, poly1, xp);

        if (!nmod_poly_equal(poly2, poly3))
        {
            flint_printf("FAIL: not equal to compose (inflation = %wu)\n", inflation);
            flint_printf("poly1:\n"); nmod_poly_print(poly1); flint_printf("\n\n");
            flint_printf("poly2:\n"); nmod_poly_print(poly2); flint_printf("\n\n");
            flint_printf("poly3:\n"); nmod_poly_print(poly3); flint_printf("\n\n");
            abort();
        }

        nmod_poly_inflate(poly1, poly1, inflation);
        if (!nmod_poly_equal(poly1, poly2))
        {
            flint_printf("FAIL: aliasing (inflation = %wu)\n", inflation);
            flint_printf("poly1:\n"); nmod_poly_print(poly1); flint_printf("\n\n");
            flint_printf("poly2:\n"); nmod_poly_print(poly2); flint_printf("\n\n");
            abort();
        }

        nmod_poly_clear(poly1);
        nmod_poly_clear(poly2);
        nmod_poly_clear(poly3);
        nmod_poly_clear(xp);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}
