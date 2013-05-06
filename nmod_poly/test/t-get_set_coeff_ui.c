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
    ulong j;
    flint_rand_t state;
    flint_randinit(state);

    printf("get/set_coeff_ui....");
    fflush(stdout);

    /* Check aliasing of a and b */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        nmod_poly_t a;
        mp_limb_t n = n_randtest_not_zero(state);
        mp_limb_t c1 = n_randtest(state), c2;
        
        j = n_randint(state, 100);

        nmod_poly_init(a, n);
        nmod_poly_randtest(a, state, n_randint(state, 100));

        nmod_poly_set_coeff_ui(a, j, c1);
        c2 = nmod_poly_get_coeff_ui(a, j);
        
        result = (c2 == c1 % n);
        if (!result)
        {
            printf("FAIL:\n");
            printf("j = %lu, c1 = %lu, c2 = %lu, n = %lu\n", j, c1, c2, a->mod.n);
            nmod_poly_print(a), printf("\n\n");
            abort();
        }

        nmod_poly_clear(a);
    }

    flint_randclear(state);

    printf("PASS\n");
    return 0;
}
