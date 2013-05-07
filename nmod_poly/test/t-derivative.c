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
#include "fmpz.h"

int
main(void)
{
    int i, j, result = 1;
    fmpz_t t;
    flint_rand_t state;
    flint_randinit(state);
    
    printf("derivative....");
    fflush(stdout);

    fmpz_init(t);

    /* Check derivative by hand */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        nmod_poly_t a, b;
        mp_limb_t n = n_randtest_not_zero(state);
        
        nmod_poly_init(a, n);
        nmod_poly_init(b, n);
        nmod_poly_randtest(a, state, n_randint(state, 100));
        
        nmod_poly_derivative(b, a);
        
        if (a->length <= 1)
            result = (b->length == 0);
        else
        {
            for (j = 1; j < a->length; j++)
            {
                fmpz_set_ui(t, nmod_poly_get_coeff_ui(a, j));
                fmpz_mul_ui(t, t, j);
                fmpz_mod_ui(t, t, n);
                result &= (fmpz_get_ui(t) == nmod_poly_get_coeff_ui(b, j - 1));
            }
        }
        
        if (!result)
        {
            printf("FAIL:\n");
            printf("a->length = %ld, n = %lu\n", a->length, a->mod.n);
            nmod_poly_print(a), printf("\n\n");
            nmod_poly_print(b), printf("\n\n");
            abort();
        }

        nmod_poly_clear(a);
        nmod_poly_clear(b);
    }

    fmpz_clear(t);

    /* Check aliasing */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        nmod_poly_t a, b;
        mp_limb_t n = n_randtest_not_zero(state);
        
        nmod_poly_init(a, n);
        nmod_poly_init(b, n);
        nmod_poly_randtest(a, state, n_randint(state, 100));
        
        nmod_poly_derivative(b, a);
        nmod_poly_derivative(a, a);
        
        result = nmod_poly_equal(a, b);
        if (!result)
        {
            printf("FAIL:\n");
            printf("a->length = %ld, n = %lu\n", a->length, a->mod.n);
            nmod_poly_print(a), printf("\n\n");
            nmod_poly_print(b), printf("\n\n");
            abort();
        }

        nmod_poly_clear(a);
        nmod_poly_clear(b);
    }

    flint_randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return 0;
}
