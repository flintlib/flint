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

    printf("make_monic....");
    fflush(stdout);

    /* Check new leading coeff = gcd old leading coeff and modulus */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        nmod_poly_t a, b;
        mp_limb_t n = n_randtest_not_zero(state);
        mp_limb_t l;

        nmod_poly_init(a, n);
        nmod_poly_init(b, n);
        
        if (n == 1) continue;
        do { nmod_poly_randtest(a, state, n_randint(state, 100) + 1); } while (a->length == 0);
        
        nmod_poly_make_monic(b, a);
        l = n_gcd(a->mod.n, a->coeffs[a->length - 1]);
        
        result = (l == b->coeffs[b->length - 1]);
        if (!result)
        {
            printf("FAIL:\n");
            printf("l = %lu, a->lead = %ld, n = %lu\n", 
                l, a->coeffs[a->length - 1], a->mod.n);
            nmod_poly_print(a), printf("\n\n");
            nmod_poly_print(b), printf("\n\n");
            abort();
        }

        nmod_poly_clear(a);
        nmod_poly_clear(b);
    }

    /* test aliasing */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        nmod_poly_t a;
        mp_limb_t n = n_randtest_not_zero(state);
        mp_limb_t l;

        nmod_poly_init(a, n);
        
        if (n == 1) continue;
        do { nmod_poly_randtest(a, state, n_randint(state, 100) + 1); } while (a->length == 0);
        
        l = n_gcd(a->mod.n, a->coeffs[a->length - 1]);
        nmod_poly_make_monic(a, a);
        
        result = (l == a->coeffs[a->length - 1]);
        if (!result)
        {
            printf("FAIL:\n");
            printf("l = %lu, a->lead = %ld, n = %lu\n", 
                l, a->coeffs[a->length - 1], a->mod.n);
            nmod_poly_print(a), printf("\n\n");
            abort();
        }

        nmod_poly_clear(a);
    }

    flint_randclear(state);

    printf("PASS\n");
    return 0;
}
