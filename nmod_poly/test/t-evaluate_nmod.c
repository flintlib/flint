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
    flint_rand_t state;
    flint_randinit(state);
    
    printf("evaluate_nmod....");
    fflush(stdout);

    /* Check evaluation at 1 gives sum of coeffs */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        nmod_poly_t a;
        mp_limb_t n = n_randtest_not_zero(state);
        mp_limb_t sum, eval;

        nmod_poly_init(a, n);
        nmod_poly_randtest(a, state, n_randint(state, 100));
        
        eval = nmod_poly_evaluate_nmod(a, 1);
        
        sum = 0;
        for (j = 0; j < a->length; j++)
           sum = n_addmod(sum, nmod_poly_get_coeff_ui(a, j), n);
        
        result = (sum == eval);
        if (!result)
        {
            printf("FAIL:\n");
            printf("a->length = %ld, n = %lu\n", a->length, a->mod.n);
            printf("sum = %lu, eval = %lu\n", sum, eval);
            nmod_poly_print(a), printf("\n\n");
            abort();
        }

        nmod_poly_clear(a);
    }

    /* Check a(c) + b(c) = (a + b)(c) */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        nmod_poly_t a, b;
        mp_limb_t n = n_randtest_not_zero(state);
        mp_limb_t eval1, eval2, c;

        nmod_poly_init(a, n);
        nmod_poly_init(b, n);
        nmod_poly_randtest(a, state, n_randint(state, 100));
        nmod_poly_randtest(b, state, n_randint(state, 100));
        
        c = n_randint(state, n);
        
        eval1 = nmod_poly_evaluate_nmod(a, c);
        eval1 = n_addmod(eval1, nmod_poly_evaluate_nmod(b, c), n);
        
        nmod_poly_add(a, a, b);
        eval2 = nmod_poly_evaluate_nmod(a, c);
        

        result = (eval1 == eval2);
        if (!result)
        {
            printf("FAIL:\n");
            printf("eval1 = %lu, eval2 = %lu\n", eval1, eval2);
            nmod_poly_print(a), printf("\n\n");
            nmod_poly_print(b), printf("\n\n");
            abort();
        }

        nmod_poly_clear(a);
        nmod_poly_clear(b);
    }

    flint_randclear(state);

    printf("PASS\n");
    return 0;
}
