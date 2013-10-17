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
    Copyright (C) 2012 Andres Goens

******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include "fq_poly.h"

int
main(void)
{
    int i, result;
    flint_rand_t state;
    flint_randinit(state);

    flint_printf("make_monic....");
    fflush(stdout);

    /* test aliasing */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {

        fmpz_t p;
        long d,len;
        fq_ctx_t ctx;

        fq_poly_t a,b;
        fmpz_init(p);
        fmpz_set_ui(p, n_randprime(state, 2 + n_randint(state, 3), 1));
        d = n_randint(state, 10) + 1;
        len = n_randint(state, 15) + 1;
        fq_ctx_init_conway(ctx, p, d, "a");
        fq_poly_init(a, ctx);
        fq_poly_init(b, ctx);

        fq_poly_randtest_not_zero(a, state, len, ctx);
        fq_poly_make_monic(b,a,ctx);
        fq_poly_make_monic(a,a,ctx);
        result = fq_poly_equal(a, b, ctx);
        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("a = "), fq_poly_print_pretty(a, "X", ctx), flint_printf("\n");
            flint_printf("b = "), fq_poly_print_pretty(a, "X", ctx), flint_printf("\n");
            
        }

        fq_poly_clear(a, ctx);
    }



    /* Check new leading coeff = 1 */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {

        fmpz_t p;
        long d,len;
        fq_ctx_t ctx;

        fq_poly_t a;

        fmpz_init(p);
        fmpz_set_ui(p, n_randprime(state, 2 + n_randint(state, 3), 1));
        d = n_randint(state, 10) + 1;
        len = n_randint(state, 15) + 1;
        fq_ctx_init_conway(ctx, p, d, "a");
        fq_poly_init(a, ctx);
        
        fq_poly_randtest_not_zero(a, state, len, ctx);
        fq_poly_make_monic(a,a,ctx);
        
        result = fq_is_one(fq_poly_lead(a, ctx), ctx);
        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("a = "), fq_poly_print_pretty(a, "X", ctx), flint_printf("\n");
            flint_printf("lead ="), fq_print_pretty(fq_poly_lead(a, ctx), ctx), flint_printf("\n");
        }

        fq_poly_clear(a, ctx);

    }
    flint_randclear(state);

    flint_printf("PASS\n");
    return 0;
}
