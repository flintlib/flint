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

    Copyright (C) 2012 Sebastian Pancratz 
    Copyright (C) 2012 Andres Goens

******************************************************************************/

#include <stdio.h>
#include <stdlib.h>

#include "fq_poly.h"

#include "ulong_extras.h"
#include "long_extras.h"

int
main(void)
{
    int i, result;
    flint_rand_t state;

    flint_printf("hamming_weight... ");
    fflush(stdout);

    flint_randinit(state);

    /* Check consistency */
        for (i = 0; i < 2000; i++)
    {
        long len;
        fq_ctx_t ctx;
	long w1,w2;
        fq_poly_t a, b;

        len = n_randint(state, 15) + 1;
        fq_ctx_randtest(ctx, state);
        fq_poly_init(a, ctx);
        fq_poly_init(b, ctx);

        fq_poly_randtest(a, state, len, ctx);
	fq_poly_set(b, a, ctx);

	w1 = fq_poly_hamming_weight(a, ctx);
	w2 = fq_poly_hamming_weight(b, ctx);

        result = (w1 == w2);
        if (!result)
        {
            flint_printf("FAIL:\n\n");
            flint_printf("a = "), fq_poly_print_pretty(a, "X", ctx), flint_printf("\n");
            flint_printf("b = "), fq_poly_print_pretty(b, "X", ctx), flint_printf("\n");
	    flint_printf("w1 = %ld \n w2 = %ld \n",w1,w2);
            abort();
        }

        fq_poly_clear(a, ctx);
        fq_poly_clear(b, ctx);


        fq_ctx_clear(ctx);
	}

    /* Check that wt(a+b) \leq wt(a) + wt(b) */
        for (i = 0; i < 2000; i++)
    {
        long len;
        fq_ctx_t ctx;
	long w1,w2,wsum;
        fq_poly_t a, b, c;

        len = n_randint(state, 15) + 1;
        fq_ctx_randtest(ctx, state);
        fq_poly_init(a, ctx);
        fq_poly_init(b, ctx);
        fq_poly_init(c, ctx);

        fq_poly_randtest(a, state, len, ctx);
        fq_poly_randtest(b, state, len, ctx);
	fq_poly_add(c,a,b,ctx);
	
	w1 = fq_poly_hamming_weight(a, ctx);
	w2 = fq_poly_hamming_weight(b, ctx);
	wsum = fq_poly_hamming_weight(c, ctx);
	
        result = (wsum <= w1+w2);
        if (!result)
        {
            flint_printf("FAIL:\n\n");
            flint_printf("a = "), fq_poly_print_pretty(a, "X", ctx), flint_printf("\n");
            flint_printf("b = "), fq_poly_print_pretty(b, "X", ctx), flint_printf("\n");
	    flint_printf("w1 = %ld \n w2 = %ld \n wsum = %ld",w1,w2,wsum);
            abort();
        }

        fq_poly_clear(a, ctx);
        fq_poly_clear(b, ctx);
        fq_poly_clear(c, ctx);


        fq_ctx_clear(ctx);
	} 

    flint_randclear(state);
    _fmpz_cleanup();
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}

