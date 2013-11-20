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

    Copyright (C) 2013 Mike Hansen

******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fq_nmod.h"
#include "long_extras.h"

int
main(void)
{
    int i, result;
    flint_rand_t state;

    printf("mul_ui....");
    fflush(stdout);

    flint_randinit(state);
    /* Check aliasing of a, b */
    for (i = 0; i < 2000; i++)
    {
        fq_nmod_ctx_t ctx;
        ulong x;
        fq_nmod_t a, b;

        fq_nmod_ctx_randtest(ctx, state);
        
        fq_nmod_init(a, ctx);
        fq_nmod_init(b, ctx);

        fq_nmod_randtest(a, state, ctx);
        x = z_randtest(state);
        fq_nmod_mul_ui(b, a, x, ctx);
        fq_nmod_mul_ui(a, a, x, ctx);

        result = (fq_nmod_equal(a, b, ctx));
        if (!result)
        {
            printf("FAIL 1:\n\n");
            printf("a = "), fq_nmod_print_pretty(a, ctx), printf("\n");
            printf("b = "), fq_nmod_print_pretty(b, ctx), printf("\n");
	    printf("x = %lu\n",x);
            abort();
        }

        fq_nmod_clear(a, ctx);
        fq_nmod_clear(b, ctx);

        fq_nmod_ctx_clear(ctx);
    }

    /* compare with direct multiplication */
    for (i = 0; i < 2000; i++)
    {
        fq_nmod_ctx_t ctx;
        ulong x;
        fq_nmod_t a, c;
	nmod_poly_t b;

        fq_nmod_ctx_randtest(ctx, state);
        
        fq_nmod_init(a, ctx);
        fq_nmod_init(c, ctx);
	nmod_poly_init(b, ctx->mod.n);

        fq_nmod_randtest(a, state, ctx);
        x = n_randint(state, ctx->mod.n);
        fq_nmod_mul_ui(c, a, x, ctx);
        nmod_poly_scalar_mul_nmod(b,a,x);

        result = (fq_nmod_equal(c, b, ctx));
        if (!result)
        {
            printf("FAIL 2:\n\n");
            printf("a = "), fq_nmod_print_pretty(a, ctx), printf("\n");
            printf("b = "), fq_nmod_print_pretty(b, ctx), printf("\n");
            printf("c = "), fq_nmod_print_pretty(c, ctx), printf("\n");
	    printf("x = %lu\n",x);
            abort();
        }

        fq_nmod_clear(a, ctx);
        fq_nmod_clear(c, ctx);
        nmod_poly_clear(b);
        fq_nmod_ctx_clear(ctx);
    }

    flint_randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return 0;
}
