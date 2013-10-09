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
#include "ulong_extras.h"
#include "fq_poly.h"

int
main(void)
{
    int i, result;
    flint_rand_t state;

    printf("mulmod....");
    fflush(stdout);

    flint_randinit(state);

    /* Check aliasing of res and a */
    for (i = 0; i < 2000; i++)
    {
        fq_ctx_t ctx;
        fq_poly_t a, b, res, f;

        fq_ctx_randtest(ctx, state);

        fq_poly_init(a);
        fq_poly_init(b);
        fq_poly_init(f);
        fq_poly_init(res);

        fq_poly_randtest(a, state, n_randint(state, 50), ctx);
        fq_poly_randtest(b, state, n_randint(state, 50), ctx);
        fq_poly_randtest_not_zero(f, state, n_randint(state, 50) + 1, ctx);

        fq_poly_mulmod(res, a, b, f, ctx);
        fq_poly_mulmod(a, a, b, f, ctx);

        result = (fq_poly_equal(res, a));
        if (!result)
        {
            printf("FAIL:\n");
            printf("a:\n"); fq_poly_print(a, ctx), printf("\n\n");
            printf("b:\n"); fq_poly_print(b, ctx), printf("\n\n");
            printf("f:\n"); fq_poly_print(f, ctx), printf("\n\n");
            printf("res:\n"); fq_poly_print(res, ctx), printf("\n\n");
            abort();
        }

        fq_poly_clear(a);
        fq_poly_clear(b);
        fq_poly_clear(f);
        fq_poly_clear(res);
        
        fq_ctx_clear(ctx);
    }

    /* Check aliasing of res and b */
    for (i = 0; i < 2000; i++)
    {
        fq_ctx_t ctx;
        fq_poly_t a, b, f, res;

        fq_ctx_randtest(ctx, state);

        fq_poly_init(a);
        fq_poly_init(b);
        fq_poly_init(f);
        fq_poly_init(res);

        fq_poly_randtest(a, state, n_randint(state, 50), ctx);
        fq_poly_randtest(b, state, n_randint(state, 50), ctx);
        fq_poly_randtest_not_zero(f, state, n_randint(state, 50) + 1, ctx);

        fq_poly_mulmod(res, a, b, f, ctx);
        fq_poly_mulmod(b, a, b, f, ctx);

        result = (fq_poly_equal(res, b));
        if (!result)
        {
            printf("FAIL:\n");
            printf("a:\n"); fq_poly_print(a, ctx), printf("\n\n");
            printf("b:\n"); fq_poly_print(b, ctx), printf("\n\n");
            printf("f:\n"); fq_poly_print(f, ctx), printf("\n\n");
            printf("res:\n"); fq_poly_print(res, ctx), printf("\n\n");
            abort();
        }

        fq_poly_clear(a);
        fq_poly_clear(b);
        fq_poly_clear(f);
        fq_poly_clear(res);

        fq_ctx_clear(ctx);
    }

    /* Check aliasing of res and f */
    for (i = 0; i < 2000; i++)
    {
        fq_ctx_t ctx;
        fq_poly_t a, b, f, res;

        fq_ctx_randtest(ctx, state);
        
        fq_poly_init(a);
        fq_poly_init(b);
        fq_poly_init(f);
        fq_poly_init(res);

        fq_poly_randtest(a, state, n_randint(state, 50), ctx);
        fq_poly_randtest(b, state, n_randint(state, 50), ctx);
        fq_poly_randtest_not_zero(f, state, n_randint(state, 50) + 1, ctx);

        fq_poly_mulmod(res, a, b, f, ctx);
        fq_poly_mulmod(f, a, b, f, ctx);

        result = (fq_poly_equal(res, f));
        if (!result)
        {
            printf("FAIL:\n");
            printf("a:\n"); fq_poly_print(a, ctx), printf("\n\n");
            printf("b:\n"); fq_poly_print(b, ctx), printf("\n\n");
            printf("f:\n"); fq_poly_print(f, ctx), printf("\n\n");
            printf("res:\n"); fq_poly_print(res, ctx), printf("\n\n");
            abort();
        }

        fq_poly_clear(a);
        fq_poly_clear(b);
        fq_poly_clear(f);
        fq_poly_clear(res);
        
        fq_ctx_clear(ctx);
    }

    /* No aliasing */
    for (i = 0; i < 1000; i++)
    {
        fq_ctx_t ctx;
        fq_poly_t a, b, res1, res2, t, f;

        fq_ctx_randtest(ctx, state);

        fq_poly_init(a);
        fq_poly_init(b);
        fq_poly_init(f);

        fq_poly_randtest(a, state, n_randint(state, 50), ctx);
        fq_poly_randtest(b, state, n_randint(state, 50), ctx);
        fq_poly_randtest_not_zero(f, state, n_randint(state, 50) + 1, ctx);

        fq_poly_init(res1);
        fq_poly_init(res2);
        fq_poly_init(t);
        fq_poly_mulmod(res1, a, b, f, ctx);
        fq_poly_mul(res2, a, b, ctx);
        fq_poly_divrem(t, res2, res2, f, ctx);

        result = (fq_poly_equal(res1, res2));
        if (!result)
        {
            printf("FAIL:\n");
            printf("a:\n"); fq_poly_print(a, ctx), printf("\n\n");
            printf("b:\n"); fq_poly_print(b, ctx), printf("\n\n");
            printf("f:\n"); fq_poly_print(f, ctx), printf("\n\n");
            printf("res1:\n"); fq_poly_print(res1, ctx), printf("\n\n");
            printf("res2:\n"); fq_poly_print(res2, ctx), printf("\n\n");
            abort();
        }

        fq_poly_clear(a);
        fq_poly_clear(b);
        fq_poly_clear(f);
        fq_poly_clear(res1);
        fq_poly_clear(res2);
        fq_poly_clear(t);
        
        fq_ctx_clear(ctx);
    }

    flint_randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return 0;
}
