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

#include "fq_poly.h"

int
main(void)
{
    int i, result;
    flint_rand_t state;

    printf("divrem_divconquer....");
    fflush(stdout);

    flint_randinit(state);

    /* Check q*b + r = a, no aliasing */
    for (i = 0; i < 5000; i++)
    {
        fq_ctx_t ctx;
        fq_poly_t a, b, q, r, t;

        fq_ctx_randtest(ctx, state);

        fq_poly_init(a, ctx);
        fq_poly_init(b, ctx);
        fq_poly_init(q, ctx);
        fq_poly_init(r, ctx);
        fq_poly_init(t, ctx);
        fq_poly_randtest(a, state, n_randint(state, 100), ctx);
        fq_poly_randtest_not_zero(b, state, n_randint(state, 100) + 1, ctx);

        fq_poly_divrem_divconquer(q, r, a, b, ctx);

        fq_poly_mul(t, q, b, ctx);
        fq_poly_add(t, t, r, ctx);

        result = (fq_poly_equal(a, t, ctx));
        if (!result)
        {
            printf("FAIL #1:\n");
            printf("a = "), fq_poly_print(a, ctx), printf("\n\n");
            printf("b = "), fq_poly_print(b, ctx), printf("\n\n");
            printf("q = "), fq_poly_print(q, ctx), printf("\n\n");
            printf("r = "), fq_poly_print(r, ctx), printf("\n\n");
            printf("t = "), fq_poly_print(t, ctx), printf("\n\n");
            abort();
        }

        fq_poly_clear(a, ctx);
        fq_poly_clear(b, ctx);
        fq_poly_clear(q, ctx);
        fq_poly_clear(r, ctx);
        fq_poly_clear(t, ctx);
        fq_ctx_clear(ctx);
    }

    /* Alias a and q, b and r */
    for (i = 0; i < 500; i++)
    {
        fq_ctx_t ctx;
        fq_poly_t a, b, q, r;

        fq_ctx_randtest(ctx, state);

        fq_poly_init(a, ctx);
        fq_poly_init(b, ctx);
        fq_poly_init(q, ctx);
        fq_poly_init(r, ctx);
        fq_poly_randtest(a, state, n_randint(state, 100), ctx);
        fq_poly_randtest_not_zero(b, state, n_randint(state, 100) + 1, ctx);

        fq_poly_divrem_divconquer(q, r, a, b, ctx);
        fq_poly_divrem_divconquer(a, b, a, b, ctx);

        result = (fq_poly_equal(q, a, ctx) && fq_poly_equal(r, b, ctx));
        if (!result)
        {
            printf("FAIL #2:\n");
            printf("a = "), fq_poly_print(a, ctx), printf("\n\n");
            printf("b = "), fq_poly_print(b, ctx), printf("\n\n");
            printf("q = "), fq_poly_print(q, ctx), printf("\n\n");
            printf("r = "), fq_poly_print(r, ctx), printf("\n\n");
            abort();
        }

        fq_poly_clear(a, ctx);
        fq_poly_clear(b, ctx);
        fq_poly_clear(q, ctx);
        fq_poly_clear(r, ctx);
        fq_ctx_clear(ctx);
    }

    /* Alias b and q, a and r */
    for (i = 0; i < 500; i++)
    {
        fq_ctx_t ctx;
        fq_poly_t a, b, q, r;

        fq_ctx_randtest(ctx, state);

        fq_poly_init(a, ctx);
        fq_poly_init(b, ctx);
        fq_poly_init(q, ctx);
        fq_poly_init(r, ctx);
        fq_poly_randtest(a, state, n_randint(state, 100), ctx);
        fq_poly_randtest_not_zero(b, state, n_randint(state, 100) + 1, ctx);

        fq_poly_divrem_divconquer(q, r, a, b, ctx);
        fq_poly_divrem_divconquer(b, a, a, b, ctx);

        result = (fq_poly_equal(q, b, ctx) && fq_poly_equal(r, a, ctx));
        if (!result)
        {
            printf("FAIL #3:\n");
            printf("a = "), fq_poly_print(a, ctx), printf("\n\n");
            printf("b = "), fq_poly_print(b, ctx), printf("\n\n");
            printf("q = "), fq_poly_print(q, ctx), printf("\n\n");
            printf("r = "), fq_poly_print(r, ctx), printf("\n\n");
            abort();
        }

        fq_poly_clear(a, ctx);
        fq_poly_clear(b, ctx);
        fq_poly_clear(q, ctx);
        fq_poly_clear(r, ctx);
        fq_ctx_clear(ctx);
    }

    flint_randclear(state);
    flint_cleanup();
    printf("PASS\n");
    return 0;
}
