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

    Copyright (C) 2011 Sebastian Pancratz
    Copyright (C) 2009 William Hart
    Copyright (C) 2013 Mike Hansen

******************************************************************************/

#include "fq_poly.h"

int
main(void)
{
    int i, result;
    flint_rand_t state;

    flint_printf("div_basecase....");
    fflush(stdout);

    flint_randinit(state);

    /* Compare to divrem_basecase */
    for (i = 0; i < 500; i++)
    {
        fq_ctx_t ctx;
        fq_poly_t a, b, q, q2, r2;

        fq_ctx_randtest(ctx, state);
        
        fq_poly_init(a);
        fq_poly_init(b);
        fq_poly_init(q);
        fq_poly_init(q2);
        fq_poly_init(r2);

        fq_poly_randtest(a, state, n_randint(state, 100), ctx);
        fq_poly_randtest_not_zero(b, state, n_randint(state, 100) + 1, ctx);

        fq_poly_div_basecase(q, a, b, ctx);
        fq_poly_divrem_basecase(q2, r2, a, b, ctx);

        result = (fq_poly_equal(q, q2));
        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("ctx = "), fq_ctx_print(ctx), flint_printf("\n\n");
            flint_printf("a = "), fq_poly_print(a, ctx), flint_printf("\n\n");
            flint_printf("b = "), fq_poly_print(b, ctx), flint_printf("\n\n");
            flint_printf("q = "), fq_poly_print(q, ctx), flint_printf("\n\n");
            flint_printf("q2 = "), fq_poly_print(q2, ctx), flint_printf("\n\n");
            flint_printf("r2 = "), fq_poly_print(r2, ctx), flint_printf("\n\n");
            abort();
        }

        fq_poly_clear(a);
        fq_poly_clear(b);
        fq_poly_clear(q);
        fq_poly_clear(q2);
        fq_poly_clear(r2);
        fq_ctx_clear(ctx);
    }

    /* Alias a and q */
    for (i = 0; i < 500; i++)
    {
        fq_ctx_t ctx;
        fq_poly_t a, b, q;

        fq_ctx_randtest(ctx, state);

        fq_poly_init(a);
        fq_poly_init(b);
        fq_poly_init(q);
        fq_poly_randtest(a, state, n_randint(state, 100), ctx);
        fq_poly_randtest_not_zero(b, state, n_randint(state, 100) + 1, ctx);

        fq_poly_div_basecase(q, a, b, ctx);
        fq_poly_div_basecase(a, a, b, ctx);

        result = (fq_poly_equal(q, a));
        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("ctx = "), fq_ctx_print(ctx), flint_printf("\n\n");
            flint_printf("a = "), fq_poly_print(a, ctx), flint_printf("\n\n");
            flint_printf("b = "), fq_poly_print(b, ctx), flint_printf("\n\n");
            flint_printf("q = "), fq_poly_print(q, ctx), flint_printf("\n\n");
            abort();
        }

        fq_poly_clear(a);
        fq_poly_clear(b);
        fq_poly_clear(q);
        fq_ctx_clear(ctx);
    }

    /* Alias b and q */
    for (i = 0; i < 500; i++)
    {
        fq_ctx_t ctx;
        fq_poly_t a, b, q;

        fq_ctx_randtest(ctx, state);

        fq_poly_init(a);
        fq_poly_init(b);
        fq_poly_init(q);
        fq_poly_randtest(a, state, n_randint(state, 100), ctx);
        fq_poly_randtest_not_zero(b, state, n_randint(state, 100) + 1, ctx);

        fq_poly_div_basecase(q, a, b, ctx);
        fq_poly_div_basecase(b, a, b, ctx);

        result = (fq_poly_equal(q, b));
        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("ctx = "), fq_ctx_print(ctx), flint_printf("\n\n");
            flint_printf("a = "), fq_poly_print(a, ctx), flint_printf("\n\n");
            flint_printf("b = "), fq_poly_print(b, ctx), flint_printf("\n\n");
            flint_printf("q = "), fq_poly_print(q, ctx), flint_printf("\n\n");
            abort();
        }

        fq_poly_clear(a);
        fq_poly_clear(b);
        fq_poly_clear(q);
        fq_ctx_clear(ctx);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return 0;
}
