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
    Copyright (C) 2013 Martin Lee

******************************************************************************/

#include "fq_poly.h"

int
main(void)
{
    int i, result;
    flint_rand_t state;
    flint_randinit(state);

    flint_printf("div_newton_preinv....");
    fflush(stdout);

    /* Check result of divrem */
    for (i = 0; i < 50 * flint_test_multiplier(); i++)
    {
        fq_poly_t a, b, binv, q, r, test;

        fq_ctx_t ctx;
        
        fq_ctx_randtest(ctx, state);

        fq_poly_init(a);
        fq_poly_init(b);
        fq_poly_init(binv);
        fq_poly_init(q);
        fq_poly_init(r);
        fq_poly_init(test);

        do
            fq_poly_randtest(b, state, n_randint(state, 200), ctx);
        while (b->length <= 2);
        fq_poly_randtest(a, state, n_randint(state, 200), ctx);
        if (a->length > 2*(b->length)-3)
            fq_poly_truncate (a, 2*(b->length)-3, ctx);

        fq_poly_reverse (binv, b, b->length, ctx);
        fq_poly_inv_series_newton (binv, binv, b->length, ctx);
        fq_poly_div_newton_preinv(q, a, b, binv, ctx);
        fq_poly_divrem (test, r, a, b, ctx);

        result = (fq_poly_equal(q, test));
        if (!result)
        {
            flint_printf("FAIL:\n");
            fq_poly_print(test, ctx), flint_printf("\n\n");
            fq_poly_print(q, ctx), flint_printf("\n\n");
            fq_poly_print(a, ctx), flint_printf("\n\n");
            fq_poly_print(b, ctx), flint_printf("\n\n");
            abort();
        }
        
        fq_poly_clear(a);
        fq_poly_clear(b);
        fq_poly_clear(binv);
        fq_poly_clear(q);
        fq_poly_clear(r);
        fq_poly_clear(test);
        fq_ctx_clear(ctx);
    }

    /* Check aliasing of a and q */
    for (i = 0; i < 50 * flint_test_multiplier(); i++)
    {
        fq_poly_t a, b, binv, q;

        fq_ctx_t ctx;

        fq_ctx_randtest(ctx, state);
        
        fq_poly_init(a);
        fq_poly_init(b);
        fq_poly_init(binv);
        fq_poly_init(q);

        do
            fq_poly_randtest(b, state, n_randint(state, 200), ctx);
        while (b->length <= 2);
        fq_poly_randtest(a, state, n_randint(state, 200), ctx);
        if (a->length > 2*(b->length)-3)
            fq_poly_truncate (a, 2*(b->length)-3, ctx);

        fq_poly_reverse (binv, b, b->length, ctx);
        fq_poly_inv_series_newton (binv, binv, b->length, ctx);

        fq_poly_div_newton_preinv(q, a, b, binv, ctx);
        fq_poly_div_newton_preinv(a, a, b, binv, ctx);

        result = (fq_poly_equal(a, q));
        if (!result)
        {
            flint_printf("FAIL:\n");
            fq_poly_print(a, ctx), flint_printf("\n\n");
            fq_poly_print(b, ctx), flint_printf("\n\n");
            fq_poly_print(q, ctx), flint_printf("\n\n");
            abort();
        }

        fq_poly_clear(a);
        fq_poly_clear(b);
        fq_poly_clear(binv);
        fq_poly_clear(q);
        fq_ctx_clear(ctx);
    }

    /* Check aliasing of b and q */
    for (i = 0; i < 50 * flint_test_multiplier(); i++)
    {
        fq_poly_t a, b, binv, q;

        fq_ctx_t ctx;

        fq_ctx_randtest(ctx, state);

        fq_poly_init(a);
        fq_poly_init(b);
        fq_poly_init(binv);
        fq_poly_init(q);
        
        do
            fq_poly_randtest(b, state, n_randint(state, 200), ctx);
        while (b->length <= 2);
        fq_poly_randtest(a, state, n_randint(state, 200), ctx);
        if (a->length > 2*(b->length)-3)
            fq_poly_truncate (a, 2*(b->length)-3, ctx);

        fq_poly_reverse (binv, b, b->length, ctx);
        fq_poly_inv_series_newton (binv, binv, b->length, ctx);

        fq_poly_div_newton_preinv(q, a, b, binv, ctx);
        fq_poly_div_newton_preinv(b, a, b, binv, ctx);

        result = (fq_poly_equal(b, q));
        if (!result)
        {
            flint_printf("FAIL:\n");
            fq_poly_print(a, ctx), flint_printf("\n\n");
            fq_poly_print(b, ctx), flint_printf("\n\n");
            fq_poly_print(q, ctx), flint_printf("\n\n");
            abort();
        }

        fq_poly_clear(a);
        fq_poly_clear(b);
        fq_poly_clear(binv);
        fq_poly_clear(q);
        fq_ctx_clear(ctx);
    }

    /* Check aliasing of binv and q */
    for (i = 0; i < 50 * flint_test_multiplier(); i++)
    {
        fq_poly_t a, b, binv, q;

        fq_ctx_t ctx;

        fq_ctx_randtest(ctx, state);

        fq_poly_init(a);
        fq_poly_init(b);
        fq_poly_init(binv);
        fq_poly_init(q);

        
        do
            fq_poly_randtest(b, state, n_randint(state, 200), ctx);
        while (b->length <= 2);
        fq_poly_randtest(a, state, n_randint(state, 200), ctx);
        if (a->length > 2*(b->length)-3)
            fq_poly_truncate (a, 2*(b->length)-3, ctx);

        fq_poly_reverse (binv, b, b->length, ctx);
        fq_poly_inv_series_newton (binv, binv, b->length, ctx);

        fq_poly_div_newton_preinv(q, a, b, binv, ctx);
        fq_poly_div_newton_preinv(binv, a, b, binv, ctx);

        result = (fq_poly_equal(binv, q));
        if (!result)
        {
            flint_printf("FAIL:\n");
            fq_poly_print(a, ctx), flint_printf("\n\n");
            fq_poly_print(b, ctx), flint_printf("\n\n");
            fq_poly_print(binv, ctx), flint_printf("\n\n");
            fq_poly_print(q, ctx), flint_printf("\n\n");
            abort();
        }

        fq_poly_clear(a);
        fq_poly_clear(b);
        fq_poly_clear(binv);
        fq_poly_clear(q);
        fq_ctx_clear(ctx);
    }

    flint_randclear(state);

    flint_printf("PASS\n");
    return 0;
}
