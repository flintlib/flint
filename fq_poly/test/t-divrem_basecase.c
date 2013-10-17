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

    flint_printf("divrem_basecase... ");
    fflush(stdout);

    flint_randinit(state);

    /* Check q*b + r == a */
    for (i = 0; i < 2000; i++)
    {
        fq_ctx_t ctx;

        fq_poly_t a, b, c, q, r;

        fq_ctx_randtest(ctx, state);
        fq_poly_init(a, ctx);
        fq_poly_init(b, ctx);
        fq_poly_init(c, ctx);
        fq_poly_init(q, ctx);
        fq_poly_init(r, ctx);

        fq_poly_randtest(         a, state, n_randint(state, 50),     ctx);
        fq_poly_randtest_not_zero(b, state, n_randint(state, 50) + 1, ctx);

        fq_poly_divrem_basecase(q, r, a, b, ctx);
        fq_poly_mul(c, q, b, ctx);
        fq_poly_add(c, c, r, ctx);

        result = (fq_poly_equal(a, c, ctx));
        if (!result)
        {
            flint_printf("FAIL:\n\n");
            flint_printf("a = "), fq_poly_print_pretty(a, "X", ctx), flint_printf("\n");
            flint_printf("b = "), fq_poly_print_pretty(b, "X", ctx), flint_printf("\n");
            flint_printf("c = "), fq_poly_print_pretty(c, "X", ctx), flint_printf("\n");
            flint_printf("q = "), fq_poly_print_pretty(q, "X", ctx), flint_printf("\n");
            flint_printf("r = "), fq_poly_print_pretty(r, "X", ctx), flint_printf("\n");
            abort();
        }

        fq_poly_clear(a, ctx);
        fq_poly_clear(b, ctx);
        fq_poly_clear(c, ctx);
        fq_poly_clear(q, ctx);
        fq_poly_clear(r, ctx);

        fq_ctx_clear(ctx);
    }

    /* Check aliasing: a and r */
    for (i = 0; i < 2000; i++)
    {
        fq_ctx_t ctx;

        fq_poly_t a, b, q, r;

        fq_ctx_randtest(ctx, state);
        fq_poly_init(a, ctx);
        fq_poly_init(b, ctx);
        fq_poly_init(q, ctx);
        fq_poly_init(r, ctx);

        fq_poly_randtest(         a, state, n_randint(state, 50),     ctx);
        fq_poly_randtest_not_zero(b, state, n_randint(state, 50) + 1, ctx);

        fq_poly_divrem_basecase(q, r, a, b, ctx);
        fq_poly_divrem_basecase(q, a, a, b, ctx);

        result = (fq_poly_equal(a, r, ctx));
        if (!result)
        {
            flint_printf("FAIL:\n\n");
            flint_printf("a = "), fq_poly_print_pretty(a, "X", ctx), flint_printf("\n");
            flint_printf("b = "), fq_poly_print_pretty(b, "X", ctx), flint_printf("\n");
            flint_printf("q = "), fq_poly_print_pretty(q, "X", ctx), flint_printf("\n");
            flint_printf("r = "), fq_poly_print_pretty(r, "X", ctx), flint_printf("\n");
            abort();
        }

        fq_poly_clear(a, ctx);
        fq_poly_clear(b, ctx);
        fq_poly_clear(q, ctx);
        fq_poly_clear(r, ctx);

        fq_ctx_clear(ctx);
    }

    /* Check aliasing: b and r */
    for (i = 0; i < 2000; i++)
    {
        fq_ctx_t ctx;

        fq_poly_t a, b, q, r;

        fq_ctx_randtest(ctx, state);
        fq_poly_init(a, ctx);
        fq_poly_init(b, ctx);
        fq_poly_init(q, ctx);
        fq_poly_init(r, ctx);

        fq_poly_randtest(         a, state, n_randint(state, 50),     ctx);
        fq_poly_randtest_not_zero(b, state, n_randint(state, 50) + 1, ctx);

        fq_poly_divrem_basecase(q, r, a, b, ctx);
        fq_poly_divrem_basecase(q, b, a, b, ctx);

        result = (fq_poly_equal(b, r, ctx));
        if (!result)
        {
            flint_printf("FAIL:\n\n");
            flint_printf("a = "), fq_poly_print_pretty(a, "X", ctx), flint_printf("\n");
            flint_printf("b = "), fq_poly_print_pretty(b, "X", ctx), flint_printf("\n");
            flint_printf("q = "), fq_poly_print_pretty(q, "X", ctx), flint_printf("\n");
            flint_printf("r = "), fq_poly_print_pretty(r, "X", ctx), flint_printf("\n");
            abort();
        }

        fq_poly_clear(a, ctx);
        fq_poly_clear(b, ctx);
        fq_poly_clear(q, ctx);
        fq_poly_clear(r, ctx);

        fq_ctx_clear(ctx);
    }

    /* Check aliasing: a and q */
    for (i = 0; i < 2000; i++)
    {
        fq_ctx_t ctx;

        fq_poly_t a, b, q, r;

        fq_ctx_randtest(ctx, state);
        fq_poly_init(a, ctx);
        fq_poly_init(b, ctx);
        fq_poly_init(q, ctx);
        fq_poly_init(r, ctx);

        fq_poly_randtest(         a, state, n_randint(state, 50),     ctx);
        fq_poly_randtest_not_zero(b, state, n_randint(state, 50) + 1, ctx);

        fq_poly_divrem_basecase(q, r, a, b, ctx);
        fq_poly_divrem_basecase(a, r, a, b, ctx);

        result = (fq_poly_equal(a, q, ctx));
        if (!result)
        {
            flint_printf("FAIL:\n\n");
            flint_printf("a = "), fq_poly_print_pretty(a, "X", ctx), flint_printf("\n");
            flint_printf("b = "), fq_poly_print_pretty(b, "X", ctx), flint_printf("\n");
            flint_printf("q = "), fq_poly_print_pretty(q, "X", ctx), flint_printf("\n");
            flint_printf("r = "), fq_poly_print_pretty(r, "X", ctx), flint_printf("\n");
            abort();
        }

        fq_poly_clear(a, ctx);
        fq_poly_clear(b, ctx);
        fq_poly_clear(q, ctx);
        fq_poly_clear(r, ctx);

        fq_ctx_clear(ctx);
    }

    /* Check aliasing: b and q */
    for (i = 0; i < 2000; i++)
    {
        fq_ctx_t ctx;

        fq_poly_t a, b, q, r;

        fq_ctx_randtest(ctx, state);
        fq_poly_init(a, ctx);
        fq_poly_init(b, ctx);
        fq_poly_init(q, ctx);
        fq_poly_init(r, ctx);

        fq_poly_randtest(         a, state, n_randint(state, 50),     ctx);
        fq_poly_randtest_not_zero(b, state, n_randint(state, 50) + 1, ctx);

        fq_poly_divrem_basecase(q, r, a, b, ctx);
        fq_poly_divrem_basecase(b, r, a, b, ctx);

        result = (fq_poly_equal(b, q, ctx));
        if (!result)
        {
            flint_printf("FAIL:\n\n");
            flint_printf("a = "), fq_poly_print_pretty(a, "X", ctx), flint_printf("\n");
            flint_printf("b = "), fq_poly_print_pretty(b, "X", ctx), flint_printf("\n");
            flint_printf("q = "), fq_poly_print_pretty(q, "X", ctx), flint_printf("\n");
            flint_printf("r = "), fq_poly_print_pretty(r, "X", ctx), flint_printf("\n");
            abort();
        }

        fq_poly_clear(a, ctx);
        fq_poly_clear(b, ctx);
        fq_poly_clear(q, ctx);
        fq_poly_clear(r, ctx);

        fq_ctx_clear(ctx);
    }

    flint_randclear(state);
    _fmpz_cleanup();
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}

