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

#include "fq.h"
#include "ulong_extras.h"
#include "long_extras.h"

int
main(void)
{
    int i, result;
    flint_rand_t state;

    flint_printf("mul... ");
    fflush(stdout);

    flint_randinit(state);

    /* Check aliasing: a = a * b */
    for (i = 0; i < 2000; i++)
    {
        fq_ctx_t ctx;
        fq_t a, b, c;

        fq_ctx_randtest(ctx, state);
        
        fq_init(a, ctx);
        fq_init(b, ctx);
        fq_init(c, ctx);

        fq_randtest(a, state, ctx);
        fq_randtest(b, state, ctx);

        fq_mul(c, a, b, ctx);
        fq_mul(a, a, b, ctx);

        result = (fq_equal(a, c, ctx));
        if (!result)
        {
            flint_printf("FAIL:\n\n");
            flint_printf("a = "), fq_print_pretty(a, ctx), flint_printf("\n");
            flint_printf("b = "), fq_print_pretty(b, ctx), flint_printf("\n");
            flint_printf("c = "), fq_print_pretty(c, ctx), flint_printf("\n");
            abort();
        }

        fq_clear(a, ctx);
        fq_clear(b, ctx);
        fq_clear(c, ctx);

        fq_ctx_clear(ctx);
    }

    /* Check aliasing: b = a * b */
    for (i = 0; i < 2000; i++)
    {
        fq_ctx_t ctx;
        fq_t a, b, c;

        fq_ctx_randtest(ctx, state);
        
        fq_init(a, ctx);
        fq_init(b, ctx);
        fq_init(c, ctx);

        fq_randtest(a, state, ctx);
        fq_randtest(b, state, ctx);

        fq_mul(c, a, b, ctx);
        fq_mul(b, a, b, ctx);

        result = (fq_equal(b, c, ctx));
        if (!result)
        {
            flint_printf("FAIL:\n\n");
            flint_printf("a = "), fq_print_pretty(a, ctx), flint_printf("\n");
            flint_printf("b = "), fq_print_pretty(b, ctx), flint_printf("\n");
            flint_printf("c = "), fq_print_pretty(c, ctx), flint_printf("\n");
            abort();
        }

        fq_clear(a, ctx);
        fq_clear(b, ctx);
        fq_clear(c, ctx);

        fq_ctx_clear(ctx);
    }

    /* Check aliasing: a = a * a */
    for (i = 0; i < 2000; i++)
    {
        fq_ctx_t ctx;

        fq_t a, c;

        fq_ctx_randtest(ctx, state);
        
        fq_init(a, ctx);
        fq_init(c, ctx);

        fq_randtest(a, state, ctx);

        fq_mul(c, a, a, ctx);
        fq_mul(a, a, a, ctx);

        result = (fq_equal(a, c, ctx));
        if (!result)
        {
            flint_printf("FAIL:\n\n");
            flint_printf("a = "), fq_print_pretty(a, ctx), flint_printf("\n");
            flint_printf("c = "), fq_print_pretty(c, ctx), flint_printf("\n");
            abort();
        }

        fq_clear(a, ctx);
        fq_clear(c, ctx);

        fq_ctx_clear(ctx);
    }

    /* Check that a * b == b * a */
    for (i = 0; i < 2000; i++)
    {
        fq_ctx_t ctx;
        fq_t a, b, c1, c2;

        fq_ctx_randtest(ctx, state);

        fq_init(a, ctx);
        fq_init(b, ctx);
        fq_init(c1, ctx);
        fq_init(c2, ctx);

        fq_randtest(a, state, ctx);
        fq_randtest(b, state, ctx);

        fq_mul(c1, a, b, ctx);
        fq_mul(c2, b, a, ctx);

        result = (fq_equal(c1, c2, ctx));
        if (!result)
        {
            flint_printf("FAIL:\n\n");
            flint_printf("a  = "), fq_print_pretty(a, ctx), flint_printf("\n");
            flint_printf("b  = "), fq_print_pretty(b, ctx), flint_printf("\n");
            flint_printf("c1 = "), fq_print_pretty(c1, ctx), flint_printf("\n");
            flint_printf("c2 = "), fq_print_pretty(c2, ctx), flint_printf("\n");
            abort();
        }

        fq_clear(a, ctx);
        fq_clear(b, ctx);
        fq_clear(c1, ctx);
        fq_clear(c2, ctx);

        fq_ctx_clear(ctx);
    }

    /* Check that (a * b) * c == a * (b * c) */
    for (i = 0; i < 2000; i++)
    {
        fq_ctx_t ctx;

        fq_t a, b, c, lhs, rhs;

        fq_ctx_randtest(ctx, state);
        
        fq_init(a, ctx);
        fq_init(b, ctx);
        fq_init(c, ctx);
        fq_init(lhs, ctx);
        fq_init(rhs, ctx);

        fq_randtest(a, state, ctx);
        fq_randtest(b, state, ctx);
        fq_randtest(c, state, ctx);

        fq_mul(lhs, a, b, ctx);
        fq_mul(lhs, lhs, c, ctx);
        fq_mul(rhs, b, c, ctx);
        fq_mul(rhs, a, rhs, ctx);

        result = (fq_equal(lhs, rhs, ctx));
        if (!result)
        {
            flint_printf("FAIL:\n\n");
            flint_printf("a   = "), fq_print_pretty(a, ctx), flint_printf("\n");
            flint_printf("b   = "), fq_print_pretty(b, ctx), flint_printf("\n");
            flint_printf("c   = "), fq_print_pretty(c, ctx), flint_printf("\n");
            flint_printf("lhs = "), fq_print_pretty(lhs, ctx), flint_printf("\n");
            flint_printf("rhs = "), fq_print_pretty(rhs, ctx), flint_printf("\n");
            abort();
        }

        fq_clear(a, ctx);
        fq_clear(b, ctx);
        fq_clear(c, ctx);
        fq_clear(lhs, ctx);
        fq_clear(rhs, ctx);

        fq_ctx_clear(ctx);
    }

    flint_randclear(state);
    _fmpz_cleanup();
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}

