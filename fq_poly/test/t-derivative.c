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

    flint_printf("derivative... ");
    fflush(stdout);

    flint_randinit(state);

    /* Check aliasing */
    for (i = 0; i < 2000; i++)
    {
        long len;
        fq_ctx_t ctx;

        fq_poly_t a, b, c;

        len = n_randint(state, 15) + 1;
        fq_ctx_randtest(ctx, state);
        fq_poly_init(a, ctx);
        fq_poly_init(b, ctx);
        fq_poly_init(c, ctx);

        fq_poly_randtest(a, state, len, ctx);
        fq_poly_set(b, a, ctx);
        fq_poly_derivative(c, b, ctx);
        fq_poly_derivative(b, b, ctx);

        result = (fq_poly_equal(b, c, ctx));
        if (!result)
        {
            flint_printf("FAIL:\n\n");
            flint_printf("a = "), fq_poly_print_pretty(a, "X", ctx), flint_printf("\n");
            flint_printf("b = "), fq_poly_print_pretty(b, "X", ctx), flint_printf("\n");
            flint_printf("c = "), fq_poly_print_pretty(c, "X", ctx), flint_printf("\n");
            abort();
        }

        fq_poly_clear(a, ctx);
        fq_poly_clear(b, ctx);
        fq_poly_clear(c, ctx);

        fq_ctx_clear(ctx);
    }

    /* Check constants have derivative zero */
    for (i = 0; i < 2000; i++)
    {
        fq_ctx_t ctx;

        fq_poly_t a, b;

        fq_ctx_randtest(ctx, state);
        fq_poly_init(a, ctx);
        fq_poly_init(b, ctx);

        fq_poly_randtest(a, state, n_randint(state, 2), ctx);
        fq_poly_derivative(b, a, ctx);

        result = (fq_poly_is_zero(b, ctx));
        if (!result)
        {
            flint_printf("FAIL:\n\n");
            flint_printf("a = "), fq_poly_print_pretty(a, "X", ctx), flint_printf("\n");
            flint_printf("b = "), fq_poly_print_pretty(b, "X", ctx), flint_printf("\n");
            abort();
        }

        fq_poly_clear(a, ctx);
        fq_poly_clear(b, ctx);

        fq_ctx_clear(ctx);
    }

    /* Check (f g)' == f' g + f g'  */
    for (i = 0; i < 2000; i++)
    {
        fq_ctx_t ctx;

        fq_poly_t a, b, c, d, lhs, rhs;

        fq_ctx_randtest(ctx, state);
        fq_poly_init(a, ctx);
        fq_poly_init(b, ctx);
        fq_poly_init(c, ctx);
        fq_poly_init(d, ctx);
        fq_poly_init(lhs, ctx);
        fq_poly_init(rhs, ctx);

        fq_poly_randtest(a, state, n_randint(state, 100), ctx);
        fq_poly_randtest(b, state, n_randint(state, 100), ctx);

        fq_poly_mul(lhs, a, b, ctx);
        fq_poly_derivative(lhs, lhs, ctx);
        fq_poly_derivative(c, a, ctx);
        fq_poly_derivative(d, b, ctx);
        fq_poly_mul(c, c, b, ctx);
        fq_poly_mul(d, a, d, ctx);
        fq_poly_add(rhs, c, d, ctx);

        result = (fq_poly_equal(lhs, rhs, ctx));
        if (!result)
        {
            flint_printf("FAIL:\n\n");
            flint_printf("a   = "), fq_poly_print_pretty(a,   "X", ctx), flint_printf("\n");
            flint_printf("b   = "), fq_poly_print_pretty(b,   "X", ctx), flint_printf("\n");
            flint_printf("c   = "), fq_poly_print_pretty(c,   "X", ctx), flint_printf("\n");
            flint_printf("d   = "), fq_poly_print_pretty(d,   "X", ctx), flint_printf("\n");
            flint_printf("lhs = "), fq_poly_print_pretty(lhs, "X", ctx), flint_printf("\n");
            flint_printf("rhs = "), fq_poly_print_pretty(rhs, "X", ctx), flint_printf("\n");
            abort();
        }

        fq_poly_clear(a, ctx);
        fq_poly_clear(b, ctx);
        fq_poly_clear(c, ctx);
        fq_poly_clear(d, ctx);
        fq_poly_clear(lhs, ctx);
        fq_poly_clear(rhs, ctx);

        fq_ctx_clear(ctx);
    }

    flint_randclear(state);
    _fmpz_cleanup();
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}

