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

    flint_printf("scalar_addmul_fq... ");
    fflush(stdout);

    flint_randinit(state);

    /* Check aliasing */
    for (i = 0; i < 2000; i++)
    {
        long len;
        fq_ctx_t ctx;

        fq_poly_t a, b, c;
        fq_t x;

        len = n_randint(state, 15) + 1;
        fq_ctx_randtest(ctx, state);

        fq_poly_init(a, ctx);
        fq_poly_init(b, ctx);
        fq_poly_init(c, ctx);
        fq_init(x, ctx);

        fq_poly_randtest(a, state, len, ctx);
        fq_poly_set(b, a, ctx);
        fq_poly_set(c, a, ctx);

        fq_poly_scalar_addmul_fq(b, a, x, ctx);
        fq_poly_scalar_addmul_fq(a, a, x, ctx);

        result = (fq_poly_equal(a, b, ctx));
        if (!result)
        {
            flint_printf("FAIL (alias):\n\n");
            flint_printf("a = "), fq_poly_print_pretty(a, "X", ctx), flint_printf("\n");
            flint_printf("b = "), fq_poly_print_pretty(b, "X", ctx), flint_printf("\n");
            flint_printf("c = "), fq_poly_print_pretty(c, "X", ctx), flint_printf("\n");
            flint_printf("x = "), fq_print_pretty(x, ctx), flint_printf("\n");
            abort();
        }

        fq_poly_clear(a, ctx);
        fq_poly_clear(b, ctx);
        fq_poly_clear(c, ctx);
        fq_clear(x, ctx);

        fq_ctx_clear(ctx);
    }

    /* Check that b += x*a is the same as c = b + x*a */
    for (i = 0; i < 2000; i++)
    {
        long len;
        fq_ctx_t ctx;

        fq_poly_t a, b, c;
        fq_t x;

        len = n_randint(state, 15) + 1;
        fq_ctx_randtest(ctx, state);

        fq_poly_init(a, ctx);
        fq_poly_init(b, ctx);
        fq_poly_init(c, ctx);
        fq_init(x, ctx);

        fq_poly_randtest(a, state, len, ctx);
        fq_poly_randtest(b, state, len, ctx);
        fq_poly_set(a, c, ctx);

        fq_poly_scalar_mul_fq(c, a, x, ctx);
        fq_poly_add(c, b, c, ctx);
        fq_poly_scalar_addmul_fq(b, a, x, ctx);

        result = (fq_poly_equal(b, c, ctx));
        if (!result)
        {
            flint_printf("FAIL:\n\n");
            flint_printf("a = "), fq_poly_print_pretty(a, "X", ctx), flint_printf("\n");
            flint_printf("b = "), fq_poly_print_pretty(b, "X", ctx), flint_printf("\n");
            flint_printf("c = "), fq_poly_print_pretty(c, "X", ctx), flint_printf("\n");
            flint_printf("x = "), fq_print_pretty(x, ctx), flint_printf("\n");
            abort();
        }

        fq_poly_clear(a, ctx);
        fq_poly_clear(b, ctx);
        fq_poly_clear(c, ctx);
        fq_clear(x, ctx);

        fq_ctx_clear(ctx);
    }

    flint_randclear(state);
    _fmpz_cleanup();
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}

