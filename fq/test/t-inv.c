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

    flint_printf("inv... ");
    fflush(stdout);

    flint_randinit(state);

    /* Check aliasing: a = ~a */
    for (i = 0; i < 2000; i++)
    {
        fq_ctx_t ctx;
        fq_t a, b, c;

        fq_ctx_randtest(ctx, state);

        fq_init(a, ctx);
        fq_init(b, ctx);
        fq_init(c, ctx);

        fq_randtest_not_zero(a, state, ctx);
        fq_set(b, a, ctx);

        fq_inv(c, b, ctx);
        fq_inv(b, b, ctx);

        result = (fq_equal(b, c, ctx));
        if (!result)
        {
            flint_printf("FAIL (aliasing):\n\n");
            flint_printf("a = "), fq_print_pretty(a, ctx), flint_printf("\n");
            flint_printf("b = "), fq_print_pretty(b, ctx), flint_printf("\n");
            flint_printf("c = "), fq_print_pretty(c, ctx), flint_printf("\n");
            fq_ctx_print(ctx);
            abort();
        }

        fq_clear(a, ctx);
        fq_clear(b, ctx);
        fq_clear(c, ctx);

        fq_ctx_clear(ctx);
    }

    /* Check a * ~a == 1 for units*/
    for (i = 0; i < 2000; i++)
    {
        fq_ctx_t ctx;
        fq_t a, b, c;

        fq_ctx_randtest(ctx, state);

        fq_init(a, ctx);
        fq_init(b, ctx);
        fq_init(c, ctx);

        fq_randtest_not_zero(a, state, ctx);

        fq_inv(b, a, ctx);
        fq_mul(c, a, b, ctx);

        result = (fq_is_one(c, ctx));
        if (!result)
        {
            flint_printf("FAIL (a * (~a) == 1):\n\n");
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

    flint_randclear(state);
    _fmpz_cleanup();
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}

