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

    flint_printf("pow... ");
    fflush(stdout);

    flint_randinit(state);

    /* Check aliasing: a = a^e */
    for (i = 0; i < 1000; i++)
    {
        fq_ctx_t ctx;
        fq_t a, b;
        fmpz_t e;

        fq_ctx_randtest(ctx, state);

        fq_init(a, ctx);
        fq_init(b, ctx);
        fmpz_init(e);

        fq_randtest(a, state, ctx);
        fmpz_randtest_unsigned(e, state, 6);

        fq_pow(b, a, e, ctx);
        fq_pow(a, a, e, ctx);

        result = (fq_equal(a, b, ctx));
        if (!result)
        {
            flint_printf("FAIL (alias):\n\n");
            flint_printf("a = "), fq_print_pretty(a, ctx), flint_printf("\n");
            flint_printf("b = "), fq_print_pretty(b, ctx), flint_printf("\n");
            abort();
        }

        fq_clear(a, ctx);
        fq_clear(b, ctx);
        fmpz_clear(e);

        fq_ctx_clear(ctx);
    }

    /* Compare with multiplication, for integral values */
    for (i = 0; i < 1000; i++)
    {
        fq_ctx_t ctx;

        fq_t a, b, c;
        fmpz_t e, f;

        fq_ctx_randtest(ctx, state);
        
        fq_init(a, ctx);
        fq_init(b, ctx);
        fq_init(c, ctx);
        fmpz_init(f);
        fmpz_init(e);

        fq_randtest(a, state, ctx);
        fmpz_randtest_unsigned(e, state, 6);

        fq_pow(b, a, e, ctx);
        fq_one(c, ctx);
        for (fmpz_one(f); fmpz_cmp(f, e) <= 0; fmpz_add_ui(f, f, 1))
        {
            fq_mul(c, c, a, ctx);
        }

        result = (fq_equal(b, c, ctx));
        if (!result)
        {
            flint_printf("FAIL (cmp with mul):\n\n");
            flint_printf("a = "), fq_print_pretty(a, ctx), flint_printf("\n");
            flint_printf("b = "), fq_print_pretty(b, ctx), flint_printf("\n");
            flint_printf("c = "), fq_print_pretty(b, ctx), flint_printf("\n");
            flint_printf("e = "), fmpz_print(e), flint_printf("\n");
            abort();
        }

        fq_clear(a, ctx);
        fq_clear(b, ctx);
        fq_clear(c, ctx);
        fmpz_clear(e);
        fmpz_clear(f);

        fq_ctx_clear(ctx);
    }

    flint_randclear(state);
    _fmpz_cleanup();
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}

