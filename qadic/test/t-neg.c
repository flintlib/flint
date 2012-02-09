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

#include "qadic.h"
#include "ulong_extras.h"
#include "long_extras.h"

int
main(void)
{
    int i, result;
    flint_rand_t state;

    printf("neg... ");
    fflush(stdout);

    flint_randinit(state);

    /* Check aliasing: a = -a */
    for (i = 0; i < 2000; i++)
    {
        fmpz_t p;
        long d, N;
        qadic_ctx_t ctx;

        qadic_t a, b, c;

        fmpz_init(p);
        fmpz_set_ui(p, n_randprime(state, 2 + n_randint(state, 3), 1));
        d = n_randint(state, 10) + 1;
        N = z_randint(state, 50) + 1;
        qadic_ctx_init_conway(ctx, p, d, N, "a", PADIC_SERIES);

        qadic_init(a);
        qadic_init(b);
        qadic_init(c);

        qadic_randtest(a, state, ctx);
        qadic_set(b, a);

        qadic_neg(c, b, ctx);
        qadic_neg(b, b, ctx);

        result = (qadic_equal(b, c));
        if (!result)
        {
            printf("FAIL:\n\n");
            printf("a = "), qadic_print_pretty(a, ctx), printf("\n");
            printf("b = "), qadic_print_pretty(b, ctx), printf("\n");
            printf("c = "), qadic_print_pretty(c, ctx), printf("\n");
            abort();
        }

        qadic_clear(a);
        qadic_clear(b);
        qadic_clear(c);

        fmpz_clear(p);
        qadic_ctx_clear(ctx);
    }

    /* Check a - b == a + (-b) */
    for (i = 0; i < 2000; i++)
    {
        fmpz_t p;
        long d, N;
        qadic_ctx_t ctx;

        qadic_t a, b, c1, c2;

        fmpz_init(p);
        fmpz_set_ui(p, n_randprime(state, 2 + n_randint(state, 3), 1));
        d = n_randint(state, 10) + 1;
        N = z_randint(state, 50) + 1;
        qadic_ctx_init_conway(ctx, p, d, N, "a", PADIC_SERIES);

        qadic_init(a);
        qadic_init(b);
        qadic_init(c1);
        qadic_init(c2);

        qadic_randtest(a, state, ctx);
        qadic_randtest(b, state, ctx);

        qadic_sub(c1, a, b, ctx);
        qadic_neg(c2, b, ctx);
        qadic_add(c2, a, c2, ctx);

        result = (qadic_equal(c1, c2));
        if (!result)
        {
            printf("FAIL:\n\n");
            printf("a = "), qadic_print_pretty(a, ctx), printf("\n");
            printf("b = "), qadic_print_pretty(b, ctx), printf("\n");
            printf("c1 = "), qadic_print_pretty(c1, ctx), printf("\n");
            printf("c2 = "), qadic_print_pretty(c2, ctx), printf("\n");
            abort();
        }

        qadic_clear(a);
        qadic_clear(b);
        qadic_clear(c1);
        qadic_clear(c2);

        fmpz_clear(p);
        qadic_ctx_clear(ctx);
    }

    flint_randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}

