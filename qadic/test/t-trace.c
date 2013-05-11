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

#include "qadic.h"
#include "ulong_extras.h"
#include "long_extras.h"

int
main(void)
{
    int i, result;
    flint_rand_t state;

    printf("trace... ");
    fflush(stdout);

    flint_randinit(state);

    /* Compare with sum of Galois conjugates */
    for (i = 0; i < 2000; i++)
    {
        fmpz_t p;
        long d, N;
        qadic_ctx_t ctx;

        qadic_t a, b, c;
        padic_t x, y;
        long j;
        int ans;

        fmpz_init_set_ui(p, n_randprime(state, 2 + n_randint(state, 3), 1));
        d = n_randint(state, 10) + 1;
        N = z_randint(state, 50) + 1;
        qadic_ctx_init_conway(ctx, p, d, FLINT_MAX(0,N-10), FLINT_MAX(0,N+10), "a", PADIC_SERIES);

        qadic_init2(a, N);
        qadic_init2(b, N);
        qadic_init2(c, N);
        padic_init2(x, N);
        padic_init2(y, N);

        qadic_randtest(a, state, ctx);

        qadic_trace(x, a, ctx);

        qadic_zero(b);
        for (j = 0; j < d; j++)
        {
            qadic_frobenius(c, a, j, ctx);
            qadic_add(b, b, c, ctx);
        }
        ans = qadic_get_padic(y, b, ctx);

        result = (ans && padic_equal(x, y));
        if (!result)
        {
            printf("FAIL:\n\n");
            printf("a = "), qadic_print_pretty(a, ctx), printf("\n");
            printf("b = "), qadic_print_pretty(b, ctx), printf("\n");
            printf("x = "), padic_print(x, &ctx->pctx), printf("\n");
            printf("y = "), padic_print(y, &ctx->pctx), printf("\n");
            for (j = 0; j < d; j++)
            {
                qadic_frobenius(c, a, j, ctx);
                printf("sigma^%ld = ", j), qadic_print_pretty(c, ctx), printf("\n");
            }
            printf("ans = %d\n", ans);
            qadic_ctx_print(ctx);
            abort();
        }

        qadic_clear(a);
        qadic_clear(b);
        qadic_clear(c);
        padic_clear(x);
        padic_clear(y);

        fmpz_clear(p);
        qadic_ctx_clear(ctx);
    }

    flint_randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}

