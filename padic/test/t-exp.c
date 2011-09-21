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

******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <mpir.h>
#include "flint.h"
#include "ulong_extras.h"
#include "padic.h"

int
main(void)
{
    int i, result;
    flint_rand_t state;

    printf("exp... ");
    fflush(stdout);

    flint_randinit(state);

/** p == 2 *******************************************************************/

    /* Check aliasing: a = exp(a) */
    for (i = 0; i < 10000; i++)
    {
        fmpz_t p;
        long N;
        padic_ctx_t ctx;

        padic_t a, b;
        int ans1, ans2;

        fmpz_init(p);
        fmpz_set_ui(p, 2);
        N = n_randint(state, 50) + 1;
        padic_ctx_init(ctx, p, N, PADIC_SERIES);

        padic_init(a, ctx);
        padic_init(b, ctx);

        padic_randtest(a, state, ctx);

        ans1 = padic_exp(b, a, ctx);
        ans2 = padic_exp(a, a, ctx);

        result = ((ans1 == ans2) && (!ans1 || padic_equal(a, b, ctx)));
        if (!result)
        {
            printf("FAIL (aliasing):\n\n");
            printf("a = "), padic_print(a, ctx), printf("\n");
            printf("b = "), padic_print(b, ctx), printf("\n");
            printf("ans1 = %d\n", ans1);
            printf("ans2 = %d\n", ans2);
            abort();
        }

        padic_clear(a, ctx);
        padic_clear(b, ctx);

        fmpz_clear(p);
        padic_ctx_clear(ctx);
    }

    /* Check aliasing: exp(a + b) == exp(a) exp(b) */
    for (i = 0; i < 10000; i++)
    {
        fmpz_t p;
        long N;
        padic_ctx_t ctx;

        padic_t a, b, c, d, e, f, g;
        int ans1, ans2, ans3;

        fmpz_init(p);
        fmpz_set_ui(p, 2);
        N = n_randint(state, 10) + 1;
        padic_ctx_init(ctx, p, N, PADIC_SERIES);

        padic_init(a, ctx);
        padic_init(b, ctx);
        padic_init(c, ctx);
        padic_init(d, ctx);
        padic_init(e, ctx);
        padic_init(f, ctx);
        padic_init(g, ctx);

        padic_randtest(a, state, ctx);
        padic_randtest(b, state, ctx);
        padic_add(c, a, b, ctx);

        ans1 = padic_exp(d, a, ctx);
        ans2 = padic_exp(e, b, ctx);
        padic_mul(f, d, e, ctx);

        ans3 = padic_exp(g, c, ctx);

        result = (!ans1 || !ans2 || (ans3 && padic_equal(f, g, ctx)));
        if (!result)
        {
            printf("FAIL (functional equation):\n\n");
            printf("a                 = "), padic_print(a, ctx), printf("\n");
            printf("b                 = "), padic_print(b, ctx), printf("\n");
            printf("c = a + b         = "), padic_print(c, ctx), printf("\n");
            printf("d = exp(a)        = "), padic_print(d, ctx), printf("\n");
            printf("e = exp(b)        = "), padic_print(e, ctx), printf("\n");
            printf("f = exp(a) exp(b) = "), padic_print(f, ctx), printf("\n");
            printf("g = exp(a + b)    = "), padic_print(g, ctx), printf("\n");
            abort();
        }

        padic_clear(a, ctx);
        padic_clear(b, ctx);
        padic_clear(c, ctx);
        padic_clear(d, ctx);
        padic_clear(e, ctx);
        padic_clear(f, ctx);
        padic_clear(g, ctx);

        fmpz_clear(p);
        padic_ctx_clear(ctx);
    }

/** p > 2 ********************************************************************/

    /* Check aliasing: a = exp(a) */
    for (i = 0; i < 10000; i++)
    {
        fmpz_t p;
        long N;
        padic_ctx_t ctx;

        padic_t a, b;
        int ans1, ans2;

        fmpz_init(p);
        fmpz_set_ui(p, n_randprime(state, 5, 1));
        N = n_randint(state, 50) + 1;
        padic_ctx_init(ctx, p, N, PADIC_SERIES);

        padic_init(a, ctx);
        padic_init(b, ctx);

        padic_randtest(a, state, ctx);

        ans1 = padic_exp(b, a, ctx);
        ans2 = padic_exp(a, a, ctx);

        result = ((ans1 == ans2) && (!ans1 || padic_equal(a, b, ctx)));
        if (!result)
        {
            printf("FAIL (aliasing):\n\n");
            printf("a = "), padic_print(a, ctx), printf("\n");
            printf("b = "), padic_print(b, ctx), printf("\n");
            printf("ans1 = %d\n", ans1);
            printf("ans2 = %d\n", ans2);
            abort();
        }

        padic_clear(a, ctx);
        padic_clear(b, ctx);

        fmpz_clear(p);
        padic_ctx_clear(ctx);
    }

    /* Check aliasing: exp(a + b) == exp(a) exp(b) */
    for (i = 0; i < 10000; i++)
    {
        fmpz_t p;
        long N;
        padic_ctx_t ctx;

        padic_t a, b, c, d, e, f, g;
        int ans1, ans2, ans3;

        fmpz_init(p);
        fmpz_set_ui(p, n_randprime(state, 5, 1));
        N = n_randint(state, 10) + 1;
        padic_ctx_init(ctx, p, N, PADIC_SERIES);

        padic_init(a, ctx);
        padic_init(b, ctx);
        padic_init(c, ctx);
        padic_init(d, ctx);
        padic_init(e, ctx);
        padic_init(f, ctx);
        padic_init(g, ctx);

        padic_randtest(a, state, ctx);
        padic_randtest(b, state, ctx);
        padic_add(c, a, b, ctx);

        ans1 = padic_exp(d, a, ctx);
        ans2 = padic_exp(e, b, ctx);
        padic_mul(f, d, e, ctx);

        ans3 = padic_exp(g, c, ctx);

        result = (!ans1 || !ans2 || (ans3 && padic_equal(f, g, ctx)));
        if (!result)
        {
            printf("FAIL (functional equation):\n\n");
            printf("a                 = "), padic_print(a, ctx), printf("\n");
            printf("b                 = "), padic_print(b, ctx), printf("\n");
            printf("c = a + b         = "), padic_print(c, ctx), printf("\n");
            printf("d = exp(a)        = "), padic_print(d, ctx), printf("\n");
            printf("e = exp(b)        = "), padic_print(e, ctx), printf("\n");
            printf("f = exp(a) exp(b) = "), padic_print(f, ctx), printf("\n");
            printf("g = exp(a + b)    = "), padic_print(g, ctx), printf("\n");
            abort();
        }

        padic_clear(a, ctx);
        padic_clear(b, ctx);
        padic_clear(c, ctx);
        padic_clear(d, ctx);
        padic_clear(e, ctx);
        padic_clear(f, ctx);
        padic_clear(g, ctx);

        fmpz_clear(p);
        padic_ctx_clear(ctx);
    }

    flint_randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}

