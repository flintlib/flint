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
#include "fmpz.h"
#include "padic.h"
#include "ulong_extras.h"
#include "long_extras.h"

int main(void)
{
    int i, result;
    flint_rand_t state;

    printf("sqrt... ");
    fflush(stdout);

    flint_randinit(state);

/* PRIME p = 2 ***************************************************************/

    /* Check aliasing: a = sqrt(a) */
    for (i = 0; i < 10000; i++)
    {
        fmpz_t p;
        long N;
        padic_ctx_t ctx;

        int ans1, ans2;
        padic_t a, d;

        fmpz_init(p);
        fmpz_set_ui(p, 2);
        N = z_randint(state, 100);
        padic_ctx_init(ctx, p, N, PADIC_SERIES);

        padic_init(a, ctx);
        padic_init(d, ctx);

        padic_randtest(a, state, ctx);

        ans1 = padic_sqrt(d, a, ctx);
        ans2 = padic_sqrt(a, a, ctx);

        result = ((ans1 == ans2) && (!ans1 || padic_equal(a, d, ctx)));
        if (!result)
        {
            printf("FAIL (aliasing):\n\n");
            printf("p = "), fmpz_print(p), printf("\n");
            printf("a = "), padic_print(a, ctx), printf("\n");
            printf("d = "), padic_print(d, ctx), printf("\n");
            abort();
        }

        padic_clear(a, ctx);
        padic_clear(d, ctx);

        fmpz_clear(p);
        padic_ctx_clear(ctx);
    }

    /* Test random elements */
    for (i = 0; i < 10000; i++)
    {
        fmpz_t p;
        long N;
        padic_ctx_t ctx;

        int ans;
        padic_t a, b, d;

        fmpz_init(p);
        fmpz_set_ui(p, 2);
        N = z_randint(state, 100);
        padic_ctx_init(ctx, p, N, PADIC_SERIES);

        padic_init(a, ctx);
        padic_init(b, ctx);
        padic_init(d, ctx);

        padic_randtest(a, state, ctx);

        ans = padic_sqrt(b, a, ctx);

        padic_mul(d, b, b, ctx);

        if (ans && padic_val(a) < 0)
        {
            padic_ctx_t ctx2;
            padic_t a2, d2;

            padic_ctx_init(ctx2, p, N + padic_val(a), PADIC_SERIES);

            padic_init(a2, ctx2);
            padic_init(d2, ctx2);
            padic_set(a2, a, ctx2);
            padic_set(d2, d, ctx2);
            padic_reduce(a2, ctx2);
            padic_reduce(d2, ctx2);

            result = (padic_equal(a2, d2, ctx2));
            if (!result)
            {
                printf("FAIL (random elements):\n\n");
                printf("a = "), padic_print(a, ctx), printf("\n");
                printf("b = "), padic_print(b, ctx), printf("\n");
                printf("d = "), padic_print(d, ctx), printf("\n");
                printf("a2 = "), padic_print(a2, ctx2), printf("\n");
                printf("d2 = "), padic_print(d2, ctx2), printf("\n");
                printf("p = "), fmpz_print(p), printf("\n");
                printf("ans = %d\n", ans);
                abort();
            }

            padic_clear(a2, ctx2);
            padic_clear(d2, ctx2);
            padic_ctx_clear(ctx2);
        }
        else
        {
            result = (!ans || padic_equal(a, d, ctx));
            if (!result)
            {
                printf("FAIL (random elements):\n\n");
                printf("a = "), padic_print(a, ctx), printf("\n");
                printf("b = "), padic_print(b, ctx), printf("\n");
                printf("d = "), padic_print(d, ctx), printf("\n");
                printf("p = "), fmpz_print(p), printf("\n");
                printf("ans = %d\n", ans);
                abort();
            }
        }

        padic_clear(a, ctx);
        padic_clear(b, ctx);
        padic_clear(d, ctx);

        fmpz_clear(p);
        padic_ctx_clear(ctx);
    }

    /* Test random squares */
    for (i = 0; i < 10000; i++)
    {
        fmpz_t p;
        long N;
        padic_ctx_t ctx;

        int ans;
        padic_t a, b, c, d;

        fmpz_init(p);
        fmpz_set_ui(p, 2);
        N = z_randint(state, 100);
        padic_ctx_init(ctx, p, N, PADIC_SERIES);

        padic_init(a, ctx);
        padic_init(b, ctx);
        padic_init(c, ctx);
        padic_init(d, ctx);

        padic_randtest(b, state, ctx);
        padic_mul(a, b, b, ctx);

        ans = padic_sqrt(c, a, ctx);

        padic_mul(d, c, c, ctx);

        if (ans && padic_val(a) < 0)
        {
            padic_ctx_t ctx2;
            padic_t a2, d2;

            padic_ctx_init(ctx2, p, N + padic_val(a), PADIC_SERIES);

            padic_init(a2, ctx2);
            padic_init(d2, ctx2);
            padic_set(a2, a, ctx2);
            padic_set(d2, d, ctx2);
            padic_reduce(a2, ctx2);
            padic_reduce(d2, ctx2);

            result = (padic_equal(a2, d2, ctx2));
            if (!result)
            {
                printf("FAIL (random elements):\n\n");
                printf("a = "), padic_print(a, ctx), printf("\n");
                printf("b = "), padic_print(b, ctx), printf("\n");
                printf("c = "), padic_print(c, ctx), printf("\n");
                printf("d = "), padic_print(d, ctx), printf("\n");
                printf("a2 = "), padic_print(a2, ctx2), printf("\n");
                printf("d2 = "), padic_print(d2, ctx2), printf("\n");
                printf("p = "), fmpz_print(p), printf("\n");
                printf("ans = %d\n", ans);
                abort();
            }

            padic_clear(a2, ctx2);
            padic_clear(d2, ctx2);
            padic_ctx_clear(ctx2);
        }
        else
        {
            result = (ans && padic_equal(a, d, ctx));
            if (!result)
            {
                printf("FAIL (random squares):\n\n");
                printf("a = "), padic_print(a, ctx), printf("\n");
                printf("b = "), padic_print(b, ctx), printf("\n");
                printf("c = "), padic_print(c, ctx), printf("\n");
                printf("d = "), padic_print(d, ctx), printf("\n");
                printf("p = "), fmpz_print(p), printf("\n");
                printf("ans = %d\n", ans);
                abort();
            }
        }

        padic_clear(a, ctx);
        padic_clear(b, ctx);
        padic_clear(c, ctx);
        padic_clear(d, ctx);

        fmpz_clear(p);
        padic_ctx_clear(ctx);
    }

/* PRIME p > 2 ***************************************************************/

    /* Check aliasing: a = sqrt(a) */
    for (i = 0; i < 10000; i++)
    {
        fmpz_t p;
        long N;
        padic_ctx_t ctx;

        int ans1, ans2;
        padic_t a, d;

        fmpz_init(p);
        fmpz_set_ui(p, n_randprime(state, 5, 1));
        N = z_randint(state, 100);
        padic_ctx_init(ctx, p, N, PADIC_SERIES);

        padic_init(a, ctx);
        padic_init(d, ctx);

        padic_randtest(a, state, ctx);

        ans1 = padic_sqrt(d, a, ctx);
        ans2 = padic_sqrt(a, a, ctx);

        result = ((ans1 == ans2) && (!ans1 || padic_equal(a, d, ctx)));
        if (!result)
        {
            printf("FAIL (aliasing):\n\n");
            printf("a = "), padic_print(a, ctx), printf("\n");
            printf("d = "), padic_print(d, ctx), printf("\n");
            abort();
        }

        padic_clear(a, ctx);
        padic_clear(d, ctx);

        fmpz_clear(p);
        padic_ctx_clear(ctx);
    }

    /* Test random elements */
    for (i = 0; i < 10000; i++)
    {
        fmpz_t p;
        long N;
        padic_ctx_t ctx;

        int ans;
        padic_t a, b, d;

        fmpz_init(p);
        fmpz_set_ui(p, n_randprime(state, 5, 1));
        N = z_randint(state, 100);
        padic_ctx_init(ctx, p, N, PADIC_SERIES);

        padic_init(a, ctx);
        padic_init(b, ctx);
        padic_init(d, ctx);

        padic_randtest(a, state, ctx);

        ans = padic_sqrt(b, a, ctx);

        padic_mul(d, b, b, ctx);

        if (ans && padic_val(a) < 0)
        {
            padic_ctx_t ctx2;
            padic_t a2, d2;

            padic_ctx_init(ctx2, p, N + padic_val(a), PADIC_SERIES);

            padic_init(a2, ctx2);
            padic_init(d2, ctx2);
            padic_set(a2, a, ctx2);
            padic_set(d2, d, ctx2);
            padic_reduce(a2, ctx2);
            padic_reduce(d2, ctx2);

            result = (padic_equal(a2, d2, ctx2));
            if (!result)
            {
                printf("FAIL (random elements):\n\n");
                printf("a = "), padic_print(a, ctx), printf("\n");
                printf("b = "), padic_print(b, ctx), printf("\n");
                printf("d = "), padic_print(d, ctx), printf("\n");
                printf("a2 = "), padic_print(a2, ctx2), printf("\n");
                printf("d2 = "), padic_print(d2, ctx2), printf("\n");
                printf("p = "), fmpz_print(p), printf("\n");
                printf("ans = %d\n", ans);
                abort();
            }

            padic_clear(a2, ctx2);
            padic_clear(d2, ctx2);
            padic_ctx_clear(ctx2);
        }
        else
        {
            result = (!ans || padic_equal(a, d, ctx));
            if (!result)
            {
                printf("FAIL (random elements):\n\n");
                printf("a = "), padic_print(a, ctx), printf("\n");
                printf("b = "), padic_print(b, ctx), printf("\n");
                printf("d = "), padic_print(d, ctx), printf("\n");
                printf("p = "), fmpz_print(p), printf("\n");
                printf("ans = %d\n", ans);
                abort();
            }
        }

        padic_clear(a, ctx);
        padic_clear(b, ctx);
        padic_clear(d, ctx);

        fmpz_clear(p);
        padic_ctx_clear(ctx);
    }

    /* Test random squares */
    for (i = 0; i < 10000; i++)
    {
        fmpz_t p;
        long N;
        padic_ctx_t ctx;

        int ans;
        padic_t a, b, c, d;

        fmpz_init(p);
        fmpz_set_ui(p, n_randprime(state, 5, 1));
        N = z_randint(state, 100);
        padic_ctx_init(ctx, p, N, PADIC_SERIES);

        padic_init(a, ctx);
        padic_init(b, ctx);
        padic_init(c, ctx);
        padic_init(d, ctx);

        padic_randtest(b, state, ctx);
        padic_mul(a, b, b, ctx);

        ans = padic_sqrt(c, a, ctx);

        padic_mul(d, c, c, ctx);

        if (ans && padic_val(a) < 0)
        {
            padic_ctx_t ctx2;
            padic_t a2, d2;

            padic_ctx_init(ctx2, p, N + padic_val(a), PADIC_SERIES);

            padic_init(a2, ctx2);
            padic_init(d2, ctx2);
            padic_set(a2, a, ctx2);
            padic_set(d2, d, ctx2);
            padic_reduce(a2, ctx2);
            padic_reduce(d2, ctx2);

            result = (padic_equal(a2, d2, ctx2));
            if (!result)
            {
                printf("FAIL (random elements):\n\n");
                printf("a = "), padic_print(a, ctx), printf("\n");
                printf("b = "), padic_print(b, ctx), printf("\n");
                printf("c = "), padic_print(c, ctx), printf("\n");
                printf("d = "), padic_print(d, ctx), printf("\n");
                printf("a2 = "), padic_print(a2, ctx2), printf("\n");
                printf("d2 = "), padic_print(d2, ctx2), printf("\n");
                printf("p = "), fmpz_print(p), printf("\n");
                printf("ans = %d\n", ans);
                abort();
            }

            padic_clear(a2, ctx2);
            padic_clear(d2, ctx2);
            padic_ctx_clear(ctx2);
        }
        else
        {
            result = (ans && padic_equal(a, d, ctx));
            if (!result)
            {
                printf("FAIL (random squares):\n\n");
                printf("a = "), padic_print(a, ctx), printf("\n");
                printf("b = "), padic_print(b, ctx), printf("\n");
                printf("c = "), padic_print(c, ctx), printf("\n");
                printf("d = "), padic_print(d, ctx), printf("\n");
                printf("p = "), fmpz_print(p), printf("\n");
                printf("ans = %d\n", ans);
                abort();
            }
        }

        padic_clear(a, ctx);
        padic_clear(b, ctx);
        padic_clear(c, ctx);
        padic_clear(d, ctx);

        fmpz_clear(p);
        padic_ctx_clear(ctx);
    }

    flint_randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}
