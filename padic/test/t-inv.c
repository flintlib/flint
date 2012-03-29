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
#include "long_extras.h"
#include "padic.h"

int
main(void)
{
    int i, result;
    flint_rand_t state;

    printf("inv... ");
    fflush(stdout);

    flint_randinit(state);

/* PRIME p = 2 ***************************************************************/

    /* Check aliasing: a = a^{-1} (mod p^N) */
    for (i = 0; i < 1000; i++)
    {
        fmpz_t p;
        long N;
        padic_ctx_t ctx;

        padic_t a, d;

        fmpz_init(p);
        fmpz_set_ui(p, 2);
        N = z_randint(state, 1500) + 1;
        padic_ctx_init(ctx, p, N, PADIC_SERIES);

        padic_init(a, ctx);
        padic_init(d, ctx);

        padic_randtest_not_zero(a, state, ctx);

        padic_inv(d, a, ctx);
        padic_inv(a, a, ctx);

        result = (padic_equal(a, d, ctx));
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

    /* Check that correct only mod p^{N} */
    for (i = 0; i < 1000; i++)
    {
        fmpz_t p;
        long N;
        padic_ctx_t ctx;

        padic_t a, b, d;
        long v;

        fmpz_init(p);
        fmpz_set_ui(p, 2);
        N = n_randint(state, 1500) + 1;
        padic_ctx_init(ctx, p, N, PADIC_SERIES);

        padic_init(a, ctx);
        padic_init(b, ctx);
        padic_init(d, ctx);

        padic_randtest_not_zero(a, state, ctx);
        v = padic_val(a);

        if (-v < N) /* Otherwise, no precision left */
        {
            padic_ctx_t ctx2;

            padic_ctx_init(ctx2, p, N - FLINT_ABS(v), PADIC_SERIES);

            padic_inv(b, a, ctx);
            padic_mul(d, a, b, ctx);

            padic_reduce(d, ctx2);

            result = (padic_is_one(d, ctx2));
            if (!result)
            {
                printf("FAIL (a * a^{-1} == 1):\n\n");
                printf("a = "), padic_print(a, ctx), printf("\n");
                printf("b = "), padic_print(b, ctx), printf("\n");
                printf("d = "), padic_print(d, ctx2), printf("\n");
                abort();
            }
            padic_ctx_clear(ctx2);
        }

        padic_clear(a, ctx);
        padic_clear(b, ctx);
        padic_clear(d, ctx);

        fmpz_clear(p);
        padic_ctx_clear(ctx);
    }

/* PRIME p > 2 ***************************************************************/

    /* Check aliasing: a = a^{-1} (mod p^N) */
    for (i = 0; i < 1000; i++)
    {
        fmpz_t p;
        long N;
        padic_ctx_t ctx;

        padic_t a, d;

        fmpz_init(p);
        fmpz_set_ui(p, n_randprime(state, 5, 1));
        N = z_randint(state, 1500) + 1;
        padic_ctx_init(ctx, p, N, PADIC_SERIES);

        padic_init(a, ctx);
        padic_init(d, ctx);

        padic_randtest_not_zero(a, state, ctx);

        padic_inv(d, a, ctx);
        padic_inv(a, a, ctx);

        result = (padic_equal(a, d, ctx));
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

    /* Check that correct only mod p^{N} */
    for (i = 0; i < 1000; i++)
    {
        fmpz_t p;
        long N;
        padic_ctx_t ctx;

        padic_t a, b, d;
        long v;

        fmpz_init(p);
        fmpz_set_ui(p, n_randprime(state, 5, 1));
        N = n_randint(state, 1500) + 1;
        padic_ctx_init(ctx, p, N, PADIC_SERIES);

        padic_init(a, ctx);
        padic_init(b, ctx);
        padic_init(d, ctx);

        padic_randtest_not_zero(a, state, ctx);
        v = padic_val(a);

        if (-v < N) /* Otherwise, no precision left */
        {
            padic_ctx_t ctx2;

            padic_ctx_init(ctx2, p, N - FLINT_ABS(v), PADIC_SERIES);

            padic_inv(b, a, ctx);
            padic_mul(d, a, b, ctx);

            padic_reduce(d, ctx2);

            result = (padic_is_one(d, ctx2));
            if (!result)
            {
                printf("FAIL (a * a^{-1} == 1):\n\n");
                printf("a = "), padic_print(a, ctx), printf("\n");
                printf("b = "), padic_print(b, ctx), printf("\n");
                printf("d = "), padic_print(d, ctx2), printf("\n");
                abort();
            }
            padic_ctx_clear(ctx2);
        }

        padic_clear(a, ctx);
        padic_clear(b, ctx);
        padic_clear(d, ctx);

        fmpz_clear(p);
        padic_ctx_clear(ctx);
    }

    flint_randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}

