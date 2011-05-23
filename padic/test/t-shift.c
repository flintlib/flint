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

    printf("shift... ");
    fflush(stdout);

    flint_randinit(state);

    /* Check aliasing */
    for (i = 0; i < 10000; i++)
    {
        fmpz_t p;
        long N;
        padic_ctx_t ctx;

        padic_t a, b, c;
        long v;

        fmpz_init(p);
        fmpz_set_ui(p, n_randprime(state, 5, 1));
        N = z_randint(state, 100);
        padic_ctx_init(ctx, p, N, PADIC_SERIES);

        padic_init(a, ctx);
        padic_init(b, ctx);
        padic_init(c, ctx);

        padic_randtest(a, state, ctx);
        v = z_randint(state, (FLINT_ABS(N) + 4) / 3);

        padic_set(b, a, ctx);
        padic_shift(c, b, v, ctx);
        padic_shift(b, b, v, ctx);

        result = (padic_equal(b, c, ctx));
        if (!result)
        {
            printf("FAIL:\n\n");
            printf("a = "), padic_print(a, ctx), printf("\n");
            printf("b = "), padic_print(b, ctx), printf("\n");
            printf("c = "), padic_print(c, ctx), printf("\n");
            abort();
        }

        padic_clear(a, ctx);
        padic_clear(b, ctx);
        padic_clear(c, ctx);

        fmpz_clear(p);
        padic_ctx_clear(ctx);
    }

    /* Check that (a * b) * c == a * (b * c), correct only mod p^{N-v} */
    for (i = 0; i < 10000; i++)
    {
        fmpz_t p;
        long N;
        padic_ctx_t ctx;

        padic_t a, b, c;
        long v, v1, v2;

        fmpz_init(p);
        fmpz_set_ui(p, n_randprime(state, 5, 1));
        N = z_randint(state, 100);
        padic_ctx_init(ctx, p, N, PADIC_SERIES);

        padic_init(a, ctx);
        padic_init(b, ctx);
        padic_init(c, ctx);

        padic_randtest(a, state, ctx);
        v1 = z_randint(state, (FLINT_ABS(N) + 4) / 3);
        v2 = z_randint(state, (FLINT_ABS(N) + 4) / 3);

        padic_shift(b, a, v1, ctx);
        padic_shift(b, b, v2, ctx);

        padic_shift(c, a, v2, ctx);
        padic_shift(c, c, v1, ctx);

        v = FLINT_MIN(v1, v2);
        v = FLINT_MIN(v, padic_val(a));
        v = FLINT_MIN(v, 0);

        if ((v >= 0) || (-v < N)) /* Otherwise, no precision left */
        {
            padic_ctx_t ctx2;

            padic_ctx_init(ctx2, p, (v >= 0) ? N : N + v, PADIC_SERIES);

            padic_normalise(b, ctx2);
            padic_normalise(c, ctx2);

            result = (padic_equal(b, c, ctx2));
            if (!result)
            {
                printf("FAIL:\n\n");
                printf("a = "), padic_print(a, ctx), printf("\n");
                printf("b = "), padic_print(b, ctx2), printf("\n");
                printf("c = "), padic_print(c, ctx2), printf("\n");
                abort();
            }
        }

        padic_clear(a, ctx);
        padic_clear(b, ctx);
        padic_clear(c, ctx);

        fmpz_clear(p);
        padic_ctx_clear(ctx);
    }

    flint_randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}

