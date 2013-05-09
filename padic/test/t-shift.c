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

    Copyright (C) 2011, 2012 Sebastian Pancratz

******************************************************************************/

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
        len_t N;
        padic_ctx_t ctx;

        padic_t a, b, c;
        len_t v;

        fmpz_init_set_ui(p, n_randtest_prime(state, 0));
        N = n_randint(state, PADIC_TEST_PREC_MAX - PADIC_TEST_PREC_MIN) 
            + PADIC_TEST_PREC_MIN;
        padic_ctx_init(ctx, p, FLINT_MAX(0, N-10), FLINT_MAX(0, N+10), PADIC_SERIES);

        padic_init2(a, N);
        padic_init2(b, N);
        padic_init2(c, N);

        padic_randtest(a, state, ctx);
        v = z_randint(state, (FLINT_ABS(N) + 4) / 3);

        padic_set(b, a, ctx);
        padic_shift(c, b, v, ctx);
        padic_shift(b, b, v, ctx);

        result = (padic_equal(b, c));
        if (!result)
        {
            printf("FAIL:\n\n");
            printf("a = "), padic_print(a, ctx), printf("\n");
            printf("b = "), padic_print(b, ctx), printf("\n");
            printf("c = "), padic_print(c, ctx), printf("\n");
            abort();
        }

        padic_clear(a);
        padic_clear(b);
        padic_clear(c);

        fmpz_clear(p);
        padic_ctx_clear(ctx);
    }

    /* Check that (a * b) * c == a * (b * c), correct only mod p^{N-v} */
    for (i = 0; i < 10000; i++)
    {
        fmpz_t p;
        len_t N;
        padic_ctx_t ctx;

        padic_t a, b, c;
        len_t v, v1, v2;

        fmpz_init_set_ui(p, n_randtest_prime(state, 0));
        N = n_randint(state, PADIC_TEST_PREC_MAX - PADIC_TEST_PREC_MIN) 
            + PADIC_TEST_PREC_MIN;
        padic_ctx_init(ctx, p, FLINT_MAX(0, N-10), FLINT_MAX(0, N+10), PADIC_SERIES);

        padic_init2(a, N);
        padic_init2(b, N);
        padic_init2(c, N);

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
            len_t N2 = (v >= 0) ? N : N + v;

            padic_prec(b) = N2;
            padic_prec(c) = N2;

            padic_reduce(b, ctx);
            padic_reduce(c, ctx);

            result = (padic_equal(b, c));
            if (!result)
            {
                printf("FAIL:\n\n");
                printf("a = "), padic_print(a, ctx), printf("\n");
                printf("b = "), padic_print(b, ctx), printf("\n");
                printf("c = "), padic_print(c, ctx), printf("\n");
                abort();
            }
        }

        padic_clear(a);
        padic_clear(b);
        padic_clear(c);

        fmpz_clear(p);
        padic_ctx_clear(ctx);
    }

    flint_randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}

