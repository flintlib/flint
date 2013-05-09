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
        len_t N;
        padic_ctx_t ctx;

        int ans1, ans2;
        padic_t a, d;

        fmpz_init_set_ui(p, 2);
        N = n_randint(state, PADIC_TEST_PREC_MAX - PADIC_TEST_PREC_MIN) 
            + PADIC_TEST_PREC_MIN;
        padic_ctx_init(ctx, p, FLINT_MAX(0, N-10), FLINT_MAX(0, N+10), PADIC_SERIES);

        padic_init2(a, N);
        padic_init2(d, N);

        padic_randtest(a, state, ctx);

        ans1 = padic_sqrt(d, a, ctx);
        ans2 = padic_sqrt(a, a, ctx);

        result = ((ans1 == ans2) && (!ans1 || padic_equal(a, d)));
        if (!result)
        {
            printf("FAIL (aliasing):\n\n");
            printf("p = "), fmpz_print(p), printf("\n");
            printf("a = "), padic_print(a, ctx), printf("\n");
            printf("d = "), padic_print(d, ctx), printf("\n");
            abort();
        }

        padic_clear(a);
        padic_clear(d);

        fmpz_clear(p);
        padic_ctx_clear(ctx);
    }

    /* Test random elements */
    for (i = 0; i < 10000; i++)
    {
        fmpz_t p;
        len_t N;
        padic_ctx_t ctx;

        int ans;
        padic_t a, b, d;

        fmpz_init_set_ui(p, 2);
        N = n_randint(state, PADIC_TEST_PREC_MAX - PADIC_TEST_PREC_MIN) 
            + PADIC_TEST_PREC_MIN;
        padic_ctx_init(ctx, p, FLINT_MAX(0, N-10), FLINT_MAX(0, N+10), PADIC_SERIES);

        padic_init2(a, N);
        padic_init2(b, N);
        padic_init2(d, N);

        padic_randtest(a, state, ctx);

        ans = padic_sqrt(b, a, ctx);

        padic_mul(d, b, b, ctx);

        if (ans && padic_val(a) < 0)
        {
            len_t N2 = N + padic_val(a);
            padic_t a2, d2;

            padic_init2(a2, N2);
            padic_init2(d2, N2);

            padic_set(a2, a, ctx);
            padic_set(d2, d, ctx);

            result = (padic_equal(a2, d2));
            if (!result)
            {
                printf("FAIL (random elements):\n\n");
                printf("a  = "), padic_print(a, ctx),  printf("\n");
                printf("b  = "), padic_print(b, ctx),  printf("\n");
                printf("d  = "), padic_print(d, ctx),  printf("\n");
                printf("a2 = "), padic_print(a2, ctx), printf("\n");
                printf("d2 = "), padic_print(d2, ctx), printf("\n");
                printf("p = "), fmpz_print(p), printf("\n");
                printf("ans = %d\n", ans);
                abort();
            }

            padic_clear(a2);
            padic_clear(d2);
        }
        else
        {
            result = (!ans || padic_equal(a, d));
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

        padic_clear(a);
        padic_clear(b);
        padic_clear(d);

        fmpz_clear(p);
        padic_ctx_clear(ctx);
    }

    /* Test random squares */
    for (i = 0; i < 10000; i++)
    {
        fmpz_t p;
        len_t N;
        padic_ctx_t ctx;

        int ans;
        padic_t a, b, c, d;

        fmpz_init_set_ui(p, 2);
        N = n_randint(state, PADIC_TEST_PREC_MAX - PADIC_TEST_PREC_MIN) 
            + PADIC_TEST_PREC_MIN;
        padic_ctx_init(ctx, p, FLINT_MAX(0, N-10), FLINT_MAX(0, N+10), PADIC_SERIES);

        padic_init2(a, N);
        padic_init2(b, N);
        padic_init2(c, N);
        padic_init2(d, N);

        padic_randtest(b, state, ctx);
        padic_mul(a, b, b, ctx);

        ans = padic_sqrt(c, a, ctx);

        padic_mul(d, c, c, ctx);

        if (ans && padic_val(a) < 0)
        {
            len_t N2 = N + padic_val(a);
            padic_t a2, d2;

            padic_init2(a2, N2);
            padic_init2(d2, N2);

            padic_set(a2, a, ctx);
            padic_set(d2, d, ctx);

            result = (padic_equal(a2, d2));
            if (!result)
            {
                printf("FAIL (random elements):\n\n");
                printf("a  = "), padic_print(a, ctx),  printf("\n");
                printf("b  = "), padic_print(b, ctx),  printf("\n");
                printf("c  = "), padic_print(c, ctx),  printf("\n");
                printf("d  = "), padic_print(d, ctx),  printf("\n");
                printf("a2 = "), padic_print(a2, ctx), printf("\n");
                printf("d2 = "), padic_print(d2, ctx), printf("\n");
                printf("p = "), fmpz_print(p), printf("\n");
                printf("ans = %d\n", ans);
                abort();
            }

            padic_clear(a2);
            padic_clear(d2);
        }
        else
        {
            result = (ans && padic_equal(a, d));
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

        padic_clear(a);
        padic_clear(b);
        padic_clear(c);
        padic_clear(d);

        fmpz_clear(p);
        padic_ctx_clear(ctx);
    }

/* PRIME p > 2 ***************************************************************/

    /* Check aliasing: a = sqrt(a) */
    for (i = 0; i < 10000; i++)
    {
        fmpz_t p;
        len_t N;
        padic_ctx_t ctx;

        int ans1, ans2;
        padic_t a, d;

        fmpz_init_set_ui(p, n_randtest_prime(state, 0));
        N = n_randint(state, PADIC_TEST_PREC_MAX - PADIC_TEST_PREC_MIN) 
            + PADIC_TEST_PREC_MIN;
        padic_ctx_init(ctx, p, FLINT_MAX(0, N-10), FLINT_MAX(0, N+10), PADIC_SERIES);

        padic_init2(a, N);
        padic_init2(d, N);

        padic_randtest(a, state, ctx);

        ans1 = padic_sqrt(d, a, ctx);
        ans2 = padic_sqrt(a, a, ctx);

        result = ((ans1 == ans2) && (!ans1 || padic_equal(a, d)));
        if (!result)
        {
            printf("FAIL (aliasing):\n\n");
            printf("a = "), padic_print(a, ctx), printf("\n");
            printf("d = "), padic_print(d, ctx), printf("\n");
            abort();
        }

        padic_clear(a);
        padic_clear(d);

        fmpz_clear(p);
        padic_ctx_clear(ctx);
    }

    /* Test random elements */
    for (i = 0; i < 10000; i++)
    {
        fmpz_t p;
        len_t N;
        padic_ctx_t ctx;

        int ans;
        padic_t a, b, d;

        fmpz_init_set_ui(p, n_randtest_prime(state, 0));
        N = n_randint(state, PADIC_TEST_PREC_MAX - PADIC_TEST_PREC_MIN) 
            + PADIC_TEST_PREC_MIN;
        padic_ctx_init(ctx, p, FLINT_MAX(0, N-10), FLINT_MAX(0, N+10), PADIC_SERIES);

        padic_init2(a, N);
        padic_init2(b, N);
        padic_init2(d, N);

        padic_randtest(a, state, ctx);

        ans = padic_sqrt(b, a, ctx);

        padic_mul(d, b, b, ctx);

        if (ans && padic_val(a) < 0)
        {
            len_t N2 = N + padic_val(a);
            padic_t a2, d2;

            padic_init2(a2, N2);
            padic_init2(d2, N2);

            padic_set(a2, a, ctx);
            padic_set(d2, d, ctx);

            result = (padic_equal(a2, d2));
            if (!result)
            {
                printf("FAIL (random elements):\n\n");
                printf("a  = "), padic_print(a, ctx),  printf("\n");
                printf("b  = "), padic_print(b, ctx),  printf("\n");
                printf("d  = "), padic_print(d, ctx),  printf("\n");
                printf("a2 = "), padic_print(a2, ctx), printf("\n");
                printf("d2 = "), padic_print(d2, ctx), printf("\n");
                printf("p = "), fmpz_print(p), printf("\n");
                printf("ans = %d\n", ans);
                abort();
            }

            padic_clear(a2);
            padic_clear(d2);
        }
        else
        {
            result = (!ans || padic_equal(a, d));
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

        padic_clear(a);
        padic_clear(b);
        padic_clear(d);

        fmpz_clear(p);
        padic_ctx_clear(ctx);
    }

    /* Test random squares */
    for (i = 0; i < 10000; i++)
    {
        fmpz_t p;
        len_t N;
        padic_ctx_t ctx;

        int ans;
        padic_t a, b, c, d;

        fmpz_init_set_ui(p, n_randtest_prime(state, 0));
        N = n_randint(state, PADIC_TEST_PREC_MAX - PADIC_TEST_PREC_MIN) 
            + PADIC_TEST_PREC_MIN;
        padic_ctx_init(ctx, p, FLINT_MAX(0, N-10), FLINT_MAX(0, N+10), PADIC_SERIES);

        padic_init2(a, N);
        padic_init2(b, N);
        padic_init2(c, N);
        padic_init2(d, N);

        padic_randtest(b, state, ctx);
        padic_mul(a, b, b, ctx);

        ans = padic_sqrt(c, a, ctx);

        padic_mul(d, c, c, ctx);

        if (ans && padic_val(a) < 0)
        {
            len_t N2 = N + padic_val(a);
            padic_t a2, d2;

            padic_init2(a2, N2);
            padic_init2(d2, N2);
            padic_set(a2, a, ctx);
            padic_set(d2, d, ctx);

            result = (padic_equal(a2, d2));
            if (!result)
            {
                printf("FAIL (random elements):\n\n");
                printf("a  = "), padic_print(a, ctx),  printf("\n");
                printf("b  = "), padic_print(b, ctx),  printf("\n");
                printf("c  = "), padic_print(c, ctx),  printf("\n");
                printf("d  = "), padic_print(d, ctx),  printf("\n");
                printf("a2 = "), padic_print(a2, ctx), printf("\n");
                printf("d2 = "), padic_print(d2, ctx), printf("\n");
                printf("p = "), fmpz_print(p), printf("\n");
                printf("ans = %d\n", ans);
                abort();
            }

            padic_clear(a2);
            padic_clear(d2);
        }
        else
        {
            result = (ans && padic_equal(a, d));
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

        padic_clear(a);
        padic_clear(b);
        padic_clear(c);
        padic_clear(d);

        fmpz_clear(p);
        padic_ctx_clear(ctx);
    }

    flint_randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}
