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

    printf("add... ");
    fflush(stdout);

    flint_randinit(state);

    /* Check aliasing: a = a + b */
    for (i = 0; i < 10000; i++)
    {
        fmpz_t p;
        len_t N;
        padic_ctx_t ctx;

        padic_t a, b, d;

        fmpz_init_set_ui(p, n_randtest_prime(state, 0));
        N = n_randint(state, PADIC_TEST_PREC_MAX - PADIC_TEST_PREC_MIN) 
            + PADIC_TEST_PREC_MIN;
        padic_ctx_init(ctx, p, FLINT_MAX(0, N-10), FLINT_MAX(0, N+10), PADIC_SERIES);

        padic_init2(a, N);
        padic_init2(b, N);
        padic_init2(d, N);

        padic_randtest(a, state, ctx);
        padic_randtest(b, state, ctx);

        padic_add(d, a, b, ctx);
        padic_add(a, a, b, ctx);

        result = (padic_equal(a, d));
        if (!result)
        {
            printf("FAIL:\n\n");
            printf("a = "), padic_print(a, ctx), printf("\n");
            printf("b = "), padic_print(b, ctx), printf("\n");
            printf("d = "), padic_print(d, ctx), printf("\n");
            abort();
        }

        padic_clear(a);
        padic_clear(b);
        padic_clear(d);

        fmpz_clear(p);
        padic_ctx_clear(ctx);
    }

    /* Check aliasing: b = a + b */
    for (i = 0; i < 10000; i++)
    {
        fmpz_t p;
        len_t N;
        padic_ctx_t ctx;

        padic_t a, b, d;

        fmpz_init_set_ui(p, n_randtest_prime(state, 0));
        N = n_randint(state, PADIC_TEST_PREC_MAX - PADIC_TEST_PREC_MIN) 
            + PADIC_TEST_PREC_MIN;
        padic_ctx_init(ctx, p, FLINT_MAX(0, N-10), FLINT_MAX(0, N+10), PADIC_SERIES);

        padic_init2(a, N);
        padic_init2(b, N);
        padic_init2(d, N);

        padic_randtest(a, state, ctx);
        padic_randtest(b, state, ctx);

        padic_add(d, a, b, ctx);
        padic_add(b, a, b, ctx);

        result = (padic_equal(b, d));
        if (!result)
        {
            printf("FAIL:\n\n");
            printf("a = "), padic_print(a, ctx), printf("\n");
            printf("b = "), padic_print(b, ctx), printf("\n");
            printf("d = "), padic_print(d, ctx), printf("\n");
            abort();
        }

        padic_clear(a);
        padic_clear(b);
        padic_clear(d);

        fmpz_clear(p);
        padic_ctx_clear(ctx);
    }

    /* Check aliasing: a = a + a */
    for (i = 0; i < 10000; i++)
    {
        fmpz_t p;
        len_t N;
        padic_ctx_t ctx;

        padic_t a, d;

        fmpz_init_set_ui(p, n_randtest_prime(state, 0));
        N = n_randint(state, PADIC_TEST_PREC_MAX - PADIC_TEST_PREC_MIN) 
            + PADIC_TEST_PREC_MIN;
        padic_ctx_init(ctx, p, FLINT_MAX(0, N-10), FLINT_MAX(0, N+10), PADIC_SERIES);

        padic_init2(a, N);
        padic_init2(d, N);

        padic_randtest(a, state, ctx);

        padic_add(d, a, a, ctx);
        padic_add(a, a, a, ctx);

        result = (padic_equal(a, d));
        if (!result)
        {
            printf("FAIL:\n\n");
            printf("a = "), padic_print(a, ctx), printf("\n");
            printf("d = "), padic_print(d, ctx), printf("\n");
            abort();
        }

        padic_clear(a);
        padic_clear(d);

        fmpz_clear(p);
        padic_ctx_clear(ctx);
    }

    /* Check that a + b == b + a */
    for (i = 0; i < 10000; i++)
    {
        fmpz_t p;
        len_t N;
        padic_ctx_t ctx;

        padic_t a, b, c, d;

        fmpz_init_set_ui(p, n_randtest_prime(state, 0));
        N = n_randint(state, PADIC_TEST_PREC_MAX - PADIC_TEST_PREC_MIN) 
            + PADIC_TEST_PREC_MIN;
        padic_ctx_init(ctx, p, FLINT_MAX(0, N-10), FLINT_MAX(0, N+10), PADIC_SERIES);

        padic_init2(a, N);
        padic_init2(b, N);
        padic_init2(c, N);
        padic_init2(d, N);

        padic_randtest(a, state, ctx);
        padic_randtest(b, state, ctx);

        padic_add(c, a, b, ctx);
        padic_add(d, b, a, ctx);

        result = (padic_equal(c, d));
        if (!result)
        {
            printf("FAIL:\n\n");
            printf("a = "), padic_print(a, ctx), printf("\n");
            printf("b = "), padic_print(b, ctx), printf("\n");
            printf("c = "), padic_print(c, ctx), printf("\n");
            printf("d = "), padic_print(d, ctx), printf("\n");
            abort();
        }

        padic_clear(a);
        padic_clear(b);
        padic_clear(c);
        padic_clear(d);

        fmpz_clear(p);
        padic_ctx_clear(ctx);
    }

    /* Check that (a + b) + c == a + (b + c) */
    for (i = 0; i < 10000; i++)
    {
        fmpz_t p;
        len_t N;
        padic_ctx_t ctx;

        padic_t a, b, c, d, e;

        fmpz_init_set_ui(p, n_randtest_prime(state, 0));
        N = n_randint(state, PADIC_TEST_PREC_MAX - PADIC_TEST_PREC_MIN) 
            + PADIC_TEST_PREC_MIN;
        padic_ctx_init(ctx, p, FLINT_MAX(0, N-10), FLINT_MAX(0, N+10), PADIC_SERIES);

        padic_init2(a, N);
        padic_init2(b, N);
        padic_init2(c, N);
        padic_init2(d, N);
        padic_init2(e, N);

        padic_randtest(a, state, ctx);
        padic_randtest(b, state, ctx);
        padic_randtest(c, state, ctx);

        padic_add(d, a, b, ctx);
        padic_add(d, d, c, ctx);

        padic_add(e, b, c, ctx);
        padic_add(e, a, e, ctx);

        result = (padic_equal(d, e));
        if (!result)
        {
            printf("FAIL:\n\n");
            printf("a = "), padic_print(a, ctx), printf("\n");
            printf("b = "), padic_print(b, ctx), printf("\n");
            printf("c = "), padic_print(c, ctx), printf("\n");
            printf("d = "), padic_print(d, ctx), printf("\n");
            printf("e = "), padic_print(e, ctx), printf("\n");
            abort();
        }

        padic_clear(a);
        padic_clear(b);
        padic_clear(c);
        padic_clear(d);
        padic_clear(e);

        fmpz_clear(p);
        padic_ctx_clear(ctx);
    }

    /* Check that a + 0 == a */
    for (i = 0; i < 10000; i++)
    {
        fmpz_t p;
        len_t N;
        padic_ctx_t ctx;

        padic_t a, b;

        fmpz_init_set_ui(p, n_randtest_prime(state, 0));
        N = n_randint(state, PADIC_TEST_PREC_MAX - PADIC_TEST_PREC_MIN) 
            + PADIC_TEST_PREC_MIN;
        padic_ctx_init(ctx, p, FLINT_MAX(0, N-10), FLINT_MAX(0, N+10), PADIC_SERIES);

        padic_init2(a, N);
        padic_init2(b, N);

        padic_randtest(a, state, ctx);
        padic_zero(b);

        padic_add(b, a, b, ctx);

        result = (padic_equal(a, b));
        if (!result)
        {
            printf("FAIL:\n\n");
            printf("a = "), padic_print(a, ctx), printf("\n");
            printf("b = "), padic_print(b, ctx), printf("\n");
            abort();
        }

        padic_clear(a);
        padic_clear(b);

        fmpz_clear(p);
        padic_ctx_clear(ctx);
    }

    flint_randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}

