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
#include "padic.h"

/*
    Set-up.  Currently we only test the logarithm for positive values of N. 
    This is important as for negative N, exp(0) is 1, which is 0 mod p^N, 
    and then log(0) does not converge.
 */
static len_t __rand_prec(flint_rand_t state, len_t i)
{
    len_t N;

    N = n_randint(state, PADIC_TEST_PREC_MAX) + 1;

    return N;
}

int
main(void)
{
    int i, result;
    flint_rand_t state;

    printf("log... ");
    fflush(stdout);

    flint_randinit(state);

/** p == 2 *******************************************************************/

    /* Check aliasing: a = log(a) */
    for (i = 0; i < 1000; i++)
    {
        fmpz_t p = {2L};
        len_t N;
        padic_ctx_t ctx;

        padic_t a, b;
        int ans1, ans2;

        N = __rand_prec(state, i);
        padic_ctx_init(ctx, p, FLINT_MAX(0, N-10), FLINT_MAX(0, N+10), PADIC_SERIES);

        padic_init2(a, N);
        padic_init2(b, N);

        padic_randtest(a, state, ctx);

        padic_one(b);
        padic_add(a, a, b, ctx);

        ans1 = padic_log(b, a, ctx);
        ans2 = padic_log(a, a, ctx);

        result = (ans1 == ans2) && (!ans1 || padic_equal(a, b));
        if (!result)
        {
            printf("FAIL (aliasing):\n\n");
            printf("a = "), padic_print(a, ctx), printf("\n");
            printf("b = "), padic_print(b, ctx), printf("\n");
            abort();
        }

        padic_clear(a);
        padic_clear(b);

        padic_ctx_clear(ctx);
    }

    /* Check: log(a) + log(b) == log(a * b) */
    for (i = 0; i < 10000; i++)
    {
        fmpz_t p = {2L};
        len_t N;
        padic_ctx_t ctx;

        padic_t a, b, c, d, e, f, g;
        int ans1, ans2, ans3;

        N = __rand_prec(state, i);
        padic_ctx_init(ctx, p, FLINT_MAX(0, N-10), FLINT_MAX(0, N+10), PADIC_SERIES);

        padic_init2(a, N);
        padic_init2(b, N);
        padic_init2(c, N);
        padic_init2(d, N);
        padic_init2(e, N);
        padic_init2(f, N);
        padic_init2(g, N);

        padic_randtest(a, state, ctx);
        padic_randtest(b, state, ctx);

        padic_one(c);
        padic_add(a, a, c, ctx);
        padic_add(b, b, c, ctx);

        padic_mul(c, a, b, ctx);

        ans1 = padic_log(d, a, ctx);
        ans2 = padic_log(e, b, ctx);
        padic_add(f, d, e, ctx);

        ans3 = padic_log(g, c, ctx);

        result = (!ans1 || !ans2 || (ans3 && padic_equal(f, g)));
        if (!result)
        {
            printf("FAIL (functional equation):\n\n");
            printf("a                   = "), padic_print(a, ctx), printf("\n");
            printf("b                   = "), padic_print(b, ctx), printf("\n");
            printf("c = a * b           = "), padic_print(c, ctx), printf("\n");
            printf("d = log(a)          = "), padic_print(d, ctx), printf("\n");
            printf("e = log(b)          = "), padic_print(e, ctx), printf("\n");
            printf("f = log(a) + log(b) = "), padic_print(f, ctx), printf("\n");
            printf("g = log(a * b)      = "), padic_print(g, ctx), printf("\n");
            abort();
        }

        padic_clear(a);
        padic_clear(b);
        padic_clear(c);
        padic_clear(d);
        padic_clear(e);
        padic_clear(f);
        padic_clear(g);

        padic_ctx_clear(ctx);
    }

    /* Check: log(exp(x)) == x */
    for (i = 0; i < 10000; i++)
    {
        fmpz_t p = {2L};
        len_t N;
        padic_ctx_t ctx;

        padic_t a, b, c;
        int ans1, ans2;

        N = __rand_prec(state, i);
        padic_ctx_init(ctx, p, FLINT_MAX(0, N-10), FLINT_MAX(0, N+10), PADIC_SERIES);

        padic_init2(a, N);
        padic_init2(b, N);
        padic_init2(c, N);

        padic_randtest(a, state, ctx);

        ans1 = padic_exp(b, a, ctx);
        if (ans1)
            ans2 = padic_log(c, b, ctx);

        result = !ans1 || (ans1 == ans2 && padic_equal(a, c));
        if (!result)
        {
            printf("FAIL (log(exp(x)) == x):\n\n");
            printf("a = "), padic_print(a, ctx), printf("\n");
            printf("b = "), padic_print(b, ctx), printf("\n");
            printf("c = "), padic_print(c, ctx), printf("\n");
            printf("ans1 = %d\n", ans1);
            printf("ans2 = %d\n", ans2);
            abort();
        }

        padic_clear(a);
        padic_clear(b);
        padic_clear(c);

        padic_ctx_clear(ctx);
    }

/** p > 2 ********************************************************************/

    /* Check aliasing: a = log(a) */
    for (i = 0; i < 1000; i++)
    {
        fmpz_t p;
        len_t N;
        padic_ctx_t ctx;

        padic_t a, b;
        int ans1, ans2;

        fmpz_init_set_ui(p, n_randtest_prime(state, 0));
        N = __rand_prec(state, i);
        padic_ctx_init(ctx, p, FLINT_MAX(0, N-10), FLINT_MAX(0, N+10), PADIC_SERIES);

        padic_init2(a, N);
        padic_init2(b, N);

        padic_randtest(a, state, ctx);

        padic_one(b);
        padic_add(a, a, b, ctx);

        ans1 = padic_log(b, a, ctx);
        ans2 = padic_log(a, a, ctx);

        result = (ans1 == ans2) && (!ans1 || padic_equal(a, b));
        if (!result)
        {
            printf("FAIL (aliasing):\n\n");
            printf("a = "), padic_print(a, ctx), printf("\n");
            printf("b = "), padic_print(b, ctx), printf("\n");
            abort();
        }

        padic_clear(a);
        padic_clear(b);

        fmpz_clear(p);
        padic_ctx_clear(ctx);
    }

    /* Check: log(a) + log(b) == log(a * b) */
    for (i = 0; i < 10000; i++)
    {
        fmpz_t p;
        len_t N;
        padic_ctx_t ctx;

        padic_t a, b, c, d, e, f, g;
        int ans1, ans2, ans3;

/*      fmpz_init_set_ui(p, n_randtest_prime(state, 0)); */
        fmpz_init_set_ui(p, n_randprime(state, 5, 1));
        N = __rand_prec(state, i);
        padic_ctx_init(ctx, p, FLINT_MAX(0, N-10), FLINT_MAX(0, N+10), PADIC_SERIES);

        padic_init2(a, N);
        padic_init2(b, N);
        padic_init2(c, N);
        padic_init2(d, N);
        padic_init2(e, N);
        padic_init2(f, N);
        padic_init2(g, N);

        padic_randtest(a, state, ctx);
        padic_randtest(b, state, ctx);

        padic_one(c);
        padic_add(a, a, c, ctx);

        padic_one(c);
        padic_add(b, b, c, ctx);

        padic_mul(c, a, b, ctx);

        ans1 = padic_log(d, a, ctx);
        ans2 = padic_log(e, b, ctx);
        padic_add(f, d, e, ctx);

        ans3 = padic_log(g, c, ctx);

        result = (!ans1 || !ans2 || (ans3 && padic_equal(f, g)));
        if (!result)
        {
            printf("FAIL (functional equation):\n\n");
            printf("a                   = "), padic_print(a, ctx), printf("\n");
            printf("b                   = "), padic_print(b, ctx), printf("\n");
            printf("c = a * b           = "), padic_print(c, ctx), printf("\n");
            printf("d = log(a)          = "), padic_print(d, ctx), printf("\n");
            printf("e = log(b)          = "), padic_print(e, ctx), printf("\n");
            printf("f = log(a) + log(b) = "), padic_print(f, ctx), printf("\n");
            printf("g = log(a * b)      = "), padic_print(g, ctx), printf("\n");
            printf("ans1 = %d\n", ans1);
            printf("ans2 = %d\n", ans2);
            printf("ans3 = %d\n", ans3);
            abort();
        }

        padic_clear(a);
        padic_clear(b);
        padic_clear(c);
        padic_clear(d);
        padic_clear(e);
        padic_clear(f);
        padic_clear(g);

        fmpz_clear(p);
        padic_ctx_clear(ctx);
    }

    /* Check: log(exp(x)) == x */
    for (i = 0; i < 10000; i++)
    {
        fmpz_t p;
        len_t N;
        padic_ctx_t ctx;

        padic_t a, b, c;
        int ans1, ans2;

        fmpz_init_set_ui(p, n_randtest_prime(state, 0));
        N = __rand_prec(state, i);
        padic_ctx_init(ctx, p, FLINT_MAX(0, N-10), FLINT_MAX(0, N+10), PADIC_SERIES);

        padic_init2(a, N);
        padic_init2(b, N);
        padic_init2(c, N);

        padic_randtest(a, state, ctx);

        ans1 = padic_exp(b, a, ctx);
        if (ans1)
            ans2 = padic_log(c, b, ctx);

        result = !ans1 || (ans1 == ans2 && padic_equal(a, c));
        if (!result)
        {
            printf("FAIL (log(exp(x)) == x):\n\n");
            printf("a = "), padic_print(a, ctx), printf("\n");
            printf("b = "), padic_print(b, ctx), printf("\n");
            printf("c = "), padic_print(c, ctx), printf("\n");
            printf("ans1 = %d\n", ans1);
            printf("ans2 = %d\n", ans2);
            abort();
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

