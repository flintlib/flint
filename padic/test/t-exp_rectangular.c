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

int
main(void)
{
    int i, result;
    flint_rand_t state;

    printf("exp_rectangular... ");
    fflush(stdout);

    flint_randinit(state);

/** p == 2 *******************************************************************/

    /* Check aliasing: a = exp(a) */
    for (i = 0; i < 10000; i++)
    {
        fmpz_t p;
        len_t N;
        padic_ctx_t ctx;

        padic_t a, b;
        int ans1, ans2;

        fmpz_init_set_ui(p, 2);

        N = n_randint(state, PADIC_TEST_PREC_MAX - PADIC_TEST_PREC_MIN) 
            + PADIC_TEST_PREC_MIN;

        padic_ctx_init(ctx, p, FLINT_MAX(0, N-10), FLINT_MAX(0, N+10), PADIC_SERIES);

        padic_init2(a, N);
        padic_init2(b, N);

        padic_randtest(a, state, ctx);

        ans1 = padic_exp_rectangular(b, a, ctx);
        ans2 = padic_exp_rectangular(a, a, ctx);

        result = ((ans1 == ans2) && (!ans1 || padic_equal(a, b)));
        if (!result)
        {
            printf("FAIL (aliasing):\n\n");
            printf("a = "), padic_print(a, ctx), printf("\n");
            printf("b = "), padic_print(b, ctx), printf("\n");
            printf("ans1 = %d\n", ans1);
            printf("ans2 = %d\n", ans2);
            abort();
        }

        padic_clear(a);
        padic_clear(b);

        fmpz_clear(p);
        padic_ctx_clear(ctx);
    }

    /* Functional equation: exp(a + b) == exp(a) exp(b) */
    for (i = 0; i < 10000; i++)
    {
        fmpz_t p;
        len_t N;
        padic_ctx_t ctx;

        padic_t a, b, c, d, e, f, g;
        int ans1, ans2, ans3;

        fmpz_init_set_ui(p, 2);

        N = n_randint(state, PADIC_TEST_PREC_MAX - PADIC_TEST_PREC_MIN) 
            + PADIC_TEST_PREC_MIN;

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
        padic_add(c, a, b, ctx);

        ans1 = padic_exp_rectangular(d, a, ctx);
        ans2 = padic_exp_rectangular(e, b, ctx);
        padic_mul(f, d, e, ctx);

        ans3 = padic_exp_rectangular(g, c, ctx);

        result = (!ans1 || !ans2 || (ans3 && padic_equal(f, g)));
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

/** p > 2 ********************************************************************/

    /* Check aliasing: a = exp(a) */
    for (i = 0; i < 10000; i++)
    {
        fmpz_t p;
        len_t N;
        padic_ctx_t ctx;

        padic_t a, b;
        int ans1, ans2;

        fmpz_init_set_ui(p, n_randtest_prime(state, 0));

        N = n_randint(state, PADIC_TEST_PREC_MAX - PADIC_TEST_PREC_MIN) 
            + PADIC_TEST_PREC_MIN;

        padic_ctx_init(ctx, p, FLINT_MAX(0, N-10), FLINT_MAX(0, N+10), PADIC_SERIES);

        padic_init2(a, N);
        padic_init2(b, N);

        padic_randtest(a, state, ctx);

        ans1 = padic_exp_rectangular(b, a, ctx);
        ans2 = padic_exp_rectangular(a, a, ctx);

        result = ((ans1 == ans2) && (!ans1 || padic_equal(a, b)));
        if (!result)
        {
            printf("FAIL (aliasing):\n\n");
            printf("a = "), padic_print(a, ctx), printf("\n");
            printf("b = "), padic_print(b, ctx), printf("\n");
            printf("ans1 = %d\n", ans1);
            printf("ans2 = %d\n", ans2);
            abort();
        }

        padic_clear(a);
        padic_clear(b);

        fmpz_clear(p);
        padic_ctx_clear(ctx);
    }

    /* Functional equation: exp(a + b) == exp(a) exp(b) */
    for (i = 0; i < 10000; i++)
    {
        fmpz_t p;
        len_t N;
        padic_ctx_t ctx;

        padic_t a, b, c, d, e, f, g;
        int ans1, ans2, ans3;

        fmpz_init_set_ui(p, n_randtest_prime(state, 0));

        N = n_randint(state, PADIC_TEST_PREC_MAX - PADIC_TEST_PREC_MIN) 
            + PADIC_TEST_PREC_MIN;

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
        padic_add(c, a, b, ctx);

        ans1 = padic_exp_rectangular(d, a, ctx);
        ans2 = padic_exp_rectangular(e, b, ctx);
        padic_mul(f, d, e, ctx);

        ans3 = padic_exp_rectangular(g, c, ctx);

        result = (!ans1 || !ans2 || (ans3 && padic_equal(f, g)));
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

    flint_randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}

