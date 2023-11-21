/*
    Copyright (C) 2012 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "ulong_extras.h"
#include "long_extras.h"
#include "qadic.h"

int
_artin_schreier_preimage(fmpz *rop, const fmpz *op, slong len,
                         const fmpz *a, const slong *j, slong lena);

TEST_FUNCTION_START(qadic_sqrt, state)
{
    int i, result;

/* PRIME p = 2 ***************************************************************/

    /* Check Artin Schreier preimages */
    for (i = 0; i < 200 * flint_test_multiplier(); i++)
    {
        fmpz_t p = {WORD(2)};
        slong d;
        qadic_ctx_t ctx;

        int ans;
        qadic_t a, b, c;

        d = n_randint(state, 10) + 1;
        qadic_ctx_init_conway(ctx, p, d, 1, 1, "X", PADIC_SERIES);

        qadic_init2(a, 1);
        qadic_init2(b, 1);
        qadic_init2(c, 1);

        qadic_randtest_val(a, state, 0, ctx);
        padic_poly_fit_length(b, d);

        ans = _artin_schreier_preimage(b->coeffs, a->coeffs, a->length,
                                       ctx->a, ctx->j, ctx->len);

        b->val = 0;
        _padic_poly_set_length(b, d);
        _padic_poly_normalise(b);

        if (ans)
        {
            qadic_mul(c, b, b, ctx);
            qadic_add(c, c, b, ctx);

            result = qadic_equal(a, c);

            if (!result)
            {
                flint_printf("FAIL (Artin Schreier preimages):\n\n");
                flint_printf("a = "), qadic_print_pretty(a, ctx), flint_printf("\n");
                flint_printf("b = "), qadic_print_pretty(b, ctx), flint_printf("\n");
                flint_printf("c = "), qadic_print_pretty(c, ctx), flint_printf("\n");
                qadic_ctx_print(ctx);
                fflush(stdout);
                flint_abort();
            }
        }

        qadic_clear(a);
        qadic_clear(b);
        qadic_clear(c);

        qadic_ctx_clear(ctx);
    }

    /* Check aliasing: a = sqrt(a) */
    for (i = 0; i < 200 * flint_test_multiplier(); i++)
    {
        fmpz_t p = {WORD(2)};
        slong d, N;
        qadic_ctx_t ctx;

        int ans1, ans2;
        qadic_t a, b, c;

        d = n_randint(state, 10) + 1;
        N = z_randint(state, 50) + 1;
        qadic_ctx_init_conway(ctx, p, d, FLINT_MAX(0,N-10), FLINT_MAX(0,N+10), "X", PADIC_SERIES);

        qadic_init2(a, N);
        qadic_init2(b, N);
        qadic_init2(c, N);

        qadic_randtest(a, state, ctx);
        qadic_set(c, a, ctx);

        ans1 = qadic_sqrt(b, a, ctx);
        ans2 = qadic_sqrt(a, a, ctx);

        result = ((ans1 == ans2) && (!ans1 || qadic_equal(a, b)));
        if (!result)
        {
            flint_printf("FAIL (aliasing):\n\n");
            flint_printf("a = "), qadic_print_pretty(a, ctx), flint_printf("\n");
            flint_printf("b = "), qadic_print_pretty(b, ctx), flint_printf("\n");
            flint_printf("c = "), qadic_print_pretty(c, ctx), flint_printf("\n");
            flint_printf("ans1,ans2 = %d,%d\n", ans1, ans2);
            qadic_ctx_print(ctx);
            fflush(stdout);
            flint_abort();
        }

        qadic_clear(a);
        qadic_clear(b);
        qadic_clear(c);

        qadic_ctx_clear(ctx);
    }

    /* Test random squares over finite fields */
    for (i = 0; i < 200 * flint_test_multiplier(); i++)
    {
        fmpz_t p = {WORD(2)};
        slong deg, N;
        qadic_ctx_t ctx;

        int ans, ans2;
        qadic_t a, b, c, c2;

        deg = n_randint(state, 10) + 1;
        N = 1;
        qadic_ctx_init_conway(ctx, p, deg, FLINT_MAX(0, N-10), FLINT_MAX(0, N+10), "X", PADIC_SERIES);

        qadic_init2(a, N);
        qadic_init2(b, N);
        qadic_init2(c, N);
        qadic_init2(c2, 0);

        qadic_randtest_val(b, state, 0, ctx);
        qadic_mul(a, b, b, ctx);

        ans = qadic_sqrt(c, a, ctx);
        ans2 = qadic_sqrt(c2, a, ctx);
        if (ans)
        {
            qadic_t d, e;
            qadic_init2(d, N + qadic_val(a)/2);
            qadic_init2(e, N + qadic_val(a)/2);

            qadic_mul(d, c, c, ctx);
            qadic_set(e, a, ctx);

            result = (qadic_equal(d, e));
            if (!result)
            {
                flint_printf("FAIL (a = b^2, c = sqrt(a), d = c^2 == a):\n\n");
                flint_printf("a = "), qadic_print_pretty(a, ctx), flint_printf("\n");
                flint_printf("b = "), qadic_print_pretty(b, ctx), flint_printf("\n");
                flint_printf("c = "), qadic_print_pretty(c, ctx), flint_printf("\n");
                flint_printf("d = "), qadic_print_pretty(d, ctx), flint_printf("\n");
                flint_printf("e = "), qadic_print_pretty(e, ctx), flint_printf("\n");
                flint_printf("ans = %d\n", ans);
                flint_printf("N = %wd\n", N);
                flint_printf("N + val(a)/2 = %wd\n", N + qadic_val(a)/2);
                qadic_ctx_print(ctx);
                fflush(stdout);
                flint_abort();
            }

            result = (ans == ans2);
            if (!result)
            {
                flint_printf("FAIL (zero output precision, random squares over finite field)\n\n");
                fflush(stdout);
                flint_abort();
            }

            qadic_clear(d);
            qadic_clear(e);
        }
        /* there is no reason for that to work */
        /*else
        {
            flint_printf("FAIL (a = b^2, c = sqrt(a), d = c^2 == a):\n\n");
            flint_printf("a = "), qadic_print_pretty(a, ctx), flint_printf("\n");
            flint_printf("b = "), qadic_print_pretty(b, ctx), flint_printf("\n");
            flint_printf("c = "), qadic_print_pretty(c, ctx), flint_printf("\n");
            flint_printf("ans = %d\n", ans);
            qadic_ctx_print(ctx);
            fflush(stdout);
            flint_abort();
            }*/

        qadic_clear(a);
        qadic_clear(b);
        qadic_clear(c);
        qadic_clear(c2);

        qadic_ctx_clear(ctx);
    }

    /* Test random elements over finite fields */
    for (i = 0; i < 200 * flint_test_multiplier(); i++)
    {
        fmpz_t p = {WORD(2)};
        slong d, N;
        qadic_ctx_t ctx;

        int ans, ans2;
        qadic_t a, b, b2;

        d = n_randint(state, 10) + 1;
        N = 1;

        qadic_ctx_init_conway(ctx, p, d, FLINT_MAX(0,N-10), FLINT_MAX(0,N+10), "X", PADIC_SERIES);

        qadic_init2(a, N);
        qadic_init2(b, N);
        qadic_init2(b2, 0);

        qadic_randtest_val(a, state, 0, ctx);

        ans = qadic_sqrt(b, a, ctx);
        ans2 = qadic_sqrt(b2, a, ctx);

        result = (ans == ans2);
        if (!result)
        {
            flint_printf("FAIL (zero output precision, random elements):\n\n");
            fflush(stdout);
            flint_abort();
        }

        if (ans)
        {
            qadic_t c, d;
            qadic_init2(c, N + qadic_val(a)/2);
            qadic_init2(d, N + qadic_val(a)/2);

            qadic_mul(c, b, b, ctx);
            qadic_set(d, a, ctx);

            result = (qadic_equal(c, d));
            if (!result)
            {
                flint_printf("FAIL (random elements):\n\n");
                flint_printf("a = "), qadic_print_pretty(a, ctx), flint_printf("\n");
                flint_printf("b = "), qadic_print_pretty(b, ctx), flint_printf("\n");
                flint_printf("c = "), qadic_print_pretty(c, ctx), flint_printf("\n");
                flint_printf("d = "), qadic_print_pretty(d, ctx), flint_printf("\n");
                flint_printf("ans = %d\n", ans);
                flint_printf("N = %wd\n", N);
                flint_printf("N + val(a)/2 = %wd\n", N + qadic_val(a)/2);
                qadic_ctx_print(ctx);
                fflush(stdout);
                flint_abort();
            }

            qadic_clear(c);
            qadic_clear(d);
        }

        qadic_clear(a);
        qadic_clear(b);
        qadic_clear(b2);

        qadic_ctx_clear(ctx);
    }

    /* Test random squares */
    for (i = 0; i < 200 * flint_test_multiplier(); i++)
    {
        fmpz_t p = {WORD(2)};
        slong deg, N, N2;
        qadic_ctx_t ctx;

        int ans, ans2;
        qadic_t a, b, c, c2;

        deg = n_randint(state, 10) + 1;
        /* N >= 3 */
        N = n_randint(state, 50) + 3;
        N2 = FLINT_MAX(n_randint(state, N), 3);
        qadic_ctx_init_conway(ctx, p, deg, FLINT_MAX(0, N-10), FLINT_MAX(0, N+10), "X", PADIC_SERIES);

        qadic_init2(a, N);
        qadic_init2(b, N);
        qadic_init2(c, N);
        qadic_init2(c2, N2);

        qadic_randtest_val(b, state, 0, ctx);
        qadic_randtest_int(b, state, ctx);
        qadic_mul(a, b, b, ctx);

        ans = qadic_sqrt(c, a, ctx);
        ans2 = qadic_sqrt(c2, a, ctx);

        result = (ans == ans2);
        if (!result)
        {
            flint_printf("FAIL (output prec lower than input prec, random squares):\n\n");
            fflush(stdout);
            flint_abort();
        }

        if (ans)
        {
            qadic_t d, e, u, v, w;
            qadic_init2(d, N + qadic_val(a)/2);
            qadic_init2(e, N + qadic_val(a)/2);
            qadic_init2(u, N);
            qadic_init2(v, N + qadic_val(a)/2);
            qadic_init2(w, N + qadic_val(a)/2);

            qadic_mul(d, c, c, ctx);
            qadic_set(e, a, ctx);

            result = (qadic_equal(d, e));
            if (!result)
            {
                flint_printf("FAIL (a = b^2, c = sqrt(a), d = c^2 == a):\n\n");
                flint_printf("a = "), qadic_print_pretty(a, ctx), flint_printf("\n");
                flint_printf("b = "), qadic_print_pretty(b, ctx), flint_printf("\n");
                flint_printf("c = "), qadic_print_pretty(c, ctx), flint_printf("\n");
                flint_printf("d = "), qadic_print_pretty(d, ctx), flint_printf("\n");
                flint_printf("e = "), qadic_print_pretty(e, ctx), flint_printf("\n");
                flint_printf("ans = %d\n", ans);
                flint_printf("N = %wd\n", N);
                flint_printf("N + val(a)/2 = %wd\n", N + qadic_val(a)/2);
                qadic_ctx_print(ctx);
                fflush(stdout);
                flint_abort();
            }

            qadic_clear(d);
            qadic_clear(e);
            qadic_clear(u);
            qadic_clear(v);
            qadic_clear(w);
        }
        else
        {
            if (N - qadic_val(a) >= 3)
            {
                flint_printf("FAIL (a = b^2, c = sqrt(a), d = c^2 == a):\n\n");
                flint_printf("a = "), qadic_print_pretty(a, ctx), flint_printf("\n");
                flint_printf("b = "), qadic_print_pretty(b, ctx), flint_printf("\n");
                flint_printf("c = "), qadic_print_pretty(c, ctx), flint_printf("\n");
                flint_printf("ans = %d\n", ans);
                flint_printf("N = %wd\n", N);
                flint_printf("N + val(a)/2 = %wd\n", N + qadic_val(a)/2);
                qadic_ctx_print(ctx);
                fflush(stdout);
                flint_abort();
            }
        }

        qadic_clear(a);
        qadic_clear(b);
        qadic_clear(c);
        qadic_clear(c2);

        qadic_ctx_clear(ctx);
    }

    /* Test random elements */
    for (i = 0; i < 200 * flint_test_multiplier(); i++)
    {
        fmpz_t p = {WORD(2)};
        slong d, N, N2;
        qadic_ctx_t ctx;

        int ans, ans2;
        qadic_t a, b, b2;

        d = n_randint(state, 10) + 1;
        N = z_randint(state, 50) + 1;
        N2 = N - n_randint(state, 10);

        qadic_ctx_init_conway(ctx, p, d, FLINT_MAX(0,N-10), FLINT_MAX(0,N+10), "X", PADIC_SERIES);

        qadic_init2(a, N);
        qadic_init2(b, N);
        qadic_init2(b2, N2);

        qadic_randtest(a, state, ctx);

        ans = qadic_sqrt(b, a, ctx);
        ans2 = qadic_sqrt(b2, a, ctx);

        result = (ans == ans2);
        if (!result)
        {
            flint_printf("FAIL (output prec lower than input prec, random elements):\n\n");
            fflush(stdout);
            flint_abort();
        }

        if (ans)
        {
            qadic_t c, d;
            qadic_init2(c, N + qadic_val(a)/2);
            qadic_init2(d, N + qadic_val(a)/2);

            qadic_mul(c, b, b, ctx);
            qadic_set(d, a, ctx);

            result = (qadic_equal(c, d));
            if (!result)
            {
                flint_printf("FAIL (random elements):\n\n");
                flint_printf("a = "), qadic_print_pretty(a, ctx), flint_printf("\n");
                flint_printf("b = "), qadic_print_pretty(b, ctx), flint_printf("\n");
                flint_printf("c = "), qadic_print_pretty(c, ctx), flint_printf("\n");
                flint_printf("d = "), qadic_print_pretty(d, ctx), flint_printf("\n");
                flint_printf("ans = %d\n", ans);
                flint_printf("N = %wd\n", N);
                flint_printf("N + val(a)/2 = %wd\n", N + qadic_val(a)/2);
                qadic_ctx_print(ctx);
                fflush(stdout);
                flint_abort();
            }

            qadic_clear(c);
            qadic_clear(d);
        }

        qadic_clear(a);
        qadic_clear(b);
        qadic_clear(b2);

        qadic_ctx_clear(ctx);
    }

/* PRIME p != 2 **************************************************************/

    /* Check aliasing: a = sqrt(a) */
    for (i = 0; i < 200 * flint_test_multiplier(); i++)
    {
        fmpz_t p;
        slong d, N, q;
        qadic_ctx_t ctx;

        int ans1, ans2;
        qadic_t a, b, c;

        q = 2;
        while (q == 2)
            q = n_randprime(state, 2 + n_randint(state, 3), 1);
        fmpz_init_set_ui(p, q);
        d = n_randint(state, 10) + 1;
        N = z_randint(state, 50) + 1;
        qadic_ctx_init_conway(ctx, p, d, FLINT_MAX(0,N-10), FLINT_MAX(0,N+10), "X", PADIC_SERIES);

        qadic_init2(a, N);
        qadic_init2(b, N);
        qadic_init2(c, N);

        qadic_randtest(a, state, ctx);
        qadic_set(c, a, ctx);

        ans1 = qadic_sqrt(b, a, ctx);
        ans2 = qadic_sqrt(a, a, ctx);

        result = ((ans1 == ans2) && (!ans1 || qadic_equal(a, b)));
        if (!result)
        {
            flint_printf("FAIL (aliasing):\n\n");
            flint_printf("a = "), qadic_print_pretty(a, ctx), flint_printf("\n");
            flint_printf("b = "), qadic_print_pretty(b, ctx), flint_printf("\n");
            flint_printf("c = "), qadic_print_pretty(c, ctx), flint_printf("\n");
            flint_printf("ans1,ans2 = %d,%d\n", ans1, ans2);
            qadic_ctx_print(ctx);
            fflush(stdout);
            flint_abort();
        }

        qadic_clear(a);
        qadic_clear(b);
        qadic_clear(c);

        fmpz_clear(p);
        qadic_ctx_clear(ctx);
    }

    /* Test random squares over finite fields */
    for (i = 0; i < 200 * flint_test_multiplier(); i++)
    {
        fmpz_t p;
        slong deg, N, q;
        qadic_ctx_t ctx;

        int ans, ans2;
        qadic_t a, b, c, c2;

        q = 2;
        while (q == 2)
            q = n_randprime(state, 2 + n_randint(state, 3), 1);
        fmpz_init_set_ui(p, q);
        deg = n_randint(state, 10) + 1;
        N = 1;
        qadic_ctx_init_conway(ctx, p, deg, FLINT_MAX(0, N-10), FLINT_MAX(0, N+10), "X", PADIC_SERIES);

        qadic_init2(a, N);
        qadic_init2(b, N);
        qadic_init2(c, N);
        qadic_init2(c2, 0);

        qadic_randtest_val(b, state, 0, ctx);
        qadic_mul(a, b, b, ctx);

        ans = qadic_sqrt(c, a, ctx);
        ans2 = qadic_sqrt(c2, a, ctx);

        result = (ans == ans2);
        if (!result)
        {
            flint_printf("FAIL (output prec lower than input prec, random squares over finite field):\n\n");
            fflush(stdout);
            flint_abort();
        }

        if (ans)
        {
            qadic_t d, e;
            qadic_init2(d, N + qadic_val(a)/2);
            qadic_init2(e, N + qadic_val(a)/2);

            qadic_mul(d, c, c, ctx);
            qadic_set(e, a, ctx);

            result = (qadic_equal(d, e));
            if (!result)
            {
                flint_printf("FAIL (a = b^2, c = sqrt(a), d = c^2 == a):\n\n");
                flint_printf("a = "), qadic_print_pretty(a, ctx), flint_printf("\n");
                flint_printf("b = "), qadic_print_pretty(b, ctx), flint_printf("\n");
                flint_printf("c = "), qadic_print_pretty(c, ctx), flint_printf("\n");
                flint_printf("d = "), qadic_print_pretty(d, ctx), flint_printf("\n");
                flint_printf("e = "), qadic_print_pretty(e, ctx), flint_printf("\n");
                flint_printf("ans = %d\n", ans);
                flint_printf("N = %wd\n", N);
                flint_printf("N + val(a)/2 = %wd\n", N + qadic_val(a)/2);
                qadic_ctx_print(ctx);
                fflush(stdout);
                flint_abort();
            }

            qadic_clear(d);
            qadic_clear(e);
        }
        else
        {
            flint_printf("FAIL (a = b^2, c = sqrt(a), d = c^2 == a):\n\n");
            flint_printf("a = "), qadic_print_pretty(a, ctx), flint_printf("\n");
            flint_printf("b = "), qadic_print_pretty(b, ctx), flint_printf("\n");
            flint_printf("c = "), qadic_print_pretty(c, ctx), flint_printf("\n");
            flint_printf("ans = %d\n", ans);
            qadic_ctx_print(ctx);
            fflush(stdout);
            flint_abort();
        }

        qadic_clear(a);
        qadic_clear(b);
        qadic_clear(c);
        qadic_clear(c2);

        fmpz_clear(p);
        qadic_ctx_clear(ctx);
    }

    /* Test random elements over finite fields */
    for (i = 0; i < 200 * flint_test_multiplier(); i++)
    {
        fmpz_t p;
        slong d, N, q;
        qadic_ctx_t ctx;

        int ans, ans2;
        qadic_t a, b, b2;

        q = 2;
        while (q == 2)
            q = n_randprime(state, 2 + n_randint(state, 3), 1);
        fmpz_init_set_ui(p, q);
        d = n_randint(state, 10) + 1;
        N = 1;

        qadic_ctx_init_conway(ctx, p, d, FLINT_MAX(0,N-10), FLINT_MAX(0,N+10), "X", PADIC_SERIES);

        qadic_init2(a, N);
        qadic_init2(b, N);
        qadic_init2(b2, 0);

        qadic_randtest_val(a, state, 0, ctx);

        ans = qadic_sqrt(b, a, ctx);
        ans2 = qadic_sqrt(b2, a, ctx);

        result = (ans == ans2);
        if (!result)
        {
            flint_printf("FAIL (output prec lower than input prec, random elements over finite field):\n\n");
            fflush(stdout);
            flint_abort();
        }

        if (ans)
        {
            qadic_t c, d;
            qadic_init2(c, N + qadic_val(a)/2);
            qadic_init2(d, N + qadic_val(a)/2);

            qadic_mul(c, b, b, ctx);
            qadic_set(d, a, ctx);

            result = (qadic_equal(c, d));
            if (!result)
            {
                flint_printf("FAIL (random elements):\n\n");
                flint_printf("a = "), qadic_print_pretty(a, ctx), flint_printf("\n");
                flint_printf("b = "), qadic_print_pretty(b, ctx), flint_printf("\n");
                flint_printf("c = "), qadic_print_pretty(c, ctx), flint_printf("\n");
                flint_printf("d = "), qadic_print_pretty(d, ctx), flint_printf("\n");
                flint_printf("ans = %d\n", ans);
                flint_printf("N = %wd\n", N);
                flint_printf("N + val(a)/2 = %wd\n", N + qadic_val(a)/2);
                qadic_ctx_print(ctx);
                fflush(stdout);
                flint_abort();
            }

            qadic_clear(c);
            qadic_clear(d);
        }

        qadic_clear(a);
        qadic_clear(b);
        qadic_clear(b2);

        fmpz_clear(p);
        qadic_ctx_clear(ctx);
    }

    /* Test random squares */
    for (i = 0; i < 200 * flint_test_multiplier(); i++)
    {
        fmpz_t p;
        slong deg, N, N2, q;
        qadic_ctx_t ctx;

        int ans, ans2;
        qadic_t a, b, c, c2;

        q = 2;
        while (q == 2)
            q = n_randprime(state, 2 + n_randint(state, 3), 1);
        fmpz_init_set_ui(p, q);
        deg = n_randint(state, 10) + 1;
        N = z_randint(state, 50) + 1;
        N2 = N - n_randint(state, 10);
        qadic_ctx_init_conway(ctx, p, deg, FLINT_MAX(0, N-10), FLINT_MAX(0, N+10), "X", PADIC_SERIES);
        qadic_init2(a, N);
        qadic_init2(b, N);
        qadic_init2(c, N);
        qadic_init2(c2, N2);

        qadic_randtest_val(b, state, 0, ctx); /* XXX */
        qadic_mul(a, b, b, ctx);

        ans = qadic_sqrt(c, a, ctx);
        ans2 = qadic_sqrt(c2, a, ctx);

        result = (ans == ans2);
        if (!result)
        {
            flint_printf("FAIL (output prec lower than input prec, random squares):\n\n");
            fflush(stdout);
            flint_abort();
        }

        if (ans)
        {
            qadic_t d, e;
            qadic_init2(d, N + qadic_val(a)/2);
            qadic_init2(e, N + qadic_val(a)/2);

            qadic_mul(d, c, c, ctx);
            qadic_set(e, a, ctx);

            result = (qadic_equal(d, e));
            if (!result)
            {
                flint_printf("FAIL (a = b^2, c = sqrt(a), d = c^2 == a):\n\n");
                flint_printf("a = "), qadic_print_pretty(a, ctx), flint_printf("\n");
                flint_printf("b = "), qadic_print_pretty(b, ctx), flint_printf("\n");
                flint_printf("c = "), qadic_print_pretty(c, ctx), flint_printf("\n");
                flint_printf("d = "), qadic_print_pretty(d, ctx), flint_printf("\n");
                flint_printf("e = "), qadic_print_pretty(e, ctx), flint_printf("\n");
                flint_printf("ans = %d\n", ans);
                flint_printf("N = %wd\n", N);
                flint_printf("N + val(a)/2 = %wd\n", N + qadic_val(a)/2);
                qadic_ctx_print(ctx);
                fflush(stdout);
                flint_abort();
            }

            qadic_clear(d);
            qadic_clear(e);
        }
        else
        {
            flint_printf("FAIL (a = b^2, c = sqrt(a), d = c^2 == a):\n\n");
            flint_printf("a = "), qadic_print_pretty(a, ctx), flint_printf("\n");
            flint_printf("b = "), qadic_print_pretty(b, ctx), flint_printf("\n");
            flint_printf("c = "), qadic_print_pretty(c, ctx), flint_printf("\n");
            flint_printf("ans = %d\n", ans);
            qadic_ctx_print(ctx);
            fflush(stdout);
            flint_abort();
        }

        qadic_clear(a);
        qadic_clear(b);
        qadic_clear(c);
        qadic_clear(c2);

        fmpz_clear(p);
        qadic_ctx_clear(ctx);
    }

    /* Test random elements */
    for (i = 0; i < 200 * flint_test_multiplier(); i++)
    {
        fmpz_t p;
        slong d, N, N2, q;
        qadic_ctx_t ctx;

        int ans, ans2;
        qadic_t a, b, b2;

        q = 2;
        while (q == 2)
            q = n_randprime(state, 2 + n_randint(state, 3), 1);
        fmpz_init_set_ui(p, q);
        d = n_randint(state, 10) + 1;
        N = z_randint(state, 50) + 1;
        N2 = N - n_randint(state, 10);

        qadic_ctx_init_conway(ctx, p, d, FLINT_MAX(0,N-10), FLINT_MAX(0,N+10), "X", PADIC_SERIES);

        qadic_init2(a, N);
        qadic_init2(b, N);
        qadic_init2(b2, N2);

        qadic_randtest(a, state, ctx);

        ans = qadic_sqrt(b, a, ctx);
        ans2 = qadic_sqrt(b2, a, ctx);

        result = (ans == ans2);
        if (!result)
        {
            flint_printf("FAIL (output prec lower than input prec, random elements):\n\n");
            fflush(stdout);
            flint_abort();
        }

        if (ans)
        {
            qadic_t c, d;
            qadic_init2(c, N + qadic_val(a)/2);
            qadic_init2(d, N + qadic_val(a)/2);

            qadic_mul(c, b, b, ctx);
            qadic_set(d, a, ctx);

            result = (qadic_equal(c, d));
            if (!result)
            {
                flint_printf("FAIL (random elements):\n\n");
                flint_printf("a = "), qadic_print_pretty(a, ctx), flint_printf("\n");
                flint_printf("b = "), qadic_print_pretty(b, ctx), flint_printf("\n");
                flint_printf("c = "), qadic_print_pretty(c, ctx), flint_printf("\n");
                flint_printf("d = "), qadic_print_pretty(d, ctx), flint_printf("\n");
                flint_printf("ans = %d\n", ans);
                flint_printf("N = %wd\n", N);
                flint_printf("N + val(a)/2 = %wd\n", N + qadic_val(a)/2);
                qadic_ctx_print(ctx);
                fflush(stdout);
                flint_abort();
            }

            qadic_clear(c);
            qadic_clear(d);
        }

        qadic_clear(a);
        qadic_clear(b);
        qadic_clear(b2);

        fmpz_clear(p);
        qadic_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
