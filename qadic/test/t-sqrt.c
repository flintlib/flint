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

    Copyright (C) 2012 Sebastian Pancratz

******************************************************************************/

#include <stdio.h>
#include <stdlib.h>

#include "qadic.h"
#include "ulong_extras.h"
#include "long_extras.h"

extern int 
_artin_schreier_preimage(fmpz *rop, const fmpz *op, long len, 
                         const fmpz *a, const long *j, long lena);

int main(void)
{
    int i, result;
    flint_rand_t state;

    printf("sqrt... ");
    fflush(stdout);

    flint_randinit(state);

/* PRIME p = 2 ***************************************************************/

    /* Check Artin Schreier preimages */
    for (i = 0; i < 1000; i++)
    {
        fmpz_t p = {2L};
        long d;
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
                printf("FAIL (Artin Schreier preimages):\n\n");
                printf("a = "), qadic_print_pretty(a, ctx), printf("\n");
                printf("b = "), qadic_print_pretty(b, ctx), printf("\n");
                printf("c = "), qadic_print_pretty(c, ctx), printf("\n");
                qadic_ctx_print(ctx);
                abort();
            }
        }

        qadic_clear(a);
        qadic_clear(b);
        qadic_clear(c);

        qadic_ctx_clear(ctx);
    }

    /* Check aliasing: a = sqrt(a) */
    for (i = 0; i < 1000; i++)
    {
        fmpz_t p = {2L};
        long d, N;
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
            printf("FAIL (aliasing):\n\n");
            printf("a = "), qadic_print_pretty(a, ctx), printf("\n");
            printf("b = "), qadic_print_pretty(b, ctx), printf("\n");
            printf("c = "), qadic_print_pretty(c, ctx), printf("\n");
            printf("ans1,ans2 = %d,%d\n", ans1, ans2);
            qadic_ctx_print(ctx);
            abort();
        }

        qadic_clear(a);
        qadic_clear(b);

        qadic_ctx_clear(ctx);
    }

    /* Test random squares over finite fields */
    for (i = 0; i < 1000; i++)
    {
        fmpz_t p = {2L};
        long deg, N;
        qadic_ctx_t ctx;

        int ans;
        qadic_t a, b, c;

        deg = n_randint(state, 10) + 1;
        N = 1;
        qadic_ctx_init_conway(ctx, p, deg, FLINT_MAX(0, N-10), FLINT_MAX(0, N+10), "X", PADIC_SERIES);

        qadic_init2(a, N);
        qadic_init2(b, N);
        qadic_init2(c, N);

        qadic_randtest_val(b, state, 0, ctx);
        qadic_mul(a, b, b, ctx);

        ans = qadic_sqrt(c, a, ctx);
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
                printf("FAIL (a = b^2, c = sqrt(a), d = c^2 == a):\n\n");
                printf("a = "), qadic_print_pretty(a, ctx), printf("\n");
                printf("b = "), qadic_print_pretty(b, ctx), printf("\n");
                printf("c = "), qadic_print_pretty(c, ctx), printf("\n");
                printf("d = "), qadic_print_pretty(d, ctx), printf("\n");
                printf("e = "), qadic_print_pretty(e, ctx), printf("\n");
                printf("ans = %d\n", ans);
                printf("N = %ld\n", N);
                printf("N + val(a)/2 = %ld\n", N + qadic_val(a)/2);
                qadic_ctx_print(ctx);
                abort();
            }

            qadic_clear(d);
            qadic_clear(e);
        }
        else
        {
            printf("FAIL (a = b^2, c = sqrt(a), d = c^2 == a):\n\n");
            printf("a = "), qadic_print_pretty(a, ctx), printf("\n");
            printf("b = "), qadic_print_pretty(b, ctx), printf("\n");
            printf("c = "), qadic_print_pretty(c, ctx), printf("\n");
            printf("ans = %d\n", ans);
            qadic_ctx_print(ctx);
            abort();
        }

        qadic_clear(a);
        qadic_clear(b);
        qadic_clear(c);

        qadic_ctx_clear(ctx);
    }

    /* Test random elements over finite fields */
    for (i = 0; i < 1000; i++)
    {
        fmpz_t p = {2L};
        long d, N;
        qadic_ctx_t ctx;

        int ans;
        qadic_t a, b;

        d = n_randint(state, 10) + 1;
        N = 1;

        qadic_ctx_init_conway(ctx, p, d, FLINT_MAX(0,N-10), FLINT_MAX(0,N+10), "X", PADIC_SERIES);

        qadic_init2(a, N);
        qadic_init2(b, N);

        qadic_randtest_val(a, state, 0, ctx);

        ans = qadic_sqrt(b, a, ctx);
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
                printf("FAIL (random elements):\n\n");
                printf("a = "), qadic_print_pretty(a, ctx), printf("\n");
                printf("b = "), qadic_print_pretty(b, ctx), printf("\n");
                printf("c = "), qadic_print_pretty(c, ctx), printf("\n");
                printf("d = "), qadic_print_pretty(d, ctx), printf("\n");
                printf("ans = %d\n", ans);
                printf("N = %ld\n", N);
                printf("N + val(a)/2 = %ld\n", N + qadic_val(a)/2);
                qadic_ctx_print(ctx);
                abort();
            }

            qadic_clear(c);
            qadic_clear(d);
        }

        qadic_clear(a);
        qadic_clear(b);

        qadic_ctx_clear(ctx);
    }

    /* Test random squares in Zq */
    for (i = 0; i < 1000; i++)
    {
        fmpz_t p = {2L};
        long deg, N;
        qadic_ctx_t ctx;

        int ans;
        qadic_t a, b, c;

        deg = n_randint(state, 10) + 1;
        N = n_randint(state, 50) + 3;  /* N >= 3 */
        qadic_ctx_init_conway(ctx, p, deg, FLINT_MAX(0, N-10), FLINT_MAX(0, N+10), "X", PADIC_SERIES);
printf("i=%d\n", i), fflush(stdout);

        qadic_init2(a, N);
        qadic_init2(b, N);
        qadic_init2(c, N);

        qadic_randtest_int(b, state, ctx);
        qadic_mul(a, b, b, ctx);

        ans = qadic_sqrt(c, a, ctx);
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
                printf("FAIL (a = b^2, c = sqrt(a), d = c^2 == a):\n\n");
                printf("a = "), qadic_print_pretty(a, ctx), printf("\n");
                printf("b = "), qadic_print_pretty(b, ctx), printf("\n");
                printf("c = "), qadic_print_pretty(c, ctx), printf("\n");
                printf("d = "), qadic_print_pretty(d, ctx), printf("\n");
                printf("e = "), qadic_print_pretty(e, ctx), printf("\n");
                printf("ans = %d\n", ans);
                printf("N = %ld\n", N);
                printf("N + val(a)/2 = %ld\n", N + qadic_val(a)/2);
                qadic_ctx_print(ctx);
                abort();
            }

            qadic_clear(d);
            qadic_clear(e);
        }
        else
        {
            printf("FAIL (a = b^2, c = sqrt(a), d = c^2 == a):\n\n");
            printf("a = "), qadic_print_pretty(a, ctx), printf("\n");
            printf("b = "), qadic_print_pretty(b, ctx), printf("\n");
            printf("c = "), qadic_print_pretty(c, ctx), printf("\n");
            printf("ans = %d\n", ans);
            qadic_ctx_print(ctx);
            abort();
        }

        qadic_clear(a);
        qadic_clear(b);
        qadic_clear(c);

        fmpz_clear(p);
        qadic_ctx_clear(ctx);
    }

#if (0)

/* PRIME (any) ***************************************************************/

    /* Check aliasing: a = sqrt(a) */
    for (i = 0; i < 1000; i++)
    {
        fmpz_t p;
        long d, N;
        qadic_ctx_t ctx;

        int ans1, ans2;
        qadic_t a, b, c;

        fmpz_init_set_ui(p, n_randint(state, 2) ? 2 : n_randprime(state, 2 + n_randint(state, 3), 1));
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
            printf("FAIL (aliasing):\n\n");
            printf("a = "), qadic_print_pretty(a, ctx), printf("\n");
            printf("b = "), qadic_print_pretty(b, ctx), printf("\n");
            printf("c = "), qadic_print_pretty(c, ctx), printf("\n");
            printf("ans1,ans2 = %d,%d\n", ans1, ans2);
            qadic_ctx_print(ctx);
            abort();
        }

        qadic_clear(a);
        qadic_clear(b);

        fmpz_clear(p);
        qadic_ctx_clear(ctx);
    }

    /* Test random squares over finite fields */
    for (i = 0; i < 1000; i++)
    {
        fmpz_t p;
        long deg, N;
        qadic_ctx_t ctx;

        int ans;
        qadic_t a, b, c;

        fmpz_init_set_ui(p, n_randint(state, 2) ? 2 : n_randprime(state, 2 + n_randint(state, 3), 1));
        deg = n_randint(state, 10) + 1;
        N = 1;
        qadic_ctx_init_conway(ctx, p, deg, FLINT_MAX(0, N-10), FLINT_MAX(0, N+10), "X", PADIC_SERIES);

        qadic_init2(a, N);
        qadic_init2(b, N);
        qadic_init2(c, N);

        qadic_randtest_val(b, state, 0, ctx);
        qadic_mul(a, b, b, ctx);

        ans = qadic_sqrt(c, a, ctx);
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
                printf("FAIL (a = b^2, c = sqrt(a), d = c^2 == a):\n\n");
                printf("a = "), qadic_print_pretty(a, ctx), printf("\n");
                printf("b = "), qadic_print_pretty(b, ctx), printf("\n");
                printf("c = "), qadic_print_pretty(c, ctx), printf("\n");
                printf("d = "), qadic_print_pretty(d, ctx), printf("\n");
                printf("e = "), qadic_print_pretty(e, ctx), printf("\n");
                printf("ans = %d\n", ans);
                printf("N = %ld\n", N);
                printf("N + val(a)/2 = %ld\n", N + qadic_val(a)/2);
                qadic_ctx_print(ctx);
                abort();
            }

            qadic_clear(d);
            qadic_clear(e);
        }
        else
        {
            printf("FAIL (a = b^2, c = sqrt(a), d = c^2 == a):\n\n");
            printf("a = "), qadic_print_pretty(a, ctx), printf("\n");
            printf("b = "), qadic_print_pretty(b, ctx), printf("\n");
            printf("c = "), qadic_print_pretty(c, ctx), printf("\n");
            printf("ans = %d\n", ans);
            qadic_ctx_print(ctx);
            abort();
        }

        qadic_clear(a);
        qadic_clear(b);
        qadic_clear(c);

        fmpz_clear(p);
        qadic_ctx_clear(ctx);
    }

    /* Test random elements over finite fields */
    for (i = 0; i < 1000; i++)
    {
        fmpz_t p;
        long d, N;
        qadic_ctx_t ctx;

        int ans;
        qadic_t a, b;

        fmpz_init_set_ui(p, n_randint(state, 2) ? 2 : n_randprime(state, 2 + n_randint(state, 3), 1));
        d = n_randint(state, 10) + 1;
        N = 1;

        qadic_ctx_init_conway(ctx, p, d, FLINT_MAX(0,N-10), FLINT_MAX(0,N+10), "X", PADIC_SERIES);

        qadic_init2(a, N);
        qadic_init2(b, N);

        qadic_randtest_val(a, state, 0, ctx);

        ans = qadic_sqrt(b, a, ctx);
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
                printf("FAIL (random elements):\n\n");
                printf("a = "), qadic_print_pretty(a, ctx), printf("\n");
                printf("b = "), qadic_print_pretty(b, ctx), printf("\n");
                printf("c = "), qadic_print_pretty(c, ctx), printf("\n");
                printf("d = "), qadic_print_pretty(d, ctx), printf("\n");
                printf("ans = %d\n", ans);
                printf("N = %ld\n", N);
                printf("N + val(a)/2 = %ld\n", N + qadic_val(a)/2);
                qadic_ctx_print(ctx);
                abort();
            }

            qadic_clear(c);
            qadic_clear(d);
        }

        qadic_clear(a);
        qadic_clear(b);

        fmpz_clear(p);
        qadic_ctx_clear(ctx);
    }

    /* Test random squares */
    for (i = 0; i < 1000000; i++)
    {
        fmpz_t p;
        long deg, N;
        qadic_ctx_t ctx;

        int ans;
        qadic_t a, b, c;

        fmpz_init_set_ui(p, n_randint(state, 2) ? 2 : n_randprime(state, 2 + n_randint(state, 3), 1));
        deg = /* n_randint(state, 10) + 1; */ 3;
        N = /* z_randint(state, 50) + 1; */ 2;
        qadic_ctx_init_conway(ctx, p, deg, FLINT_MAX(0, N-10), FLINT_MAX(0, N+10), "X", PADIC_SERIES);
printf("i=%d\n", i), fflush(stdout);
        qadic_init2(a, N);
        qadic_init2(b, N);
        qadic_init2(c, N);

        qadic_randtest_val(b, state, 0, ctx); /* XXX */
        qadic_mul(a, b, b, ctx);

        ans = qadic_sqrt(c, a, ctx);
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
                printf("FAIL (a = b^2, c = sqrt(a), d = c^2 == a):\n\n");
                printf("a = "), qadic_print_pretty(a, ctx), printf("\n");
                printf("b = "), qadic_print_pretty(b, ctx), printf("\n");
                printf("c = "), qadic_print_pretty(c, ctx), printf("\n");
                printf("d = "), qadic_print_pretty(d, ctx), printf("\n");
                printf("e = "), qadic_print_pretty(e, ctx), printf("\n");
                printf("ans = %d\n", ans);
                printf("N = %ld\n", N);
                printf("N + val(a)/2 = %ld\n", N + qadic_val(a)/2);
                qadic_ctx_print(ctx);
                abort();
            }

            qadic_clear(d);
            qadic_clear(e);
        }
        else
        {
            printf("FAIL (a = b^2, c = sqrt(a), d = c^2 == a):\n\n");
            printf("a = "), qadic_print_pretty(a, ctx), printf("\n");
            printf("b = "), qadic_print_pretty(b, ctx), printf("\n");
            printf("c = "), qadic_print_pretty(c, ctx), printf("\n");
            printf("ans = %d\n", ans);
            qadic_ctx_print(ctx);
            abort();
        }

        qadic_clear(a);
        qadic_clear(b);
        qadic_clear(c);

        fmpz_clear(p);
        qadic_ctx_clear(ctx);
    }

    /* Test random elements */
    for (i = 0; i < 100000*0; i++)
    {
        fmpz_t p;
        long d, N;
        qadic_ctx_t ctx;

        int ans;
        qadic_t a, b;

        fmpz_init_set_ui(p, n_randint(state, 2) ? 2 : n_randprime(state, 2 + n_randint(state, 3), 1));
        d = n_randint(state, 10) + 1;
        N = z_randint(state, 50) + 1;

        qadic_ctx_init_conway(ctx, p, d, FLINT_MAX(0,N-10), FLINT_MAX(0,N+10), "X", PADIC_SERIES);

        qadic_init2(a, N);
        qadic_init2(b, N);

        qadic_randtest(a, state, ctx);

        ans = qadic_sqrt(b, a, ctx);
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
                printf("FAIL (random elements):\n\n");
                printf("a = "), qadic_print_pretty(a, ctx), printf("\n");
                printf("b = "), qadic_print_pretty(b, ctx), printf("\n");
                printf("c = "), qadic_print_pretty(c, ctx), printf("\n");
                printf("d = "), qadic_print_pretty(d, ctx), printf("\n");
                printf("ans = %d\n", ans);
                printf("N = %ld\n", N);
                printf("N + val(a)/2 = %ld\n", N + qadic_val(a)/2);
                qadic_ctx_print(ctx);
                abort();
            }

            qadic_clear(c);
            qadic_clear(d);
        }

        qadic_clear(a);
        qadic_clear(b);

        fmpz_clear(p);
        qadic_ctx_clear(ctx);
    }

#endif

    flint_randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}

