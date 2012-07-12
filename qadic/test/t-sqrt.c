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
        qadic_ctx_init_conway(ctx, p, d, 1, "X", PADIC_SERIES);

        qadic_init(a);
        qadic_init(b);
        qadic_init(c);

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
        const fmpz_t p = {2L};
        long d, N;
        qadic_ctx_t ctx;

        int ans1, ans2;
        qadic_t a, b;

        d = n_randint(state, 10) + 1;
        N = z_randint(state, 50) + 1;
        qadic_ctx_init_conway(ctx, p, d, N, "X", PADIC_SERIES);

        qadic_init(a);
        qadic_init(b);

        qadic_randtest(a, state, ctx);

        ans1 = qadic_sqrt(b, a, ctx);
        ans2 = qadic_sqrt(a, a, ctx);

        result = ((ans1 == ans2) && (!ans1 || qadic_equal(a, b)));
        if (!result)
        {
            printf("FAIL (aliasing):\n\n");
            printf("a = "), qadic_print_pretty(a, ctx), printf("\n");
            printf("b = "), qadic_print_pretty(b, ctx), printf("\n");
            qadic_ctx_print(ctx);
            abort();
        }

        qadic_clear(a);
        qadic_clear(b);

        qadic_ctx_clear(ctx);
    }

/* PRIME p > 2 ***************************************************************/

    /* Check aliasing: a = sqrt(a) */
    for (i = 0; i < 1000; i++)
    {
        fmpz_t p;
        long d, N;
        qadic_ctx_t ctx;

        int ans1, ans2;
        qadic_t a, b;

        fmpz_init(p);
        fmpz_set_ui(p, n_randprime(state, 3 + n_randint(state, 3), 1));
        d = n_randint(state, 10) + 1;
        N = z_randint(state, 50) + 1;
        qadic_ctx_init_conway(ctx, p, d, N, "X", PADIC_SERIES);

        qadic_init(a);
        qadic_init(b);

        qadic_randtest(a, state, ctx);

        ans1 = qadic_sqrt(b, a, ctx);
        ans2 = qadic_sqrt(a, a, ctx);

        result = ((ans1 == ans2) && (!ans1 || qadic_equal(a, b)));
        if (!result)
        {
            printf("FAIL (aliasing):\n\n");
            printf("a = "), qadic_print_pretty(a, ctx), printf("\n");
            printf("b = "), qadic_print_pretty(b, ctx), printf("\n");
            qadic_ctx_print(ctx);
            abort();
        }

        qadic_clear(a);
        qadic_clear(b);

        fmpz_clear(p);
        qadic_ctx_clear(ctx);
    }

    /* Test random elements */
    for (i = 0; i < 1000; i++)
    {
        fmpz_t p;
        long d, N;
        qadic_ctx_t ctx;

        int ans;
        qadic_t a, b, c;

        fmpz_init(p);
        fmpz_set_ui(p, n_randprime(state, 3 + n_randint(state, 3), 1));
        d = n_randint(state, 10) + 1;
        N = z_randint(state, 50) + 1;
        qadic_ctx_init_conway(ctx, p, d, N, "X", PADIC_SERIES);

        qadic_init(a);
        qadic_init(b);
        qadic_init(c);

        qadic_randtest(a, state, ctx);

        ans = qadic_sqrt(b, a, ctx);

        qadic_mul(c, b, b, ctx);

        if (ans && a->val < 0)
        {
            qadic_t a2, c2;

            qadic_init(a2);
            qadic_init(c2);
            qadic_scalar_mod_ppow(a2, a, N + a->val, ctx);
            qadic_scalar_mod_ppow(c2, c, N + a->val, ctx);

            result = (qadic_equal(a2, c2));
            if (!result)
            {
                printf("FAIL (random elements):\n\n");
                printf("a  = "), qadic_print_pretty(a, ctx), printf("\n");
                printf("b  = "), qadic_print_pretty(b, ctx), printf("\n");
                printf("c  = "), qadic_print_pretty(c, ctx), printf("\n");
                printf("a2 = "), qadic_print_pretty(a2, ctx), printf("\n");
                printf("c2 = "), qadic_print_pretty(c2, ctx), printf("\n");
                printf("ans = %d\n", ans);
                qadic_ctx_print(ctx);
                abort();
            }

            qadic_clear(a2);
            qadic_clear(c2);
        }
        else
        {
            result = (!ans || qadic_equal(a, c));
            if (!result)
            {
                printf("FAIL (random elements):\n\n");
                printf("a = "), qadic_print_pretty(a, ctx), printf("\n");
                printf("b = "), qadic_print_pretty(b, ctx), printf("\n");
                printf("c = "), qadic_print_pretty(c, ctx), printf("\n");
                printf("ans = %d\n", ans);
                qadic_ctx_print(ctx);
                abort();
            }
        }

        qadic_clear(a);
        qadic_clear(b);
        qadic_clear(c);

        fmpz_clear(p);
        qadic_ctx_clear(ctx);
    }

    /* Test random squares */
    for (i = 0; i < 1000; i++)
    {
        fmpz_t p;
        long deg, N;
        qadic_ctx_t ctx;

        int ans;
        qadic_t a, b, c, d;

        fmpz_init(p);
        fmpz_set_ui(p, n_randprime(state, 3 + n_randint(state, 3), 1));
        deg = n_randint(state, 10) + 1;
        N = z_randint(state, 50) + 1;
        qadic_ctx_init_conway(ctx, p, deg, N, "X", PADIC_SERIES);

        qadic_init(a);
        qadic_init(b);
        qadic_init(c);
        qadic_init(d);

        qadic_randtest(b, state, ctx);
        qadic_mul(a, b, b, ctx);

        ans = qadic_sqrt(c, a, ctx);

        qadic_mul(d, c, c, ctx);

        if (ans && a->val < 0)
        {
            qadic_t a2, d2;

            qadic_init(a2);
            qadic_init(d2);
            qadic_scalar_mod_ppow(a2, a, N + a->val, ctx);
            qadic_scalar_mod_ppow(d2, d, N + a->val, ctx);

            result = (qadic_equal(a2, d2));
            if (!result)
            {
                printf("FAIL (a = b^2, c = sqrt(a), d = c^2 == a):\n\n");
                printf("a  = "), qadic_print_pretty(a, ctx), printf("\n");
                printf("b  = "), qadic_print_pretty(b, ctx), printf("\n");
                printf("c  = "), qadic_print_pretty(c, ctx), printf("\n");
                printf("d  = "), qadic_print_pretty(d, ctx), printf("\n");
                printf("a2 = "), qadic_print_pretty(a2, ctx), printf("\n");
                printf("d2 = "), qadic_print_pretty(d2, ctx), printf("\n");
                printf("p  = "), fmpz_print(p), printf("\n");
                printf("ans = %d\n", ans);
                qadic_ctx_print(ctx);
                abort();
            }

            qadic_clear(a2);
            qadic_clear(d2);
        }
        else
        {
            result = (ans && qadic_equal(a, d));
            if (!result)
            {
                printf("FAIL (a = b^2, c = sqrt(a), d = c^2 == a):\n\n");
                printf("a = "), qadic_print_pretty(a, ctx), printf("\n");
                printf("b = "), qadic_print_pretty(b, ctx), printf("\n");
                printf("c = "), qadic_print_pretty(c, ctx), printf("\n");
                printf("d = "), qadic_print_pretty(d, ctx), printf("\n");
                printf("ans = %d\n", ans);
                qadic_ctx_print(ctx);
                abort();
            }
        }

        qadic_clear(a);
        qadic_clear(b);
        qadic_clear(c);
        qadic_clear(d);

        fmpz_clear(p);
        qadic_ctx_clear(ctx);
    }

    flint_randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}

