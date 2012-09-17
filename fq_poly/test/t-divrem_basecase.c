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

#include "fq_poly.h"

#include "ulong_extras.h"
#include "long_extras.h"

int
main(void)
{
    int i, result;
    flint_rand_t state;

    printf("divrem_basecase... ");
    fflush(stdout);

    flint_randinit(state);

    /* Check q*b + r == a */
    for (i = 0; i < 2000; i++)
    {
        fmpz_t p;
        long d;
        fq_ctx_t ctx;

        fq_poly_t a, b, c, q, r;

        fmpz_init(p);
        fmpz_set_ui(p, n_randprime(state, 2 + n_randint(state, 3), 1));
        d   = n_randint(state, 10) + 1;
        fq_ctx_init_conway(ctx, p, d, "a");
        fq_poly_init(a);
        fq_poly_init(b);
        fq_poly_init(c);
        fq_poly_init(q);
        fq_poly_init(r);

        fq_poly_randtest(         a, state, n_randint(state, 50),     ctx);
        fq_poly_randtest_not_zero(b, state, n_randint(state, 50) + 1, ctx);

        fq_poly_divrem_basecase(q, r, a, b, ctx);
        fq_poly_mul(c, q, b, ctx);
        fq_poly_add(c, c, r, ctx);

        result = (fq_poly_equal(a, c));
        if (!result)
        {
            printf("FAIL:\n\n");
            printf("a = "), fq_poly_print_pretty(a, "X", ctx), printf("\n");
            printf("b = "), fq_poly_print_pretty(b, "X", ctx), printf("\n");
            printf("c = "), fq_poly_print_pretty(c, "X", ctx), printf("\n");
            printf("q = "), fq_poly_print_pretty(q, "X", ctx), printf("\n");
            printf("r = "), fq_poly_print_pretty(r, "X", ctx), printf("\n");
            abort();
        }

        fq_poly_clear(a);
        fq_poly_clear(b);
        fq_poly_clear(c);
        fq_poly_clear(q);
        fq_poly_clear(r);

        fmpz_clear(p);
        fq_ctx_clear(ctx);
    }

    /* Check aliasing: a and r */
    for (i = 0; i < 2000; i++)
    {
        fmpz_t p;
        long d;
        fq_ctx_t ctx;

        fq_poly_t a, b, q, r;

        fmpz_init(p);
        fmpz_set_ui(p, n_randprime(state, 2 + n_randint(state, 3), 1));
        d   = n_randint(state, 10) + 1;
        fq_ctx_init_conway(ctx, p, d, "a");
        fq_poly_init(a);
        fq_poly_init(b);
        fq_poly_init(q);
        fq_poly_init(r);

        fq_poly_randtest(         a, state, n_randint(state, 50),     ctx);
        fq_poly_randtest_not_zero(b, state, n_randint(state, 50) + 1, ctx);

        fq_poly_divrem_basecase(q, r, a, b, ctx);
        fq_poly_divrem_basecase(q, a, a, b, ctx);

        result = (fq_poly_equal(a, r));
        if (!result)
        {
            printf("FAIL:\n\n");
            printf("a = "), fq_poly_print_pretty(a, "X", ctx), printf("\n");
            printf("b = "), fq_poly_print_pretty(b, "X", ctx), printf("\n");
            printf("q = "), fq_poly_print_pretty(q, "X", ctx), printf("\n");
            printf("r = "), fq_poly_print_pretty(r, "X", ctx), printf("\n");
            abort();
        }

        fq_poly_clear(a);
        fq_poly_clear(b);
        fq_poly_clear(q);
        fq_poly_clear(r);

        fmpz_clear(p);
        fq_ctx_clear(ctx);
    }

    /* Check aliasing: b and r */
    for (i = 0; i < 2000; i++)
    {
        fmpz_t p;
        long d;
        fq_ctx_t ctx;

        fq_poly_t a, b, q, r;

        fmpz_init(p);
        fmpz_set_ui(p, n_randprime(state, 2 + n_randint(state, 3), 1));
        d   = n_randint(state, 10) + 1;
        fq_ctx_init_conway(ctx, p, d, "a");
        fq_poly_init(a);
        fq_poly_init(b);
        fq_poly_init(q);
        fq_poly_init(r);

        fq_poly_randtest(         a, state, n_randint(state, 50),     ctx);
        fq_poly_randtest_not_zero(b, state, n_randint(state, 50) + 1, ctx);

        fq_poly_divrem_basecase(q, r, a, b, ctx);
        fq_poly_divrem_basecase(q, b, a, b, ctx);

        result = (fq_poly_equal(b, r));
        if (!result)
        {
            printf("FAIL:\n\n");
            printf("a = "), fq_poly_print_pretty(a, "X", ctx), printf("\n");
            printf("b = "), fq_poly_print_pretty(b, "X", ctx), printf("\n");
            printf("q = "), fq_poly_print_pretty(q, "X", ctx), printf("\n");
            printf("r = "), fq_poly_print_pretty(r, "X", ctx), printf("\n");
            abort();
        }

        fq_poly_clear(a);
        fq_poly_clear(b);
        fq_poly_clear(q);
        fq_poly_clear(r);

        fmpz_clear(p);
        fq_ctx_clear(ctx);
    }

    /* Check aliasing: a and q */
    for (i = 0; i < 2000; i++)
    {
        fmpz_t p;
        long d;
        fq_ctx_t ctx;

        fq_poly_t a, b, q, r;

        fmpz_init(p);
        fmpz_set_ui(p, n_randprime(state, 2 + n_randint(state, 3), 1));
        d   = n_randint(state, 10) + 1;
        fq_ctx_init_conway(ctx, p, d, "a");
        fq_poly_init(a);
        fq_poly_init(b);
        fq_poly_init(q);
        fq_poly_init(r);

        fq_poly_randtest(         a, state, n_randint(state, 50),     ctx);
        fq_poly_randtest_not_zero(b, state, n_randint(state, 50) + 1, ctx);

        fq_poly_divrem_basecase(q, r, a, b, ctx);
        fq_poly_divrem_basecase(a, r, a, b, ctx);

        result = (fq_poly_equal(a, q));
        if (!result)
        {
            printf("FAIL:\n\n");
            printf("a = "), fq_poly_print_pretty(a, "X", ctx), printf("\n");
            printf("b = "), fq_poly_print_pretty(b, "X", ctx), printf("\n");
            printf("q = "), fq_poly_print_pretty(q, "X", ctx), printf("\n");
            printf("r = "), fq_poly_print_pretty(r, "X", ctx), printf("\n");
            abort();
        }

        fq_poly_clear(a);
        fq_poly_clear(b);
        fq_poly_clear(q);
        fq_poly_clear(r);

        fmpz_clear(p);
        fq_ctx_clear(ctx);
    }

    /* Check aliasing: b and q */
    for (i = 0; i < 2000; i++)
    {
        fmpz_t p;
        long d;
        fq_ctx_t ctx;

        fq_poly_t a, b, q, r;

        fmpz_init(p);
        fmpz_set_ui(p, n_randprime(state, 2 + n_randint(state, 3), 1));
        d   = n_randint(state, 10) + 1;
        fq_ctx_init_conway(ctx, p, d, "a");
        fq_poly_init(a);
        fq_poly_init(b);
        fq_poly_init(q);
        fq_poly_init(r);

        fq_poly_randtest(         a, state, n_randint(state, 50),     ctx);
        fq_poly_randtest_not_zero(b, state, n_randint(state, 50) + 1, ctx);

        fq_poly_divrem_basecase(q, r, a, b, ctx);
        fq_poly_divrem_basecase(b, r, a, b, ctx);

        result = (fq_poly_equal(b, q));
        if (!result)
        {
            printf("FAIL:\n\n");
            printf("a = "), fq_poly_print_pretty(a, "X", ctx), printf("\n");
            printf("b = "), fq_poly_print_pretty(b, "X", ctx), printf("\n");
            printf("q = "), fq_poly_print_pretty(q, "X", ctx), printf("\n");
            printf("r = "), fq_poly_print_pretty(r, "X", ctx), printf("\n");
            abort();
        }

        fq_poly_clear(a);
        fq_poly_clear(b);
        fq_poly_clear(q);
        fq_poly_clear(r);

        fmpz_clear(p);
        fq_ctx_clear(ctx);
    }

    flint_randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}

