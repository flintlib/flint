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

#include <stdio.h>
#include <stdlib.h>
#include <mpir.h>
#include "flint.h"
#include "ulong_extras.h"
#include "padic.h"

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
    for (i = 0; i < 10000; i++)
    {
        fmpz_t p = {2L};
        long N;
        padic_ctx_t ctx;

        padic_t a, b;

        N = n_randint(state, 50) + 1;
        padic_ctx_init(ctx, p, N, PADIC_SERIES);

        _padic_init(a);
        _padic_init(b);

        padic_randtest(a, state, ctx);

        if (!_padic_is_zero(a) && padic_val(a) < 2)
        {
            padic_val(a) = 2;
            padic_reduce(a, ctx);
        }

        _padic_one(b);
        _padic_add(a, a, b, ctx);

        padic_log(b, a, ctx);
        padic_log(a, a, ctx);

        result = (padic_equal(a, b, ctx));
        if (!result)
        {
            printf("FAIL (aliasing):\n\n");
            printf("a = "), padic_print(a, ctx), printf("\n");
            printf("b = "), padic_print(b, ctx), printf("\n");
            abort();
        }

        _padic_clear(a);
        _padic_clear(b);

        padic_ctx_clear(ctx);
    }

    /* Check: log(a) + log(b) == log(a * b) */
    for (i = 0; i < 10000; i++)
    {
        fmpz_t p = {2L};
        long N;
        padic_ctx_t ctx;

        padic_t a, b, c, d, e, f, g;

        N = n_randint(state, 50) + 1;
        padic_ctx_init(ctx, p, N, PADIC_SERIES);

        _padic_init(a);
        _padic_init(b);
        _padic_init(c);
        _padic_init(d);
        _padic_init(e);
        _padic_init(f);
        _padic_init(g);

        padic_randtest(a, state, ctx);
        padic_randtest(b, state, ctx);

        if (!_padic_is_zero(a) && padic_val(a) < 2)
        {
            padic_val(a) = 2;
            padic_reduce(a, ctx);
        }
        if (!_padic_is_zero(b) && padic_val(b) < 2)
        {
            padic_val(b) = 2;
            padic_reduce(b, ctx);
        }

        _padic_one(c);
        _padic_add(a, a, c, ctx);
        _padic_add(b, b, c, ctx);

        padic_mul(c, a, b, ctx);

        padic_log(d, a, ctx);
        padic_log(e, b, ctx);
        padic_add(f, d, e, ctx);

        padic_log(g, c, ctx);

        result = (padic_equal(f, g, ctx));
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

        _padic_clear(a);
        _padic_clear(b);
        _padic_clear(c);
        _padic_clear(d);
        _padic_clear(e);
        _padic_clear(f);
        _padic_clear(g);

        padic_ctx_clear(ctx);
    }

    /* Check: log(exp(x)) == x */
    for (i = 0; i < 10000; i++)
    {
        fmpz_t p = {2L};
        long N;
        padic_ctx_t ctx;

        padic_t a, b, c;

        N = n_randint(state, 50) + 1;
        padic_ctx_init(ctx, p, N, PADIC_SERIES);

        _padic_init(a);
        _padic_init(b);
        _padic_init(c);

        padic_randtest(a, state, ctx);
        if (!_padic_is_zero(a) && padic_val(a) < 2)
        {
            padic_val(a) = 2;
            padic_reduce(a, ctx);
        }

        padic_exp(b, a, ctx);
        padic_log(c, b, ctx);

        result = (padic_equal(a, c, ctx));
        if (!result)
        {
            printf("FAIL (log(exp(x)) == x):\n\n");
            printf("a = "), padic_print(a, ctx), printf("\n");
            printf("b = "), padic_print(b, ctx), printf("\n");
            printf("c = "), padic_print(c, ctx), printf("\n");
            abort();
        }

        _padic_clear(a);
        _padic_clear(b);
        _padic_clear(c);

        padic_ctx_clear(ctx);
    }

/** p > 2 ********************************************************************/

    /* Check aliasing: a = log(a) */
    for (i = 0; i < 10000; i++)
    {
        fmpz_t p;
        long pbits, N;
        padic_ctx_t ctx;

        padic_t a, b;

        fmpz_init(p);
        pbits = n_randint(state, 4) + 2;
        fmpz_set_ui(p, n_randprime(state, pbits, 1));
        N = n_randint(state, 50) + 1;
        padic_ctx_init(ctx, p, N, PADIC_SERIES);

        _padic_init(a);
        _padic_init(b);

        padic_randtest(a, state, ctx);

        if (!_padic_is_zero(a) && (padic_val(a) < 1 || (*p == 2L && padic_val(a) < 2)))
        {
            padic_val(a) = (*p == 2L) + 1;
            padic_reduce(a, ctx);
        }

        _padic_one(b);
        _padic_add(a, a, b, ctx);

        padic_log(b, a, ctx);
        padic_log(a, a, ctx);

        result = (padic_equal(a, b, ctx));
        if (!result)
        {
            printf("FAIL (aliasing):\n\n");
            printf("a = "), padic_print(a, ctx), printf("\n");
            printf("b = "), padic_print(b, ctx), printf("\n");
            abort();
        }

        _padic_clear(a);
        _padic_clear(b);

        fmpz_clear(p);
        padic_ctx_clear(ctx);
    }

    /* Check: log(a) + log(b) == log(a * b) */
    for (i = 0; i < 10000; i++)
    {
        fmpz_t p;
        long pbits, N;
        padic_ctx_t ctx;

        padic_t a, b, c, d, e, f, g;

        fmpz_init(p);
        pbits = n_randint(state, 4) + 2;
        fmpz_set_ui(p, n_randprime(state, pbits, 1));
        N = n_randint(state, 50) + 1;
        padic_ctx_init(ctx, p, N, PADIC_SERIES);

        _padic_init(a);
        _padic_init(b);
        _padic_init(c);
        _padic_init(d);
        _padic_init(e);
        _padic_init(f);
        _padic_init(g);

        padic_randtest(a, state, ctx);

        if (!_padic_is_zero(a) && (padic_val(a) < 1 || (*p == 2L && padic_val(a) < 2)))
        {
            padic_val(a) = (*p == 2L) + 1;
            padic_reduce(a, ctx);
        }

        _padic_one(c);
        _padic_add(a, a, c, ctx);

        padic_randtest(b, state, ctx);

        if (!_padic_is_zero(b) && (padic_val(b) < 1 || (*p == 2L && padic_val(b) < 2)))
        {
            padic_val(b) = (*p == 2L) + 1;
            padic_reduce(b, ctx);
        }

        _padic_one(c);
        _padic_add(b, b, c, ctx);

        padic_mul(c, a, b, ctx);

        padic_log(d, a, ctx);
        padic_log(e, b, ctx);
        padic_add(f, d, e, ctx);

        padic_log(g, c, ctx);

        result = (padic_equal(f, g, ctx));
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

        _padic_clear(a);
        _padic_clear(b);
        _padic_clear(c);
        _padic_clear(d);
        _padic_clear(e);
        _padic_clear(f);
        _padic_clear(g);

        fmpz_clear(p);
        padic_ctx_clear(ctx);
    }

    /* Check: log(exp(x)) == x */
    for (i = 0; i < 10000; i++)
    {
        fmpz_t p;
        long pbits, N;
        padic_ctx_t ctx;

        padic_t a, b, c;

        fmpz_init(p);
        pbits = n_randint(state, 4) + 2;
        fmpz_set_ui(p, n_randprime(state, pbits, 1));
        N = n_randint(state, 50) + 1;
        padic_ctx_init(ctx, p, N, PADIC_SERIES);

        _padic_init(a);
        _padic_init(b);
        _padic_init(c);

        padic_randtest(a, state, ctx);

        if (!_padic_is_zero(a) && (padic_val(a) < 1 || (*p == 2L && padic_val(a) < 2)))
        {
            padic_val(a) = (*p == 2L) + 1;
            padic_reduce(a, ctx);
        }

        padic_exp(b, a, ctx);
        padic_log(c, b, ctx);

        result = (padic_equal(a, c, ctx));
        if (!result)
        {
            printf("FAIL (log(exp(x)) == x):\n\n");
            printf("a = "), padic_print(a, ctx), printf("\n");
            printf("b = "), padic_print(b, ctx), printf("\n");
            printf("c = "), padic_print(c, ctx), printf("\n");
            abort();
        }

        _padic_clear(a);
        _padic_clear(b);
        _padic_clear(c);

        fmpz_clear(p);
        padic_ctx_clear(ctx);
    }

    flint_randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}

