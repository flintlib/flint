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

    Copyright (C) 2011 Sebastian Pancratz

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

/** p > 2 ********************************************************************/

    /* Check aliasing: a = log(a) */
    for (i = 0; i < 10000; i++)
    {
        fmpz_t p;
        long N;
        padic_ctx_t ctx;

        padic_t a, b;

        fmpz_init(p);
        fmpz_set_ui(p, n_randprime(state, 5, 1));
        N = n_randint(state, 50) + 1;
        padic_ctx_init(ctx, p, N, PADIC_SERIES);

        _padic_init(a);
        _padic_init(b);

        padic_randtest_not_zero(a, state, ctx);
        if (padic_val(a) < 1)
            padic_val(a) = 1;
        _padic_one(b);
        _padic_add(a, a, b, ctx);

        padic_log(b, a, ctx);
        padic_log(a, a, ctx);

        result = (padic_equal(a, b, ctx));
        if (!result)
        {
            printf("FAIL (aliasing):\n\n");
            printf("a = "), _padic_print(a, ctx), printf("\n");
            printf("b = "), _padic_print(b, ctx), printf("\n");
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
        long N;
        padic_ctx_t ctx;

        padic_t a, b, c, d, e, f, g;

        fmpz_init(p);
        fmpz_set_ui(p, n_randprime(state, 5, 1));
        N = n_randint(state, 50) + 1;
        padic_ctx_init(ctx, p, N, PADIC_SERIES);

        _padic_init(a);
        _padic_init(b);
        _padic_init(c);
        _padic_init(d);
        _padic_init(e);
        _padic_init(f);
        _padic_init(g);

        padic_randtest_not_zero(a, state, ctx);
        if (padic_val(a) < 1) padic_val(a) = 1;
        _padic_one(c);
        _padic_add(a, a, c, ctx);

        padic_randtest_not_zero(b, state, ctx);
        if (padic_val(b) < 1) padic_val(b) = 1;
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

    flint_randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}

