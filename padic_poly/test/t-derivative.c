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

    Copyright (C) 2009 William Hart
    Copyright (C) 2012 Sebastian Pancratz

******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include "padic_poly.h"
#include "ulong_extras.h"
#include "long_extras.h"

int
main(void)
{
    int i, result;
    flint_rand_t state;

    padic_ctx_t ctx;
    fmpz_t p;
    long N;

    printf("derivative... ");
    fflush(stdout);

    flint_randinit(state);

    /* Aliasing */
    for (i = 0; i < 10000; i++)
    {
        padic_poly_t a, b, c;

        fmpz_init(p);
        fmpz_set_ui(p, n_randprime(state, 5, 1));
        N = z_randint(state, 50);
        padic_ctx_init(ctx, p, N, PADIC_SERIES);

        padic_poly_init(a);
        padic_poly_init(b);
        padic_poly_init(c);

        padic_poly_randtest(a, state, n_randint(state, 100), ctx);
        padic_poly_set(b, a);

        padic_poly_derivative(c, b, ctx);
        padic_poly_derivative(b, b, ctx);

        result = (padic_poly_equal(b, c));
        if (!result)
        {
            printf("FAIL (alias):\n");
            printf("a = "), padic_poly_print(a, ctx), printf("\n\n");
            printf("b = "), padic_poly_print(b, ctx), printf("\n\n");
            printf("c = "), padic_poly_print(c, ctx), printf("\n\n");
            abort();
        }

        padic_poly_clear(a);
        padic_poly_clear(b);
        padic_poly_clear(c);

        padic_ctx_clear(ctx);
        fmpz_clear(p);
    }

    /* Compare with derivative over QQ */
    for (i = 0; i < 10000; i++)
    {
        padic_poly_t a, b, c;
        fmpq_poly_t aQQ, bQQ;

        fmpz_init(p);
        fmpz_set_ui(p, n_randprime(state, 5, 1));
        N = z_randint(state, 50);
        padic_ctx_init(ctx, p, N, PADIC_SERIES);

        padic_poly_init(a);
        padic_poly_init(b);
        padic_poly_init(c);

        padic_poly_randtest(a, state, n_randint(state, 100), ctx);

        fmpq_poly_init(aQQ);
        fmpq_poly_init(bQQ);

        padic_poly_derivative(b, a, ctx);

        padic_poly_get_fmpq_poly(aQQ, a, ctx);
        fmpq_poly_derivative(bQQ, aQQ);
        padic_poly_set_fmpq_poly(c, bQQ, ctx);

        result = (padic_poly_equal(b, c));
        if (!result)
        {
            printf("FAIL (cmp with QQ):\n");
            printf("a = "), padic_poly_print(a, ctx), printf("\n\n");
            printf("b = "), padic_poly_print(b, ctx), printf("\n\n");
            printf("c = "), padic_poly_print(c, ctx), printf("\n\n");
            printf("aQQ = "), fmpq_poly_print(aQQ), printf("\n\n");
            printf("bQQ = "), fmpq_poly_print(bQQ), printf("\n\n");
            abort();
        }

        padic_poly_clear(a);
        padic_poly_clear(b);
        padic_poly_clear(c);
        fmpq_poly_clear(aQQ);
        fmpq_poly_clear(bQQ);

        padic_ctx_clear(ctx);
        fmpz_clear(p);
    }

    flint_randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}

