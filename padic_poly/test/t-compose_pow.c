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
#include "long_extras.h"
#include "ulong_extras.h"
#include "padic_poly.h"

int
main(void)
{
    int i, result;
    flint_rand_t state;

    padic_ctx_t ctx;
    fmpz_t p;
    long N;

    printf("compose_pow... ");
    fflush(stdout);

    flint_randinit(state);

    /* Aliasing */
    for (i = 0; i < 1000; i++)
    {
        padic_poly_t a, b, c;
        long k;

        fmpz_init(p);
        fmpz_set_ui(p, n_randprime(state, 5, 1));
        N = z_randint(state, 50);
        padic_ctx_init(ctx, p, N, PADIC_SERIES);

        padic_poly_init(a);
        padic_poly_init(b);
        padic_poly_init(c);

        padic_poly_randtest(a, state, n_randint(state, 100), ctx);
        padic_poly_set(b, a);
        k = n_randint(state, 20) + 1;

        padic_poly_compose_pow(c, b, k, ctx);
        padic_poly_compose_pow(b, b, k, ctx);

        result = (padic_poly_equal(b, c));
        if (!result)
        {
            printf("FAIL:\n");
            padic_poly_print(a, ctx), printf("\n\n");
            padic_poly_print(b, ctx), printf("\n\n");
            padic_poly_print(c, ctx), printf("\n\n");
            abort();
        }

        padic_poly_clear(a);
        padic_poly_clear(b);
        padic_poly_clear(c);

        padic_ctx_clear(ctx);
        fmpz_clear(p);
    }

    /* Compare with usual composition */
    for (i = 0; i < 1000; i++)
    {
        padic_poly_t f, g, h1, h2;
        long k;
        padic_t one;

        fmpz_init(p);
        fmpz_set_ui(p, n_randprime(state, 5, 1));
        N = z_randint(state, 50);
        padic_ctx_init(ctx, p, N, PADIC_VAL_UNIT);

        padic_poly_init(f);
        padic_poly_init(g);
        padic_poly_init(h1);
        padic_poly_init(h2);

        padic_poly_randtest(f, state, n_randint(state, 40), ctx);
        k = n_randint(state, 20) + 1;

        padic_poly_compose_pow(h1, f, k, ctx);

        _padic_init(one);
        _padic_one(one);
        padic_poly_fit_length(g, k + 1);
        padic_poly_set_coeff_padic(g, k, one, ctx);
        _padic_poly_set_length(g, k + 1);
        _padic_clear(one);
        padic_poly_compose(h2, f, g, ctx);

        result = (padic_poly_equal(h1, h2));
        if (!result)
        {
            printf("FAIL (cmp with composition):\n");
            printf("f  = "), padic_poly_print(f, ctx), printf("\n\n");
            printf("g  = "), padic_poly_print(g, ctx), printf("\n\n");
            printf("h1 = "), padic_poly_print(h1, ctx), printf("\n\n");
            printf("h2 = "), padic_poly_print(h2, ctx), printf("\n\n");
            printf("p = "), fmpz_print(p), printf("\n\n");
            printf("N = %ld\n\n", ctx->N);
            abort();
        }

        padic_poly_clear(f);
        padic_poly_clear(g);
        padic_poly_clear(h1);
        padic_poly_clear(h2);

        padic_ctx_clear(ctx);
        fmpz_clear(p);
    }

    flint_randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}
