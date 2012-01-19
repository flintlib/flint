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

    printf("shift_left/right... ");
    fflush(stdout);

    flint_randinit(state);

    /* Aliasing for left shift */
    for (i = 0; i < 10000; i++)
    {
        padic_poly_t a, b, c;
        long shift = n_randint(state, 100);

        fmpz_init(p);
        fmpz_set_ui(p, n_randprime(state, 5, 1));
        N = z_randint(state, 50);
        padic_ctx_init(ctx, p, N, PADIC_SERIES);

        padic_poly_init(a);
        padic_poly_init(b);
        padic_poly_init(c);

        padic_poly_randtest(a, state, n_randint(state, 100), ctx);
        padic_poly_set(b, a);

        padic_poly_shift_left(c, b, shift);
        padic_poly_shift_left(b, b, shift);

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

    /* Aliasing for shift right */
    for (i = 0; i < 10000; i++)
    {
        padic_poly_t a, b, c;
        long shift = n_randint(state, 100);

        fmpz_init(p);
        fmpz_set_ui(p, n_randprime(state, 5, 1));
        N = z_randint(state, 50);
        padic_ctx_init(ctx, p, N, PADIC_SERIES);

        padic_poly_init(a);
        padic_poly_init(b);
        padic_poly_init(c);

        padic_poly_randtest(a, state, n_randint(state, 100), ctx);
        padic_poly_set(b, a);

        padic_poly_shift_right(c, b, shift, ctx);
        padic_poly_shift_right(b, b, shift, ctx);

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

    /* Check shift left then right does nothing */
    for (i = 0; i < 10000; i++)
    {
        padic_poly_t a, b, c;
        long shift = n_randint(state, 100);

        fmpz_init(p);
        fmpz_set_ui(p, n_randprime(state, 5, 1));
        N = z_randint(state, 50);
        padic_ctx_init(ctx, p, N, PADIC_SERIES);

        padic_poly_init(a);
        padic_poly_init(b);
        padic_poly_init(c);

        padic_poly_randtest(a, state, n_randint(state, 100), ctx);

        padic_poly_shift_left(b, a, shift);
        padic_poly_shift_right(c, b, shift, ctx);

        result = (padic_poly_equal(a, c));
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

    flint_randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}
