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
    Copyright (C) 2012 Andres Goens

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

    printf("sqr_KS... ");
    fflush(stdout);

    flint_randinit(state);

    /* Check aliasing: a = a * a */
    for (i = 0; i < 2000; i++)
    {
        fmpz_t p;
        long d, len;
        fq_ctx_t ctx;

        fq_poly_t a, c;

        fmpz_init(p);
        fmpz_set_ui(p, n_randprime(state, 2 + n_randint(state, 3), 1));
        d   = n_randint(state, 10) + 1;
        len = n_randint(state, 50) + 1;
        fq_ctx_init_conway(ctx, p, d, "a");
        fq_poly_init(a);
        fq_poly_init(c);

        fq_poly_randtest(a, state, len, ctx);

        fq_poly_sqr_KS(c, a, ctx);
        fq_poly_sqr_KS(a, a, ctx);

        result = (fq_poly_equal(a, c));
        if (!result)
        {
            printf("FAIL (aliasing):\n\n");
            printf("a = "), fq_poly_print_pretty(a, "X", ctx), printf("\n");
            printf("c = "), fq_poly_print_pretty(c, "X", ctx), printf("\n");
            abort();
        }

        fq_poly_clear(a);
        fq_poly_clear(c);

        fmpz_clear(p);
        fq_ctx_clear(ctx);
    }

    /* Check that a^2 + a^2 == a * (a + a) */
    for (i = 0; i < 2000; i++)
    {
        fmpz_t p;
        long deg, len;
        fq_ctx_t ctx;

        fq_poly_t a, b, c, d;

        fmpz_init(p);
        fmpz_set_ui(p, n_randprime(state, 2 + n_randint(state, 3), 1));
        deg = n_randint(state, 10) + 1;
        len = n_randint(state, 50) + 1;
        fq_ctx_init_conway(ctx, p, deg, "a");
        fq_poly_init(a);
        fq_poly_init(b);
        fq_poly_init(c);
        fq_poly_init(d);

        fq_poly_randtest(a, state, len, ctx);

        fq_poly_sqr_KS(b, a, ctx);
        fq_poly_add(c, b, b, ctx);

        fq_poly_add(d, a, a, ctx);
        fq_poly_mul(d, a, d, ctx);

        result = (fq_poly_equal(c, d));
        if (!result)
        {
            printf("FAIL (a^2 + a^2 == a(a + a)):\n\n");
            printf("a = "), fq_poly_print_pretty(a, "X", ctx), printf("\n");
            printf("b = "), fq_poly_print_pretty(b, "X", ctx), printf("\n");
            printf("c = "), fq_poly_print_pretty(c, "X", ctx), printf("\n");
            printf("d = "), fq_poly_print_pretty(d, "X", ctx), printf("\n");
            abort();
        }

        fq_poly_clear(a);
        fq_poly_clear(b);
        fq_poly_clear(c);
        fq_poly_clear(d);

        fmpz_clear(p);
        fq_ctx_clear(ctx);
    }

    /* Compare mul() */
    for (i = 0; i < 2000; i++)
    {
        fmpz_t p;
        long deg, len;
        fq_ctx_t ctx;

        fq_poly_t a, b, c;

        fmpz_init(p);
        fmpz_set_ui(p, n_randprime(state, 2 + n_randint(state, 3), 1));
        deg = n_randint(state, 10) + 1;
        len = n_randint(state, 50) + 1;
        fq_ctx_init_conway(ctx, p, deg, "a");
        fq_poly_init(a);
        fq_poly_init(b);
        fq_poly_init(c);

        fq_poly_randtest(a, state, len, ctx);

        fq_poly_mul(b, a, a, ctx);
        fq_poly_sqr_KS(c, a, ctx);

        result = (fq_poly_equal(b, c));
        if (!result)
        {
            printf("FAIL (cmp with mul):\n\n");
            printf("a = "), fq_poly_print_pretty(a, "X", ctx), printf("\n");
            printf("b = "), fq_poly_print_pretty(b, "X", ctx), printf("\n");
            printf("c = "), fq_poly_print_pretty(c, "X", ctx), printf("\n");
            abort();
        }

        fq_poly_clear(a);
        fq_poly_clear(b);
        fq_poly_clear(c);

        fmpz_clear(p);
        fq_ctx_clear(ctx);
    }

    flint_randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}

