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

    printf("compose_horner... ");
    fflush(stdout);

    flint_randinit(state);

    /* Check aliasing of the first argument */
    for (i = 0; i < 50; i++)
    {
        fmpz_t p;
        long d, len;
        fq_ctx_t ctx;

        fq_poly_t f, g, h;

        fmpz_init(p);
        fmpz_set_ui(p, n_randprime(state, 2 + n_randint(state, 3), 1));
        d   = n_randint(state, 10) + 1;
        len = n_randint(state, 15) + 1;
        fq_ctx_init_conway(ctx, p, d, "a");
        fq_poly_init(f);
        fq_poly_init(g);
        fq_poly_init(h);

        fq_poly_randtest(f, state, n_randint(state, 40), ctx);
        fq_poly_randtest(g, state, n_randint(state, 20), ctx);

        fq_poly_compose_horner(h, f, g, ctx);
        fq_poly_compose_horner(f, f, g, ctx);

        result = (fq_poly_equal(f, h));
        if (!result)
        {
            printf("FAIL:\n\n");
            printf("f = "), fq_poly_print_pretty(f, "X", ctx), printf("\n");
            printf("g = "), fq_poly_print_pretty(g, "X", ctx), printf("\n");
            printf("h = "), fq_poly_print_pretty(h, "X", ctx), printf("\n");
            abort();
        }

        fq_poly_clear(f);
        fq_poly_clear(g);
        fq_poly_clear(h);

        fmpz_clear(p);
        fq_ctx_clear(ctx);
    }

    /* Check aliasing of the second argument */
    for (i = 0; i < 50; i++)
    {
        fmpz_t p;
        long d, len;
        fq_ctx_t ctx;

        fq_poly_t f, g, h;

        fmpz_init(p);
        fmpz_set_ui(p, n_randprime(state, 2 + n_randint(state, 3), 1));
        d   = n_randint(state, 10) + 1;
        len = n_randint(state, 15) + 1;
        fq_ctx_init_conway(ctx, p, d, "a");
        fq_poly_init(f);
        fq_poly_init(g);
        fq_poly_init(h);

        fq_poly_randtest(f, state, n_randint(state, 40), ctx);
        fq_poly_randtest(g, state, n_randint(state, 20), ctx);

        fq_poly_compose_horner(h, f, g, ctx);
        fq_poly_compose_horner(g, f, g, ctx);

        result = (fq_poly_equal(g, h));
        if (!result)
        {
            printf("FAIL:\n\n");
            printf("f = "), fq_poly_print_pretty(f, "X", ctx), printf("\n");
            printf("g = "), fq_poly_print_pretty(g, "X", ctx), printf("\n");
            printf("h = "), fq_poly_print_pretty(h, "X", ctx), printf("\n");
            abort();
        }

        fq_poly_clear(f);
        fq_poly_clear(g);
        fq_poly_clear(h);

        fmpz_clear(p);
        fq_ctx_clear(ctx);
    }

    /* Compare with the naive method */
    for (i = 0; i < 50; i++)
    {
        fmpz_t p;
        long d, len;
        fq_ctx_t ctx;

        fq_poly_t f, g, h, s, t;
        long k;

        fmpz_init(p);
        fmpz_set_ui(p, n_randprime(state, 2 + n_randint(state, 3), 1));
        d   = n_randint(state, 10) + 1;
        len = n_randint(state, 15) + 1;
        fq_ctx_init_conway(ctx, p, d, "a");
        fq_poly_init(f);
        fq_poly_init(g);
        fq_poly_init(h);
        fq_poly_init(s);
        fq_poly_init(t);

        fq_poly_randtest(g, state, n_randint(state, 40), ctx);
        fq_poly_randtest(h, state, n_randint(state, 20), ctx);

        fq_poly_one(t);
        for (k = 0; k < fq_poly_length(g); k++)
        {
            fq_poly_scalar_addmul_fq(s, t, g->coeffs + k, ctx);
            fq_poly_mul(t, t, h, ctx);
        }

        fq_poly_compose_horner(f, g, h, ctx);

        result = (fq_poly_equal(f, s));
        if (!result)
        {
            printf("FAIL:\n\n");
            printf("f = "), fq_poly_print_pretty(f, "X", ctx), printf("\n");
            printf("g = "), fq_poly_print_pretty(g, "X", ctx), printf("\n");
            printf("h = "), fq_poly_print_pretty(h, "X", ctx), printf("\n");
            printf("s = "), fq_poly_print_pretty(s, "X", ctx), printf("\n");
            printf("t = "), fq_poly_print_pretty(t, "X", ctx), printf("\n");
            abort();
        }

        fq_poly_clear(f);
        fq_poly_clear(g);
        fq_poly_clear(h);
        fq_poly_clear(s);
        fq_poly_clear(t);

        fmpz_clear(p);
        fq_ctx_clear(ctx);
    }

    flint_randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}

