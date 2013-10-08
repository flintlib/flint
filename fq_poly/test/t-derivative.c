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

    printf("derivative... ");
    fflush(stdout);

    flint_randinit(state);

    /* Check aliasing */
    for (i = 0; i < 2000; i++)
    {
        fmpz_t p;
        long d,len;
        fq_ctx_t ctx;

        fq_poly_t a, b, c;

        fmpz_init(p);
        fmpz_set_ui(p, n_randprime(state, 2 + n_randint(state, 3), 1));
        d = n_randint(state, 10) + 1;
        len = n_randint(state, 15) + 1;
        fq_ctx_init_conway(ctx, p, d, "a");
        fq_poly_init(a);
        fq_poly_init(b);
        fq_poly_init(c);

        fq_poly_randtest(a, state, len, ctx);
        fq_poly_set(b, a);
        fq_poly_derivative(c, b, ctx);
        fq_poly_derivative(b, b, ctx);

        result = (fq_poly_equal(b, c));
        if (!result)
        {
            printf("FAIL:\n\n");
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

    /* Check constants have derivative zero */
    for (i = 0; i < 2000; i++)
    {
        fmpz_t p;
        long d;
        fq_ctx_t ctx;

        fq_poly_t a, b;

        fmpz_init(p);
        fmpz_set_ui(p, n_randprime(state, 2 + n_randint(state, 3), 1));
        d = n_randint(state, 10) + 1;
        fq_ctx_init_conway(ctx, p, d, "a");
        fq_poly_init(a);
        fq_poly_init(b);

        fq_poly_randtest(a, state, n_randint(state, 2), ctx);
        fq_poly_derivative(b, a, ctx);

        result = (fq_poly_is_zero(b));
        if (!result)
        {
            printf("FAIL:\n\n");
            printf("a = "), fq_poly_print_pretty(a, "X", ctx), printf("\n");
            printf("b = "), fq_poly_print_pretty(b, "X", ctx), printf("\n");
            abort();
        }

        fq_poly_clear(a);
        fq_poly_clear(b);

        fmpz_clear(p);
        fq_ctx_clear(ctx);
    }

    /* Check (f g)' == f' g + f g'  */
    for (i = 0; i < 2000; i++)
    {
        fmpz_t p;
        long deg;
        fq_ctx_t ctx;

        fq_poly_t a, b, c, d, lhs, rhs;

        fmpz_init(p);
        fmpz_set_ui(p, n_randprime(state, 2 + n_randint(state, 3), 1));
        deg = n_randint(state, 10) + 1;
        fq_ctx_init_conway(ctx, p, deg, "a");
        fq_poly_init(a);
        fq_poly_init(b);
        fq_poly_init(c);
        fq_poly_init(d);
        fq_poly_init(lhs);
        fq_poly_init(rhs);

        fq_poly_randtest(a, state, n_randint(state, 100), ctx);
        fq_poly_randtest(b, state, n_randint(state, 100), ctx);

        fq_poly_mul(lhs, a, b, ctx);
        fq_poly_derivative(lhs, lhs, ctx);
        fq_poly_derivative(c, a, ctx);
        fq_poly_derivative(d, b, ctx);
        fq_poly_mul(c, c, b, ctx);
        fq_poly_mul(d, a, d, ctx);
        fq_poly_add(rhs, c, d, ctx);

        result = (fq_poly_equal(lhs, rhs));
        if (!result)
        {
            printf("FAIL:\n\n");
            printf("a   = "), fq_poly_print_pretty(a,   "X", ctx), printf("\n");
            printf("b   = "), fq_poly_print_pretty(b,   "X", ctx), printf("\n");
            printf("c   = "), fq_poly_print_pretty(c,   "X", ctx), printf("\n");
            printf("d   = "), fq_poly_print_pretty(d,   "X", ctx), printf("\n");
            printf("lhs = "), fq_poly_print_pretty(lhs, "X", ctx), printf("\n");
            printf("rhs = "), fq_poly_print_pretty(rhs, "X", ctx), printf("\n");
            abort();
        }

        fq_poly_clear(a);
        fq_poly_clear(b);
        fq_poly_clear(c);
        fq_poly_clear(d);
        fq_poly_clear(lhs);
        fq_poly_clear(rhs);

        fmpz_clear(p);
        fq_ctx_clear(ctx);
    }

    flint_randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}

