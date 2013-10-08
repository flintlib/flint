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

    printf("pow... ");
    fflush(stdout);

    flint_randinit(state);

    /* Check aliasing  */
    for (i = 0; i < 200; i++)
    {
        fmpz_t p;
        long d, len;
        fq_ctx_t ctx;

        fq_poly_t a, b, c;
        ulong exp;

        fmpz_init(p);
        fmpz_set_ui(p, n_randprime(state, 2 + n_randint(state, 3), 1));
        d = n_randint(state, 10) + 1;
        len = n_randint(state, 15) + 1;
        fq_ctx_init_conway(ctx, p, d, "a");
        fq_poly_init(a);
        fq_poly_init(b);
        fq_poly_init(c);

        fq_poly_randtest(a, state, len, ctx);
        exp = n_randtest(state) % 20UL;
        fq_poly_set(b, a);
        fq_poly_pow(c, b, exp, ctx);
        fq_poly_pow(b, b, exp, ctx);

        result = (fq_poly_equal(b, c));
        if (!result)
        {
            printf("FAIL:\n\n");
            printf("a = "), fq_poly_print_pretty(a, "X", ctx), printf("\n");
            printf("b = "), fq_poly_print_pretty(b, "X", ctx), printf("\n");
            printf("c = "), fq_poly_print_pretty(c, "X", ctx), printf("\n");
            printf("exp = %lu\n", exp);
            abort();
        }

        fq_poly_clear(a);
        fq_poly_clear(b);
        fq_poly_clear(c);

        fmpz_clear(p);
        fq_ctx_clear(ctx);
    }

    /* Compare with repeated multiplications by the base */
    for (i = 0; i < 1000; i++)
    {
        fmpz_t p;
        long d, len;
        fq_ctx_t ctx;

        fq_poly_t a, b, c;
        ulong exp;

        fmpz_init(p);
        fmpz_set_ui(p, n_randprime(state, 2 + n_randint(state, 3), 1));
        d   = n_randint(state, 10) + 1;
        len = n_randint(state, 15) + 1;
        fq_ctx_init_conway(ctx, p, d, "a");
        fq_poly_init(a);
        fq_poly_init(b);
        fq_poly_init(c);

        fq_poly_randtest(b, state, len, ctx);
        exp = n_randtest(state) % 20UL;

        fq_poly_pow(a, b, exp, ctx);

        if (exp == 0)
        {
            fq_poly_one(c);
        }
        else
        {
            long j;

            fq_poly_set(c, b);
            for (j = 1; j < exp; j++)
                fq_poly_mul(c, c, b, ctx);
        }

        result = (fq_poly_equal(a, c));
        if (!result)
        {
            printf("FAIL:\n\n");
            printf("a = "), fq_poly_print_pretty(a, "X", ctx), printf("\n");
            printf("b = "), fq_poly_print_pretty(b, "X", ctx), printf("\n");
            printf("c = "), fq_poly_print_pretty(c, "X", ctx), printf("\n");
            printf("exp = %lu\n", exp);
            fq_ctx_print(ctx);
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

