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

#include "fq.h"
#include "fq_poly.h"

#include "ulong_extras.h"
#include "long_extras.h"

int
main(void)
{
    int i, result;
    flint_rand_t state;

    printf("neg... ");
    fflush(stdout);

    flint_randinit(state);

    /* Check aliasing: a = -a */
    for (i = 0; i < 2000; i++)
    {
        fmpz_t p;
        long d,len;
        fq_ctx_t ctx;

        fq_poly_t a, b;

        fmpz_init(p);
        fmpz_set_ui(p, n_randprime(state, 2 + n_randint(state, 3), 1));
        d = n_randint(state, 10) + 1;
        len = n_randint(state, 15) + 1;
        fq_ctx_init_conway(ctx, p, d,"a", PADIC_SERIES);
        fq_poly_init(a,ctx);
        fq_poly_init(b,ctx);

        fq_poly_randtest(a, ctx, state,len);

        fq_poly_neg(b, a);
        fq_poly_neg(a, a);

        result = (fq_poly_equal(a, b));
        if (!result)
        {
            printf("FAIL:\n\n");
            /*
            printf("a = "), fq_print_pretty(a, ctx), printf("\n");
            printf("b = "), fq_print_pretty(b, ctx), printf("\n");
            printf("c = "), fq_print_pretty(c, ctx), printf("\n");
            */
            abort();
        }

        fq_poly_clear(a);
        fq_poly_clear(b);

        fmpz_clear(p);
        fq_ctx_clear(ctx);
    }


    /* Check that -(-a) == a */
    for (i = 0; i < 2000; i++)
    {
        fmpz_t p;
        long d,len;
        fq_ctx_t ctx;

        fq_poly_t a, b;

        fmpz_init(p);
        fmpz_set_ui(p, n_randprime(state, 2 + n_randint(state, 3), 1));
        d = n_randint(state, 10) + 1;
        len = n_randint(state, 15) + 1;
        fq_ctx_init_conway(ctx, p, d,"a", PADIC_SERIES);
        fq_poly_init(a,ctx);
        fq_poly_init(b,ctx);

        fq_poly_randtest(a, ctx, state,len);

        fq_poly_neg(b, a);
        fq_poly_neg(b, b);

        result = (fq_poly_equal(a , b));
        if (!result)
        {
            printf("FAIL:\n\n");
            /*
            printf("a = "), fq_print_pretty(a, ctx), printf("\n");
            printf("b = "), fq_print_pretty(b, ctx), printf("\n");
            printf("c = "), fq_print_pretty(c, ctx), printf("\n");
            */
            abort();
        }

        fq_poly_clear(a);
        fq_poly_clear(b);

        fmpz_clear(p);
        fq_ctx_clear(ctx);
    }


    /* Check that (a + (-a)) == 0 */
    for (i = 0; i < 2000; i++)
    {
        fmpz_t p;
        long d,len;
        fq_ctx_t ctx;

        fq_poly_t a, b;

        fmpz_init(p);
        fmpz_set_ui(p, n_randprime(state, 2 + n_randint(state, 3), 1));
        d = n_randint(state, 10) + 1;
        len = n_randint(state, 15) + 1;
        fq_ctx_init_conway(ctx, p, d,"a", PADIC_SERIES);
        fq_poly_init(a,ctx);
        fq_poly_init(b,ctx);

        fq_poly_randtest(a, ctx, state,len);

        fq_poly_neg(b, a);
        fq_poly_add(a, a, b);


        result = (fq_poly_is_zero(a));
        if (!result)
        {
            printf("FAIL:\n\n");
            /*
            printf("a = "), fq_print_pretty(a, ctx), printf("\n");
            printf("b = "), fq_print_pretty(b, ctx), printf("\n");
            printf("c = "), fq_print_pretty(c, ctx), printf("\n");
            */
            abort();
        }

        fq_poly_clear(a);
        fq_poly_clear(b);

        fmpz_clear(p);
        fq_ctx_clear(ctx);
    }


    flint_randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}

