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

    Copyright (C) 2012 Sebastian Pancratz, 2012 Andres Goens

******************************************************************************/

#include <stdio.h>

#include "fq.h"
#include "ulong_extras.h"
#include "long_extras.h"

int
main(void)
{
    int i, result;
    flint_rand_t state;

    printf("trace... ");
    fflush(stdout);

    flint_randinit(state);

    /* Compare with sum of Galois conjugates */
    for (i = 0; i < 2000; i++)
    {
        fmpz_t p;
        long d;
        fq_ctx_t ctx;

        fq_t a, b, c, x;
        long j;

        fmpz_init(p);
        fmpz_set_ui(p, n_randprime(state, 2 + n_randint(state, 3), 1));
        d = n_randint(state, 10) + 1;
        fq_ctx_init_conway(ctx, p, d, "a", PADIC_TERSE);

        fq_init(a);
        fq_init(b);
        fq_init(c);
        fq_init(x);

        fq_randtest(a, state, ctx);

        fq_trace(x, a, ctx);

        fq_zero(b);
        for (j = 0; j < d; j++)
        {
            fq_frobenius(c, a, j, ctx);
            fq_add(b, b, c, ctx);
        }

        result = fq_equal(x, b);
        if (!result)
        {
            printf("FAIL:\n\n");
            printf("a = "), fq_print_pretty(a, ctx), printf("\n");
            printf("b = "), fq_print_pretty(b, ctx), printf("\n");
            printf("x = "), fq_print_pretty(x,ctx), printf("\n");
            for (j = 0; j < d; j++)
            {
                fq_frobenius(c, a, j, ctx);
                printf("sigma^%ld = ", j), fq_print_pretty(c, ctx), printf("\n");
            }
            abort();
        }

        fq_clear(a);
        fq_clear(b);
        fq_clear(c);
        fq_clear(x);


        fmpz_clear(p);
        fq_ctx_clear(ctx);
    }

    flint_randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}

