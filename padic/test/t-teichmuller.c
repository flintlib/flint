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

    printf("teichmuller... ");
    fflush(stdout);

    flint_randinit(state);

    /* Check aliasing */
    for (i = 0; i < 10000; i++)
    {
        fmpz_t p;
        long N;
        padic_ctx_t ctx;

        padic_t a, b, c;

        fmpz_init(p);
        fmpz_set_ui(p, n_randprime(state, 5, 1));
        N = n_randint(state, 100) + 1;
        padic_ctx_init(ctx, p, N, PADIC_SERIES);

        padic_init(a, ctx);
        padic_init(b, ctx);
        padic_init(c, ctx);

        padic_randtest_not_zero(a, state, ctx);
        padic_val(a) = FLINT_ABS(padic_val(a));
        padic_set(b, a, ctx);

        padic_teichmuller(c, b, ctx);
        padic_teichmuller(b, b, ctx);

        result = (padic_equal(b, c, ctx));
        if (!result)
        {
            printf("FAIL:\n\n");
            printf("a = "), padic_print(a, ctx), printf("\n");
            printf("b = "), padic_print(b, ctx), printf("\n");
            printf("c = "), padic_print(c, ctx), printf("\n");
            abort();
        }

        padic_clear(a, ctx);
        padic_clear(b, ctx);
        padic_clear(c, ctx);

        fmpz_clear(p);
        padic_ctx_clear(ctx);
    }

    /* Check correctness */
    for (i = 0; i < 10000; i++)
    {
        fmpz_t p;
        long N;
        padic_ctx_t ctx;

        padic_t a, b, c;

        fmpz_init(p);
        fmpz_set_ui(p, n_randprime(state, 5, 1));
        N = n_randint(state, 100) + 1;
        padic_ctx_init(ctx, p, N, PADIC_SERIES);

        padic_init(a, ctx);
        padic_init(b, ctx);
        padic_init(c, ctx);

        padic_randtest(a, state, ctx);
        padic_val(a) = FLINT_ABS(padic_val(a));

        padic_teichmuller(b, a, ctx);

        padic_pow_si(c, b, *p, ctx);  /* XXX.  Assumes p is small */
        padic_sub(c, c, b, ctx);

        result = (padic_is_zero(c, ctx));
        if (!result)
        {
            printf("FAIL:\n\n");
            printf("a = "), padic_print(a, ctx), printf("\n");
            printf("b = "), padic_print(b, ctx), printf("\n");
            printf("c = "), padic_print(c, ctx), printf("\n");
            abort();
        }

        padic_clear(a, ctx);
        padic_clear(b, ctx);
        padic_clear(c, ctx);

        fmpz_clear(p);
        padic_ctx_clear(ctx);
    }

    flint_randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}

