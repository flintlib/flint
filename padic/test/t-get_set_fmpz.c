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

    printf("get_set_fmpz... ");
    fflush(stdout);

    flint_randinit(state);

    /* Check that Zp(QQ(x)) == x. */
    for (i = 0; i < 10000; i++)
    {
        fmpz_t p;
        long N;
        padic_ctx_t ctx;

        padic_t a, b;
        fmpz_t c;

        fmpz_init(p);
        fmpz_set_ui(p, n_randprime(state, 5, 1));
        N = n_randint(state, 50) + 1;
        padic_ctx_init(ctx, p, N, PADIC_SERIES);

        padic_init(a, ctx);
        padic_init(b, ctx);
        fmpz_init(c);

        padic_randtest(a, state, ctx);
        if (padic_val(a) < 0)
        {
            padic_val(a) = - padic_val(a);
            padic_normalise(a, ctx);
        }

        padic_get_fmpz(c, a, ctx);
        padic_set_fmpz(b, c, ctx);

        result = (padic_equal(a, b, ctx));
        if (!result)
        {
            printf("FAIL:\n\n");
            printf("a = "), padic_print(a, ctx), printf("\n");
            printf("b = "), padic_print(b, ctx), printf("\n");
            printf("c = "), fmpz_print(c), printf("\n");
            abort();
        }

        padic_clear(a, ctx);
        padic_clear(b, ctx);
        fmpz_clear(c);

        fmpz_clear(p);
        padic_ctx_clear(ctx);
    }

    flint_randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}

