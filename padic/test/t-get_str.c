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
#include <string.h>
#include <mpir.h>
#include "flint.h"
#include "ulong_extras.h"
#include "long_extras.h"
#include "fmpq.h"
#include "padic.h"

int
main(void)
{
    int i, result;
    flint_rand_t state;

    printf("get_str... ");
    fflush(stdout);

    flint_randinit(state);

    for (i = 0; i < 10000; i++)
    {
        fmpz_t p;
        long N;
        padic_ctx_t ctx;

        padic_t x;
        fmpq_t y;

        char *s, *t;

        fmpz_init(p);
        fmpz_set_ui(p, n_randprime(state, 5, 1));
        N = z_randint(state, 50) + 1;
        padic_ctx_init(ctx, p, N, PADIC_TERSE);

        padic_init(x, ctx);
        fmpq_init(y);

        padic_randtest(x, state, ctx);
        _padic_get_fmpq(y, x, ctx);

        s = _padic_get_str(NULL, x, ctx);
        t = fmpq_get_str(NULL, 10, y);

        result = strcmp(s, t) == 0;
        if (!result)
        {
            printf("FAIL:\n\n");
            printf("x = "), padic_print(x, ctx), printf("\n");
            printf("y = "), fmpq_clear(y), printf("\n");
            abort();
        }

        free(s);
        free(t);
        _padic_clear(x);
        fmpq_clear(y);
        padic_ctx_clear(ctx);
    }

    flint_randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}

