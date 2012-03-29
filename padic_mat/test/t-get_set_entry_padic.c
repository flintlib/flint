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
#include "long_extras.h"
#include "padic.h"
#include "padic_mat.h"

int
main(void)
{
    int i, result;
    flint_rand_t state;

    printf("get/ set_entry_padic... ");
    fflush(stdout);

    flint_randinit(state);

    for (i = 0; i < 10000; i++)
    {
        fmpz_t p;
        long N;
        padic_ctx_t ctx;
        long m, n;

        padic_mat_t a;
        padic_t x, y;
        long r, c;

        fmpz_init(p);
        fmpz_set_ui(p, n_randprime(state, 5, 1));
        N = z_randint(state, 50);
        padic_ctx_init(ctx, p, N, PADIC_VAL_UNIT);

        m = n_randint(state, 20) + 1;
        n = n_randint(state, 20) + 1;

        padic_mat_init(a, m, n);
        _padic_init(x);
        _padic_init(y);

        padic_mat_randtest(a, state, ctx);
        padic_randtest_not_zero(x, state, ctx);

        r = n_randint(state, m);
        c = n_randint(state, n);

        padic_mat_set_entry_padic(a, r, c, x, ctx);
        padic_mat_get_entry_padic(y, a, r, c, ctx);

        result = (_padic_equal(x, y) && _padic_mat_is_canonical(a, p));
        if (!result)
        {
            printf("FAIL:\n\n");
            printf("a = "), padic_mat_print(a, ctx), printf("\n");
            printf("x = "), padic_print(x, ctx), printf("\n");
            printf("y = "), padic_print(y, ctx), printf("\n");
            abort();
        }

        padic_mat_clear(a);
        _padic_clear(x);
        _padic_clear(y);

        fmpz_clear(p);
        padic_ctx_clear(ctx);
    }

    flint_randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}

