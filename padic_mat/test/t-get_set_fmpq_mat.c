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

    Copyright (C) 2011, 2012 Sebastian Pancratz

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

    printf("get/ set_fmpq_mat... ");
    fflush(stdout);

    flint_randinit(state);

    /* Qp -> QQ -> Qp */
    for (i = 0; i < 1000; i++)
    {
        fmpz_t p;
        long N;
        padic_ctx_t ctx;
        long m, n;

        padic_mat_t a, c;
        fmpq_mat_t b;

        fmpz_init(p);
        fmpz_set_ui(p, n_randprime(state, 5, 1));
        N = z_randint(state, 50);
        padic_ctx_init(ctx, p, N, PADIC_SERIES);

        m = n_randint(state, 10);
        n = n_randint(state, 10);

        padic_mat_init(a, m, n);
        padic_mat_init(c, m, n);
        fmpq_mat_init(b, m, n);

        padic_mat_randtest(a, state, ctx);
        padic_mat_get_fmpq_mat(b, a, ctx);
        padic_mat_set_fmpq_mat(c, b, ctx);

        result = (padic_mat_equal(a, c) && _padic_mat_is_canonical(a, p));

        if (!result)
        {
            printf("FAIL:\n\n");
            printf("a = "), padic_mat_print(a, ctx), printf("\n");
            printf("c = "), padic_mat_print(c, ctx), printf("\n");
            printf("b = "), fmpq_mat_print(b), printf("\n");
            abort();
        }

        padic_mat_clear(a);
        padic_mat_clear(c);
        fmpq_mat_clear(b);

        fmpz_clear(p);
        padic_ctx_clear(ctx);
    }

    flint_randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}

