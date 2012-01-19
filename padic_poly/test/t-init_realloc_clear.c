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
#include <mpir.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_poly.h"
#include "long_extras.h"
#include "ulong_extras.h"
#include "padic_poly.h"

int
main(void)
{
    int i;
    flint_rand_t state;

    printf("init/init2/realloc/clear... ");
    fflush(stdout);

    flint_randinit(state);

    for (i = 0; i < 10000; i++)
    {
        padic_poly_t a;

        padic_poly_init2(a, n_randint(state, 100));
        padic_poly_clear(a);
    }

    for (i = 0; i < 10000; i++)
    {
        fmpz_t p;
        padic_poly_t a;

        fmpz_init(p);
        fmpz_set_ui(p, 7);

        padic_poly_init2(a, n_randint(state, 100));
        padic_poly_realloc(a, n_randint(state, 100), p);
        padic_poly_clear(a);

        fmpz_clear(p);
    }

    for (i = 0; i < 10000; i++)
    {
        padic_ctx_t ctx;
        fmpz_t p;
        long N;

        padic_poly_t a;

        fmpz_init(p);
        fmpz_set_ui(p, n_randprime(state, 5, 1));
        N = z_randint(state, 50);
        padic_ctx_init(ctx, p, N, PADIC_SERIES);

        padic_poly_init(a);
        padic_poly_randtest(a, state, n_randint(state, 100), ctx);
        padic_poly_clear(a);

        padic_ctx_clear(ctx);
        fmpz_clear(p);
    }

    flint_randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}

