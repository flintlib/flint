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

#include <string.h>
#include "ulong_extras.h"
#include "long_extras.h"
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
        len_t N;
        padic_ctx_t ctx;

        padic_t x;
        fmpq_t y;

        char *s, *t;

        fmpz_init_set_ui(p, n_randtest_prime(state, 0));
        N = n_randint(state, PADIC_TEST_PREC_MAX - PADIC_TEST_PREC_MIN) 
            + PADIC_TEST_PREC_MIN;
        padic_ctx_init(ctx, p, FLINT_MAX(0, N-10), FLINT_MAX(0, N+10), PADIC_TERSE);

        padic_init2(x, N);
        fmpq_init(y);

        padic_randtest(x, state, ctx);
        padic_get_fmpq(y, x, ctx);

        s = padic_get_str(NULL, x, ctx);
        t = fmpq_get_str(NULL, 10, y);

        result = strcmp(s, t) == 0;
        if (!result)
        {
            printf("FAIL:\n\n");
            printf("x = "), padic_print(x, ctx), printf("\n");
            printf("y = "), fmpq_print(y), printf("\n");
            abort();
        }

        free(s);
        free(t);
        padic_clear(x);
        fmpq_clear(y);
        padic_ctx_clear(ctx);
        fmpz_clear(p);
    }

    flint_randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}

