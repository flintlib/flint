/*
    Copyright (C) 2011 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "ulong_extras.h"
#include "long_extras.h"
#include "padic.h"
#include "padic_mat.h"

int
main(void)
{
    int i, result;

    fmpz_t p;
    slong N;
    padic_ctx_t ctx;
    slong m, n;

    FLINT_TEST_INIT(state);

    flint_printf("get/ set_entry_padic... ");
    fflush(stdout);    

    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        padic_mat_t a;
        padic_t x, y;
        slong r, c;

        fmpz_init_set_ui(p, n_randtest_prime(state, 0));
        N = n_randint(state, PADIC_TEST_PREC_MAX - PADIC_TEST_PREC_MIN) 
            + PADIC_TEST_PREC_MIN;
        padic_ctx_init(ctx, p, FLINT_MAX(0, N-10), FLINT_MAX(0, N+10), PADIC_SERIES);

        m = n_randint(state, 20) + 1;
        n = n_randint(state, 20) + 1;

        padic_mat_init2(a, m, n, N);
        padic_init2(x, N);
        padic_init2(y, N);

        padic_mat_randtest(a, state, ctx);
        padic_randtest_not_zero(x, state, ctx);

        r = n_randint(state, m);
        c = n_randint(state, n);

        padic_mat_set_entry_padic(a, r, c, x, ctx);
        padic_mat_get_entry_padic(y, a, r, c, ctx);

        result = (padic_equal(x, y) && padic_mat_is_reduced(a, ctx));
        if (!result)
        {
            flint_printf("FAIL:\n\n");
            flint_printf("a = "), padic_mat_print(a, ctx), flint_printf("\n");
            flint_printf("x = "), padic_print(x, ctx), flint_printf("\n");
            flint_printf("y = "), padic_print(y, ctx), flint_printf("\n");
            abort();
        }

        padic_mat_clear(a);
        padic_clear(x);
        padic_clear(y);

        fmpz_clear(p);
        padic_ctx_clear(ctx);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}

