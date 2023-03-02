/*
    Copyright (C) 2011, 2012 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <string.h>
#include "ulong_extras.h"
#include "long_extras.h"
#include "padic.h"

int
main(void)
{
    int i, result;
    FLINT_TEST_INIT(state);

    flint_printf("get_str... ");
    fflush(stdout);

    

    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpz_t p;
        slong N;
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
            flint_printf("FAIL:\n\n");
            flint_printf("x = "), padic_print(x, ctx), flint_printf("\n");
            flint_printf("y = "), fmpq_print(y), flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }

        flint_free(s);
        flint_free(t);
        padic_clear(x);
        fmpq_clear(y);
        padic_ctx_clear(ctx);
        fmpz_clear(p);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}

