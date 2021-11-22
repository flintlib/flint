/*
    Copyright (C) 2011, 2012, 2013 Sebastian Pancratz

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

    flint_printf("get/ set_fmpq_mat... ");
    fflush(stdout);

    /* Qp -> QQ -> Qp */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        padic_mat_t a, c;
        fmpq_mat_t b;

        fmpz_init_set_ui(p, n_randtest_prime(state, 0));
        N = n_randint(state, PADIC_TEST_PREC_MAX - PADIC_TEST_PREC_MIN) 
            + PADIC_TEST_PREC_MIN;
        padic_ctx_init(ctx, p, FLINT_MAX(0, N-10), FLINT_MAX(0, N+10), PADIC_SERIES);

        m = n_randint(state, 10);
        n = n_randint(state, 10);

        padic_mat_init2(a, m, n, N);
        padic_mat_init2(c, m, n, N);
        fmpq_mat_init(b, m, n);

        padic_mat_randtest(a, state, ctx);
        padic_mat_get_fmpq_mat(b, a, ctx);
        padic_mat_set_fmpq_mat(c, b, ctx);

        result = (padic_mat_equal(a, c) && padic_mat_is_reduced(a, ctx));

        if (!result)
        {
            flint_printf("FAIL:\n\n");
            flint_printf("a = "), padic_mat_print(a, ctx), flint_printf("\n");
            flint_printf("c = "), padic_mat_print(c, ctx), flint_printf("\n");
            flint_printf("b = "), fmpq_mat_print(b), flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }

        padic_mat_clear(a);
        padic_mat_clear(c);
        fmpq_mat_clear(b);

        fmpz_clear(p);
        padic_ctx_clear(ctx);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}

