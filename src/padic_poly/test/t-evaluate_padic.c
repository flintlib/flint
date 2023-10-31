/*
    Copyright (C) 2011, 2012 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "ulong_extras.h"
#include "fmpq.h"
#include "fmpq_poly.h"
#include "padic_poly.h"

TEST_FUNCTION_START(padic_poly_evaluate_padic, state)
{
    int i, result;

    padic_ctx_t ctx;
    fmpz_t p;
    slong N;

    /* Compare with the computation over QQ */
    for (i = 0; i < 20 * flint_test_multiplier(); i++)
    {
        padic_poly_t f;
        fmpq_poly_t fQQ;
        padic_t a, y, z;
        fmpq_t aQQ, yQQ;

        fmpz_init_set_ui(p, n_randtest_prime(state, 0));
        N = n_randint(state, PADIC_TEST_PREC_MAX - PADIC_TEST_PREC_MIN)
            + PADIC_TEST_PREC_MIN;
        padic_ctx_init(ctx, p, FLINT_MAX(0, N-10), FLINT_MAX(0, N+10), PADIC_SERIES);

        padic_poly_init2(f, 0, N);
        fmpq_poly_init(fQQ);
        padic_init2(a, N);
        padic_init2(y, N);
        padic_init2(z, N);
        fmpq_init(aQQ);
        fmpq_init(yQQ);

        padic_poly_randtest(f, state, n_randint(state, 80), ctx);
        padic_randtest(a, state, ctx);

        padic_poly_get_fmpq_poly(fQQ, f, ctx);
        padic_get_fmpq(aQQ, a, ctx);

        padic_poly_evaluate_padic(y, f, a, ctx);
        fmpq_poly_evaluate_fmpq(yQQ, fQQ, aQQ);

        padic_set_fmpq(z, yQQ, ctx);

        if (padic_val(a) >= 0)
        {
            result = (padic_equal(y, z));
            if (!result)
            {
                flint_printf("FAIL (cmp with QQ):\n");
                flint_printf("f = "), padic_poly_print(f, ctx), flint_printf("\n\n");
                flint_printf("a = "), padic_print(a, ctx), flint_printf("\n\n");
                flint_printf("y = "), padic_print(y, ctx), flint_printf("\n\n");
                flint_printf("z = "), padic_print(z, ctx), flint_printf("\n\n");
                fflush(stdout);
                flint_abort();
            }
        }
        else
        {
            slong N2 = N + (f->length - 1) * padic_val(a);
            padic_t y2, z2;

            padic_init2(y2, N2);
            padic_init2(z2, N2);

            padic_set(y2, y, ctx);
            padic_set(z2, z, ctx);

            result = (padic_equal(y2, z2));
            if (!result)
            {
                flint_printf("FAIL (cmp with QQ):\n");
                flint_printf("f  = "), padic_poly_print(f, ctx), flint_printf("\n\n");
                flint_printf("a  = "), padic_print(a,  ctx), flint_printf("\n\n");
                flint_printf("y  = "), padic_print(y,  ctx), flint_printf("\n\n");
                flint_printf("z  = "), padic_print(z,  ctx), flint_printf("\n\n");
                flint_printf("y2 = "), padic_print(y2, ctx), flint_printf("\n\n");
                flint_printf("z2 = "), padic_print(z2, ctx), flint_printf("\n\n");
                fflush(stdout);
                flint_abort();
            }

            padic_clear(y2);
            padic_clear(z2);
        }

        padic_poly_clear(f);
        fmpq_poly_clear(fQQ);
        padic_clear(a);
        padic_clear(y);
        padic_clear(z);
        fmpq_clear(aQQ);
        fmpq_clear(yQQ);

        padic_ctx_clear(ctx);
        fmpz_clear(p);
    }

    TEST_FUNCTION_END(state);
}
