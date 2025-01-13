/*
    Copyright (C) 2023 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpq.h"
#include "test_helpers.h"
#include "ulong_extras.h"
#include "gr_poly.h"

int
test_mul_scalar(flint_rand_t state, int which)
{
    gr_ctx_t ctx;
    gr_ptr c;
    ulong c_ui = 0;
    slong c_si = 0;
    fmpz_t c_fmpz;
    fmpq_t c_fmpq;
    gr_poly_t A, B, C;
    int status = GR_SUCCESS;

    gr_ctx_init_random(ctx, state);

    GR_TMP_INIT(c, ctx);

    gr_poly_init(A, ctx);
    gr_poly_init(B, ctx);
    gr_poly_init(C, ctx);

    ulong n = 10;
    GR_MUST_SUCCEED(gr_poly_randtest(A, state, 1 + n_randint(state, n), ctx));

    ulong max_scalar = 1 << 8;
    switch (which)
    {
        case 0:
            c_ui = n_randint(state, max_scalar);
            status |= gr_set_ui(c, c_ui, ctx);
            status |= gr_poly_mul_scalar_ui(C, A, c_ui, ctx);
            break;
        case 1:
            c_si = n_randint(state, max_scalar);
            if (c_si & 1)
                c_si = - c_si;
            status |= gr_set_si(c, c_si, ctx);
            status |= gr_poly_mul_scalar_si(C, A, c_si, ctx);
            break;
        case 2:
            fmpz_init(c_fmpz);
            fmpz_randtest(c_fmpz, state, 8);
            status |= gr_set_fmpz(c,c_fmpz, ctx);
            status |= gr_poly_mul_scalar_fmpz(C, A, c_fmpz, ctx);
            fmpz_clear(c_fmpz);
            break;
        case 3:
            fmpq_init(c_fmpq);
            fmpq_randtest(c_fmpq, state, 8);
            status |= gr_set_fmpq(c, c_fmpq, ctx);
            status |= gr_poly_mul_scalar_fmpq(C, A, c_fmpq, ctx);
            fmpq_clear(c_fmpq);
            break;

        default:
            flint_abort();
    }

    /* Do the same multiplication but with the gr_ptr type scalar */
    status |= gr_poly_mul_scalar(B, A, c, ctx);

    if (status == GR_SUCCESS && gr_poly_equal(B, C, ctx) == T_FALSE)
    {
        flint_printf("FAIL\n\n");
        flint_printf("which = %d\n\n", which);
        gr_ctx_println(ctx);
        flint_printf("A = "); gr_poly_print(A, ctx); flint_printf("\n\n");
        flint_printf("B = "); gr_poly_print(B, ctx); flint_printf("\n\n");
        flint_printf("C = "); gr_poly_print(C, ctx); flint_printf("\n\n");
        flint_printf("c = "); gr_print(c, ctx); flint_printf("\n\n");
        flint_abort();
    }

    GR_TMP_CLEAR(c, ctx);

    gr_poly_clear(A, ctx);
    gr_poly_clear(B, ctx);
    gr_poly_clear(C, ctx);

    gr_ctx_clear(ctx);

    return status;
}

TEST_FUNCTION_START(gr_poly_mul_scalar, state)
{
    slong iter;

    for (iter = 0; iter < 1000; iter++)
    {
        test_mul_scalar(state, n_randint(state, 4));
    }

    TEST_FUNCTION_END(state);
}
