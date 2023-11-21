/*
    Copyright (C) 2023 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "ulong_extras.h"
#include "fmpz.h"
#include "gr_vec.h"
#include "gr_poly.h"

TEST_FUNCTION_START(gr_poly_roots_other, state)
{
    slong iter;

    for (iter = 0; iter < 10000; iter++)
    {
        int status;
        gr_ctx_t ctx, ctx_other, ZZ;
        gr_poly_t f, f_other, g, h, q, r;
        gr_vec_t roots, mult;
        slong i, j;

        gr_ctx_init_random(ctx, state);

        if (n_randint(state, 2))
            gr_ctx_init_fmpz(ctx_other);
        else
            gr_ctx_init_random(ctx_other, state);

        gr_ctx_init_fmpz(ZZ);

        gr_poly_init(f, ctx);
        gr_poly_init(f_other, ctx_other);

        gr_vec_init(roots, 0, ctx);
        gr_vec_init(mult, 0, ZZ);

        status = GR_SUCCESS;

        status |= gr_poly_randtest(f_other, state, 1 + n_randint(state, 6), ctx_other);
        status |= gr_poly_roots_other(roots, mult, f_other, ctx_other, 0, ctx);

        if (status == GR_SUCCESS)
        {
            status |= gr_poly_set_gr_poly_other(f, f_other, ctx_other, ctx);

            if (status == GR_SUCCESS)
            {
                gr_poly_init(g, ctx);
                gr_poly_init(h, ctx);
                gr_poly_init(q, ctx);
                gr_poly_init(r, ctx);

                status |= gr_poly_one(g, ctx);

                for (i = 0; i < roots->length; i++)
                {
                    /* (x-r) */
                    status |= gr_poly_set_scalar(h, gr_vec_entry_ptr(roots, i, ctx), ctx);
                    status |= gr_poly_neg(h, h, ctx);
                    status |= gr_poly_set_coeff_si(h, 1, 1, ctx);

                    for (j = 0; j < fmpz_get_si(gr_vec_entry_ptr(mult, i, ZZ)); j++)
                        status |= gr_poly_mul(g, g, h, ctx);
                }

                status |= gr_poly_divrem(q, r, f, g, ctx);

                if (status == GR_SUCCESS && gr_poly_is_zero(r, ctx) == T_FALSE)
                {
                    flint_printf("FAIL: product of factors does not divide polynomial\n\n");
                    flint_printf("f = "); gr_poly_print(f, ctx); flint_printf("\n");
                    flint_printf("g = "); gr_poly_print(g, ctx); flint_printf("\n");
                    flint_abort();
                }

                gr_poly_clear(g, ctx);
                gr_poly_clear(h, ctx);
                gr_poly_clear(q, ctx);
                gr_poly_clear(r, ctx);
            }
        }

        gr_poly_clear(f, ctx);
        gr_poly_clear(f_other, ctx_other);
        gr_vec_clear(roots, ctx);
        gr_vec_clear(mult, ZZ);

        gr_ctx_clear(ctx);
        gr_ctx_clear(ctx_other);
        gr_ctx_clear(ZZ);
    }

    TEST_FUNCTION_END(state);
}
