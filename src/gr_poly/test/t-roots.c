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

TEST_FUNCTION_START(gr_poly_roots, state)
{
    slong iter;

    for (iter = 0; iter < 1000 * flint_test_multiplier(); iter++)
    {
        int status;
        gr_ctx_t ctx, ZZ;
        gr_poly_t f, g, h, q, r;
        gr_vec_t roots, mult;
        slong i, j, expected_mult;
        int init_from_roots;

        gr_ctx_init_random(ctx, state);
        gr_ctx_init_fmpz(ZZ);

        gr_poly_init(f, ctx);
        gr_vec_init(roots, 0, ctx);
        gr_vec_init(mult, 0, ZZ);

        status = GR_SUCCESS;

        init_from_roots = n_randint(state, 2);
        expected_mult = n_randint(state, 4);

        if (init_from_roots)
        {
            gr_ptr c = gr_heap_init(ctx);

            status |= gr_randtest_not_zero(c, state, ctx);

            if (status == GR_SUCCESS)
            {
                gr_poly_init(h, ctx);

                status |= gr_poly_set_scalar(f, c, ctx);

                for (i = 0; i < expected_mult; i++)
                {
                    status |= gr_poly_set_scalar(h, c, ctx);
                    status |= gr_poly_neg(h, h, ctx);
                    status |= gr_poly_set_coeff_si(h, 1, 1, ctx);
                    status |= gr_poly_mul(f, f, h, ctx);

                    if (n_randint(state, 3) == 0)
                    {
                        status |= gr_poly_mul(f, f, h, ctx);
                        expected_mult++;
                    }
                }

                gr_poly_clear(h, ctx);
            }
            else
            {
                init_from_roots = 0;
            }

            gr_heap_clear(c, ctx);
        }

        if (!init_from_roots)
            status |= gr_poly_randtest(f, state, 1 + n_randint(state, 6), ctx);

        status |= gr_poly_roots(roots, mult, f, 0, ctx);

        if (status == GR_SUCCESS)
        {
            gr_poly_init(g, ctx);
            gr_poly_init(h, ctx);
            gr_poly_init(q, ctx);
            gr_poly_init(r, ctx);

            status |= gr_poly_one(g, ctx);

            if (gr_ctx_is_integral_domain(ctx) == T_TRUE)
            {
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

                if (status == GR_SUCCESS && init_from_roots && g->length != f->length)
                {
                    flint_printf("FAIL: wrong multiplicity for polynomial that splits\n\n");
                    flint_printf("f = "); gr_poly_print(f, ctx); flint_printf("\n");
                    flint_printf("g = "); gr_poly_print(g, ctx); flint_printf("\n");
                    flint_abort();
                }
            }
            else
            {
                for (i = 0; i < roots->length; i++)
                {
                    gr_ptr c = gr_heap_init(ctx);

                    status |= gr_poly_evaluate(c, f, gr_vec_entry_ptr(roots, i, ctx), ctx);

                    if (status == GR_SUCCESS && gr_is_zero(c, ctx) == T_FALSE)
                    {
                        flint_printf("FAIL: f(r) != 0\n\n");
                        flint_printf("f = "); gr_poly_print(f, ctx); flint_printf("\n");
                        flint_printf("r = "); gr_print(gr_vec_entry_ptr(roots, i, ctx), ctx); flint_printf("\n");
                        flint_printf("r = "); gr_print(c, ctx); flint_printf("\n");
                        flint_abort();
                    }

                    gr_heap_clear(c, ctx);
                }
            }

            gr_poly_clear(g, ctx);
            gr_poly_clear(h, ctx);
            gr_poly_clear(q, ctx);
            gr_poly_clear(r, ctx);
        }

        gr_poly_clear(f, ctx);
        gr_vec_clear(roots, ctx);
        gr_vec_clear(mult, ZZ);

        gr_ctx_clear(ctx);
        gr_ctx_clear(ZZ);
    }

    TEST_FUNCTION_END(state);
}
