/*
    Copyright (C) 2025 Ricardo Buring

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "gr_mat.h"
#include "gr_poly.h"

FLINT_DLL extern gr_static_method_table _ca_methods;

TEST_GR_FUNCTION_START(gr_mat_gr_poly_solve_lode_newton, state, count_success, count_unable, count_domain)
{
    slong iter;

    for (iter = 0; iter < 1000; iter++)
    {
        gr_ctx_t ctx, poly_ctx;
        gr_mat_t A_numerator, Y, Y0, AY, err, tmp_mat;
        gr_poly_t A_denominator, tmp_poly;
        slong n, den_len, sol_len, i, j;
        int ic_satisfied, ode_satisfied, status = GR_SUCCESS;

        n = n_randint(state, 6);
        den_len = 1 + n_randint(state, 6);
        sol_len = 1 + n_randint(state, 16);

        gr_ctx_init_random(ctx, state);
        /* Hack: avoid CA because slow */
        while (ctx->methods == _ca_methods || gr_ctx_is_finite_characteristic(ctx) != T_FALSE || gr_ctx_is_field(ctx) != T_TRUE)
        {
            gr_ctx_clear(ctx);
            gr_ctx_init_random(ctx, state);
        }
        gr_ctx_init_gr_poly(poly_ctx, ctx);

        gr_mat_init(A_numerator, n, n, poly_ctx);
        gr_poly_init(A_denominator, ctx);
        gr_mat_init(Y, n, n, poly_ctx);
        gr_mat_init(Y0, n, n, ctx);
        gr_mat_init(AY, n, n, poly_ctx);
        gr_mat_init(err, n, n, poly_ctx);
        gr_mat_init(tmp_mat, n, n, poly_ctx);

        if (iter % 2 == 0)
        {
            status |= gr_mat_randtest(A_numerator, state, poly_ctx);
            status |= gr_poly_randtest(A_denominator, state, den_len, ctx);
            while (gr_poly_length(A_denominator, ctx) == 0 || gr_is_invertible(gr_poly_coeff_srcptr(A_denominator, 0, ctx), ctx) != T_TRUE)
            {
                status |= gr_poly_randtest(A_denominator, state, den_len, ctx);
            }
        }
        else
        {
            /* TODO: Make tmp_poly a gr_ore_poly_t instead. */
            gr_poly_init(tmp_poly, poly_ctx);
            gr_poly_randtest(tmp_poly, state, n + 1, poly_ctx);
            while (gr_poly_length(tmp_poly, poly_ctx) != n + 1
                     || gr_poly_length(gr_poly_coeff_srcptr(tmp_poly, n, poly_ctx), ctx) == 0
                     || gr_is_invertible(gr_poly_coeff_srcptr(gr_poly_coeff_srcptr(tmp_poly, n, poly_ctx), 0, ctx), ctx) != T_TRUE)
            {
                status |= gr_poly_randtest(tmp_poly, state, n + 1, poly_ctx);
            }
            status |= gr_mat_companion_fraction(A_numerator, A_denominator, tmp_poly, poly_ctx);
            gr_poly_clear(tmp_poly, poly_ctx);
        }

        /* status |= gr_mat_one(Y0, ctx); */
        status |= gr_mat_randrank(Y0, state, n, ctx);

        status |= gr_mat_gr_poly_solve_lode_newton(Y, A_numerator, A_denominator, Y0, sol_len, poly_ctx, poly_ctx);

        if (status == GR_SUCCESS)
        {
            /* Check initial condition is satisfied. */
            for (i = 0; i < n; i++)
            {
                for (j = 0; j < n; j++)
                {
                    gr_ptr entry_tmp = gr_mat_entry_ptr(tmp_mat, i, j, poly_ctx);
                    status |= gr_poly_truncate(entry_tmp, gr_mat_entry_ptr(Y, i, j, poly_ctx), 1, ctx);
                    status |= gr_poly_sub_scalar(entry_tmp, entry_tmp, gr_mat_entry_ptr(Y0, i, j, ctx), ctx);
                }
            }
            ic_satisfied = gr_mat_is_zero(tmp_mat, poly_ctx) != T_FALSE;

            /* Check ODE is satisfied. */
            status |= gr_mat_mul(AY, A_numerator, Y, poly_ctx);
            for (i = 0; i < n; i++)
            {
                for (j = 0; j < n; j++)
                {
                    gr_ptr entry_err = gr_mat_entry_ptr(err, i, j, poly_ctx);
                    gr_ptr entry_Y = gr_mat_entry_ptr(Y, i, j, poly_ctx);
                    gr_ptr entry_AY = gr_mat_entry_ptr(AY, i, j, poly_ctx);
                    status |= gr_poly_derivative(entry_err, entry_Y, ctx);
                    status |= gr_poly_mullow(entry_err, entry_err, A_denominator, sol_len - 1, ctx);
                    status |= gr_poly_truncate(entry_AY, entry_AY, sol_len - 1, ctx);
                    status |= gr_poly_sub(entry_err, entry_err, entry_AY, ctx);
                }
            }
            ode_satisfied = gr_mat_is_zero(err, poly_ctx) != T_FALSE;

            if (status == GR_SUCCESS && (!ic_satisfied || !ode_satisfied))
            {
                flint_printf("FAIL\n");
                gr_ctx_println(ctx);
                flint_printf("ic_satisfied = %d, ode_satisfied = %d\n", ic_satisfied, ode_satisfied);
                flint_printf("A_numerator = \n"); gr_mat_print(A_numerator, poly_ctx); flint_printf("\n\n");
                flint_printf("A_denominator = \n"); gr_poly_print(A_denominator, ctx); flint_printf("\n\n");
                flint_printf("Y0 = \n"); gr_mat_print(Y0, ctx); flint_printf("\n\n");
                flint_printf("sol_len = %ld\n", sol_len);
                flint_printf("Y = \n"); gr_mat_print(Y, poly_ctx); flint_printf("\n\n");
                flint_printf("AY = \n"); gr_mat_print(AY, poly_ctx); flint_printf("\n\n");
                flint_printf("err = \n"); gr_mat_print(err, poly_ctx); flint_printf("\n\n");
                flint_abort();
            }
        }

        count_success += (status == GR_SUCCESS);
        count_domain += ((status & GR_DOMAIN) != 0);
        count_unable += ((status & GR_UNABLE) != 0);

        gr_mat_clear(err, poly_ctx);
        gr_mat_clear(AY, poly_ctx);
        gr_mat_clear(Y0, ctx);
        gr_mat_clear(Y, poly_ctx);
        gr_poly_clear(A_denominator, ctx);
        gr_mat_clear(A_numerator, poly_ctx);
        gr_ctx_clear(poly_ctx);
        gr_ctx_clear(ctx);
    }

    TEST_GR_FUNCTION_END(state, count_success, count_unable, count_domain);
}
