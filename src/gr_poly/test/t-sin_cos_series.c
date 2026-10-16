/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "ulong_extras.h"
#include "gr_poly.h"

FLINT_DLL extern gr_static_method_table _ca_methods;

TEST_FUNCTION_START(gr_poly_sin_cos_series, state)
{
    slong iter;

    for (iter = 0; iter < 100 * flint_test_multiplier(); iter++)
    {
        gr_ctx_t ctx;
        slong n;
        gr_poly_t a, Sref, Cref, S, C;
        int status = GR_SUCCESS;
        int times_pi, aliasing, function;

        gr_ctx_init_random_commutative_ring(ctx, state);

        gr_poly_init(a, ctx);
        gr_poly_init(Sref, ctx);
        gr_poly_init(Cref, ctx);
        gr_poly_init(S, ctx);
        gr_poly_init(C, ctx);

        if (ctx->methods == _ca_methods)
            n = 1 + n_randint(state, 4);
        else if (gr_ctx_is_finite(ctx) == T_TRUE)
            n = 1 + n_randint(state, 100);
        else
            n = 1 + n_randint(state, 30);

        GR_MUST_SUCCEED(gr_poly_randtest(a, state, n + 2, ctx));
        GR_MUST_SUCCEED(gr_poly_randtest(S, state, n + 2, ctx));
        GR_MUST_SUCCEED(gr_poly_randtest(C, state, n + 2, ctx));

        if (n_randint(state, 4) != 0)
            status |= gr_poly_set_coeff_si(a, 0, 0, ctx);

        times_pi = n_randint(state, 2);
        function = n_randint(state, 4);    

        status |= gr_poly_sin_cos_series_basecase(Sref, Cref, a, n + n_randint(state, 2), times_pi, ctx);
        status |= gr_poly_truncate(Sref, Sref, n, ctx);
        status |= gr_poly_truncate(Cref, Cref, n, ctx);

        aliasing = n_randint(state, 3);

        if (function == 0)
        {
            if (aliasing == 0)
            {
                status |= gr_poly_set(S, a, ctx);
                status |= gr_poly_sin_cos_series_newton(S, C, S, n, n_randint(state, 20), times_pi, ctx);
            }
            else if (aliasing == 1)
            {
                status |= gr_poly_set(C, a, ctx);
                status |= gr_poly_sin_cos_series_newton(S, C, C, n, n_randint(state, 20), times_pi, ctx);
            }
            else
            {
                status |= gr_poly_sin_cos_series_newton(S, C, a, n, n_randint(state, 20), times_pi, ctx);
            }
        }
        else if (function == 1)
        {
            if (aliasing == 0)
            {
                status |= gr_poly_set(S, a, ctx);
                status |= gr_poly_sin_cos_series_tangent(S, C, S, n, times_pi, ctx);
            }
            else if (aliasing == 1)
            {
                status |= gr_poly_set(C, a, ctx);
                status |= gr_poly_sin_cos_series_tangent(S, C, C, n, times_pi, ctx);
            }
            else
            {
                status |= gr_poly_sin_cos_series_tangent(S, C, a, n, times_pi, ctx);
            }
        }
        else if (function == 2)
        {
            if (aliasing == 0)
            {
                status |= gr_poly_set(S, a, ctx);
                if (times_pi)
                    status |= gr_poly_sin_cos_pi_series(S, C, S, n, ctx);
                else
                    status |= gr_poly_sin_cos_series(S, C, S, n, ctx);
            }
            else if (aliasing == 1)
            {
                status |= gr_poly_set(C, a, ctx);
                if (times_pi)
                    status |= gr_poly_sin_cos_pi_series(S, C, C, n, ctx);
                else
                    status |= gr_poly_sin_cos_series(S, C, C, n, ctx);
            }
            else
            {
                if (times_pi)
                    status |= gr_poly_sin_cos_pi_series(S, C, a, n, ctx);
                else
                    status |= gr_poly_sin_cos_series(S, C, a, n, ctx);
            }
        }
        else
        {
            if (aliasing == 0)
            {
                status |= gr_poly_set(S, a, ctx);
                status |= gr_poly_set(C, a, ctx);
                if (times_pi)
                {
                    status |= gr_poly_sin_pi_series(S, S, n, ctx);
                    status |= gr_poly_cos_pi_series(C, C, n, ctx);
                }
                else
                {
                    status |= gr_poly_sin_series(S, S, n, ctx);
                    status |= gr_poly_cos_series(C, C, n, ctx);
                }
            }
            else
            {
                if (times_pi)
                {
                    status |= gr_poly_sin_pi_series(S, a, n, ctx);
                    status |= gr_poly_cos_pi_series(C, a, n, ctx);
                }
                else
                {
                    status |= gr_poly_sin_series(S, a, n, ctx);
                    status |= gr_poly_cos_series(C, a, n, ctx);
                }
            }
        }

        if (status == GR_SUCCESS && (gr_poly_equal(S, Sref, ctx) == T_FALSE ||
                gr_poly_equal(C, Cref, ctx) == T_FALSE))
        {
            flint_printf("FAIL: gr_poly_sin_cos_series\n\n");
            gr_ctx_println(ctx);
            flint_printf("times_pi = %d, function = %d, aliasing = %d\n", times_pi, function, aliasing);
            flint_printf("n = %wd\n", n);
            flint_printf("a = "); gr_poly_print(a, ctx); flint_printf("\n\n");
            flint_printf("Sref = "); gr_poly_print(Sref, ctx); flint_printf("\n\n");
            flint_printf("S    = "); gr_poly_print(S, ctx); flint_printf("\n\n");
            flint_printf("Cref = "); gr_poly_print(Cref, ctx); flint_printf("\n\n");
            flint_printf("C    = "); gr_poly_print(C, ctx); flint_printf("\n\n");
            flint_abort();
        }

        gr_poly_clear(a, ctx);
        gr_poly_clear(Sref, ctx);
        gr_poly_clear(Cref, ctx);
        gr_poly_clear(S, ctx);
        gr_poly_clear(C, ctx);

        gr_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}

