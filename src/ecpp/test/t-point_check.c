/*
    Copyright (C) 2026 Fredrik Johansson
    Copyright (C) 2026 Maël Hostettler
    Developed using Claude Fable 5.1 and Opus 5.5

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz.h"
#include "fmpz_mod.h"
#include "gr.h"
#include "gr_ec.h"
#include "ecpp.h"

/* compare with k P and q (k P) computed directly with gr_ec */
TEST_FUNCTION_START(ecpp_point_check, state)
{
    slong iter;

    for (iter = 0; iter < 200 * flint_test_multiplier(); iter++)
    {
        fmpz_t n, a, b, x, y, k, q, t;
        fmpz_mod_ctx_t ctx;
        gr_ctx_t R;
        gr_ec_ctx_t E;
        gr_ec_jac_point_t P, A, B;
        gr_ptr u, v;
        int expected, result;

        fmpz_init(n); fmpz_init(a); fmpz_init(b); fmpz_init(x);
        fmpz_init(y); fmpz_init(k); fmpz_init(q); fmpz_init(t);

        fmpz_randprime(n, state, 10 + n_randint(state, 200), 0);
        fmpz_mod_ctx_init(ctx, n);

        /* a random curve through a random point, with j = 0 often */
        if (n_randint(state, 4) == 0)
            fmpz_zero(a);
        else
            fmpz_randm(a, state, n);
        fmpz_randm(x, state, n);
        fmpz_randm(y, state, n);
        fmpz_mod_mul(b, y, y, ctx);
        fmpz_mod_mul(t, x, x, ctx);
        fmpz_mod_add(t, t, a, ctx);
        fmpz_mod_mul(t, t, x, ctx);
        fmpz_mod_sub(b, b, t, ctx);

        /* small multipliers so that both outcomes occur */
        fmpz_randtest_unsigned(k, state, 1 + n_randint(state, 60));
        fmpz_randtest_unsigned(q, state, 1 + n_randint(state, 60));

        ecpp_gr_ctx_init(R, ctx);
        GR_MUST_SUCCEED(gr_ctx_set_is_field(R, T_TRUE));
        u = gr_heap_init(R);
        v = gr_heap_init(R);
        GR_MUST_SUCCEED(gr_set_fmpz(u, a, R));
        GR_MUST_SUCCEED(gr_set_fmpz(v, b, R));

        if (gr_ec_ctx_init_short_weierstrass(E, R, u, v) == GR_SUCCESS)
        {
            gr_ec_jac_point_init(P, E);
            gr_ec_jac_point_init(A, E);
            gr_ec_jac_point_init(B, E);

            GR_MUST_SUCCEED(gr_set_fmpz(u, x, R));
            GR_MUST_SUCCEED(gr_set_fmpz(v, y, R));
            GR_MUST_SUCCEED(gr_ec_jac_point_set_affine(P, u, v, E));
            GR_MUST_SUCCEED(gr_ec_jac_point_mul_fmpz(A, P, k, E));
            GR_MUST_SUCCEED(gr_ec_jac_point_mul_fmpz(B, A, q, E));

            if (gr_ec_jac_point_is_inf(A, E) == T_TRUE)
                expected = ECPP_POINT_ZERO;
            else if (gr_ec_jac_point_is_inf(B, E) == T_TRUE)
                expected = ECPP_POINT_OK;
            else
                expected = ECPP_POINT_WRONG_ORDER;

            result = ecpp_point_check(x, y, a, b, k, q, ctx);

            if (result != expected)
            {
                flint_printf("FAIL: got %d, expected %d\n", result, expected);
                flint_printf("n = "); fmpz_print(n);
                flint_printf("\na = "); fmpz_print(a);
                flint_printf("\nb = "); fmpz_print(b);
                flint_printf("\nk = "); fmpz_print(k);
                flint_printf("\nq = "); fmpz_print(q);
                flint_printf("\n");
                flint_abort();
            }

            gr_ec_jac_point_clear(P, E);
            gr_ec_jac_point_clear(A, E);
            gr_ec_jac_point_clear(B, E);
            gr_ec_ctx_clear(E);
        }

        gr_heap_clear(u, R);
        gr_heap_clear(v, R);
        gr_ctx_clear(R);
        fmpz_mod_ctx_clear(ctx);

        fmpz_clear(n); fmpz_clear(a); fmpz_clear(b); fmpz_clear(x);
        fmpz_clear(y); fmpz_clear(k); fmpz_clear(q); fmpz_clear(t);
    }

    TEST_FUNCTION_END(state);
}
