/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "gr_mpoly.h"

/* verify A == sum Q[w]*B[w] + R */
static void
_check_identity(const char * name, const gr_mpoly_vec_t Q, const gr_mpoly_vec_t B,
    const gr_mpoly_t R, const gr_mpoly_t A, gr_ctx_t cctx, gr_mpoly_ctx_t ctx)
{
    gr_mpoly_t t, u;
    slong w, len = B->length;
    int status = GR_SUCCESS;

    gr_mpoly_init(t, ctx);
    gr_mpoly_init(u, ctx);

    status |= gr_mpoly_set(t, R, ctx);
    for (w = 0; w < len; w++)
    {
        status |= gr_mpoly_mul(u, gr_mpoly_vec_entry_srcptr(Q, w, ctx),
                                  gr_mpoly_vec_entry_srcptr(B, w, ctx), ctx);
        status |= gr_mpoly_add(t, t, u, ctx);
    }

    if (status == GR_SUCCESS && gr_mpoly_equal(t, A, ctx) == T_FALSE)
    {
        flint_printf("FAIL: %s: A != sum Q[w]*B[w] + R\n", name);
        gr_ctx_println(cctx);
        flint_printf("len = %wd\n", len);
        fflush(stdout);
        flint_abort();
    }

    gr_mpoly_clear(t, ctx);
    gr_mpoly_clear(u, ctx);
}

TEST_FUNCTION_START(gr_mpoly_divrem_ideal, state)
{
    slong iter;

    for (iter = 0; iter < 1000 * flint_test_multiplier(); iter++)
    {
        gr_ctx_t cctx;
        gr_mpoly_ctx_t ctx;
        gr_mpoly_t A, R;
        gr_mpoly_vec_t B, Q;
        slong w, len, alen, elen;
        flint_bitcnt_t ebits;
        int status, big;

        big = (n_randint(state, 4) == 0);

        gr_ctx_t ctx1;
        gr_ctx_init_nmod8(ctx1, n_randprime(state, 8, 1));

        switch (n_randint(state, 5))
        {
            case 0: gr_ctx_init_fmpz(cctx); break;
            case 1: gr_ctx_init_fmpq(cctx); break;
            case 2: gr_ctx_init_random_finite_field(cctx, state); break;
            case 3: gr_ctx_init_nmod(cctx, n_randtest_prime(state, 1)); break;
            default: gr_ctx_init_debug(cctx, ctx1, n_randint(state, 2), n_randint(state, 2) ? 0.0 : 0.001);
        }

        gr_mpoly_ctx_init_rand(ctx, state, cctx, big ? 3 : 10);

        len = 1 + n_randint(state, 4);
        gr_mpoly_init(A, ctx);
        gr_mpoly_init(R, ctx);
        gr_mpoly_vec_init(B, len, ctx);
        gr_mpoly_vec_init(Q, len, ctx);

        alen = 1 + n_randint(state, 10);
        elen = 1 + n_randint(state, 6);
        ebits = 2 + n_randint(state, 6);

        status = GR_SUCCESS;
        status |= gr_mpoly_randtest_bits(A, state, alen, ebits, ctx);

        for (w = 0; w < len; w++)
        {
            gr_mpoly_struct * Bw = gr_mpoly_vec_entry_ptr(B, w, ctx);
            status |= gr_mpoly_randtest_bits(Bw, state, elen, ebits, ctx);
            if (gr_mpoly_is_zero(Bw, ctx) != T_FALSE)
                status |= gr_mpoly_one(Bw, ctx);
        }

        if (status == GR_SUCCESS)
        {
            status = gr_mpoly_divrem_ideal(Q, R, A, B, ctx);

            /* field-like */
            switch (n_randint(state, 4) * 0)
            {
                case 0:
                    status = gr_mpoly_divrem_ideal(Q, R, A, B, ctx);
                    break;
                case 1:
                    status = gr_mpoly_vec_set(Q, B, ctx);
                    status |= gr_mpoly_divrem_ideal(Q, R, A, Q, ctx);
                    break;
                case 2:
                    status = gr_mpoly_set(R, A, ctx);
                    status |= gr_mpoly_divrem_ideal(Q, R, R, B, ctx);
                    break;
                default:
                    status = gr_mpoly_set(R, A, ctx);
                    status = gr_mpoly_vec_set(Q, B, ctx);
                    status |= gr_mpoly_divrem_ideal(Q, R, R, Q, ctx);
                    break;
            }

            if (status == GR_SUCCESS)
            {
                for (w = 0; w < len; w++)
                    gr_mpoly_assert_canonical(gr_mpoly_vec_entry_ptr(Q, w, ctx), ctx);
                gr_mpoly_assert_canonical(R, ctx);
                _check_identity("field", Q, B, R, A, cctx, ctx);
            }

            /* allowing non-unit leading coefficients (euclidean) */
            status = gr_mpoly_divrem_ideal_weak(Q, R, A, B, ctx);
            if (status == GR_SUCCESS)
            {
                for (w = 0; w < len; w++)
                    gr_mpoly_assert_canonical(gr_mpoly_vec_entry_ptr(Q, w, ctx), ctx);
                gr_mpoly_assert_canonical(R, ctx);
                _check_identity("weak", Q, B, R, A, cctx, ctx);
            }
        }

        gr_mpoly_clear(A, ctx);
        gr_mpoly_clear(R, ctx);
        gr_mpoly_vec_clear(B, ctx);
        gr_mpoly_vec_clear(Q, ctx);

        gr_mpoly_ctx_clear(ctx);
        gr_ctx_clear(cctx);
        gr_ctx_clear(ctx1);
    }

    TEST_FUNCTION_END(state);
}
