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

TEST_FUNCTION_START(gr_mpoly_quasidivrem_ideal, state)
{
    slong iter;

    for (iter = 0; iter < 1000 * flint_test_multiplier(); iter++)
    {
        gr_ctx_t cctx;
        gr_mpoly_ctx_t ctx;
        gr_mpoly_t A, R, t, u;
        gr_mpoly_vec_t B, Q;
        gr_ptr scale;
        slong w, len, alen, elen;
        flint_bitcnt_t ebits;
        int status;

        /* fmpz exercises the scaling; the others are fields (scale == 1) */
        switch (n_randint(state, 4))
        {
            case 0: gr_ctx_init_fmpz(cctx); break;
            case 1: gr_ctx_init_fmpq(cctx); break;
            case 2: gr_ctx_init_random_finite_field(cctx, state); break;
            default: gr_ctx_init_nmod(cctx, n_randtest_prime(state, 1)); break;
        }

        gr_mpoly_ctx_init_rand(ctx, state, cctx, 8);

        len = 1 + n_randint(state, 4);
        gr_mpoly_init(A, ctx);
        gr_mpoly_init(R, ctx);
        gr_mpoly_init(t, ctx);
        gr_mpoly_init(u, ctx);
        gr_mpoly_vec_init(B, len, ctx);
        gr_mpoly_vec_init(Q, len, ctx);
        GR_TMP_INIT(scale, cctx);

        alen = 1 + n_randint(state, 8);
        elen = 1 + n_randint(state, 5);
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
            switch (n_randint(state, 4))
            {
                case 0:
                    status = gr_mpoly_quasidivrem_ideal(scale, Q, R, A, B, ctx);
                    break;
                case 1:
                    status = gr_mpoly_vec_set(Q, B, ctx);
                    status |= gr_mpoly_quasidivrem_ideal(scale, Q, R, A, Q, ctx);
                    break;
                case 2:
                    status = gr_mpoly_set(R, A, ctx);
                    status |= gr_mpoly_quasidivrem_ideal(scale, Q, R, R, B, ctx);
                    break;
                default:
                    status = gr_mpoly_set(R, A, ctx);
                    status = gr_mpoly_vec_set(Q, B, ctx);
                    status |= gr_mpoly_quasidivrem_ideal(scale, Q, R, R, Q, ctx);
                    break;
            }

            if (status == GR_SUCCESS)
            {
                for (w = 0; w < len; w++)
                    gr_mpoly_assert_canonical(gr_mpoly_vec_entry_ptr(Q, w, ctx), ctx);
                gr_mpoly_assert_canonical(R, ctx);

                /* check scale*A == sum Q[w]*B[w] + R */
                status = gr_mpoly_set(t, R, ctx);
                for (w = 0; w < len; w++)
                {
                    status |= gr_mpoly_mul(u, gr_mpoly_vec_entry_srcptr(Q, w, ctx),
                                              gr_mpoly_vec_entry_srcptr(B, w, ctx), ctx);
                    status |= gr_mpoly_add(t, t, u, ctx);
                }
                status |= gr_mpoly_mul_scalar(u, A, scale, ctx);

                if (status == GR_SUCCESS && gr_mpoly_equal(t, u, ctx) == T_FALSE)
                {
                    flint_printf("FAIL: scale*A != sum Q[w]*B[w] + R\n");
                    gr_ctx_println(cctx);
                    flint_printf("len = %wd\n", len);
                    flint_printf("scale = "); gr_println(scale, cctx);
                    fflush(stdout);
                    flint_abort();
                }
            }
        }

        gr_mpoly_clear(A, ctx);
        gr_mpoly_clear(R, ctx);
        gr_mpoly_clear(t, ctx);
        gr_mpoly_clear(u, ctx);
        gr_mpoly_vec_clear(B, ctx);
        gr_mpoly_vec_clear(Q, ctx);
        GR_TMP_CLEAR(scale, cctx);

        gr_mpoly_ctx_clear(ctx);
        gr_ctx_clear(cctx);
    }

    TEST_FUNCTION_END(state);
}
