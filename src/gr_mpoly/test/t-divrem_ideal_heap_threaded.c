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

TEST_FUNCTION_START(gr_mpoly_divrem_ideal_heap_threaded, state)
{
    slong iter;

    for (iter = 0; iter < 1000 * flint_test_multiplier(); iter++)
    {
        gr_ctx_t cctx;
        gr_mpoly_ctx_t ctx;
        gr_mpoly_t A, R1, R2, recon, t;
        gr_mpoly_vec_t B, Q1, Q2;
        slong lenA, lenB, len, w;
        flint_bitcnt_t exp_bits;
        int status, s1, s2, nonfield;
        int big = (n_randint(state, 8) == 0);

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

        gr_mpoly_ctx_init_rand(ctx, state, cctx, 12);

        gr_mpoly_init(A, ctx);
        gr_mpoly_init(R1, ctx);
        gr_mpoly_init(R2, ctx);
        gr_mpoly_init(recon, ctx);
        gr_mpoly_init(t, ctx);
        gr_mpoly_vec_init(B, 0, ctx);
        gr_mpoly_vec_init(Q1, 0, ctx);
        gr_mpoly_vec_init(Q2, 0, ctx);

        flint_set_num_threads(1 + n_randint(state, 4));

        len = 2 + n_randint(state, 3);
        gr_mpoly_vec_set_length(B, len, ctx);

        /*
            Keep the exponent bit width modest even in the "big" case: a
            sparse dividend with a very wide exponent range reduced by a
            divisor whose leading term has small degree in some variable
            can legitimately need a huge number of quotient terms (the
            multivariate analogue of dividing x^(10^6) by (x - 1)) -- this
            is a real worst case of the algorithm, not a bug, but it is not
            something a routine regression test should try to exercise.
        */
        lenA = n_randint(state, big ? 150 : 20) + 1;
        lenB = n_randint(state, big ? 12 : 8) + 3;
        exp_bits = n_randint(state, big ? 5 : 3) + 2;
        nonfield = n_randint(state, 2);

        status = GR_SUCCESS;
        status |= gr_mpoly_randtest_bits(A, state, lenA, exp_bits, ctx);
        for (w = 0; w < len; w++)
        {
            status |= gr_mpoly_randtest_bits(gr_mpoly_vec_entry_ptr(B, w, ctx),
                                              state, lenB, exp_bits, ctx);
            if (gr_mpoly_is_zero(gr_mpoly_vec_entry_ptr(B, w, ctx), ctx) != T_FALSE)
                status |= gr_mpoly_one(gr_mpoly_vec_entry_ptr(B, w, ctx), ctx);
        }

        if (status == GR_SUCCESS)
        {
            s1 = nonfield ? gr_mpoly_divrem_ideal_weak(Q1, R1, A, B, ctx)
                          : gr_mpoly_divrem_ideal(Q1, R1, A, B, ctx);

            /* exercise aliasing of the threaded call's own Q2/R2 with its
               own A/B arguments (Q2 == B, or R2 == A); B and A themselves
               are never mutated here (Q2/R2 receive copies beforehand),
               so the later comparisons against Q1/R1/B/A below remain
               valid regardless of which case is taken */
            switch (n_randint(state, 4))
            {
                case 0:
                    s2 = nonfield ? gr_mpoly_divrem_ideal_weak_heap_threaded(Q2, R2, A, B, ctx)
                                  : gr_mpoly_divrem_ideal_heap_threaded(Q2, R2, A, B, ctx);
                    break;
                case 1:
                    status |= gr_mpoly_vec_set(Q2, B, ctx);
                    s2 = nonfield ? gr_mpoly_divrem_ideal_weak_heap_threaded(Q2, R2, A, Q2, ctx)
                                  : gr_mpoly_divrem_ideal_heap_threaded(Q2, R2, A, Q2, ctx);
                    break;
                case 2:
                    status |= gr_mpoly_set(R2, A, ctx);
                    s2 = nonfield ? gr_mpoly_divrem_ideal_weak_heap_threaded(Q2, R2, R2, B, ctx)
                                  : gr_mpoly_divrem_ideal_heap_threaded(Q2, R2, R2, B, ctx);
                    break;
                default:
                    status |= gr_mpoly_set(R2, A, ctx);
                    status |= gr_mpoly_vec_set(Q2, B, ctx);
                    s2 = nonfield ? gr_mpoly_divrem_ideal_weak_heap_threaded(Q2, R2, R2, Q2, ctx)
                                  : gr_mpoly_divrem_ideal_heap_threaded(Q2, R2, R2, Q2, ctx);
                    break;
            }

            if (cctx->which_ring != GR_CTX_DEBUG && s1 != s2)
            {
                flint_printf("FAIL: status mismatch serial=%d threaded=%d\n", s1, s2);
                flint_printf("iter = %wd\n", iter);
                gr_ctx_println(cctx);
                fflush(stdout);
                flint_abort();
            }

            if (s1 == GR_SUCCESS && s2 == GR_SUCCESS)
            {
                gr_mpoly_assert_canonical(R1, ctx);
                gr_mpoly_assert_canonical(R2, ctx);

                if (gr_mpoly_equal(R1, R2, ctx) == T_FALSE)
                {
                    flint_printf("FAIL: remainder mismatch\n");
                    flint_printf("iter = %wd\n", iter);
                    gr_ctx_println(cctx);
                    fflush(stdout);
                    flint_abort();
                }

                for (w = 0; w < len; w++)
                {
                    gr_mpoly_struct * q1w = gr_mpoly_vec_entry_ptr(Q1, w, ctx);
                    gr_mpoly_struct * q2w = gr_mpoly_vec_entry_ptr(Q2, w, ctx);
                    gr_mpoly_assert_canonical(q1w, ctx);
                    gr_mpoly_assert_canonical(q2w, ctx);
                    if (gr_mpoly_equal(q1w, q2w, ctx) == T_FALSE)
                    {
                        flint_printf("FAIL: quotient[%wd] mismatch\n", w);
                        flint_printf("iter = %wd\n", iter);
                        gr_ctx_println(cctx);
                        fflush(stdout);
                        flint_abort();
                    }
                }

                /* verify the defining identity: A == sum_w Q1[w]*B[w] + R1 */
                status = gr_mpoly_zero(recon, ctx);
                for (w = 0; w < len; w++)
                {
                    status |= gr_mpoly_mul(t, gr_mpoly_vec_entry_ptr(Q1, w, ctx),
                                               gr_mpoly_vec_entry_ptr(B, w, ctx), ctx);
                    status |= gr_mpoly_add(recon, recon, t, ctx);
                }
                status |= gr_mpoly_add(recon, recon, R1, ctx);

                if (status == GR_SUCCESS && gr_mpoly_equal(recon, A, ctx) == T_FALSE)
                {
                    flint_printf("FAIL: sum(Q_i*B_i) + R != A\n");
                    flint_printf("iter = %wd\n", iter);
                    gr_ctx_println(cctx);
                    fflush(stdout);
                    flint_abort();
                }
            }
        }

        gr_mpoly_clear(A, ctx);
        gr_mpoly_clear(R1, ctx);
        gr_mpoly_clear(R2, ctx);
        gr_mpoly_clear(recon, ctx);
        gr_mpoly_clear(t, ctx);
        gr_mpoly_vec_clear(B, ctx);
        gr_mpoly_vec_clear(Q1, ctx);
        gr_mpoly_vec_clear(Q2, ctx);
        gr_mpoly_ctx_clear(ctx);
        gr_ctx_clear(cctx);
        gr_ctx_clear(ctx1);
    }

    TEST_FUNCTION_END(state);
}
