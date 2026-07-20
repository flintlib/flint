/*
    Copyright (C) 2026 Fredrik Johansson
    Developed using Claude Fable 5

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "thread_pool.h"
#include "thread_support.h"
#include "ulong_extras.h"
#include "gr.h"
#include "gr_vec.h"
#include "gr_dft.h"

TEST_FUNCTION_START(gr_dft_inverse, state)
{
    slong iter;

    for (iter = 0; iter < 200 * flint_test_multiplier(); iter++)
    {
        gr_ctx_t ctx, rctx;
        int which = n_randint(state, 3);
        int have_real = 0, karatsuba = 0;
        int status = GR_SUCCESS;
        int alg, flags;
        int use_threads, attach_threads;
        thread_pool_handle * handles = NULL;
        slong num_handles = 0;
        ulong n, i, q = 0;
        slong sz;
        gr_ptr x, y, z, w;
        gr_dft_pre_t P;

        if (which == 0)
        {
            /* arbitrary n; pick a prime q = 1 mod n so that Z/qZ
               contains n-th roots of unity */
            ulong k;
            n = 1 + n_randint(state, 300);
            for (k = 1; ; k++)
            {
                q = k * n + 1;
                if (q >= 3 && n_is_prime(q))
                    break;
            }
            gr_ctx_init_nmod(ctx, q);
        }
        else if (which == 1)
        {
            slong prec = 64 + n_randint(state, 128);
            gr_ctx_init_complex_acb(ctx, prec);
            gr_ctx_init_real_arb(rctx, prec);
            have_real = 1;

            switch (n_randint(state, 4))
            {
                case 0:
                    /* prime lengths exercising Bluestein */
                    {
                        static const ulong bl_primes[] =
                            { 23, 29, 31, 37, 41, 101, 127 };
                        n = bl_primes[n_randint(state, 7)];
                    }
                    break;
                case 1:
                    n = UWORD(1) << n_randint(state, 8);
                    break;
                default:
                    n = 1 + n_randint(state, 200);
                    break;
            }
        }
        else
        {
            gr_ctx_init_complex_ca(ctx);
            n = 1 + n_randint(state, 6);
        }

        sz = ctx->sizeof_elem;

        /* pick an algorithm valid for this n and ring; Bluestein is
           restricted to rings supporting gr_dft_default_root for the
           power-of-two convolution length */
        {
            int algs[8];
            int nalgs = 0, pow2 = ((n & (n - 1)) == 0), ncop;
            ulong m, pp;

            for (ncop = 0, m = n, pp = 2; m > 1; pp++)
            {
                if (m % pp == 0)
                {
                    ncop++;
                    while (m % pp == 0)
                        m /= pp;
                }
            }

            algs[nalgs++] = GR_DFT_ALG_AUTO;
            algs[nalgs++] = GR_DFT_ALG_NAIVE;
            if (pow2)
            {
                algs[nalgs++] = GR_DFT_ALG_CT;
                algs[nalgs++] = GR_DFT_ALG_BAILEY;
                algs[nalgs++] = GR_DFT_ALG_SPLIT;
            }
            if (n >= 2)
                algs[nalgs++] = GR_DFT_ALG_MIXED;
            if (ncop >= 2)
                algs[nalgs++] = GR_DFT_ALG_PFA;
            if (n >= 3 && n % 2 == 1 && which == 1)
                algs[nalgs++] = GR_DFT_ALG_BLUESTEIN;

            alg = algs[n_randint(state, nalgs)];
        }
        flags = n_randint(state, 2) ? GR_DFT_SCRAMBLED : 0;
        karatsuba = have_real && n_randint(state, 2);

        /* thread settings are applied for all rings; the library only
           actually uses threads when gr_ctx_is_threadsafe certifies the
           ring (so e.g. ca plans run serially regardless) */
        use_threads = n_randint(state, 2);
        attach_threads = use_threads && n_randint(state, 2);
        flint_set_num_threads(use_threads ? 2 + n_randint(state, 2) : 1);

        w = gr_heap_init(ctx);

        if (which == 0)
        {
            /* w = g^((q-1)/n) mod q for a primitive root g mod q */
            status |= gr_set_ui(w, n_primitive_root_prime(q), ctx);
            status |= gr_pow_ui(w, w, (q - 1) / n, ctx);
        }
        else
        {
            status |= gr_dft_default_root(w, n, ctx);
        }

        if (karatsuba)
            status |= gr_dft_precomp_init_karatsuba(P, n, alg, flags, rctx, ctx);
        else
            status |= gr_dft_precomp_init_root(P, w, n, alg, flags, ctx);

        if (status == GR_SUCCESS && use_threads)
        {
            /* force threaded code paths even at small n */
            gr_dft_precomp_set_serial_block(P, 1 + (slong) n_randint(state, 16));

            if (attach_threads)
            {
                num_handles = flint_request_threads(&handles, 3);
                gr_dft_precomp_set_threads(P, handles, num_handles);
            }
        }

        x = gr_heap_init_vec(n, ctx);
        y = gr_heap_init_vec(n, ctx);
        z = gr_heap_init_vec(n, ctx);

        status |= _gr_vec_randtest(x, state, n, ctx);

        status |= gr_dft_precomp(y, x, P, ctx);

        if (n_randint(state, 2))
        {
            status |= gr_dft_inverse_precomp(z, y, P, ctx);
        }
        else
        {
            status |= _gr_vec_set(z, y, n, ctx);
            status |= gr_dft_inverse_precomp(z, z, P, ctx);
        }

        if (status == GR_SUCCESS)
        {
            for (i = 0; i < n; i++)
            {
                truth_t eq = gr_equal(GR_ENTRY(z, i, sz),
                        GR_ENTRY(x, i, sz), ctx);

                if (eq == T_FALSE || (which != 1 && eq != T_TRUE))
                {
                    flint_printf("FAIL\n\n");
                    gr_ctx_println(ctx);
                    flint_printf("n = %wu, alg = %d, flags = %d, karatsuba = %d, i = %wu\n",
                            n, alg, flags, karatsuba, i);
                    flint_printf("x = "); _gr_vec_print(x, n, ctx); flint_printf("\n");
                    flint_printf("y = "); _gr_vec_print(y, n, ctx); flint_printf("\n");
                    flint_printf("z = "); _gr_vec_print(z, n, ctx); flint_printf("\n");
                    flint_abort();
                }
            }
        }

        gr_heap_clear_vec(x, n, ctx);
        gr_heap_clear_vec(y, n, ctx);
        gr_heap_clear_vec(z, n, ctx);
        gr_heap_clear(w, ctx);
        gr_dft_precomp_clear(P);
        if (handles != NULL)
            flint_give_back_threads(handles, num_handles);
        flint_set_num_threads(1);
        gr_ctx_clear(ctx);
        if (have_real)
            gr_ctx_clear(rctx);
    }


    /* Deterministic inverse sweep over the explicit power-of-two
       algorithms (including Bailey, which the random draws above can
       miss), serially and with the threaded paths forced: forward
       then inverse must reproduce the input. */
    {
        /* ring 0: acb; ring 1: acb at 2048 bits with the complex
           Karatsuba tables (serial layout, DIT butterflies through the
           class dispatch); ring 2: nmod (depth-1 and blocked-recursion
           inverse paths) */
        static const struct { int alg; ulong n; int threads; int ring; }
        sweep[] = {
            { GR_DFT_ALG_CT,        16, 1, 0 },
            { GR_DFT_ALG_CT,        64, 3, 0 },
            { GR_DFT_ALG_CT,        16, 1, 1 },
            { GR_DFT_ALG_CT,         2, 1, 2 },
            { GR_DFT_ALG_CT,     32768, 1, 2 },
            { GR_DFT_ALG_BAILEY,    16, 1, 0 },
            { GR_DFT_ALG_BAILEY,    64, 3, 0 },
            { GR_DFT_ALG_SPLIT,     16, 1, 0 },
            { GR_DFT_ALG_SPLIT,     64, 4, 0 },
            { GR_DFT_ALG_SPLIT,     16, 1, 1 },
        };
        slong si;
        gr_ctx_t actx, kctx, krctx, nctx;
        ulong q = 998244353;

        gr_ctx_init_complex_acb(actx, 100);
        gr_ctx_init_complex_acb(kctx, 2048);
        gr_ctx_init_real_arb(krctx, 2048);
        gr_ctx_init_nmod(nctx, q);

        for (si = 0; si < (slong) (sizeof(sweep) / sizeof(sweep[0])); si++)
        {
            ulong n = sweep[si].n;
            gr_ctx_struct * ctx = (sweep[si].ring == 2) ? nctx :
                    (sweep[si].ring == 1) ? kctx : actx;
            slong sz = ctx->sizeof_elem;
            gr_dft_pre_t P;
            gr_ptr x, y, z;
            ulong i;
            int status = GR_SUCCESS;

            flint_set_num_threads(sweep[si].threads);

            if (sweep[si].ring == 2)
            {
                gr_ptr w = gr_heap_init(ctx);
                status |= gr_set_ui(w, n_primitive_root_prime(q), ctx);
                status |= gr_pow_ui(w, w, (q - 1) / n, ctx);
                status |= gr_dft_precomp_init_root(P, w, n,
                        sweep[si].alg, 0, ctx);
                gr_heap_clear(w, ctx);
            }
            else if (sweep[si].ring == 1)
            {
                status |= gr_dft_precomp_init_karatsuba(P, n, sweep[si].alg,
                        0, krctx, ctx);
            }
            else
            {
                status |= gr_dft_precomp_init(P, n, sweep[si].alg, 0, ctx);
            }
            if (sweep[si].threads > 1)
                gr_dft_precomp_set_serial_block(P, 1);

            x = gr_heap_init_vec(n, ctx);
            y = gr_heap_init_vec(n, ctx);
            z = gr_heap_init_vec(n, ctx);

            status |= _gr_vec_randtest(x, state, n, ctx);
            status |= gr_dft_precomp(y, x, P, ctx);
            status |= gr_dft_inverse_precomp(z, y, P, ctx);

            if (status != GR_SUCCESS)
                TEST_FUNCTION_FAIL("inverse sweep status, alg = %d, "
                        "n = %wu\n", sweep[si].alg, n);

            for (i = 0; i < n; i++)
            {
                truth_t eq = gr_equal(GR_ENTRY(z, i, sz),
                        GR_ENTRY(x, i, sz), ctx);
                if (eq == T_FALSE || (sweep[si].ring == 2 && eq != T_TRUE))
                    TEST_FUNCTION_FAIL("inverse sweep roundtrip, alg = %d, "
                            "n = %wu, threads = %d, i = %wu\n",
                            sweep[si].alg, n, sweep[si].threads, i);
            }

            gr_heap_clear_vec(x, n, ctx);
            gr_heap_clear_vec(y, n, ctx);
            gr_heap_clear_vec(z, n, ctx);
            gr_dft_precomp_clear(P);
            flint_set_num_threads(1);
        }

        gr_ctx_clear(actx);
        gr_ctx_clear(kctx);
        gr_ctx_clear(krctx);
        gr_ctx_clear(nctx);
    }

    TEST_FUNCTION_END(state);
}
