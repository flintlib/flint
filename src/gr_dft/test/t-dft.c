/*
    Copyright (C) 2026 Fredrik Johansson
    Developed using Claude Fable 5

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include "test_helpers.h"
#include "thread_pool.h"
#include "thread_support.h"
#include "ulong_extras.h"
#include "gr.h"
#include "gr_vec.h"
#include "gr_dft.h"

/* which: 0 = Z/pZ with an NTT prime, 1 = CC (acb), 2 = CC (ca, exact) */

TEST_FUNCTION_START(gr_dft, state)
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
        gr_ptr x, y1, y2, w;
        ulong * perm;
        gr_dft_pre_t P, Pref;

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

        status |= gr_dft_precomp_init_root(Pref, w, n, GR_DFT_ALG_NAIVE, 0, ctx);

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
        y1 = gr_heap_init_vec(n, ctx);
        y2 = gr_heap_init_vec(n, ctx);
        perm = flint_malloc(n * sizeof(ulong));

        status |= _gr_vec_randtest(x, state, n, ctx);

        status |= gr_dft_precomp(y1, x, Pref, ctx);

        if (n_randint(state, 2))
        {
            status |= gr_dft_precomp(y2, x, P, ctx);
        }
        else
        {
            /* aliased input and output */
            status |= _gr_vec_set(y2, x, n, ctx);
            status |= gr_dft_precomp(y2, y2, P, ctx);
        }

        gr_dft_precomp_output_perm(perm, P);

        if (status == GR_SUCCESS)
        {
            for (i = 0; i < n; i++)
            {
                truth_t eq = gr_equal(GR_ENTRY(y2, i, sz),
                        GR_ENTRY(y1, perm[i], sz), ctx);

                if (eq == T_FALSE || (which != 1 && eq != T_TRUE))
                {
                    flint_printf("FAIL\n\n");
                    gr_ctx_println(ctx);
                    flint_printf("n = %wu, alg = %d, flags = %d, karatsuba = %d, i = %wu\n",
                            n, alg, flags, karatsuba, i);
                    flint_printf("x = "); _gr_vec_print(x, n, ctx); flint_printf("\n");
                    flint_printf("y1 = "); _gr_vec_print(y1, n, ctx); flint_printf("\n");
                    flint_printf("y2 = "); _gr_vec_print(y2, n, ctx); flint_printf("\n");
                    flint_abort();
                }
            }
        }

        gr_heap_clear_vec(x, n, ctx);
        gr_heap_clear_vec(y1, n, ctx);
        gr_heap_clear_vec(y2, n, ctx);
        gr_heap_clear(w, ctx);
        flint_free(perm);
        gr_dft_precomp_clear(P);
        gr_dft_precomp_clear(Pref);
        if (handles != NULL)
            flint_give_back_threads(handles, num_handles);
        flint_set_num_threads(1);
        gr_ctx_clear(ctx);
        if (have_real)
            gr_ctx_clear(rctx);
    }


    /* Deterministic algorithm sweep: every explicit algorithm at fixed
       sizes, in natural and scrambled order where supported, serially
       and with the threaded code paths forced, so that coverage of
       each algorithm does not depend on the random draws above.  Also
       exercises the plain gr_dft_precomp_init constructor, the
       one-shot gr_dft / gr_dft_inverse interface, and the plan
       printer. */
    {
        /* ring 0: acb; ring 1: acb at 2048 bits with the complex
           Karatsuba tables (which only arise above the precision
           cutoff and keep the serial twiddle layout with class
           dispatch); ring 2: nmod, length 2 (the only depth-1 case,
           covering the plain batch butterflies used when no class or
           packed tables exist) */
        static const struct { int alg; ulong n; int flags; int threads;
                int ring; }
        sweep[] = {
            { GR_DFT_ALG_CT,        16, 0,                 1, 0 },
            { GR_DFT_ALG_CT,        16, GR_DFT_SCRAMBLED,  1, 0 },
            { GR_DFT_ALG_CT,        64, 0,                 3, 0 },
            { GR_DFT_ALG_CT,        64, GR_DFT_SCRAMBLED,  3, 0 },
            { GR_DFT_ALG_CT,        16, 0,                 1, 1 },
            { GR_DFT_ALG_CT,        16, GR_DFT_SCRAMBLED,  1, 1 },
            { GR_DFT_ALG_CT,         2, 0,                 1, 2 },
            { GR_DFT_ALG_CT,         2, GR_DFT_SCRAMBLED,  1, 2 },
            { GR_DFT_ALG_CT,     32768, 0,                 1, 2 },
            { GR_DFT_ALG_BAILEY,    16, 0,                 1, 0 },
            { GR_DFT_ALG_BAILEY,    16, GR_DFT_SCRAMBLED,  1, 0 },
            { GR_DFT_ALG_BAILEY,    64, 0,                 3, 0 },
            { GR_DFT_ALG_BAILEY,    64, GR_DFT_SCRAMBLED,  3, 0 },
            { GR_DFT_ALG_SPLIT,     16, 0,                 1, 0 },
            { GR_DFT_ALG_SPLIT,     64, 0,                 4, 0 },
            { GR_DFT_ALG_SPLIT,     16, 0,                 1, 1 },
            { GR_DFT_ALG_NAIVE,      7, GR_DFT_SCRAMBLED,  1, 0 },
            { GR_DFT_ALG_MIXED,     48, 0,                 1, 0 },
            { GR_DFT_ALG_PFA,       15, 0,                 1, 0 },
            { GR_DFT_ALG_BLUESTEIN, 23, 0,                 1, 0 },
            { GR_DFT_ALG_BLUESTEIN, 23, 0,                 1, 1 },
        };
        slong si;
        gr_ctx_t actx, kctx, krctx, nctx;
        FILE * devnull = tmpfile();
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
            gr_dft_pre_t P, Pref;
            gr_ptr x, y1, y2;
            ulong * perm;
            ulong i;
            int status = GR_SUCCESS;
            /* a naive reference costs n^2 multiplications, so large
               lengths (which exercise the blocked recursion rather
               than the arithmetic) are verified by an inverse
               roundtrip instead */
            int roundtrip = (n > 4096);

            flint_set_num_threads(sweep[si].threads);

            if (sweep[si].ring == 2)
            {
                /* explicit primitive root of unity mod q */
                gr_ptr w = gr_heap_init(ctx);
                status |= gr_set_ui(w, n_primitive_root_prime(q), ctx);
                status |= gr_pow_ui(w, w, (q - 1) / n, ctx);
                if (!roundtrip)
                    status |= gr_dft_precomp_init_root(Pref, w, n,
                            GR_DFT_ALG_NAIVE, 0, ctx);
                status |= gr_dft_precomp_init_root(P, w, n,
                        sweep[si].alg, sweep[si].flags, ctx);
                gr_heap_clear(w, ctx);
            }
            else if (sweep[si].ring == 1)
            {
                status |= gr_dft_precomp_init(Pref, n, GR_DFT_ALG_NAIVE,
                        0, ctx);
                status |= gr_dft_precomp_init_karatsuba(P, n, sweep[si].alg,
                        sweep[si].flags, krctx, ctx);
            }
            else
            {
                status |= gr_dft_precomp_init(Pref, n, GR_DFT_ALG_NAIVE,
                        0, ctx);
                status |= gr_dft_precomp_init(P, n, sweep[si].alg,
                        sweep[si].flags, ctx);
            }
            if (sweep[si].threads > 1)
                gr_dft_precomp_set_serial_block(P, 1);

            if (devnull != NULL)
                gr_dft_precomp_fprint(devnull, P);

            x = gr_heap_init_vec(n, ctx);
            y1 = gr_heap_init_vec(n, ctx);
            y2 = gr_heap_init_vec(n, ctx);
            perm = flint_malloc(n * sizeof(ulong));

            status |= _gr_vec_randtest(x, state, n, ctx);
            status |= gr_dft_precomp(y2, x, P, ctx);
            if (roundtrip)
                status |= gr_dft_inverse_precomp(y1, y2, P, ctx);
            else
                status |= gr_dft_precomp(y1, x, Pref, ctx);
            gr_dft_precomp_output_perm(perm, P);

            if (status != GR_SUCCESS)
                TEST_FUNCTION_FAIL("sweep status, alg = %d, n = %wu, "
                        "flags = %d\n", sweep[si].alg, n, sweep[si].flags);

            for (i = 0; i < n; i++)
            {
                truth_t eq = roundtrip
                    ? gr_equal(GR_ENTRY(y1, i, sz), GR_ENTRY(x, i, sz), ctx)
                    : gr_equal(GR_ENTRY(y2, i, sz),
                            GR_ENTRY(y1, perm[i], sz), ctx);
                if (eq == T_FALSE || (sweep[si].ring == 2 && eq != T_TRUE))
                    TEST_FUNCTION_FAIL("sweep mismatch, alg = %d, n = %wu, "
                            "flags = %d, threads = %d, i = %wu\n",
                            sweep[si].alg, n, sweep[si].flags,
                            sweep[si].threads, i);
            }

            /* one-shot interface against the same reference */
            if (sweep[si].flags == 0 && sweep[si].threads == 1 &&
                sweep[si].ring == 0)
            {
                status = gr_dft(y2, x, n, ctx);
                for (i = 0; i < n; i++)
                    if (status != GR_SUCCESS ||
                        gr_equal(GR_ENTRY(y2, i, sz),
                            GR_ENTRY(y1, i, sz), ctx) == T_FALSE)
                        TEST_FUNCTION_FAIL("one-shot dft, n = %wu\n", n);

                status = gr_dft_inverse(y2, y1, n, ctx);
                for (i = 0; i < n; i++)
                    if (status != GR_SUCCESS ||
                        gr_equal(GR_ENTRY(y2, i, sz),
                            GR_ENTRY(x, i, sz), ctx) == T_FALSE)
                        TEST_FUNCTION_FAIL("one-shot roundtrip, n = %wu\n", n);
            }

            gr_heap_clear_vec(x, n, ctx);
            gr_heap_clear_vec(y1, n, ctx);
            gr_heap_clear_vec(y2, n, ctx);
            flint_free(perm);
            gr_dft_precomp_clear(P);
            if (!roundtrip)
                gr_dft_precomp_clear(Pref);
            flint_set_num_threads(1);
        }

        /* invalid algorithm/length combinations are rejected */
        {
            gr_dft_pre_t P;

            if (gr_dft_precomp_init(P, 3, GR_DFT_ALG_CT, 0, actx)
                    != GR_DOMAIN ||
                gr_dft_precomp_init(P, 8, GR_DFT_ALG_PFA, 0, actx)
                    != GR_DOMAIN ||
                gr_dft_precomp_init(P, 8, GR_DFT_ALG_BLUESTEIN, 0, actx)
                    != GR_DOMAIN)
                TEST_FUNCTION_FAIL("invalid plan accepted\n");
        }

        if (devnull != NULL)
            fclose(devnull);
        gr_ctx_clear(actx);
        gr_ctx_clear(kctx);
        gr_ctx_clear(krctx);
        gr_ctx_clear(nctx);
    }

    TEST_FUNCTION_END(state);
}
