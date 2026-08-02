/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fft_small.h"
#include "gr.h"
#include "gr_poly.h"
#include "nmod_poly.h"
#include "nmod_poly_mat.h"

TEST_FUNCTION_START(nmod_poly_mat_mulmid_fft_small, state)

    slong n_engaged = 0, n_refused = 0;
{
#if FLINT_HAVE_FFT_SMALL
    for (slong iter = 0; iter < 100 * flint_test_multiplier(); iter++)
    {
        nmod_poly_mat_t A, B, C, D;
        nmod_poly_t t;
        ulong n;
        slong ar, k, bc, i, j;
        ulong maxlen, zn_max;
        slong zl, zh;
        int large = n_randint(state, 10) == 0;
        int square, alias, success;
        ulong save_limit = flint_fft_small_max_transformed_ring_size;

        flint_set_num_threads(1 + n_randint(state, 8));

        /* the occasional large iterations cross the internal threading
           thresholds of the transform and export paths, so that both
           entry-level and intra-transform threading run */
        ar = n_randint(state, large ? 3 : 5);
        k  = n_randint(state, large ? 3 : 5);
        bc = n_randint(state, large ? 3 : 5);
        square = (iter % 5 == 0);
        if (square)
        { k = ar; bc = ar; }
        /* a starved storage budget forces the blocked wrapper */
        if (iter % 4 == 3)
            flint_fft_small_max_transformed_ring_size = UWORD(1) << 18;
        maxlen = large ? 1500 + n_randint(state, 3000)
                       : 1 + n_randint(state, 300);

        /* random moduli, sometimes an fft prime for the matching-prime
           plan, sometimes an fft-friendly prime for the direct branch */
        switch (n_randint(state, 4))
        {
            case 0:  n = UWORD(0x0003f00000000001); break;
            case 1:  n = UWORD(1) + (UWORD(3) << 29); break;
            default: n = n_randtest_not_zero(state);
                     n += (n == 1);
                     break;
        }

        nmod_poly_mat_init(A, ar, k, n);
        nmod_poly_mat_init(B, k, bc, n);
        nmod_poly_mat_init(C, ar, bc, n);
        nmod_poly_mat_init(D, ar, bc, n);
        nmod_poly_init(t, n);

        nmod_poly_mat_randtest(A, state, maxlen);
        nmod_poly_mat_randtest(B, state, maxlen);

        /* window, occasionally reaching past the product length */
        zn_max = 2*maxlen + 2;
        zl = (slong) n_randint(state, zn_max);
        zh = zl + 1 + (slong) n_randint(state, zn_max - (ulong) zl);

        /* aliasing C with A requires matching dimensions, so k == bc;
           the reference must be computed before the operand is clobbered */
        alias = (k == bc) && n_randint(state, 8) == 0;

        nmod_poly_mat_mul(D, A, square ? A : B);

        if (alias)
            success = nmod_poly_mat_mulmid_fft_small(A, A,
                    square ? A : B, zl, zh);
        else
            success = nmod_poly_mat_mulmid_fft_small(C, A,
                    square ? A : B, zl, zh);

        if (!success)
        {
            /* the driver may decline when its profitability or memory
               model judges the transformed representation not worth it
               (small workloads, single-prime or tiny moduli); the
               dispatcher falls back to other algorithms then, so a
               refusal leaves nothing to verify here */
            n_refused++;
            flint_fft_small_max_transformed_ring_size = save_limit;
            nmod_poly_mat_clear(A); nmod_poly_mat_clear(B);
            nmod_poly_mat_clear(C); nmod_poly_mat_clear(D);
            nmod_poly_clear(t);
            continue;
        }
        n_engaged++;

        for (i = 0; i < ar; i++)
        for (j = 0; j < bc; j++)
        {
            nmod_poly_struct* got = alias ? nmod_poly_mat_entry(A, i, j)
                                          : nmod_poly_mat_entry(C, i, j);

            nmod_poly_shift_right(t, nmod_poly_mat_entry(D, i, j), zl);
            nmod_poly_truncate(t, zh - zl);

            if (!nmod_poly_equal(t, got))
            {
                TEST_FUNCTION_FAIL(
                        "entry (%wd, %wd) of %wd x %wd x %wd\n"
                        "modulus = %wu, maxlen = %wu\n"
                        "zl = %wd, zh = %wd, alias = %d\n"
                        "got %{nmod_poly}\nexpected %{nmod_poly}\n",
                        i, j, ar, k, bc, n, maxlen, zl, zh, alias,
                        got, t);
            }
        }

        flint_fft_small_max_transformed_ring_size = save_limit;
        nmod_poly_mat_clear(A);
        nmod_poly_mat_clear(B);
        nmod_poly_mat_clear(C);
        nmod_poly_mat_clear(D);
        nmod_poly_clear(t);
    }
    /* the mixed sizes and moduli must exercise the engaged path; without
       fft_small the driver never runs and there is nothing to require */
    if (n_engaged == 0)
        TEST_FUNCTION_FAIL("driver never engaged (%wd refusals)\n", n_refused);
#endif

    #if FLINT_HAVE_FFT_SMALL
    /* ring compliance for the transformed nmod polynomial contexts:
       operations beyond the provisioning answer GR_UNABLE, which the
       framework accepts */
    {
        slong it, n_engaged = 0;
        for (it = 0; it < 4; it++)
        {
            gr_ctx_t base, tctx;
            /* the first modulus and size are deterministic so the
               engagement assertion cannot fail on unlucky draws */
            ulong n = (it % 2) ? UWORD(0x0003f00000000001)
                               : n_randtest_prime(state, 0);
            /* forced: only implementation bounds decline, so small
               unprofitable sizes -- where bugs are easiest to catch --
               are testable */
            gr_transformed_poly_workload_t wl = {{ 0, 0, 0, 0, 0, 1 }};
            gr_ctx_init_nmod(base, n);
            /* the terms bound must cover the convolution accumulation
               of a full product (min operand length per coefficient) */
            slong Nc = (it == 0) ? 8
                                 : 2 + (slong) n_randint(state, 40);
            if (_gr_nmod_ctx_init_transformed_poly_repr(tctx, base,
                    Nc, Nc, wl) == GR_SUCCESS)
            {
                gr_test_ring(tctx, 20, 0);
                gr_ctx_clear(tctx);
                n_engaged++;
            }
            gr_ctx_clear(base);
        }
        if (n_engaged == 0)
            TEST_FUNCTION_FAIL("forced tpoly context never engaged\n");
    }

    /* directed exercise of the conversion-out family: negated
       accumulations (the biased path), the destructive get, and the
       windowed destructive get, against nmod_poly reference */
    {
        slong it;
        for (it = 0; it < 8; it++)
        {
            gr_ctx_t base, tctx;
            gr_transformed_poly_workload_t wl = {{ 0, 0, 0, 0, 0, 1 }};
            ulong n = n_randtest_prime(state, 0);
            slong N = 2 + n_randint(state, 50);
            /* products must fit the element bound N */
            slong half = FLINT_MAX(1, N / 2);
            slong la = 1 + n_randint(state, half), lb = 1 + n_randint(state, half);
            slong zl, zh;

            gr_ctx_init_nmod(base, n);
            if (_gr_nmod_ctx_init_transformed_poly_repr(tctx, base, N,
                    N, wl) != GR_SUCCESS)
            { gr_ctx_clear(base); continue; }
            {
                nmod_poly_t a, b, c, d, ref, got;
                gr_ptr xa, xb, xc, xd, acc;
                slong len = 0;
                nmod_poly_init(a, n); nmod_poly_init(b, n);
                nmod_poly_init(c, n); nmod_poly_init(d, n);
                nmod_poly_init(ref, n); nmod_poly_init(got, n);
                nmod_poly_randtest(a, state, la);
                nmod_poly_randtest(b, state, lb);
                nmod_poly_randtest(c, state, 1 + n_randint(state, half));
                nmod_poly_randtest(d, state, 1 + n_randint(state, half));
                /* nonzero operands: the conversions accept length 0,
                   but the reference bookkeeping below is simpler
                   without it */
                nmod_poly_set_coeff_ui(a, 0, 1);
                nmod_poly_set_coeff_ui(b, 0, 1);
                nmod_poly_set_coeff_ui(c, 0, 1);
                nmod_poly_set_coeff_ui(d, 0, 1);
                xa = gr_heap_init(tctx); xb = gr_heap_init(tctx);
                xc = gr_heap_init(tctx); xd = gr_heap_init(tctx);
                acc = gr_heap_init(tctx);
                GR_MUST_SUCCEED(_gr_set_gr_poly(xa, a->coeffs, a->length, base, tctx));
                GR_MUST_SUCCEED(_gr_set_gr_poly(xb, b->coeffs, b->length, base, tctx));
                GR_MUST_SUCCEED(_gr_set_gr_poly(xc, c->coeffs, c->length, base, tctx));
                GR_MUST_SUCCEED(_gr_set_gr_poly(xd, d->coeffs, d->length, base, tctx));
                /* acc = a b - c d: submul drives the negated path */
                GR_MUST_SUCCEED(gr_mul(acc, xa, xb, tctx));
                GR_MUST_SUCCEED(gr_submul(acc, xc, xd, tctx));
                nmod_poly_mul(ref, a, b);
                nmod_poly_mul(got, c, d);
                nmod_poly_sub(ref, ref, got);
                /* windowed destructive get of [zl, zh) */
                zh = N;
                zl = n_randint(state, zh);
                nmod_poly_fit_length(got, zh - zl);
                GR_MUST_SUCCEED(_gr_nmod_tpoly_get_gr_poly_window_destructive(
                        got->coeffs, acc, zl, zh, tctx));
                _nmod_poly_set_length(got, zh - zl);
                _nmod_poly_normalise(got);
                nmod_poly_shift_right(ref, ref, zl);
                nmod_poly_truncate(ref, zh - zl);
                if (!nmod_poly_equal(got, ref))
                    TEST_FUNCTION_FAIL("windowed negs get: N=%wd n=%wu "
                            "zl=%wd zh=%wd\n", N, n, zl, zh);
                /* aliased accumulation, which the domain bookkeeping
                   must route through a temporary: acc += acc * xb and
                   acc -= xa * acc */
                GR_MUST_SUCCEED(gr_mul(acc, xa, xb, tctx));
                {
                    int st1 = gr_addmul(acc, acc, xb, tctx);
                    int st2 = gr_submul(acc, xa, acc, tctx);
                    if (st1 == GR_DOMAIN || st2 == GR_DOMAIN)
                        TEST_FUNCTION_FAIL("aliased addmul/submul "
                                "returned GR_DOMAIN\n");
                }
                gr_heap_clear(acc, tctx);
                acc = gr_heap_init(tctx);

                /* full destructive get on a fresh negated accumulation;
                   the windowed destructive get consumed acc, whose
                   contract allows only clearing afterwards */
                gr_heap_clear(acc, tctx);
                acc = gr_heap_init(tctx);
                GR_MUST_SUCCEED(gr_mul(acc, xa, xb, tctx));
                GR_MUST_SUCCEED(gr_submul(acc, xc, xd, tctx));
                nmod_poly_mul(ref, a, b);
                nmod_poly_mul(got, c, d);
                nmod_poly_sub(ref, ref, got);
                nmod_poly_fit_length(got, N);
                GR_MUST_SUCCEED(_gr_get_gr_poly_destructive(got->coeffs, &len, acc, base, tctx));
                _nmod_poly_set_length(got, len);
                _nmod_poly_normalise(got);
                if (!nmod_poly_equal(got, ref))
                    TEST_FUNCTION_FAIL("destructive negs get: N=%wd n=%wu\n",
                            N, n);
                gr_heap_clear(xa, tctx); gr_heap_clear(xb, tctx);
                gr_heap_clear(xc, tctx); gr_heap_clear(xd, tctx);
                gr_heap_clear(acc, tctx);
                nmod_poly_clear(a); nmod_poly_clear(b);
                nmod_poly_clear(c); nmod_poly_clear(d);
                nmod_poly_clear(ref); nmod_poly_clear(got);
            }
            gr_ctx_clear(tctx);
            gr_ctx_clear(base);
        }
    }
#endif

    TEST_FUNCTION_END(state);
}
