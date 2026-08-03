/*
    Copyright (C) 2025 The FLINT team

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "nmod.h"
#include "nmod_vec.h"
#include "gr.h"
#include "gr_poly.h"
#include "fft_small.h"

/* Direct test of the transformed nmod polynomial ring: accumulate
   K products in transform space against an nmod_vec reference, with
   the workload forced so that small, unprofitable shapes -- in
   particular sub-block transform lengths, which reach the ring only
   through production gates otherwise -- are exercised, along with
   moduli beyond 50 bits and windowed exports (the export applies the
   deferred normalizer, so every comparison checks it). */
TEST_FUNCTION_START(gr_transformed_nmod_poly, state)
{
    for (slong iter = 0; iter < 200 * flint_test_multiplier(); iter++)
    {
        nmod_t mod;
        ulong modbits = 2 + n_randint(state, 61);
        int large = n_randint(state, 6) == 0;
        slong K = 1 + n_randint(state, 3);
        slong an = 1 + n_randint(state, large ? 600 : 90);
        slong bn = 1 + n_randint(state, large ? 600 : 90);
        slong zn = an + bn - 1;
        slong zl = n_randint(state, zn);
        slong zh = zl + 1 + n_randint(state, zn - zl);
        ulong terms;
        gr_ctx_t ctx, tctx;
        gr_transformed_poly_workload_t wl;
        nn_ptr a, b, t, ref, got;
        gr_ptr A, B, T, Z;
        slong k, glen;
        int status = GR_SUCCESS;

        nmod_init(&mod, n_randbits(state, modbits));

        terms = (ulong) K * (ulong) FLINT_MIN(an, bn) + 2;

        wl->num_inputs = 2 * K;
        wl->num_muls = K;
        wl->num_outputs = 1;
        wl->num_live = 0;
        wl->mem_limit = 0;
        wl->force = 1;

        gr_ctx_init_nmod(ctx, mod.n);

        if (gr_ctx_init_gr_poly_transformed_repr(tctx, ctx, zn,
                        (slong) terms, wl) != GR_SUCCESS)
        {
            /* implementation bounds (e.g. the prime product cannot
               cover the declared capacity) */
            gr_ctx_clear(ctx);
            continue;
        }

        a = _nmod_vec_init(K * an);
        b = _nmod_vec_init(K * bn);
        t = _nmod_vec_init(zn);
        ref = _nmod_vec_init(zn);
        got = _nmod_vec_init(zn);

        for (k = 0; k < K * an; k++) a[k] = n_randint(state, mod.n);
        for (k = 0; k < K * bn; k++) b[k] = n_randint(state, mod.n);

        /* reference: ref = sum_k a_k * b_k as plain convolutions */
        _nmod_vec_zero(ref, zn);
        for (k = 0; k < K; k++)
        {
            slong i, j;
            _nmod_vec_zero(t, zn);
            for (i = 0; i < an; i++)
                for (j = 0; j < bn; j++)
                    t[i + j] = nmod_add(t[i + j],
                        nmod_mul(a[k*an + i], b[k*bn + j], mod), mod);
            _nmod_vec_add(ref, ref, t, zn, mod);
        }

        A = gr_heap_init(tctx);
        B = gr_heap_init(tctx);
        T = gr_heap_init(tctx);
        Z = gr_heap_init(tctx);

        status |= gr_zero(Z, tctx);
        for (k = 0; k < K; k++)
        {
            status |= _gr_set_gr_poly(A, GR_ENTRY(a, k*an, sizeof(ulong)),
                                      an, ctx, tctx);
            status |= _gr_set_gr_poly(B, GR_ENTRY(b, k*bn, sizeof(ulong)),
                                      bn, ctx, tctx);
            status |= gr_mul(T, A, B, tctx);
            status |= gr_add(Z, Z, T, tctx);
        }

        if (status != GR_SUCCESS)
            TEST_FUNCTION_FAIL("ring ops declined: modbits %wu an %wd "
                "bn %wd K %wd\n", modbits, an, bn, K);

        /* full export */
        glen = 0;
        status |= _gr_get_gr_poly(got, &glen, Z, ctx, tctx);
        if (status != GR_SUCCESS || glen > zn)
            TEST_FUNCTION_FAIL("export declined: modbits %wu an %wd "
                "bn %wd K %wd\n", modbits, an, bn, K);
        for (k = glen; k < zn; k++) got[k] = 0;
        if (!_nmod_vec_equal(got, ref, zn))
            TEST_FUNCTION_FAIL("full export mismatch: modbits %wu an %wd "
                "bn %wd K %wd zn %wd\n", modbits, an, bn, K, zn);

        /* windowed export */
        _nmod_vec_zero(got, zn);
        status |= _gr_get_gr_poly_window(got, Z, zl, zh, ctx, tctx);
        if (status != GR_SUCCESS)
            TEST_FUNCTION_FAIL("window declined: modbits %wu an %wd "
                "bn %wd zl %wd zh %wd\n", modbits, an, bn, zl, zh);
        if (!_nmod_vec_equal(got, ref + zl, zh - zl))
            TEST_FUNCTION_FAIL("window mismatch: modbits %wu an %wd "
                "bn %wd K %wd zl %wd zh %wd\n",
                modbits, an, bn, K, zl, zh);

        gr_heap_clear(A, tctx);
        gr_heap_clear(B, tctx);
        gr_heap_clear(T, tctx);
        gr_heap_clear(Z, tctx);
        _nmod_vec_clear(a);
        _nmod_vec_clear(b);
        _nmod_vec_clear(t);
        _nmod_vec_clear(ref);
        _nmod_vec_clear(got);
        gr_ctx_clear(tctx);
        gr_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
