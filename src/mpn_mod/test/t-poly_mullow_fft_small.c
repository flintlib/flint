/*
    Copyright (C) 2024 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "mpn_mod.h"
#include "gr_poly.h"

TEST_FUNCTION_START(mpn_mod_poly_mullow_fft_small, state)
{
    gr_ctx_t ctx;
    slong iter;

    for (iter = 0; iter < 100 * flint_test_multiplier(); iter++)
    {
        gr_ctx_init_mpn_mod_randtest(ctx, state);
        flint_set_num_threads(1 + n_randint(state, 4));
        _gr_poly_test_mullow((gr_method_poly_binary_trunc_op) _mpn_mod_poly_mullow_fft_small, NULL, state, 10, 50, ctx);
        flint_set_num_threads(1);
        gr_ctx_clear(ctx);
    }

    for (iter = 0; iter < 100 * flint_test_multiplier(); iter++)
    {
        gr_ctx_init_mpn_mod_randtest(ctx, state);
        flint_set_num_threads(1 + n_randint(state, 4));
        _gr_poly_test_mulmid((gr_method_poly_binary_trunc2_op) _mpn_mod_poly_mulmid_fft_small, NULL, state, 10, 50, ctx);
        flint_set_num_threads(1);
        gr_ctx_clear(ctx);
    }

    for (iter = 0; iter < 10 * flint_test_multiplier(); iter++)
    {
        gr_ctx_init_mpn_mod_randtest(ctx, state);
        flint_set_num_threads(1 + n_randint(state, 4));
        _gr_poly_test_mullow((gr_method_poly_binary_trunc_op) _mpn_mod_poly_mullow_fft_small,
                            (gr_method_poly_binary_trunc_op) _mpn_mod_poly_mullow_KS, state, 10, 5000, ctx);
        flint_set_num_threads(1);
        gr_ctx_clear(ctx);
    }

    for (iter = 0; iter < 10 * flint_test_multiplier(); iter++)
    {
        gr_ctx_init_mpn_mod_randtest(ctx, state);
        flint_set_num_threads(1 + n_randint(state, 4));
        _gr_poly_test_mulmid((gr_method_poly_binary_trunc2_op) _mpn_mod_poly_mulmid_fft_small,
                            (gr_method_poly_binary_trunc2_op) _mpn_mod_poly_mulmid_KS, state, 10, 2000, ctx);
        flint_set_num_threads(1);
        gr_ctx_clear(ctx);
    }

    #if FLINT_HAVE_FFT_SMALL
    /* ring compliance for the transformed mpn_mod polynomial contexts:
       forced initialization so small unprofitable sizes -- where bugs
       are easiest to catch -- are testable */
    {
        slong it, n_engaged = 0;
        for (it = 0; it < 4; it++)
        {
            gr_ctx_t base, tctx;
            fmpz_t m;
            gr_transformed_poly_workload_t wl = {{ 0, 0, 0, 0, 0, 1 }};
            slong nb, Nc;

            /* the first attempt is deterministic and within the crt
               capacity on any fft_small platform, so the engagement
               assertion below cannot fail on unlucky draws: moduli
               whose coefficient products exceed the prime set refuse
               even when forced (a genuine implementation bound); the
               last attempt deliberately reaches beyond it to keep the
               refusal path exercised, without asserting refusal, since
               a larger prime table may admit it */
            if (it == 0)
            { nb = 2 * FLINT_BITS - 3; Nc = 6; }
            else if (it == 3)
            { nb = 5 * FLINT_BITS; Nc = 4; }
            else
            { nb = FLINT_BITS * (2 + (slong) n_randint(state, 2));
              Nc = 2 + (slong) n_randint(state, 24); }

            fmpz_init(m);
            fmpz_randprime(m, state, nb, 0);
            if (gr_ctx_init_mpn_mod(base, m) == GR_SUCCESS)
            {
                if (_gr_mpn_mod_ctx_init_transformed_poly_repr(tctx, base,
                        Nc, Nc, wl) == GR_SUCCESS)
                {
                    gr_test_ring(tctx, 15, 0);
                    gr_ctx_clear(tctx);
                    n_engaged++;
                }
                gr_ctx_clear(base);
            }
            fmpz_clear(m);
        }
        if (n_engaged == 0)
            TEST_FUNCTION_FAIL("forced tpoly context never engaged\n");
    }

    /* the unforced constructor exercises the profitability model in
       both directions: a heavy declared workload accepts, the default
       declines at unprofitable sizes */
    {
        gr_ctx_t base, tctx;
        fmpz_t m;
        gr_transformed_poly_workload_t heavy =
            {{ 64, 100000, 64, 0, 0, 0 }};
        fmpz_init(m);
        fmpz_randprime(m, state, 3 * FLINT_BITS, 0);
        if (gr_ctx_init_mpn_mod(base, m) == GR_SUCCESS)
        {
            if (_gr_mpn_mod_ctx_init_transformed_poly_repr(tctx, base,
                    64, 64, heavy) == GR_SUCCESS)
                gr_ctx_clear(tctx);
            if (_gr_mpn_mod_ctx_init_transformed_poly_repr(tctx, base,
                    4, 2, NULL) == GR_SUCCESS)
                gr_ctx_clear(tctx);
            gr_ctx_clear(base);
        }
        fmpz_clear(m);
    }

    /* directed exercise: negated accumulation and the destructive
       conversion out, against a gr_poly reference over the base */
    {
        slong it;
        for (it = 0; it < 6; it++)
        {
            gr_ctx_t base, tctx;
            fmpz_t m;
            gr_transformed_poly_workload_t wl = {{ 0, 0, 0, 0, 0, 1 }};
            slong N = 4 + (slong) n_randint(state, 20);
            slong half = FLINT_MAX(1, N / 2);
            fmpz_init(m);
            fmpz_randprime(m, state,
                    FLINT_BITS * (2 + (slong) n_randint(state, 3)), 0);
            if (gr_ctx_init_mpn_mod(base, m) != GR_SUCCESS)
            { fmpz_clear(m); continue; }
            if (_gr_mpn_mod_ctx_init_transformed_poly_repr(tctx, base, N,
                    N, wl) != GR_SUCCESS)
            { gr_ctx_clear(base); fmpz_clear(m); continue; }
            {
                gr_poly_t a, b, c, d, ref, t2, got;
                gr_ptr xa, xb, xc, xd, acc;
                slong len = 0;
                gr_poly_init(a, base); gr_poly_init(b, base);
                gr_poly_init(c, base); gr_poly_init(d, base);
                gr_poly_init(ref, base); gr_poly_init(t2, base);
                gr_poly_init(got, base);
                GR_MUST_SUCCEED(gr_poly_randtest(a, state, half, base));
                GR_MUST_SUCCEED(gr_poly_randtest(b, state, half, base));
                GR_MUST_SUCCEED(gr_poly_randtest(c, state, half, base));
                GR_MUST_SUCCEED(gr_poly_randtest(d, state, half, base));
                GR_MUST_SUCCEED(gr_poly_set_coeff_si(a, 0, 1, base));
                GR_MUST_SUCCEED(gr_poly_set_coeff_si(b, 0, 1, base));
                GR_MUST_SUCCEED(gr_poly_set_coeff_si(c, 0, 1, base));
                GR_MUST_SUCCEED(gr_poly_set_coeff_si(d, 0, 1, base));
                xa = gr_heap_init(tctx); xb = gr_heap_init(tctx);
                xc = gr_heap_init(tctx); xd = gr_heap_init(tctx);
                acc = gr_heap_init(tctx);
                GR_MUST_SUCCEED(_gr_set_gr_poly(xa, a->coeffs, a->length,
                        base, tctx));
                GR_MUST_SUCCEED(_gr_set_gr_poly(xb, b->coeffs, b->length,
                        base, tctx));
                GR_MUST_SUCCEED(_gr_set_gr_poly(xc, c->coeffs, c->length,
                        base, tctx));
                GR_MUST_SUCCEED(_gr_set_gr_poly(xd, d->coeffs, d->length,
                        base, tctx));
                GR_MUST_SUCCEED(gr_mul(acc, xa, xb, tctx));
                GR_MUST_SUCCEED(gr_submul(acc, xc, xd, tctx));
                /* aliased accumulation must route through a temporary
                   internally; UNABLE is a legitimate answer at tight
                   provisioning, a crash or GR_DOMAIN is not */
                {
                    int st1 = gr_addmul(acc, acc, xb, tctx);
                    (void) st1;
                }
                gr_heap_clear(acc, tctx);
                acc = gr_heap_init(tctx);
                GR_MUST_SUCCEED(gr_mul(acc, xa, xb, tctx));
                GR_MUST_SUCCEED(gr_submul(acc, xc, xd, tctx));
                GR_MUST_SUCCEED(gr_poly_mul(ref, a, b, base));
                GR_MUST_SUCCEED(gr_poly_mul(t2, c, d, base));
                GR_MUST_SUCCEED(gr_poly_sub(ref, ref, t2, base));
                gr_poly_fit_length(got, N, base);
                GR_MUST_SUCCEED(_gr_get_gr_poly_destructive(got->coeffs,
                        &len, acc, base, tctx));
                _gr_poly_set_length(got, len, base);
                _gr_poly_normalise(got, base);
                if (gr_poly_equal(got, ref, base) != T_TRUE)
                    TEST_FUNCTION_FAIL("mpn_mod negs destructive get: "
                            "N=%wd\n", N);
                gr_heap_clear(xa, tctx); gr_heap_clear(xb, tctx);
                gr_heap_clear(xc, tctx); gr_heap_clear(xd, tctx);
                gr_heap_clear(acc, tctx);
                gr_poly_clear(a, base); gr_poly_clear(b, base);
                gr_poly_clear(c, base); gr_poly_clear(d, base);
                gr_poly_clear(ref, base); gr_poly_clear(t2, base);
                gr_poly_clear(got, base);
            }
            gr_ctx_clear(tctx);
            gr_ctx_clear(base);
            fmpz_clear(m);
        }
    }
#endif

    TEST_FUNCTION_END(state);
}
