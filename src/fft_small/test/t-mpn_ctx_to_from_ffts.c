/*
    Copyright (C) 2026 Fredrik Johansson

    Developed using Claude Fable 5

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <string.h>
#include <math.h>
#include "test_helpers.h"
#include "ulong_extras.h"
#include "fft_small.h"

/* Stress the templated digit-splitting input passes directly: for
   every instantiated (np, bits) pair, extract digits from random
   mpn data through the vectorized pass and compare each residue
   against an independently computed digit modulo each prime. */

static ulong
ref_digit(nn_srcptr x, slong n, slong bits, slong i, nmod_t mod)
{
    /* value of the bits-wide window at bit offset bits*i, mod p */
    slong lo = bits * i, hi = lo + bits;
    ulong r = 0, pow = 1;
    slong t;

    for (t = lo; t < hi && t < FLINT_BITS * n; t++)
    {
        if ((x[t / FLINT_BITS] >> (t % FLINT_BITS)) & 1)
            r = nmod_add(r, pow, mod);
        pow = nmod_add(pow, pow, mod);
    }
    return r;
}

TEST_FUNCTION_START(mpn_ctx_to_from_ffts, state)
{
    static const slong cfg[][2] = {
        { 3, 64 }, { 4, 64 }, { 4, 84 }, { 4, 88 }, { 4, 92 },
        { 5, 84 }, { 5, 88 }, { 5, 92 },
        { 5, 112 }, { 5, 116 }, { 5, 120 },
        { 6, 112 }, { 6, 116 }, { 6, 120 },
        { 6, 126 }, { 6, 128 }, { 6, 136 }, { 6, 140 }, { 6, 144 },
        { 7, 126 }, { 7, 128 }, { 7, 136 }, { 7, 140 }, { 7, 144 },
        { 7, 160 }, { 7, 164 }, { 7, 168 },
        { 8, 160 }, { 8, 164 }, { 8, 168 },
        { 8, 184 }, { 8, 188 }, { 8, 192 }
    };
    mpn_ctx_struct * R = get_default_mpn_ctx();
    slong iter;

    for (iter = 0; iter < 40 * flint_test_multiplier(); iter++)
    {
        slong ci = n_randint(state, sizeof(cfg) / sizeof(cfg[0]));
        slong np = cfg[ci][0], bits = cfg[ci][1];
        slong m = WORD(32) << n_randint(state, 3);
        slong n = (bits * m) / FLINT_BITS + 1, i, k;
        to_ffts_func tf = _mpn_ctx_to_ffts_func((ulong) np,
                                                (ulong) bits);
        ulong stop_easy, stop_hard, align;
        slong depth = FLINT_BIT_COUNT(m) - 1, dsz;
        nn_ptr x;
        double * buf;
        const vec4d * tp;

        if (tf == NULL)
            TEST_FUNCTION_FAIL("missing pass np = %wd, bits = %wd\n",
                               np, bits);

        for (k = 0; k < np; k++)
            sd_fft_ctx_fit_depth(R->ffts + k, depth);
        dsz = (slong) sd_fft_ctx_data_size(depth);

        x = flint_malloc(n * sizeof(ulong));
        for (i = 0; i < n; i++)
            x[i] = n_randlimb(state);
        if (n_randint(state, 4) == 0)
            x[n - 1] = 0;

        buf = flint_aligned_alloc(FLINT_FFT_SMALL_ALIGNMENT,
                  n_round_up(np * dsz * sizeof(double),
                             FLINT_FFT_SMALL_ALIGNMENT));

        align = bits % 8 == 0 ? 4 : bits % 4 == 0 ? 8 : 16;
        stop_easy = ((FLINT_BITS*(ulong)(n - 1) - 33)/bits)
                        & -(ulong) align;
        stop_easy = FLINT_MIN(stop_easy, (ulong) m & -(ulong) align);
        stop_hard = FLINT_MIN((ulong) m,
                        (FLINT_BITS*(ulong)(n - 1) + bits - 1)/bits);
        tp = R->vec_two_pow_tab[(np + 3)/4 - 1];

        tf(R->ffts, buf, dsz, x, n - 1, (ulong) m, tp,
           0, stop_easy, stop_easy, stop_hard);

        for (k = 0; k < np; k++)
        {
            nmod_t mod;
            nmod_init(&mod, R->ffts[k].p);
            for (i = 0; i < (slong) stop_hard; i++)
            {
                double v = (buf + k * dsz)[i];
                slong iv = (slong) v;
                ulong got, want;

                if (v != floor(v))
                    TEST_FUNCTION_FAIL("nonintegral residue "
                        "np = %wd bits = %wd k = %wd i = %wd\n",
                        np, bits, k, i);
                got = iv < 0 ? (ulong)(iv + (slong) mod.n)
                             : (ulong) iv;
                got = n_mod2_preinv(got, mod.n, mod.ninv);
                want = ref_digit(x, n - 1, bits, i, mod);
                if (got != want)
                    TEST_FUNCTION_FAIL("np = %wd bits = %wd k = %wd "
                        "i = %wd (m = %wd)\n", np, bits, k, i, m);
            }
        }
        flint_aligned_free(buf);
        flint_free(x);
    }

    TEST_FUNCTION_END(state);
}
