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
#include "test_helpers.h"
#include "ulong_extras.h"
#include "fmpz.h"
#include "mpn_extras.h"
#include "fft_small.h"

/* reference: z = x*y mod 2^N+1 by a plain product and two folds,
   v = f0 - f1 + f2 with f0, f1 the low/next n limbs and f2 the top
   two limbs of the (2n+2)-limb full product */
static void
negmul_ref(nn_ptr z, nn_srcptr x, nn_srcptr y, slong n)
{
    nn_ptr full = flint_malloc((2 * n + 2) * sizeof(ulong));
    ulong cy, br;

    flint_mpn_mul(full, x, n + 1, y, n + 1);
    cy = mpn_add(z, full, n, full + 2 * n, 2);
    br = mpn_sub_n(z, z, full + n, n);
    z[n] = 0;
    if (cy > br)
    {
        if (mpn_sub_1(z, z, n, 1))
        {
            memset(z, 0, n * sizeof(ulong));
            z[n] = 1;
        }
    }
    else if (br > cy)
    {
        if (mpn_add_1(z, z, n, 1))
        {
            memset(z, 0, n * sizeof(ulong));
            z[n] = 1;
        }
    }
    flint_free(full);
}

TEST_FUNCTION_START(sd_fft_mpn_mulmod_2expp1, state)
{
    /* every instantiated wide digit size, plus small-b digit
       shapes; np is chosen internally by the context */
    static const slong cfg[] = {
        64, 84, 88, 92, 112, 116, 120, 126, 128, 136, 140, 144,
        160, 164, 168, 184, 188, 192, 34, 40, 46, 50
    };
    slong iter;

    for (iter = 0; iter < 300 * flint_test_multiplier(); iter++)
    {
        slong ci = n_randint(state, sizeof(cfg) / sizeof(cfg[0]));
        slong b = cfg[ci];
        slong m = WORD(16) << n_randint(state, 4);
        slong N = m * b, n = N / FLINT_BITS;
        sd_fft_mpn_mulmod_2expp1_ctx_t C;
        nn_ptr x, y, z, r;
        slong i;

        if (N % FLINT_BITS != 0)
            continue;

        sd_fft_mpn_mulmod_2expp1_ctx_init(C, get_default_mpn_ctx(), N,
                (b > 50 && n_randint(state, 4) == 0) ? 0 : m);
        if (C->m != m)   /* default splitting may differ; adopt it */
        {
            sd_fft_mpn_mulmod_2expp1_ctx_clear(C);
            sd_fft_mpn_mulmod_2expp1_ctx_init(C, get_default_mpn_ctx(),
                                              N, m);
        }
        sd_fft_mpn_mulmod_2expp1_scratch_t S;
        sd_fft_mpn_mulmod_2expp1_scratch_init(S, C);

        /* the context must pick the least workable prime count */
        if (C->np > 3)
        {
            fmpz_t P;
            slong k;
            fmpz_init_set_ui(P, UWORD(1));
            for (k = 0; k < C->np - 1; k++)
                fmpz_mul_ui(P, P, get_default_mpn_ctx()->ffts[k].p);
            if ((slong) fmpz_bits(P) - 1 >= 2 * b
                    + FLINT_BIT_COUNT(m) - 1 + 1
                && (b <= 50
                    || _mpn_ctx_to_ffts_func((ulong)(C->np - 1),
                                             (ulong) b) != NULL))
                TEST_FUNCTION_FAIL("np = %wd not minimal for b = %wd, "
                                   "m = %wd\n", C->np, b, m);
            fmpz_clear(P);
        }

        x = flint_malloc((n + 1) * sizeof(ulong));
        y = flint_malloc((n + 1) * sizeof(ulong));
        z = flint_malloc((n + 1) * sizeof(ulong));
        r = flint_malloc((n + 1) * sizeof(ulong));

        for (i = 0; i < n; i++)
        {
            x[i] = n_randlimb(state);
            y[i] = n_randlimb(state);
        }
        x[n] = y[n] = 0;

        switch (n_randint(state, 9))
        {
            case 8:   /* unreduced: top limb set, low limbs kept */
                x[n] = 1;
                break;
            case 0:   /* x = 2^N */
                memset(x, 0, n * sizeof(ulong));
                x[n] = 1;
                break;
            case 1:   /* y = 2^N */
                memset(y, 0, n * sizeof(ulong));
                y[n] = 1;
                break;
            case 2:   /* squaring */
                memcpy(y, x, (n + 1) * sizeof(ulong));
                break;
            case 3:   /* sparse operand */
                memset(x, 0, n * sizeof(ulong));
                x[n_randint(state, n)] = n_randlimb(state);
                break;
            default:
                break;
        }

        negmul_ref(r, x, y, n);
        sd_fft_mpn_mulmod_2expp1(C, z, x, y, S);

        if (z[n] != r[n] || mpn_cmp(z, r, n) != 0)
            TEST_FUNCTION_FAIL("b = %wd, m = %wd, N = %wd, np = %wd\n",
                               b, m, N, C->np);

        /* in-place: z aliases x, as the convolution pointwise does */
        memcpy(z, x, (n + 1) * sizeof(ulong));
        sd_fft_mpn_mulmod_2expp1(C, z, z, y, S);
        if (z[n] != r[n] || mpn_cmp(z, r, n) != 0)
            TEST_FUNCTION_FAIL("aliased: b = %wd, m = %wd\n",
                               b, m);

        flint_free(x);
        flint_free(y);
        flint_free(z);
        flint_free(r);
        sd_fft_mpn_mulmod_2expp1_scratch_clear(S);
        sd_fft_mpn_mulmod_2expp1_ctx_clear(C);
    }

    /* second-family shapes: depths where the family prime count
       runs out of capacity and np is bumped (e.g. b = 92 at large
       depth needs five primes) */
    {
        static const slong cfg2[][2] = {
            { 92, 4096 }, { 128, 8192 }, { 168, 2048 }
        };
        slong ci;
        for (ci = 0; ci < 3; ci++)
        {
            slong b = cfg2[ci][0], m = cfg2[ci][1];
            slong N = m * b, n = N / FLINT_BITS, i;
            sd_fft_mpn_mulmod_2expp1_ctx_t C;
            nn_ptr x, y, z, r;

            sd_fft_mpn_mulmod_2expp1_ctx_init(C, get_default_mpn_ctx(), N,
                (b > 50 && n_randint(state, 4) == 0) ? 0 : m);
        if (C->m != m)   /* default splitting may differ; adopt it */
        {
            sd_fft_mpn_mulmod_2expp1_ctx_clear(C);
            sd_fft_mpn_mulmod_2expp1_ctx_init(C, get_default_mpn_ctx(),
                                              N, m);
        }
        sd_fft_mpn_mulmod_2expp1_scratch_t S;
        sd_fft_mpn_mulmod_2expp1_scratch_init(S, C);
            x = flint_malloc((n + 1) * sizeof(ulong));
            y = flint_malloc((n + 1) * sizeof(ulong));
            z = flint_malloc((n + 1) * sizeof(ulong));
            r = flint_malloc((n + 1) * sizeof(ulong));
            for (i = 0; i < n; i++)
            {
                x[i] = n_randlimb(state);
                y[i] = n_randlimb(state);
            }
            x[n] = y[n] = 0;
            negmul_ref(r, x, y, n);
            sd_fft_mpn_mulmod_2expp1(C, z, x, y, S);
            if (z[n] != r[n] || mpn_cmp(z, r, n) != 0)
                TEST_FUNCTION_FAIL("deep: b = %wd, m = %wd, np = %wd\n",
                                   b, m, C->np);
            flint_free(x);
            flint_free(y);
            flint_free(z);
            flint_free(r);
            sd_fft_mpn_mulmod_2expp1_scratch_clear(S);
            sd_fft_mpn_mulmod_2expp1_ctx_clear(C);
        }
    }

    TEST_FUNCTION_END(state);
}
