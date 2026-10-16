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
#include "mpn_extras.h"
#include "thread_support.h"
#include "fft.h"
#include "fft_small.h"

/* fft_convolution_sd_fft against classic fft_convolution on
   identical inputs, across transform depths, thread counts (to
   exercise the threaded pointwise helpers and both sides of the
   fused-versus-flat dispatch), truncation, and squaring. */

TEST_FUNCTION_START(fft_convolution_sd_fft, state)
{
#if FLINT_HAVE_FFT_SMALL
    slong depth, ti, iter;
    static const int threads[] = { 1, 2, 4, 10 };
    slong save_threads = flint_get_num_threads();

    for (depth = 3; depth <= 8; depth++)
    for (ti = 0; ti < 4; ti++)
    for (iter = 0; iter < 2; iter++)
    {
        mp_size_t n = (UWORD(1) << depth);
        mp_size_t n1 = (UWORD(1) << (depth/2));
        slong b, mm, N = 0;
        mp_size_t limbs, size, trunc, i, len1, len2;
        mp_limb_t * ptr;
        mp_limb_t ** ii, ** jj, ** ii2, ** jj2;
        mp_limb_t ** t1, ** t2, ** s1, ** tt;
        sd_fft_mpn_mulmod_2expp1_ctx_t C;
        slong nt;
        int squaring;

        /* an engine-compatible ring that the outer transform
           accepts: N a multiple of n with an even w */
        for (b = 64; N == 0 && b <= 128; b += 64)
            for (mm = 16; mm <= 512; mm *= 2)
                if ((mm * b) % n == 0 && ((mm * b) / n) % 2 == 0
                    && mm * b >= 24 * n)
                {
                    N = mm * b;
                    break;
                }
        if (N == 0)
            continue;

        flint_set_num_threads(threads[ti]);
        nt = flint_get_num_threads();

        limbs = N / FLINT_BITS;
        size = limbs + 1;
        trunc = 2*n + n_randint(state, 2*n) + 1;
        trunc = 2*n1*((trunc + 2*n1 - 1)/(2*n1));
        squaring = (n_randint(state, 4) == 0);
        len1 = 1 + n_randint(state, trunc);
        len2 = squaring ? len1 : 1 + n_randint(state, trunc);

        sd_fft_mpn_mulmod_2expp1_ctx_init(C, get_default_mpn_ctx(),
                                          N, mm);

        ii = flint_malloc((4*(n + n*size) + 5*size*nt)*sizeof(mp_limb_t)
                          + (4*n + 5*nt)*sizeof(mp_limb_t *));
        for (i = 0, ptr = (mp_limb_t *)(ii + 4*n + 5*nt); i < 4*n;
             i++, ptr += size)
        {
            ii[i] = ptr;
            if (i < len1)
                random_fermat(ii[i], state, limbs);
            else
                flint_mpn_zero(ii[i], size);
            mpn_normmod_2expp1(ii[i], limbs);
        }
        t1 = (mp_limb_t **)(ii + 4*n);
        t2 = t1 + nt;
        s1 = t2 + nt;
        tt = s1 + nt;
        for (i = 0; i < nt; i++)
        {
            t1[i] = ptr; ptr += size;
            t2[i] = ptr; ptr += size;
            s1[i] = ptr; ptr += size;
            tt[i] = ptr; ptr += 2*size;
        }

        if (squaring)
            jj = ii;
        else
        {
            jj = flint_malloc(4*(n + n*size)*sizeof(mp_limb_t)
                              + 4*n*sizeof(mp_limb_t *));
            for (i = 0, ptr = (mp_limb_t *)(jj + 4*n); i < 4*n;
                 i++, ptr += size)
            {
                jj[i] = ptr;
                if (i < len2)
                    random_fermat(jj[i], state, limbs);
                else
                    flint_mpn_zero(jj[i], size);
                mpn_normmod_2expp1(jj[i], limbs);
            }
        }

        /* classic copies */
        ii2 = flint_malloc(4*(n + n*size)*sizeof(mp_limb_t)
                           + 4*n*sizeof(mp_limb_t *));
        for (i = 0, ptr = (mp_limb_t *)(ii2 + 4*n); i < 4*n;
             i++, ptr += size)
        {
            ii2[i] = ptr;
            flint_mpn_copyi(ii2[i], ii[i], size);
        }
        if (squaring)
            jj2 = ii2;
        else
        {
            jj2 = flint_malloc(4*(n + n*size)*sizeof(mp_limb_t)
                               + 4*n*sizeof(mp_limb_t *));
            for (i = 0, ptr = (mp_limb_t *)(jj2 + 4*n); i < 4*n;
                 i++, ptr += size)
            {
                jj2[i] = ptr;
                flint_mpn_copyi(jj2[i], jj[i], size);
            }
        }

        if (n_randint(state, 4) == 0)   /* NULL context = classic */
            fft_convolution_sd_fft(NULL, ii, jj, depth, limbs, trunc,
                                   t1, t2, s1, tt);
        else
            fft_convolution_sd_fft(C, ii, jj, depth, limbs, trunc,
                                   t1, t2, s1, tt);
        fft_convolution(ii2, jj2, depth, limbs, trunc,
                        t1, t2, s1, tt);

        for (i = 0; i < trunc; i++)
        {
            mpn_normmod_2expp1(ii[i], limbs);
            mpn_normmod_2expp1(ii2[i], limbs);
            if (mpn_cmp(ii[i], ii2[i], size) != 0)
                TEST_FUNCTION_FAIL("depth = %wd, N = %wd, m = %wd, "
                    "threads = %wd, trunc = %wd, i = %wd, sq = %d\n",
                    depth, N, mm, nt, (slong) trunc, (slong) i,
                    squaring);
        }

        sd_fft_mpn_mulmod_2expp1_ctx_clear(C);
        flint_free(ii);
        if (!squaring)
            flint_free(jj);
        flint_free(ii2);
        if (!squaring)
            flint_free(jj2);
    }

    flint_set_num_threads(save_threads);

    /* chooser edges: a huge tight size walks the lattice through
       divisibility/parity skips and capacity-infeasible digit
       sizes (np < 0) before settling; the result must be a valid
       ring at least as large as the request */
    {
        slong mo, r;
        r = fft_adjust_limbs_sd_fft(1000000, 16, &mo);
        if (r < 1000000 || (mo != 0 && (64*r) % mo != 0))
            TEST_FUNCTION_FAIL("chooser edge: r = %wd, m = %wd\n",
                               r, mo);
        r = fft_adjust_limbs_sd_fft(300, 1024, &mo);  /* divisibility
                                skips at a large outer length */
        if (r < 300)
            TEST_FUNCTION_FAIL("chooser edge 2: r = %wd\n", r);
    }
#endif

    TEST_FUNCTION_END(state);
}
