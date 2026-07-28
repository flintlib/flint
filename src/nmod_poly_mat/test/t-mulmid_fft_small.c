/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
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
        int alias, success;

        flint_set_num_threads(1 + n_randint(state, 8));

        /* the occasional large iterations cross the internal threading
           thresholds of the transform and export paths, so that both
           entry-level and intra-transform threading run */
        ar = n_randint(state, large ? 3 : 5);
        k  = n_randint(state, large ? 3 : 5);
        bc = n_randint(state, large ? 3 : 5);
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

        nmod_poly_mat_mul(D, A, B);

        if (alias)
            success = nmod_poly_mat_mulmid_fft_small(A, A, B, zl, zh);
        else
            success = nmod_poly_mat_mulmid_fft_small(C, A, B, zl, zh);

        if (!success)
        {
            /* the driver may decline when its profitability or memory
               model judges the transformed representation not worth it
               (small workloads, single-prime or tiny moduli); the
               dispatcher falls back to other algorithms then, so a
               refusal leaves nothing to verify here */
            n_refused++;
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

    TEST_FUNCTION_END(state);
}
