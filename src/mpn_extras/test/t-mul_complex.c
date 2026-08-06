/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz.h"
#include "mpn_extras.h"

static void
_set_signed(fmpz_t f, nn_srcptr a, slong n, int s)
{
    /* a returned length may legitimately be zero: ar = ai with equal
       signs makes the real part of a square vanish */
    if (n <= 0)
    {
        fmpz_zero(f);
        return;
    }

    fmpz_set_ui_array(f, a, n);
    if (s)
        fmpz_neg(f, f);
}

/* Error bound of the high variants, measured against the exact value:
   one flint_mpn_mulhigh_n read as n limbs errs by (-1 - eps, +eps) ulp
   of its lowest limb, eps = (n + 4)/2^64 -- the floor drop of its
   returned limb is one-sided in [0, 1) and the documented n + 2 ulp of
   that limb contributes eps either way. Composing: two products per
   output classically gives |error| < 2 + 2 eps; the Karatsuba mulhigh
   combination c3 - (c1 + c2) gives < 2 + 3 eps; the transformed path
   errs within (-1.5, +0.5). Since the outputs are integers they differ
   from the floor of the exact value by at most 3 in the lowest returned
   limb, which is what is asserted here (observed worst case: 2).
   Returns 0 on success, 1 if the magnitude is off by more than tol,
   2 if the sign is wrong. */
static int
_check_high(const fmpz_t ref, nn_srcptr z, slong n, int sgn, slong nn,
            slong tol)
{
    fmpz_t want, got, diff;
    int res = 0;

    fmpz_init(want);
    fmpz_init(got);
    fmpz_init(diff);

    fmpz_abs(want, ref);
    fmpz_fdiv_q_2exp(want, want, FLINT_BITS * nn);
    fmpz_set_ui_array(got, z, n);
    fmpz_sub(diff, got, want);
    fmpz_abs(diff, diff);

    if (fmpz_cmp_ui(diff, tol) > 0)
        res = 1;
    else if (!fmpz_is_zero(want) && !fmpz_is_zero(got) &&
                sgn != (fmpz_sgn(ref) < 0))
        res = 2;

    fmpz_clear(want);
    fmpz_clear(got);
    fmpz_clear(diff);

    return res;
}

/* Complex multiplication of limb arrays with separate sign bits, on both
   the classical and the transformed (fft_small) paths: sizes straddle
   the internal cutoff, and zero components exercise the canonical zero. */
TEST_FUNCTION_START(flint_mpn_mul_complex, state)
{
    /* some iterations force the transformed paths regardless of size,
       covering the square pointwise ops and the signed exports */
    slong save_mul_cutoff = flint_mpn_mul_complex_fft_cutoff;
    slong save_sqr_cutoff = flint_mpn_sqr_complex_fft_cutoff;

    for (slong iter = 0; iter < 40 * flint_test_multiplier(); iter++)
    {
        if (iter % 3 == 1)
        {
            flint_mpn_mul_complex_fft_cutoff = 1;
            flint_mpn_sqr_complex_fft_cutoff = 1;
        }
        else
        {
            flint_mpn_mul_complex_fft_cutoff = save_mul_cutoff;
            flint_mpn_sqr_complex_fft_cutoff = save_sqr_cutoff;
        }

        /* the square path: under the forced cutoff this reaches the
           transformed square ops and their signed exports; the
           Karatsuba square is the reference */
        if (iter % 3 == 1)
        {
            slong sn = 1 + n_randint(state, 300);
            nn_ptr sa = flint_malloc(sn * sizeof(ulong));
            nn_ptr sb = flint_malloc(sn * sizeof(ulong));
            nn_ptr z1 = flint_malloc((2 * sn + 2) * sizeof(ulong));
            nn_ptr z2 = flint_malloc((2 * sn + 2) * sizeof(ulong));
            nn_ptr z3 = flint_malloc((2 * sn + 2) * sizeof(ulong));
            nn_ptr z4 = flint_malloc((2 * sn + 2) * sizeof(ulong));
            slong l1, l2, l3, l4, q;
            int s1 = (int) n_randint(state, 2), s2 = (int) n_randint(state, 2);
            for (q = 0; q < sn; q++)
            { sa[q] = n_randtest(state); sb[q] = n_randtest(state); }
            if (flint_mpn_sqr_complex_fft_small(z1, &l1, z2, &l2,
                    sa, sn, s1, sb, sn, s2))
            {
                flint_mpn_sqr_complex_karatsuba(z3, &l3, z4, &l4,
                        sa, sn, s1, sb, sn, s2);
                if (l1 != l3 || l2 != l4
                    || (l1 != 0 && mpn_cmp(z1, z3, FLINT_ABS(l1)) != 0)
                    || (l2 != 0 && mpn_cmp(z2, z4, FLINT_ABS(l2)) != 0))
                    TEST_FUNCTION_FAIL("transformed square vs karatsuba: "
                            "sn=%wd s=%d%d\n", sn, s1, s2);
            }
            flint_free(sa); flint_free(sb);
            flint_free(z1); flint_free(z2); flint_free(z3); flint_free(z4);
        }


        /* the transformed path is off by default, so drive it explicitly
           on some iterations; sizes straddle the Karatsuba threshold */
        int fftpath = (iter % 4 == 0);
        slong n = fftpath ? 1024 + n_randint(state, 300)
                          : 1 + n_randint(state, 200);
        slong arn, ain, brn, bin, wfull, wsqr;
        int cancel_a, cancel_b;

        flint_mpn_mul_complex_fft_cutoff = fftpath ? 1 : WORD_MAX;
        flint_mpn_sqr_complex_fft_cutoff = fftpath ? 1 : WORD_MAX;

        /* the transformed conversions thread from a size threshold */
        flint_set_num_threads(1 + n_randint(state, 4));

        /* independent part lengths; every third iteration is
           deliberately lopsided so the schoolbook selection is covered
           alongside the Karatsuba and transformed ones */
        arn = n; ain = n; brn = n; bin = n;
        if (iter % 3 == 0)
        {
            arn = 1 + n_randint(state, n);
            ain = 1 + n_randint(state, n);
            brn = 1 + n_randint(state, n);
            bin = 1 + n_randint(state, n);
        }

        /* sometimes make the Karatsuba auxiliary sums cancel: with
           ai = ar and opposite signs, ar + ai collapses, which shrinks
           the intermediate lengths -- where a width taken from the wrong
           buffer shows up. Decided here, before the output widths are
           computed from the lengths. */
        cancel_a = (n_randint(state, 5) == 0);
        cancel_b = (n_randint(state, 5) == 0);
        if (cancel_a)
            arn = ain = n;
        if (cancel_b)
            brn = bin = n;
        nn_ptr ar, ai, br, bi, zr, zi, hr, hi;
        fmpz_t fa, fb, fc, fd, refr, refi, got;
        int sar, sai, sbr, sbi, shr, shi;
        slong lzr, lzi;
        slong j;

        ar = flint_malloc(n * sizeof(ulong));
        ai = flint_malloc(n * sizeof(ulong));
        br = flint_malloc(n * sizeof(ulong));
        bi = flint_malloc(n * sizeof(ulong));
        wfull = FLINT_MAX(arn, ain) + FLINT_MAX(brn, bin) + 1;
        wsqr = 2 * FLINT_MAX(arn, ain) + 1;
        zr = flint_malloc(FLINT_MAX(wfull, wsqr) * sizeof(ulong));
        zi = flint_malloc(FLINT_MAX(wfull, wsqr) * sizeof(ulong));
        hr = flint_malloc((n + 1) * sizeof(ulong));
        hi = flint_malloc((n + 1) * sizeof(ulong));

        for (j = 0; j < n; j++)
        {
            ar[j] = n_randlimb(state);
            ai[j] = n_randlimb(state);
            br[j] = n_randlimb(state);
            bi[j] = n_randlimb(state);
        }

        if (cancel_a)
            flint_mpn_copyi(ai, ar, n);
        if (cancel_b)
            flint_mpn_copyi(bi, br, n);

        /* zero components must not disturb the signs */
        if (n_randint(state, 6) == 0)
            flint_mpn_zero(ai, n);
        if (n_randint(state, 6) == 0)
            flint_mpn_zero(br, n);

        sar = n_randint(state, 2);
        sai = n_randint(state, 2);
        sbr = n_randint(state, 2);
        sbi = n_randint(state, 2);

        fmpz_init(fa);
        fmpz_init(fb);
        fmpz_init(fc);
        fmpz_init(fd);
        fmpz_init(refr);
        fmpz_init(refi);
        fmpz_init(got);

        _set_signed(fa, ar, arn, sar);
        _set_signed(fb, ai, ain, sai);
        _set_signed(fc, br, brn, sbr);
        _set_signed(fd, bi, bin, sbi);

        /* (ar + i ai)(br + i bi) */
        fmpz_mul(refr, fa, fc);
        fmpz_submul(refr, fb, fd);
        fmpz_mul(refi, fa, fd);
        fmpz_addmul(refi, fb, fc);

        /* exercise the three public methods directly as well as the
           dispatcher, whatever the shape */
        if (iter % 4 == 1)
            flint_mpn_mul_complex_classical(zr, &lzr, zi, &lzi, ar, arn, sar,
                    ai, ain, sai, br, brn, sbr, bi, bin, sbi);
        else if (iter % 4 == 2)
            flint_mpn_mul_complex_karatsuba(zr, &lzr, zi, &lzi, ar, arn, sar,
                    ai, ain, sai, br, brn, sbr, bi, bin, sbi);
        else if (iter % 4 == 3 &&
                 flint_mpn_mul_complex_fft_small(zr, &lzr, zi, &lzi, ar, arn,
                     sar, ai, ain, sai, br, brn, sbr, bi, bin, sbi))
            ;
        else
            flint_mpn_mul_complex(zr, &lzr, zi, &lzi, ar, arn, sar, ai, ain,
                                  sai, br, brn, sbr, bi, bin, sbi);

        if (FLINT_ABS(lzr) > wfull || FLINT_ABS(lzi) > wfull)
            TEST_FUNCTION_FAIL("mul: returned length exceeds the documented "
                               "capacity (%wd, %wd against %wd)\n",
                               lzr, lzi, wfull);
        if ((lzr != 0 && zr[FLINT_ABS(lzr) - 1] == 0) ||
            (lzi != 0 && zi[FLINT_ABS(lzi) - 1] == 0))
            TEST_FUNCTION_FAIL("mul: returned length not normalized\n");

        _set_signed(got, zr, FLINT_ABS(lzr), lzr < 0);
        if (!fmpz_equal(got, refr))
            TEST_FUNCTION_FAIL("mul: real part wrong (arn = %wd, ain = %wd, "
                               "brn = %wd, bin = %wd)\n", arn, ain, brn, bin);
        _set_signed(got, zi, FLINT_ABS(lzi), lzi < 0);
        if (!fmpz_equal(got, refi))
            TEST_FUNCTION_FAIL("mul: imaginary part wrong (arn = %wd, "
                               "ain = %wd, brn = %wd, bin = %wd)\n",
                               arn, ain, brn, bin);

        /* the high variants take one length for all four parts */
        _set_signed(fa, ar, n, sar);
        _set_signed(fb, ai, n, sai);
        _set_signed(fc, br, n, sbr);
        _set_signed(fd, bi, n, sbi);
        fmpz_mul(refr, fa, fc);
        fmpz_submul(refr, fb, fd);
        fmpz_mul(refi, fa, fd);
        fmpz_addmul(refi, fb, fc);

        flint_mpn_mulhigh_n_complex(hr, &shr, hi, &shi, ar, sar, ai, sai,
                                    br, sbr, bi, sbi, n);
        if (_check_high(refr, hr, n + 1, shr, n, 3))
            TEST_FUNCTION_FAIL("mulhigh: real part wrong (n = %wd)\n", n);
        if (_check_high(refi, hi, n + 1, shi, n, 3))
            TEST_FUNCTION_FAIL("mulhigh: imaginary part wrong (n = %wd)\n", n);

        /* (ar + i ai)^2, per-part lengths */
        _set_signed(fa, ar, arn, sar);
        _set_signed(fb, ai, ain, sai);
        fmpz_mul(refr, fa, fa);
        fmpz_submul(refr, fb, fb);
        fmpz_mul(refi, fa, fb);
        fmpz_mul_2exp(refi, refi, 1);

        if (iter % 4 == 1)
            flint_mpn_sqr_complex_classical(zr, &lzr, zi, &lzi,
                    ar, arn, sar, ai, ain, sai);
        else if (iter % 4 == 2)
            flint_mpn_sqr_complex_karatsuba(zr, &lzr, zi, &lzi,
                    ar, arn, sar, ai, ain, sai);
        else if (iter % 4 == 3 &&
                 flint_mpn_sqr_complex_fft_small(zr, &lzr, zi, &lzi,
                     ar, arn, sar, ai, ain, sai))
            ;
        else
            flint_mpn_sqr_complex(zr, &lzr, zi, &lzi, ar, arn, sar,
                                  ai, ain, sai);

        _set_signed(got, zr, FLINT_ABS(lzr), lzr < 0);
        if (!fmpz_equal(got, refr))
            TEST_FUNCTION_FAIL("sqr: real part wrong (arn = %wd, ain = %wd)\n",
                               arn, ain);
        _set_signed(got, zi, FLINT_ABS(lzi), lzi < 0);
        if (!fmpz_equal(got, refi))
            TEST_FUNCTION_FAIL("sqr: imaginary part wrong (arn = %wd, "
                               "ain = %wd)\n", arn, ain);

        /* the high square again takes one length */
        _set_signed(fa, ar, n, sar);
        _set_signed(fb, ai, n, sai);
        fmpz_mul(refr, fa, fa);
        fmpz_submul(refr, fb, fb);
        fmpz_mul(refi, fa, fb);
        fmpz_mul_2exp(refi, refi, 1);

        flint_mpn_sqrhigh_n_complex(hr, &shr, hi, &shi, ar, sar, ai, sai, n);
        if (_check_high(refr, hr, n + 1, shr, n, 3))
            TEST_FUNCTION_FAIL("sqrhigh: real part wrong (n = %wd)\n", n);
        if (_check_high(refi, hi, n + 1, shi, n, 3))
            TEST_FUNCTION_FAIL("sqrhigh: imaginary part wrong (n = %wd)\n", n);

        flint_free(ar);
        flint_free(ai);
        flint_free(br);
        flint_free(bi);
        flint_free(zr);
        flint_free(zi);
        flint_free(hr);
        flint_free(hi);
        fmpz_clear(fa);
        fmpz_clear(fb);
        fmpz_clear(fc);
        fmpz_clear(fd);
        fmpz_clear(refr);
        fmpz_clear(refi);
        fmpz_clear(got);
    }

    flint_mpn_mul_complex_fft_cutoff = 1024;
    flint_mpn_sqr_complex_fft_cutoff = 4096;

    flint_mpn_mul_complex_fft_cutoff = save_mul_cutoff;
    flint_mpn_sqr_complex_fft_cutoff = save_sqr_cutoff;
    TEST_FUNCTION_END(state);
}
