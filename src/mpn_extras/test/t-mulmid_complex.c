/*
    Copyright (C) 2026 Fredrik Johansson

    Developed using Claude Opus 5.5

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "mpn_extras.h"
#include "ulong_extras.h"
#if FLINT_HAVE_FFT_SMALL
#include "fft_small.h"
#endif

/* limbs [lo, hi) of the two's complement of a signed length result */
static void
_window_mid(nn_ptr r, mp_size_t lo, mp_size_t hi, nn_srcptr z, slong zlen)
{
    mp_size_t zn = FLINT_ABS(zlen);
    nn_ptr t = flint_malloc(hi * sizeof(ulong));

    flint_mpn_zero(t, hi);
    if (zn > 0)
        flint_mpn_copyi(t, z, FLINT_MIN(zn, hi));
    if (zlen < 0)
        mpn_neg(t, t, hi);
    flint_mpn_copyi(r, t + lo, hi - lo);
    flint_free(t);
}

/* is |computed - exact| (mod B^w, centered) <= (bhi B + blo)? */
static int
_within_mid(nn_srcptr c, nn_srcptr e, mp_size_t w, ulong bhi, ulong blo)
{
    nn_ptr d = flint_malloc((w + 2) * sizeof(ulong));
    mp_size_t i;
    int ok;

    mpn_sub_n(d, c, e, w);
    if (d[w - 1] >> (FLINT_BITS - 1))
        mpn_neg(d, d, w);

    ok = 1;
    for (i = 2; i < w; i++)
        ok = ok && (d[i] == 0);
    if (ok)
    {
        ulong dhi = (w >= 2) ? d[1] : 0;
        ok = (dhi < bhi) || (dhi == bhi && d[0] <= blo);
    }

    flint_free(d);
    return ok;
}

static void
_rand_part_mid(nn_ptr x, mp_size_t * xn, mp_size_t maxn, flint_rand_t state)
{
    mp_size_t i;

    *xn = 1 + n_randint(state, maxn);

    switch (n_randint(state, 6))
    {
        case 0:
            flint_mpn_zero(x, *xn);
            break;
        case 1:
            flint_mpn_zero(x, *xn);
            x[0] = n_randtest(state);
            break;
        case 2:
            for (i = 0; i < *xn; i++)
                x[i] = ~UWORD(0);
            break;
        default:
            flint_mpn_rrandom(x, state, *xn);
    }
}

TEST_FUNCTION_START(flint_mpn_mulmid_complex, state)
{
    slong iter;
    slong save_sparse_cutoff = flint_mpn_mulmid_complex_sparse_fft_cutoff;

    for (iter = 0; iter < 2000 * flint_test_multiplier(); iter++)
    {
        mp_size_t maxn, arn, ain, brn, bin, zn, lo, hi, w, m;
        int ar_sgn, ai_sgn, br_sgn, bi_sgn, squaring, alg, ok;
        nn_ptr ar, ai, br, bi, zr, zi, rr, ri, er, ei;
        slong zr_len, zi_len;
        ulong bhi, blo;

        maxn = n_randint(state, 8) ? 40 : (n_randint(state, 3) ? 800 : 4000);
        squaring = n_randint(state, 4) == 0;

        ar = flint_malloc(4 * maxn * sizeof(ulong));
        ai = ar + maxn;
        br = ai + maxn;
        bi = br + maxn;

        _rand_part_mid(ar, &arn, maxn, state);
        _rand_part_mid(ai, &ain, maxn, state);
        _rand_part_mid(br, &brn, maxn, state);
        _rand_part_mid(bi, &bin, maxn, state);
        ar_sgn = n_randint(state, 2);
        ai_sgn = n_randint(state, 2);
        br_sgn = n_randint(state, 2);
        bi_sgn = n_randint(state, 2);

        /* the Newton shape: a twice as long as b */
        if (!squaring && n_randint(state, 3) == 0)
        {
            brn = bin = FLINT_MAX(1, FLINT_MIN(arn, ain) / 2);
            flint_mpn_rrandom(br, state, brn);
            flint_mpn_rrandom(bi, state, bin);
        }

        if (squaring)
        {
            br = ar; brn = arn; br_sgn = ar_sgn;
            bi = ai; bin = ain; bi_sgn = ai_sgn;
        }

        zn = FLINT_MAX(arn, ain) + FLINT_MAX(brn, bin) + 1;

        switch (n_randint(state, 3))
        {
            case 0:
                /* the middle third, the Newton window */
                lo = FLINT_MIN(FLINT_MAX(brn, bin), zn - 1);
                hi = FLINT_MIN(FLINT_MAX(arn, ain) + 1, zn + 1);
                if (hi <= lo)
                    hi = lo + 1;
                break;
            case 1:
                lo = 0;
                hi = 1 + n_randint(state, zn + 1);
                break;
            default:
                lo = n_randint(state, zn + 1);
                hi = lo + 1 + n_randint(state, zn + 1 - lo + 1);
        }
        w = hi - lo;

        zr = flint_malloc(2 * zn * sizeof(ulong));
        zi = zr + zn;
        er = flint_malloc(4 * w * sizeof(ulong));
        ei = er + w;
        rr = ei + w;
        ri = rr + w;

        if (squaring)
            flint_mpn_sqr_complex(zr, &zr_len, zi, &zi_len, ar, arn, ar_sgn, ai, ain, ai_sgn);
        else
            flint_mpn_mul_complex_classical(zr, &zr_len, zi, &zi_len, ar, arn, ar_sgn,
                                  ai, ain, ai_sgn, br, brn, br_sgn, bi, bin, bi_sgn);
        _window_mid(er, lo, hi, zr, zr_len);
        _window_mid(ei, lo, hi, zi, zi_len);

        /* the zero-part transformed path from any size in the dispatcher */
        flint_mpn_mulmid_complex_sparse_fft_cutoff = n_randint(state, 2) ? 1 : save_sparse_cutoff;

        /* 0 dispatcher, 1 classical, 2 Karatsuba, 3 transformed */
        alg = n_randint(state, 4);
        ok = 1;

        if (squaring)
        {
            switch (alg)
            {
                case 0: flint_mpn_sqrmid_complex(rr, ri, lo, hi, ar, arn, ar_sgn, ai, ain, ai_sgn); break;
                case 1: flint_mpn_sqrmid_complex_classical(rr, ri, lo, hi, ar, arn, ar_sgn, ai, ain, ai_sgn); break;
                case 2: flint_mpn_sqrmid_complex_karatsuba(rr, ri, lo, hi, ar, arn, ar_sgn, ai, ain, ai_sgn); break;
                default: ok = flint_mpn_sqrmid_complex_fft_small(rr, ri, lo, hi, ar, arn, ar_sgn, ai, ain, ai_sgn);
            }
        }
        else
        {
            switch (alg)
            {
                case 0: flint_mpn_mulmid_complex(rr, ri, lo, hi, ar, arn, ar_sgn, ai, ain, ai_sgn,
                                    br, brn, br_sgn, bi, bin, bi_sgn); break;
                case 1: flint_mpn_mulmid_complex_classical(rr, ri, lo, hi, ar, arn, ar_sgn, ai, ain, ai_sgn,
                                    br, brn, br_sgn, bi, bin, bi_sgn); break;
                case 2: flint_mpn_mulmid_complex_karatsuba(rr, ri, lo, hi, ar, arn, ar_sgn, ai, ain, ai_sgn,
                                    br, brn, br_sgn, bi, bin, bi_sgn); break;
                default: ok = flint_mpn_mulmid_complex_fft_small(rr, ri, lo, hi, ar, arn, ar_sgn, ai, ain, ai_sgn,
                                    br, brn, br_sgn, bi, bin, bi_sgn);
            }
        }

        /* the documented error bounds: exact for lo = 0, one unit on the
           transformed path, 4 (m + 1) B units otherwise */
        m = FLINT_MIN(FLINT_MAX(arn, ain), FLINT_MAX(brn, bin));
        if (lo == 0)
            bhi = 0, blo = 0;
        else if (alg == 3)
            bhi = 0, blo = 1;
        else
            bhi = 4 * (m + 1), blo = 0;

        if (ok && !(_within_mid(rr, er, w, bhi, blo) && _within_mid(ri, ei, w, bhi, blo)))
        {
            flint_printf("FAIL: flint_mpn_%s_complex, algorithm %d\n",
                         squaring ? "sqrmid" : "mulmid", alg);
            flint_printf("lo = %wd, hi = %wd, arn = %wd, ain = %wd, brn = %wd, bin = %wd, signs %d %d %d %d\n",
                         lo, hi, arn, ain, brn, bin, ar_sgn, ai_sgn, br_sgn, bi_sgn);
            flint_printf("real %s, imaginary %s\n",
                         _within_mid(rr, er, w, bhi, blo) ? "ok" : "wrong",
                         _within_mid(ri, ei, w, bhi, blo) ? "ok" : "wrong");
            fflush(stdout);
            flint_abort();
        }

        flint_free(ar);
        flint_free(zr);
        flint_free(er);
    }

#if FLINT_HAVE_FFT_SMALL
    /* The transformed path at the edge of its cyclic window: shapes
       whose wrapped coefficients land just below the first coefficient
       reaching the window (zlen - 2^d slightly above c_lo would corrupt
       the window); random windows almost never land there. */
    for (iter = 0; iter < 10 * flint_test_multiplier(); iter++)
    {
        mpn_ctx_struct * R = get_default_mpn_ctx();
        fft_small_plan_t P;
        mp_size_t an, bn, lo, hi, w, zn;
        ulong bits, nn, alen, blen, zlen, d, wl, k, c_target;
        nn_ptr a, zr, zi, rr, ri, er, ei;
        slong zr_len, zi_len;
        int sg[4], j;

        an = 300 + n_randint(state, 3000);
        bn = an / 2 + n_randint(state, 3);

        if (!fft_small_plan_init_mpn(P, R, an, bn, 2, 1))
            continue;
        bits = P->bits;
        nn = (P->crts + P->np - 1)->coeff_len;
        fft_small_plan_clear(P);

        alen = n_cdiv(FLINT_BITS * an, bits);
        blen = n_cdiv(FLINT_BITS * bn, bits);
        zlen = alen + blen - 1;
        d = n_flog2(zlen);
        wl = n_pow2(d);
        k = 1 + n_randint(state, 40);

        if (d < LG_BLK_SZ || n_round_up(alen, BLK_SZ) > wl || zlen - wl <= k)
            continue;

        c_target = zlen - wl - k;
        lo = (mp_size_t) ((c_target * bits) / FLINT_BITS + nn);
        hi = (mp_size_t) ((wl * bits) / FLINT_BITS) - 1;
        if (hi <= lo + 1 || hi > an + bn)
            continue;
        w = hi - lo;
        zn = an + bn + 1;

        a = flint_malloc(4 * an * sizeof(ulong));
        flint_mpn_rrandom(a, state, 4 * an);
        zr = flint_malloc((2 * zn + 4 * w) * sizeof(ulong));
        zi = zr + zn;
        rr = zi + zn;
        ri = rr + w;
        er = ri + w;
        ei = er + w;

        for (j = 0; j < 4; j++)
            sg[j] = n_randint(state, 2);

        flint_mpn_mul_complex(zr, &zr_len, zi, &zi_len, a, an, sg[0], a + an, an, sg[1],
                              a + 2 * an, bn, sg[2], a + 3 * an, bn, sg[3]);
        _window_mid(er, lo, hi, zr, zr_len);
        _window_mid(ei, lo, hi, zi, zi_len);

        if (flint_mpn_mulmid_complex_fft_small(rr, ri, lo, hi, a, an, sg[0], a + an, an, sg[1],
                              a + 2 * an, bn, sg[2], a + 3 * an, bn, sg[3])
            && !(_within_mid(rr, er, w, 0, 1) && _within_mid(ri, ei, w, 0, 1)))
        {
            flint_printf("FAIL: flint_mpn_mulmid_complex_fft_small at the cyclic edge\n");
            flint_printf("an = %wd, bn = %wd, lo = %wd, hi = %wd\n", an, bn, lo, hi);
            fflush(stdout);
            flint_abort();
        }

        flint_free(a);
        flint_free(zr);
    }
#endif

    flint_mpn_mulmid_complex_sparse_fft_cutoff = save_sparse_cutoff;

    TEST_FUNCTION_END(state);
}
