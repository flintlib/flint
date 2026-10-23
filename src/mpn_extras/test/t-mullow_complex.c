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

/* the two's-complement residue mod B^n of a signed length result */
static void
_residue(nn_ptr r, mp_size_t n, nn_srcptr z, slong zlen)
{
    mp_size_t zn = FLINT_ABS(zlen);

    flint_mpn_zero(r, n);
    if (zn > 0)
        flint_mpn_copyi(r, z, FLINT_MIN(zn, n));
    if (zlen < 0)
        mpn_neg(r, r, n);
}

static void
_rand_part(nn_ptr x, mp_size_t * xn, mp_size_t maxn, flint_rand_t state)
{
    mp_size_t i;

    *xn = 1 + n_randint(state, maxn);

    switch (n_randint(state, 6))
    {
        case 0:
            /* zero */
            flint_mpn_zero(x, *xn);
            break;
        case 1:
            /* much shorter than its partner, with insignificant limbs */
            flint_mpn_zero(x, *xn);
            x[0] = n_randtest(state);
            break;
        case 2:
            /* all ones: extreme carries */
            for (i = 0; i < *xn; i++)
                x[i] = ~UWORD(0);
            break;
        default:
            flint_mpn_rrandom(x, state, *xn);
    }
}

TEST_FUNCTION_START(flint_mpn_mullow_complex, state)
{
    slong iter;
    slong save_sparse_cutoff = flint_mpn_mul_complex_sparse_fft_cutoff;

    for (iter = 0; iter < 3000 * flint_test_multiplier(); iter++)
    {
        mp_size_t maxn, arn, ain, brn, bin, n, zn;
        int ar_sgn, ai_sgn, br_sgn, bi_sgn, squaring, alg, ok;
        nn_ptr ar, ai, br, bi, zr, zi, rr, ri, er, ei;
        slong zr_len, zi_len;

        maxn = n_randint(state, 10) ? 40 : (n_randint(state, 4) ? 600 : 3000);
        squaring = n_randint(state, 4) == 0;

        ar = flint_malloc(4 * maxn * sizeof(ulong));
        ai = ar + maxn;
        br = ai + maxn;
        bi = br + maxn;

        _rand_part(ar, &arn, maxn, state);
        _rand_part(ai, &ain, maxn, state);
        _rand_part(br, &brn, maxn, state);
        _rand_part(bi, &bin, maxn, state);
        ar_sgn = n_randint(state, 2);
        ai_sgn = n_randint(state, 2);
        br_sgn = n_randint(state, 2);
        bi_sgn = n_randint(state, 2);

        if (squaring)
        {
            br = ar; brn = arn; br_sgn = ar_sgn;
            bi = ai; bin = ain; bi_sgn = ai_sgn;
        }

        /* windows inside the product and beyond it */
        zn = FLINT_MAX(arn, ain) + FLINT_MAX(brn, bin) + 1;
        n = 1 + n_randint(state, zn + 2);

        zr = flint_malloc(2 * zn * sizeof(ulong));
        zi = zr + zn;
        er = flint_malloc(4 * n * sizeof(ulong));
        ei = er + n;
        rr = ei + n;
        ri = rr + n;

        if (squaring)
            flint_mpn_sqr_complex(zr, &zr_len, zi, &zi_len, ar, arn, ar_sgn, ai, ain, ai_sgn);
        else
            flint_mpn_mul_complex_classical(zr, &zr_len, zi, &zi_len, ar, arn, ar_sgn,
                                  ai, ain, ai_sgn, br, brn, br_sgn, bi, bin, bi_sgn);
        _residue(er, n, zr, zr_len);
        _residue(ei, n, zi, zi_len);

        /* the zero-part transformed path from any size in the dispatcher */
        flint_mpn_mul_complex_sparse_fft_cutoff = n_randint(state, 2) ? 1 : save_sparse_cutoff;

        /* 0 dispatcher, 1 classical, 2 Karatsuba, 3 transformed */
        alg = n_randint(state, 4);
        ok = 1;

        if (squaring)
        {
            switch (alg)
            {
                case 0: flint_mpn_sqrlow_complex(rr, ri, n, ar, arn, ar_sgn, ai, ain, ai_sgn); break;
                case 1: flint_mpn_sqrlow_complex_classical(rr, ri, n, ar, arn, ar_sgn, ai, ain, ai_sgn); break;
                case 2: flint_mpn_sqrlow_complex_karatsuba(rr, ri, n, ar, arn, ar_sgn, ai, ain, ai_sgn); break;
                default: ok = flint_mpn_sqrlow_complex_fft_small(rr, ri, n, ar, arn, ar_sgn, ai, ain, ai_sgn);
            }
        }
        else
        {
            switch (alg)
            {
                case 0: flint_mpn_mullow_complex(rr, ri, n, ar, arn, ar_sgn, ai, ain, ai_sgn,
                                    br, brn, br_sgn, bi, bin, bi_sgn); break;
                case 1: flint_mpn_mullow_complex_classical(rr, ri, n, ar, arn, ar_sgn, ai, ain, ai_sgn,
                                    br, brn, br_sgn, bi, bin, bi_sgn); break;
                case 2: flint_mpn_mullow_complex_karatsuba(rr, ri, n, ar, arn, ar_sgn, ai, ain, ai_sgn,
                                    br, brn, br_sgn, bi, bin, bi_sgn); break;
                default: ok = flint_mpn_mullow_complex_fft_small(rr, ri, n, ar, arn, ar_sgn, ai, ain, ai_sgn,
                                    br, brn, br_sgn, bi, bin, bi_sgn);
            }
        }

        /* the transformed path may decline (no fft_small) */
        if (ok && (mpn_cmp(rr, er, n) != 0 || mpn_cmp(ri, ei, n) != 0))
        {
            flint_printf("FAIL: flint_mpn_%s_complex, algorithm %d\n",
                         squaring ? "sqrlow" : "mullow", alg);
            flint_printf("n = %wd, arn = %wd, ain = %wd, brn = %wd, bin = %wd, signs %d %d %d %d\n",
                         n, arn, ain, brn, bin, ar_sgn, ai_sgn, br_sgn, bi_sgn);
            flint_printf("real %s, imaginary %s\n",
                         mpn_cmp(rr, er, n) ? "wrong" : "ok",
                         mpn_cmp(ri, ei, n) ? "wrong" : "ok");
            fflush(stdout);
            flint_abort();
        }

        flint_free(ar);
        flint_free(zr);
        flint_free(er);
    }

    flint_mpn_mul_complex_sparse_fft_cutoff = save_sparse_cutoff;

    TEST_FUNCTION_END(state);
}
