/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "mpn_extras.h"
#include "arb.h"
#include "fixed.h"

/* Adversarial stress inputs for the bitwise reductions: for exp,
   sparse sums of the tabulated logarithms L_i = log(1 + 2^-i)
   (perturbed by a few ulp), so that the greedy subtraction runs into
   exact matches, single-limb model ties, near-zero residuals and --
   when an index is duplicated -- the extra correction step at i = r;
   for log, sparse products of the factors 1 + 2^-i built by exact
   shift-and-adds, so that the reduction must rediscover the factor
   set, driving the deficit through zero and near-tie states. */

static double
ulp_error_vs_arb(nn_srcptr res, slong rn, const arb_t exact, slong n)
{
    arb_t va, delta;
    fmpz_t t;
    mag_t mb;
    double u;

    arb_init(va);
    arb_init(delta);
    fmpz_init(t);
    mag_init(mb);

    fmpz_set_ui_array(t, res, rn);
    arb_set_fmpz(va, t);
    arb_mul_2exp_si(va, va, -FLINT_BITS * n);
    arb_sub(delta, va, exact, ARF_PREC_EXACT);
    arb_get_mag(mb, delta);
    mag_mul_2exp_si(mb, mb, FLINT_BITS * n);
    u = mag_get_d(mb);

    arb_clear(va);
    arb_clear(delta);
    fmpz_clear(t);
    mag_clear(mb);

    return u;
}

TEST_FUNCTION_START(fixed_bitwise_rs_stress, state)
{
    slong iter;

    for (iter = 0; iter < 200 * flint_test_multiplier(); iter++)
    {
        slong nmin = (FLINT_BITS == 64) ? 1 : 2;
        slong n = nmin + n_randint(state, (iter % 5 == 0) ? 40 : 10);
        slong imax = FLINT_MIN(FLINT_BITS * n - 16, 448);
        int r = (iter % 4 == 0) ? 0
                : 32 + (int) n_randint(state, FLINT_MAX(1, imax - 32));
        int reff = (r == 0) ? 32 : r;
        slong j, k, nsum;
        ulong x[52], res[53], cy;
        arb_t xa, exact;
        fmpz_t f;
        slong prec = FLINT_BITS * n + 128;
        double u;

        /* make sure the table covers the indices we want to use */
        _fixed_exp_logs_ensure(n, FLINT_MIN((slong) reff + 2, imax));

        /* ---- exp: x = sum of tabulated L_i, sparse and possibly
           with repeats (every fifth iteration duplicates the index
           at the reduction limit, forcing the extra step) ---- */
        flint_mpn_zero(x, n);
        nsum = 1 + n_randint(state, 10);
        for (j = 0; j < nsum; j++)
        {
            slong i = 1 + n_randint(state, FLINT_MIN((slong) reff + 2,
                imax));
            nn_srcptr Li = _fixed_exp_logs_entry(i, n);

            cy = mpn_add_n(x, x, Li, n);
            if (cy)
                mpn_sub_n(x, x, Li, n);    /* would exceed 1: revert */
        }
        if (iter % 5 == 0 && reff <= imax)
        {
            /* duplicate L_r (twice on top of whatever is there) */
            nn_srcptr Lr = _fixed_exp_logs_entry(
                FLINT_MIN((slong) reff, _fixed_exp_logs_max_index()), n);

            for (k = 0; k < 2; k++)
            {
                cy = mpn_add_n(x, x, Lr, n);
                if (cy)
                    mpn_sub_n(x, x, Lr, n);
            }
        }
        /* perturb by a few ulp in either direction */
        k = (slong) n_randint(state, 7) - 3;
        if (k > 0)
        {
            cy = mpn_add_1(x, x, n, (ulong) k);
            if (cy)
                mpn_sub_1(x, x, n, (ulong) k);
        }
        else if (k < 0)
        {
            cy = mpn_sub_1(x, x, n, (ulong) (-k));
            if (cy)
                mpn_add_1(x, x, n, (ulong) (-k));
        }

        fixed_exp_bitwise_rs(res, x, n, r);

        arb_init(xa);
        arb_init(exact);
        fmpz_init(f);
        fmpz_set_ui_array(f, x, n);
        arb_set_fmpz(xa, f);
        arb_mul_2exp_si(xa, xa, -FLINT_BITS * n);
        arb_exp(exact, xa, prec);
        u = ulp_error_vs_arb(res, n + 1, exact, n);
        if (u > (double) FIXED_EXP_BITWISE_RS_MAX_ERR(n, r))
            TEST_FUNCTION_FAIL("exp: n = %wd, r = %d, ulp = %f\n",
                n, r, u);

        /* ---- log: X = product of factors 1 + 2^-i staying below 2,
           built by exact shift-and-adds (every fifth iteration
           attempts a duplicated factor at the reduction limit) ---- */
        {
            ulong P[53], sh[53];
            slong wn = n + 1, tries;
            int rl = (iter % 4 == 0) ? 0
                     : 32 + (int) n_randint(state, FLINT_MAX(1, imax - 32));

            flint_mpn_zero(P, n);
            P[n] = 1;

            tries = 2 + (slong) n_randint(state, 12);
            for (j = 0; j < tries; j++)
            {
                slong i;

                if (iter % 5 == 0 && j < 2)
                    i = FLINT_MAX(1, FLINT_MIN((rl == 0) ? 32 : rl,
                        imax));
                else
                    i = 1 + n_randint(state, imax);

                if (i >= FLINT_BITS * wn)
                    continue;
                mpn_rshift(sh, P, wn, 1);     /* dummy init */
                {
                    slong q = i / FLINT_BITS, b = i % FLINT_BITS;

                    flint_mpn_zero(sh, wn);
                    if (b)
                        mpn_rshift(sh, P + q, wn - q, (int) b);
                    else
                        flint_mpn_copyi(sh, P + q, wn - q);
                }
                cy = mpn_add(P, P, wn, sh, n);   /* low n limbs of v */
                if (cy || P[n] > 1)
                {
                    mpn_sub(P, P, wn, sh, n);    /* crossed 2: revert */
                }
            }

            /* perturb the fraction by a few ulp */
            k = (slong) n_randint(state, 7) - 3;
            if (k > 0)
            {
                cy = mpn_add_1(P, P, n, (ulong) k);
                if (cy)
                    mpn_sub_1(P, P, n, (ulong) k);
            }
            else if (k < 0)
            {
                cy = mpn_sub_1(P, P, n, (ulong) (-k));
                if (cy)
                    mpn_add_1(P, P, n, (ulong) (-k));
            }

            fixed_log1p_bitwise_rs(res, P, n, rl);

            fmpz_set_ui_array(f, P, n);
            arb_set_fmpz(xa, f);
            arb_mul_2exp_si(xa, xa, -FLINT_BITS * n);
            arb_add_ui(xa, xa, 1, prec);
            arb_log(exact, xa, prec);
            u = ulp_error_vs_arb(res, n, exact, n);
            if (u > (double) FIXED_LOG1P_BITWISE_RS_MAX_ERR(n, rl))
                TEST_FUNCTION_FAIL("log: n = %wd, r = %d, ulp = %f\n",
                    n, rl, u);
        }

        /* ---- sin/cos: x = sum of rotation angles 2 atan(2^-i),
           sparse with possible repeats: the table stores the
           halves, so sum the entries and double ---- */
        if (FLINT_BITS == 64 || n >= 2)
        {
            ulong ys[53], yc[53];
            int rt = (iter % 4 == 0) ? 0
                     : 32 + (int) n_randint(state,
                        FLINT_MAX(1, imax - 32));
            int rteff = (rt == 0) ? 32 : rt;

            _fixed_atans_ensure(n, FLINT_MIN((slong) rteff + 2, imax));

            flint_mpn_zero(x, n);
            nsum = 1 + n_randint(state, 8);
            for (j = 0; j < nsum; j++)
            {
                slong i = 1 + n_randint(state,
                    FLINT_MIN((slong) rteff + 2, imax));
                nn_srcptr Ai = _fixed_atans_entry(i, n);

                cy = mpn_add_n(x, x, Ai, n);
                if (cy)
                    mpn_sub_n(x, x, Ai, n);
            }
            /* double to turn the summed half-angles into rotation
               angles (skip if it would overflow the fraction) */
            if (x[n - 1] >> (FLINT_BITS - 1) == 0)
                mpn_lshift(x, x, n, 1);

            k = (slong) n_randint(state, 7) - 3;
            if (k > 0)
            {
                cy = mpn_add_1(x, x, n, (ulong) k);
                if (cy)
                    mpn_sub_1(x, x, n, (ulong) k);
            }
            else if (k < 0)
            {
                cy = mpn_sub_1(x, x, n, (ulong) (-k));
                if (cy)
                    mpn_add_1(x, x, n, (ulong) (-k));
            }

            fixed_sin_cos_bitwise_rs(ys, yc, x, n, rt);

            fmpz_set_ui_array(f, x, n);
            arb_set_fmpz(xa, f);
            arb_mul_2exp_si(xa, xa, -FLINT_BITS * n);
            arb_sin(exact, xa, prec);
            u = ulp_error_vs_arb(ys, n + 1, exact, n);
            if (u > (double) FIXED_SIN_COS_BITWISE_RS_MAX_ERR(n, rt))
                TEST_FUNCTION_FAIL("sin: n = %wd, r = %d, ulp = %f\n",
                    n, rt, u);
            arb_cos(exact, xa, prec);
            u = ulp_error_vs_arb(yc, n + 1, exact, n);
            if (u > (double) FIXED_SIN_COS_BITWISE_RS_MAX_ERR(n, rt))
                TEST_FUNCTION_FAIL("cos: n = %wd, r = %d, ulp = %f\n",
                    n, rt, u);
        }

        /* ---- atan: t = Y/X of a sparse rotation product
           prod (1 + i 2^-i) applied to (1, 0), kept below the
           diagonal (Y < X, i.e. angle < pi/4), so the vectoring
           must rediscover the factor set ---- */
        {
            ulong WX[54], WY[54], va2[53], vb2[53], nd2[107], ya[53];
            slong wn2 = n + 1, tries;
            int rt = (iter % 4 == 0) ? 0
                     : 32 + (int) n_randint(state,
                        FLINT_MAX(1, imax - 32));

            flint_mpn_zero(WX, n);
            WX[n] = 1;
            flint_mpn_zero(WY, wn2);

            tries = 1 + (slong) n_randint(state, 10);
            for (j = 0; j < tries; j++)
            {
                slong i = 1 + n_randint(state, imax);
                slong qq = i / FLINT_BITS, b = i % FLINT_BITS;

                if (qq >= wn2)
                    continue;
                flint_mpn_zero(va2, n);
                flint_mpn_zero(vb2, n);
                if (b)
                {
                    mpn_rshift(va2, WX + qq, wn2 - qq, (int) b);
                    mpn_rshift(vb2, WY + qq, wn2 - qq, (int) b);
                }
                else
                {
                    flint_mpn_copyi(va2, WX + qq, wn2 - qq);
                    flint_mpn_copyi(vb2, WY + qq, wn2 - qq);
                }
                mpn_sub(WX, WX, wn2, vb2, n);
                mpn_add(WY, WY, wn2, va2, n);
                if (WY[n] != 0 || mpn_cmp(WY, WX, wn2) >= 0)
                {
                    /* crossed the diagonal: revert */
                    mpn_sub(WY, WY, wn2, va2, n);
                    mpn_add(WX, WX, wn2, vb2, n);
                }
            }

            /* t = floor(WY 2^(b n) / WX): WY < WX so t < 1.  In
               the ROTATION direction X decreases (|W| grows into
               the imaginary part), so WX < 1 is possible and its
               top limb may be zero: normalize the divisor length,
               since mpn_tdiv_qr requires a nonzero leading limb.
               (The vectoring loop inside fixed_atan_bitwise_rs has
               no such issue: there X only grows from 1.) */
            {
                slong dn = wn2;
                ulong qbuf[108];

                while (dn > 1 && WX[dn - 1] == 0)
                    dn--;

                flint_mpn_zero(nd2, n);
                flint_mpn_copyi(nd2 + n, WY, wn2);
                mpn_tdiv_qr(qbuf, nd2, 0, nd2, n + wn2, WX, dn);
                flint_mpn_zero(x, n);
                flint_mpn_copyi(x, qbuf, n);
            }

            k = (slong) n_randint(state, 7) - 3;
            if (k > 0)
            {
                cy = mpn_add_1(x, x, n, (ulong) k);
                if (cy)
                    mpn_sub_1(x, x, n, (ulong) k);
            }
            else if (k < 0)
            {
                cy = mpn_sub_1(x, x, n, (ulong) (-k));
                if (cy)
                    mpn_add_1(x, x, n, (ulong) (-k));
            }

            fixed_atan_bitwise_rs(ya, x, n, rt);

            fmpz_set_ui_array(f, x, n);
            arb_set_fmpz(xa, f);
            arb_mul_2exp_si(xa, xa, -FLINT_BITS * n);
            arb_atan(exact, xa, prec);
            u = ulp_error_vs_arb(ya, n, exact, n);
            if (u > (double) FIXED_ATAN_BITWISE_RS_MAX_ERR(n, rt))
                TEST_FUNCTION_FAIL("atan: n = %wd, r = %d, ulp = %f\n",
                    n, rt, u);
        }

        arb_clear(xa);
        arb_clear(exact);
        fmpz_clear(f);
    }

    /* cache lifecycle: the loop above has warmed both cached
       reduction tables; query the angle-table extent, free the
       tables through the public explicit-clear entry points, and
       verify that fresh calls rebuild them and still compute
       correctly; finish with flint_cleanup to run the registered
       thread-exit path over freshly rebuilt tables as well */
    {
        ulong x[3], y[4];
        arb_t xa, exact;
        fmpz_t f;
        double u;
        slong n = 3, prec = FLINT_BITS * n + 128;

        if (_fixed_atans_max_index() < 1)
            TEST_FUNCTION_FAIL("angle table not warmed\n");

        _fixed_exp_logs_clear();
        _fixed_atans_clear();

        arb_init(xa); arb_init(exact); fmpz_init(f);
        flint_mpn_rrandom(x, state, n);
        x[n - 1] >>= 1;

        fixed_exp_bitwise_rs(y, x, n, 64);
        fmpz_set_ui_array(f, x, n);
        arb_set_fmpz(xa, f);
        arb_mul_2exp_si(xa, xa, -FLINT_BITS * n);
        arb_exp(exact, xa, prec);
        u = ulp_error_vs_arb(y, n + 1, exact, n);
        if (u > (double) FIXED_EXP_BITWISE_RS_MAX_ERR(n, 64))
            TEST_FUNCTION_FAIL("exp after cleanup: ulp = %f\n", u);

        fixed_tan_bitwise_rs(y, x, n, 64);
        arb_tan(exact, xa, prec);
        u = ulp_error_vs_arb(y, n + 1, exact, n);
        if (u > (double) FIXED_TAN_BITWISE_RS_MAX_ERR(n, 64))
            TEST_FUNCTION_FAIL("tan after cleanup: ulp = %f\n", u);

        arb_clear(xa); arb_clear(exact); fmpz_clear(f);

        flint_cleanup();
    }

    TEST_FUNCTION_END(state);
}
