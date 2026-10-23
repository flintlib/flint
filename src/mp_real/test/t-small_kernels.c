/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "arb.h"
#include "mp_real.h"
#include "mp_real/impl.h"

/* The table-driven kernels at two and three limbs (small.h) against arb,
   within their stated error bounds, with arguments at the table
   boundaries and, for atan(Y / X), ratios close to multiples of 2^-9;
   the paired exp(v), exp(-v); and the three-limb divisions against
   mpn_tdiv_qr. */

#if FLINT_BITS == 64

/* r = (x, len) B^-frac */
static void
_small_get_arf(arf_t r, nn_srcptr x, slong len, slong frac)
{
    fmpz_t t;
    fmpz_init(t);
    fmpz_set_ui_array(t, x, len);
    arf_set_fmpz(r, t);
    arf_mul_2exp_si(r, r, -FLINT_BITS * frac);
    fmpz_clear(t);
}

/* |ref - (y, len) B^-n| <= err B^-n */
static int
_small_check(const arb_t ref, nn_srcptr y, slong len, slong n, ulong err)
{
    arb_t t;
    arf_t f;
    mag_t m;
    int ok;

    arb_init(t);
    arf_init(f);
    mag_init(m);
    _small_get_arf(f, y, len, n);
    arb_sub_arf(t, ref, f, 4 * FLINT_BITS * n);
    arb_mul_2exp_si(t, t, FLINT_BITS * n);
    arb_get_mag(m, t);
    ok = (mag_cmp_2exp_si(m, 62) < 0) && (mag_get_d(m) <= (double) err);
    arb_clear(t);
    arf_clear(f);
    mag_clear(m);
    return ok;
}

static void
_small_randtest(nn_ptr v, slong n, flint_rand_t state)
{
    slong k;
    int r = n_randint(state, 6);

    for (k = 0; k < n; k++)
        v[k] = n_randtest(state);

    if (r == 0)
    {
        /* close to a table boundary j 2^-7, 2^-9, 2^-14 or 2^-21 */
        static const int bits[4] = { 7, 9, 14, 21 };
        int b = bits[n_randint(state, 4)];
        ulong mask = (UWORD(1) << (FLINT_BITS - b)) - 1;
        if (n_randint(state, 2))
        {
            v[n - 1] &= ~mask;
            for (k = 0; k < n - 1; k++)
                v[k] = n_randint(state, 4);
        }
        else
        {
            v[n - 1] |= mask;
            for (k = 0; k < n - 1; k++)
                v[k] = UWORD_MAX - n_randint(state, 4);
        }
    }
    else if (r == 1)
    {
        /* small */
        v[n - 1] >>= n_randint(state, FLINT_BITS);
    }
}

#endif

TEST_FUNCTION_START(mp_real_small_kernels, state)
{
#if FLINT_BITS == 64
    slong iter;

    for (iter = 0; iter < 5000 * flint_test_multiplier(); iter++)
    {
        ulong v[3], x[3], y[5], y2[5], err;
        slong n = 1 + n_randint(state, 3), prec, k;
        int which = n_randint(state, 7);
        arb_t a, b, ref;

        /* the one-limb kernels: exp, log1p, sin_cos and atan */
        if (which >= 4 && n == 1)
            n = 2;

        prec = 4 * FLINT_BITS * n;
        arb_init(a);
        arb_init(b);
        arb_init(ref);

        _small_randtest(v, n, state);
        _small_get_arf(arb_midref(a), v, n, n);

        if (which == 0)
        {
            if (n == 1)
                _mp_real_small_exp_1_mpn(y, &err, v);
            else if (n == 2)
                _mp_real_small_exp_2_mpn(y, &err, v);
            else
                _mp_real_small_exp_3_mpn(y, &err, v);
            arb_exp(ref, a, prec);
            if (!_small_check(ref, y, n + 1, n, err))
                TEST_FUNCTION_FAIL("exp, n = %wd, v = %{ulong*}\n", n, v, n);
        }
        else if (which == 1)
        {
            if (n == 1)
                _mp_real_small_log1p_1_mpn(y, &err, v);
            else if (n == 2)
                _mp_real_small_log1p_2_mpn(y, &err, v);
            else
                _mp_real_small_log1p_3_mpn(y, &err, v);
            arb_log1p(ref, a, prec);
            if (!_small_check(ref, y, n, n, err))
                TEST_FUNCTION_FAIL("log1p, n = %wd, v = %{ulong*}\n", n, v, n);
        }
        else if (which == 2)
        {
            if (n == 1)
                _mp_real_small_sin_cos_1_mpn(y, y2, &err, v);
            else if (n == 2)
                _mp_real_small_sin_cos_2_mpn(y, y2, &err, v);
            else
                _mp_real_small_sin_cos_3_mpn(y, y2, &err, v);
            arb_sin_cos(ref, b, a, prec);
            if (!_small_check(ref, y, n + 1, n, err) || !_small_check(b, y2, n + 1, n, err))
                TEST_FUNCTION_FAIL("sin_cos, n = %wd, v = %{ulong*}\n", n, v, n);
        }
        else if (which == 3)
        {
            if (n == 1)
                _mp_real_small_atan_1_mpn(y, &err, v);
            else if (n == 2)
                _mp_real_small_atan_2_mpn(y, &err, v);
            else
                _mp_real_small_atan_3_mpn(y, &err, v);
            arb_atan(ref, a, prec);
            if (!_small_check(ref, y, n, n, err))
                TEST_FUNCTION_FAIL("atan, n = %wd, v = %{ulong*}\n", n, v, n);
        }
        else if (which == 6)
        {
            /* exp(v) and exp(-v) together */
            if (n == 2)
                _mp_real_small_exp_pm_2(y + 2, y + 1, y, y2 + 2, y2 + 1, y2, v[1], v[0]);
            else
                _mp_real_small_exp_pm_3(y + 3, y + 2, y + 1, y, y2 + 3, y2 + 2, y2 + 1, y2, v[2], v[1], v[0]);
            err = (n == 2) ? MP_REAL_SMALL_EXP_ERR_2 : MP_REAL_SMALL_EXP_ERR_3;
            arb_exp(ref, a, prec);
            arb_neg(b, a);
            arb_exp(b, b, prec);
            if (!_small_check(ref, y, n + 1, n, err) || !_small_check(b, y2, n + 1, n, err))
                TEST_FUNCTION_FAIL("exp_pm, n = %wd, v = %{ulong*}\n", n, v, n);
        }
        else if (which == 5)
        {
            /* the 3-limb division: floor(N B^3 / D) exactly for N < D */
            ulong N[6], D[3], Q[4], R[3], w2, w1, w0;

            for (k = 0; k < 3; k++)
            {
                D[k] = n_randtest(state);
                N[k + 3] = n_randtest(state);
                N[k] = 0;
            }
            D[2] |= UWORD(1) << (FLINT_BITS - 1);
            if (n_randint(state, 4) == 0)
            {
                /* N = D - small */
                N[3] = D[0]; N[4] = D[1]; N[5] = D[2];
                mpn_sub_1(N + 3, N + 3, 3, 1 + n_randint(state, 3));
            }
            if (mpn_cmp(N + 3, D, 3) >= 0)
                mpn_rshift(N + 3, N + 3, 3, 1);

            _mp_real_small_divq_3(&w2, &w1, &w0, N[5], N[4], N[3], D[2], D[1], D[0]);
            mpn_tdiv_qr(Q, R, 0, N, 6, D, 3);
            if (Q[3] != 0 || Q[2] != w2 || Q[1] != w1 || Q[0] != w0)
                TEST_FUNCTION_FAIL("divq_3, N = %{ulong*}, D = %{ulong*}\n", N + 3, 3, D, 3);

            /* _mp_real_divq_4_3z: floor((y1, y0, 0, 0) / S), s2 != 0 */
            {
                ulong S[3], Y[4], P[4], q[2], t, c;
                int tries;

                S[0] = n_randtest(state) | 1;
                S[1] = n_randtest(state);
                S[2] = n_randtest_not_zero(state);
                Y[0] = n_randtest(state);
                Y[1] = n_randtest(state);

                /* (y1, y0) B^2 = (t - 1) S B + (S - c) after normalization
                   (t S = c mod B, t S < B^3): the second quotient limb is
                   B - 1 with a partial remainder equal to D in its top
                   two limbs, the special case of the 3/2 division */
                for (tries = (n_randint(state, 2) ? 0 : 20); tries < 20; tries++)
                {
                    S[2] = n_randtest_not_zero(state);
                    c = 1 + n_randint(state, 3);
                    if (c > S[0])
                        break;
                    t = c * n_binvert(S[0]);
                    if (t <= 1)
                        continue;
                    P[3] = mpn_mul_1(P, S, 3, t);
                    mpn_sub_1(P, P, 4, c);
                    if (P[3] == 0 && P[0] == 0)
                    {
                        Y[0] = P[1];
                        Y[1] = P[2];
                        break;
                    }
                }

                _mp_real_divq_4_3z(q, Y[1], Y[0], S[2], S[1], S[0]);
                N[0] = N[1] = 0;
                N[2] = Y[0];
                N[3] = Y[1];
                mpn_tdiv_qr(Q, R, 0, N, 4, S, 3);
                if (Q[1] != q[1] || Q[0] != q[0])
                    TEST_FUNCTION_FAIL("divq_4_3z, y = %{ulong*}, S = %{ulong*}\n", Y, 2, S, 3);
            }
        }
        else
        {
            /* atan(Y / X), X >= 1/2, Y <= X */
            for (k = 0; k < n; k++)
                x[k] = n_randtest(state);
            x[n - 1] |= UWORD(1) << (FLINT_BITS - 1);

            switch (n_randint(state, 4))
            {
                case 0:
                    /* Y = X */
                    for (k = 0; k < n; k++)
                        v[k] = x[k];
                    break;
                case 1:
                {
                    /* Y close to j X / 512 */
                    ulong c[4];
                    c[n] = mpn_mul_1(c, x, n, n_randint(state, 513));
                    mpn_rshift(c, c, n + 1, 9);
                    for (k = 0; k < n; k++)
                        v[k] = c[k];
                    if (n_randint(state, 2))
                        mpn_add_1(v, v, n, n_randint(state, 3));
                    else
                        mpn_sub_1(v, v, n, n_randint(state, 3));
                    break;
                }
                default:
                    break;
            }

            if (mpn_cmp(v, x, n) > 0)
            {
                for (k = 0; k < n; k++)
                    v[k] = x[k] - (k == 0 ? n_randint(state, 2) : 0);
                if (mpn_cmp(v, x, n) > 0)
                    v[0] = x[0];
            }

            if (n == 2)
            {
                _mp_real_small_atan_ratio_2(y + 1, y, v[1], v[0], x[1], x[0]);
                err = MP_REAL_SMALL_ATAN_ERR_2;
            }
            else
            {
                _mp_real_small_atan_ratio_3(y + 2, y + 1, y, v[2], v[1], v[0], x[2], x[1], x[0]);
                err = MP_REAL_SMALL_ATAN_ERR_3;
            }

            _small_get_arf(arb_midref(a), v, n, n);
            _small_get_arf(arb_midref(b), x, n, n);
            arb_div(ref, a, b, prec);
            arb_atan(ref, ref, prec);
            if (!_small_check(ref, y, n, n, err))
                TEST_FUNCTION_FAIL("atan ratio, n = %wd, Y = %{ulong*}, X = %{ulong*}\n", n, v, n, x, n);
        }

        arb_clear(a);
        arb_clear(b);
        arb_clear(ref);
    }

    TEST_FUNCTION_END(state);
#else
    TEST_FUNCTION_END_SKIPPED(state);
#endif
}
