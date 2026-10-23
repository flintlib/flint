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
#include "fmpz_vec.h"
#include "arb.h"
#include "acb.h"
#include "dfloat.h"

/* The sums over j of the Platt multi-evaluation: _dfloat_platt_smk
   overlaps a direct evaluation in arb, with a small radius, and the
   results with a small table bound M (where most large j are split as
   j = x y with x, y <= M, or computed directly) overlap those with the
   default table; _dfloat_platt_smk_dd is the same without the factor
   exp(-i t0 log sqrt(pi)). */

/* the bucket boundaries, as get_smk_points in platt_multieval.c:
   ceil(exp(pi (2m - 1) / B) / sqrt(pi)) */
static void
_t_smk_points(fmpz * res, slong A, slong B)
{
    slong m, N = A * B, prec = 64;
    arb_t x, u, v;

    arb_init(x);
    arb_init(u);
    arb_init(v);
    for (m = 0; m < N; m++)
    {
        while (1)
        {
            arb_const_pi(u, prec);
            arb_div_si(u, u, B, prec);
            arb_const_sqrt_pi(v, prec);
            arb_inv(v, v, prec);
            arb_set_si(x, 2 * m - 1);
            arb_mul(x, x, u, prec);
            arb_exp(x, x, prec);
            arb_mul(x, x, v, prec);
            arb_ceil(x, x, prec);
            if (arb_get_unique_fmpz(res + m, x))
                break;
            prec *= 2;
        }
    }
    arb_clear(x);
    arb_clear(u);
    arb_clear(v);
}

/* S[k N + m] = sum over the j <= J of bucket m of z_j b_j^k, z_j =
   j^(-1/2) exp(-i t0 log(j sqrt(pi))), b_j = log(j sqrt(pi))/(2 pi) - m/B */
static void
_t_smk_ref(acb_ptr S, const fmpz * pts, const arb_t t0, slong A, slong B,
    ulong J, slong K, slong prec)
{
    slong N = A * B, m = 0, k;
    ulong j;
    arb_t L, ls, tp, b, r, ph;
    acb_t z, w;

    arb_init(L);
    arb_init(ls);
    arb_init(tp);
    arb_init(b);
    arb_init(r);
    arb_init(ph);
    acb_init(z);
    acb_init(w);

    arb_const_sqrt_pi(ls, prec);
    arb_log(ls, ls, prec);
    arb_const_pi(tp, prec);
    arb_mul_2exp_si(tp, tp, 1);
    _acb_vec_zero(S, K * N);

    for (j = 1; j <= J; j++)
    {
        while (m + 1 < N && fmpz_cmp_ui(pts + m + 1, j) <= 0)
            m++;
        arb_log_ui(L, j, prec);
        arb_add(L, L, ls, prec);
        arb_mul(ph, t0, L, prec);
        arb_neg(ph, ph);
        arb_sin_cos(acb_imagref(z), acb_realref(z), ph, prec);
        arb_rsqrt_ui(r, j, prec);
        acb_mul_arb(z, z, r, prec);
        arb_div(b, L, tp, prec);
        arb_set_si(r, m);
        arb_div_si(r, r, B, prec);
        arb_sub(b, b, r, prec);
        acb_set(w, z);
        for (k = 0; k < K; k++)
        {
            acb_add(S + k * N + m, S + k * N + m, w, prec);
            acb_mul_arb(w, w, b, prec);
        }
    }

    arb_clear(L);
    arb_clear(ls);
    arb_clear(tp);
    arb_clear(b);
    arb_clear(r);
    arb_clear(ph);
    acb_clear(z);
    acb_clear(w);
}

TEST_FUNCTION_START(platt_smk, state)
{
    slong iter;

    DFLOAT_TEST_SKIP_IF_UNSUPPORTED(state);

    for (iter = 0; iter < 20 * flint_test_multiplier(); iter++)
    {
        slong A, B, N, K, prec, i;
        ulong J, M;
        int pn, ref, ok, ok2, ok3;
        fmpz * pts;
        fmpz_t T;
        arb_t t0;
        acb_ptr S1, S2, S3;
        double * S5, W;

        /* (A >= 2: the buckets cover j up to exp(2 pi A)/sqrt(pi) > J,
           so that |b_j| <= 1/(2B), and the radius is small) */
        A = 2 + n_randint(state, 3);
        B = 1 + n_randint(state, 32);
        N = A * B;
        K = 1 + n_randint(state, 20);
        pn = 3 + n_randint(state, 2);
        ref = (iter % 2 == 0);
        J = ref ? 1 + n_randint(state, 1500) : 1 + n_randint(state, 30000);
        /* a table bound from sqrt(J) (the smallest allowed) to J/2 */
        M = 1 + n_randint(state, J / 2 + 1);

        fmpz_init(T);
        arb_init(t0);
        fmpz_randbits_unsigned(T, state, 1 + n_randint(state, 45));
        arb_set_fmpz(t0, T);
        if (n_randint(state, 2))
        {
            arb_set_si(t0, n_randint(state, 16));
            arb_mul_2exp_si(t0, t0, -4);
            arb_add_fmpz(t0, t0, T, 128);
        }
        prec = 128 + fmpz_bits(T);

        pts = _fmpz_vec_init(N);
        _t_smk_points(pts, A, B);
        S1 = _acb_vec_init(K * N);
        S2 = _acb_vec_init(K * N);
        S3 = _acb_vec_init(K * N);
        S5 = flint_malloc(sizeof(double) * 5 * K * N);

        ok = _dfloat_platt_smk(S1, pts, t0, A, B, J, K, pn, 0, prec);
        ok2 = _dfloat_platt_smk(S2, pts, t0, A, B, J, K, pn, M, prec);
        ok3 = _dfloat_platt_smk_dd(S5, pts, t0, A, B, J, K, pn, M);

        if (!ok || !ok2 || !ok3)
        {
            flint_printf("FAIL: unsupported\n");
            flint_abort();
        }

        if (ref)
            _t_smk_ref(S3, pts, t0, A, B, J, K, prec);

        /* the radius, against W >= the sum of j^(-1/2) */
        W = 2.0 + 2.0 * sqrt((double) J);

        for (i = 0; i < K * N; i++)
        {
            mag_t rad;
            int bad = !acb_overlaps(S1 + i, S2 + i);

            if (ref)
                bad = bad || !acb_overlaps(S1 + i, S3 + i);

            /* triple-double and double-double moments (k < 14), then
               double moments */
            mag_init(rad);
            mag_max(rad, arb_radref(acb_realref(S1 + i)), arb_radref(acb_imagref(S1 + i)));
            bad = bad || mag_cmp_2exp_si(rad, ((i / N < 14) ? -90 : -40) + (slong) log2(W)) > 0;
            mag_max(rad, arb_radref(acb_realref(S2 + i)), arb_radref(acb_imagref(S2 + i)));
            bad = bad || mag_cmp_2exp_si(rad, ((i / N < 14) ? -90 : -40) + (slong) log2(W)) > 0;
            mag_clear(rad);

            if (bad)
            {
                flint_printf("FAIL: platt_smk (entry k = %wd, m = %wd)\n", i / N, i % N);
                flint_printf("A = %wd, B = %wd, K = %wd, J = %wu, M = %wu, pn = %d\n", A, B, K, J, M, pn);
                flint_printf("t0 = "); arb_printd(t0, 30); flint_printf("\n");
                flint_printf("S1 = "); acb_printd(S1 + i, 30); flint_printf("\n");
                flint_printf("S2 = "); acb_printd(S2 + i, 30); flint_printf("\n");
                if (ref)
                {
                    flint_printf("S3 = "); acb_printd(S3 + i, 30); flint_printf("\n");
                }
                flint_abort();
            }
        }

        /* the compact version: S2 = c S5, c = exp(-i t0 log sqrt(pi)) */
        {
            acb_t c, x;
            arb_t w;
            mag_t r;

            acb_init(c);
            acb_init(x);
            arb_init(w);
            mag_init(r);
            arb_const_sqrt_pi(w, prec);
            arb_log(w, w, prec);
            arb_mul(w, w, t0, prec);
            arb_neg(w, w);
            arb_sin_cos(acb_imagref(c), acb_realref(c), w, prec);

            for (i = 0; i < K * N; i++)
            {
                const double * e = S5 + 5 * i;

                arb_set_d(acb_realref(x), e[0]);
                arb_set_d(w, e[1]);
                arb_add(acb_realref(x), acb_realref(x), w, prec);
                arb_set_d(acb_imagref(x), e[2]);
                arb_set_d(w, e[3]);
                arb_add(acb_imagref(x), acb_imagref(x), w, prec);
                mag_set_d(r, e[4]);
                arb_add_error_mag(acb_realref(x), r);
                arb_add_error_mag(acb_imagref(x), r);
                acb_mul(x, x, c, prec);

                if (!acb_overlaps(x, S2 + i))
                {
                    flint_printf("FAIL: platt_smk_dd (entry %wd)\n", i);
                    flint_printf("x = "); acb_printd(x, 30); flint_printf("\n");
                    flint_printf("S2 = "); acb_printd(S2 + i, 30); flint_printf("\n");
                    flint_abort();
                }
            }

            acb_clear(c);
            acb_clear(x);
            arb_clear(w);
            mag_clear(r);
        }

        /* only P = 3 and 4 */
        if (_dfloat_platt_smk_dd(S5, pts, t0, A, B, J, K, 2, 0) != 0)
        {
            flint_printf("FAIL: P = 2\n");
            flint_abort();
        }

        fmpz_clear(T);
        arb_clear(t0);
        _fmpz_vec_clear(pts, N);
        _acb_vec_clear(S1, K * N);
        _acb_vec_clear(S2, K * N);
        _acb_vec_clear(S3, K * N);
        flint_free(S5);
    }

    TEST_FUNCTION_END(state);
}
