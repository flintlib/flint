/*
    Copyright (C) 2021 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb_hypgeom.h"
#include "bernoulli.h"

void
arb_hypgeom_gamma_stirling_sum_improved(arb_t s, const arb_t z, slong N, slong K, slong prec)
{
    arb_t b, t, zinv, w, u, S2, S3, S4;
    mag_t zinv_mag, err;
    slong n, term_mag, term_prec;
    slong * term_mags;
    slong i, j, k, M;
    slong * Mk;
    arb_ptr upow, ukpow;
    slong kodd, kpow_exp, wp;
    fmpz_t kpow;
    slong m;

    if (N <= 1)
    {
        arb_zero(s);
        return;
    }

    if (N == 2)
    {
        arb_mul_ui(s, z, 12, prec);
        arb_inv(s, s, prec);
        return;
    }

    if (K == 0)
    {
        if (prec <= 128)
            K = 1;
        else if (prec <= 1024)
            K = 2;
        else
        {
            K = 4 + 0.1 * sqrt(FLINT_MAX(prec - 4096, 0));
            K = FLINT_MIN(K, 100);
        }
    }

    arb_init(b);
    arb_init(t);
    arb_init(zinv);
    arb_init(w);
    arb_init(u);
    arb_init(S2);
    arb_init(S3);
    arb_init(S4);
    mag_init(zinv_mag);
    mag_init(err);

    arb_inv(zinv, z, prec);
    arb_mul(w, zinv, zinv, prec);

    arb_get_mag(zinv_mag, zinv);

    term_mags = flint_malloc(sizeof(ulong) * (N + 1));

    /* Avoid possible overflow. (This case is not interesting to handle well,
       but we need to handle it.) */
    if (mag_cmp_2exp_si(zinv_mag, 10) > 0)
        K = 1;

    /* Avoid possible underflow. */
    if (mag_cmp_2exp_si(zinv_mag, -100000) < 0)
        mag_set_ui_2exp_si(zinv_mag, 1, -100000);

    _arb_hypgeom_gamma_stirling_term_bounds(term_mags, zinv_mag, N + 1);

    /* We will assume below that term_mags is nonincreasing */
    for (n = N - 2; n >= 1; n--)
    {
        if (term_mags[n] < term_mags[n + 1])
            term_mags[n] = term_mags[n + 1];
    }

    Mk = NULL;

    if (K > 1)
    {
        Mk = flint_malloc(sizeof(slong) * (K + 1));
        Mk[0] = Mk[1] = N;

        for (k = 2; k <= K; k++)
        {
            double log2_k;

            Mk[k] = N;

            log2_k = log(k) * (1.0 / 0.693147180559945309);

            while (Mk[k] > 2)
            {
                slong err, Mnew;

                Mnew = Mk[k] - 1;
                err = term_mags[Mnew] - log2_k * (2 * Mnew) + FLINT_BIT_COUNT(N - Mnew);

                if (err < -prec)
                    Mk[k] = Mnew;
                else
                    break;
            }
        }

        /* We will assume that Mk is nonincreasing. */
        for (k = 2; k <= K; k++)
        {
            if (Mk[k] > Mk[k - 1])
                Mk[k] = Mk[k - 1];
        }

        while (K >= 2 && Mk[K] == Mk[K - 1])
            K--;

        M = Mk[K];

        mag_hurwitz_zeta_uiui(err, 2 * M, K);
        mag_mul_ui(err, err, N - M);
        mag_mul_2exp_si(err, err, term_mags[M]);

        for (k = 1; k < K; k++)
        {
            mag_t t;
            mag_init(t);
            mag_set_ui_lower(t, k);
            mag_inv(t, t);
            mag_pow_ui(t, t, 2 * Mk[k]);
            mag_mul_ui(t, t, N - Mk[k]);
            mag_mul_2exp_si(t, t, term_mags[Mk[k]]);
            mag_add(err, err, t);
            mag_clear(t);
        }
    }
    else
    {
        M = N;
    }

    m = sqrt(N - M);
    m = FLINT_MAX(m, 1);

    /* S3 precision */
    wp = prec + term_mags[M];
    wp = FLINT_MIN(wp, prec);
    wp = FLINT_MAX(wp, 10);

    arb_zero(S3);

    /* todo: could avoid one mul + div with precomputation here */
    /* u = -1 / (2 pi z)^2 */
    arb_const_pi(u, wp);
    arb_mul(u, u, z, wp);
    arb_mul_2exp_si(u, u, 1);
    arb_mul(u, u, u, wp);
    arb_inv(u, u, wp);
    arb_neg(u, u);

    fmpz_init(kpow);
    upow = _arb_vec_init(m + 1);
    ukpow = _arb_vec_init(m + 1);

    _arb_vec_set_powers(upow, u, m + 1, wp);

    for (kodd = 1; kodd < K; kodd += 2)
    {
        for (k = kodd; k < K; k *= 2)
        {
            if (k == 1)
            {
                _arb_vec_set(ukpow, upow, m + 1);
                fmpz_one(kpow);
                kpow_exp = 0;
            }
            else if (k == kodd)
            {
                arb_set(ukpow + 0, upow + 0);
                for (j = 1; j <= FLINT_MIN(Mk[k] - M - 1, m); j++)
                {
                    if (j == 1)
                        fmpz_set_ui(kpow, k * k);
                    else
                        fmpz_mul_ui(kpow, kpow, k * k);
                    arb_div_fmpz(ukpow + j, upow + j, kpow, wp);
                }

                /* set kpow = k^(2M) */
                fmpz_ui_pow_ui(kpow, k * k, M);
                kpow_exp = 0;
            }
            else
            {
                /* compute x / k^(2j) given x / (k/2)^(2j) */
                for (j = 1; j <= FLINT_MIN(Mk[k] - M - 1, m); j++)
                    arb_mul_2exp_si(ukpow + j, ukpow + j, -2 * j);
                kpow_exp += 2 * M;
            }

            arb_zero(S4);

            for (n = Mk[k] - 1; n >= M; n--)
            {
                i = n - M;

                term_prec = prec + term_mags[n];
                term_prec = FLINT_MIN(wp, prec);
                term_prec = FLINT_MAX(term_prec, 10);

                /* note: set_round makes small difference here */
                arb_fma_ui(S4, S4, (2 * n) * (2 * n - 1), ukpow + i % m, term_prec);

                if (i != 0 && i % m == 0)
                {
                    /* note: set_round makes small difference here */
                    arb_mul(S4, S4, ukpow + m, term_prec);
                }
            }

            /* divide by k^(2M) */
            if (k != 1)
            {
                if (!fmpz_is_one(kpow))
                    arb_div_fmpz(S4, S4, kpow, wp);

                arb_mul_2exp_si(S4, S4, -kpow_exp);
            }

            arb_add(S3, S3, S4, wp);
        }
    }

    /* multiply by -2 (2M-2)! u^M z */
    fmpz_fac_ui(kpow, 2 * M - 2);
    arb_mul_fmpz(S3, S3, kpow, wp);
    arb_pow_ui(u, u, M, wp);
    arb_mul(S3, S3, u, wp);
    arb_set_round(t, z, wp);
    arb_mul(S3, S3, t, wp);
    arb_mul_2exp_si(S3, S3, 1);
    arb_neg(S3, S3);

    arb_add_error_mag(S3, err);

    _arb_vec_clear(upow, m + 1);
    _arb_vec_clear(ukpow, m + 1);
    fmpz_clear(kpow);

    arb_zero(S2);

    m = sqrt(M);
    upow = _arb_vec_init(m + 1);
    _arb_vec_set_powers(upow, w, m + 1, prec);

    BERNOULLI_ENSURE_CACHED(2 * M - 2);

    {
        fmpz_t d, e, f, g, h;
        fmpz q[4];

        fmpz_init(d);
        fmpz_init(e);
        fmpz_init(f);
        fmpz_init(g);
        fmpz_init(h);

        fmpz_init(q);
        fmpz_init(q + 1);
        fmpz_init(q + 2);
        fmpz_init(q + 3);

        for (n = M - 1; n >= 1; n--)
        {
            i = n - 1;

            if (i >= 4 && i % m >= 3 && prec >= 512 && prec <= 4096)
            {
                term_mag = term_mags[n - 3];
                term_prec = prec + term_mag;
                term_prec = FLINT_MIN(term_prec, prec);
                term_prec = FLINT_MAX(term_prec, 10);

                fmpz_mul_ui(d, fmpq_denref(bernoulli_cache + 2 * (n - 0)), 2 * (n - 0) * (2 * (n - 0) - 1));
                fmpz_mul_ui(e, fmpq_denref(bernoulli_cache + 2 * (n - 1)), 2 * (n - 1) * (2 * (n - 1) - 1));
                fmpz_mul_ui(g, fmpq_denref(bernoulli_cache + 2 * (n - 2)), 2 * (n - 2) * (2 * (n - 2) - 1));
                fmpz_mul_ui(h, fmpq_denref(bernoulli_cache + 2 * (n - 3)), 2 * (n - 3) * (2 * (n - 3) - 1));

                /* q3 = egh   q2 = dgh   q1 = deh   q0 = deg   d = degh */
                fmpz_mul(q + 3, e, g);
                fmpz_mul(q + 0, q + 3, d);
                fmpz_mul(q + 3, q + 3, h);
                fmpz_mul(q + 2, d, h);
                fmpz_mul(q + 1, q + 2, e);
                fmpz_mul(q + 2, q + 2, g);
                fmpz_mul(d, q + 3, d);

                fmpz_mul(q + 3, q + 3, fmpq_numref(bernoulli_cache + 2 * (n - 0)));
                fmpz_mul(q + 2, q + 2, fmpq_numref(bernoulli_cache + 2 * (n - 1)));
                fmpz_mul(q + 1, q + 1, fmpq_numref(bernoulli_cache + 2 * (n - 2)));
                fmpz_mul(q + 0, q + 0, fmpq_numref(bernoulli_cache + 2 * (n - 3)));

                arb_dot_fmpz(t, NULL, 0, upow + i % m - 3, 1, q, 1, 4, term_prec);
                arb_div_fmpz(t, t, d, term_prec);
                arb_add(S2, S2, t, term_prec);

                n -= 3;
                i -= 3;
            }
            else if (i >= 3 && i % m >= 2)
            {
                term_mag = term_mags[n - 2];
                term_prec = prec + term_mag;
                term_prec = FLINT_MIN(term_prec, prec);
                term_prec = FLINT_MAX(term_prec, 10);

                fmpz_mul_ui(d, fmpq_denref(bernoulli_cache + 2 * (n - 0)), 2 * (n - 0) * (2 * (n - 0) - 1));
                fmpz_mul_ui(e, fmpq_denref(bernoulli_cache + 2 * (n - 1)), 2 * (n - 1) * (2 * (n - 1) - 1));
                fmpz_mul_ui(g, fmpq_denref(bernoulli_cache + 2 * (n - 2)), 2 * (n - 2) * (2 * (n - 2) - 1));

                fmpz_mul(q + 2, e, g);
                fmpz_mul(q + 2, q + 2, fmpq_numref(bernoulli_cache + 2 * (n - 0)));
                fmpz_mul(q + 1, d, g);
                fmpz_mul(q + 1, q + 1, fmpq_numref(bernoulli_cache + 2 * (n - 1)));
                fmpz_mul(q + 0, d, e);
                fmpz_mul(d, q + 0, g);
                fmpz_mul(q + 0, q + 0, fmpq_numref(bernoulli_cache + 2 * (n - 2)));

                arb_dot_fmpz(t, NULL, 0, upow + i % m - 2, 1, q, 1, 3, term_prec);
                arb_div_fmpz(t, t, d, term_prec);
                arb_add(S2, S2, t, term_prec);

                n -= 2;
                i -= 2;
            }
            else if (i >= 1 && i % m >= 1)
            {
                term_mag = term_mags[n - 1];
                term_prec = prec + term_mag;
                term_prec = FLINT_MIN(term_prec, prec);
                term_prec = FLINT_MAX(term_prec, 10);

                fmpz_mul_ui(d, fmpq_denref(bernoulli_cache + 2 * (n - 0)), 2 * (n - 0) * (2 * (n - 0) - 1));
                fmpz_mul_ui(e, fmpq_denref(bernoulli_cache + 2 * (n - 1)), 2 * (n - 1) * (2 * (n - 1) - 1));

                fmpz_mul(f, fmpq_numref(bernoulli_cache + 2 * (n - 0)), e);
                arb_set_round(u, upow + i % m, term_prec);
                arb_mul_fmpz(t, u, f, term_prec);                             /* todo: output-sensitive mul */
                fmpz_mul(f, fmpq_numref(bernoulli_cache + 2 * (n - 1)), d);
                arb_set_round(u, upow + i % m - 1, term_prec);
                arb_mul_fmpz(u, u, f, term_prec);                             /* todo: output-sensitive mul */
                arb_add(t, t, u, term_prec);

                fmpz_mul(d, d, e);
                arb_div_fmpz(t, t, d, term_prec);
                arb_add(S2, S2, t, term_prec);

                n--;
                i--;
            }
            else
            {
                term_mag = term_mags[n];
                term_prec = prec + term_mag;
                term_prec = FLINT_MIN(term_prec, prec);
                term_prec = FLINT_MAX(term_prec, 10);

                arb_set_round(u, upow + i % m, term_prec);
                arb_mul_fmpz(t, u, fmpq_numref(bernoulli_cache + 2 * n), term_prec);   /* todo: output-sensitive mul */
                fmpz_mul_ui(d, fmpq_denref(bernoulli_cache + 2 * n), 2 * n * (2 * n - 1));
                arb_div_fmpz(t, t, d, term_prec);

                arb_add(S2, S2, t, term_prec);
            }

            if (i != 0 && i % m == 0)
            {
                arb_set_round(u, upow + m, term_prec);   /* todo: output-sensitive mul */
                arb_mul(S2, S2, u, term_prec);
            }
        }

        fmpz_clear(q);
        fmpz_clear(q + 1);
        fmpz_clear(q + 2);
        fmpz_clear(q + 3);

        fmpz_clear(d);
        fmpz_clear(e);
        fmpz_clear(f);
        fmpz_clear(g);
        fmpz_clear(h);
    }

    _arb_vec_clear(upow, m + 1);

    arb_mul(S2, S2, zinv, prec);
    arb_add(s, S2, S3, prec);

    flint_free(term_mags);

    if (Mk != NULL)
        flint_free(Mk);

    arb_clear(b);
    arb_clear(t);
    arb_clear(zinv);
    arb_clear(w);
    arb_clear(u);
    arb_clear(S2);
    arb_clear(S3);
    arb_clear(S4);
    mag_clear(zinv_mag);
    mag_clear(err);
}

