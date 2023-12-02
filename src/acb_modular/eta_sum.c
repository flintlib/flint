/*
    Copyright (C) 2014 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb.h"
#include "acb_modular.h"

static const int pentagonal_best_m[] = {
  2, 5, 7, 11, 13, 17, 19, 23, 25, 35,
  55, 65, 77, 91, 119, 133, 143, 175, 275, 325,
  385, 455, 595, 665, 715, 935, 1001, 1309, 1463, 1547,
  1729, 1925, 2275, 2975, 3325, 3575, 4675, 5005, 6545, 7315,
  7735, 8645, 10465, 11305, 12155, 13585, 16445, 17017, 19019, 23023,
  24871, 25025, 32725, 36575, 38675, 43225, 52325, 56525, 60775, 67925,
  82225, 85085, 95095, 115115, 124355, 145145, 146965, 168245, 177905, 198835,
  224315, 230945, 279565, 312455, 323323, 391391, 425425, 475475, 575575, 621775,
  725725, 734825, 841225, 889525, 994175, 0
};

static const int pentagonal_best_m_residues[] = {
  2, 3, 4, 6, 7, 9, 10, 12, 11, 12,
  18, 21, 24, 28, 36, 40, 42, 44, 66, 77,
  72, 84, 108, 120, 126, 162, 168, 216, 240, 252,
  280, 264, 308, 396, 440, 462, 594, 504, 648, 720,
  756, 840, 1008, 1080, 1134, 1260, 1512, 1512, 1680, 2016,
  2160, 1848, 2376, 2640, 2772, 3080, 3696, 3960, 4158, 4620,
  5544, 4536, 5040, 6048, 6480, 7560, 7560, 8640, 9072, 10080,
  11340, 11340, 13608, 15120, 15120, 18144, 16632, 18480, 22176, 23760,
  27720, 27720, 31680, 33264, 36960, 0
};

slong acb_modular_rs_optimal_m(const int * best_ms, const int * num_residues, slong N);

#define PENTAGONAL(N) ((((N)+2)/2) * ((3*(N)+5)/2)/2)

void
_acb_modular_mul(acb_t z, acb_t tmp1, acb_t tmp2, const acb_t x, const acb_t y, slong wprec, slong prec)
{
    if (prec <= 1024)
    {
        acb_mul(z, x, y, wprec);
    }
    else if (x == y)
    {
        acb_set_round(tmp1, x, wprec);
        acb_mul(z, tmp1, tmp1, wprec);
    }
    else
    {
        acb_set_round(tmp1, x, wprec);
        acb_set_round(tmp2, y, wprec);
        acb_mul(z, tmp1, tmp2, wprec);
    }
}

void
_acb_modular_eta_sum_basecase(acb_t eta, const acb_t q, double log2q_approx, slong N, slong prec)
{
    slong e, e1, e2, k, k1, k2, num, term_prec;
    slong *exponents, *aindex, *bindex;
    acb_ptr qpow;
    acb_t tmp1, tmp2;
    double log2term_approx;

    if (N <= 5)
    {
        if (N <= 1)
        {
            acb_set_ui(eta, N != 0);
        }
        else if (N == 2)
        {
            acb_sub_ui(eta, q, 1, prec);
            acb_neg(eta, eta);
        }
        else
        {
            acb_mul(eta, q, q, prec);
            acb_add(eta, eta, q, prec);
            acb_neg(eta, eta);
            acb_add_ui(eta, eta, 1, prec);
        }
        return;
    }

    num = 1;
    while (PENTAGONAL(num) < N)
        num++;

    acb_init(tmp1);
    acb_init(tmp2);

    exponents = flint_malloc(sizeof(slong) * 3 * num);
    aindex = exponents + num;
    bindex = aindex + num;

    qpow = _acb_vec_init(num);

    acb_modular_addseq_eta(exponents, aindex, bindex, num);
    acb_set_round(qpow + 0, q, prec);

    acb_zero(eta);

    for (k = 0; k < num; k++)
    {
        e = exponents[k];

        log2term_approx = e * log2q_approx;
        term_prec = FLINT_MIN(FLINT_MAX(prec + log2term_approx + 16.0, 16.0), prec);

        if (k > 0)
        {
            k1 = aindex[k];
            k2 = bindex[k];
            e1 = exponents[k1];
            e2 = exponents[k2];

            if (e == e1 + e2)
            {
                _acb_modular_mul(qpow + k, tmp1, tmp2, qpow + k1, qpow + k2, term_prec, prec);
            }
            else if (e == 2 * e1 + e2)
            {
                _acb_modular_mul(qpow + k, tmp1, tmp2, qpow + k1, qpow + k1, term_prec, prec);
                _acb_modular_mul(qpow + k, tmp1, tmp2, qpow + k, qpow + k2, term_prec, prec);
            }
            else
            {
                flint_throw(FLINT_ERROR, "exponent not in addition sequence!\n");
            }
        }

        if (k % 4 <= 1)
            acb_sub(eta, eta, qpow + k, prec);
        else
            acb_add(eta, eta, qpow + k, prec);
    }

    acb_add_ui(eta, eta, 1, prec);

    flint_free(exponents);
    _acb_vec_clear(qpow, num);
    acb_clear(tmp1);
    acb_clear(tmp2);
}

void
_acb_modular_eta_sum_rs(acb_t eta, const acb_t q, double log2q_approx, slong N, slong prec)
{
    slong * tab;
    slong k, term_prec, i, e, eprev;
    slong m, num_pentagonal;
    double log2term_approx;
    acb_ptr qpow;
    acb_t tmp1, tmp2;

    acb_init(tmp1);
    acb_init(tmp2);

    /* choose rectangular splitting parameters */
    m = acb_modular_rs_optimal_m(pentagonal_best_m, pentagonal_best_m_residues, N);

    /* build addition sequence */
    tab = flint_calloc(m + 1, sizeof(slong));

    for (k = 0; PENTAGONAL(k) < N; k++)
        tab[PENTAGONAL(k) % m] = -1;
    num_pentagonal = k;
    tab[m] = -1;

    /* compute powers in addition sequence */
    qpow = _acb_vec_init(m + 1);
    acb_modular_fill_addseq(tab, m + 1);

    for (k = 0; k < m + 1; k++)
    {
        if (k == 0)
        {
            acb_one(qpow + k);
        }
        else if (k == 1)
        {
            acb_set_round(qpow + k, q, prec);
        }
        else if (tab[k] != 0)
        {
            log2term_approx = k * log2q_approx;
            term_prec = FLINT_MIN(FLINT_MAX(prec + log2term_approx + 16.0, 16.0), prec);
            _acb_modular_mul(qpow + k, tmp1, tmp2, qpow + tab[k], qpow + k - tab[k], term_prec, prec);
        }
    }

    /* compute eta */
    acb_zero(eta);
    term_prec = prec;

    for (k = num_pentagonal - 1; k >= 0; k--)
    {
        e = PENTAGONAL(k);  /* exponent */
        eprev = PENTAGONAL(k+1);

        log2term_approx = e * log2q_approx;
        term_prec = FLINT_MIN(FLINT_MAX(prec + log2term_approx + 16.0, 16.0), prec);

        /* giant steps */
        for (i = e / m; i < eprev / m; i++)
        {
            if (!acb_is_zero(eta))
                _acb_modular_mul(eta, tmp1, tmp2, eta, qpow + m, term_prec, prec);
        }

        if (k % 4 <= 1)
            acb_sub(eta, eta, qpow + (e % m), prec);
        else
            acb_add(eta, eta, qpow + (e % m), prec);
    }

    acb_add_ui(eta, eta, 1, prec);

    acb_clear(tmp1);
    acb_clear(tmp2);

    _acb_vec_clear(qpow, m + 1);
    flint_free(tab);
}

void
acb_modular_eta_sum(acb_t eta, const acb_t q, slong prec)
{
    mag_t err, qmag;
    double log2q_approx;
    int q_is_real;
    slong N;

    mag_init(err);
    mag_init(qmag);

    q_is_real = arb_is_zero(acb_imagref(q));

    acb_get_mag(qmag, q);
    log2q_approx = mag_get_d_log2_approx(qmag);

    if (log2q_approx >= 0.0)
    {
        N = 1;
        mag_inf(err);
    }
    else  /* Pick N and compute error bound */
    {
        N = 0;
        while (0.05 * N * N < prec)
        {
            if (log2q_approx * PENTAGONAL(N) < -prec - 2)
                break;
            N++;
        }
        N = PENTAGONAL(N);

        mag_geom_series(err, qmag, N);
        if (mag_is_inf(err))
            N = 1;
    }

    if (N < 400)
        _acb_modular_eta_sum_basecase(eta, q, log2q_approx, N, prec);
    else
        _acb_modular_eta_sum_rs(eta, q, log2q_approx, N, prec);

    if (q_is_real)
        arb_add_error_mag(acb_realref(eta), err);
    else
        acb_add_error_mag(eta, err);

    mag_clear(err);
    mag_clear(qmag);
}

