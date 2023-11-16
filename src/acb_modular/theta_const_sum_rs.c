/*
    Copyright (C) 2015 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb.h"
#include "acb_modular.h"

static const int square_best_m[] = {
  2, 3, 4, 8, 12, 16, 32, 48, 80, 96,
  112, 144, 240, 288, 336, 480, 560, 576, 720, 1008,
  1440, 1680, 2016, 2640, 2880, 3600, 4032, 5040, 7920, 9360,
  10080, 15840, 18480, 20160, 25200, 31680, 37440, 39600, 44352, 50400,
  55440, 65520, 85680, 95760, 102960, 110880, 131040, 171360, 191520, 205920,
  221760, 262080, 277200, 327600, 383040, 411840, 514800, 554400, 655200, 720720,
  942480, 0
};

static const int square_best_m_residues[] = {
  2, 2, 2, 3, 4, 4, 7, 8, 12, 14,
  16, 16, 24, 28, 32, 42, 48, 48, 48, 64,
  84, 96, 112, 144, 144, 176, 192, 192, 288, 336,
  336, 504, 576, 576, 704, 864, 1008, 1056, 1152, 1232,
  1152, 1344, 1728, 1920, 2016, 2016, 2352, 3024, 3360, 3528,
  3456, 4032, 4224, 4928, 5760, 6048, 7392, 7392, 8624, 8064,
  10368, 0
};

static const int trigonal_best_m[] = {
  2, 6, 10, 14, 18, 30, 42, 66, 70, 90,
  126, 198, 210, 330, 390, 450, 630, 990, 1170, 1386,
  1638, 2142, 2310, 2730, 3150, 4950, 5850, 6930, 8190, 10710,
  11970, 12870, 16830, 18018, 23562, 26334, 27846, 30030, 34650, 40950,
  53550, 59850, 64350, 84150, 90090, 117810, 131670, 139230, 155610, 188370,
  203490, 218790, 244530, 270270, 296010, 306306, 342342, 414414, 447678, 450450,
  589050, 658350, 696150, 778050, 941850, 0
};

static const int trigonal_best_m_residues[] = {
  1, 2, 3, 4, 4, 6, 8, 12, 12, 12,
  16, 24, 24, 36, 42, 44, 48, 72, 84, 96,
  112, 144, 144, 168, 176, 264, 308, 288, 336, 432,
  480, 504, 648, 672, 864, 960, 1008, 1008, 1056, 1232,
  1584, 1760, 1848, 2376, 2016, 2592, 2880, 3024, 3360, 4032,
  4320, 4536, 5040, 5544, 6048, 6048, 6720, 8064, 8640, 7392,
  9504, 10560, 11088, 12320, 14784, 0,
};

slong
acb_modular_rs_optimal_m(const int * best_ms, const int * num_residues, slong N)
{
    slong i, m, cost, best_i, best_m, best_cost;

    best_i = 0;
    best_m = best_ms[0];
    best_cost = WORD_MAX;

    for (i = 0; (m = best_ms[i]) != 0; i++)
    {
        cost = N / m + num_residues[i];

        if (i == 0 || cost < best_cost)
        {
            best_i = i;
            best_cost = cost;
            best_m = m;
        }
    }

    /* flint_printf("N = %wd, best_m = %wd, best_cost = %wd, s(m) = %d\n",
        N, best_m, best_cost, num_residues[best_i]); */
    i = best_i;

    return best_m;
}

void
acb_modular_theta_const_sum_rs(acb_t theta2, acb_t theta3, acb_t theta4,
    const acb_t q, slong N, slong prec)
{
    slong * tab;
    slong k, term_prec, i, e, eprev;
    slong M, m2, m3, num_square, num_trigonal;
    double log2q_approx, log2term_approx;
    acb_ptr qpow;
    acb_t tmp1, tmp2;
    mag_t qmag;

    mag_init(qmag);
    acb_get_mag(qmag, q);
    log2q_approx = mag_get_d_log2_approx(qmag);
    mag_clear(qmag);

    acb_init(tmp1);
    acb_init(tmp2);

    /* choose rectangular splitting parameters */
    m2 = acb_modular_rs_optimal_m(trigonal_best_m, trigonal_best_m_residues, N);
    m3 = acb_modular_rs_optimal_m(square_best_m, square_best_m_residues, N);
    M = FLINT_MAX(m2, m3) + 1;

    /* build addition sequence */
    tab = flint_calloc(M, sizeof(slong));

    for (k = 0; k*(k+1) < N; k++)
        tab[(k*(k+1)) % m2] = -1;
    num_trigonal = k;

    for (k = 0; k*k < N; k++)
        tab[(k*k) % m3] = -1;
    num_square = k;

    tab[m2] = -1;
    tab[m3] = -1;

    /* compute powers in addition sequence */
    qpow = _acb_vec_init(M);
    acb_modular_fill_addseq(tab, M);

    for (k = 0; k < M; k++)
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

    /* compute theta2 */
    acb_zero(theta2);
    term_prec = prec;

    for (k = num_trigonal - 1; k >= 0; k--)
    {
        e = k * (k + 1);  /* exponent */
        eprev = (k + 1) * (k + 2);

        log2term_approx = e * log2q_approx;
        term_prec = FLINT_MIN(FLINT_MAX(prec + log2term_approx + 16.0, 16.0), prec);

        /* giant steps */
        for (i = e / m2; i < eprev / m2; i++)
        {
            if (!acb_is_zero(theta2))
                _acb_modular_mul(theta2, tmp1, tmp2, theta2, qpow + m2, term_prec, prec);
        }

        acb_add(theta2, theta2, qpow + (e % m2), prec);
    }

    acb_mul_2exp_si(theta2, theta2, 1);

    /* compute theta3, theta4 */
    acb_zero(theta3);
    acb_zero(theta4);
    term_prec = prec;

    for (k = num_square - 1; k >= 0; k--)
    {
        e = k * k;  /* exponent */
        eprev = (k + 1) * (k + 1);

        log2term_approx = e * log2q_approx;
        term_prec = FLINT_MIN(FLINT_MAX(prec + log2term_approx + 16.0, 16.0), prec);

        /* giant steps */
        for (i = e / m3; i < eprev / m3; i++)
        {
            if (!acb_is_zero(theta3))
                _acb_modular_mul(theta3, tmp1, tmp2, theta3, qpow + m3, term_prec, prec);

            if (!acb_is_zero(theta4))
                _acb_modular_mul(theta4, tmp1, tmp2, theta4, qpow + m3, term_prec, prec);
        }

        if (k == 0)
        {
            acb_mul_2exp_si(theta3, theta3, 1);
            acb_mul_2exp_si(theta4, theta4, 1);
        }

        acb_add(theta3, theta3, qpow + (e % m3), prec);

        if (k % 2 == 0)
            acb_add(theta4, theta4, qpow + (e % m3), prec);
        else
            acb_sub(theta4, theta4, qpow + (e % m3), prec);
    }

    acb_clear(tmp1);
    acb_clear(tmp2);

    _acb_vec_clear(qpow, M);
    flint_free(tab);
}

