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

void
acb_modular_theta_const_sum_basecase(acb_t theta2, acb_t theta3, acb_t theta4,
    const acb_t q, slong N, slong prec)
{
    slong * tab;
    slong k, term_prec;
    double log2q_approx, log2term_approx;
    mag_t qmag;
    acb_ptr qpow;
    acb_t s1, s2, s3, t1, t2;

    if (N < 2)
    {
        acb_set_ui(theta2, 2 * (N > 0));
        acb_set_ui(theta3, N > 0);
        acb_set(theta4, theta3);
        return;
    }

    if (N < 25)
    {
        acb_t q1, q2, q4, q8, q16;

        acb_init(q1);
        acb_init(q2);
        acb_init(q4);
        acb_init(q8);
        acb_init(q16);

        acb_set_round(q1, q, prec);

        if (N > 2) acb_mul(q2, q1, q1, prec);
        if (N > 4) acb_mul(q4, q2, q2, prec);
        if (N > 9) acb_mul(q8, q4, q4, prec);
        if (N > 16) acb_mul(q16, q8, q8, prec);

        /* theta2 = 2 + 2q^2 + 2q^4 [2q^2 + 2q^8 + 2q^16] */
        if (N > 6)
        {
            if (N > 12)
            {
                acb_add(theta2, q2, q8, prec);
                if (N > 20)
                    acb_add(theta2, theta2, q16, prec);
                acb_mul(theta2, theta2, q4, prec);
            }
            else
            {
                acb_mul(theta2, q2, q4, prec);
            }
            acb_add(theta2, theta2, q2, prec);
            acb_add_ui(theta2, theta2, 1, prec);
        }
        else if (N > 2)
            acb_add_ui(theta2, q2, 1, prec);
        else
            acb_one(theta2);

        acb_mul_2exp_si(theta2, theta2, 1);

        /* theta3 = [1 + 2q^4 + 2q^16] + [2q + 2q^9] */
        /* theta4 = [1 + 2q^4 + 2q^16] - [2q + 2q^9] */
        if (N > 4)
        {
            if (N > 16)
                acb_add(q4, q4, q16, prec);
            acb_mul_2exp_si(q4, q4, 1);
            acb_add_ui(q4, q4, 1, prec);

            if (N > 9)
                acb_addmul(q1, q1, q8, prec);
            acb_mul_2exp_si(q1, q1, 1);

            acb_add(theta3, q4, q1, prec);
            acb_sub(theta4, q4, q1, prec);
        }
        else
        {
            acb_mul_2exp_si(q1, q1, 1);
            acb_add_ui(theta3, q1, 1, prec);
            acb_sub_ui(theta4, q1, 1, prec);
            acb_neg(theta4, theta4);
        }

        acb_clear(q1);
        acb_clear(q2);
        acb_clear(q4);
        acb_clear(q8);
        acb_clear(q16);

        return;
    }

    mag_init(qmag);
    acb_init(s1);
    acb_init(s2);
    acb_init(s3);
    acb_init(t1);
    acb_init(t2);

    tab = flint_calloc(N, sizeof(slong));
    qpow = _acb_vec_init(N);

    for (k = 0; k*(k+1) < N; k++) tab[k*(k+1)] = -1;
    for (k = 0; 4*k*k < N; k++) tab[4*k*k] = -1;
    for (k = 0; 4*k*(k+1) + 1 < N; k++) tab[4*k*(k+1)] = -1;
    if (N > 0) tab[0] = 0;
    if (N > 1) tab[1] = 1;

    acb_modular_fill_addseq(tab, N);

    acb_get_mag(qmag, q);
    log2q_approx = mag_get_d_log2_approx(qmag);

    for (k = 0; k < N; k++)
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
            _acb_modular_mul(qpow + k, t1, t2, qpow + tab[k], qpow + k - tab[k], term_prec, prec);
        }
    }

    for (k = 0; k*(k+1) < N; k++) acb_add(s1, s1, qpow + k*(k+1), prec);
    for (k = 1; 4*k*k < N; k++) acb_add(s2, s2, qpow + 4*k*k, prec);
    for (k = 0; 4*k*(k+1) + 1 < N; k++) acb_add(s3, s3, qpow + 4*k*(k+1), prec);

    /*
    theta2 = 2 + 2q^2 + 2q^6 + 2q^12 + 2q^20 + 2q^30 + ...
    theta3 = 1 + 2 (q^4 + q^16 + ...) + 2q (1 + q^8 + q^24 + ...)
    theta4 = 1 + 2 (q^4 + q^16 + ...) - 2q (1 + q^8 + q^24 + ...)
    */
    acb_mul(s3, s3, q, prec);
    acb_mul_2exp_si(s3, s3, 1);
    acb_mul_2exp_si(s2, s2, 1);
    acb_add(theta3, s2, s3, prec);
    acb_sub(theta4, s2, s3, prec);
    acb_add_ui(theta3, theta3, 1, prec);
    acb_add_ui(theta4, theta4, 1, prec);
    acb_mul_2exp_si(theta2, s1, 1);

    _acb_vec_clear(qpow, N);
    flint_free(tab);

    acb_clear(s1);
    acb_clear(s2);
    acb_clear(s3);
    acb_clear(t1);
    acb_clear(t2);
    mag_clear(qmag);
}

