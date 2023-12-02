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

void
acb_modular_theta_sum(acb_ptr theta1,
                          acb_ptr theta2,
                          acb_ptr theta3,
                          acb_ptr theta4,
    const acb_t w, int w_is_unit, const acb_t q, slong len, slong prec)
{
    mag_t qmag, wmag, vmag;
    mag_ptr err;
    double log2q_approx, log2w_approx, log2term_approx;
    slong e, e1, e2, k, k1, k2, r, n, N, WN, term_prec;
    slong *exponents, *aindex, *bindex;
    acb_ptr qpow, wpow, vpow;
    acb_t tmp1, tmp2, v;
    int q_is_real, w_is_one;

    q_is_real = arb_is_zero(acb_imagref(q));
    w_is_one = acb_is_one(w);

    if (w_is_one && len == 1)
    {
        acb_modular_theta_const_sum(theta2, theta3, theta4, q, prec);
        acb_zero(theta1);
        return;
    }

    mag_init(qmag);
    mag_init(wmag);
    mag_init(vmag);

    acb_init(tmp1);
    acb_init(tmp2);
    acb_init(v);
    err = _mag_vec_init(len);

    if (w_is_one)
        acb_one(v);
    else if (w_is_unit)
        acb_conj(v, w);
    else
        acb_inv(v, w, prec);

    acb_get_mag(qmag, q);
    log2q_approx = mag_get_d_log2_approx(qmag);

    if (w_is_unit)
    {
        mag_one(wmag);
        mag_one(vmag);
        log2w_approx = 0.0;
    }
    else
    {
        acb_get_mag(wmag, w);
        acb_get_mag(vmag, v);
        mag_max(wmag, wmag, vmag);
        log2w_approx = mag_get_d_log2_approx(wmag);
    }

    if (log2q_approx >= 0.0)
    {
        N = 1;
        for (r = 0; r < len; r++)
            mag_inf(err + r);
    }
    else  /* Pick N and compute error bound */
    {
        mag_t den, cmag, dmag;

        mag_init(den);
        mag_init(cmag);
        mag_init(dmag);

        N = 1;

        while (0.05 * N * N < prec)
        {
            log2term_approx = log2q_approx * ((N+2)*(N+2)/4) + (N+2)*log2w_approx;

            if (log2term_approx < -prec - 2)
                break;

            N++;
        }

        if (len == 1)
        {
            if (w_is_unit)
            {
                mag_one(den);
                mag_sub_lower(den, den, qmag);  /* 1 - |q| is good enough */
            }
            else  /* denominator: 1 - |q|^(floor((N+1)/2)+1) * max(|w|,1/|w|) */
            {
                mag_pow_ui(err, qmag, (N + 1) / 2 + 1);
                mag_mul(err, err, wmag);
                mag_one(den);
                mag_sub_lower(den, den, err);
            }

            /* no convergence */
            if (mag_is_zero(den))
            {
                N = 1;
                mag_inf(err);
            }
            else if (w_is_unit)
            {
                mag_pow_ui(err, qmag, ((N + 2) * (N + 2)) / 4);
                mag_div(err, err, den);
                mag_mul_2exp_si(err, err, 1);
            }
            else
            {
                mag_pow_ui(err, qmag, ((N + 2) * (N + 2)) / 4);
                mag_pow_ui(vmag, wmag, N + 2);
                mag_mul(err, err, vmag);
                mag_div(err, err, den);
                mag_mul_2exp_si(err, err, 1);
            }
        }
        else
        {
            /* numerator: 2 |q|^E * max(|w|,|v|)^(N+2) * (N+2)^r */
            mag_pow_ui(err, qmag, ((N + 2) * (N + 2)) / 4);

            if (!w_is_one)
            {
                mag_pow_ui(vmag, wmag, N + 2);
                mag_mul(err, err, vmag);
            }

            mag_mul_2exp_si(err, err, 1);

            for (r = 1; r < len; r++)
                mag_mul_ui(err + r, err + r - 1, N + 2);

            /* den: 1 - |q|^floor((N+1)/2+1) * max(|w|,|v|) * exp(r/(N+2)) */
            mag_pow_ui(cmag, qmag, (N + 1) / 2 + 1);
            mag_mul(cmag, cmag, wmag);

            for (r = 0; r < len; r++)
            {
                mag_set_ui(dmag, r);
                mag_div_ui(dmag, dmag, N + 2);
                mag_exp(dmag, dmag);
                mag_mul(dmag, cmag, dmag);
                mag_one(den);
                mag_sub_lower(den, den, dmag);

                if (mag_is_zero(den))
                    mag_inf(err + r);
                else
                    mag_div(err + r, err + r, den);
            }
        }

        /* don't do work if we can't determine the zeroth derivative */
        if (mag_is_inf(err))
            N = 1;

        mag_clear(den);
        mag_clear(cmag);
        mag_clear(dmag);
    }

    exponents = flint_malloc(sizeof(slong) * 3 * N);
    aindex = exponents + N;
    bindex = aindex + N;

    qpow = _acb_vec_init(N);

    acb_modular_addseq_theta(exponents, aindex, bindex, N);
    acb_set_round(qpow + 0, q, prec);

    _acb_vec_zero(theta1, len);
    _acb_vec_zero(theta2, len);
    _acb_vec_zero(theta3, len);
    _acb_vec_zero(theta4, len);

    WN = (N + 3) / 2;

    /* compute powers of w^2 and 1/w^2 */
    /* todo: conjugates... */
    if (!w_is_one)
    {
        wpow = _acb_vec_init(WN);
        vpow = _acb_vec_init(WN + 1);

        acb_mul(tmp1, w, w, prec);
        acb_mul(tmp2, v, v, prec);

        _acb_vec_set_powers(wpow, tmp1, WN, prec);
        _acb_vec_set_powers(vpow, tmp2, WN + 1, prec);
    }
    else
    {
        wpow = vpow = NULL;
    }

    for (k = 0; k < N; k++)
    {
        e = exponents[k];

        log2term_approx = e * log2q_approx + (k+2) * log2w_approx;
        term_prec = FLINT_MIN(FLINT_MAX(prec + log2term_approx + 16.0, 16.0), prec);

        if (k > 0)
        {
            k1 = aindex[k];
            k2 = bindex[k];
            e1 = exponents[k1];
            e2 = exponents[k2];

            if (e == e1 + e2)
            {
                _acb_modular_mul(qpow + k, tmp1, tmp2,
                    qpow + k1, qpow + k2, term_prec, prec);
            }
            else if (e == 2 * e1 + e2)
            {
                _acb_modular_mul(qpow + k, tmp1, tmp2,
                    qpow + k1, qpow + k1, term_prec, prec);
                _acb_modular_mul(qpow + k, tmp1, tmp2,
                    qpow + k, qpow + k2, term_prec, prec);
            }
            else
            {
                flint_throw(FLINT_ERROR, "exponent not in addition sequence!\n");
            }
        }

        if (w_is_one && len == 1)
        {
            if (k % 2 == 0)
            {
                acb_add(theta3, theta3, qpow + k, prec);

                if (k % 4 == 0)
                    acb_sub(theta4, theta4, qpow + k, prec);
                else
                    acb_add(theta4, theta4, qpow + k, prec);
            }
            else
            {
                acb_add(theta2, theta2, qpow + k, prec);
            }
        }
        else
        {
            n = k / 2 + 1;

            if (k % 2 == 0)
            {
                acb_ptr term;

                if (w_is_one)
                {
                    acb_mul_2exp_si(tmp1, qpow + k, 1);
                    acb_zero(tmp2);
                }
                else
                {
                    /* tmp1 = w^(2n) + v^(2n) ~= 2 cos(2n) */
                    acb_add(tmp1, wpow + n, vpow + n, term_prec);
                    acb_mul(tmp1, qpow + k, tmp1, term_prec);

                    /* tmp2 = w^(2n) - v^(2n) ~= 2 sin(2n) */
                    if (len > 1)
                    {
                        acb_sub(tmp2, wpow + n, vpow + n, term_prec);
                        acb_mul(tmp2, qpow + k, tmp2, term_prec);
                    }
                }

                /* compute all the derivatives */
                for (r = 0; r < len; r++)
                {
                    term = (r % 2 == 0) ? tmp1 : tmp2;

                    if (r == 1)
                        acb_mul_ui(term, term, 2 * n, term_prec);
                    else if (r > 1)
                        acb_mul_ui(term, term, 4 * n * n, term_prec);

                    acb_add(theta3 + r, theta3 + r, term, prec);

                    if (k % 4 == 0)
                        acb_sub(theta4 + r, theta4 + r, term, prec);
                    else
                        acb_add(theta4 + r, theta4 + r, term, prec);
                }
            }
            else
            {
                acb_ptr term;

                if (w_is_one)
                {
                    acb_mul_2exp_si(tmp1, qpow + k, 1);
                    acb_zero(tmp2);
                }
                else
                {
                    /* tmp1 = w^(2n) + v^(2n+2) ~= 2 cos(2n+1) / w */
                    acb_add(tmp1, wpow + n, vpow + n + 1, term_prec);
                    acb_mul(tmp1, qpow + k, tmp1, term_prec);

                    /* tmp2 = w^(2n) - v^(2n+2) ~= 2 sin(2n+1) / w */
                    acb_sub(tmp2, wpow + n, vpow + n + 1, term_prec);
                    acb_mul(tmp2, qpow + k, tmp2, term_prec);
                }

                /* compute all the derivatives */
                for (r = 0; r < len; r++)
                {
                    if (r > 0)
                    {
                        acb_mul_ui(tmp1, tmp1, 2 * n + 1, term_prec);
                        acb_mul_ui(tmp2, tmp2, 2 * n + 1, term_prec);
                    }

                    term = (r % 2 == 0) ? tmp2 : tmp1;

                    if (k % 4 == 1)
                        acb_sub(theta1 + r, theta1 + r, term, prec);
                    else
                        acb_add(theta1 + r, theta1 + r, term, prec);

                    term = (r % 2 == 0) ? tmp1 : tmp2;

                    acb_add(theta2 + r, theta2 + r, term, prec);
                }
            }
        }
    }

    if (w_is_one && len == 1)
    {
        acb_mul_2exp_si(theta2, theta2, 1);
        acb_mul_2exp_si(theta3, theta3, 1);
        acb_mul_2exp_si(theta4, theta4, 1);
    }

    /* theta1: w * sum + 2 sin */
    /* theta2: w * sum + 2 cos */

    if (!w_is_one)
    {
        _acb_vec_scalar_mul(theta1, theta1, len, w, prec);
        _acb_vec_scalar_mul(theta2, theta2, len, w, prec);

        acb_add(tmp1, w, v, prec);
        acb_sub(tmp2, w, v, prec);
    }
    else
    {
        acb_set_ui(tmp1, 2);
        acb_zero(tmp2);
    }

    for (r = 0; r < len; r++)
    {
        acb_add(theta1 + r, theta1 + r, (r % 2 == 0) ? tmp2 : tmp1, prec);
        acb_add(theta2 + r, theta2 + r, (r % 2 == 0) ? tmp1 : tmp2, prec);
    }

    /*
    Coefficient r in the z-expansion gains a factor: pi^r / r!
    times a sign:

        + 2 cos   =  +1 * (exp + 1/exp)
        - 2 sin   =  +i * (exp - 1/exp)
        - 2 cos   =  -1 * (exp + 1/exp)
        + 2 sin   =  -i * (exp - 1/exp)
        ...
    */

    acb_mul_onei(theta1, theta1);
    acb_neg(theta1, theta1);

    for (r = 1; r < len; r++)
    {
        if (r % 4 == 0)
        {
            acb_mul_onei(theta1 + r, theta1 + r);
            acb_neg(theta1 + r, theta1 + r);
        }
        else if (r % 4 == 1)
        {
            acb_mul_onei(theta2 + r, theta2 + r);
            acb_mul_onei(theta3 + r, theta3 + r);
            acb_mul_onei(theta4 + r, theta4 + r);
        }
        else if (r % 4 == 2)
        {
            acb_mul_onei(theta1 + r, theta1 + r);

            acb_neg(theta2 + r, theta2 + r);
            acb_neg(theta3 + r, theta3 + r);
            acb_neg(theta4 + r, theta4 + r);
        }
        else
        {
            acb_neg(theta1 + r, theta1 + r);

            acb_mul_onei(theta2 + r, theta2 + r);
            acb_mul_onei(theta3 + r, theta3 + r);
            acb_mul_onei(theta4 + r, theta4 + r);

            acb_neg(theta2 + r, theta2 + r);
            acb_neg(theta3 + r, theta3 + r);
            acb_neg(theta4 + r, theta4 + r);
        }
    }

    /* Add error bound. Note that this must be done after the
       rearrangements above, and before scaling by pi^r / r! below. */
    for (r = 0; r < len; r++)
    {
        if (q_is_real && w_is_unit)    /* result must be real */
        {
            arb_add_error_mag(acb_realref(theta1 + r), err + r);
            arb_add_error_mag(acb_realref(theta2 + r), err + r);
            arb_add_error_mag(acb_realref(theta3 + r), err + r);
            arb_add_error_mag(acb_realref(theta4 + r), err + r);
        }
        else
        {
            acb_add_error_mag(theta1 + r, err + r);
            acb_add_error_mag(theta2 + r, err + r);
            acb_add_error_mag(theta3 + r, err + r);
            acb_add_error_mag(theta4 + r, err + r);
        }
    }

    if (len > 1)
    {
        arb_t c, d;

        arb_init(c);
        arb_init(d);

        arb_const_pi(c, prec);
        arb_set(d, c);

        for (r = 1; r < len; r++)
        {
            acb_mul_arb(theta1 + r, theta1 + r, d, prec);
            acb_mul_arb(theta2 + r, theta2 + r, d, prec);
            acb_mul_arb(theta3 + r, theta3 + r, d, prec);
            acb_mul_arb(theta4 + r, theta4 + r, d, prec);

            if (r + 1 < len)
            {
                arb_mul(d, d, c, prec);
                arb_div_ui(d, d, r + 1, prec);
            }
        }

        arb_clear(c);
        arb_clear(d);
    }

    acb_add_ui(theta3, theta3, 1, prec);
    acb_add_ui(theta4, theta4, 1, prec);

    if (!w_is_one)
    {
        _acb_vec_clear(wpow, WN);
        _acb_vec_clear(vpow, WN + 1);
    }

    flint_free(exponents);
    _acb_vec_clear(qpow, N);
    acb_clear(tmp1);
    acb_clear(tmp2);
    acb_clear(v);
    mag_clear(qmag);
    mag_clear(wmag);
    mag_clear(vmag);
    _mag_vec_clear(err, len);
}

