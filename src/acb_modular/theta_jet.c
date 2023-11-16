/*
    Copyright (C) 2017 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpq.h"
#include "acb_modular.h"
#include "acb_poly.h"

static void
_acb_vec_mul_4th_root(acb_ptr y, acb_srcptr x, slong len, int r, slong prec)
{
    slong k;
    r &= 7;

    if (r == 0)
    {
        _acb_vec_set(y, x, len);
    }
    else if (r == 4)
    {
        _acb_vec_neg(y, x, len);
    }
    else if (r == 2)
    {
        for (k = 0; k < len; k++)
            acb_mul_onei(y + k, x + k);
    }
    else if (r == 6)
    {
        for (k = 0; k < len; k++)
        {
            acb_mul_onei(y + k, x + k);
            acb_neg(y + k, y + k);
        }
    }
    else
    {
        fmpq_t t;
        acb_t w;
        fmpq_init(t);
        acb_init(w);
        fmpq_set_si(t, r, 4);
        arb_sin_cos_pi_fmpq(acb_imagref(w), acb_realref(w), t, prec);
        _acb_vec_scalar_mul(y, x, len, w, prec);
        fmpq_clear(t);
        acb_clear(w);
    }
}

void
acb_modular_theta_jet(acb_ptr theta1, acb_ptr theta2,
    acb_ptr theta3, acb_ptr theta4, const acb_t z, const acb_t tau,
    slong len, slong prec)
{
    fmpq_t t;
    psl2z_t g;
    arf_t one_minus_eps;
    acb_t z_prime, tau_prime, q, q4, w, A;
    acb_ptr B;
    acb_ptr thetas[4];
    int w_is_unit, R[4], S[4], C, rescale;
    slong k;

    if (len == 0)
        return;

    if (len == 1)
    {
        acb_modular_theta(theta1, theta2, theta3, theta4, z, tau, prec);
        return;
    }

    psl2z_init(g);
    arf_init(one_minus_eps);
    acb_init(tau_prime);

    /* move tau to the fundamental domain */
    arf_set_ui_2exp_si(one_minus_eps, 63, -6);
    acb_modular_fundamental_domain_approx(tau_prime, g, tau,
        one_minus_eps, prec);

    if (psl2z_is_one(g) &&
            arf_cmpabs_2exp_si(arb_midref(acb_imagref(z)), 4) <= 0)
    {
        acb_modular_theta_jet_notransform(theta1, theta2, theta3, theta4,
            z, tau, len, prec);
    }
    else
    {
        fmpq_init(t);
        acb_init(z_prime);
        acb_init(q);
        acb_init(q4);
        acb_init(w);
        acb_init(A);
        B = _acb_vec_init(len);
        thetas[0] = _acb_vec_init(len);
        thetas[1] = _acb_vec_init(len);
        thetas[2] = _acb_vec_init(len);
        thetas[3] = _acb_vec_init(len);

        /* compute transformation parameters */
        acb_modular_theta_transform(R, S, &C, g);

        if (C == 0)
        {
            acb_set(z_prime, z);
            acb_one(A);
            rescale = 0;
        }
        else
        {
            rescale = 1;

            /* B = 1/(c*tau+d) (temporarily) */
            acb_mul_fmpz(B, tau, &g->c, prec);
            acb_add_fmpz(B, B, &g->d, prec);
            acb_inv(B, B, prec);

            /* z' = -z/(c*tau+d) */
            acb_mul(z_prime, z, B, prec);
            acb_neg(z_prime, z_prime);

            /* A = sqrt(i/(c*tau+d)) */
            acb_mul_onei(A, B);
            acb_sqrt(A, A, prec);

            /* B = c/(c*tau+d) */
            acb_mul_fmpz(B, B, &g->c, prec);

            /* B[2] = -c/(c*tau+d) */
            if (len >= 3)
                acb_neg(B + 2, B);

            if (len >= 2)
            {
                /* B[1] = -2*z*c/(c*tau+d) */
                acb_mul(B + 1, B, z, prec);
                acb_mul_2exp_si(B + 1, B + 1, 1);
                acb_neg(B + 1, B + 1);
            }

            acb_mul(B, z_prime, z, prec);
            acb_mul_fmpz(B, B, &g->c, prec);

            /* we will have B = exp(-pi i c (z+x)^2/(c*tau+d))
               after computing the exponential later */
        }

        /* reduce z_prime modulo tau_prime if the imaginary part is large */
        if (arf_cmpabs_2exp_si(arb_midref(acb_imagref(z_prime)), 4) > 0)
        {
            arb_t nn;
            arb_init(nn);
            arf_div(arb_midref(nn), arb_midref(acb_imagref(z_prime)),
                arb_midref(acb_imagref(tau_prime)), prec, ARF_RND_DOWN);
            arf_mul_2exp_si(arb_midref(nn), arb_midref(nn), 1);
            arf_add_ui(arb_midref(nn), arb_midref(nn), 1, prec, ARF_RND_DOWN);
            arf_mul_2exp_si(arb_midref(nn), arb_midref(nn), -1);
            arf_floor(arb_midref(nn), arb_midref(nn));

            /* transform z_prime further */
            acb_submul_arb(z_prime, tau_prime, nn, prec);

            /* add -tau n^2 - 2n(z+x)' to B */
            arb_mul_2exp_si(nn, nn, 1);
            acb_submul_arb(B, z_prime, nn, prec);
            if (len >= 2)
            {
                acb_t u;
                acb_init(u);

                /* the x picks up a factor -1/(tau*c+d) */
                if (rescale)
                {
                    acb_mul_fmpz(u, tau, &g->c, prec);
                    acb_add_fmpz(u, u, &g->d, prec);
                    acb_inv(u, u, prec);
                    acb_neg(u, u);
                    acb_mul_arb(u, u, nn, prec);
                    acb_sub(B + 1, B + 1, u, prec);
                }
                else
                {
                    acb_sub_arb(B + 1, B + 1, nn, prec);
                }

                acb_clear(u);
            }
            arb_mul_2exp_si(nn, nn, -1);
            arb_sqr(nn, nn, prec);
            acb_submul_arb(B, tau_prime, nn, prec);

            /* theta1, theta4 pick up factors (-1)^n */
            if (!arf_is_int_2exp_si(arb_midref(nn), 1))
            {
                int i;
                for (i = 0; i < 4; i++)
                {
                    if (S[i] == 0 || S[i] == 3)
                        R[i] += 4;
                }
            }

            C = 1;

            arb_clear(nn);
        }

        if (C != 0)
            _acb_poly_exp_pi_i_series(B, B, FLINT_MIN(len, 3), len, prec);

        /* compute q_{1/4}, q */
        acb_mul_2exp_si(q4, tau_prime, -2);
        acb_exp_pi_i(q4, q4, prec);
        acb_pow_ui(q, q4, 4, prec);

        /* compute w */
        acb_exp_pi_i(w, z_prime, prec);
        w_is_unit = arb_is_zero(acb_imagref(z_prime));

        /* evaluate theta functions of transformed variables */
        acb_modular_theta_sum(thetas[0], thetas[1], thetas[2], thetas[3],
            w, w_is_unit, q, len, prec);

        /* correct for change of variables */
        if (rescale)
        {
            /* [-1/(tau*c+d)]]^k */
            acb_mul_fmpz(z_prime, tau, &g->c, prec);
            acb_add_fmpz(z_prime, z_prime, &g->d, prec);
            acb_inv(z_prime, z_prime, prec);
            acb_neg(z_prime, z_prime);

            acb_set(w, z_prime);

            for (k = 1; k < len; k++)
            {
                acb_mul(thetas[0] + k, thetas[0] + k, w, prec);
                acb_mul(thetas[1] + k, thetas[1] + k, w, prec);
                acb_mul(thetas[2] + k, thetas[2] + k, w, prec);
                acb_mul(thetas[3] + k, thetas[3] + k, w, prec);
                acb_mul(w, w, z_prime, prec);
            }
        }

        /* todo: fuse */
        _acb_vec_scalar_mul(thetas[0], thetas[0], len, q4, prec);
        _acb_vec_scalar_mul(thetas[1], thetas[1], len, q4, prec);

        /* multiply by roots of unity */
        _acb_vec_mul_4th_root(theta1, thetas[S[0]], len, R[0], prec);
        _acb_vec_mul_4th_root(theta2, thetas[S[1]], len, R[1], prec);
        _acb_vec_mul_4th_root(theta3, thetas[S[2]], len, R[2], prec);
        _acb_vec_mul_4th_root(theta4, thetas[S[3]], len, R[3], prec);

        if (C != 0)
        {
            _acb_vec_scalar_mul(B, B, len, A, prec);

            _acb_poly_mullow(thetas[0], theta1, len, B, len, len, prec);
            _acb_poly_mullow(thetas[1], theta2, len, B, len, len, prec);
            _acb_poly_mullow(thetas[2], theta3, len, B, len, len, prec);
            _acb_poly_mullow(thetas[3], theta4, len, B, len, len, prec);

            for (k = 0; k < len; k++) acb_swap(theta1 + k, thetas[0] + k);
            for (k = 0; k < len; k++) acb_swap(theta2 + k, thetas[1] + k);
            for (k = 0; k < len; k++) acb_swap(theta3 + k, thetas[2] + k);
            for (k = 0; k < len; k++) acb_swap(theta4 + k, thetas[3] + k);
        }

        fmpq_clear(t);
        acb_clear(z_prime);
        acb_clear(q);
        acb_clear(q4);
        acb_clear(w);
        acb_clear(A);
        _acb_vec_clear(B, len);
        _acb_vec_clear(thetas[0], len);
        _acb_vec_clear(thetas[1], len);
        _acb_vec_clear(thetas[2], len);
        _acb_vec_clear(thetas[3], len);
    }

    psl2z_clear(g);
    arf_clear(one_minus_eps);
    acb_clear(tau_prime);
}

