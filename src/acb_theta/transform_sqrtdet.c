/*
    Copyright (C) 2023 Jean Kieffer

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb_poly.h"
#include "acb_mat.h"
#include "acb_theta.h"

static void
acb_theta_sqrt_branch(acb_t res, const acb_t x, acb_srcptr rts_neg, slong nb_neg,
    acb_srcptr rts_pos, slong nb_pos, const acb_t sqrt_lead, slong prec)
{
    acb_t s, t;
    slong k;

    acb_init(s);
    acb_init(t);

    acb_set(s, sqrt_lead);
    for (k = 0; k < nb_neg; k++)
    {
        acb_sub(t, x, &rts_neg[k], prec);
        acb_sqrt_analytic(t, t, 1, prec);
        acb_mul(s, s, t, prec);
    }
    for (k = 0; k < nb_pos; k++)
    {
        acb_sub(t, &rts_pos[k], x, prec);
        acb_sqrt_analytic(t, t, 1, prec);
        acb_mul(s, s, t, prec);
    }

    acb_set(res, s);

    acb_clear(s);
    acb_clear(t);
}

void
acb_theta_transform_sqrtdet(acb_t res, const acb_mat_t tau, slong prec)
{
    slong g = acb_mat_nrows(tau);
    flint_rand_t state;
    acb_mat_t A, B, C;
    acb_poly_t pol, h;
    acb_ptr rts, rts_neg, rts_pos;
    acb_t z, rt, mu;
    arb_t x;
    slong k, j, nb_neg, nb_pos;
    int success = 0;

    flint_randinit(state);
    acb_mat_init(A, g, g);
    acb_mat_init(B, g, g);
    acb_mat_init(C, g, g);
    acb_poly_init(pol);
    acb_poly_init(h);
    rts = _acb_vec_init(g);
    rts_neg = _acb_vec_init(g);
    rts_pos = _acb_vec_init(g);
    acb_init(z);
    acb_init(rt);
    acb_init(mu);
    arb_init(x);

    /* Choose a purely imaginary matrix A and compute pol s.t. pol(-1) = det(A)
       and pol(1) = det(tau): pol(t) is

       det(A + (t+1)/2 (tau - A)) = det(A) det(I - (t+1)/2 (I - A^{-1}tau))

       We want to get the g roots of this polynomial to compute the branch
       cuts. This can fail e.g. when det(tau - A) = 0, so pick A at random
       until the roots can be found */
    for (k = 0; (k < 100) && !success; k++)
    {
        acb_mat_onei(A);
        for (j = 0; j < g; j++)
        {
            arb_urandom(x, state, prec);
            arb_add(acb_imagref(acb_mat_entry(A, j, j)),
                acb_imagref(acb_mat_entry(A, j, j)), x, prec);
        }
        acb_mat_inv(B, A, prec);
        acb_mat_mul(B, B, tau, prec);
        acb_mat_one(C);
        acb_mat_sub(C, C, B, prec);

        /* Get reverse of charpoly */
        acb_mat_charpoly(h, C, prec);
        acb_poly_zero(pol);
        for (j = 0; j <= g; j++)
        {
            acb_poly_get_coeff_acb(z, h, j);
            acb_poly_set_coeff_acb(pol, g - j, z);
        }
        acb_poly_one(h);
        acb_poly_set_coeff_si(h, 1, 1);
        acb_poly_scalar_mul_2exp_si(h, h, -1);
        acb_poly_compose(pol, pol, h, prec);

        success = (acb_poly_find_roots(rts, pol, NULL, 0, prec) == g);

        /* Check that no root intersects the [-1,1] segment */
        for (j = 0; (j < g) && success; j++)
        {
            if (arb_contains_zero(acb_imagref(&rts[j])))
            {
                arb_abs(x, acb_realref(&rts[j]));
                arb_sub_si(x, x, 1, prec);
                success = arb_is_positive(x);
            }
        }
    }

    if (success)
    {
        /* Partition the roots between positive & negative real parts to
           compute branch for sqrt(pol) */
        nb_neg = 0;
        nb_pos = 0;
        for (k = 0; k < g; k++)
        {
            if (arb_is_negative(acb_realref(&rts[k])))
            {
                acb_set(&rts_neg[nb_neg], &rts[k]);
                nb_neg++;
            }
            else
            {
                acb_set(&rts_pos[nb_pos], &rts[k]);
                nb_pos++;
            }
        }
        acb_mat_det(rt, A, prec);
        acb_mul(rt, rt, acb_poly_get_coeff_ptr(pol, g), prec);
        acb_sqrts(rt, z, rt, prec);

        /* Set mu to +-1 such that mu*sqrt_branch gives the correct value at A,
           i.e. i^(g/2) * something positive */
        acb_mat_det(mu, A, prec);
        acb_mul_i_pow_si(mu, mu, -g);
        acb_sqrt(mu, mu, prec);
        acb_set_si(z, g);
        acb_mul_2exp_si(z, z, -2);
        acb_exp_pi_i(z, z, prec);
        acb_mul(mu, mu, z, prec);
        acb_set_si(z, -1);
        acb_theta_sqrt_branch(z, z, rts_neg, nb_neg, rts_pos, nb_pos, rt, prec);
        acb_div(mu, mu, z, prec);

        /* Compute square root branch at z=1 to get sqrtdet */
        acb_set_si(z, 1);
        acb_theta_sqrt_branch(rt, z, rts_neg, nb_neg, rts_pos, nb_pos, rt, prec);
        acb_mul(rt, rt, mu, prec);
        acb_mat_det(res, tau, prec);
        acb_theta_agm_sqrt(res, res, rt, 1, prec);
    }
    else
    {
        acb_mat_det(res, tau, prec);
        acb_sqrts(res, z, res, prec);
        acb_union(res, res, z, prec);
    }

    flint_randclear(state);
    acb_mat_clear(A);
    acb_mat_clear(B);
    acb_mat_clear(C);
    acb_poly_clear(pol);
    acb_poly_clear(h);
    _acb_vec_clear(rts, g);
    _acb_vec_clear(rts_pos, g);
    _acb_vec_clear(rts_neg, g);
    acb_clear(z);
    acb_clear(rt);
    acb_clear(mu);
    arb_clear(x);
}
