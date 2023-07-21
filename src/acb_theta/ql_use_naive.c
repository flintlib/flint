/*
    Copyright (C) 2023 Jean Kieffer

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "acb_theta.h"

int
acb_theta_ql_use_naive(acb_ptr r, acb_srcptr t, acb_srcptr z, arb_srcptr dist,
    const acb_mat_t tau, slong d, slong prec, acb_theta_ql_worker_t worker_d)
{
    slong g = acb_mat_nrows(tau);
    slong n = 1 << g;
    slong nb_a = 1 << (g - d);
    slong nb_th = 1 << d;
    slong nb_t = (_acb_vec_is_zero(t, g) ? 1 : 3);
    arb_mat_t Yinv, cho, cho1;
    acb_mat_t tau0, star, tau1;
    arb_ptr offset, z_offset, new_dist;
    acb_ptr v, w, new_z, new_th;
    arf_t eps, R2;
    arb_t max_dist, x;
    acb_t c, f;
    acb_theta_eld_t E;
    slong* pts;
    slong newprec, fullprec;
    ulong a;
    slong j, k, l;
    int res = 1;

    if (d == 0)
    {
        res = worker_d(r, t, z, dist, tau, prec);
        return res;
    }

    arb_mat_init(Yinv, g, g);
    arb_mat_init(cho, g, g);
    arb_mat_init(cho1, g - d, g - d);
    acb_mat_init(tau0, d, d);
    acb_mat_init(star, d, g - d);
    acb_mat_init(tau1, g - d, g - d);
    offset = _arb_vec_init(g - d);
    z_offset = _arb_vec_init(g);
    new_dist = _arb_vec_init(nb_th);
    v = _acb_vec_init(g - d);
    w = _acb_vec_init(g - d);
    new_z = _acb_vec_init(d);
    new_th = _acb_vec_init(nb_th * nb_t);
    arf_init(R2);
    arb_init(max_dist);
    arb_init(x);
    acb_init(c);
    acb_init(f);

    acb_theta_ql_blocks(tau0, star, tau1, tau, d);
    acb_theta_eld_cho(cho, tau, prec);
    acb_theta_eld_cho(cho1, tau1, prec);
    
    acb_mat_get_imag(Yinv, tau);
    arb_mat_inv(Yinv, Yinv, prec);
    _acb_vec_get_imag(z_offset, z, g);
    arb_mat_vector_mul_col(z_offset, Yinv, z_offset, prec);

    _acb_vec_zero(r, n * nb_t);
    for (a = 0; a < nb_a; a++)
    {
        /* Get R2 */
        arb_zero(max_dist);
        for (k = a; k < n; k += nb_a)
        {
            arb_max(max_dist, max_dist, &dist[k], prec);
        }
        fullprec = prec + acb_theta_ql_addprec(max_dist);
        arf_one(eps);
        arf_mul_2exp_si(eps, eps, -fullprec);
        acb_theta_naive_radius(R2, cho, 0, eps, prec);

        /* Get offset */
        acb_theta_char_get_arb(offset, a, g - d);
        _arb_vec_add(offset, offset, z_offset + d, g - d, prec);
        
        /* Make ellipsoid and list points */
        acb_theta_eld_init(E, g - d, g - d);
        acb_theta_eld_fill(E, cho1, R2, offset, prec);
        pts = flint_malloc(acb_theta_eld_nb_pts(E) * (g - d) * sizeof(slong));
        acb_theta_eld_points(pts, E);        

        /* Compute th_rec at each point using worker and sum */
        for (k = 0; (k < acb_theta_eld_nb_pts(E)) && res; k++)
        {
            /* Set v to pt + a1/2 */
            acb_theta_char_get_acb(v, a, g - d);
            for (j = 0; j < g - d; j++)
            {
                acb_add_si(&v[j], &v[j], pts[k * (g - d) + j], prec);
            }
            
            /* Get new_z and cofactor at 0 */
            acb_mat_vector_mul_col(new_z, star, v, prec);
            _acb_vec_add(new_z, new_z, z, d, prec);
            acb_dot(f, NULL, 0, v, 1, z + d, 1, g - d, prec);
            acb_mul_2exp_si(f, f, 1);
            acb_mat_vector_mul_col(w, tau1, v, prec);
            acb_dot(f, f, 0, w, 1, v, 1, d, prec);

            /* Get new distances and relative precision */
            acb_theta_ql_sqr_dists_a(new_dist, new_z, tau0, prec);
            acb_theta_ql_sqr_dist_pt(max_dist, offset, cho1, pts + k * (g - d), prec);
            newprec = prec;
            for (j = 0; j < nb_th; j++)
            {
                arb_sub(x, &dist[a + j * nb_a], max_dist, prec);
                arb_sub(x, x, &new_dist[j], prec);
                newprec = FLINT_MIN(newprec, acb_theta_ql_addprec(x)); /* <= prec */
                newprec = FLINT_MAX(newprec, ACB_THETA_ELD_DEFAULT_PREC);
            }            
            
            /* Call worker */
            res = worker_d(new_th, t, new_z, new_dist, tau0, newprec);

            /* Rescale to set r; cofactor depends on t */
            for (l = 0; l < nb_t; l++)
            {
                acb_dot(c, NULL, 0, v, 1, t + d, 1, g - d, prec);
                acb_mul_2exp_si(c, c, 1);
                acb_add(c, c, f, prec);
                acb_exp_pi_i(c, c, prec);
                _acb_vec_scalar_mul(new_th + l * nb_th, new_th + l * nb_th,
                    nb_th, c, prec);
                for (j = 0; j < nb_th; j++)
                {
                    acb_add(&r[l * n + j * nb_a + a], &r[l * n + j * nb_a + a],
                        &new_th[j], fullprec);
                }
            }
        }

        /* Add error */
        for (j = 0; j < nb_t * nb_th; j++)
        {
            for (l = 0; l < nb_t; l++)
            {
                acb_add_error_arf(&r[l * n + j * nb_a + a], eps);
            }
        }

        acb_theta_eld_clear(E);
        flint_free(pts);
    }
    
    arb_mat_clear(Yinv);
    arb_mat_clear(cho);
    arb_mat_clear(cho1);
    acb_mat_clear(tau0);
    acb_mat_clear(star);
    acb_mat_clear(tau1);
    _arb_vec_clear(offset, g - d);
    _arb_vec_clear(z_offset, g);
    _arb_vec_clear(new_dist, nb_th);
    _acb_vec_clear(v, g - d);
    _acb_vec_clear(w, g - d);
    _acb_vec_clear(new_z, d);
    _acb_vec_clear(new_th, nb_th * nb_t);
    arf_clear(R2);
    arb_clear(max_dist);
    arb_clear(x);
    acb_clear(c);
    acb_clear(f);
    return res;
}
