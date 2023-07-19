/*
    Copyright (C) 2023 Jean Kieffer

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "acb_theta.h"

static void
acb_theta_ql_tree_init_rec(acb_theta_ql_tree_t T, slong* cuts, slong nb_cuts,
    slong g, acb_srcptr z, const arb_mat_t cho, const acb_mat_t tau, slong prec)
{
    slong d = (nb_cuts == 0 ? 0 : cuts[nb_cuts - 1]);
    slong n = 1 << (g - d);
    arb_mat_t Y;
    arb_mat_t Yinv;
    acb_mat_t star;
    acb_ptr new_z, col;
    arf_t eps, R2;
    arb_ptr offset;
    slong* points;
    slong tot;
    ulong a;
    slong j, k;
    
    acb_theta_ql_tree_ambient_dim(T) = arb_mat_nrows(cho);
    acb_theta_ql_tree_dim(T) = d;
    acb_theta_ql_tree_nb_steps(T) = acb_theta_ql_new_nb_steps(cho, d, prec);
    acb_theta_ql_tree_z(T) = _acb_vec_init(g - d);
    _acb_vec_set(acb_theta_ql_tree_z(T), z, g - d);
    
    if (nb_cuts == 0)
    {
        T->eld = NULL;
        T->nb_children = NULL;
        T->index_children = NULL;
        T->children = NULL;
        acb_theta_ql_tree_total_children(T) = 0;
        return;
    }
    
    arb_mat_init(Y, g - d, g - d);
    arb_mat_init(Yinv, g - d, g - d);
    acb_mat_init(star, g - d, d);
    new_z = _acb_vec_init(g - d);
    col = _acb_vec_init(d);
    arf_init(eps);
    arf_init(R2);
    offset = _arb_vec_init(d);
    T->eld = flint_malloc(n * sizeof(struct acb_theta_ql_tree_struct));
    T->nb_children = flint_malloc((n + 1) * sizeof(slong));
    T->index_children = flint_malloc(n * sizeof(slong));

    /* Set matrices */
    for (j = d; j < g; j++)
    {
        for (k = d; k < g; k++)
        {
            arb_set(arb_mat_entry(Y, j - d, k - d), arb_mat_entry(cho, j, k));
        }
    }
    arb_mat_inv(Yinv, Y, prec);
    arb_mat_transpose(Yinv, Yinv);
    for (j = 0; j < d; j++)
    {
        for (k = d; k < g; k++)
        {
            acb_set(acb_mat_entry(star, j, k), acb_mat_entry(tau, j, k));
        }
    }

    /* Gather data for ellipsoids */
    _acb_vec_get_imag(offset, z, g - d);
    arb_mat_vector_mul_col(offset, Yinv, offset, prec);    
    arb_mat_scalar_mul_2exp_si(Y, Y, acb_theta_ql_tree_nb_steps(T));
    arf_one(eps);
    arf_mul_2exp_si(eps, eps, -prec); /* This is wrong; need the distance */
    acb_theta_naive_radius(R2, Y, 0, eps, ACB_THETA_ELD_DEFAULT_PREC);
    
    /* Make ellipsoids */
    for (a = 0; a < n; a++)
    {
        acb_theta_eld_init(acb_theta_ql_tree_eld(T, a), g - d, g - d);
        acb_theta_eld_fill(acb_theta_ql_tree_eld(T, a), Y, R2, offset, NULL, a,
            ACB_THETA_ELD_DEFAULT_PREC);
        acb_theta_ql_tree_nb_children(T, a) = acb_theta_eld_nb_pts(acb_theta_ql_tree_eld(T, a));
    }

    /* Count children, list points */
    acb_theta_ql_tree_index_children(T, 0) = 0;
    for (a = 0; a < n + 1; a++)
    {
        acb_theta_ql_tree_index_children(T, a) = acb_theta_ql_tree_index_children(T, a - 1)
            + acb_theta_ql_tree_nb_children(T, a - 1);
    }
    tot = acb_theta_ql_tree_index_children(T, n);
    acb_theta_ql_tree_total_children(T) = tot;
    points = flint_malloc(tot * (g - d) * sizeof(slong));
    for (a = 0; a < n; a++)
    {
        acb_theta_eld_points(points + acb_theta_ql_tree_index_children(T, a) * (g - d),
            acb_theta_ql_tree_eld(T, a));
    }
    
    /* Recursive calls */
    T->children = flint_malloc(tot * sizeof(struct acb_theta_ql_tree_struct));
    for (k = 0; k < tot; k++)
    {
        for (j = 0; j < g - d; j++)
        {
            acb_set_si(&col[j], points[k * (g - d) + j]);
        }
        acb_mat_vector_mul_col(new_z, star, col, prec);
        _acb_vec_add(new_z, new_z, z, g - d, prec);
        acb_theta_ql_tree_init_rec(T, cuts, nb_cuts - 1, d, new_z, cho, tau, prec);
    }
    
    arb_mat_clear(Y);
    arb_mat_clear(Yinv);
    acb_mat_clear(star);
    _acb_vec_clear(new_z, g - d);
    _acb_vec_clear(col, d);
    arf_clear(eps);
    arf_clear(R2);
    _arb_vec_clear(offset, d);
    flint_free(points);    
}

void
acb_theta_ql_tree_init(acb_theta_ql_tree_t T, acb_srcptr z,
    const acb_mat_t tau, slong prec)
{
    slong g = acb_mat_nrows(tau);
    slong* cuts;
    slong nb_cuts;
    arb_mat_t cho;

    arb_mat_init(cho, g, g);
    cuts = flint_malloc(g * sizeof(slong));

    acb_mat_get_imag(cho, tau);
    arb_mat_cho(cho, cho, prec);
    arb_mat_transpose(cho, cho);

    nb_cuts = acb_theta_ql_cuts(cuts, cho, prec);
    acb_theta_ql_tree_init_rec(T, cuts, nb_cuts, g, z, cho, tau, prec);
    
    arb_mat_clear(cho);
    flint_free(cuts);
}
