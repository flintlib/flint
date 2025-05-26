/*
    Copyright (C) 2023 Jean Kieffer

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb.h"
#include "arb_mat.h"
#include "acb.h"
#include "acb_mat.h"
#include "acb_theta.h"

/* Should not be used in tests */
static void
acb_theta_dist_unif(arb_t d, const arb_mat_t cho, slong prec)
{
    slong g = arb_mat_nrows(cho);
    arb_ptr v;
    slong k;

    v = _arb_vec_init(g);

    for (k = 0; k < g; k++)
    {
        arb_zero_pm_one(&v[k]);
        arb_mul_2exp_si(&v[k], &v[k], -1);
    }
    arb_mat_vector_mul_col(v, cho, v, prec);
    arb_dot(d, NULL, 0, v, 1, v, 1, g, prec);

    _arb_vec_clear(v, g);
}

static void
acb_theta_dist_pt(arb_t d, arb_srcptr v, const arb_mat_t cho, const slong * n, slong prec)
{
    slong g = arb_mat_nrows(cho);
    arb_ptr w;
    slong k;

    w = _arb_vec_init(g);

    for (k = 0; k < g; k++)
    {
        arb_set_si(&w[k], n[k]);
    }
    arb_mat_vector_mul_col(w, cho, w, prec);
    _arb_vec_add(w, w, v, g, prec);
    arb_dot(d, NULL, 0, w, 1, w, 1, g, prec);

    _arb_vec_clear(w, g);
}

static void
acb_theta_dist_ubound(arf_t u, arb_srcptr v, const arb_mat_t cho, slong prec)
{
    slong g = acb_mat_nrows(cho);
    slong nb = 1 << g;
    arb_mat_t choinv;
    arb_ptr x;
    slong * approx;
    slong * pt;
    arb_t d;
    arf_t b;
    slong j, k;
    int r = 1;

    arb_mat_init(choinv, g, g);
    x = _arb_vec_init(g);
    approx = flint_malloc(2 * g * sizeof(slong));
    pt = flint_malloc(g * sizeof(slong));
    arb_init(d);
    arf_init(b);

    arb_mat_one(choinv);
    arb_mat_solve_triu(choinv, cho, choinv, 0, prec);
    arb_mat_vector_mul_col(x, choinv, v, prec);
    r = _arb_vec_is_finite(x, g);

    for (k = 0; (k < g) && r; k++)
    {
        r = (arf_cmpabs_2exp_si(arb_midref(&x[k]), 30) <= 0);
        if (r)
        {
            approx[2 * k] = - arf_get_si(arb_midref(&x[k]), ARF_RND_FLOOR);
            approx[2 * k + 1] = - arf_get_si(arb_midref(&x[k]), ARF_RND_CEIL);
        }
    }

    arf_pos_inf(u);
    if (r)
    {
        for (k = 0; k < nb; k++)
        {
            for (j = 0; j < g; j++)
            {
                pt[j] = approx[2 * j + (k & (1 << j) ? 0 : 1)];
            }
            acb_theta_dist_pt(d, v, cho, pt, prec);
            arb_get_ubound_arf(b, d, prec);
            arf_min(u, u, b);
        }
    }

    arb_mat_clear(choinv);
    _arb_vec_clear(x, g);
    flint_free(approx);
    flint_free(pt);
    arb_clear(d);
    arf_clear(b);
}

static void
acb_theta_dist_lat(arb_t d, arb_srcptr v, const arb_mat_t cho, slong prec)
{
    slong g = arb_mat_nrows(cho);
    acb_theta_eld_t E;
    slong nb;
    slong * pts;
    arf_t u;
    arb_t x;
    slong k;
    int b;

    acb_theta_eld_init(E, g, g);
    arf_init(u);
    arb_init(x);

    acb_theta_dist_ubound(u, v, cho, prec);
    b = acb_theta_eld_set(E, cho, u, v);

    if (b)
    {
        nb = acb_theta_eld_nb_pts(E);
        pts = flint_malloc(nb * g * sizeof(slong));
        acb_theta_eld_points(pts, E);

        arb_pos_inf(d);
        for (k = 0; k < nb; k++)
        {
            acb_theta_dist_pt(x, v, cho, pts + k * g, prec);
            arb_min(d, d, x, prec);
        }

        flint_free(pts);
    }
    else
    {
        /* Should not happen in tests */
        acb_theta_dist_unif(d, cho, prec);
    }
    arb_nonnegative_part(d, d);

    acb_theta_eld_clear(E);
    arf_clear(u);
    arb_clear(x);
}

void
acb_theta_eld_distances(arb_ptr ds, acb_srcptr zs, slong nb,
    const acb_mat_t tau, slong prec)
{
    slong g = acb_mat_nrows(tau);
    slong n = 1 << g;
    arb_mat_t yinv, cho;
    arb_ptr v, w;
    ulong a;
    slong j;

    arb_mat_init(yinv, g, g);
    arb_mat_init(cho, g, g);
    v = _arb_vec_init(g);
    w = _arb_vec_init(g);

    acb_siegel_cho_yinv(cho, yinv, tau, prec);

    for (j = 0; j < nb; j++)
    {
        _acb_vec_get_imag(v, zs + j * g, g);
        arb_mat_vector_mul_col(v, yinv, v, prec);

        for (a = 0; a < n; a++)
        {
            acb_theta_char_get_arb(w, a, g);
            _arb_vec_add(w, v, w, g, prec);
            arb_mat_vector_mul_col(w, cho, w, prec);
            acb_theta_dist_lat(&ds[j * n + a], w, cho, prec);
        }
    }

    arb_mat_clear(yinv);
    arb_mat_clear(cho);
    _arb_vec_clear(v, g);
    _arb_vec_clear(w, g);
}
