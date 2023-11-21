/*
    Copyright (C) 2023 Jean Kieffer

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb_mat.h"
#include "acb_theta.h"

#define ACB_THETA_ELD_MAX_PTS 1000000
#define ACB_THETA_ELD_MAX_ERR 100

static void
slong_vec_max(slong * r, slong * v1, slong * v2, slong d)
{
    slong k;
    for (k = 0; k < d; k++)
    {
        r[k] = FLINT_MAX(v1[k], v2[k]);
    }
}

static int
arf_get_si_safe(slong * m, const arf_t x, arf_rnd_t rnd)
{
    if (!arf_is_finite(x))
    {
        return 0;
    }
    else if (arf_cmpabs_2exp_si(x, FLINT_BITS - 4) > 0)
    {
        return 0;
    }
    else
    {
        *m = arf_get_si(x, rnd);
        return 1;
    }
}

static int
acb_theta_eld_interval(slong * min, slong * mid, slong * max, const arb_t ctr, const arf_t rad)
{
    slong lp = ACB_THETA_LOW_PREC;
    slong e;
    arb_t y;
    arf_t b;
    int res;

    arb_init(y);
    arf_init(b);

    arf_set_mag(b, arb_radref(ctr));
    res = arf_get_si_safe(&e, b, ARF_RND_NEAR);
    if (res)
    {
        res = (e <= ACB_THETA_ELD_MAX_ERR);
    }

    res = res && arf_get_si_safe(mid, arb_midref(ctr), ARF_RND_NEAR);

    arb_set_arf(y, rad);
    arb_add(y, ctr, y, lp);
    arb_get_ubound_arf(b, y, lp);
    res = res && arf_get_si_safe(max, b, ARF_RND_FLOOR);

    arb_set_arf(y, rad);
    arb_sub(y, ctr, y, lp);
    arb_get_lbound_arf(b, y, lp);
    res = res && arf_get_si_safe(min, b, ARF_RND_CEIL);

    arb_clear(y);
    arf_clear(b);
    return res;
}

static void
acb_theta_eld_next_R2(arf_t next_R2, const arf_t R2, const arb_t gamma, const arb_t v, slong k)
{
    slong lp = ACB_THETA_LOW_PREC;
    arb_t x;
    arf_t sub;

    arb_init(x);
    arf_init(sub);

    /* Set next_R2 to R2 - (v + gamma*k)^2 */
    arb_mul_si(x, gamma, k, lp);
    arb_add(x, x, v, lp);
    arb_sqr(x, x, lp);

    arb_get_lbound_arf(sub, x, lp);
    arf_sub(next_R2, R2, sub, lp, ARF_RND_CEIL);

    arb_clear(x);
    arf_clear(sub);
}

static void
acb_theta_eld_init_children(acb_theta_eld_t E, slong nr, slong nl)
{
    slong d = acb_theta_eld_dim(E);
    slong g = acb_theta_eld_ambient_dim(E);
    slong k;

    if (nr > 0) /* should always be the case */
    {
        E->rchildren = flint_malloc(nr * sizeof(struct acb_theta_eld_struct));
        acb_theta_eld_nr(E) = nr;
        for (k = 0; k < nr; k++)
        {
            acb_theta_eld_init(acb_theta_eld_rchild(E, k), d - 1, g);
        }
    }
    if (nl > 0)
    {
        E->lchildren = flint_malloc(nl * sizeof(struct acb_theta_eld_struct));
        acb_theta_eld_nl(E) = nl;
        for (k = 0; k < nl; k++)
        {
            acb_theta_eld_init(acb_theta_eld_lchild(E, k), d - 1, g);
        }
    }
}

static int
acb_theta_eld_init_interval(acb_theta_eld_t E, const arb_mat_t C,
    const arf_t R2, arb_srcptr v, slong * last_coords)
{
    slong min, mid, max;
    slong d = acb_theta_eld_dim(E);
    slong g = acb_theta_eld_ambient_dim(E);
    slong lp = ACB_THETA_LOW_PREC;
    slong k;
    arb_t x, ctr;
    arf_t rad;
    int res;

    arb_init(x);
    arb_init(ctr);
    arf_init(rad);

    for (k = 0; k < g - d; k++)
    {
        E->last_coords[k] = last_coords[k];
    }

    if (arf_cmp_si(R2, 0) < 0)
    {
        arf_zero(rad);
    }
    else
    {
        arb_set_arf(x, R2);
        arb_sqrt(x, x, lp);
        arb_div(x, x, arb_mat_entry(C, d - 1, d - 1), lp);
        arb_get_ubound_arf(rad, x, lp);
    }

    arb_div(ctr, &v[d - 1], arb_mat_entry(C, d - 1, d - 1), lp);
    arb_neg(ctr, ctr);
    res = acb_theta_eld_interval(&min, &mid, &max, ctr, rad);

    if (res)
    {
        acb_theta_eld_min(E) = min;
        acb_theta_eld_mid(E) = mid;
        acb_theta_eld_max(E) = max;
    }

    arb_clear(x);
    arb_clear(ctr);
    arf_clear(rad);
    return res;
}

/* Main recursive function in dimension d>1 */

static int
acb_theta_eld_set_rec(acb_theta_eld_t E, const arb_mat_t C,
    const arf_t R2, arb_srcptr v, slong * last_coords)
{
    slong d = acb_theta_eld_dim(E);
    slong g = acb_theta_eld_ambient_dim(E);
    slong lp = ACB_THETA_LOW_PREC;
    slong min, mid, max, k;
    arf_t next_R2;
    slong *next_coords;
    arb_ptr v_diff;
    arb_ptr v_mid;
    arb_ptr next_v;
    slong c;
    slong nr, nl;
    int res;

    res = acb_theta_eld_init_interval(E, C, R2, v, last_coords);
    if (!res)
    {
        return 0;
    }
    min = acb_theta_eld_min(E);
    mid = acb_theta_eld_mid(E);
    max = acb_theta_eld_max(E);

    /* Induction only if d > 1 and min <= max */
    if (min > max)
    {
        acb_theta_eld_nb_pts(E) = 0;
        if (d == 1)
        {
            acb_theta_eld_nb_border(E) = 2;
        }
        else
        {
            acb_theta_eld_nb_border(E) = 0;
        }
        for (k = 0; k < d; k++)
        {
            acb_theta_eld_box(E, k) = 0;
        }
        return 1;
    }
    else if (d == 1)
    {
        acb_theta_eld_nb_pts(E) = max - min + 1;
        acb_theta_eld_nb_border(E) = 2;
        acb_theta_eld_box(E, 0) = FLINT_MAX(max, -min);
        return (acb_theta_eld_nb_pts(E) <= ACB_THETA_ELD_MAX_PTS);
    }

    /* Begin main function */
    arf_init(next_R2);
    next_coords = flint_malloc((g - d + 1) * sizeof(slong));
    v_diff = _arb_vec_init(d - 1);
    v_mid = _arb_vec_init(d - 1);
    next_v = _arb_vec_init(d - 1);

    /* Initialize children */
    nr = max - mid + 1;
    nl = mid - min;
    acb_theta_eld_init_children(E, nr, nl);

    /* Set v_mid, v_diff */
    for (k = 0; k < d - 1; k++)
    {
        arb_set(&v_diff[k], arb_mat_entry(C, k, d - 1));
        arb_mul_si(&v_mid[k], &v_diff[k], mid, lp);
    }
    _arb_vec_add(v_mid, v_mid, v, d - 1, lp);
    for (k = 0; k < g - d; k++)
    {
        next_coords[k + 1] = last_coords[k];
    }

    /* Set children recursively */
    acb_theta_eld_nb_pts(E) = 0;
    acb_theta_eld_nb_border(E) = 0;
    acb_theta_eld_box(E, d - 1) = FLINT_MAX(max, -min);
    for (k = 0; k < d - 1; k++)
    {
        acb_theta_eld_box(E, k) = 0;
    }

    /* Right loop */
    _arb_vec_set(next_v, v_mid, d - 1);
    for (k = 0; (k < nr) && res; k++)
    {
        if (k > 0)
        {
            _arb_vec_add(next_v, next_v, v_diff, d - 1, lp);
        }

        c = mid + k;
        acb_theta_eld_next_R2(next_R2, R2, arb_mat_entry(C, d - 1, d - 1), &v[d - 1], c);
        next_coords[0] = c;
        res = acb_theta_eld_set_rec(acb_theta_eld_rchild(E, k), C, next_R2,
            next_v, next_coords);
        if (res)
        {
            acb_theta_eld_nb_pts(E) += acb_theta_eld_nb_pts(acb_theta_eld_rchild(E, k));
            acb_theta_eld_nb_border(E) += acb_theta_eld_nb_border(acb_theta_eld_rchild(E, k));
            slong_vec_max(E->box, E->box, acb_theta_eld_rchild(E, k)->box, d - 1);
            res = (acb_theta_eld_nb_pts(E) <= ACB_THETA_ELD_MAX_PTS);
        }
    }

    /* Left loop */
    _arb_vec_set(next_v, v_mid, d - 1);
    for (k = 0; (k < nl) && res; k++)
    {
        _arb_vec_sub(next_v, next_v, v_diff, d - 1, lp);

        c = mid - (k + 1);
        acb_theta_eld_next_R2(next_R2, R2, arb_mat_entry(C, d - 1, d - 1), &v[d - 1], c);
        next_coords[0] = c;
        res = acb_theta_eld_set_rec(acb_theta_eld_lchild(E, k), C, next_R2,
            next_v, next_coords);

        if (res) /* we expect this always holds */
        {
            acb_theta_eld_nb_pts(E) += acb_theta_eld_nb_pts(acb_theta_eld_lchild(E, k));
            acb_theta_eld_nb_border(E) += acb_theta_eld_nb_border(acb_theta_eld_lchild(E, k));
            slong_vec_max(E->box, E->box, acb_theta_eld_lchild(E, k)->box, d - 1);
            res = (acb_theta_eld_nb_pts(E) <= ACB_THETA_ELD_MAX_PTS);
        }
    }

    arf_clear(next_R2);
    flint_free(next_coords);
    _arb_vec_clear(v_diff, d - 1);
    _arb_vec_clear(v_mid, d - 1);
    _arb_vec_clear(next_v, d - 1);
    return res;
}

int
acb_theta_eld_set(acb_theta_eld_t E, const arb_mat_t C, const arf_t R2, arb_srcptr v)
{
    slong d = acb_theta_eld_dim(E);
    slong g = acb_theta_eld_ambient_dim(E);

    acb_theta_eld_clear(E);
    acb_theta_eld_init(E, d, g);
    return acb_theta_eld_set_rec(E, C, R2, v, NULL);
}
