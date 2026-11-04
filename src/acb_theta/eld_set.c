/*
    Copyright (C) 2023 Jean Kieffer

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb.h"
#include "double_interval.h"
#include "arb_mat.h"
#include "acb_theta.h"

#define ACB_THETA_ELD_MAX_PTS 10000000
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
    slong d = E->dim;
    slong g = E->ambient_dim;
    slong k;

    if (nr > 0) /* should always be the case */
    {
        E->rchildren = flint_malloc(nr * sizeof(struct acb_theta_eld_struct));
        E->nr = nr;
        for (k = 0; k < nr; k++)
        {
            acb_theta_eld_init(&E->rchildren[k], d - 1, g);
        }
    }
    if (nl > 0)
    {
        E->lchildren = flint_malloc(nl * sizeof(struct acb_theta_eld_struct));
        E->nl = nl;
        for (k = 0; k < nl; k++)
        {
            acb_theta_eld_init(&E->lchildren[k], d - 1, g);
        }
    }
}

static int
acb_theta_eld_init_interval(acb_theta_eld_t E, const arb_mat_t C,
    const arf_t R2, arb_srcptr v, slong * last_coords)
{
    slong min, mid, max;
    slong d = E->dim;
    slong g = E->ambient_dim;
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
        E->min = min;
        E->mid = mid;
        E->max = max;
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
    slong d = E->dim;
    slong g = E->ambient_dim;
    slong min, mid, max;
    slong lp = ACB_THETA_LOW_PREC;
    arf_t next_R2;
    slong *next_coords;
    arb_ptr v_diff;
    arb_ptr v_mid;
    arb_ptr next_v;
    slong c, nr, nl, k;
    int res;

    res = acb_theta_eld_init_interval(E, C, R2, v, last_coords);
    if (!res)
    {
        return 0;
    }
    min = E->min;
    mid = E->mid;
    max = E->max;

    /* Induction only if d > 1 and min <= max */
    if (min > max)
    {
        E->nb_pts = 0;
        if (d == 1)
        {
            E->nb_border = 2;
        }
        else
        {
            E->nb_border = 0;
        }
        for (k = 0; k < d; k++)
        {
            E->box[k] = 0;
        }
        return 1;
    }
    else if (d == 1)
    {
        E->nb_pts = max - min + 1;
        E->nb_border = 2;
        E->box[0] = FLINT_MAX(max, -min);
        return (E->nb_pts <= ACB_THETA_ELD_MAX_PTS);
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
    E->nb_pts = 0;
    E->nb_border = 0;
    E->box[d - 1] = FLINT_MAX(max, -min);
    for (k = 0; k < d - 1; k++)
    {
        E->box[k] = 0;
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
        res = acb_theta_eld_set_rec(&E->rchildren[k], C, next_R2, next_v, next_coords);
        if (res)
        {
            E->nb_pts += (&E->rchildren[k])->nb_pts;
            E->nb_border += (&E->rchildren[k])->nb_border;
            slong_vec_max(E->box, E->box, (&E->rchildren[k])->box, d - 1);
            res = ((E->nb_pts) <= ACB_THETA_ELD_MAX_PTS);
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
        res = acb_theta_eld_set_rec(&E->lchildren[k], C, next_R2, next_v, next_coords);

        if (res)
        {
            E->nb_pts += (&E->lchildren[k])->nb_pts;
            E->nb_border += (&E->lchildren[k])->nb_border;
            slong_vec_max(E->box, E->box, (&E->lchildren[k])->box, d - 1);
            res = ((E->nb_pts) <= ACB_THETA_ELD_MAX_PTS);
        }
    }

    arf_clear(next_R2);
    flint_free(next_coords);
    _arb_vec_clear(v_diff, d - 1);
    _arb_vec_clear(v_mid, d - 1);
    _arb_vec_clear(next_v, d - 1);
    return res;
}

/* The same recursion in double interval arithmetic (outward rounded), used
   when all the data converts to finite double intervals of moderate size.
   It encloses the same ellipsoid, possibly with a few more points on its
   border, and avoids the cost of arb arithmetic at every node, which
   dominates at low precision (the enumeration is otherwise more expensive
   than the theta summation itself). */

#define ELD_D_MAX_COORD ((FLINT_BITS == 64) ? 0x1p56 : 0x1p28)

static int
acb_theta_eld_interval_d(slong * min, slong * mid, slong * max, di_t ctr, double rad)
{
    double lo, hi, m;

    if (!(ctr.b - ctr.a <= 2.0 * ACB_THETA_ELD_MAX_ERR))
        return 0;
    m = 0.5 * (ctr.a + ctr.b);
    lo = _di_below(ctr.a - rad);
    hi = _di_above(ctr.b + rad);
    if (!(lo > -ELD_D_MAX_COORD && hi < ELD_D_MAX_COORD))
        return 0;
    *mid = (slong) nearbyint(m);
    *min = (slong) ceil(lo);
    *max = (slong) floor(hi);
    return 1;
}

static int
acb_theta_eld_set_rec_d(acb_theta_eld_t E, const di_t * C, double R2,
    const di_t * v, slong * last_coords)
{
    slong d = E->dim;
    slong g = E->ambient_dim;
    slong min, mid, max, nr, nl, c, k, j;
    di_t gam = C[(d - 1) * g + (d - 1)];
    di_t ctr, x;
    double rad;
    slong * next_coords;
    di_t * next_v;
    int res = 1;

    for (k = 0; k < g - d; k++)
        E->last_coords[k] = last_coords[k];

    /* interval for the last coordinate: |gam n + v_{d-1}| <= sqrt(R2) */
    if (!(gam.a > 0.0))
        return 0;
    rad = (R2 <= 0.0) ? 0.0 : _di_above(_di_above(sqrt(R2)) / gam.a);
    ctr = di_neg(di_fast_div(v[d - 1], gam));
    if (!acb_theta_eld_interval_d(&min, &mid, &max, ctr, rad))
        return 0;
    E->min = min;
    E->mid = mid;
    E->max = max;

    if (min > max)
    {
        E->nb_pts = 0;
        E->nb_border = (d == 1) ? 2 : 0;
        for (k = 0; k < d; k++)
            E->box[k] = 0;
        return 1;
    }
    else if (d == 1)
    {
        E->nb_pts = max - min + 1;
        E->nb_border = 2;
        E->box[0] = FLINT_MAX(max, -min);
        return (E->nb_pts <= ACB_THETA_ELD_MAX_PTS);
    }

    nr = max - mid + 1;
    nl = mid - min;
    acb_theta_eld_init_children(E, nr, nl);

    next_coords = flint_malloc((g - d + 1) * sizeof(slong));
    next_v = flint_malloc((d - 1) * sizeof(di_t));
    for (k = 0; k < g - d; k++)
        next_coords[k + 1] = last_coords[k];

    E->nb_pts = 0;
    E->nb_border = 0;
    E->box[d - 1] = FLINT_MAX(max, -min);
    for (k = 0; k < d - 1; k++)
        E->box[k] = 0;

    for (j = 0; (j < nr + nl) && res; j++)
    {
        struct acb_theta_eld_struct * child;
        double next_R2;

        if (j < nr)
        {
            c = mid + j;
            child = &E->rchildren[j];
        }
        else
        {
            c = mid - (j - nr + 1);
            child = &E->lchildren[j - nr];
        }

        /* next_R2 = R2 - (v_{d-1} + gam c)^2, rounded up */
        x = di_fast_add(v[d - 1], di_fast_mul_d(gam, (double) c));
        next_R2 = _di_above(R2 - di_fast_sqr(x).a);
        /* next_v = v + c C[0..d-2][d-1] */
        for (k = 0; k < d - 1; k++)
            next_v[k] = di_fast_add(v[k], di_fast_mul_d(C[k * g + (d - 1)], (double) c));
        next_coords[0] = c;

        res = acb_theta_eld_set_rec_d(child, C, next_R2, next_v, next_coords);
        if (res)
        {
            E->nb_pts += child->nb_pts;
            E->nb_border += child->nb_border;
            slong_vec_max(E->box, E->box, child->box, d - 1);
            res = (E->nb_pts <= ACB_THETA_ELD_MAX_PTS);
        }
    }

    flint_free(next_coords);
    flint_free(next_v);
    return res;
}

static int
acb_theta_eld_set_d(acb_theta_eld_t E, const arb_mat_t C, const arf_t R2, arb_srcptr v)
{
    slong g = E->ambient_dim;
    slong i, j;
    di_t * Cd, * vd;
    double r2;
    int res = 1;

    if (!arf_is_finite(R2) || arf_cmpabs_2exp_si(R2, 900) > 0)
        return 0;
    r2 = arf_get_d(R2, ARF_RND_UP);

    Cd = flint_malloc(g * g * sizeof(di_t));
    vd = flint_malloc(g * sizeof(di_t));
    for (i = 0; i < g && res; i++)
    {
        for (j = 0; j < g && res; j++)
        {
            Cd[i * g + j] = arb_get_di(arb_mat_entry(C, i, j));
            res = (fabs(Cd[i * g + j].a) < 1e150 && fabs(Cd[i * g + j].b) < 1e150);
        }
        vd[i] = arb_get_di(v + i);
        res = res && (fabs(vd[i].a) < 1e150 && fabs(vd[i].b) < 1e150);
    }

    if (res)
        res = acb_theta_eld_set_rec_d(E, Cd, r2, vd, NULL);

    flint_free(Cd);
    flint_free(vd);
    return res;
}

int
acb_theta_eld_set(acb_theta_eld_t E, const arb_mat_t C, const arf_t R2, arb_srcptr v)
{
    slong d = E->dim;
    slong g = E->ambient_dim;

    acb_theta_eld_clear(E);
    acb_theta_eld_init(E, d, g);
    if (d == g && acb_theta_eld_set_d(E, C, R2, v))
        return 1;

    acb_theta_eld_clear(E);
    acb_theta_eld_init(E, d, g);
    return acb_theta_eld_set_rec(E, C, R2, v, NULL);
}
