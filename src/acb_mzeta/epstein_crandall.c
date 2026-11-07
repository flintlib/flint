/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdlib.h>
#include "arb_mat.h"
#include "acb_theta.h"
#include "impl.h"

/* Crandall's representation (Buchheit, Busse, Gutendorf, Thm. 2.11):

     Z = (pi/lam^2)^(nu/2) / Gamma(nu/2) [ sum_{z in L} G_nu((z - x)/lam) e(-y.z)
         + lam^d / V sum_{k in L*} G_{d-nu}(lam (k + y)) e(-x.(k + y)) ]

   with G_nu(u) = t^(-nu/2) Gamma(nu/2, t), t = pi |u|^2.  Both sums run over
   the points of ellipsoids |C n + v|^2 <= T (acb_theta_eld), T from the
   theta tail bound with |G| <= e^-t / (t - v).  With an exact Gram matrix
   the quadratic form is exact and the points are grouped by its value
   (shells), so G is evaluated once per shell; +-n are paired when the
   center is a lattice point.  Each G is computed to a relative precision
   tapered by its magnitude e^-t. */

typedef struct
{
    arb_ptr q;
    acb_ptr ph;
    slong len, alloc;
}
shells_struct;

static void
shells_push(shells_struct * S, const arb_t q, const acb_t ph)
{
    if (S->len == S->alloc)
    {
        slong na = FLINT_MAX(16, 2 * S->alloc), i;
        arb_ptr nq = _arb_vec_init(na);
        acb_ptr np = _acb_vec_init(na);
        for (i = 0; i < S->len; i++)
        {
            arb_swap(nq + i, S->q + i);
            acb_swap(np + i, S->ph + i);
        }
        if (S->alloc)
        {
            _arb_vec_clear(S->q, S->alloc);
            _acb_vec_clear(S->ph, S->alloc);
        }
        S->q = nq;
        S->ph = np;
        S->alloc = na;
    }
    arb_set(S->q + S->len, q);
    acb_set(S->ph + S->len, ph);
    S->len++;
}

typedef struct
{
    arf_srcptr key;
    slong idx;
}
sort_entry;

static int
cmp_entry(const void * a, const void * b)
{
    return arf_cmp(((const sort_entry *) a)->key, ((const sort_entry *) b)->key);
}

/* merge the entries with equal exact q (if all q are exact) */
static void
shells_merge(shells_struct * S, slong prec)
{
    slong i, n = S->len, k;
    sort_entry * idx;
    arb_ptr nq;
    acb_ptr np;

    for (i = 0; i < n; i++)
        if (!arb_is_exact(S->q + i))
            return;
    if (n < 2)
        return;

    idx = flint_malloc(n * sizeof(sort_entry));
    for (i = 0; i < n; i++)
    {
        idx[i].key = arb_midref(S->q + i);
        idx[i].idx = i;
    }
    qsort(idx, n, sizeof(sort_entry), cmp_entry);

    nq = _arb_vec_init(n);
    np = _acb_vec_init(n);
    k = 0;
    for (i = 0; i < n; i++)
    {
        if (k > 0 && arb_equal(nq + k - 1, S->q + idx[i].idx))
            acb_add(np + k - 1, np + k - 1, S->ph + idx[i].idx, prec);
        else
        {
            arb_set(nq + k, S->q + idx[i].idx);
            acb_set(np + k, S->ph + idx[i].idx);
            k++;
        }
    }
    _arb_vec_clear(S->q, S->alloc);
    _acb_vec_clear(S->ph, S->alloc);
    S->q = nq;
    S->ph = np;
    S->alloc = n;
    S->len = k;
    flint_free(idx);
}

/* v^T M v */
static void
quad_form(arb_t res, const arb_mat_t M, arb_srcptr v, slong d, slong prec)
{
    slong i, j;
    arb_t s;
    arb_init(s);
    arb_zero(res);
    for (i = 0; i < d; i++)
    {
        arb_zero(s);
        for (j = i + 1; j < d; j++)
            arb_addmul(s, arb_mat_entry(M, i, j), v + j, prec);
        arb_mul_2exp_si(s, s, 1);
        arb_addmul(s, arb_mat_entry(M, i, i), v + i, prec);
        arb_addmul(res, s, v + i, prec);
    }
    arb_clear(s);
}

/* sum_n G_{2a}(s (n + c)^T M (n + c)) e(-p.(n + pofs)) over the points of
   |C (n + c)|^2 <= T, plus the tail bound as a radius */
static int
lattice_sum(acb_t res, const _acb_mzeta_epstein_problem_struct * P, int side,
    const acb_t a, const acb_t ga, const arb_mat_t C, const arb_mat_t M, const arb_t s,
    arb_srcptr c, arb_srcptr p, arb_srcptr pofs, const arb_t K, double v,
    slong bits, slong wp)
{
    slong d = arb_mat_nrows(C);
    slong i, j, n, nb;
    slong * pts;
    arb_ptr vc, nc;
    arb_t q, tq, ang, tmp;
    acb_t ph, g;
    arf_t R2;
    mag_t tail;
    double T;
    int pair, p_zero, ok = 1;
    shells_struct S;
    acb_theta_eld_t E;

    vc = _arb_vec_init(d);
    nc = _arb_vec_init(d);
    arb_init(q); arb_init(tq); arb_init(ang); arb_init(tmp);
    acb_init(ph); acb_init(g);
    arf_init(R2);
    mag_init(tail);
    acb_theta_eld_init(E, d, d);
    S.q = NULL; S.ph = NULL; S.len = S.alloc = 0;

    /* +-n give the same t when c = 0 (the omitted points are then
       symmetric, except the regularisation point) */
    pair = _arb_vec_is_zero(c, d) && !(side == 1 && P->reg);
    p_zero = _arb_vec_is_zero(p, d) && _arb_vec_is_zero(pofs, d);

    T = _acb_mzeta_epstein_choose_T(K, d, v, bits);
    if (!(T < 1e9))
    {
        ok = 0;
        goto cleanup;
    }
    _acb_mzeta_epstein_tail_bound(tail, K, d, T, v);
    arf_set_d(R2, T);

    for (i = 0; i < d; i++)
        for (j = i; j < d; j++)
            arb_addmul(vc + i, arb_mat_entry(C, i, j), c + j, wp);

    if (!acb_theta_eld_set(E, C, R2, vc))
    {
        ok = 0;
        goto cleanup;
    }

    nb = acb_theta_eld_nb_pts(E);
    pts = flint_malloc(FLINT_MAX(nb, 1) * d * sizeof(slong));
    acb_theta_eld_points(pts, E);

    for (n = 0; n < nb; n++)
    {
        slong * pt = pts + n * d;
        int iszero = 1, sign = 0;

        for (i = 0; i < d; i++)
        {
            if (pt[i] != 0)
            {
                iszero = 0;
                if (sign == 0)
                    sign = (pt[i] > 0) ? 1 : -1;
            }
        }
        if (_acb_mzeta_epstein_skip(P, side, pt))
            continue;
        if (pair && sign < 0)
            continue;

        for (i = 0; i < d; i++)
            arb_add_si(nc + i, c + i, pt[i], wp);
        quad_form(q, M, nc, d, wp);

        if (p_zero)
            acb_one(ph);
        else
        {
            arb_zero(ang);
            for (i = 0; i < d; i++)
            {
                arb_add_si(tmp, pofs + i, pt[i], wp);
                arb_addmul(ang, p + i, tmp, wp);
            }
            arb_mul_2exp_si(ang, ang, 1);
            if (pair && !iszero)
            {
                /* pofs = 0 here: e(-p.n) + e(p.n) = 2 cos(2 pi p.n) */
                arb_cos_pi(acb_realref(ph), ang, wp);
                arb_zero(acb_imagref(ph));
            }
            else
            {
                arb_neg(ang, ang);
                arb_sin_cos_pi(acb_imagref(ph), acb_realref(ph), ang, wp);
            }
        }
        if (pair && !iszero)
            acb_mul_2exp_si(ph, ph, 1);

        shells_push(&S, q, ph);
    }
    flint_free(pts);

    shells_merge(&S, wp);

    acb_zero(res);
    for (n = 0; n < S.len; n++)
    {
        slong tp;
        double tm;
        if (acb_is_zero(S.ph + n))
            continue;
        arb_mul(tq, S.q + n, s, wp);
        tm = arf_get_d(arb_midref(tq), ARF_RND_DOWN);
        if (!(tm >= 0.0 && tm < 1e15))
            tm = 0.0;
        tp = FLINT_MAX(wp - (slong) (tm * 1.4426950408889634), 32);
        _acb_mzeta_gamma_upper_scaled(g, a, tq, ga, tp);
        acb_addmul(res, g, S.ph + n, wp);
    }
    acb_add_error_mag(res, tail);

cleanup:
    _arb_vec_clear(vc, d);
    _arb_vec_clear(nc, d);
    arb_clear(q); arb_clear(tq); arb_clear(ang); arb_clear(tmp);
    acb_clear(ph); acb_clear(g);
    arf_clear(R2);
    mag_clear(tail);
    acb_theta_eld_clear(E);
    if (S.alloc)
    {
        _arb_vec_clear(S.q, S.alloc);
        _acb_vec_clear(S.ph, S.alloc);
    }
    return ok;
}

int
acb_mzeta_epstein_crandall(acb_t res, const acb_t nu, const arb_mat_t A,
    arb_srcptr x, arb_srcptr y, int flags, acb_mzeta_cache_t cache, slong prec)
{
    slong d = arb_mat_nrows(A);
    slong wp = prec + 16 + 2 * FLINT_BIT_COUNT(prec);
    const _acb_mzeta_lattice_struct * L;
    _acb_mzeta_epstein_problem_struct P;
    acb_mzeta_cache_t tmp_cache;
    acb_t Sd, Sr, ga, gb;
    arb_t s;
    arb_ptr cneg, zero;
    int ok;

    if (d < 1 || arb_mat_ncols(A) != d)
    {
        acb_indeterminate(res);
        return 0;
    }

    if (cache == NULL)
    {
        acb_mzeta_cache_init(tmp_cache);
        ok = acb_mzeta_epstein_crandall(res, nu, A, x, y, flags, tmp_cache, prec);
        acb_mzeta_cache_clear(tmp_cache);
        return ok;
    }

    L = _acb_mzeta_lattice_get(cache, A, wp);
    if (L == NULL)
    {
        acb_indeterminate(res);
        return 1;
    }

    if (!_acb_mzeta_epstein_problem_init(&P, L, nu, x, y, flags, wp))
    {
        _acb_mzeta_epstein_problem_clear(&P);
        acb_indeterminate(res);
        return 1;
    }

    acb_init(Sd); acb_init(Sr); acb_init(ga); acb_init(gb);
    arb_init(s);
    cneg = _arb_vec_init(d);
    zero = _arb_vec_init(d);

    acb_gamma(ga, P.a, wp + 32);
    acb_gamma(gb, P.b, wp + 32);

    /* direct: t = (pi/lam^2) (n - u0)^T G (n - u0), phase e(-w0.n) */
    _arb_vec_neg(cneg, P.u0, d);
    arb_const_pi(s, wp);
    arb_div(s, s, L->lam2, wp);
    ok = lattice_sum(Sd, &P, 0, P.a, ga, L->Cd, L->G, s, cneg, P.w0, zero,
        L->Kd, P.v0, prec + 10, wp);

    /* reciprocal: t = pi lam^2 (m + w0)^T G^-1 (m + w0), phase e(-u0.(m + w0)) */
    arb_const_pi(s, wp);
    arb_mul(s, s, L->lam2, wp);
    ok = ok && lattice_sum(Sr, &P, 1, P.b, gb, L->Cr, L->Ginv, s, P.w0, P.u0, P.w0,
        L->Kr, P.v1, prec + 10, wp);

    if (ok)
        _acb_mzeta_epstein_finish(res, Sd, Sr, L, &P, prec, wp);
    else
        acb_indeterminate(res);

    acb_clear(Sd); acb_clear(Sr); acb_clear(ga); acb_clear(gb);
    arb_clear(s);
    _arb_vec_clear(cneg, d);
    _arb_vec_clear(zero, d);
    _acb_mzeta_epstein_problem_clear(&P);
    return 1;
}
