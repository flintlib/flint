/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_vec.h"
#include "arb_mat.h"
#include "acb_hypgeom.h"
#include "impl.h"

/* 1: all exactly zero, 0: some entry certainly nonzero, -1: undecided */
static int
vec_zero_status(arb_srcptr v, slong d)
{
    slong i;
    int exact = 1;
    for (i = 0; i < d; i++)
    {
        if (!arb_contains_zero(v + i))
            return 0;
        if (!arb_is_zero(v + i))
            exact = 0;
    }
    return exact ? 1 : -1;
}

/* v = v0 + m with m the rounded midpoint */
static int
split_int(arb_ptr v0, fmpz * m, arb_srcptr v, slong d, slong prec)
{
    slong i;
    for (i = 0; i < d; i++)
    {
        if (!arb_is_finite(v + i) || mag_cmp_2exp_si(arb_radref(v + i), -2) > 0
            || arf_cmpabs_2exp_si(arb_midref(v + i), FLINT_BITS - 4) > 0)
            return 0;
        arf_get_fmpz(m + i, arb_midref(v + i), ARF_RND_NEAR);
        arb_sub_fmpz(v0 + i, v + i, m + i, prec);
    }
    return 1;
}

int
_acb_mzeta_epstein_problem_init(_acb_mzeta_epstein_problem_struct * P,
    const _acb_mzeta_lattice_struct * L, const acb_t nu, arb_srcptr x, arb_srcptr y,
    int flags, slong prec)
{
    slong d = L->d, i, j;
    arb_ptr u;
    int ok;

    P->d = d;
    P->u0 = _arb_vec_init(d);
    P->w0 = _arb_vec_init(d);
    P->w = _arb_vec_init(d);
    P->m = _fmpz_vec_init(d);
    P->l = _fmpz_vec_init(d);
    P->lpt = flint_malloc(d * sizeof(slong));
    acb_init(P->nu);
    acb_init(P->a);
    acb_init(P->b);
    u = _arb_vec_init(d);

    acb_set(P->nu, nu);
    acb_mul_2exp_si(P->a, nu, -1);
    acb_set_si(P->b, d);
    acb_sub(P->b, P->b, nu, prec);
    acb_mul_2exp_si(P->b, P->b, -1);
    P->real = arb_is_zero(acb_imagref(nu));

    /* u = A^-1 x, w = A^T y */
    for (i = 0; i < d; i++)
    {
        for (j = 0; j < d; j++)
        {
            arb_addmul(u + i, arb_mat_entry(L->Ainv, i, j), x + j, prec);
            arb_addmul(P->w + i, arb_mat_entry(L->A, j, i), y + j, prec);
        }
    }

    ok = split_int(P->u0, P->m, u, d, prec) && split_int(P->w0, P->l, P->w, d, prec)
        && acb_is_finite(nu);

    if (ok)
    {
        P->xs = vec_zero_status(P->u0, d);
        P->ys = vec_zero_status(P->w0, d);
        /* a shift within the radius of a lattice point: the function is
           discontinuous (or singular) there */
        ok = (P->xs >= 0 && P->ys >= 0);
    }

    if (ok)
    {
        /* v = an upper bound for max(0, Re(a) - 1) over the ball (the
           parameter of the tail bound |G| <= e^-t / (t - v)), as a double
           rounded upward */
        arf_t ub;
        arf_init(ub);
        ok = arb_is_finite(acb_realref(P->a)) && arb_is_finite(acb_realref(P->b))
            && arf_cmpabs_2exp_si(arb_midref(acb_realref(P->a)), 20) < 0
            && arf_cmpabs_2exp_si(arb_midref(acb_realref(P->b)), 20) < 0
            && mag_cmp_2exp_si(arb_radref(acb_realref(P->a)), 20) < 0;
        if (ok)
        {
            arb_get_ubound_arf(ub, acb_realref(P->a), 64);
            arf_sub_ui(ub, ub, 1, 64, ARF_RND_UP);
            P->v0 = (arf_sgn(ub) > 0) ? arf_get_d(ub, ARF_RND_UP) : 0.0;
            arb_get_ubound_arf(ub, acb_realref(P->b), 64);
            arf_sub_ui(ub, ub, 1, 64, ARF_RND_UP);
            P->v1 = (arf_sgn(ub) > 0) ? arf_get_d(ub, ARF_RND_UP) : 0.0;
        }
        arf_clear(ub);

        /* the regularised function replaces the reciprocal term of y itself
           (m = l); at y = 0 it equals the singular term of Z */
        P->regflag = (flags & ACB_MZETA_EPSTEIN_REGULARIZED) != 0;
        P->reg = (flags & ACB_MZETA_EPSTEIN_REGULARIZED) && !(P->ys == 1 && _fmpz_vec_is_zero(P->l, d));
        for (i = 0; i < d; i++)
            P->lpt[i] = fmpz_get_si(P->l + i);
    }

    _arb_vec_clear(u, d);
    return ok;
}

void
_acb_mzeta_epstein_problem_clear(_acb_mzeta_epstein_problem_struct * P)
{
    slong d = P->d;
    _arb_vec_clear(P->u0, d);
    _arb_vec_clear(P->w0, d);
    _arb_vec_clear(P->w, d);
    _fmpz_vec_clear(P->m, d);
    _fmpz_vec_clear(P->l, d);
    flint_free(P->lpt);
    acb_clear(P->nu);
    acb_clear(P->a);
    acb_clear(P->b);
}

int
_acb_mzeta_epstein_skip(const _acb_mzeta_epstein_problem_struct * P, int side, const slong * n)
{
    slong i, d = P->d;
    int zero = 1, isl;

    for (i = 0; i < d; i++)
        zero = zero && (n[i] == 0);

    if (side == 0)
        return zero && P->xs == 1;

    if (zero && P->ys == 1)
        return 1;

    if (P->reg)
    {
        isl = 1;
        for (i = 0; i < d; i++)
            isl = isl && (n[i] == P->lpt[i]);
        return isl;
    }

    return 0;
}

double
_acb_mzeta_epstein_choose_T(const arb_t K, slong d, double v, slong bits)
{
    double T, target, lo;
    slong i;

    /* a parameter only: the bound itself is _acb_mzeta_epstein_tail_bound */
    if (!arb_is_finite(K) || arf_cmpabs_2exp_si(arb_midref(K), 1000) > 0)
        return INFINITY;
    target = bits * 0.6931471805599453 + log(FLINT_MAX(arf_get_d(arb_midref(K), ARF_RND_UP), 1.0) * 1.000001);
    lo = FLINT_MAX(4.0, v + 1.0);
    T = FLINT_MAX(lo, target);
    for (i = 0; i < 40; i++)
        T = FLINT_MAX(lo, target + 0.5 * (d - 1) * log(T) - log(T - v));
    return ceil(T * 1.0001 + 0.01);
}

/* The tail bound K T^((d-1)/2) e^-T / (T - v), valid for T >= 4 (the
   condition of the theta tail bound) and T > v; +inf otherwise.  T and v
   are doubles (exact), the bound is computed in ball arithmetic. */
void
_acb_mzeta_epstein_tail_bound(mag_t res, const arb_t K, slong d, double T, double v)
{
    arb_t t, u, w;
    slong prec = 64;

    if (!(T >= 4.0 && T > v && T < 1e15 && v >= 0.0) || !arb_is_finite(K))
    {
        mag_inf(res);
        return;
    }

    arb_init(t);
    arb_init(u);
    arb_init(w);

    arb_set_d(t, T);
    arb_log(u, t, prec);
    arb_mul_ui(u, u, d - 1, prec);
    arb_mul_2exp_si(u, u, -1);
    arb_sub(u, u, t, prec);
    arb_exp(u, u, prec);
    arb_mul(u, u, K, prec);
    arb_set_d(w, v);
    arb_sub(t, t, w, prec);
    if (arb_is_positive(t))
    {
        arb_div(u, u, t, prec);
        arb_get_mag(res, u);
    }
    else
        mag_inf(res);

    arb_clear(t);
    arb_clear(u);
    arb_clear(w);
}

void
_acb_mzeta_epstein_finish(acb_t res, const acb_t Sd0, const acb_t Sr0,
    const _acb_mzeta_lattice_struct * L, const _acb_mzeta_epstein_problem_struct * P,
    slong prec, slong wp)
{
    slong d = L->d, i;
    acb_t Sd, Sr, rg, t, u;
    arb_t c, pi, ang, tmp;

    acb_init(Sd); acb_init(Sr); acb_init(rg); acb_init(t); acb_init(u);
    arb_init(c); arb_init(pi); arb_init(ang); arb_init(tmp);

    acb_set(Sd, Sd0);
    acb_set(Sr, Sr0);
    arb_const_pi(pi, wp);

    /* regularisation: + (-g_{d-nu}(t_l)) e(-u0.w), t_l = pi lam^2 |A^-T w|^2,
       -g_b(t) = -Gamma(b) gamma*(b, t) */
    if (P->reg)
    {
        arb_t q;
        arb_init(q);
        arb_zero(q);
        for (i = 0; i < d; i++)
        {
            slong j;
            arb_zero(tmp);
            for (j = i + 1; j < d; j++)
                arb_addmul(tmp, arb_mat_entry(L->Ginv, i, j), P->w + j, wp);
            arb_mul_2exp_si(tmp, tmp, 1);
            arb_addmul(tmp, arb_mat_entry(L->Ginv, i, i), P->w + i, wp);
            arb_addmul(q, tmp, P->w + i, wp);
        }
        arb_mul(q, q, pi, wp);
        arb_mul(q, q, L->lam2, wp);
        acb_set_arb(t, q);
        acb_hypgeom_gamma_lower(u, P->b, t, 2, wp);
        acb_gamma(t, P->b, wp);
        acb_mul(u, u, t, wp);
        arb_zero(ang);
        for (i = 0; i < d; i++)
            arb_addmul(ang, P->u0 + i, P->w + i, wp);
        arb_mul_2exp_si(ang, ang, 1);
        arb_neg(ang, ang);
        arb_sin_cos_pi(acb_imagref(t), acb_realref(t), ang, wp);
        acb_submul(Sr, u, t, wp);
        arb_clear(q);
    }

    /* c = lam^d / V */
    arb_pow_ui(c, L->lam, d, wp);
    arb_div(c, c, L->vol, wp);

    acb_rgamma(rg, P->a, wp);
    acb_mul_arb(Sr, Sr, c, wp);
    acb_add(Sd, Sd, Sr, wp);
    acb_mul(Sd, Sd, rg, wp);

    /* x in L: G_nu(0) / Gamma(nu/2) = -1 / Gamma(nu/2 + 1) */
    if (P->xs == 1)
    {
        acb_add_ui(t, P->a, 1, wp);
        acb_rgamma(t, t, wp);
        acb_sub(Sd, Sd, t, wp);
    }

    /* y in L*: c G_{d-nu}(0) / Gamma(nu/2) = c (2 / (nu - d)) / Gamma(nu/2) */
    if (P->ys == 1)
    {
        acb_sub_si(t, P->nu, d, wp);
        acb_inv(t, t, wp);
        acb_mul_2exp_si(t, t, 1);
        acb_mul(t, t, rg, wp);
        acb_mul_arb(t, t, c, wp);
        acb_add(Sd, Sd, t, wp);
    }

    /* (pi / lam^2)^(nu/2) */
    arb_div(c, pi, L->lam2, wp);
    acb_set_arb(t, c);
    acb_pow(t, t, P->a, wp);
    acb_mul(Sd, Sd, t, wp);

    /* translation phase e(-w.m); for the regularised function also the
       factor e(x.y) = e(u.w): together e(u0.w) */
    arb_zero(ang);
    if (P->regflag)
    {
        for (i = 0; i < d; i++)
            arb_addmul(ang, P->u0 + i, P->w + i, wp);
    }
    else
    {
        for (i = 0; i < d; i++)
        {
            arb_mul_fmpz(tmp, P->w + i, P->m + i, wp);
            arb_sub(ang, ang, tmp, wp);
        }
    }
    if (!arb_is_zero(ang))
    {
        arb_mul_2exp_si(ang, ang, 1);
        arb_sin_cos_pi(acb_imagref(t), acb_realref(t), ang, wp);
        acb_mul(Sd, Sd, t, wp);
    }

    acb_set_round(res, Sd, prec);

    acb_clear(Sd); acb_clear(Sr); acb_clear(rg); acb_clear(t); acb_clear(u);
    arb_clear(c); arb_clear(pi); arb_clear(ang); arb_clear(tmp);
}
