/*
    Copyright (C) 2026 Fredrik Johansson
    Developed using Claude Fable 5.2

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/*
    Rigorous quadrature of the Euler-type integrals of the confluent
    and Gauss hypergeometric functions on a real interval,

        I = int_A^B exp(phi(t)) dt,
        phi(t) = lambda t + sum_{j < m} mu_j log(al_j + be_j t),

    with real parameters (balls), al_j + be_j t > 0 on (A, B), A = 0 and
    B finite or infinite, and mu_j >= 0 for a term which vanishes at an
    end of the interval.  (The points of the cells are doubles, which
    are dense at 0 only: an integral over (0, 1) with a term vanishing
    at 1 is computed as two integrals over (0, 1/2), in t and in 1 - t.)

    The interval is divided into cells.  On a cell [m - h, m + h], with
    d_j = m + al_j / be_j (the distance from m to the singularity of the
    term j),

        phi(m + u) - phi(m) = phi'(m) u
            + sum_j mu_j (log(1 + u / d_j) - u / d_j),

    in which only the terms of the second order and beyond remain of
    the logarithms: this gives bounds for the integrand on the cell and
    on the Bernstein ellipses around it which do not suffer from the
    cancellation between the terms of phi (see el_bound_est below).  A
    cell on which the bound, times the length, is below the tolerance
    is not integrated.  On the others, the Gauss-Legendre rule with n
    nodes has an error of at most (64 / 15) M h / ((rho - 1) rho^(2 n - 1))
    with the bound M on the ellipse with the parameter
    rho = X + sqrt(X^2 - 1); the number of nodes follows from this,
    minimised over X, and a cell is bisected if its two halves take
    fewer nodes together (or if no rule will do).  The subdivision
    starts from the part of the interval around each peak of the
    integrand outside which the integrand is negligible.

    At the end 0 where a term j vanishes (al_j = 0), the cell [0, l] is
    bounded by l (be_j l)^mu_j sup exp(phi without that term) if this is
    negligible, and otherwise integrated term by term with the power
    series of the other factor (el_end_series).  The tail beyond T on
    an infinite interval is at most exp(phi(T)) / kappa if
    phi' <= -kappa < 0 there.

    The integrand is evaluated as exp(phi(t) - c) with c close to the
    maximum of phi, with the precision raised by the number of bits of
    the terms of phi (which cancel), in vectorized dfloat arithmetic
    with N <= 4 doubles or in arb arithmetic, where the precision is
    lowered at the nodes where the integrand is small.

    What is rigorous and what is not.  The cells, the rules, the
    reference value c and the tolerances are chosen from double
    precision estimates; every bound which enters the radius of the
    result is computed in arb and mag arithmetic for the choices made,
    with the input balls as parameters, so that the result is an
    enclosure whatever the choices are.
*/

#include <math.h>
#include <stdlib.h>
#include "arb.h"
#include "arb_hypgeom.h"
#include "arb_hypgeom/impl.h"
#include "acb_calc/impl.h"
#include "dfloat.h"

#define EL_MAXTERMS 3
#define EL_MAXPEAKS 8
#define EL_LN2 0.69314718055994530942
/* a cell is integrated if its half-length is at most this fraction of
   the distance of its midpoint to the nearest singularity */
#define EL_THETA 0.75

typedef struct
{
    int m;
    arb_srcptr lam;
    arb_srcptr mu[EL_MAXTERMS], al[EL_MAXTERMS], be[EL_MAXTERMS];
    double dlam, dmu[EL_MAXTERMS], dal[EL_MAXTERMS], dbe[EL_MAXTERMS];
    double B;                /* the interval is (0, B); B can be infinite */
    int sing0;               /* the term vanishing at 0, or -1 */
    int unit[EL_MAXTERMS];   /* be_j = 1 or -1 (or 0 for another number) */
}
el_kernel;

#define N 1
#define DB(name) d1b_##name
#define DBV(name) _d1b_##name
#define ELQ(name) elq_##name##_d1
#include "exp_log_quad.inc"
#undef N
#undef DB
#undef DBV
#undef ELQ

#define N 2
#define DB(name) d2b_##name
#define DBV(name) _d2b_##name
#define ELQ(name) elq_##name##_d2
#include "exp_log_quad.inc"
#undef N
#undef DB
#undef DBV
#undef ELQ

#define N 3
#define DB(name) d3b_##name
#define DBV(name) _d3b_##name
#define ELQ(name) elq_##name##_d3
#include "exp_log_quad.inc"
#undef N
#undef DB
#undef DBV
#undef ELQ

#define N 4
#define DB(name) d4b_##name
#define DBV(name) _d4b_##name
#define ELQ(name) elq_##name##_d4
#include "exp_log_quad.inc"
#undef N
#undef DB
#undef DBV
#undef ELQ

/* Heuristic: phi, phi', phi'' at t in double precision, without the
   term skip */
static void
el_phi(double * phi, double * d1, double * d2, const el_kernel * K, double t, int skip)
{
    double p = K->dlam * t, q = K->dlam, r = 0.0, u, v;
    int j;

    for (j = 0; j < K->m; j++)
    {
        if (j == skip || K->dmu[j] == 0.0)
            continue;
        u = K->dal[j] + K->dbe[j] * t;
        v = K->dbe[j] / u;
        p += K->dmu[j] * log(u);
        q += K->dmu[j] * v;
        r -= K->dmu[j] * v * v;
    }

    *phi = p; *d1 = q; *d2 = r;
}

/* the distance from t to the nearest singularity */
static double
el_dist(const el_kernel * K, double t)
{
    double d = INFINITY;
    int j;

    for (j = 0; j < K->m; j++)
        if (K->dmu[j] != 0.0 && K->dbe[j] != 0.0)
            d = fmin(d, fabs(t + K->dal[j] / K->dbe[j]));
    return d;
}

/* -log(1 - x) - x */
static double
el_lrem(double x)
{
    if (!(x < 1.0))
        return INFINITY;
    if (x < 0.01)
        return 0.5 * x * x / (1.0 - x);
    return -log1p(-x) - x;
}

/* -(lambda + sum_{mu_j be_j > 0} mu_j be_j / (al_j + be_j t)): phi' <= -kappa
   beyond t if this is positive, since the positive terms decrease and
   the others are negative (be_j > 0) */
static double
el_kappa(const el_kernel * K, double t)
{
    double s = -K->dlam, v;
    int j;

    for (j = 0; j < K->m; j++)
    {
        v = K->dmu[j] * K->dbe[j] / (K->dal[j] + K->dbe[j] * t);
        if (v > 0.0)
            s -= v;
    }
    return s;
}

/* The bounds.  With d_j = m + al_j / be_j (the distance from m to the
   singularity of the term j) and zeta = u / d_j,

       phi(m + u) - phi(m) = phi'(m) u + sum_j mu_j r(u / d_j),
       r(zeta) = log(1 + zeta) - zeta = -zeta^2 / 2 + zeta^3 / 3 - ...

   On the ellipse with the semi-axes R = h X (real) and h sqrt(X^2 - 1)
   around m (the Bernstein ellipse with the parameter
   rho = X + sqrt(X^2 - 1) of the cell [m - h, m + h]), with
   A = R / |d_j| < 1 and q = 1 - 1 / X^2: |zeta| <= A and
   Re(-zeta^2 / 2) <= q A^2 / 2 (at the top of the ellipse), so that

       Re r(zeta) <= min(L(A), q A^2 / 2 + A^3 / (3 (1 - A))),
       -Re r(zeta) <= L(A),      L(A) = -log(1 - A) - A = sum_{k >= 2} A^k / k,

   the first for the terms with mu_j >= 0 and the second for the others
   (the bound L(A) holds for both), and Re(phi'(m) u) <= |phi'(m)| R.
   On the cell itself (zeta real, |zeta| <= A <= 1; q < 0 below), r <= 0,
   and the terms with mu_j >= 0 do not contribute at all. */

/* an estimate of the bound for sup Re(phi - c) on the ellipse (the
   cell itself if q < 0) */
static double
el_bound_est(const el_kernel * K, double m, double R, double q, double cd)
{
    double phi, d1, d2, s, A, v;
    int j;

    el_phi(&phi, &d1, &d2, K, m, -1);
    s = phi - cd + fabs(d1) * R;
    for (j = 0; j < K->m; j++)
    {
        if (K->dmu[j] != 0.0 && K->dbe[j] != 0.0)
        {
            A = R / fabs(m + K->dal[j] / K->dbe[j]);
            v = el_lrem(A);
            if (K->dmu[j] > 0.0)
            {
                if (q < 0.0)
                    v = (A <= 1.0) ? 0.0 : INFINITY;
                else if (A < 1.0)
                    v = fmin(v, 0.5 * q * A * A + A * A * A / (3.0 * (1.0 - A)));
            }
            s += fabs(K->dmu[j]) * v;
        }
    }
    return s;
}

/* Rigorous: phi and phi' at m without the term skip */
static void
el_phi_arb(arb_t phi, arb_t d1, const el_kernel * K, const arb_t m, int skip, slong prec)
{
    arb_t u, v;
    int j;

    arb_init(u);
    arb_init(v);

    arb_mul(phi, K->lam, m, prec);
    arb_set(d1, K->lam);

    for (j = 0; j < K->m; j++)
    {
        if (j == skip || arb_is_zero(K->mu[j]))
            continue;
        arb_mul(u, K->be[j], m, prec);
        arb_add(u, u, K->al[j], prec);
        arb_div(v, K->be[j], u, prec);
        arb_addmul(d1, K->mu[j], v, prec);
        arb_log(u, u, prec);
        arb_addmul(phi, K->mu[j], u, prec);
    }

    arb_clear(u);
    arb_clear(v);
}

/* Rigorous: x >= R / |d_j| (infinite if d_j may be zero) */
static void
el_ratio(mag_t x, const el_kernel * K, int j, const arb_t m, const mag_t R, slong prec)
{
    arb_t u;
    mag_t d;

    /* (a constant term) */
    if (arb_is_zero(K->be[j]))
    {
        mag_zero(x);
        return;
    }

    arb_init(u);
    mag_init(d);

    /* |al_j + be_j m| / |be_j| */
    arb_mul(u, K->be[j], m, prec);
    arb_add(u, u, K->al[j], prec);
    arb_div(u, u, K->be[j], prec);
    arb_get_mag_lower(d, u);
    if (mag_is_zero(d))
        mag_inf(x);
    else
        mag_div(x, R, d);

    arb_clear(u);
    mag_clear(d);
}

/* Rigorous: res >= -log(1 - x) - x = sum_{k >= 2} x^k / k, which is at
   most x^2 / (2 (1 - x)) */
static void
el_lrem_mag(mag_t res, const mag_t x)
{
    if (mag_cmp_2exp_si(x, 0) >= 0)
    {
        mag_inf(res);
    }
    else if (mag_cmp_2exp_si(x, -2) < 0)
    {
        mag_t t;
        mag_init(t);
        mag_one(t);
        mag_sub_lower(t, t, x);
        mag_mul(res, x, x);
        mag_div(res, res, t);
        mag_mul_2exp_si(res, res, -1);
        mag_clear(t);
    }
    else
    {
        /* (increasing in x: at the upper bound) */
        arb_t t, u;
        arb_init(t);
        arb_init(u);
        arf_set_mag(arb_midref(u), x);
        arb_neg(t, u);
        arb_log1p(t, t, 64);
        arb_neg(t, t);
        arb_sub(t, t, u, 64);
        arb_get_mag(res, t);
        arb_clear(t);
        arb_clear(u);
    }
}

/* Rigorous: M >= sup |exp(phi(t) - c)| over the ellipse around m with
   the real semi-axis R and the parameter q (over [m - R, m + R] if
   q < 0), without the term skip */
static void
el_bound(mag_t M, const el_kernel * K, const arb_t c, const arb_t m, const mag_t R,
    double q, int skip, slong prec)
{
    arb_t phi, d1;
    mag_t t, x, y, z;
    int j;

    arb_init(phi);
    arb_init(d1);
    mag_init(t);
    mag_init(x);
    mag_init(y);
    mag_init(z);

    el_phi_arb(phi, d1, K, m, skip, prec);
    arb_sub(phi, phi, c, prec);

    arb_get_mag(t, d1);
    mag_mul(t, t, R);

    for (j = 0; j < K->m; j++)
    {
        if (j == skip || arb_is_zero(K->mu[j]))
            continue;

        /* x = A */
        el_ratio(x, K, j, m, R, prec);

        if (arb_is_nonnegative(K->mu[j]) && q < 0.0)
        {
            if (mag_cmp_2exp_si(x, 0) > 0)
                mag_inf(t);
            continue;
        }

        el_lrem_mag(y, x);

        if (arb_is_nonnegative(K->mu[j]) && mag_cmp_2exp_si(x, 0) < 0)
        {
            /* q A^2 / 2 + A^3 / (3 (1 - A)) */
            mag_one(z);
            mag_sub_lower(z, z, x);
            mag_mul_ui_lower(z, z, 3);
            mag_div(z, x, z);
            mag_set_d(M, 0.5 * q);
            mag_add(z, z, M);
            mag_mul(z, z, x);
            mag_mul(z, z, x);
            mag_min(y, y, z);
        }

        arb_get_mag(M, K->mu[j]);
        mag_addmul(t, y, M);
    }

    arb_add_error_mag(phi, t);
    arb_exp(phi, phi, prec);
    arb_get_mag(M, phi);
    if (!arb_is_finite(phi))
        mag_inf(M);

    arb_clear(phi);
    arb_clear(d1);
    mag_clear(t);
    mag_clear(x);
    mag_clear(y);
    mag_clear(z);
}

typedef struct
{
    const el_kernel * K;
    arb_t c;                 /* the reference value (exact) */
    double cd;
    slong prec;
    slong bprec;             /* for the bounds */
    double logtol;           /* log of the tolerance for a cell, relative
                                to exp(c) */
    double magd;             /* the size of the terms of phi */
    double peak[EL_MAXPEAKS];
    int npeaks;
    int imax;                /* the largest rule */
    slong budget;            /* of nodes */
    arb_t sum;
    mag_t err;
    void * scratch;
}
el_ctx;

/* The rule gi on the cell [m - h, m + h] (m and h exact), where the
   integrand is about exp(lmax) at most.  Returns 0 on failure. */
static int
el_cell_value(arb_t res, el_ctx * ctx, const arb_t m, const arb_t h, double md, double hd,
    int gi, double lmax)
{
    const el_kernel * K = ctx->K;
    double bits, lmag = log2(ctx->magd);
    slong n = _acb_calc_gl_steps[gi], k, wp;
    int j, N, ok = 1;

    /* the bits of the arithmetic: of the result, and of the terms of
       the exponent (which cancel); fewer where the integrand is small */
    bits = ctx->prec + 10.0 + lmag;
    bits = fmax(bits + fmin(lmax / EL_LN2 + 2.0, 0.0), lmag + 14.0);
    /* (not where the integrand relative to exp(c) is too small for a
       double: at a high precision, a cell far below the peak) */
    N = (dfloat_is_supported() && bits <= 4 * 53 && gi < ACB_CALC_GL_DFLOAT_STEPS
            && lmax / EL_LN2 > -900.0)
        ? (int) ceil(fmax(bits / 53.0, 1.0)) : 0;

    if (N != 0)
    {
        const void * xs, * ws;

        _acb_calc_gl_nodes_dfloat(&xs, &ws, N, gi);
        if (N == 1) ok = elq_direct_d1(res, K, ctx->c, m, h, xs, ws, n, ctx->scratch);
        else if (N == 2) ok = elq_direct_d2(res, K, ctx->c, m, h, xs, ws, n, ctx->scratch);
        else if (N == 3) ok = elq_direct_d3(res, K, ctx->c, m, h, xs, ws, n, ctx->scratch);
        else ok = elq_direct_d4(res, K, ctx->c, m, h, xs, ws, n, ctx->scratch);
    }
    else
    {
        /* With arb arithmetic, at a precision which is lowered per
           node: the exponent as the dot product
           lambda t + sum_j mu_j log(al_j + be_j t). */
        arb_srcptr xs, ws;
        arb_struct xv[EL_MAXTERMS + 1], yv[EL_MAXTERMS + 1];
        arb_t t, v, u[EL_MAXTERMS];
        int idx[EL_MAXTERMS], nx = 0;
        slong kk, p;

        arb_init(t);
        arb_init(v);
        for (j = 0; j < K->m; j++)
        {
            arb_init(u[j]);
            if (!arb_is_zero(K->mu[j]))
                idx[nx++] = j;
        }

        xv[0] = *K->lam;
        for (j = 0; j < nx; j++)
        {
            xv[j + 1] = *K->mu[idx[j]];
            yv[j + 1] = *u[j];
        }

        wp = (slong) (ctx->prec + 10.0 + lmag) + 8;
        _acb_calc_gl_nodes(&xs, &ws, gi, wp);
        arb_zero(res);

        for (k = 0; k < n; k++)
        {
            double xd, pd, d1, d2;

            kk = (2 * k < n) ? k : n - 1 - k;
            xd = arf_get_d(arb_midref(xs + kk), ARF_RND_NEAR);
            if (2 * k >= n)
                xd = -xd;
            el_phi(&pd, &d1, &d2, K, md + hd * xd, -1);
            p = wp + (slong) fmax(fmin((pd - ctx->cd) / EL_LN2 + 2.0, 0.0), -1e9);
            p = FLINT_MAX(p, (slong) (lmag + 14.0));

            /* t = m +- h x */
            arb_mul(t, xs + kk, h, p);
            if (2 * k >= n)
                arb_sub(t, m, t, p);
            else
                arb_add(t, m, t, p);

            for (j = 0; j < nx; j++)
            {
                int i = idx[j];

                /* al_i + be_i t */
                if (K->unit[i] == 1)
                    arb_add(u[j], K->al[i], t, p);
                else if (K->unit[i] == -1)
                    arb_sub(u[j], K->al[i], t, p);
                else
                {
                    arb_mul(u[j], K->be[i], t, p);
                    arb_add(u[j], u[j], K->al[i], p);
                }
                arb_log(u[j], u[j], p);
                yv[j + 1] = *u[j];
            }

            yv[0] = *t;
            arb_dot(v, ctx->c, 1, xv, 1, yv, 1, nx + 1, p);
            arb_neg(v, v);
            arb_exp(v, v, p);
            arb_addmul(res, v, ws + kk, wp);
        }
        arb_mul(res, res, h, wp);

        arb_clear(t);
        arb_clear(v);
        for (j = 0; j < K->m; j++)
            arb_clear(u[j]);
    }

    return ok && arb_is_finite(res);
}

/* The cell [0, l] at the end where the term s vanishes (al_s = 0,
   be_s > 0), by the power series of the other factor: with

       g(t) = exp(lambda t + sum_{j != s} mu_j log(al_j + be_j t))
            = sum_k g_k t^k,

   int_0^l (be_s t)^mu_s g(t) dt = be_s^mu_s sum_k g_k l^(k + mu_s + 1) / (k + mu_s + 1).

   Since D g' = N g with the polynomials D = prod_{j != s} (al_j + be_j t)
   and N = lambda D + sum_j mu_j be_j prod_{i != j, s} (al_i + be_i t), the
   coefficients satisfy

       D_0 (k + 1) g_{k+1} = sum_i N_i g_{k-i} - sum_{i >= 1} D_i (k - i + 1) g_{k-i+1}.

   By Cauchy's inequality |g_k| <= M r^-k with M >= sup |g| on |t| <= r,

       M = exp(|lambda| r) prod_j (al_j +- |be_j| r)^mu_j

   (the sign of mu_j; r below the distance to the other singularities),
   so that with r = 4 l the terms from K on add up to at most
   M 4^-K (4 / 3) / (K + mu_s + 1), times the same factor.  Returns 0 if
   this does not apply. */
static int
el_end_series(el_ctx * ctx, const arb_t l, double ld)
{
    const el_kernel * K = ctx->K;
    int s = K->sing0, idx[EL_MAXTERMS], no = 0, i, j, ok = 1;
    slong k, nterms, wp;
    double logM, logg, pref, r = 4.0 * ld;
    arb_t D[3], Nn[3], t, u, sum, pw, lk, M;
    arb_ptr g;
    mag_t tail;

    for (j = 0; j < K->m; j++)
    {
        if (j == s || K->dmu[j] == 0.0)
            continue;
        /* (the other terms: analytic and positive on the disk) */
        if (!(K->dal[j] - fabs(K->dbe[j]) * r > 0.5 * K->dal[j]) || !arb_is_positive(K->al[j]))
            return 0;
        idx[no++] = j;
    }
    if (no > 2 || !arb_is_zero(K->al[s]) || !arb_is_positive(K->be[s]))
        return 0;

    /* the number of terms, from estimates (and logg: of g(0)) */
    logM = fabs(K->dlam) * r;
    logg = 0.0;
    for (i = 0; i < no; i++)
    {
        j = idx[i];
        logM += K->dmu[j] * log(K->dal[j] + ((K->dmu[j] > 0.0) ? 1.0 : -1.0) * fabs(K->dbe[j]) * r);
        logg += K->dmu[j] * log(K->dal[j]);
    }
    pref = (K->dmu[s] + 1.0) * log(ld) + K->dmu[s] * log(K->dbe[s]) - ctx->cd;
    nterms = (slong) ((logM + pref - ctx->logtol) / log(4.0)) + 3;
    if (nterms < 2)
        nterms = 2;
    /* (the terms can be larger than the sum by the growth of g on the
       disk) */
    if (nterms > 20000 || !(logM - logg < 60.0) || !(logM - logg >= 0.0))
        return 0;

    wp = ctx->prec + 30 + 2 * FLINT_BIT_COUNT(nterms) + (slong) ((logM - logg) / EL_LN2);

    for (i = 0; i < 3; i++)
    {
        arb_init(D[i]);
        arb_init(Nn[i]);
    }
    arb_init(t); arb_init(u); arb_init(sum); arb_init(pw); arb_init(lk); arb_init(M);
    mag_init(tail);
    g = _arb_vec_init(nterms + 1);

    /* D and N */
    arb_one(D[0]);
    for (i = 0; i < no; i++)
    {
        j = idx[i];
        /* D *= (al_j + be_j t) */
        arb_mul(D[2], D[1], K->be[j], wp);
        arb_mul(D[1], D[1], K->al[j], wp);
        arb_addmul(D[1], D[0], K->be[j], wp);
        arb_mul(D[0], D[0], K->al[j], wp);
    }
    for (i = 0; i < 3; i++)
        arb_mul(Nn[i], D[i], K->lam, wp);
    if (no == 1)
    {
        arb_addmul(Nn[0], K->mu[idx[0]], K->be[idx[0]], wp);
    }
    else if (no == 2)
    {
        for (i = 0; i < 2; i++)
        {
            int j1 = idx[i], j2 = idx[1 - i];
            /* mu_j1 be_j1 (al_j2 + be_j2 t) */
            arb_mul(t, K->mu[j1], K->be[j1], wp);
            arb_addmul(Nn[0], t, K->al[j2], wp);
            arb_addmul(Nn[1], t, K->be[j2], wp);
        }
    }

    /* g_0 = prod al_j^mu_j, relative to exp(c) */
    arb_neg(t, ctx->c);
    for (i = 0; i < no; i++)
    {
        arb_log(u, K->al[idx[i]], wp);
        arb_addmul(t, u, K->mu[idx[i]], wp);
    }
    arb_exp(g, t, wp);

    for (k = 0; k < nterms; k++)
    {
        arb_zero(t);
        for (i = 0; i <= 2 && i <= k; i++)
            arb_addmul(t, Nn[i], g + k - i, wp);
        for (i = 1; i <= 2 && i <= k; i++)
        {
            arb_mul_ui(u, g + k - i + 1, k - i + 1, wp);
            arb_submul(t, D[i], u, wp);
        }
        arb_div(t, t, D[0], wp);
        arb_div_ui(g + k + 1, t, k + 1, wp);
    }

    /* sum_k g_k l^k / (k + mu_s + 1) */
    arb_one(lk);
    for (k = 0; k < nterms; k++)
    {
        arb_add_ui(t, K->mu[s], k + 1, wp);
        arb_div(t, g + k, t, wp);
        arb_addmul(sum, t, lk, wp);
        arb_mul(lk, lk, l, wp);
    }

    /* M, relative to exp(c) */
    arb_mul_2exp_si(u, l, 2);
    arb_abs(t, K->lam);
    arb_mul(M, t, u, wp);
    arb_sub(M, M, ctx->c, wp);
    for (i = 0; i < no; i++)
    {
        j = idx[i];
        /* al_j -+ |be_j| r: positive (the term is analytic on the
           disk), and the larger or the smaller according to the sign
           of mu_j */
        arb_abs(t, K->be[j]);
        arb_mul(t, t, u, wp);
        arb_sub(pw, K->al[j], t, wp);
        if (!arb_is_positive(pw))
            ok = 0;
        if (arb_is_nonnegative(K->mu[j]))
            arb_add(t, K->al[j], t, wp);
        else if (arb_is_negative(K->mu[j]))
            arb_set(t, pw);
        else
            ok = 0;
        arb_log(t, t, wp);
        arb_addmul(M, t, K->mu[j], wp);
    }
    arb_exp(M, M, wp);
    /* 4^-K (4 / 3) / (K + mu_s + 1) <= 4^-K (4 / 3) / K */
    arb_mul_2exp_si(M, M, -2 * nterms + 2);
    arb_div_ui(M, M, 3 * nterms, wp);
    arb_get_mag(tail, M);
    arb_add_error_mag(sum, tail);

    /* times be_s^mu_s l^(mu_s + 1) */
    arb_add_ui(t, K->mu[s], 1, wp);
    arb_pow(pw, l, t, wp);
    arb_mul(sum, sum, pw, wp);
    arb_pow(pw, K->be[s], K->mu[s], wp);
    arb_mul(sum, sum, pw, wp);

    ok = ok && arb_is_finite(sum);
    if (ok)
        arb_add(ctx->sum, ctx->sum, sum, ctx->prec + 30);

    for (i = 0; i < 3; i++)
    {
        arb_clear(D[i]);
        arb_clear(Nn[i]);
    }
    arb_clear(t); arb_clear(u); arb_clear(sum); arb_clear(pw); arb_clear(lk); arb_clear(M);
    mag_clear(tail);
    _arb_vec_clear(g, nterms + 1);
    return ok;
}

/* The number of nodes for the cell [m - h, m + h] (estimates): of the
   smallest rule whose error bound, minimised over the ellipse parameter
   X, is below the tolerance; 0 for a cell which is negligible (if this
   is allowed), and infinity if there is no rule. */
static double
el_cell_plan(int * rule, double * bestX, const el_ctx * ctx, double md, double hd, int negl)
{
    static const double Xs[] = { 1.25, 1.5, 2.0, 3.0, 4.0, 6.0, 8.0, 12.0, 24.0,
        1e2, 1e3, 1e4, 1e6, 1e8 };
    const el_kernel * K = ctx->K;
    double dc = el_dist(K, md), bestn = INFINITY;
    int ix, gi;

    if (negl && hd < dc
        && log(2.0 * hd) + el_bound_est(K, md, hd, -1.0, ctx->cd) <= ctx->logtol - 1.0)
        return 0.0;

    if (!(hd <= EL_THETA * dc))
        return INFINITY;

    for (ix = 0; ix < (int) (sizeof(Xs) / sizeof(double)); ix++)
    {
        double X = Xs[ix], rho = X + sqrt(X * X - 1.0), Rd = hd * X, need;

        if (!(Rd < 0.98 * dc))
            break;
        need = log(64.0 / 15.0) + el_bound_est(K, md, Rd, 1.0 - 1.0 / (X * X), ctx->cd) + log(hd)
            - log(rho - 1.0) - ctx->logtol + 0.5;
        need = 0.5 * (need / log(rho) + 1.0);
        if (need < bestn)
        {
            bestn = need;
            *bestX = X;
        }
    }

    for (gi = 0; gi <= ctx->imax; gi++)
        if (_acb_calc_gl_steps[gi] >= bestn)
            break;

    if (gi > ctx->imax)
        return INFINITY;

    *rule = gi;
    return (double) _acb_calc_gl_steps[gi];
}

/* The cell [u, v] (exact numbers): bounded, integrated or bisected.
   Returns 0 on failure. */
static int
el_cell(el_ctx * ctx, const arf_t u, const arf_t v, int depth)
{
    const el_kernel * K = ctx->K;
    double ud = arf_get_d(u, ARF_RND_NEAR), vd = arf_get_d(v, ARF_RND_NEAR);
    double md = 0.5 * (ud + vd), hd = 0.5 * (vd - ud), phi, d1, d2, sup, dc;
    int s0 = (K->sing0 >= 0 && arf_is_zero(u) && !arb_is_zero(K->mu[K->sing0]));
    int i, ok = 1, done = 0, negl = 1;
    arb_t m, h;
    mag_t R, M;

    /* (also: cells which doubles do not resolve) */
    ctx->budget -= 10;   /* (the work, in nodes) */
    if (depth > 4000 || ctx->budget < 0 || !(hd > 1e-13 * vd))
        return 0;

    arb_init(m);
    arb_init(h);
    mag_init(R);
    mag_init(M);

    /* (exactly) */
    arf_add(arb_midref(m), u, v, ARF_PREC_EXACT, ARF_RND_DOWN);
    arf_mul_2exp_si(arb_midref(m), arb_midref(m), -1);
    arf_sub(arb_midref(h), v, u, ARF_PREC_EXACT, ARF_RND_DOWN);
    arf_mul_2exp_si(arb_midref(h), arb_midref(h), -1);
    arf_get_mag(R, arb_midref(h));

    if (s0)
    {
        /* The cell at the end where the term j vanishes: its factor is
           at most (|be_j| v)^mu_j there. */
        int j = K->sing0;
        double len = 2.0 * hd, est;

        sup = -INFINITY;
        for (i = 0; i < 3; i++)
        {
            el_phi(&phi, &d1, &d2, K, (i == 0) ? ud : (i == 1) ? md : vd, j);
            sup = fmax(sup, phi);
        }
        est = log(len) + K->dmu[j] * log(fabs(K->dbe[j]) * len) + sup - ctx->cd;

        if (est <= ctx->logtol - 2.0)
        {
            arb_t f;
            mag_t t;

            arb_init(f);
            mag_init(t);

            el_bound(M, K, ctx->c, m, R, -1.0, j, ctx->bprec);
            arb_mul_2exp_si(f, h, 1);
            arb_mul(f, f, K->be[j], ctx->bprec);
            arb_abs(f, f);
            arb_pow(f, f, K->mu[j], ctx->bprec);
            arb_get_mag(t, f);
            mag_mul(M, M, t);
            mag_mul(M, M, R);
            mag_mul_2exp_si(M, M, 1);
            {
                arb_t e;
                arb_init(e);
                arb_set_d(e, ctx->logtol);
                arb_exp(e, e, 64);
                arb_get_mag_lower(t, e);
                arb_clear(e);
            }
            if (arb_is_finite(f) && mag_cmp(M, t) <= 0)
            {
                mag_add(ctx->err, ctx->err, M);
                done = 1;
            }

            arb_clear(f);
            mag_clear(t);
        }

        /* not negligible: by the power series at the end, if the cell
           is short enough for it */
        if (!done)
        {
            double q = fabs(K->dlam);
            int i2;

            for (i2 = 0; i2 < K->m; i2++)
                if (i2 != j && K->dmu[i2] != 0.0)
                    q += fabs(K->dmu[i2] * K->dbe[i2] / K->dal[i2]);

            if (q * len <= 2.0)
            {
                arb_t l;
                arb_init(l);
                arf_set(arb_midref(l), v);
                done = el_end_series(ctx, l, vd);
                arb_clear(l);
            }
        }
    }
    else
    {
        sup = -INFINITY;
        for (i = 0; i < 3; i++)
        {
            el_phi(&phi, &d1, &d2, K, (i == 0) ? ud : (i == 1) ? md : vd, -1);
            sup = fmax(sup, phi);
        }
        for (i = 0; i < ctx->npeaks; i++)
        {
            if (ctx->peak[i] >= ud && ctx->peak[i] <= vd)
            {
                el_phi(&phi, &d1, &d2, K, ctx->peak[i], -1);
                sup = fmax(sup, phi);
            }
        }
        sup -= ctx->cd;

        /* negligible: bounded on the cell */
        if (el_cell_plan(&i, &dc, ctx, md, hd, 1) == 0.0)
        {
            negl = 0;
            mag_t t;
            arb_t e;

            mag_init(t);
            arb_init(e);
            el_bound(M, K, ctx->c, m, R, -1.0, -1, ctx->bprec);
            mag_mul(M, M, R);
            mag_mul_2exp_si(M, M, 1);
            arb_set_d(e, ctx->logtol);
            arb_exp(e, e, 64);
            arb_get_mag_lower(t, e);
            if (mag_cmp(M, t) <= 0)
            {
                mag_add(ctx->err, ctx->err, M);
                done = 1;
            }
            mag_clear(t);
            arb_clear(e);
        }

        if (!done)
        {
            double bestX = 0.0, cw, cl, cr, X2;
            int gi = 0, g2;

            /* the cell as a whole if this does not take more nodes than
               its two halves (which are not taken for negligible if the
               bound for the cell has just failed where the estimate did
               not: wide parameters) */
            cw = el_cell_plan(&gi, &bestX, ctx, md, hd, 0);
            if (cw < INFINITY)
            {
                cl = el_cell_plan(&g2, &X2, ctx, md - 0.5 * hd, 0.5 * hd, negl);
                cr = el_cell_plan(&g2, &X2, ctx, md + 0.5 * hd, 0.5 * hd, negl);
            }
            else
                cl = cr = 0.0;

            if (cw <= cl + cr)
            {
                arb_t val;
                mag_t e, RX;

                arb_init(val);
                mag_init(e);
                mag_init(RX);

                ctx->budget -= _acb_calc_gl_steps[gi];
                ok = el_cell_value(val, ctx, m, h, md, hd, gi, sup);

                if (ok)
                {
                    mag_set_d(RX, bestX);
                    mag_mul(RX, RX, R);
                    el_bound(M, K, ctx->c, m, RX, 1.0 - 1.0 / (bestX * bestX), -1, ctx->bprec);
                    _acb_calc_gl_error_bound(e, M, R, bestX, _acb_calc_gl_steps[gi]);
                    mag_add(ctx->err, ctx->err, e);
                    arb_add(ctx->sum, ctx->sum, val, ctx->prec + 30);
                    ok = mag_is_finite(e);
                }
                done = 1;

                arb_clear(val);
                mag_clear(e);
                mag_clear(RX);
            }
        }
    }

    if (!done)
    {
        ok = el_cell(ctx, u, arb_midref(m), depth + 1)
            && el_cell(ctx, arb_midref(m), v, depth + 1);
    }

    arb_clear(m);
    arb_clear(h);
    mag_clear(R);
    mag_clear(M);
    return ok;
}

static int
el_cmp_d(const void * a, const void * b)
{
    double x = *(const double *) a, y = *(const double *) b;
    return (x > y) - (x < y);
}

/* res = int_0^B exp(lambda t + sum_j mu_j log(al_j + be_j t)) dt, where
   sing0 is the index of the term which vanishes at 0 (or -1), and B is
   a double or infinite (no term may vanish at B).  The tolerance is
   2^-prec relative to the integral, or to exp(logscale) if that is
   larger (for an integral which is a part of a sum; -inf otherwise).
   Returns 0 (without a result) if the integral is not computed:
   parameters which are not finite or too large for the double
   precision planning, mu_j which is not nonnegative at 0, or an
   integrand whose peak is not found. */
int
_arb_hypgeom_integrate_exp_log(arb_t res, const arb_t lambda, arb_srcptr mu,
    arb_srcptr al, arb_srcptr be, int m, int sing0, double B, double logscale, slong prec)
{
    int infinite = !(B < INFINITY);
    el_kernel K;
    el_ctx ctx;
    double grid[512], pts[1024], phi, d1, d2, best = -INFINITY, tbest = 0.0;
    double logI = -INFINITY, scale, T = 1.0;
    double pk[EL_MAXPEAKS], pw[EL_MAXPEAKS];
    int j, k, ng = 0, np = 0, npk = 0, ok = 1;
    arf_t u, v;

    if (m > EL_MAXTERMS || !arb_is_finite(lambda)
        || arf_cmpabs_2exp_si(arb_midref(lambda), 40) >= 0)
        return 0;

    K.m = m;
    K.lam = lambda;
    K.dlam = arf_get_d(arb_midref(lambda), ARF_RND_NEAR);
    K.B = B;
    K.sing0 = sing0;
    scale = 1.0;
    for (j = 0; j < m; j++)
    {
        if (!arb_is_finite(mu + j) || !arb_is_finite(al + j) || !arb_is_finite(be + j)
            || arf_cmpabs_2exp_si(arb_midref(mu + j), 40) >= 0
            || arf_cmpabs_2exp_si(arb_midref(al + j), 40) >= 0
            || arf_cmpabs_2exp_si(arb_midref(be + j), 40) >= 0)
            return 0;
        K.mu[j] = mu + j; K.al[j] = al + j; K.be[j] = be + j;
        K.dmu[j] = arf_get_d(arb_midref(mu + j), ARF_RND_NEAR);
        K.dal[j] = arf_get_d(arb_midref(al + j), ARF_RND_NEAR);
        K.dbe[j] = arf_get_d(arb_midref(be + j), ARF_RND_NEAR);
        K.unit[j] = arb_is_one(be + j) ? 1 : arb_equal_si(be + j, -1) ? -1 : 0;
        scale += fabs(K.dmu[j]);
        if (j == sing0 && (!arb_is_nonnegative(mu + j) || !arb_is_zero(al + j)))
            return 0;
    }

    if (infinite && !(K.dlam < 0.0))
        return 0;

    /* The peaks: the local maxima of phi, located on a grid which is
       geometric towards the ends, and refined by bisection */
    if (!infinite)
    {
        for (k = 6; k <= 70; k++)
            grid[ng++] = ldexp(B, -k);
        for (k = 1; k <= 31; k++)
            grid[ng++] = B * k / 32.0;
    }
    else
    {
        scale /= -K.dlam;
        for (k = -60; k <= 60; k++)
            grid[ng++] = ldexp(scale, k);
    }
    qsort(grid, ng, sizeof(double), el_cmp_d);

    for (k = 0; k <= ng; k++)
    {
        double a = (k == 0) ? 0.0 : grid[k - 1], b = (k == ng) ? K.B : grid[k];
        double da, db, t0;
        int it;

        /* the sign of phi' at the two ends (at an end of the interval:
           +-inf if a term vanishes there) */
        if (k == 0 && sing0 >= 0 && K.dmu[sing0] > 0.0)
            da = 1.0;
        else
            el_phi(&phi, &da, &d2, &K, a, -1);
        if (k == ng && infinite)
            db = -1.0;
        else
            el_phi(&phi, &db, &d2, &K, b, -1);

        if (k == 0 && da <= 0.0)
        {
            /* (decreasing from the left end) */
            t0 = 0.0;
        }
        else if (k == ng && db > 0.0)
        {
            t0 = K.B;
        }
        else if (da > 0.0 && db <= 0.0 && b < INFINITY)
        {
            for (it = 0; it < 80 && b - a > 1e-14 * b; it++)
            {
                t0 = 0.5 * (a + b);
                el_phi(&phi, &d1, &d2, &K, t0, -1);
                if (d1 > 0.0)
                    a = t0;
                else
                    b = t0;
            }
            t0 = 0.5 * (a + b);
        }
        else
            continue;

        if (npk >= EL_MAXPEAKS)
            return 0;

        el_phi(&phi, &d1, &d2, &K, t0, -1);
        if (!(phi > -1e300) || !(phi < 1e300))
            continue;

        pk[npk] = t0;
        /* the width: of a peak inside, or of the decay from an end */
        if (t0 == 0.0 || t0 == K.B)
            pw[npk] = 1.0 / (fabs(d1) + 1e-300);
        else
            pw[npk] = 1.0 / sqrt(fabs(d2) + 1e-300);
        pw[npk] = fmin(pw[npk], infinite ? scale : B);

        if (phi > best)
        {
            best = phi;
            tbest = t0;
        }
        npk++;
    }

    if (npk == 0 || !(best > -1e300))
        return 0;

    for (k = 0; k < npk; k++)
    {
        el_phi(&phi, &d1, &d2, &K, pk[k], -1);
        phi += log(2.5 * pw[k]) - best;
        logI = (logI == -INFINITY) ? phi : fmax(logI, phi) + log1p(exp(-fabs(logI - phi)));
    }

    /* the context */
    ctx.K = &K;
    ctx.prec = prec;
    ctx.cd = best;
    ctx.npeaks = npk;
    for (k = 0; k < npk; k++)
        ctx.peak[k] = pk[k];
    /* (relative to exp(c)) */
    if (logscale - best > logI)
        logI = logscale - best;
    ctx.logtol = logI - (prec + 5.0) * EL_LN2 - log(32.0);
    /* (rules with up to prec / 2 + 60 nodes) */
    for (ctx.imax = 0; ctx.imax + 1 < ACB_CALC_GL_STEPS
        && _acb_calc_gl_steps[ctx.imax + 1] <= 0.5 * prec + 60.0; ctx.imax++) ;
    ctx.budget = 20000 + 200 * prec;

    /* the size of the terms of phi where the integrand is not negligible */
    ctx.magd = fabs(best) + 2.0;
    for (k = 0; k < npk; k++)
    {
        double t = pk[k], s = fabs(K.dlam * t);
        for (j = 0; j < m; j++)
            if (K.dmu[j] != 0.0)
                s += fabs(K.dmu[j]) * (fabs(log(fabs(K.dal[j] + K.dbe[j] * t) + 1e-300)) + 1.0);
        if (s < 1e300)
            ctx.magd = fmax(ctx.magd, s + fabs(best));
    }
    ctx.bprec = 96 + (slong) log2(ctx.magd);

    /* The points of the initial subdivision: the ends of the part
       around each peak outside which the integrand is negligible
       (located by doubling steps from the peak, then by bisection). */
    pts[np++] = 0.0;
    for (k = 0; k < npk; k++)
    {
        int s;

        for (s = -1; s <= 1; s += 2)
        {
            double lo = pk[k], hi = pk[k], step = 2.0 * pw[k], lim = ctx.logtol - 2.0 + best;
            int found = 0;

            for (j = 0; j < 200; j++)
            {
                hi = pk[k] + s * step;
                if (!(hi > 0.0) || !(hi < K.B))
                    break;
                el_phi(&phi, &d1, &d2, &K, hi, -1);
                if (phi < lim)
                {
                    found = 1;
                    break;
                }
                lo = hi;
                step *= 2.0;
            }

            if (found)
            {
                for (j = 0; j < 6; j++)
                {
                    double mid = 0.5 * (lo + hi);
                    el_phi(&phi, &d1, &d2, &K, mid, -1);
                    if (phi < lim)
                        hi = mid;
                    else
                        lo = mid;
                }
                pts[np++] = hi;
                if (infinite && s == 1)
                    T = fmax(T, hi);
            }
            else if (infinite && s == 1)
                return 0;
        }
    }
    if (infinite)
    {
        /* far enough for the bound for the tail */
        for (j = 0; j < 200; j++)
        {
            double kap = el_kappa(&K, T);
            el_phi(&phi, &d1, &d2, &K, T, -1);
            if (kap > 0.0 && phi - best - log(kap) < ctx.logtol - 3.0)
                break;
            pts[np++] = T;
            T *= 2.0;
        }
        if (j == 200)
            return 0;
        pts[np++] = T;
    }
    else
        pts[np++] = B;

    qsort(pts, np, sizeof(double), el_cmp_d);

    /* (no cells which are tiny compared to their neighbours: a peak at
       a point of the subdivision, up to rounding) */
    {
        double end = pts[np - 1];
        int n2 = 1;

        for (k = 1; k < np; k++)
        {
            if (pts[k] - pts[n2 - 1] > 1e-9 * end && end - pts[k] > 1e-9 * end)
                pts[n2++] = pts[k];
        }
        pts[n2++] = end;
        np = n2;
    }

    arb_init(ctx.c);
    arb_init(ctx.sum);
    mag_init(ctx.err);
    arf_init(u);
    arf_init(v);
    ctx.scratch = flint_malloc(ACB_CALC_GL_DFLOAT_SCRATCH(4));

    /* c = phi at the highest peak, rounded (an exact number) */
    {
        arb_t t, d;
        arb_init(t);
        arb_init(d);
        arb_set_d(t, tbest);
        if (tbest == 0.0 || tbest == K.B)
            arb_set_d(ctx.c, best);
        else
        {
            el_phi_arb(ctx.c, d, &K, t, -1, ctx.bprec);
            arb_get_mid_arb(ctx.c, ctx.c);
        }
        if (!arb_is_finite(ctx.c))
            ok = 0;
        ctx.cd = arf_get_d(arb_midref(ctx.c), ARF_RND_NEAR);
        arb_clear(t);
        arb_clear(d);
    }

    for (k = 0; k + 1 < np && ok; k++)
    {
        if (!(pts[k + 1] > pts[k]))
            continue;
        arf_set_d(u, pts[k]);
        arf_set_d(v, pts[k + 1]);
        ok = el_cell(&ctx, u, v, 0);
    }

    /* the tail: phi' <= lambda + sum_{mu_j be_j > 0} mu_j be_j / (al_j + be_j T)
       = -kappa for t >= T (the terms decrease), and the tail is at most
       exp(phi(T) - c) / kappa */
    if (ok && infinite)
    {
        arb_t t, kap, x, y;
        mag_t e;

        arb_init(t); arb_init(kap); arb_init(x); arb_init(y);
        mag_init(e);

        arb_set_d(t, T);
        arb_neg(kap, lambda);
        for (j = 0; j < m; j++)
        {
            if (arb_is_zero(mu + j))
                continue;
            /* (be_j > 0 on an infinite interval) */
            if (!arb_is_positive(be + j) || !arb_is_nonnegative(al + j))
                ok = 0;
            arb_mul(x, be + j, t, ctx.bprec);
            arb_add(x, x, al + j, ctx.bprec);
            arb_div(x, be + j, x, ctx.bprec);
            arb_mul(x, x, mu + j, ctx.bprec);
            /* max(x, 0) */
            if (!arb_is_negative(x))
            {
                arb_get_ubound_arf(arb_midref(y), x, ctx.bprec);
                mag_zero(arb_radref(y));
                if (arf_sgn(arb_midref(y)) > 0)
                    arb_sub(kap, kap, y, ctx.bprec);
            }
        }

        if (!arb_is_positive(kap))
            ok = 0;
        else
        {
            el_phi_arb(x, y, &K, t, -1, ctx.bprec);
            arb_sub(x, x, ctx.c, ctx.bprec);
            arb_exp(x, x, ctx.bprec);
            arb_div(x, x, kap, ctx.bprec);
            arb_get_mag(e, x);
            mag_add(ctx.err, ctx.err, e);
            ok = ok && arb_is_finite(x);
        }

        arb_clear(t); arb_clear(kap); arb_clear(x); arb_clear(y);
        mag_clear(e);
    }

    if (ok)
    {
        arb_t t;
        arb_init(t);
        arb_add_error_mag(ctx.sum, ctx.err);
        arb_exp(t, ctx.c, prec + 30);
        arb_mul(res, ctx.sum, t, prec + 30);
        ok = arb_is_finite(res);
        arb_clear(t);
    }

    flint_free(ctx.scratch);
    arb_clear(ctx.c);
    arb_clear(ctx.sum);
    mag_clear(ctx.err);
    arf_clear(u);
    arf_clear(v);
    return ok;
}

/* log |x|, approximately (from above, within log 2), for any exponent:
   arf_abs_bound_lt_2exp_si is clamped to the range of an slong, which
   integrals with large parameters exceed on 32-bit systems (with
   parameters of 2^36, 1F1 has an integral of 2^(-1.3 10^11)) */
static double
el_log_abs(const arb_t x)
{
    if (arf_is_zero(arb_midref(x)) || !arb_is_finite(x))
        return -INFINITY;
    return EL_LN2 * fmpz_get_d(ARF_EXPREF(arb_midref(x)));
}

/* an estimate of max phi on (0, B) from a few points */
static double
el_quick_max(const arb_t lambda, arb_srcptr mu, arb_srcptr al, arb_srcptr be, int m, double B)
{
    double best = -INFINITY, t, p;
    int j, k;

    for (k = 1; k <= 32; k++)
    {
        t = B * k / 32.0;
        p = arf_get_d(arb_midref(lambda), ARF_RND_NEAR) * t;
        for (j = 0; j < m; j++)
            if (!arb_is_zero(mu + j))
                p += arf_get_d(arb_midref(mu + j), ARF_RND_NEAR)
                    * log(arf_get_d(arb_midref(al + j), ARF_RND_NEAR)
                        + arf_get_d(arb_midref(be + j), ARF_RND_NEAR) * t);
        if (p > best)
            best = p;
    }
    return best;
}

/* The integral over (0, 1) as the sum of the integrals over (0, 1/2)
   of the kernel (lambda, mu, al, be) in t and of the kernel
   (lambda2, mu2, al2, be2) in 1 - t, the latter times exp(shift); the
   term 0 of each vanishes at 0.  The half with the larger integrand
   is computed first, and sets the tolerance of the other. */
static int
el_two_halves(arb_t res, const arb_t lambda, arb_srcptr mu, arb_srcptr al, arb_srcptr be,
    const arb_t lambda2, arb_srcptr mu2, arb_srcptr al2, arb_srcptr be2, const arb_t shift,
    int m, slong prec)
{
    arb_t r1, r2, e;
    double sh = arb_is_zero(shift) ? 0.0 : arf_get_d(arb_midref(shift), ARF_RND_NEAR);
    int ok, first2;

    arb_init(r1);
    arb_init(r2);
    arb_init(e);

    first2 = el_quick_max(lambda2, mu2, al2, be2, m, 0.5) + sh
        > el_quick_max(lambda, mu, al, be, m, 0.5);

    if (!first2)
    {
        ok = _arb_hypgeom_integrate_exp_log(r1, lambda, mu, al, be, m, 0, 0.5, -INFINITY, prec)
            && _arb_hypgeom_integrate_exp_log(r2, lambda2, mu2, al2, be2, m, 0, 0.5,
                el_log_abs(r1) - sh, prec);
    }
    else
    {
        ok = _arb_hypgeom_integrate_exp_log(r2, lambda2, mu2, al2, be2, m, 0, 0.5, -INFINITY, prec)
            && _arb_hypgeom_integrate_exp_log(r1, lambda, mu, al, be, m, 0, 0.5,
                el_log_abs(r2) + sh, prec);
    }

    if (ok)
    {
        if (!arb_is_zero(shift))
        {
            arb_exp(e, shift, prec + 30);
            arb_mul(r2, r2, e, prec + 30);
        }
        arb_add(res, r1, r2, prec + 30);
        ok = arb_is_finite(res);
    }

    arb_clear(r1);
    arb_clear(r2);
    arb_clear(e);
    return ok;
}

/* int_0^1 e^(z t) t^a1 (1 - t)^ba1 dt */
int
_arb_hypgeom_1f1_integral(arb_t res, const arb_t a1, const arb_t ba1, const arb_t z, slong prec)
{
    arb_struct mu[2], mu2[2];
    arb_t al[2], be[2], mz;
    int ok;

    arb_init(al[0]); arb_init(al[1]); arb_init(be[0]); arb_init(be[1]); arb_init(mz);

    /* t^a1 (1 - t)^ba1, and in s = 1 - t: e^z e^(-z s) s^ba1 (1 - s)^a1 */
    arb_one(be[0]);
    arb_one(al[1]);
    arb_set_si(be[1], -1);
    mu[0] = *a1; mu[1] = *ba1;
    mu2[0] = *ba1; mu2[1] = *a1;
    arb_neg(mz, z);

    ok = el_two_halves(res, z, mu, al[0], be[0], mz, mu2, al[0], be[0], z, 2, prec);

    arb_clear(al[0]); arb_clear(al[1]); arb_clear(be[0]); arb_clear(be[1]); arb_clear(mz);
    return ok;
}

/* int_0^inf e^(-z t) t^a1 (1 + t)^ba1 dt */
int
_arb_hypgeom_u_integral(arb_t res, const arb_t a1, const arb_t ba1, const arb_t z, slong prec)
{
    arb_struct mu[2];
    arb_t al[2], be[2], mz;
    int ok;

    arb_init(al[0]); arb_init(al[1]); arb_init(be[0]); arb_init(be[1]); arb_init(mz);

    arb_one(be[0]);
    arb_one(al[1]);
    arb_one(be[1]);
    mu[0] = *a1; mu[1] = *ba1;
    arb_neg(mz, z);

    ok = _arb_hypgeom_integrate_exp_log(res, mz, mu, al[0], be[0], 2, 0, INFINITY, -INFINITY, prec);

    arb_clear(al[0]); arb_clear(al[1]); arb_clear(be[0]); arb_clear(be[1]); arb_clear(mz);
    return ok;
}

/* int_0^1 t^b1 (1 - t)^cb1 (1 - z t)^nega dt */
int
_arb_hypgeom_2f1_integral(arb_t res, const arb_t b1, const arb_t cb1, const arb_t nega,
    const arb_t z, slong prec)
{
    arb_struct mu[3], mu2[3];
    arb_t al[3], be[3], al2[3], be2[3], zero;
    int ok, j;

    for (j = 0; j < 3; j++)
    {
        arb_init(al[j]); arb_init(be[j]); arb_init(al2[j]); arb_init(be2[j]);
    }
    arb_init(zero);

    /* t^b1 (1 - t)^cb1 (1 - z t)^nega, and in s = 1 - t:
       s^cb1 (1 - s)^b1 ((1 - z) + z s)^nega */
    arb_one(be[0]);
    arb_one(al[1]);
    arb_set_si(be[1], -1);
    arb_one(al[2]);
    arb_neg(be[2], z);
    mu[0] = *b1; mu[1] = *cb1; mu[2] = *nega;

    arb_one(be2[0]);
    arb_one(al2[1]);
    arb_set_si(be2[1], -1);
    arb_sub_ui(al2[2], z, 1, prec + 60);
    arb_neg(al2[2], al2[2]);
    arb_set(be2[2], z);
    mu2[0] = *cb1; mu2[1] = *b1; mu2[2] = *nega;

    ok = el_two_halves(res, zero, mu, al[0], be[0], zero, mu2, al2[0], be2[0], zero, 3, prec);

    for (j = 0; j < 3; j++)
    {
        arb_clear(al[j]); arb_clear(be[j]); arb_clear(al2[j]); arb_clear(be2[j]);
    }
    arb_clear(zero);
    return ok;
}
