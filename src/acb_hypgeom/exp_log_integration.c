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
    Rigorous quadrature of

        I = int_P exp(phi(t)) dt,
        phi(t) = k + lambda t + sum_{j < m} mu_j Log(al_j + be_j t),

    along a polygon P in the complex plane, with complex parameters
    (balls) and principal logarithms.  This is the complex counterpart
    of arb_hypgeom/exp_log_integration.c, which see; the differences
    are the following.

    The polygon.  The vertices are exact numbers.  The first and the
    last vertex can be points where a term vanishes (al_j + be_j t = 0:
    a branch point, where Re(mu_j) >= 0 is required).  The polygon must
    not meet the branch cuts of the logarithms, the rays on which
    al_j + be_j t is a negative number: this is verified on every cell
    as a part of its bound (the bound is infinite otherwise), so that
    the integrand is the analytic function which the formula defines on
    a neighbourhood of P.  It is up to the caller to relate this to the
    integral which is wanted (see the second part of this file: the
    cuts are rotated away from the path by replacing (al_j, be_j) by
    (c al_j, c be_j) and adjusting k).

    Bounds.  On the chord [m - h, m + h], with zeta = u / h for a point
    m + u and d_j = m + al_j / be_j,

        phi(m + u) - phi(m) = p zeta + q zeta^2 + sum_j mu_j r3(u / d_j),
        p = phi'(m) h,  q = phi''(m) h^2 / 2,
        r3(x) = log(1 + x) - x + x^2 / 2,  |r3(x)| <= |x|^3 / (3 (1 - |x|)).

    On the Bernstein ellipse with the parameter rho = X + Y,
    Y = sqrt(X^2 - 1), zeta = X cos(s) + i Y sin(s), and

        Re(p zeta) <= sqrt((X Re p)^2 + (Y Im p)^2),
        Re(q zeta^2) <= Re(q) / 2 + sqrt(|q|^2 (rho^4 + rho^-4) + 2 Re(q^2)) / 4

    (the maxima of the two trigonometric polynomials), |u| <= |h| X.  On
    a contour of steepest descent through a saddle, p = 0 and q < 0,
    and the second bound is |q| (X^2 - 1): the growth of the Gaussian
    off its axis.  On the chord itself, Re(p zeta + q zeta^2) <= |Re p|
    + max(Re q, 0).  The bound of the first order
    Re(phi(m + u) - phi(m)) <= |p| X + sum_j |mu_j| L(|u / d_j|),
    L(A) = -log(1 - A) - A, is used where it is smaller.

    The identity above holds for the principal logarithms where the
    disk |u| <= R around m does not meet a cut: R is less than |w| / |be_j|
    if Re(w) > 0 and than |Im(w)| / |be_j| otherwise, w = al_j + be_j m.

    Evaluation.  In dfloat arithmetic, in the local form
    exp((phi(m) - c) + lambda u + sum_j mu_j log(1 + u / d_j)) with
    phi(m) - c computed once per cell in arb arithmetic: the terms
    which cancel are then of the size of the variation of the exponent
    over the cell and not of phi itself.  In arb arithmetic, directly.

    End cells.  At a vertex E where the term s vanishes, the cell
    [E, E + tau] is bounded if it is negligible, and integrated term by
    term otherwise: with g(tau) = exp(phi without the term s) =
    sum_k g_k tau^k,

        int_0^tau1 (be_s tau)^mu_s g(tau) dtau
            = (be_s tau1)^mu_s tau1 sum_k g_k tau1^k / (k + mu_s + 1).
*/

#include <math.h>
#include <string.h>
#include "acb.h"
#include "acb_hypgeom.h"
#include "acb_hypgeom/impl.h"
#include "acb_calc/impl.h"
#include "dfloat.h"

#define CEL_MAXTERMS 3
#define CEL_LN2 0.69314718055994530942
/* the disk around an ellipse stays away from the singularities and the
   cuts by this factor (estimates) */
#define CEL_MARGIN 1.02

#define N 1
#define DB(name) d1b_##name
#define DBV(name) _d1b_##name
#define CELQ(name) celq_##name##_d1
#include "exp_log_quad.inc"
#undef N
#undef DB
#undef DBV
#undef CELQ

#define N 2
#define DB(name) d2b_##name
#define DBV(name) _d2b_##name
#define CELQ(name) celq_##name##_d2
#include "exp_log_quad.inc"
#undef N
#undef DB
#undef DBV
#undef CELQ

#define N 3
#define DB(name) d3b_##name
#define DBV(name) _d3b_##name
#define CELQ(name) celq_##name##_d3
#include "exp_log_quad.inc"
#undef N
#undef DB
#undef DBV
#undef CELQ

#define N 4
#define DB(name) d4b_##name
#define DBV(name) _d4b_##name
#define CELQ(name) celq_##name##_d4
#include "exp_log_quad.inc"
#undef N
#undef DB
#undef DBV
#undef CELQ

/* complex doubles */
typedef struct { double re, im; } cel_cd;

static inline cel_cd cel_mk(double re, double im) { cel_cd z; z.re = re; z.im = im; return z; }
static inline cel_cd cel_add(cel_cd a, cel_cd b) { return cel_mk(a.re + b.re, a.im + b.im); }
static inline cel_cd cel_sub(cel_cd a, cel_cd b) { return cel_mk(a.re - b.re, a.im - b.im); }
static inline cel_cd cel_mul(cel_cd a, cel_cd b) { return cel_mk(a.re * b.re - a.im * b.im, a.re * b.im + a.im * b.re); }
static inline cel_cd cel_scale(cel_cd a, double s) { return cel_mk(a.re * s, a.im * s); }
static inline double cel_abs(cel_cd a) { return hypot(a.re, a.im); }
static inline cel_cd cel_div(cel_cd a, cel_cd b)
{
    double d = b.re * b.re + b.im * b.im;
    return cel_mk((a.re * b.re + a.im * b.im) / d, (a.im * b.re - a.re * b.im) / d);
}
static inline cel_cd cel_log(cel_cd a) { return cel_mk(log(hypot(a.re, a.im)), atan2(a.im, a.re)); }

static cel_cd
cel_get(const acb_t x)
{
    return cel_mk(arf_get_d(arb_midref(acb_realref(x)), ARF_RND_NEAR),
        arf_get_d(arb_midref(acb_imagref(x)), ARF_RND_NEAR));
}

typedef struct
{
    int m;
    acb_srcptr lam, konst;
    acb_srcptr mu, al, be;
    cel_cd dlam, dkonst, dmu[CEL_MAXTERMS], dal[CEL_MAXTERMS], dbe[CEL_MAXTERMS];
    int end0, end1;          /* the term vanishing at the first or the
                                last vertex, or -1 */
    acb_srcptr E0, E1;       /* these vertices */
    /* The stable form, for large parameters: with the origin at (or
       close to) a saddle point, all al_j nonzero, rq_j = be_j / al_j
       and zeta_j = rq_j t,

           phi(t) - phi(0) = g1 t + sum_j mu_j (log(1 + zeta_j) - zeta_j),
           phi'(t) = g1 - sum_j mu_j rq_j zeta_j / (1 + zeta_j),

       g1 = phi'(0) = lambda + sum_j mu_j rq_j, where the terms of the
       first order, which cancel, have been summed in ball arithmetic
       (g1), and phi0 = phi(0) minus the reference value. */
    int stable;
    cel_cd g1, phi0, rq[CEL_MAXTERMS];
}
cel_kernel;

/* log(1 + x) - x */
static cel_cd
cel_l1pmx(cel_cd x)
{
    double ax = cel_abs(x);

    if (ax < 0.1)
    {
        /* -x^2 / 2 + x^3 / 3 - ... */
        cel_cd s = cel_mk(0.0, 0.0), pw = cel_mul(x, x);
        int k;

        for (k = 2; k < 40; k++)
        {
            cel_cd term = cel_scale(pw, ((k & 1) ? 1.0 : -1.0) / k);
            s = cel_add(s, term);
            if (cel_abs(term) < 1e-18 * cel_abs(s))
                break;
            pw = cel_mul(pw, x);
        }
        return s;
    }
    else
    {
        return cel_mk(0.5 * log1p(2.0 * x.re + x.re * x.re + x.im * x.im) - x.re,
            atan2(x.im, 1.0 + x.re) - x.im);
    }
}

/* Log(al_j + be_j t) - Log(al_j) - zeta_j with principal logarithms */
static cel_cd
cel_lpm(const cel_kernel * K, int j, cel_cd t)
{
    cel_cd z = cel_mul(K->rq[j], t), r = cel_l1pmx(z);

    if (cel_abs(z) >= 0.1)
    {
        cel_cd w = cel_mul(K->dal[j], cel_mk(1.0 + z.re, z.im));
        double d = atan2(w.im, w.re) - atan2(K->dal[j].im, K->dal[j].re) - atan2(z.im, 1.0 + z.re);
        r.im += 6.283185307179586477 * floor(d / 6.283185307179586477 + 0.5);
    }
    return r;
}

/* Heuristic: phi (minus the reference value, in the stable form), phi'
   and phi'' at t in double precision, without the term skip */
static void
cel_phi(cel_cd * phi, cel_cd * d1, cel_cd * d2, const cel_kernel * K, cel_cd t, int skip)
{
    cel_cd p = cel_add(K->dkonst, cel_mul(K->dlam, t)), q = K->dlam, r = cel_mk(0.0, 0.0), w, v;
    int j;

    if (K->stable)
    {
        p = cel_add(K->phi0, cel_mul(K->g1, t));
        q = K->g1;
        for (j = 0; j < K->m; j++)
        {
            cel_cd z, u, a;

            if (K->dmu[j].re == 0.0 && K->dmu[j].im == 0.0)
                continue;
            a = cel_mul(K->dmu[j], K->rq[j]);
            if (j == skip)
            {
                q = cel_sub(q, a);
                continue;
            }
            z = cel_mul(K->rq[j], t);
            u = cel_mk(1.0 + z.re, z.im);
            p = cel_add(p, cel_mul(K->dmu[j], cel_lpm(K, j, t)));
            q = cel_sub(q, cel_div(cel_mul(a, z), u));
            r = cel_sub(r, cel_div(cel_mul(a, K->rq[j]), cel_mul(u, u)));
        }
        *phi = p; *d1 = q; *d2 = r;
        return;
    }

    for (j = 0; j < K->m; j++)
    {
        if (j == skip || (K->dmu[j].re == 0.0 && K->dmu[j].im == 0.0))
            continue;
        w = cel_add(K->dal[j], cel_mul(K->dbe[j], t));
        v = cel_div(K->dbe[j], w);
        p = cel_add(p, cel_mul(K->dmu[j], cel_log(w)));
        v = cel_mul(K->dmu[j], v);
        q = cel_add(q, v);
        r = cel_sub(r, cel_mul(v, cel_div(K->dbe[j], w)));
    }

    *phi = p; *d1 = q; *d2 = r;
}

/* -log(1 - x) - x */
static double
cel_lrem(double x)
{
    if (!(x < 1.0))
        return INFINITY;
    if (x < 0.01)
        return 0.5 * x * x / (1.0 - x);
    return -log1p(-x) - x;
}

/* The quantities of the estimates for a chord [m - h, m + h] */
typedef struct
{
    double phire, ah;
    cel_cd p, q;
    double amu[CEL_MAXTERMS], dist[CEL_MAXTERMS], cut[CEL_MAXTERMS];
    int nt;
}
cel_local;

static void
cel_local_init(cel_local * W, const cel_kernel * K, cel_cd m, cel_cd h, int skip)
{
    cel_cd phi, d1, d2, w;
    int j;

    /* (the value with all the terms) */
    cel_phi(&phi, &d1, &d2, K, m, -1);
    W->phire = phi.re;
    if (skip >= 0)
        cel_phi(&w, &d1, &d2, K, m, skip);
    W->ah = cel_abs(h);
    W->p = cel_mul(d1, h);
    W->q = cel_scale(cel_mul(d2, cel_mul(h, h)), 0.5);
    W->nt = 0;

    for (j = 0; j < K->m; j++)
    {
        double amu = cel_abs(K->dmu[j]), abe = cel_abs(K->dbe[j]);

        if (j == skip || amu == 0.0 || abe == 0.0)
            continue;
        w = cel_add(K->dal[j], cel_mul(K->dbe[j], m));
        W->amu[W->nt] = amu;
        W->dist[W->nt] = cel_abs(w) / abe;
        W->cut[W->nt] = (w.re > 0.0) ? W->dist[W->nt] : fabs(w.im) / abe;
        W->nt++;
    }
}

/* An estimate of the bound for sup Re(phi - c) on the ellipse with the
   parameter X around the chord (on the chord itself if X <= 1);
   infinite if the disk of radius |h| X meets a singularity or a cut. */
static double
cel_local_est(const cel_local * W, double X, double cre)
{
    double R, S1 = 0.0, S3 = 0.0, B1, B2, A;
    cel_cd p = W->p, q = W->q;
    int j, chord = !(X > 1.0);

    if (chord)
        X = 1.0;
    R = W->ah * X;

    for (j = 0; j < W->nt; j++)
    {
        if (!(R * (chord ? 1.0 : CEL_MARGIN) < W->cut[j]))
            return INFINITY;
        A = R / W->dist[j];
        S1 += W->amu[j] * cel_lrem(A);
        S3 += W->amu[j] * A * A * A / (3.0 * (1.0 - A));
    }

    if (chord)
    {
        /* the maximum of Re(p) x + Re(q) x^2 on [-1, 1] */
        double a = fabs(p.re), b = q.re;

        B1 = a + S1;
        if (b >= 0.0)
            B2 = a + b;
        else if (a >= -2.0 * b)
            B2 = a + b;
        else
            B2 = -0.25 * a * a / b;
        B2 += S3;
    }
    else
    {
        double Y = sqrt(X * X - 1.0), rho = X + Y, r4 = rho * rho * rho * rho;
        double aq2 = q.re * q.re + q.im * q.im;

        B1 = cel_abs(p) * X + S1;
        B2 = 0.5 * q.re + hypot(X * p.re, Y * p.im)
            + 0.25 * sqrt(fmax(aq2 * (r4 + 1.0 / r4) + 2.0 * (q.re * q.re - q.im * q.im), 0.0)) + S3;
    }

    return W->phire - cre + fmin(B1, B2);
}

static double
cel_bound_est(const cel_kernel * K, cel_cd m, cel_cd h, double X, double cre, int skip)
{
    cel_local W;
    cel_local_init(&W, K, m, h, skip);
    return cel_local_est(&W, X, cre);
}

/* Rigorous: phi, phi' and phi'' (if not NULL) at m without the term skip */
static void
cel_phi_acb(acb_t phi, acb_t d1, acb_t d2, const cel_kernel * K, const acb_t m, int skip, slong prec)
{
    acb_t w, v;
    int j;

    acb_init(w);
    acb_init(v);

    acb_mul(phi, K->lam, m, prec);
    acb_add(phi, phi, K->konst, prec);
    if (d1 != NULL)
        acb_set(d1, K->lam);
    if (d2 != NULL)
        acb_zero(d2);

    for (j = 0; j < K->m; j++)
    {
        if (j == skip || acb_is_zero(K->mu + j))
            continue;
        acb_mul(w, K->be + j, m, prec);
        acb_add(w, w, K->al + j, prec);
        if (d1 != NULL)
        {
            acb_div(v, K->be + j, w, prec);
            acb_mul(v, v, K->mu + j, prec);
            acb_add(d1, d1, v, prec);
            if (d2 != NULL)
            {
                acb_mul(v, v, K->be + j, prec);
                acb_div(v, v, w, prec);
                acb_sub(d2, d2, v, prec);
            }
        }
        acb_log(w, w, prec);
        acb_addmul(phi, K->mu + j, w, prec);
    }

    acb_clear(w);
    acb_clear(v);
}

/* Rigorous: A >= R |be_j| / |w| for w = al_j + be_j m, infinite if the
   disk of radius R around m may meet the singularity or the cut of
   the term j (a constant term has neither). */
static void
cel_ratio(mag_t A, const cel_kernel * K, int j, const acb_t m, const mag_t R, slong prec)
{
    acb_t w;
    mag_t b, d, e;

    if (acb_is_zero(K->be + j))
    {
        mag_zero(A);
        return;
    }

    acb_init(w);
    mag_init(b);
    mag_init(d);
    mag_init(e);

    acb_mul(w, K->be + j, m, prec);
    acb_add(w, w, K->al + j, prec);
    acb_get_mag(b, K->be + j);
    acb_get_mag_lower(d, w);

    /* e = R |be_j| */
    mag_mul(e, R, b);

    if (mag_is_zero(d))
        mag_inf(A);
    else
    {
        mag_div(A, e, d);

        /* the cut: the disk around w with the radius e does not reach
           the negative numbers if Re(w) > 0 and e < |w|, or if
           e < |Im(w)| */
        if (!arb_is_positive(acb_realref(w)))
        {
            arb_get_mag_lower(d, acb_imagref(w));
            if (mag_cmp(e, d) >= 0)
                mag_inf(A);
        }
        else if (mag_cmp_2exp_si(A, 0) >= 0)
            mag_inf(A);
    }

    acb_clear(w);
    mag_clear(b);
    mag_clear(d);
    mag_clear(e);
}

/* Rigorous: res >= -log(1 - x) - x */
static void
cel_lrem_mag(mag_t res, const mag_t x)
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

/* Rigorous: M >= sup |exp(phi(t) - c)| on the ellipse with the
   parameter X around the chord [m - h, m + h] (on the chord if
   X <= 1), without the term skip. */
static void
cel_bound(mag_t M, const cel_kernel * K, const acb_t c, const acb_t m, const acb_t h,
    double X, int skip, slong prec)
{
    acb_t phi, d1, d2, p, q;
    arb_t B1, B2, t, u, ax, ay;
    mag_t R, A, S1, S3, y, z, amu;
    int j, chord = !(X > 1.0);

    if (chord)
        X = 1.0;

    acb_init(phi); acb_init(d1); acb_init(d2); acb_init(p); acb_init(q);
    arb_init(B1); arb_init(B2); arb_init(t); arb_init(u); arb_init(ax); arb_init(ay);
    mag_init(R); mag_init(A); mag_init(S1); mag_init(S3); mag_init(y); mag_init(z); mag_init(amu);

    cel_phi_acb(phi, d1, d2, K, m, skip, prec);
    acb_sub(phi, phi, c, prec);
    acb_mul(p, d1, h, prec);
    acb_mul(q, d2, h, prec);
    acb_mul(q, q, h, prec);
    acb_mul_2exp_si(q, q, -1);

    /* R = |h| X */
    acb_get_mag(R, h);
    mag_set_d(y, X);
    mag_mul(R, R, y);

    for (j = 0; j < K->m; j++)
    {
        if (j == skip || acb_is_zero(K->mu + j))
            continue;

        cel_ratio(A, K, j, m, R, prec);
        acb_get_mag(amu, K->mu + j);

        cel_lrem_mag(y, A);
        mag_addmul(S1, y, amu);

        /* A^3 / (3 (1 - A)) */
        if (mag_cmp_2exp_si(A, 0) >= 0)
            mag_inf(z);
        else
        {
            mag_one(z);
            mag_sub_lower(z, z, A);
            mag_mul_ui_lower(z, z, 3);
            mag_div(z, A, z);
            mag_mul(z, z, A);
            mag_mul(z, z, A);
        }
        mag_addmul(S3, z, amu);
    }

    if (chord)
    {
        /* B1 = |Re p| + S1; B2 = S3 + the maximum of a x + b x^2 on
           [-1, 1], a = |Re p|, b = Re q: a + b if b >= 0 or a >= 2 |b|,
           and at most min(a, a^2 / (4 |b|)) (the maximum without the
           constraint) otherwise */
        arb_abs(B1, acb_realref(p));
        if (arb_is_nonnegative(acb_realref(q)))
            arb_add(B2, B1, acb_realref(q), prec);
        else if (arb_is_negative(acb_realref(q)))
        {
            arb_mul_2exp_si(t, acb_realref(q), 1);
            arb_add(t, t, B1, prec);
            if (arb_is_nonnegative(t))
                arb_add(B2, B1, acb_realref(q), prec);
            else
            {
                arb_mul(t, B1, B1, prec);
                arb_div(t, t, acb_realref(q), prec);
                arb_mul_2exp_si(t, t, -2);
                arb_neg(t, t);
                arb_min(B2, B1, t, prec);
            }
        }
        else
            arb_set(B2, B1);
    }
    else
    {
        /* ax = X, ay = Y */
        arb_set_d(ax, X);
        arb_mul(ay, ax, ax, prec);
        arb_sub_ui(ay, ay, 1, prec);
        arb_sqrt(ay, ay, prec);

        /* B1 = |p| X */
        acb_abs(B1, p, prec);
        arb_mul(B1, B1, ax, prec);

        /* B2 = Re(q) / 2 + sqrt((X Re p)^2 + (Y Im p)^2)
                + sqrt(|q|^2 (rho^4 + rho^-4) + 2 Re(q^2)) / 4 */
        arb_mul(t, ax, acb_realref(p), prec);
        arb_mul(u, ay, acb_imagref(p), prec);
        arb_hypot(B2, t, u, prec);
        arb_mul_2exp_si(t, acb_realref(q), -1);
        arb_add(B2, B2, t, prec);

        arb_add(t, ax, ay, prec);
        arb_pow_ui(t, t, 4, prec);
        arb_inv(u, t, prec);
        arb_add(t, t, u, prec);          /* rho^4 + rho^-4 */
        arb_mul(u, acb_realref(q), acb_realref(q), prec);
        arb_mul(ax, acb_imagref(q), acb_imagref(q), prec);
        arb_add(ay, u, ax, prec);        /* |q|^2 */
        arb_sub(u, u, ax, prec);         /* Re(q^2) */
        arb_mul(t, t, ay, prec);
        arb_mul_2exp_si(u, u, 1);
        arb_add(t, t, u, prec);
        /* (nonnegative; rounding) */
        arb_zero(u);
        arb_max(t, t, u, prec);
        arb_sqrt(t, t, prec);
        arb_mul_2exp_si(t, t, -2);
        arb_add(B2, B2, t, prec);
    }

    arb_add_error_mag(B1, S1);
    arb_add_error_mag(B2, S3);
    /* (upper bounds) */
    arb_get_ubound_arf(arb_midref(B1), B1, prec);
    mag_zero(arb_radref(B1));
    arb_get_ubound_arf(arb_midref(B2), B2, prec);
    mag_zero(arb_radref(B2));
    arb_min(B1, B1, B2, prec);

    arb_add(t, acb_realref(phi), B1, prec);
    arb_exp(t, t, prec);
    arb_get_mag(M, t);
    if (!arb_is_finite(t))
        mag_inf(M);

    acb_clear(phi); acb_clear(d1); acb_clear(d2); acb_clear(p); acb_clear(q);
    arb_clear(B1); arb_clear(B2); arb_clear(t); arb_clear(u); arb_clear(ax); arb_clear(ay);
    mag_clear(R); mag_clear(A); mag_clear(S1); mag_clear(S3); mag_clear(y); mag_clear(z); mag_clear(amu);
}

typedef struct
{
    const cel_kernel * K;
    acb_srcptr c;            /* the reference value (exact) */
    acb_t kc;                /* konst - c */
    cel_cd cd;
    slong prec;
    slong bprec;             /* for the bounds */
    double logtol;           /* log of the tolerance for a cell, relative
                                to |exp(c)| */
    double magd;             /* the size of the terms of phi */
    int imax;                /* the largest rule */
    slong budget;            /* of nodes */
    acb_t sum;
    mag_t err;
    void * scratch;
}
cel_ctx;

/* The rule gi on the chord [m - h, m + h] (m and h exact), where the
   integrand is about exp(lmax) at most.  Returns 0 on failure. */
static int
cel_cell_value(acb_t res, cel_ctx * ctx, const acb_t m, const acb_t h, cel_cd md, cel_cd hd,
    int gi, double lmax)
{
    const cel_kernel * K = ctx->K;
    double bits, lmag, lmag2, Amax, Aj[CEL_MAXTERMS], ah = cel_abs(hd), taper = fmin(lmax / CEL_LN2 + 2.0, 0.0);
    slong n = _acb_calc_gl_steps[gi], k, wp, nt = 0;
    int j, N, ok = 1, idx[CEL_MAXTERMS], nx = 0, poly = 0;
    cel_cd w;

    for (j = 0; j < K->m; j++)
        if (!acb_is_zero(K->mu + j))
            idx[nx++] = j;

    /* the size of the terms of the exponent in the local form, and in
       the form where the terms of the first order are summed */
    lmag = cel_abs(K->dlam) * ah + 2.0 + fabs(lmax);
    lmag2 = 2.0 + fabs(lmax);
    Amax = 0.0;
    for (j = 0; j < nx; j++)
    {
        double A;
        if (K->stable)
        {
            cel_cd z = cel_mul(K->rq[idx[j]], md);
            A = ah * cel_abs(K->rq[idx[j]]) / hypot(1.0 + z.re, z.im);
        }
        else
        {
            w = cel_add(K->dal[idx[j]], cel_mul(K->dbe[idx[j]], md));
            A = ah * cel_abs(K->dbe[idx[j]]) / cel_abs(w);
        }
        Aj[j] = A;
        Amax = fmax(Amax, A);
        /* (the logarithm of 1 + q x has an absolute error) */
        lmag += cel_abs(K->dmu[idx[j]]) * (1.0 + A / (1.0 - fmin(A, 0.99)));
        lmag2 += cel_abs(K->dmu[idx[j]]) * A * A / (1.0 - fmin(A, 0.99));
    }

    bits = ctx->prec + 10.0 + log2(lmag);
    bits = fmax(bits + taper, log2(lmag) + 14.0);
    /* (not where the integrand relative to exp(c) is too small for a
       double: at a high precision, a cell far below the highest point) */
    N = (dfloat_is_supported() && bits <= 4 * 53 && gi < ACB_CALC_GL_DFLOAT_STEPS
            && taper > -900.0)
        ? FLINT_MAX(1, (int) ceil(bits / 53.0)) : 0;

    /* The exponent as a polynomial in x, exp(sum_k c_k x^k) with
       c_0 = phi(m) - c, c_1 = phi'(m) h and
       c_k = sum_j mu_j (-1)^(k+1) q_j^k / k, q_j = be_j h / (al_j + be_j m),
       for large parameters: the coefficients are of the size of the
       exponent, where the terms of the local form cancel.  The rest
       after the terms to x^nt is at most
       sum_j |mu_j| |q_j|^(nt+1) / ((nt + 1) (1 - |q_j|)). */
    if (K->stable && dfloat_is_supported() && gi < ACB_CALC_GL_DFLOAT_STEPS && Amax <= 0.3
        && taper > -900.0)
    {
        cel_cd pd, d1, d2;
        double bits2, tb, tol = ldexp(1.0, -(int) FLINT_MIN(ctx->prec + 12, 1000));
        int N2;

        cel_phi(&pd, &d1, &d2, K, md, -1);
        lmag2 += cel_abs(d1) * ah;
        bits2 = ctx->prec + 10.0 + log2(lmag2);
        bits2 = fmax(bits2 + taper, log2(lmag2) + 14.0);
        N2 = (bits2 <= 4 * 53) ? FLINT_MAX(1, (int) ceil(bits2 / 53.0)) : 0;

        if (N2 != 0 && (N == 0 || N2 < N))
        {
            for (nt = 2; nt < 48; nt++)
            {
                tb = 0.0;
                for (j = 0; j < nx; j++)
                    tb += cel_abs(K->dmu[idx[j]]) * pow(Aj[j], nt + 1.0) / ((nt + 1.0) * (1.0 - Aj[j]));
                if (tb < tol)
                    break;
            }
            if (nt < 48)
            {
                poly = 1;
                N = N2;
            }
        }
    }

    if (poly)
    {
        acb_ptr e;
        acb_t q, pw, d1;
        mag_t tr, xm, tm;
        slong hp = 53 * N + 10 + (slong) log2(ctx->magd) + 10, lp = 53 * N + 20;

        e = _acb_vec_init(nt + 1);
        acb_init(q); acb_init(pw); acb_init(d1);
        mag_init(tr); mag_init(xm); mag_init(tm);

        cel_phi_acb(e, d1, NULL, K, m, -1, hp);
        acb_sub(e, e, ctx->c, hp);
        acb_mul(e + 1, d1, h, hp);

        for (j = 0; j < nx; j++)
        {
            acb_mul(q, K->be + idx[j], m, hp);
            acb_add(q, q, K->al + idx[j], hp);
            acb_div(q, K->be + idx[j], q, lp);
            acb_mul(q, q, h, lp);

            acb_mul(pw, q, K->mu + idx[j], lp);
            for (k = 2; k <= nt; k++)
            {
                /* (-1)^(k+1) mu_j q^k / k */
                acb_mul(pw, pw, q, lp);
                acb_neg(pw, pw);
                acb_div_ui(d1, pw, k, lp);
                acb_add(e + k, e + k, d1, lp);
            }

            acb_get_mag(xm, q);
            if (mag_cmp_2exp_si(xm, 0) >= 0)
                ok = 0;
            else
            {
                mag_one(tm);
                mag_sub_lower(tm, tm, xm);
                mag_pow_ui(xm, xm, nt + 1);
                mag_div(xm, xm, tm);
                mag_div_ui(xm, xm, nt + 1);
                acb_get_mag(tm, K->mu + idx[j]);
                mag_mul(xm, xm, tm);
                mag_add(tr, tr, xm);
            }
        }
        acb_add_error_mag(e, tr);

        if (ok)
        {
            _acb_calc_gl_exp_poly_dfloat(res, e, nt + 1, -1.0, 0.0, 1.0, 0.0, gi, N, ctx->scratch);
            acb_mul(res, res, h, 53 * N + 10);
        }

        _acb_vec_clear(e, nt + 1);
        acb_clear(q); acb_clear(pw); acb_clear(d1);
        mag_clear(tr); mag_clear(xm); mag_clear(tm);
    }
    else if (N != 0)
    {
        const void * xs, * ws;
        acb_t E, lh;
        acb_ptr q, mu;
        slong hp = 53 * N + 10 + (slong) log2(ctx->magd) + 10;

        acb_init(E);
        acb_init(lh);
        q = _acb_vec_init(CEL_MAXTERMS);
        mu = _acb_vec_init(CEL_MAXTERMS);

        /* E = phi(m) - c, and q_j = be_j h / (al_j + be_j m) */
        cel_phi_acb(E, NULL, NULL, K, m, -1, hp);
        acb_sub(E, E, ctx->c, hp);
        acb_mul(lh, K->lam, h, hp);
        for (j = 0; j < nx; j++)
        {
            acb_mul(q + j, K->be + idx[j], m, hp);
            acb_add(q + j, q + j, K->al + idx[j], hp);
            acb_div(q + j, K->be + idx[j], q + j, hp);
            acb_mul(q + j, q + j, h, hp);
            acb_set(mu + j, K->mu + idx[j]);
        }

        _acb_calc_gl_nodes_dfloat(&xs, &ws, N, gi);
        if (N == 1) ok = celq_direct_d1(res, E, lh, mu, q, nx, xs, ws, n, ctx->scratch);
        else if (N == 2) ok = celq_direct_d2(res, E, lh, mu, q, nx, xs, ws, n, ctx->scratch);
        else if (N == 3) ok = celq_direct_d3(res, E, lh, mu, q, nx, xs, ws, n, ctx->scratch);
        else ok = celq_direct_d4(res, E, lh, mu, q, nx, xs, ws, n, ctx->scratch);
        acb_mul(res, res, h, 53 * N + 10);

        acb_clear(E);
        acb_clear(lh);
        _acb_vec_clear(q, CEL_MAXTERMS);
        _acb_vec_clear(mu, CEL_MAXTERMS);
    }
    else
    {
        /* With arb arithmetic, at a precision which is lowered per
           node: the exponent as the dot product
           (k - c) + lambda t + sum_j mu_j log(al_j + be_j t). */
        arb_srcptr xs, ws;
        acb_struct xv[CEL_MAXTERMS + 1], yv[CEL_MAXTERMS + 1];
        acb_t t, v, u[CEL_MAXTERMS];
        slong kk, p;
        double gmag = log2(ctx->magd);

        acb_init(t);
        acb_init(v);
        for (j = 0; j < CEL_MAXTERMS; j++)
            acb_init(u[j]);

        xv[0] = *K->lam;
        for (j = 0; j < nx; j++)
            xv[j + 1] = K->mu[idx[j]];

        wp = (slong) (ctx->prec + 10.0 + gmag) + 8;
        _acb_calc_gl_nodes(&xs, &ws, gi, wp);
        acb_zero(res);

        for (k = 0; k < n; k++)
        {
            double xd;
            cel_cd pd, d1, d2;

            kk = (2 * k < n) ? k : n - 1 - k;
            xd = arf_get_d(arb_midref(xs + kk), ARF_RND_NEAR);
            if (2 * k >= n)
                xd = -xd;
            cel_phi(&pd, &d1, &d2, K, cel_add(md, cel_scale(hd, xd)), -1);
            p = wp + (slong) fmax(fmin((pd.re - ctx->cd.re) / CEL_LN2 + 2.0, 0.0), -1e9);
            p = FLINT_MAX(p, (slong) (gmag + 14.0));

            /* t = m +- h x */
            acb_mul_arb(t, h, xs + kk, p);
            if (2 * k >= n)
                acb_sub(t, m, t, p);
            else
                acb_add(t, m, t, p);

            for (j = 0; j < nx; j++)
            {
                acb_mul(u[j], K->be + idx[j], t, p);
                acb_add(u[j], u[j], K->al + idx[j], p);
                acb_log(u[j], u[j], p);
                yv[j + 1] = *u[j];
            }

            yv[0] = *t;
            acb_dot(v, ctx->kc, 0, xv, 1, yv, 1, nx + 1, p);
            acb_exp(v, v, p);
            acb_addmul_arb(res, v, ws + kk, wp);
        }
        acb_mul(res, res, h, wp);

        acb_clear(t);
        acb_clear(v);
        for (j = 0; j < CEL_MAXTERMS; j++)
            acb_clear(u[j]);
    }

    return ok && acb_is_finite(res);
}

/* The end cell from the vertex E where the term s vanishes to
   E + tau1 (tau1 exact), by the power series of the other factor
   g(tau) = exp(phi(E + tau) without the term s) = sum_k g_k tau^k.
   Since D g' = N g with the polynomials D = prod_{j != s} (a_j + be_j tau),
   a_j = al_j + be_j E, and N = lambda D + sum_j mu_j be_j prod_{i != j, s}
   (a_i + be_i tau),

       D_0 (k + 1) g_{k+1} = sum_i N_i g_{k-i} - sum_{i >= 1} D_i (k - i + 1) g_{k-i+1}.

   On the disk |tau| <= r = 4 |tau1| (which must not meet the
   singularities and the cuts of the other terms), with g1 the
   logarithmic derivative of g at 0 (where lambda and the terms
   mu_j be_j / a_j can cancel),

       |g(tau)| <= M = |g_0| exp(|g1| r + sum_j |mu_j| L(A_j)),
       A_j = |be_j| r / |a_j|,  L(A) = -log(1 - A) - A,

   and |g_k| <= M r^-k, so that the terms from K on add up to at most
   M 4^-K (4 / 3) / K, times the common factor (Re(mu_s) >= 0).  The
   integral is added to the sum (negated if rev: the cell [E + tau1, E]).
   Returns 0 if this does not apply. */
static int
cel_end_series(cel_ctx * ctx, int s, const acb_t E, const acb_t tau1, int rev)
{
    const cel_kernel * K = ctx->K;
    int idx[CEL_MAXTERMS], no = 0, i, j, ok = 1;
    slong k, nterms, wp, wp2;
    double grow, pref, r, ld;
    cel_cd Ed = cel_get(E), td = cel_get(tau1), w;
    acb_t D[3], Nn[3], a[CEL_MAXTERMS], t, u, sum, pw, lk, e0;
    arb_t M;
    acb_ptr g;
    mag_t tail, R, A, y;

    ld = cel_abs(td);
    r = 4.0 * ld;

    /* estimates: the growth of g on the disk */
    {
        cel_cd pd, d1, d2;
        cel_phi(&pd, &d1, &d2, K, Ed, s);
        grow = cel_abs(d1) * r;
    }
    for (j = 0; j < K->m; j++)
    {
        double amu = cel_abs(K->dmu[j]), abe = cel_abs(K->dbe[j]), cut;

        if (j == s || amu == 0.0)
            continue;
        w = cel_add(K->dal[j], cel_mul(K->dbe[j], Ed));
        cut = (w.re > 0.0) ? cel_abs(w) : fabs(w.im);
        if (!(abe * r < 0.5 * cut))
            return 0;
        grow += amu * cel_lrem(abe * r / cel_abs(w));
        idx[no++] = j;
    }
    if (no > 2 || !(grow < 60.0))
        return 0;

    for (i = 0; i < 3; i++)
    {
        acb_init(D[i]);
        acb_init(Nn[i]);
    }
    for (i = 0; i < CEL_MAXTERMS; i++)
        acb_init(a[i]);
    acb_init(t); acb_init(u); acb_init(sum); acb_init(pw); acb_init(lk); acb_init(e0);
    arb_init(M);
    mag_init(tail); mag_init(R); mag_init(A); mag_init(y);

    /* The number of terms, from the size of the first one (in ball
       arithmetic: the terms of the exponent cancel):
       e0 = k + lambda E + sum_j mu_j log(a_j) - c, and
       mu_s log(be_s tau1) */
    wp = ctx->bprec + 30;
    acb_mul(e0, K->lam, E, wp);
    acb_add(e0, e0, ctx->kc, wp);
    for (i = 0; i < no; i++)
    {
        j = idx[i];
        acb_mul(a[i], K->be + j, E, wp);
        acb_add(a[i], a[i], K->al + j, wp);
        acb_log(u, a[i], wp);
        acb_addmul(e0, u, K->mu + j, wp);
    }
    acb_mul(t, K->be + s, tau1, wp);
    acb_log(t, t, wp);
    acb_mul(t, t, K->mu + s, wp);
    pref = log(ld) + arf_get_d(arb_midref(acb_realref(t)), ARF_RND_NEAR)
        + arf_get_d(arb_midref(acb_realref(e0)), ARF_RND_NEAR);

    nterms = 2;
    if (pref - ctx->logtol + grow < 1e6)
        nterms = (slong) ((grow + pref - ctx->logtol) / log(4.0)) + 3;
    if (nterms < 2)
        nterms = 2;
    if (nterms > 20000 || !acb_is_finite(e0) || !(pref < 1e300))
    {
        for (i = 0; i < 3; i++)
        {
            acb_clear(D[i]);
            acb_clear(Nn[i]);
        }
        for (i = 0; i < CEL_MAXTERMS; i++)
            acb_clear(a[i]);
        acb_clear(t); acb_clear(u); acb_clear(sum); acb_clear(pw); acb_clear(lk); acb_clear(e0);
        arb_clear(M);
        mag_clear(tail); mag_clear(R); mag_clear(A); mag_clear(y);
        return 0;
    }

    wp = ctx->prec + 30 + 2 * FLINT_BIT_COUNT(nterms) + (slong) fmin(grow / CEL_LN2, 1e9);
    g = _acb_vec_init(nterms + 1);

    /* a_j, D and N (with the precision of the parameters: lambda and
       the mu_j be_j / a_j can cancel in N) */
    wp2 = wp + ctx->bprec;
    acb_one(D[0]);
    for (i = 0; i < no; i++)
    {
        j = idx[i];
        acb_mul(a[i], K->be + j, E, wp2);
        acb_add(a[i], a[i], K->al + j, wp2);
        /* D *= (a_j + be_j tau) */
        acb_mul(D[2], D[1], K->be + j, wp2);
        acb_mul(D[1], D[1], a[i], wp2);
        acb_addmul(D[1], D[0], K->be + j, wp2);
        acb_mul(D[0], D[0], a[i], wp2);
    }
    for (i = 0; i < 3; i++)
        acb_mul(Nn[i], D[i], K->lam, wp2);
    if (no == 1)
    {
        acb_addmul(Nn[0], K->mu + idx[0], K->be + idx[0], wp2);
    }
    else if (no == 2)
    {
        for (i = 0; i < 2; i++)
        {
            int j1 = idx[i], j2 = idx[1 - i];
            /* mu_j1 be_j1 (a_j2 + be_j2 tau) */
            acb_mul(t, K->mu + j1, K->be + j1, wp2);
            acb_addmul(Nn[0], t, a[1 - i], wp2);
            acb_addmul(Nn[1], t, K->be + j2, wp2);
        }
    }

    /* g_0 = exp(k + lambda E + sum_j mu_j log(a_j) - c) */
    wp2 = wp + ctx->bprec;
    acb_mul(t, K->lam, E, wp2);
    acb_add(t, t, ctx->kc, wp2);
    for (i = 0; i < no; i++)
    {
        acb_mul(u, K->be + idx[i], E, wp2);
        acb_add(u, u, K->al + idx[i], wp2);
        acb_log(u, u, wp2);
        acb_addmul(t, u, K->mu + idx[i], wp2);
    }
    acb_exp(g, t, wp);

    for (k = 0; k < nterms; k++)
    {
        acb_zero(t);
        for (i = 0; i <= 2 && i <= k; i++)
            acb_addmul(t, Nn[i], g + k - i, wp);
        for (i = 1; i <= 2 && i <= k; i++)
        {
            acb_mul_ui(u, g + k - i + 1, k - i + 1, wp);
            acb_submul(t, D[i], u, wp);
        }
        acb_div(t, t, D[0], wp);
        acb_div_ui(g + k + 1, t, k + 1, wp);
    }

    /* sum_k g_k tau1^k / (k + mu_s + 1) */
    acb_one(lk);
    for (k = 0; k < nterms; k++)
    {
        acb_add_ui(t, K->mu + s, k + 1, wp);
        acb_div(t, g + k, t, wp);
        acb_addmul(sum, t, lk, wp);
        acb_mul(lk, lk, tau1, wp);
    }

    /* M, and the tail */
    acb_get_mag(R, tau1);
    mag_mul_2exp_si(R, R, 2);
    cel_phi_acb(t, u, NULL, K, E, s, wp2);
    acb_get_mag(y, u);
    mag_mul(tail, y, R);
    for (i = 0; i < no; i++)
    {
        j = idx[i];
        /* A_j (infinite if the disk meets the singularity or the cut) */
        cel_ratio(A, K, j, E, R, wp);
        if (mag_cmp_2exp_si(A, 0) >= 0)
            ok = 0;
        else
        {
            cel_lrem_mag(y, A);
            acb_get_mag(A, K->mu + j);
            mag_addmul(tail, y, A);
        }
    }
    if (!arb_is_nonnegative(acb_realref(K->mu + s)))
        ok = 0;

    /* M 4^-K (4 / 3) / K */
    arf_set_mag(arb_midref(M), tail);
    mag_zero(arb_radref(M));
    arb_exp(M, M, 64);
    acb_abs(acb_realref(t), g, 64);
    arb_mul(M, M, acb_realref(t), 64);
    arb_mul_2exp_si(M, M, -2 * nterms + 2);
    arb_div_ui(M, M, 3 * nterms, 64);
    arb_get_mag(tail, M);
    acb_add_error_mag(sum, tail);

    /* times tau1 (be_s tau1)^mu_s; the chord must not lie on the cut */
    acb_mul(t, K->be + s, tau1, wp);
    if (!arb_is_positive(acb_realref(t)) && arb_contains_zero(acb_imagref(t)))
        ok = 0;
    acb_log(pw, t, wp2);
    acb_mul(pw, pw, K->mu + s, wp2);
    acb_exp(pw, pw, wp);
    acb_mul(sum, sum, pw, wp);
    acb_mul(sum, sum, tau1, wp);

    ok = ok && acb_is_finite(sum);
    if (ok)
    {
        if (rev)
            acb_sub(ctx->sum, ctx->sum, sum, ctx->prec + 30);
        else
            acb_add(ctx->sum, ctx->sum, sum, ctx->prec + 30);
    }

    for (i = 0; i < 3; i++)
    {
        acb_clear(D[i]);
        acb_clear(Nn[i]);
    }
    for (i = 0; i < CEL_MAXTERMS; i++)
        acb_clear(a[i]);
    acb_clear(t); acb_clear(u); acb_clear(sum); acb_clear(pw); acb_clear(lk); acb_clear(e0);
    arb_clear(M);
    mag_clear(tail); mag_clear(R); mag_clear(A); mag_clear(y);
    _acb_vec_clear(g, nterms + 1);
    return ok;
}

/* The number of nodes for the chord [m - h, m + h] (estimates): of the
   smallest rule whose error bound, minimised over the ellipse parameter
   X, is below the tolerance; 0 for a cell which is negligible (if this
   is allowed), and infinity if there is no rule. */
static double
cel_cell_plan(int * rule, double * bestX, const cel_ctx * ctx, cel_cd md, cel_cd hd, int negl)
{
    static const double Xs[] = { 1.25, 1.5, 2.0, 3.0, 4.0, 6.0, 8.0, 12.0, 24.0,
        1e2, 1e3, 1e4, 1e6, 1e8 };
    cel_local W;
    double bestn = INFINITY, ah = cel_abs(hd), est;
    int ix, gi;

    cel_local_init(&W, ctx->K, md, hd, -1);

    if (negl)
    {
        est = cel_local_est(&W, 0.0, ctx->cd.re);
        if (log(2.0 * ah) + est <= ctx->logtol - 1.0)
            return 0.0;
    }

    for (ix = 0; ix < (int) (sizeof(Xs) / sizeof(double)); ix++)
    {
        double X = Xs[ix], rho = X + sqrt(X * X - 1.0), need;

        est = cel_local_est(&W, X, ctx->cd.re);
        if (!(est < INFINITY))
            break;
        need = log(64.0 / 15.0) + est + log(ah) - log(rho - 1.0) - ctx->logtol + 0.5;
        need = 0.5 * (need / log(rho) + 1.0);
        if (need < bestn)
        {
            bestn = need;
            *bestX = X;
        }
        else if (need > 1.5 * bestn + 2.0)
            break;
    }

    for (gi = 0; gi <= ctx->imax; gi++)
        if (_acb_calc_gl_steps[gi] >= bestn)
            break;

    if (gi > ctx->imax)
        return INFINITY;

    *rule = gi;
    return (double) _acb_calc_gl_steps[gi];
}

/* is tol >= M, for the tolerance exp(logtol)? */
static int
cel_below_tol(const mag_t M, double logtol)
{
    arb_t e;
    mag_t t;
    int r;

    arb_init(e);
    mag_init(t);
    arb_set_d(e, logtol);
    arb_exp(e, e, 64);
    arb_get_mag_lower(t, e);
    r = (mag_cmp(M, t) <= 0);
    arb_clear(e);
    mag_clear(t);
    return r;
}

/* The cell [u, v] (exact numbers): bounded, integrated or bisected.
   Returns 0 on failure. */
static int
cel_cell(cel_ctx * ctx, const acb_t u, const acb_t v, int depth)
{
    const cel_kernel * K = ctx->K;
    cel_cd ud = cel_get(u), vd = cel_get(v), md, hd, phi, d1, d2;
    double sup, ah;
    int s0 = (K->end0 >= 0 && acb_equal(u, K->E0) && !acb_is_zero(K->mu + K->end0));
    int s1 = (K->end1 >= 0 && acb_equal(v, K->E1) && !acb_is_zero(K->mu + K->end1));
    int i, ok = 1, done = 0, negl = 1;
    acb_t m, h;
    mag_t R, M;

    md = cel_scale(cel_add(ud, vd), 0.5);
    hd = cel_scale(cel_sub(vd, ud), 0.5);
    ah = cel_abs(hd);

    /* (also: cells which doubles do not resolve) */
    ctx->budget -= 10;   /* (the work, in nodes) */
    if (depth > 4000 || ctx->budget < 0 || !(ah > 1e-13 * (cel_abs(ud) + cel_abs(vd))))
        return 0;

    acb_init(m);
    acb_init(h);
    mag_init(R);
    mag_init(M);

    /* (exactly) */
    arf_add(arb_midref(acb_realref(m)), arb_midref(acb_realref(u)), arb_midref(acb_realref(v)), ARF_PREC_EXACT, ARF_RND_DOWN);
    arf_add(arb_midref(acb_imagref(m)), arb_midref(acb_imagref(u)), arb_midref(acb_imagref(v)), ARF_PREC_EXACT, ARF_RND_DOWN);
    acb_mul_2exp_si(m, m, -1);
    arf_sub(arb_midref(acb_realref(h)), arb_midref(acb_realref(v)), arb_midref(acb_realref(u)), ARF_PREC_EXACT, ARF_RND_DOWN);
    arf_sub(arb_midref(acb_imagref(h)), arb_midref(acb_imagref(v)), arb_midref(acb_imagref(u)), ARF_PREC_EXACT, ARF_RND_DOWN);
    acb_mul_2exp_si(h, h, -1);
    acb_get_mag(R, h);

    if (s0 && s1)
    {
        /* (both ends: bisected) */
    }
    else if (s0 || s1)
    {
        /* The cell at an end E where the term j vanishes: with
           tau1 = the other end - E, its factor is at most
           |(be_j tau1)^mu_j| there. */
        int j = s0 ? K->end0 : K->end1;
        acb_srcptr E = s0 ? K->E0 : K->E1;
        cel_cd Ed = cel_get(E);
        double len = 2.0 * ah, est, q;

        /* (be_j tau1 = 2 (al_j + be_j m)) */
        est = cel_bound_est(K, md, hd, 0.0, ctx->cd.re, j);
        est += log(len) + K->dmu[j].re * CEL_LN2;

        if (est <= ctx->logtol - 2.0)
        {
            acb_t f, t1;
            mag_t t;

            acb_init(f);
            acb_init(t1);
            mag_init(t);

            cel_bound(M, K, ctx->c, m, h, 0.0, j, ctx->bprec);
            /* f = (be_j tau1)^mu_j */
            acb_mul_2exp_si(t1, h, 1);
            if (!s0)
                acb_neg(t1, t1);
            acb_mul(f, t1, K->be + j, ctx->bprec);
            if (!arb_is_positive(acb_realref(f)) && arb_contains_zero(acb_imagref(f)))
                acb_indeterminate(f);
            acb_pow(f, f, K->mu + j, ctx->bprec);
            acb_get_mag(t, f);
            mag_mul(M, M, t);
            mag_mul(M, M, R);
            mag_mul_2exp_si(M, M, 1);
            if (acb_is_finite(f) && arb_is_nonnegative(acb_realref(K->mu + j))
                && cel_below_tol(M, ctx->logtol))
            {
                mag_add(ctx->err, ctx->err, M);
                done = 1;
            }

            acb_clear(f);
            acb_clear(t1);
            mag_clear(t);
        }

        /* not negligible: by the power series at the end, if the cell
           is short enough for it */
        if (!done)
        {
            /* (the growth of the regular factor on the disk of radius 4 len) */
            cel_phi(&phi, &d1, &d2, K, Ed, j);
            q = cel_abs(d1) * 4.0 * len;
            for (i = 0; i < K->m; i++)
                if (i != j && cel_abs(K->dmu[i]) != 0.0)
                    q += cel_abs(K->dmu[i]) * cel_lrem(fmin(4.0 * len * cel_abs(K->dbe[i])
                        / cel_abs(cel_add(K->dal[i], cel_mul(K->dbe[i], Ed))), 2.0));

            if (q <= 8.0)
            {
                acb_t t1;
                acb_init(t1);
                acb_mul_2exp_si(t1, h, 1);
                if (!s0)
                    acb_neg(t1, t1);
                done = cel_end_series(ctx, j, E, t1, !s0);
                acb_clear(t1);
            }
        }
    }
    else
    {
        double bestX = 0.0, cw, cl, cr, X2;
        int gi = 0, g2;

        sup = -INFINITY;
        for (i = 0; i < 5; i++)
        {
            cel_phi(&phi, &d1, &d2, K, cel_add(ud, cel_scale(cel_sub(vd, ud), 0.25 * i)), -1);
            sup = fmax(sup, phi.re);
        }
        sup -= ctx->cd.re;

        /* negligible: bounded on the cell */
        if (cel_cell_plan(&gi, &bestX, ctx, md, hd, 1) == 0.0)
        {
            negl = 0;
            cel_bound(M, K, ctx->c, m, h, 0.0, -1, ctx->bprec);
            mag_mul(M, M, R);
            mag_mul_2exp_si(M, M, 1);
            if (cel_below_tol(M, ctx->logtol))
            {
                mag_add(ctx->err, ctx->err, M);
                done = 1;
            }
        }

        if (!done)
        {
            /* the cell as a whole if this does not take more nodes than
               its two halves (which are not taken for negligible if the
               bound for the cell has just failed where the estimate did
               not: wide parameters) */
            cw = cel_cell_plan(&gi, &bestX, ctx, md, hd, 0);
            if (cw < INFINITY)
            {
                cl = cel_cell_plan(&g2, &X2, ctx, cel_sub(md, cel_scale(hd, 0.5)), cel_scale(hd, 0.5), negl);
                cr = cel_cell_plan(&g2, &X2, ctx, cel_add(md, cel_scale(hd, 0.5)), cel_scale(hd, 0.5), negl);
            }
            else
                cl = cr = 0.0;

            if (cw <= cl + cr)
            {
                acb_t val;
                mag_t e;

                acb_init(val);
                mag_init(e);

                ctx->budget -= _acb_calc_gl_steps[gi];
                ok = cel_cell_value(val, ctx, m, h, md, hd, gi, sup);

                if (ok)
                {
                    cel_bound(M, K, ctx->c, m, h, bestX, -1, ctx->bprec);
                    _acb_calc_gl_error_bound(e, M, R, bestX, _acb_calc_gl_steps[gi]);
                    mag_add(ctx->err, ctx->err, e);
                    acb_add(ctx->sum, ctx->sum, val, ctx->prec + 30);
                    ok = mag_is_finite(e);
                }
                done = 1;

                acb_clear(val);
                mag_clear(e);
            }
        }
    }

    if (!done)
    {
        ok = cel_cell(ctx, u, m, depth + 1) && cel_cell(ctx, m, v, depth + 1);
    }

    acb_clear(m);
    acb_clear(h);
    mag_clear(R);
    mag_clear(M);
    return ok;
}

/* res = the integral of exp(k + lambda t + sum_j mu_j Log(al_j + be_j t))
   along the polygon with the np vertices pts (exact numbers), where
   the term end0 vanishes at the first vertex and the term end1 at the
   last one (or -1).  c is a reference value for the exponent (an exact
   number close to its largest value on the polygon), and logtol the
   logarithm of the tolerance for a cell relative to |exp(c)|.  Returns
   0 if the integral is not computed. */
int
_acb_hypgeom_exp_log_path(acb_t res, const acb_t lambda, acb_srcptr mu, acb_srcptr al,
    acb_srcptr be, const acb_t konst, const acb_t g1, int m, int end0, int end1, acb_srcptr pts,
    slong np, const acb_t c, double logtol, slong prec)
{
    cel_kernel K;
    cel_ctx ctx;
    int j, ok = 1;
    slong k;
    acb_t t;

    if (m > CEL_MAXTERMS || np < 2 || !acb_is_finite(lambda) || !acb_is_finite(konst)
        || !acb_is_exact(c) || !acb_is_finite(c))
        return 0;

    K.m = m;
    K.lam = lambda;
    K.konst = konst;
    K.mu = mu;
    K.al = al;
    K.be = be;
    K.dlam = cel_get(lambda);
    K.dkonst = cel_get(konst);
    K.end0 = end0;
    K.end1 = end1;
    K.E0 = pts;
    K.E1 = pts + np - 1;
    K.stable = (g1 != NULL);

    acb_init(t);

    for (j = 0; j < m && ok; j++)
    {
        if (!acb_is_finite(mu + j) || !acb_is_finite(al + j) || !acb_is_finite(be + j))
            ok = 0;
        K.dmu[j] = cel_get(mu + j);
        K.dal[j] = cel_get(al + j);
        K.dbe[j] = cel_get(be + j);

        /* a term at its end: zero there, exactly, and Re(mu_j) >= 0 */
        if (ok && (j == end0 || j == end1))
        {
            acb_mul(t, be + j, (j == end0) ? K.E0 : K.E1, 4096);
            acb_add(t, t, al + j, 4096);
            if (!acb_is_zero(t) || !arb_is_nonnegative(acb_realref(mu + j)))
                ok = 0;
        }
    }

    for (k = 0; k < np && ok; k++)
        if (!acb_is_exact(pts + k) || !acb_is_finite(pts + k))
            ok = 0;

    if (!ok)
    {
        acb_clear(t);
        return 0;
    }

    /* the stable form */
    if (K.stable)
    {
        K.g1 = cel_get(g1);
        K.phi0 = cel_mk(0.0, 0.0);
        for (j = 0; j < m; j++)
        {
            if (acb_is_zero(mu + j))
                K.rq[j] = cel_mk(0.0, 0.0);
            else if (cel_abs(K.dal[j]) == 0.0 || !(cel_abs(K.g1) < 1e300))
                K.stable = 0;
            else
                K.rq[j] = cel_div(K.dbe[j], K.dal[j]);
        }
    }

    ctx.K = &K;
    ctx.c = c;
    ctx.cd = cel_get(c);
    ctx.prec = prec;
    ctx.logtol = logtol;
    for (ctx.imax = 0; ctx.imax + 1 < ACB_CALC_GL_STEPS
        && _acb_calc_gl_steps[ctx.imax + 1] <= 0.5 * prec + 60.0; ctx.imax++) ;
    ctx.budget = 20000 + 200 * prec + 100 * np;

    /* the size of the terms of phi along the path */
    ctx.magd = cel_abs(ctx.cd) + cel_abs(K.dkonst) + 2.0;
    for (k = 0; k < np; k++)
    {
        cel_cd p = cel_get(pts + k), w;
        double s = cel_abs(cel_mul(K.dlam, p));

        for (j = 0; j < m; j++)
        {
            w = cel_add(K.dal[j], cel_mul(K.dbe[j], p));
            if (cel_abs(K.dmu[j]) != 0.0 && cel_abs(w) != 0.0)
                s += cel_abs(K.dmu[j]) * (fabs(log(cel_abs(w))) + 4.0);
        }
        if (s < 1e300)
            ctx.magd = fmax(ctx.magd, s + cel_abs(ctx.cd));
    }
    ctx.bprec = 96 + (slong) log2(ctx.magd);

    acb_init(ctx.kc);
    acb_init(ctx.sum);
    mag_init(ctx.err);
    ctx.scratch = flint_malloc(ACB_CALC_GL_DFLOAT_SCRATCH(4));

    acb_sub(ctx.kc, konst, c, prec + ctx.bprec);

    /* (in the stable form the values are relative to c) */
    if (K.stable)
    {
        acb_t z0;
        acb_init(z0);
        cel_phi_acb(t, NULL, NULL, &K, z0, -1, ctx.bprec + 30);
        acb_clear(z0);
        acb_sub(t, t, c, ctx.bprec + 30);
        K.phi0 = cel_get(t);
        ctx.cd = cel_mk(0.0, 0.0);
        if (!acb_is_finite(t))
            ok = 0;
    }

    for (k = 0; k + 1 < np && ok; k++)
    {
        if (acb_equal(pts + k, pts + k + 1))
            continue;
        ok = cel_cell(&ctx, pts + k, pts + k + 1, 0);
    }

    if (ok)
    {
        acb_add_error_mag(ctx.sum, ctx.err);
        /* (c can be large) */
        acb_exp(t, c, prec + 30 + ctx.bprec);
        acb_mul(res, ctx.sum, t, prec + 30);
        ok = acb_is_finite(res);
    }

    flint_free(ctx.scratch);
    acb_clear(ctx.kc);
    acb_clear(ctx.sum);
    mag_clear(ctx.err);
    acb_clear(t);
    return ok;
}

/* ------------------------------------------------------------------ */
/*
    The path.  For the integral over the segment [0, 1], or over a ray
    from 0 to infinity, of exp(phi) with the term end0 vanishing at 0
    (and the term end1 at 1), a polygon between the same end points is
    chosen heuristically in double precision:

    - The saddle points are the roots of phi' (of a polynomial of
      degree two at most here).  From each saddle point the two lines
      of steepest descent of Re(phi) are followed numerically until
      they have run into an end point, into another zero of the
      integrand (a point where a term with Re(mu_j) > 0 vanishes) or
      out to infinity in a direction of decay, far enough for the
      integrand to be negligible.  An end point where the integrand
      does not vanish (mu_j = 0) is itself the start of a line of
      descent.
    - The saddle points are the edges of a graph on these valleys; the
      path follows a chain of edges from one end to the other, with
      connections where the integrand is negligible in the valleys
      between.
    - The points of the lines are merged into longer chords where the
      integrand stays below its values at the ends of the chord and
      does not oscillate.

    Frames.  With parameters of the size N the saddle points have the
    width N^(-1/2), and phi is a sum of terms of the size N: neither
    the path nor the exponent on it can be resolved in doubles in the
    original coordinates.  Each saddle point (and each end point from
    which a line starts) gets a frame: the kernel translated to the
    base point (an exact number close to the saddle point) in the
    stable form

        phi(base + t) - phi(base) = g1 t + sum_j mu_j (log(1 + zeta_j) - zeta_j),
        zeta_j = be_j t / (al_j + be_j base),

    with g1 = phi'(base) computed in ball arithmetic.  The lines from
    the base point are followed in the coordinate t, and the piece of
    the path made of them is integrated with the translated kernel
    (al_j + be_j base for al_j, computed exactly, and k + lambda base
    for k).  The pieces are joined at vertices which are exact numbers
    in both frames.  The height of a route (for the choice between
    routes, and for what is negligible) is compared across frames with
    phi(base) in ball arithmetic.

    The logarithms are followed continuously along the lines (the
    lines can cross the principal cuts).  For the quadrature, each
    term gets a cut which avoids the path and the original contour:
    with c a Gaussian integer,

        log(al_j + be_j t) = Log(c (al_j + be_j t)) + gamma,

    for the branch of the logarithm on the left which is the principal
    one on the original contour, on the plane cut along the ray on
    which c (al_j + be_j t) is negative, where gamma is the constant
    Log(w) - Log(c w) for a point w of the contour.  If these rays
    avoid the polygon and the contour, the integrand with the
    logarithms on the right is analytic in a simply connected domain
    containing both and equal to the original integrand on the
    contour, so that the two integrals agree.  That the rays avoid the
    polygon is verified by the quadrature on every cell, and that they
    avoid the contour here, in ball arithmetic; nothing else about the
    path needs to be right.
*/

#define CEL_TRMAX 300
#define CEL_PMAX 1400
#define CEL_PI 3.14159265358979323846

/* the valleys (bits: a line can lead to several, through the region
   where the integrand is negligible) */
#define CEL_LE0 1
#define CEL_LE1 2
#define CEL_LINF 4
#define CEL_LV(j) (8 << (j))    /* the zero of the term j */

typedef struct
{
    cel_cd t[CEL_TRMAX];
    cel_cd phi[CEL_TRMAX];   /* followed continuously */
    int n;
    int label;
    int eidx[2];             /* the point from which the chord to an end
                                point is taken */
}
cel_trace;

typedef struct
{
    const cel_kernel * K;      /* with the principal logarithms */
    cel_cd tp[CEL_MAXTERMS];   /* the singular points */
    int sing[CEL_MAXTERMS];    /* the term has one */
    int infinite;
    int inf_valley;
    cel_cd E0, E1;
    cel_cd d0;                 /* the direction of the ray */
    double D;                  /* the integrand is negligible this far
                                  below a saddle point */
    cel_cd tilt[CEL_MAXTERMS]; /* the singular points are moved by this
                                  much for the planning (see cel_integrate) */
}
cel_plan;

static inline cel_cd cel_conj(cel_cd a) { return cel_mk(a.re, -a.im); }
static inline int cel_is0(cel_cd a) { return a.re == 0.0 && a.im == 0.0; }

/* The planning is done in the stable form of the kernel, with the
   origin at a saddle point (or at an end point or a point of the
   contour).  Along a line, the functions

       L_j(t) = log(1 + zeta_j) - zeta_j,   zeta_j = rq_j t,

   are followed continuously (the lines can go around the singular
   points): from t0 to t1 along the segment, with x = (zeta1 - zeta0) / (1 + zeta0),

       L_j(t1) - L_j(t0) = (log(1 + x) - x) - (zeta1 - zeta0) zeta0 / (1 + zeta0),

   and phi(t) - phi(0) = g1 t + sum_j mu_j L_j(t). */
static void
cel_logs_step(cel_cd * L, const cel_kernel * K, cel_cd t0, cel_cd t1)
{
    int j;

    for (j = 0; j < K->m; j++)
    {
        cel_cd z0, dz, u0;

        if (cel_is0(K->dmu[j]))
            continue;
        z0 = cel_mul(K->rq[j], t0);
        dz = cel_mul(K->rq[j], cel_sub(t1, t0));
        u0 = cel_mk(1.0 + z0.re, z0.im);
        if (cel_is0(u0))
            continue;
        L[j] = cel_add(L[j], cel_sub(cel_l1pmx(cel_div(dz, u0)), cel_div(cel_mul(dz, z0), u0)));
    }
}

static cel_cd
cel_phi_logs(const cel_kernel * K, cel_cd t, const cel_cd * L)
{
    cel_cd p = cel_mul(K->g1, t);
    int j;

    for (j = 0; j < K->m; j++)
        if (!cel_is0(K->dmu[j]))
            p = cel_add(p, cel_mul(K->dmu[j], L[j]));
    return p;
}

/* the value at x, continued along the segment from t with the logarithms L */
static cel_cd
cel_phi_from(const cel_kernel * K, cel_cd t, const cel_cd * L, cel_cd x)
{
    cel_cd L2[CEL_MAXTERMS];
    int j;

    for (j = 0; j < K->m; j++)
        L2[j] = L[j];
    cel_logs_step(L2, K, t, x);
    return cel_phi_logs(K, x, L2);
}

/* the distance from t to the nearest singular point */
static double
cel_sing_dist(const cel_plan * P, cel_cd t)
{
    double d = INFINITY;
    int j;

    for (j = 0; j < P->K->m; j++)
        if (P->sing[j])
            d = fmin(d, cel_abs(cel_sub(t, P->tp[j])));
    return d;
}

#define CEL_CHORD_SAMPLES 16

/* Does Re(phi) descend on the chord from t (with the logarithms L and
   the value phit) to the zero E of a term, or stay below lim? */
static int
cel_chord_descends(const cel_plan * P, cel_cd t, const cel_cd * L, cel_cd phit, cel_cd E, int je, double lim)
{
    /* (Samples at the scales of the chord, and at smaller ones: close
       to a saddle point with a large phi'' a chord which leaves the
       line of descent at an angle of more than pi / 4 rises over a
       short part of its length, by more than the line has descended.) */
    static const double fr[CEL_CHORD_SAMPLES] = { 0.2 / 1024, 0.2 / 512, 0.2 / 256, 0.2 / 128,
        0.2 / 64, 0.2 / 32, 0.2 / 16, 0.2 / 8, 0.2 / 4, 0.2 / 2, 0.15, 0.2, 0.4, 0.6, 0.8, 0.95 };
    double prev = phit.re, rhump = 0.0, hump = 0.0;
    int i, j;

    /* Close to the end point, with tau = t - E = r e, |e| = 1:
       Re(phi) = Re(mu log(be tau)) + Re(g(E)) + Re(g' e) r + ..., with g
       the exponent without the term which vanishes; this has a maximum
       at r = -Re(mu) / Re(g' e) (where another saddle point lies close
       to the end point), with the value
       Re(g(E)) + Re(mu (log(be e r) - 1)): one more sample.  In the
       stable form, g(E) - phi(0) = p + mu - mu log(al), with p the value
       of the form without the logarithm of the term. */
    {
        cel_cd pd, d1, d2, e = cel_sub(t, E), mu = P->K->dmu[je];
        double len = cel_abs(e), re, r;

        cel_phi(&pd, &d1, &d2, P->K, E, je);
        e = cel_scale(e, 1.0 / len);
        re = d1.re * e.re - d1.im * e.im;
        if (re < 0.0 && mu.re > 0.0)
        {
            r = -mu.re / re;
            if (r < 0.05 * len)
            {
                hump = pd.re + cel_mul(mu, cel_log(cel_scale(cel_mul(P->K->rq[je], e), r))).re;
                rhump = 1.0;
            }
        }
    }

    /* (not close by another singular point) */
    for (j = 0; j < P->K->m; j++)
    {
        cel_cd v = cel_sub(P->tp[j], t), d = cel_sub(E, t);
        double x, dist;

        if (!P->sing[j] || cel_abs(cel_sub(P->tp[j], E)) <= 1e-12 * (1.0 + cel_abs(E)))
            continue;
        x = (v.re * d.re + v.im * d.im) / (d.re * d.re + d.im * d.im);
        x = fmin(fmax(x, 0.0), 1.0);
        dist = cel_abs(cel_sub(v, cel_scale(d, x)));
        if (dist < 0.3 * cel_abs(v))
            return 0;
    }

    for (i = 0; i < CEL_CHORD_SAMPLES + (rhump != 0.0); i++)
    {
        cel_cd x = cel_add(t, cel_scale(cel_sub(E, t), fr[i < CEL_CHORD_SAMPLES ? i : 0]));
        cel_cd p = cel_phi_from(P->K, t, L, x);
        if (i == CEL_CHORD_SAMPLES)
            p.re = hump;
        if (!(p.re <= prev + 0.05) && !(p.re <= lim))
            return 0;
        prev = fmin(prev, p.re);
    }
    return 1;
}

/* The line of steepest descent from the origin of the frame: a saddle
   point with the first step given, or an ordinary point (step = 0). */
static void
cel_trace_run(cel_trace * tr, const cel_plan * P, cel_cd step)
{
    const cel_kernel * K = P->K;
    cel_cd L[CEL_MAXTERMS], phi, d1, d2, g, tn, pn, dummy, t = cel_mk(0.0, 0.0);
    double h, hprev, dist, top;
    int j, it, ee, negl, nneg = 0;

    for (j = 0; j < K->m; j++)
        L[j] = cel_mk(0.0, 0.0);

    tr->n = 0;
    tr->label = 0;
    tr->eidx[0] = tr->eidx[1] = 0;
    phi = cel_phi_logs(K, t, L);
    top = phi.re;
    tr->t[tr->n] = t;
    tr->phi[tr->n++] = phi;
    hprev = cel_abs(step);

    for (it = 0; it < CEL_TRMAX - 3; it++)
    {
        cel_phi(&dummy, &d1, &d2, K, t, -1);
        dist = cel_sing_dist(P, t);

        if (it == 0 && hprev != 0.0)
        {
            tn = cel_add(t, step);
        }
        else
        {
            if (cel_is0(d1))
                return;
            h = 0.3 * dist;
            /* (where the integrand is negligible: steps of the size of
               the distance to the singular points, to reach the valley
               from a narrow saddle) */
            if (!(top - phi.re >= P->D))
            {
                if (!cel_is0(d2))
                    h = fmin(h, 0.3 * cel_abs(d1) / cel_abs(d2));
                if (hprev != 0.0)
                    h = fmin(h, 2.0 * hprev);
                else
                    h = fmin(h, 1.0 / cel_abs(d1));
            }
            else if (!(h < 1e300))
                h = 8.0 * hprev;
            if (!(h > 1e-14 * cel_abs(t)) || !(h < 1e300))
                return;

            /* (the midpoint rule) */
            g = cel_scale(cel_conj(d1), -1.0 / cel_abs(d1));
            tn = cel_add(t, cel_scale(g, 0.5 * h));
            cel_phi(&dummy, &d1, &d2, K, tn, -1);
            if (cel_is0(d1))
                return;
            g = cel_scale(cel_conj(d1), -1.0 / cel_abs(d1));
            tn = cel_add(t, cel_scale(g, h));
            hprev = h;
        }

        cel_logs_step(L, K, t, tn);
        pn = cel_phi_logs(K, tn, L);
        t = tn;
        phi = pn;
        tr->t[tr->n] = t;
        tr->phi[tr->n++] = phi;

        if (!(phi.re < 1e300) || !(cel_abs(t) < 1e300))
            return;

        cel_phi(&dummy, &d1, &d2, K, t, -1);
        g = cel_scale(cel_conj(d1), -1.0 / (cel_abs(d1) + 1e-300));

        /* into an end point: its term dominates the gradient, or the
           chord to it descends */
        /* (the nearer one first) */
        for (ee = 0; ee < 2; ee++)
        {
            int e0 = (!P->infinite && cel_abs(cel_sub(t, P->E1)) < cel_abs(cel_sub(t, P->E0)));
            int e = e0 ? 1 - ee : ee;
            int je = (e == 0) ? K->end0 : K->end1;
            cel_cd E = (e == 0) ? P->E0 : P->E1, w, v, to;
            double rest, resta, own, al;

            if (je < 0 || (e == 1 && P->infinite) || cel_is0(K->dmu[je]))
                continue;

            w = cel_add(K->dal[je], cel_mul(K->dbe[je], t));
            v = cel_div(cel_mul(K->dmu[je], K->dbe[je]), w);
            rest = cel_abs(cel_sub(d1, v));
            own = K->dmu[je].re * cel_abs(K->dbe[je]) / cel_abs(w);
            /* (the other terms, which can cancel) */
            resta = cel_abs(K->dlam);
            for (j = 0; j < K->m; j++)
                if (j != je && P->sing[j])
                    resta += cel_abs(K->dmu[j]) * cel_abs(K->dbe[j])
                        / cel_abs(cel_add(K->dal[j], cel_mul(K->dbe[j], t)));
            to = cel_sub(E, t);

            /* (the line heads for the end point; or around it, where
               the exponent is close to imaginary) */
            al = (g.re * to.re + g.im * to.im) / cel_abs(to);
            negl = (top - phi.re >= P->D);
            (void) rest;
            if (tr->label & (e == 0 ? CEL_LE0 : CEL_LE1))
                continue;
            if (own >= 4.0 * resta
                || ((al > 0.5 || (al > -0.05 && cel_abs(v) >= 10.0 * resta))
                    && cel_chord_descends(P, t, L, phi, E, je, -INFINITY)))
            {
                tr->label |= (e == 0) ? CEL_LE0 : CEL_LE1;
                tr->eidx[e] = tr->n - 1;
                if (!negl)
                    return;
            }
            else if (negl && al > 0.0 && cel_chord_descends(P, t, L, phi, E, je, top - P->D))
            {
                /* (through the negligible; the chord starts downhill, not
                   back over the saddle point) */
                tr->label |= (e == 0) ? CEL_LE0 : CEL_LE1;
                tr->eidx[e] = tr->n - 1;
            }
        }

        /* (enough: both end points, or some steps into the negligible) */
        if (!P->infinite && (tr->label & CEL_LE0) && (tr->label & CEL_LE1))
            return;
        if (top - phi.re >= P->D && ++nneg > 100)
            return;

        if (top - phi.re >= P->D)
        {
            /* into another zero */
            for (j = 0; j < K->m; j++)
            {
                cel_cd w, v;
                (void) v;

                if (!P->sing[j] || j == K->end0 || j == K->end1 || !(K->dmu[j].re > 0.0))
                    continue;
                double resta = cel_abs(K->dlam);
                int i2;

                w = cel_add(K->dal[j], cel_mul(K->dbe[j], t));
                v = cel_div(cel_mul(K->dmu[j], K->dbe[j]), w);
                for (i2 = 0; i2 < K->m; i2++)
                    if (i2 != j && P->sing[i2])
                        resta += cel_abs(K->dmu[i2]) * cel_abs(K->dbe[i2])
                            / cel_abs(cel_add(K->dal[i2], cel_mul(K->dbe[i2], t)));
                /* (and no other singular point lies nearer, in the way) */
                if (K->dmu[j].re * cel_abs(K->dbe[j]) / cel_abs(w) >= 4.0 * resta
                    && cel_abs(cel_sub(t, P->tp[j])) <= 1.0001 * cel_sing_dist(P, t))
                {
                    tr->label |= CEL_LV(j);
                    return;
                }
            }

            /* out to infinity */
            if (P->inf_valley)
            {
                double s = 0.0, far = 0.0;

                for (j = 0; j < K->m; j++)
                {
                    if (cel_is0(K->dmu[j]) || !P->sing[j])
                        continue;
                    s += cel_abs(K->dmu[j]) * cel_abs(K->dbe[j])
                        / cel_abs(cel_add(K->dal[j], cel_mul(K->dbe[j], t)));
                    far = fmax(far, cel_abs(P->tp[j]));
                }

                if (cel_is0(K->dlam) ? (cel_abs(t) >= 4.0 * far) : (cel_abs(K->dlam) >= 4.0 * s))
                {
                    tr->label |= CEL_LINF;
                    return;
                }
            }
        }
    }
}

/* does the ray from o in the direction d meet the segment [a, b]? */
static int
cel_ray_meets(cel_cd o, cel_cd d, cel_cd a, cel_cd b)
{
    cel_cd e = cel_sub(b, a), f = cel_sub(a, o);
    double den = d.re * e.im - d.im * e.re, x, y, tol;

    tol = 1e-9 * (cel_abs(d) * cel_abs(e) + 1e-300);
    if (fabs(den) <= tol)
    {
        /* parallel: on the same line? */
        double cr = f.re * d.im - f.im * d.re;
        if (fabs(cr) > 1e-9 * cel_abs(f) * cel_abs(d))
            return 0;
        return (f.re * d.re + f.im * d.im) > 0.0
            || ((b.re - o.re) * d.re + (b.im - o.im) * d.im) > 0.0;
    }
    /* o + x d = a + y e */
    x = (f.re * e.im - f.im * e.re) / den;
    y = (f.re * d.im - f.im * d.re) / den;
    return x > 1e-9 && y >= -1e-9 && y <= 1.0 + 1e-9;
}

/* the Gaussian integers c for the cuts, in the order of preference */
static const int cel_rot[8][2] = { { 1, 0 }, { 0, 1 }, { 0, -1 }, { 1, 1 }, { 1, -1 },
    { -1, 1 }, { -1, -1 }, { -1, 0 } };

/* Rigorous: does al + be t stay off the negative numbers and zero for
   t = t0 + s d, 0 <= s <= 1 (0 <= s if ray)?  If from_zero, al + be t0
   is exactly zero, and s = 0 is excluded; if to_zero, al + be (t0 + d)
   is exactly zero, and s = 1 is excluded. */
static int
cel_avoids_cut(const acb_t al, const acb_t be, const acb_t t0, const acb_t d, int ray,
    int from_zero, int to_zero, slong prec)
{
    acb_t w0, wd, w1;
    arb_t s, x;
    int ok = 0;

    acb_init(w0); acb_init(wd); acb_init(w1);
    arb_init(s); arb_init(x);

    acb_mul(wd, be, d, prec);
    if (from_zero)
    {
        /* s wd, s > 0 */
        ok = arb_is_positive(acb_realref(wd)) || !arb_contains_zero(acb_imagref(wd));
    }
    else
    {
        acb_mul(w0, be, t0, prec);
        acb_add(w0, w0, al, prec);
        acb_add(w1, w0, wd, prec);

        /* (1 - s) w0, s < 1 */
        if (to_zero)
            ok = arb_is_positive(acb_realref(w0)) || !arb_contains_zero(acb_imagref(w0));
        else if (arb_is_positive(acb_imagref(w0))
            && (ray ? arb_is_nonnegative(acb_imagref(wd)) : arb_is_positive(acb_imagref(w1))))
            ok = 1;
        else if (arb_is_negative(acb_imagref(w0))
            && (ray ? arb_is_nonpositive(acb_imagref(wd)) : arb_is_negative(acb_imagref(w1))))
            ok = 1;
        else if (arb_is_positive(acb_realref(w0))
            && (ray ? arb_is_nonnegative(acb_realref(wd)) : arb_is_positive(acb_realref(w1))))
            ok = 1;
        else if (!arb_contains_zero(acb_imagref(wd)))
        {
            /* the crossing of the real line, at s = -Im(w0) / Im(wd) */
            arb_div(s, acb_imagref(w0), acb_imagref(wd), prec);
            arb_neg(s, s);
            arb_mul(x, s, acb_realref(wd), prec);
            arb_add(x, x, acb_realref(w0), prec);
            ok = arb_is_positive(x);
        }
    }

    acb_clear(w0); acb_clear(wd); acb_clear(w1);
    arb_clear(s); arb_clear(x);
    return ok;
}

/* A frame: the kernel translated to a base point (a saddle point, an
   end point where the integrand does not vanish, or a point of the
   contour), t = base + delta, in the stable form in the variable
   delta.  The lines from the base point and the part of the path
   around it are planned in these coordinates, which resolve the
   neighbourhood of a saddle point whatever the size of the
   parameters. */
#define CEL_NFRAMES 5
#define CEL_PIECEMAX 700

typedef struct
{
    acb_t base;               /* an exact number */
    acb_ptr al;               /* al_j + be_j base */
    acb_t g1;                 /* phi'(base) */
    cel_kernel K;
    cel_plan P;
    cel_cd babs;              /* the base point as a double */
}
cel_frame;

/* A piece of the path: the vertices around one base point, in the
   coordinates of its frame */
typedef struct
{
    int f;
    cel_cd t[CEL_PIECEMAX];
    cel_cd phi[CEL_PIECEMAX];             /* phi - phi(base), followed continuously */
    cel_cd L[CEL_PIECEMAX][CEL_MAXTERMS];
    char fix[CEL_PIECEMAX];   /* the vertex is not merged away */
    int n;
    int sidx;                 /* the index of the base point */
    int saddle;               /* which is a saddle point */
    int junction;             /* the first vertex is the last one of the
                                 piece before */
    int e0, e1;               /* the first (last) vertex is an end point of
                                 the integral */
}
cel_piece;

static int
cel_piece_push(cel_piece * pc, cel_cd t)
{
    if (pc->n >= CEL_PIECEMAX)
        return 0;
    if (pc->n > 0 && pc->t[pc->n - 1].re == t.re && pc->t[pc->n - 1].im == t.im)
        return 1;
    pc->fix[pc->n] = 0;
    pc->t[pc->n++] = t;
    return 1;
}

/* the frame at base (which is copied), for the kernel K0 with the
   settings P0; returns 0 if a term vanishes at the base point */
static int
cel_frame_init(cel_frame * F, const cel_kernel * K0, const cel_plan * P0, const acb_t base, slong wp)
{
    acb_t t;
    int j, ok = 1, m = K0->m;

    acb_init(t);
    acb_init(F->base);
    acb_init(F->g1);
    F->al = _acb_vec_init(m);

    acb_set(F->base, base);
    F->babs = cel_get(base);
    acb_set(F->g1, K0->lam);

    F->K = *K0;
    F->K.al = F->al;
    F->K.stable = 1;
    F->K.phi0 = cel_mk(0.0, 0.0);
    F->P = *P0;
    F->P.K = &F->K;

    for (j = 0; j < m; j++)
    {
        acb_mul(F->al + j, K0->be + j, base, wp);
        acb_add(F->al + j, F->al + j, K0->al + j, wp);
        F->K.dal[j] = cel_get(F->al + j);
        F->K.rq[j] = cel_mk(0.0, 0.0);
        F->P.tp[j] = cel_mk(0.0, 0.0);

        if (cel_is0(K0->dmu[j]))
            continue;
        if (cel_is0(F->K.dal[j]) || !acb_is_finite(F->al + j))
        {
            ok = 0;
            continue;
        }
        acb_div(t, K0->be + j, F->al + j, wp);
        acb_addmul(F->g1, t, K0->mu + j, wp);
        F->K.rq[j] = cel_div(K0->dbe[j], F->K.dal[j]);
        if (P0->sing[j])
            F->P.tp[j] = cel_scale(cel_div(cel_mk(1.0, 0.0), F->K.rq[j]), -1.0);
    }
    F->K.g1 = cel_get(F->g1);

    /* (the tilted kernel: the singular point of the term j at
       tp_j + tilt_j, which changes rq_j and phi'(base) = g1) */
    for (j = 0; j < m; j++)
    {
        if (P0->sing[j] && !cel_is0(P0->tilt[j]) && !cel_is0(F->P.tp[j]))
        {
            cel_cd rq = cel_scale(cel_div(cel_mk(1.0, 0.0), cel_add(F->P.tp[j], P0->tilt[j])), -1.0);

            F->K.g1 = cel_add(F->K.g1, cel_mul(K0->dmu[j], cel_sub(rq, F->K.rq[j])));
            F->K.rq[j] = rq;
            F->P.tp[j] = cel_add(F->P.tp[j], P0->tilt[j]);
        }
    }

    /* the end points, in the coordinates of the frame */
    acb_neg(t, base);
    F->P.E0 = cel_get(t);
    acb_add_ui(t, t, 1, wp);
    F->P.E1 = cel_get(t);
    if (P0->sing[K0->end0])
        F->P.tp[K0->end0] = F->P.E0;
    if (K0->end1 >= 0 && P0->sing[K0->end1])
        F->P.tp[K0->end1] = F->P.E1;

    ok = ok && acb_is_finite(F->g1);
    acb_clear(t);
    return ok;
}

/* the size of the geometry around the base point: the distances to
   the end points and the singular points */
static double
cel_frame_scale(const cel_frame * F, int infinite)
{
    double s = cel_abs(F->P.E0);
    int j;

    if (!infinite)
        s = fmax(s, cel_abs(F->P.E1));
    for (j = 0; j < F->K.m; j++)
        if (F->P.sing[j])
            s = fmax(s, cel_abs(F->P.tp[j]));
    return s;
}

static void
cel_frame_clear(cel_frame * F, int m)
{
    acb_clear(F->base);
    acb_clear(F->g1);
    _acb_vec_clear(F->al, m);
}

/* the values of phi - phi(base) and the functions L_j at the vertices
   of a piece, followed from the base point in both directions */
static void
cel_piece_track(cel_piece * pc, const cel_frame * F, int sing0, int sing1)
{
    const cel_kernel * K = &F->K;
    cel_cd L[CEL_MAXTERMS];
    int j, k, dir;

    for (dir = 1; dir >= -1; dir -= 2)
    {
        for (j = 0; j < K->m; j++)
            L[j] = cel_mk(0.0, 0.0);

        for (k = pc->sidx; k >= 0 && k < pc->n; k += dir)
        {
            if ((k == 0 && pc->e0 && sing0) || (k == pc->n - 1 && pc->e1 && sing1))
            {
                pc->phi[k] = cel_mk(-INFINITY, 0.0);
                continue;
            }
            if (k != pc->sidx)
                cel_logs_step(L, K, pc->t[k - dir], pc->t[k]);
            for (j = 0; j < K->m; j++)
                pc->L[k][j] = L[j];
            pc->phi[k] = cel_phi_logs(K, pc->t[k], L);
        }
    }
}

/* The saddle points: the roots of
   lambda prod_j w_j + sum_j mu_j be_j prod_{i != j} w_i, w_j = al_j + be_j t,
   over the terms tix with a singular point (a polynomial of degree two
   at most), in ball arithmetic, polished by Newton's method for phi'
   and rounded to exact numbers.  Returns 0 if the degree is larger. */
static int
cel_saddles(acb_ptr sad, int * ns, const cel_kernel * K, const int * tix, int nterm, slong wp)
{
    acb_t q[4], r[4], poly[4], t, u, d1, d2;
    int i, j, a, b, deg, ok = 1, k;
    double big = 0.0;

    for (i = 0; i < 4; i++)
    {
        acb_init(q[i]);
        acb_init(r[i]);
        acb_init(poly[i]);
    }
    acb_init(t); acb_init(u); acb_init(d1); acb_init(d2);
    *ns = 0;

    /* q = prod w_j */
    acb_one(q[0]);
    deg = 0;
    for (a = 0; a < nterm; a++)
    {
        j = tix[a];
        acb_zero(q[deg + 1]);
        for (i = deg + 1; i >= 1; i--)
        {
            acb_mul(q[i], q[i], K->al + j, wp);
            acb_addmul(q[i], q[i - 1], K->be + j, wp);
        }
        acb_mul(q[0], q[0], K->al + j, wp);
        deg++;
    }
    for (i = 0; i <= deg; i++)
        acb_mul(poly[i], K->lam, q[i], wp);

    for (a = 0; a < nterm; a++)
    {
        int dr = 0;

        acb_mul(r[0], K->mu + tix[a], K->be + tix[a], wp);
        for (b = 0; b < nterm; b++)
        {
            if (b == a)
                continue;
            j = tix[b];
            acb_zero(r[dr + 1]);
            for (i = dr + 1; i >= 1; i--)
            {
                acb_mul(r[i], r[i], K->al + j, wp);
                acb_addmul(r[i], r[i - 1], K->be + j, wp);
            }
            acb_mul(r[0], r[0], K->al + j, wp);
            dr++;
        }
        for (i = 0; i <= dr; i++)
            acb_add(poly[i], poly[i], r[i], wp);
    }

    /* The saddle points of the midpoints of the parameters: the path
       is a matter of heuristics, and the square root of a ball is not
       close to that of its midpoint near the cut (for conjugate saddle
       points of real parameters with complex radii). */
    for (i = 0; i <= deg; i++)
        acb_get_mid(poly[i], poly[i]);

    /* the degree */
    for (i = 0; i <= deg; i++)
        big = fmax(big, cel_abs(cel_get(poly[i])));
    while (deg > 0 && cel_abs(cel_get(poly[deg])) <= 1e-13 * big)
        deg--;

    if (deg > 2)
        ok = 0;
    else if (deg == 2)
    {
        /* q = -(b + sgn sqrt(b^2 - 4 a c)) / 2, the roots q / a and c / q */
        cel_cd bd, sd;

        acb_mul(t, poly[1], poly[1], wp);
        acb_mul(u, poly[2], poly[0], wp);
        acb_mul_2exp_si(u, u, 2);
        acb_sub(t, t, u, wp);
        acb_sqrt(t, t, wp);
        bd = cel_get(poly[1]);
        sd = cel_get(t);
        /* (scaled: the product could overflow) */
        if ((bd.re / (cel_abs(bd) + 1e-300)) * (sd.re / (cel_abs(sd) + 1e-300))
            + (bd.im / (cel_abs(bd) + 1e-300)) * (sd.im / (cel_abs(sd) + 1e-300)) < 0.0)
            acb_neg(t, t);
        acb_add(t, t, poly[1], wp);
        acb_mul_2exp_si(t, t, -1);
        acb_neg(t, t);
        if (acb_contains_zero(t))
        {
            acb_zero(sad + (*ns)++);
        }
        else
        {
            acb_div(sad + (*ns)++, t, poly[2], wp);
            acb_div(sad + (*ns)++, poly[0], t, wp);
        }
    }
    else if (deg == 1)
    {
        acb_div(sad, poly[0], poly[1], wp);
        acb_neg(sad, sad);
        *ns = 1;
    }

    for (k = 0; k < *ns && ok; k++)
    {
        acb_get_mid(sad + k, sad + k);
        for (i = 0; i < 3; i++)
        {
            cel_phi_acb(t, d1, d2, K, sad + k, -1, wp);
            acb_get_mid(d1, d1);
            acb_get_mid(d2, d2);
            if (acb_contains_zero(d2) || !acb_is_finite(d1))
                break;
            acb_div(t, d1, d2, wp);
            /* (a step which is small compared to the distance to the
               singular points) */
            {
                double st = cel_abs(cel_get(t)), di = INFINITY;
                for (a = 0; a < nterm; a++)
                {
                    acb_mul(u, K->be + tix[a], sad + k, wp);
                    acb_add(u, u, K->al + tix[a], wp);
                    di = fmin(di, cel_abs(cel_get(u)) / cel_abs(K->dbe[tix[a]]));
                }
                if (!(st < 0.1 * di))
                    break;
            }
            acb_sub(sad + k, sad + k, t, wp);
            acb_get_mid(sad + k, sad + k);
        }
        if (!acb_is_finite(sad + k))
            ok = 0;
    }

    for (i = 0; i < 4; i++)
    {
        acb_clear(q[i]);
        acb_clear(r[i]);
        acb_clear(poly[i]);
    }
    acb_clear(t); acb_clear(u); acb_clear(d1); acb_clear(d2);
    return ok;
}

/* Points on the chord from p to the end point E, for a narrow saddle
   point far from the end point: the point from where the integrand is
   negligible, and points beyond it at distances s from the base point
   which grow like (s^2 d)^(1/3), d the distance to the singular points.
   (On a cell [s, s'] of the Gaussian part the bound on the chord is
   sharp up to the terms of the third order, mu (s' / d)^3, which are to
   stay below the drop mu (s / d)^2 at its near end; with the cells of a
   bisection there would be log2 of the ratio of the scales.)  Returns
   the number of points (0 if the chord needs none). */
static int
cel_chord_split(cel_cd * x, int max, const cel_frame * Fr, cel_cd p, cel_cd E)
{
    cel_cd phi, d1, d2, v = cel_sub(E, p);
    double len = cel_abs(v), f, s, d, sn;
    int k, n = 0;

    if (!(len > 0.0))
        return 0;
    f = fmax(cel_abs(p), 1e-300) / len;
    for (k = 0; k < 60 && f < 0.25; k++, f *= 2.0)
    {
        x[0] = cel_add(p, cel_scale(v, f));
        cel_phi(&phi, &d1, &d2, &Fr->K, x[0], -1);
        if (phi.re <= -Fr->P.D - 3.0)
        {
            n = 1;
            break;
        }
    }
    if (n == 0)
        return 0;

    d = cel_sing_dist(&Fr->P, cel_mk(0.0, 0.0));
    s = cel_abs(x[0]);
    while (n < max)
    {
        sn = 0.5 * cbrt(s * s * d);
        if (!(sn > 2.0 * s))
            break;
        f = (sn - cel_abs(p)) / len;
        if (!(f < 0.25))
            break;
        x[n++] = cel_add(p, cel_scale(v, f));
        s = sn;
    }
    return n;
}

/* The chords of a piece: its vertices are merged where the chord stays
   on the same side of the singular points as the line, and below its
   ends without oscillating (or where everything is below neg, the
   level of the negligible, relative to the base point).  keep gets
   the indices of the vertices which remain; returns their number. */
static int
cel_piece_chords(int * keep, const cel_piece * pc, const cel_frame * Fr, double neg)
{
    const cel_kernel * K = &Fr->K;
    int a = 0, b, lastv = pc->n - 1, nk = 0, j;

    keep[nk++] = 0;
    while (a < lastv)
    {
        b = a + 1;
        /* (the chords at the ends, from the zeros, are kept) */
        while (b < lastv && pc->phi[a].re > -INFINITY && pc->phi[b + 1].re > -INFINITY && !pc->fix[b])
        {
            int good = 1, s;
            cel_cd pa = pc->t[a], pb = pc->t[b + 1];
            double hi = fmax(pc->phi[a].re, pc->phi[b + 1].re);

            for (j = 0; j < K->m && good; j++)
            {
                cel_cd za, zb, ua, ub;
                if (!Fr->P.sing[j])
                    continue;
                za = cel_mul(K->rq[j], pa);
                zb = cel_mul(K->rq[j], pb);
                ua = cel_mk(1.0 + za.re, za.im);
                ub = cel_mk(1.0 + zb.re, zb.im);
                /* (the change of the argument along the line, and along the chord) */
                if (fabs((pc->L[b + 1][j].im + zb.im) - (pc->L[a][j].im + za.im)
                    - cel_log(cel_div(ub, ua)).im) > 1e-6)
                    good = 0;
                /* (not close to a singular point) */
                if (cel_abs(cel_sub(pb, pa)) > 1.2 * cel_abs(cel_sub(cel_scale(cel_add(pa, pb), 0.5), Fr->P.tp[j])))
                    good = 0;
            }

            for (s = 1; s < 4 && good; s++)
            {
                cel_cd x = cel_add(pa, cel_scale(cel_sub(pb, pa), 0.25 * s));
                cel_cd p = cel_phi_from(K, pa, pc->L[a], x);
                double lin = pc->phi[a].im + 0.25 * s * (pc->phi[b + 1].im - pc->phi[a].im);

                if (p.re <= neg && hi <= neg)
                    continue;
                if (!(p.re <= hi + 0.3) || fabs(p.im - lin) > 1.0
                    || fabs(pc->phi[b + 1].im - pc->phi[a].im) > 2.0)
                    good = 0;
            }

            if (!good)
                break;
            b++;
        }
        keep[nk++] = b;
        a = b;
    }
    return nk;
}

/* the vertex k of a piece as an exact number, in the coordinates of
   its frame */
static void
cel_piece_vertex(acb_t v, const cel_piece * pc, int k, const cel_frame * Fr, const acb_t J, int infinite)
{
    if (k == 0 && pc->junction)
        acb_set(v, J);
    else if (k == 0 && pc->e0)
        acb_neg(v, Fr->base);
    else if (k == pc->n - 1 && pc->e1 && !infinite)
    {
        acb_neg(v, Fr->base);
        acb_add_ui(v, v, 1, 16384);
    }
    else
        acb_set_d_d(v, pc->t[k].re, pc->t[k].im);
}

static int
cel_integrate(acb_t res, const acb_t dir, const acb_t lambda, acb_srcptr mu,
    acb_srcptr al, acb_srcptr be, int m, int end0, int end1, int infinite, double logscale,
    slong prec, int level, const acb_struct * via, int tilted)
{
    cel_kernel K;
    cel_plan P0;
    cel_frame F[CEL_NFRAMES];
    cel_trace * tr = NULL;
    cel_piece * pcs = NULL, * bpc = NULL;
    cel_cd tpabs[CEL_MAXTERMS], * ab = NULL, dummy, d1, d2;
    int nf = 0, ntr = 0, npc = 0, nbest = 0, ns = 0, fdirect = -1;
    int j, k, i, ok = 1, nterm = 0, tix[CEL_MAXTERMS], sing0, sing1;
    int elab[4][2], etr[4][2], esad[4], efr[4], ne = 0;   /* the edges */
    int cand[9][3], cdir[9][3], clab[9][3], clen[9], ncand = 0, ci, have = 0, target, direct = 0, bdirect = 0;
    acb_t zero, one, konst, bkonst, c, t, u, total;
    acb_ptr ral, rbe, bal, bbe, sad, Jv, bJ, H, bH, tal;
    arb_t score, bscore, vmax, x;
    slong wp0, redo = 0, bits = 0;
    double logI = 0.0, toplen = 1.0;

    if (m > CEL_MAXTERMS || end0 < 0 || end0 >= m || (!infinite && (end1 < 0 || end1 >= m || end1 == end0))
        || !acb_is_finite(lambda))
        return 0;
    for (j = 0; j < m; j++)
        if (!acb_is_finite(mu + j) || !acb_is_finite(al + j) || !acb_is_finite(be + j))
            return 0;
    if (!acb_is_zero(al + end0) || !arb_is_nonnegative(acb_realref(mu + end0)))
        return 0;

    acb_init(zero); acb_init(one); acb_init(konst); acb_init(bkonst); acb_init(c); acb_init(t);
    acb_init(u); acb_init(total);
    arb_init(score); arb_init(bscore); arb_init(vmax); arb_init(x);
    ral = _acb_vec_init(m);
    rbe = _acb_vec_init(m);
    bal = _acb_vec_init(m);
    bbe = _acb_vec_init(m);
    tal = _acb_vec_init(m);
    sad = _acb_vec_init(2);
    Jv = _acb_vec_init(3);
    bJ = _acb_vec_init(3);
    H = _acb_vec_init(3);
    bH = _acb_vec_init(3);
    acb_one(one);

    if (!infinite)
    {
        acb_add(t, al + end1, be + end1, 4096);
        if (!acb_is_zero(t) || !arb_is_nonnegative(acb_realref(mu + end1)))
            ok = 0;
    }

    /* the kernel with the principal logarithms, in doubles */
    K.m = m;
    K.lam = lambda;
    K.konst = zero;
    K.mu = mu;
    K.al = al;
    K.be = be;
    K.dlam = cel_get(lambda);
    K.dkonst = cel_mk(0.0, 0.0);
    K.stable = 0;
    K.end0 = end0;
    K.end1 = infinite ? -1 : end1;
    K.E0 = zero;
    K.E1 = one;
    P0.K = &K;
    P0.infinite = infinite;
    P0.E0 = cel_mk(0.0, 0.0);
    P0.E1 = cel_mk(1.0, 0.0);
    P0.d0 = cel_mk(1.0, 0.0);
    for (j = 0; j < m; j++)
    {
        mag_t mg;

        K.dmu[j] = cel_get(mu + j);
        K.dal[j] = cel_get(al + j);
        K.dbe[j] = cel_get(be + j);
        P0.sing[j] = !cel_is0(K.dmu[j]) && !cel_is0(K.dbe[j]);
        tpabs[j] = P0.sing[j] ? cel_scale(cel_div(K.dal[j], K.dbe[j]), -1.0) : cel_mk(0.0, 0.0);
        P0.tp[j] = tpabs[j];
        P0.tilt[j] = cel_mk(0.0, 0.0);
        if (P0.sing[j])
            tix[nterm++] = j;
        if (!(cel_abs(K.dmu[j]) < 1e290) || !(cel_abs(K.dal[j]) < 1e290) || !(cel_abs(K.dbe[j]) < 1e290))
            ok = 0;

        /* (the bits of the parameters; numbers too large for doubles
           are rejected above) */
        mag_init(mg);
        acb_get_mag(mg, mu + j);
        if (mag_cmp_2exp_si(mg, 0) > 0 && mag_cmp_2exp_si(mg, 1100) < 0)
            bits = FLINT_MAX(bits, MAG_EXP(mg));
        acb_get_mag(mg, be + j);
        if (mag_cmp_2exp_si(mg, 0) > 0 && mag_cmp_2exp_si(mg, 1100) < 0)
            bits = FLINT_MAX(bits, MAG_EXP(mg));
        acb_get_mag(mg, lambda);
        if (mag_cmp_2exp_si(mg, 0) > 0 && mag_cmp_2exp_si(mg, 1100) < 0)
            bits = FLINT_MAX(bits, MAG_EXP(mg));
        mag_clear(mg);
    }
    if (!(cel_abs(K.dlam) < 1e290))
        ok = 0;

    /* (balls with huge radii: the vertices of the path are exact
       numbers computed from these) */
    {
        mag_t mg;

        mag_init(mg);
        acb_get_mag(mg, lambda);
        if (mag_cmp_2exp_si(mg, 1000) >= 0)
            ok = 0;
        for (j = 0; j < m; j++)
        {
            acb_get_mag(mg, mu + j);
            if (mag_cmp_2exp_si(mg, 1000) >= 0)
                ok = 0;
            acb_get_mag(mg, al + j);
            if (mag_cmp_2exp_si(mg, 1000) >= 0)
                ok = 0;
            acb_get_mag(mg, be + j);
            if (mag_cmp_2exp_si(mg, 1000) >= 0)
                ok = 0;
        }
        mag_clear(mg);
    }

    /* The tilted plan, for a polygon by way of via.  With real data and
       a singular point on the segment [0, 1] the line of descent from
       one real saddle point runs along the real axis into another one
       (the lines of the second lead around the singular point to the
       end), and is not a path.  For the planning, which is heuristic,
       the singular points between the end points are moved away from
       the side of via, by so little that the exponent changes by a
       phase of the order of one: the line then passes the second
       saddle point on the side of via and turns to the end point.
       The path is integrated with the kernel as it is. */
    if (tilted && via != NULL && ok)
    {
        double side = cel_get(via).im;

        for (j = 0; j < m; j++)
        {
            if (P0.sing[j] && j != end0 && j != end1 && side != 0.0 && K.dmu[j].re < 0.0)
            {
                double d = fmin(cel_abs(tpabs[j]), cel_abs(cel_sub(tpabs[j], cel_mk(1.0, 0.0))));
                double eta = fmin(0.05, 2.0 / (cel_abs(K.dmu[j]) + 1e-300));

                P0.tilt[j] = cel_mk(0.0, (side < 0.0) ? eta * d : -eta * d);
                P0.tp[j] = cel_add(tpabs[j], P0.tilt[j]);
            }
        }
    }

    wp0 = 192 + 2 * bits;
    sing0 = P0.sing[end0];
    sing1 = (!infinite && P0.sing[end1]);

    /* infinity: a valley in the directions where exp(lambda t) decays,
       or all around if lambda = 0 and the sum of the mu_j has a
       negative real part */
    if (!cel_is0(K.dlam))
        P0.inf_valley = 1;
    else
    {
        double s = 0.0;
        for (j = 0; j < m; j++)
            s += K.dmu[j].re;
        P0.inf_valley = (s < 0.0);
    }

    if (infinite)
    {
        /* the direction of the ray: a Gaussian integer, where
           exp(lambda t) decays */
        P0.d0 = cel_get(dir);
        if (cel_is0(K.dlam) || !acb_is_exact(dir) || !acb_is_finite(dir) || cel_is0(P0.d0)
            || fabs(P0.d0.re) > 2.0 || fabs(P0.d0.im) > 2.0
            || !(-(K.dlam.re / cel_abs(K.dlam)) * P0.d0.re + (K.dlam.im / cel_abs(K.dlam)) * P0.d0.im
                > 0.2 * cel_abs(P0.d0)))
            ok = 0;
    }

    P0.D = (prec + 20.0) * CEL_LN2 + 6.0;

    if (!ok)
        goto cleanup;

    tr = flint_malloc(sizeof(cel_trace) * 6);
    pcs = flint_malloc(sizeof(cel_piece) * 3);
    bpc = flint_malloc(sizeof(cel_piece) * 3);
    ab = flint_malloc(sizeof(cel_cd) * (3 * CEL_PIECEMAX + 8));

    ok = cel_saddles(sad, &ns, &K, tix, nterm, wp0);

    /* the frames and the lines of the saddle points */
    for (k = 0; k < ns && ok; k++)
    {
        cel_frame * Fr = F + nf;
        cel_cd dirn, st, o = cel_mk(0.0, 0.0);
        double a2, ang, w;

        ok = cel_frame_init(Fr, &K, &P0, sad + k, wp0);
        nf++;
        if (!ok)
            break;

        cel_phi(&dummy, &d1, &d2, &Fr->K, o, -1);
        a2 = cel_abs(d2);
        if (!(a2 > 0.0) || !(a2 < 1e300))
        {
            ok = 0;
            break;
        }
        /* u^2 phi'' < 0: u = i / sqrt(phi'') */
        ang = 0.5 * (CEL_PI - atan2(d2.im, d2.re));
        dirn = cel_mk(cos(ang), sin(ang));
        w = fmin(0.5 / sqrt(a2), 0.3 * cel_sing_dist(&Fr->P, o));
        st = cel_scale(dirn, w);
        /* (negligible: on a path of the length of the geometry, compared
           to the width of the saddle point) */
        Fr->P.D += log(fmax(cel_frame_scale(Fr, infinite) * sqrt(a2), 1.0));

        cel_trace_run(tr + ntr, &Fr->P, st);
        cel_trace_run(tr + ntr + 1, &Fr->P, cel_scale(st, -1.0));
        elab[ne][0] = tr[ntr].label;
        elab[ne][1] = tr[ntr + 1].label;
        etr[ne][0] = ntr;
        etr[ne][1] = ntr + 1;
        esad[ne] = 1;
        efr[ne] = nf - 1;
        ne++;
        ntr += 2;
    }

    /* the lines from the end points where the integrand does not vanish */
    for (i = 0; i < 2 && ok; i++)
    {
        cel_frame * Fr = F + nf;
        int je = (i == 0) ? end0 : end1;

        if ((i == 1 && infinite) || !cel_is0(K.dmu[je]))
            continue;
        ok = cel_frame_init(Fr, &K, &P0, (i == 0) ? zero : one, wp0);
        nf++;
        if (!ok)
            break;
        cel_phi(&dummy, &d1, &d2, &Fr->K, cel_mk(0.0, 0.0), -1);
        Fr->P.D += log(fmax(cel_frame_scale(Fr, infinite) * cel_abs(d1), 1.0));
        cel_trace_run(tr + ntr, &Fr->P, cel_mk(0.0, 0.0));
        elab[ne][0] = (i == 0) ? CEL_LE0 : CEL_LE1;
        elab[ne][1] = tr[ntr].label;
        etr[ne][0] = -1;
        etr[ne][1] = ntr;
        esad[ne] = 0;
        efr[ne] = nf - 1;
        ne++;
        ntr++;
    }

    if (!ok)
        goto cleanup;

    /* The routes: the chains of at most three edges from one end to the
       other through valleys between (infinity or zeros of the
       integrand), the shorter ones first, and last the contour itself */
    target = infinite ? CEL_LINF : CEL_LE1;
    {
        int a, b, cc, da, db, dc, len, v, v2, mid = ~(CEL_LE0 | CEL_LE1 | target);

        for (len = 1; len <= 3; len++)
        for (a = 0; a < ne; a++) for (da = 0; da < 2; da++)
        {
            if (!(elab[a][da] & CEL_LE0))
                continue;
            if ((elab[a][1 - da] & target) && len == 1 && ncand < 8)
            {
                cand[ncand][0] = a; cdir[ncand][0] = da; clab[ncand][0] = CEL_LE0; clen[ncand++] = 1;
            }
            for (v = 4; v < (8 << CEL_MAXTERMS); v <<= 1)
            {
                if (!(elab[a][1 - da] & v & mid))
                    continue;
                for (b = 0; b < ne; b++) for (db = 0; db < 2; db++)
                {
                    if (b == a || !(elab[b][db] & v))
                        continue;
                    if ((elab[b][1 - db] & target) && len == 2 && ncand < 8)
                    {
                        cand[ncand][0] = a; cdir[ncand][0] = da; clab[ncand][0] = CEL_LE0;
                        cand[ncand][1] = b; cdir[ncand][1] = db; clab[ncand][1] = v;
                        clen[ncand++] = 2;
                    }
                    for (v2 = 4; v2 < (8 << CEL_MAXTERMS); v2 <<= 1)
                    {
                        if (!(elab[b][1 - db] & v2 & mid) || v2 == v)
                            continue;
                        for (cc = 0; cc < ne; cc++) for (dc = 0; dc < 2; dc++)
                        {
                            if (cc == a || cc == b || !(elab[cc][dc] & v2) || !(elab[cc][1 - dc] & target))
                                continue;
                            if (len == 3 && ncand < 8)
                            {
                                cand[ncand][0] = a; cdir[ncand][0] = da; clab[ncand][0] = CEL_LE0;
                                cand[ncand][1] = b; cdir[ncand][1] = db; clab[ncand][1] = v;
                                cand[ncand][2] = cc; cdir[ncand][2] = dc; clab[ncand][2] = v2;
                                clen[ncand++] = 3;
                            }
                        }
                    }
                }
            }
        }
        clen[ncand++] = 0;
    }

    /* Of the routes for which cuts are found (the others pass a
       singular point on the wrong side), the lowest one is taken; the
       contour itself only if there is no other, and the parameters are
       not large. */
    for (ci = 0; ci < ncand; ci++)
    {
    int pass, narc = 0, found = 0, nroute = clen[ci];

    direct = (nroute == 0);
    if (direct && (have || bits > 20 || via != NULL))
        break;

    for (pass = 0; pass < 2 && !found; pass++)
    {
        int okc = 1, na = 0;

        npc = 0;
        narc = 0;

        if (direct)
        {
            /* The contour itself, which will do if the integrand does
               not oscillate much on it (small parameters, where the
               lines of descent can wind around the end points): in a
               frame at a point of the contour */
            cel_piece * pc = pcs;
            cel_frame * Fr;

            pc->n = 0;
            pc->saddle = 0;
            pc->junction = 0;
            pc->e0 = 1;
            pc->e1 = !infinite;

            if (fdirect < 0)
            {
                if (infinite)
                {
                    double sc = 1.0 / (cel_abs(K.dlam) + 1e-300);
                    for (j = 0; j < m; j++)
                        if (P0.sing[j] && j != end0)
                            sc = fmin(sc, cel_abs(tpabs[j]));
                    acb_set_d(t, ldexp(1.0, (int) floor(log2(0.125 * sc / cel_abs(P0.d0)))));
                    acb_mul(t, t, dir, wp0);
                }
                else
                {
                    acb_one(t);
                    acb_mul_2exp_si(t, t, -1);
                }
                fdirect = nf;
                okc = cel_frame_init(F + nf, &K, &P0, t, wp0);
                nf++;
            }
            pc->f = fdirect;
            Fr = F + fdirect;

            if (!okc)
            {
            }
            else if (!infinite)
            {
                for (i = 0; i <= 16 && okc; i++)
                    okc = cel_piece_push(pc, (i == 0) ? Fr->P.E0 : (i == 16) ? Fr->P.E1
                        : cel_mk(i / 16.0 - 0.5, 0.0));
                pc->sidx = 8;
            }
            else
            {
                cel_cd L[CEL_MAXTERMS], xq, xp = cel_mk(0.0, 0.0), p;
                double r = 1.0, hi = 0.0;

                for (j = 0; j < m; j++)
                    L[j] = cel_mk(0.0, 0.0);
                okc = cel_piece_push(pc, Fr->P.E0) && cel_piece_push(pc, xp);
                pc->sidx = 1;
                for (i = 0; i < 300 && okc; i++)
                {
                    /* (the points base (1.5^i), i = 1, 2, ...) */
                    r *= 1.5;
                    xq = cel_scale(Fr->babs, r - 1.0);
                    cel_logs_step(L, &Fr->K, xp, xq);
                    p = cel_phi_logs(&Fr->K, xq, L);
                    okc = cel_piece_push(pc, xq);
                    xp = xq;
                    hi = fmax(hi, p.re);
                    if (p.re <= hi - P0.D && cel_abs(K.dlam) * cel_abs(Fr->babs) * r
                        > 4.0 * cel_abs(K.dmu[0]) + 4.0 * cel_abs(K.dmu[m - 1]))
                        break;
                }
                if (i == 300)
                    okc = 0;
                toplen = cel_abs(Fr->babs) * r;
            }
            npc = 1;
        }

        /* The pieces: the lines of each edge of the route, from the
           valley where it is entered to the base point and on to the
           next valley (an edge from an end point is its one line), in
           the frame of the edge.  The first and the last vertex of the
           path are the end points, reached by chords from the lines. */
        for (i = 0; i < nroute && okc; i++)
        {
            int e = cand[ci][i], da = cdir[ci][i], lab = clab[ci][i];
            int ta = etr[e][da], tb = etr[e][1 - da];
            int rev = (!esad[e] && da == 1);   /* the line of an end point, backwards */
            int last = (i == nroute - 1), kk, nxs;
            cel_cd xs[24];
            cel_piece * pc = pcs + npc;
            cel_frame * Fr = F + efr[e];

            if (rev)
            {
                ta = etr[e][1];
                tb = -1;
            }

            pc->f = efr[e];
            pc->n = 0;
            pc->sidx = -1;
            pc->saddle = esad[e];
            pc->junction = (i > 0);
            pc->e0 = (i == 0);
            pc->e1 = (last && !infinite);

            if (i == 0 && ta >= 0)
            {
                okc = cel_piece_push(pc, Fr->P.E0);
                /* (and the points of the chord from where the integrand
                   is negligible) */
                nxs = cel_chord_split(xs, 24, Fr, tr[ta].t[tr[ta].eidx[0]], Fr->P.E0);
                for (k = nxs - 1; k >= 0 && okc; k--)
                {
                    okc = cel_piece_push(pc, xs[k]);
                    pc->fix[pc->n - 1] = 1;
                }
            }

            if (i > 0)
            {
                /* from the last vertex of the piece before (as an exact
                   number), through the valley between: far out, or
                   around a zero */
                cel_piece * pv = pcs + npc - 1;
                cel_cd p1, p2 = tr[ta].t[tr[ta].n - 1];

                acb_set_d_d(Jv + npc, pv->t[pv->n - 1].re, pv->t[pv->n - 1].im);
                acb_add(Jv + npc, Jv + npc, F[pv->f].base, 4 * wp0 + 8192);
                acb_sub(Jv + npc, Jv + npc, Fr->base, 4 * wp0 + 8192);
                p1 = cel_get(Jv + npc);
                okc = cel_piece_push(pc, p1);

                if (lab == CEL_LINF)
                {
                    cel_cd a1 = cel_add(Fr->babs, p1), a2 = cel_add(Fr->babs, p2);
                    double r1 = cel_abs(a1), r2 = cel_abs(a2);

                    cel_cd q1, q2, qs;
                    int rises = 0, s3;

                    /* (unless the chord stays as low as its ends) */
                    cel_phi(&q1, &d1, &d2, &Fr->K, p1, -1);
                    cel_phi(&q2, &d1, &d2, &Fr->K, p2, -1);
                    for (s3 = 1; s3 < 8; s3++)
                    {
                        cel_phi(&qs, &d1, &d2, &Fr->K, cel_add(p1, cel_scale(cel_sub(p2, p1), 0.125 * s3)), -1);
                        if (!(qs.re <= fmax(q1.re, q2.re) + 2.0))
                            rises = 1;
                    }

                    if (rises && (a1.re * a2.re + a1.im * a2.im) < 0.5 * r1 * r2)
                    {
                        /* by way of the direction of the fastest decay */
                        cel_cd dd = cel_is0(K.dlam) ? cel_add(a1, a2)
                            : cel_scale(cel_conj(K.dlam), -1.0 / cel_abs(K.dlam));
                        if (cel_abs(dd) > 0.0)
                            okc = okc && cel_piece_push(pc, cel_sub(cel_scale(dd,
                                2.0 * fmax(r1, r2) / cel_abs(dd)), Fr->babs));
                    }
                }
                else
                {
                    /* the short way around the zero, or the long way if
                       no cuts are found for that */
                    cel_cd z0 = cel_mk(0.0, 0.0), v1, v2;
                    double a1, a2, da2, r1, r2;
                    int steps, s2;

                    for (j = 0; j < m; j++)
                        if (lab == CEL_LV(j))
                            z0 = Fr->P.tp[j];
                    v1 = cel_sub(p1, z0);
                    v2 = cel_sub(p2, z0);
                    a1 = atan2(v1.im, v1.re);
                    a2 = atan2(v2.im, v2.re);
                    r1 = cel_abs(v1);
                    r2 = cel_abs(v2);
                    narc++;
                    da2 = a2 - a1;
                    while (da2 > CEL_PI) da2 -= 2.0 * CEL_PI;
                    while (da2 < -CEL_PI) da2 += 2.0 * CEL_PI;
                    if (pass == 1)
                        da2 += (da2 > 0.0) ? -2.0 * CEL_PI : 2.0 * CEL_PI;
                    steps = (int) (fabs(da2) / 0.5) + 1;
                    for (s2 = 1; s2 < steps && okc; s2++)
                    {
                        double f = (double) s2 / steps, r = r1 * pow(r2 / r1, f), an = a1 + f * da2;
                        okc = cel_piece_push(pc, cel_add(z0, cel_mk(r * cos(an), r * sin(an))));
                    }
                }
            }

            if (ta >= 0)
            {
                /* towards the base point (from the point of the chord
                   to the end point, on the first piece) */
                kk = (i == 0) ? tr[ta].eidx[0] : tr[ta].n - 1;
                for (k = kk; k >= 1 && okc; k--)
                    okc = cel_piece_push(pc, tr[ta].t[k]);
            }
            if (okc)
            {
                okc = cel_piece_push(pc, cel_mk(0.0, 0.0));
                pc->sidx = pc->n - 1;
            }
            if (tb >= 0)
            {
                kk = (last && !infinite) ? tr[tb].eidx[1] : tr[tb].n - 1;
                for (k = 1; k <= kk && okc; k++)
                    okc = cel_piece_push(pc, tr[tb].t[k]);
            }
            if (last && !infinite && tb >= 0 && okc)
            {
                nxs = cel_chord_split(xs, 24, Fr, tr[tb].t[tr[tb].eidx[1]], Fr->P.E1);
                for (k = 0; k < nxs && okc; k++)
                {
                    okc = cel_piece_push(pc, xs[k]);
                    pc->fix[pc->n - 1] = 1;
                }
                okc = okc && cel_piece_push(pc, Fr->P.E1);
            }
            npc++;
        }

        if (!okc || npc == 0)
            break;

        /* the values along the pieces, and the path in absolute
           coordinates (coarsely, for the choice of the cuts) */
        for (i = 0; i < npc; i++)
        {
            cel_piece * pc = pcs + i;

            if (pc->n < 2 || pc->sidx < 0)
                okc = 0;
            else
                cel_piece_track(pc, F + pc->f, sing0, sing1);
            for (k = 0; k < pc->n; k++)
                ab[na++] = cel_add(F[pc->f].babs, pc->t[k]);
        }
        if (!okc)
            break;
        ab[0] = P0.E0;
        if (!infinite)
            ab[na - 1] = P0.E1;

        /* The cuts: for each term, the ray which avoids the path and
           the contour with the widest margin */
        {
            cel_cd far = infinite ? cel_scale(P0.d0, 1e3 * (1.0 + cel_abs(ab[na - 1]))) : P0.E1;
            cel_cd dvia = (via != NULL) ? cel_get(via) : cel_mk(0.0, 0.0);
            int allok = 1;

            for (j = 0; j < m && allok; j++)
            {
                double bestcl = 0.12;
                int bi = -1, a2i;

                acb_set(ral + j, al + j);
                acb_set(rbe + j, be + j);
                if (!P0.sing[j])
                    continue;

                for (a2i = 0; a2i < 8; a2i++)
                {
                    cel_cd cc = cel_mk(cel_rot[a2i][0], cel_rot[a2i][1]);
                    /* the ray: t_j - x / (c be_j), x > 0 */
                    cel_cd dd = cel_scale(cel_div(cel_mk(1.0, 0.0), cel_mul(cc, K.dbe[j])), -1.0);
                    double cl = INFINITY, an;
                    int meets = 0;

                    dd = cel_scale(dd, 1.0 / cel_abs(dd));
                    for (k = 0; k + 1 < na && !meets; k++)
                        meets = cel_ray_meets(tpabs[j], dd, ab[k], ab[k + 1]);
                    if (via == NULL)
                        meets = meets || cel_ray_meets(tpabs[j], dd, P0.E0, far);
                    else
                        meets = meets || cel_ray_meets(tpabs[j], dd, P0.E0, dvia)
                            || cel_ray_meets(tpabs[j], dd, dvia, far);
                    if (infinite)
                        meets = meets || cel_ray_meets(tpabs[j], dd, ab[na - 1], cel_add(ab[na - 1], far));
                    if (meets)
                        continue;

                    /* the margin: the smallest angle to a vertex */
                    for (k = 0; k <= na + 1; k++)
                    {
                        cel_cd v = (k < na) ? ab[k] : (k == na) ? far
                            : (via != NULL) ? dvia : cel_scale(cel_add(P0.E0, far), 0.5);
                        v = cel_sub(v, tpabs[j]);
                        if (cel_abs(v) <= 1e-14 * (1.0 + cel_abs(tpabs[j])))
                            continue;
                        an = fabs(atan2(v.re * dd.im - v.im * dd.re, v.re * dd.re + v.im * dd.im));
                        cl = fmin(cl, an);
                    }
                    if (cl > bestcl + 0.05)
                    {
                        bestcl = cl;
                        bi = a2i;
                    }
                }

                if (bi < 0)
                    allok = 0;
                else
                {
                    acb_set_si_si(t, cel_rot[bi][0], cel_rot[bi][1]);
                    acb_mul(ral + j, al + j, t, 16384);
                    acb_mul(rbe + j, be + j, t, 16384);
                }
            }
            found = allok;
        }

        if (found)
        {
            /* The height of the path with the logarithms of the new
               cuts: at the base points in ball arithmetic, with
               k = sum_j mu_j gamma_j, gamma_j = Log(w) - Log(c w) at a
               point of the contour */
            cel_kernel KR = K;

            acb_zero(konst);
            if (infinite)
                acb_set(u, dir);
            else if (via != NULL)
                acb_set(u, via);
            else
            {
                acb_one(u);
                acb_mul_2exp_si(u, u, -1);
            }
            for (j = 0; j < m; j++)
            {
                if (acb_is_zero(mu + j))
                    continue;
                acb_mul(t, be + j, u, wp0);
                acb_add(t, t, al + j, wp0);
                acb_log(c, t, wp0);
                acb_mul(t, rbe + j, u, wp0);
                acb_add(t, t, ral + j, wp0);
                acb_log(t, t, wp0);
                acb_sub(c, c, t, wp0);
                acb_addmul(konst, c, mu + j, wp0);
            }

            KR.al = ral;
            KR.be = rbe;
            KR.konst = konst;
            for (i = 0; i < npc; i++)
            {
                double mx = -INFINITY;

                cel_phi_acb(H + i, NULL, NULL, &KR, F[pcs[i].f].base, -1, wp0);
                for (k = 0; k < pcs[i].n; k++)
                    mx = fmax(mx, pcs[i].phi[k].re);
                arb_set_d(x, mx);
                arb_add(x, x, acb_realref(H + i), wp0);
                if (i == 0)
                    arb_set(score, x);
                else
                    arb_max(score, score, x, wp0);
            }

            /* (lower by a margin: the first of equals) */
            arb_sub(x, bscore, score, wp0);
            if (arb_is_finite(score) && (!have
                || arf_cmp_d(arb_midref(x), 0.01) > 0))
            {
                have = 1;
                arb_set(bscore, score);
                bdirect = direct;
                nbest = npc;
                memcpy(bpc, pcs, sizeof(cel_piece) * npc);
                _acb_vec_set(bJ, Jv, 3);
                _acb_vec_set(bH, H, 3);
                _acb_vec_set(bal, ral, m);
                _acb_vec_set(bbe, rbe, m);
                acb_set(bkonst, konst);
            }
        }

        if (narc == 0)
            break;
    }
    }

    if (!have)
    {
        ok = 0;
        goto cleanup;
    }
    direct = bdirect;

    /* Rigorous: the contour avoids the principal and the new cuts, and
       k for these */
    {
        slong wp = wp0 + prec;
        acb_t d, half, w, lw;

        acb_init(d); acb_init(half); acb_init(w); acb_init(lw);

        if (infinite)
        {
            acb_set(d, dir);
            acb_set(half, dir);
        }
        else if (via != NULL)
        {
            /* (the two chords of the polygon: d from 0, and 1 - via,
               an exact number, from via) */
            acb_set(d, via);
            acb_set(half, via);
            acb_sub(w, one, via, ARF_PREC_EXACT);
        }
        else
        {
            acb_one(d);
            acb_one(half);
            acb_mul_2exp_si(half, half, -1);
        }

        acb_zero(konst);
        for (j = 0; j < m && ok; j++)
        {
            if (acb_is_zero(mu + j))
                continue;
            if (via != NULL)
            {
                acb_sub(w, one, via, ARF_PREC_EXACT);
                ok = cel_avoids_cut(al + j, be + j, zero, d, 0, j == end0, 0, wp)
                    && cel_avoids_cut(bal + j, bbe + j, zero, d, 0, j == end0, 0, wp)
                    && cel_avoids_cut(al + j, be + j, via, w, 0, 0, j == end1, wp)
                    && cel_avoids_cut(bal + j, bbe + j, via, w, 0, 0, j == end1, wp);
            }
            else
            ok = cel_avoids_cut(al + j, be + j, zero, d, infinite, j == end0, !infinite && j == end1, wp)
                && cel_avoids_cut(bal + j, bbe + j, zero, d, infinite, j == end0, !infinite && j == end1, wp);
            acb_mul(w, be + j, half, wp);
            acb_add(w, w, al + j, wp);
            acb_log(lw, w, wp);
            acb_mul(w, bbe + j, half, wp);
            acb_add(w, w, bal + j, wp);
            acb_log(w, w, wp);
            acb_sub(lw, lw, w, wp);
            acb_addmul(konst, lw, mu + j, wp);
        }

        acb_clear(d); acb_clear(half); acb_clear(w); acb_clear(lw);
    }

    if (!ok)
        goto cleanup;

    /* The pieces are integrated, each in its frame: with t = base + delta,
       the kernel has al_j + be_j base for al_j and k + lambda base for k. */
    {
        acb_ptr * pts = flint_malloc(sizeof(acb_ptr) * nbest);
        slong * npts = flint_malloc(sizeof(slong) * nbest);
        int * keep = flint_malloc(sizeof(int) * CEL_PIECEMAX);
        acb_ptr cs = _acb_vec_init(nbest), tk = _acb_vec_init(nbest);
        arb_ptr vI = _arb_vec_init(nbest);
        cel_kernel KT = K;
        slong wp = wp0 + prec + 64;
        double logIrel, s = 0.0;

        KT.al = tal;
        KT.be = bbe;

        for (i = 0; i < nbest && ok; i++)
        {
            cel_piece * pc = bpc + i;
            cel_frame * Fr = F + pc->f;
            double neg, mx = -INFINITY, corr;
            int nk, kt = pc->sidx;

            /* the chords */
            arb_sub(x, bscore, acb_realref(bH + i), wp0);
            neg = arf_get_d(arb_midref(x), ARF_RND_NEAR) - Fr->P.D;
            nk = cel_piece_chords(keep, pc, Fr, neg);
            pts[i] = _acb_vec_init(nk);
            npts[i] = nk;
            for (k = 0; k < nk; k++)
                cel_piece_vertex(pts[i] + k, pc, keep[k], Fr, bJ + i, infinite);

            /* the reference value: the exponent at the highest vertex */
            for (k = 0; k < pc->n; k++)
                if (pc->phi[k].re > mx)
                {
                    mx = pc->phi[k].re;
                    kt = k;
                }
            /* (exactly, for exact parameters: a term vanishes exactly
               at its end point) */
            for (j = 0; j < m; j++)
            {
                acb_mul(tal + j, bbe + j, Fr->base, 4 * wp0 + 8192);
                acb_add(tal + j, tal + j, bal + j, 4 * wp0 + 8192);
            }
            acb_mul(tk + i, lambda, Fr->base, wp);
            acb_add(tk + i, tk + i, konst, wp);
            KT.konst = tk + i;
            cel_piece_vertex(t, pc, kt, Fr, bJ + i, infinite);
            cel_phi_acb(cs + i, NULL, NULL, &KT, t, -1, wp0);
            acb_get_mid(cs + i, cs + i);
            if (!acb_is_finite(cs + i))
                ok = 0;

            /* the size of the integral: from a saddle point, or from
               the descent at an end point */
            cel_phi(&dummy, &d1, &d2, &Fr->K, pc->t[kt], -1);
            if (direct)
                corr = log(infinite ? toplen : 1.0);
            else if (pc->saddle && kt == pc->sidx)
                corr = 0.5 * log(2.0 * CEL_PI / (cel_abs(d2) + 1e-300));
            else
                corr = -log(cel_abs(d1) + 1e-300);
            arb_set_d(x, corr);
            arb_add(vI + i, acb_realref(cs + i), x, wp0);
            if (i == 0 || arf_cmp(arb_midref(vI + i), arb_midref(vmax)) > 0)
                arb_set(vmax, vI + i);
        }

        for (i = 0; i < nbest && ok; i++)
        {
            arb_sub(x, vI + i, vmax, wp0);
            s += exp(arf_get_d(arb_midref(x), ARF_RND_NEAR));
        }
        logIrel = log(s);
        if (ok && logscale > -INFINITY)
        {
            arb_set_d(x, logscale);
            arb_sub(x, x, vmax, wp0);
            if (arf_get_d(arb_midref(x), ARF_RND_NEAR) > logIrel)
                logIrel = arf_get_d(arb_midref(x), ARF_RND_NEAR);
        }
        logI = arf_get_d(arb_midref(vmax), ARF_RND_NEAR) + logIrel;

        acb_zero(total);
        for (i = 0; i < nbest && ok; i++)
        {
            cel_piece * pc = bpc + i;
            cel_frame * Fr = F + pc->f;
            double logtol;

            /* (exactly, for exact parameters: a term vanishes exactly
               at its end point) */
            for (j = 0; j < m; j++)
            {
                acb_mul(tal + j, bbe + j, Fr->base, 4 * wp0 + 8192);
                acb_add(tal + j, tal + j, bal + j, 4 * wp0 + 8192);
            }
            arb_sub(x, vmax, acb_realref(cs + i), wp0);
            logtol = arf_get_d(arb_midref(x), ARF_RND_NEAR) + logIrel - (prec + 5.0) * CEL_LN2 - log(32.0);

            ok = _acb_hypgeom_exp_log_path(t, lambda, mu, tal, bbe, tk + i, Fr->g1, m,
                (pc->e0 && i == 0) ? end0 : -1, (pc->e1 && i == nbest - 1 && !infinite) ? end1 : -1,
                pts[i], npts[i], cs + i, logtol, prec);
            if (ok)
                acb_add(total, total, t, prec + 30);
        }

        /* The tail on the ray from the last vertex T in the direction d:
           the cuts are not met, and with |w_j| increasing along the ray,
           Re(d phi') <= Re(d lambda) + sum_j |d mu_j be_j| / |w_j(T)| = -kappa,
           so that the tail is at most |d| |exp(phi(T))| / kappa. */
        if (ok && infinite)
        {
            acb_t w, v;
            arb_t kap;
            mag_t e;
            acb_srcptr T = pts[nbest - 1] + npts[nbest - 1] - 1;

            acb_init(w); acb_init(v);
            arb_init(kap);
            mag_init(e);

            KT.konst = tk + nbest - 1;
            acb_mul(v, lambda, dir, wp);
            arb_neg(kap, acb_realref(v));
            for (j = 0; j < m && ok; j++)
            {
                if (acb_is_zero(mu + j))
                    continue;
                ok = cel_avoids_cut(tal + j, bbe + j, T, dir, 1, 0, 0, wp);
                acb_mul(w, bbe + j, T, wp);
                acb_add(w, w, tal + j, wp);
                /* |w| does not decrease: Re(conj(w) be d) >= 0 */
                acb_mul(v, bbe + j, dir, wp);
                acb_conj(t, w);
                acb_mul(t, t, v, wp);
                if (!arb_is_nonnegative(acb_realref(t)))
                    ok = 0;
                acb_mul(v, v, mu + j, wp);
                acb_div(v, v, w, wp);
                acb_abs(x, v, wp);
                arb_sub(kap, kap, x, wp);
            }

            if (!ok || !arb_is_positive(kap))
                ok = 0;
            else
            {
                cel_phi_acb(v, NULL, NULL, &KT, T, -1, wp);
                arb_exp(x, acb_realref(v), wp);
                arb_div(x, x, kap, wp);
                acb_abs(kap, dir, wp);
                arb_mul(x, x, kap, wp);
                arb_get_mag(e, x);
                acb_add_error_mag(total, e);
                ok = acb_is_finite(total);
            }

            acb_clear(w); acb_clear(v);
            arb_clear(kap);
            mag_clear(e);
        }

        for (i = 0; i < nbest; i++)
            if (i < nbest && pts[i] != NULL)
                _acb_vec_clear(pts[i], npts[i]);
        flint_free(pts);
        flint_free(npts);
        flint_free(keep);
        _acb_vec_clear(cs, nbest);
        _acb_vec_clear(tk, nbest);
        _arb_vec_clear(vI, nbest);
    }

    if (ok)
        acb_set(res, total);

    /* On the contour itself the integrand can cancel: if the integral
       is smaller than was assumed for the tolerance, once more with a
       higher precision. */
    if (ok && direct)
    {
        mag_t lo;
        double lost;

        mag_init(lo);
        acb_get_mag_lower(lo, res);
        if (mag_is_zero(lo))
            ok = 0;
        else
        {
            /* (log2 |res| from the exponent, which mag_get_d_log2_approx
               clamps to the range of a small fmpz, 2^30 on 32-bit
               systems) */
            lost = logI - (fmpz_get_d(MAG_EXPREF(lo))
                + log2(MAG_MAN(lo) * ldexp(1.0, -MAG_BITS))) * CEL_LN2;
            if (lost > 4.0 && level == 0)
            {
                if (lost < 0.5 * prec + 100.0)
                    redo = prec + (slong) (lost / CEL_LN2) + 10;
                else
                    ok = 0;
            }
        }
        mag_clear(lo);
    }

cleanup:
    if (tr != NULL) flint_free(tr);
    if (pcs != NULL) flint_free(pcs);
    if (bpc != NULL) flint_free(bpc);
    if (ab != NULL) flint_free(ab);
    for (i = 0; i < nf; i++)
        cel_frame_clear(F + i, m);
    acb_clear(zero); acb_clear(one); acb_clear(konst); acb_clear(bkonst); acb_clear(c); acb_clear(t);
    acb_clear(u); acb_clear(total);
    arb_clear(score); arb_clear(bscore); arb_clear(vmax); arb_clear(x);
    _acb_vec_clear(ral, m);
    _acb_vec_clear(rbe, m);
    _acb_vec_clear(bal, m);
    _acb_vec_clear(bbe, m);
    _acb_vec_clear(tal, m);
    _acb_vec_clear(sad, 2);
    _acb_vec_clear(Jv, 3);
    _acb_vec_clear(bJ, 3);
    _acb_vec_clear(H, 3);
    _acb_vec_clear(bH, 3);

    if (ok && redo != 0)
        return cel_integrate(res, dir, lambda, mu, al, be, m, end0, end1, infinite, logscale, redo, 1, via, tilted);

    return ok;
}

/* res = the integral of exp(lambda t + sum_j mu_j Log(al_j + be_j t))
   over [0, 1], or over the ray from 0 in the direction dir (a Gaussian
   integer, in a direction of decay of exp(lambda t)) if infinite,
   where the term end0 vanishes at 0 and
   the term end1 at 1 (-1 if infinite).  The logarithms must be
   continuous on the contour (this is checked).  The tolerance is
   2^-prec relative to the integral, or to exp(logscale) if that is
   larger.  Returns 0 if the integral is not computed. */
int
_acb_hypgeom_integrate_exp_log(acb_t res, const acb_t dir, const acb_t lambda, acb_srcptr mu,
    acb_srcptr al, acb_srcptr be, int m, int end0, int end1, int infinite, double logscale,
    slong prec)
{
    return cel_integrate(res, dir, lambda, mu, al, be, m, end0, end1, infinite, logscale, prec, 0, NULL, 0);
}

/* The same over the polygon from 0 by way of the exact number via to
   1: the integral over any contour from 0 to 1 which can be deformed
   into this polygon without crossing a point where a term vanishes
   (one of them can lie on the segment [0, 1] itself, or beyond it as
   seen from via, so that the integral is the continuation of the one
   over [0, 1] to the other side of that point).  The logarithms must
   be continuous on the polygon.  If no path is found, the path is
   planned once more with the singular points moved away from the side
   of via (see cel_integrate), or in the other order if tilted is
   nonzero (for a singular point on the segment or close to it with
   real data, where the first plan would fail or follow rounding
   errors). */
int
_acb_hypgeom_integrate_exp_log_via(acb_t res, const acb_t via, const acb_t lambda, acb_srcptr mu,
    acb_srcptr al, acb_srcptr be, int m, int end0, int end1, double logscale, int tilted,
    slong prec)
{
    if (!acb_is_exact(via) || !acb_is_finite(via) || acb_is_zero(via) || acb_is_one(via))
        return 0;

    tilted = (tilted != 0);

    return cel_integrate(res, via, lambda, mu, al, be, m, end0, end1, 0, logscale, prec, 0, via, tilted)
        || cel_integrate(res, via, lambda, mu, al, be, m, end0, end1, 0, logscale, prec, 0, via, !tilted);
}
