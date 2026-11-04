/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb_mat.h"
#include "acb_mat.h"
#include "arb_hypgeom.h"
#include "acb_theta.h"
#include "thread_support.h"
#include "impl.h"

/* The two Crandall sums as Mellin integrals of theta functions.

   With G_a(t) = int_1^oo w^(a-1) e^(-t w) dw,

     S = sum'_n e(-p.(n + pofs)) G_a(t_n) = int_1^oo w^(a-1) Theta(w) dw,
     Theta(w) = sum'_n e(-p.(n + pofs)) e^(-w t_n),  t_n = (n - sg)^T Q (n - sg),

   (sum' omitting the excluded points of _acb_mzeta_epstein_skip), and

     sum_n e(-p.n) e^(-w t_n) = e^(-w sg^T Q sg) theta(z, tau),
     tau = i w Q / pi,  z = -p - i w Q sg / pi

   (theta = acb_theta_one with zero characteristic).  The cost is one
   theta value per quadrature node, independent of the number of lattice
   points below the truncation, which grows like P^(d/2) for P bits.

   Quadrature: w = e^y, f(y) = e^(a y) Theta(e^y), Gauss-Legendre on
   y-intervals of length 2h.  f is analytic for |Im y| < pi/2 (Re w > 0),
   and on the Bernstein ellipse E_rho of an interval with semi-minor axis
   theta_max < pi/2, Re w >= r_min = e^(y_c - h (rho + 1/rho)/2) cos(theta_max),
   |Theta(w)| <= Phi(Re w),

     Phi(r) = sum'_n e^(-r t_n) <= min(Phi_all(r), e^(-(r - r1) t_min) Phi_all(r1)),
     Phi_all(r) = prod_j (1 + sqrt(pi / r) / C_jj)

   (C the upper Cholesky factor of Q, summing one coordinate at a time;
   t_min a lower bound for the t_n of the non-excluded points, r1 <= r),
   so the error of n nodes is at most (64/15) h M rho^(-2n) / (rho^2 - 1)
   with M the bound of |f| on the ellipse.  The tail beyond W = e^Y is
   bounded by Phi_all(W/2) e^(W t_min/2) int_W^oo w^(Re a - 1) e^(-w t_min) dw.
   The node counts and the working precisions taper with these bounds. */

#define THETA_MAX 1.4
#define HALF_LEN 1.5

/* Phi_all(r) = prod_j (1 + sqrt(pi / r) / C_jj) as a mag, r > 0 */
static void
phi_all(mag_t res, const arb_mat_t C, const arb_t r)
{
    slong d = arb_mat_nrows(C), j;
    arb_t t, u;
    arb_init(t);
    arb_init(u);
    arb_const_pi(t, 30);
    arb_div(t, t, r, 30);
    arb_sqrt(t, t, 30);
    arb_one(u);
    for (j = 0; j < d; j++)
    {
        arb_t v;
        arb_init(v);
        arb_div(v, t, arb_mat_entry(C, j, j), 30);
        arb_add_ui(v, v, 1, 30);
        arb_mul(u, u, v, 30);
        arb_clear(v);
    }
    arb_get_mag(res, u);
    arb_clear(t);
    arb_clear(u);
}

/* Phi(r) <= min(Phi_all(r), e^(-(r/2) tmin) Phi_all(r/2)) */
static void
phi_bound(mag_t res, const arb_mat_t C, const arb_t tmin, const arb_t r)
{
    arb_t r2, e;
    mag_t m;
    arb_init(r2);
    arb_init(e);
    mag_init(m);
    phi_all(res, C, r);
    arb_mul_2exp_si(r2, r, -1);
    phi_all(m, C, r2);
    arb_mul(e, r2, tmin, 30);
    arb_neg(e, e);
    arb_exp(e, e, 30);
    {
        mag_t em;
        mag_init(em);
        arb_get_mag(em, e);
        mag_mul(m, m, em);
        mag_clear(em);
    }
    if (mag_cmp(m, res) < 0)
        mag_swap(res, m);
    arb_clear(r2);
    arb_clear(e);
    mag_clear(m);
}

typedef struct
{
    const _acb_mzeta_epstein_problem_struct * P;
    int side;
    slong d;
    acb_srcptr a;
    arb_mat_struct * Q;     /* the form */
    arb_mat_struct * C;     /* its upper Cholesky factor */
    arb_ptr sg, p, pofs;
    arb_t sQs;              /* sg^T Q sg */
    arb_ptr Qsg;
    slong nexcl;
    slong * excl;           /* excluded points */
}
theta_side_struct;

/* Theta(w) for real w (ball), at precision prec */
static void
theta_eval(acb_t res, const theta_side_struct * S, const arb_t w, slong prec)
{
    slong d = S->d, i, j;
    acb_mat_t tau;
    acb_ptr z;
    arb_t c, t;
    acb_t ph, e;

    acb_mat_init(tau, d, d);
    z = _acb_vec_init(d);
    arb_init(c);
    arb_init(t);
    acb_init(ph);
    acb_init(e);

    arb_const_pi(c, prec);
    arb_div(c, w, c, prec);     /* w / pi */
    for (i = 0; i < d; i++)
        for (j = 0; j < d; j++)
            arb_mul(acb_imagref(acb_mat_entry(tau, i, j)), arb_mat_entry(S->Q, i, j), c, prec);
    for (i = 0; i < d; i++)
    {
        arb_neg(acb_realref(z + i), S->p + i);
        arb_mul(acb_imagref(z + i), S->Qsg + i, c, prec);
        arb_neg(acb_imagref(z + i), acb_imagref(z + i));
    }
    acb_theta_one(res, z, tau, 0, prec);
    arb_mul(t, w, S->sQs, prec);
    arb_neg(t, t);
    arb_exp(t, t, prec);
    acb_mul_arb(res, res, t, prec);

    /* the excluded points: e(-p.n) e^(-w t_n) */
    for (i = 0; i < S->nexcl; i++)
    {
        const slong * n = S->excl + i * d;
        arb_ptr v = _arb_vec_init(d);
        arb_t q;
        arb_init(q);
        for (j = 0; j < d; j++)
            arb_sub_si(v + j, S->sg + j, n[j], prec);
        /* t_n = v^T Q v, phase e(-p.n) */
        arb_zero(q);
        for (j = 0; j < d; j++)
        {
            slong k;
            arb_t s;
            arb_init(s);
            for (k = 0; k < d; k++)
                arb_addmul(s, arb_mat_entry(S->Q, j, k), v + k, prec);
            arb_addmul(q, s, v + j, prec);
            arb_clear(s);
        }
        arb_mul(q, q, w, prec);
        arb_neg(q, q);
        arb_exp(t, q, prec);
        arb_zero(q);
        for (j = 0; j < d; j++)
            arb_submul_si(q, S->p + j, n[j], prec);
        arb_mul_2exp_si(q, q, 1);
        arb_sin_cos_pi(acb_imagref(ph), acb_realref(ph), q, prec);
        acb_mul_arb(ph, ph, t, prec);
        acb_sub(res, res, ph, prec);
        arb_clear(q);
        _arb_vec_clear(v, d);
    }

    /* the constant phase e(-p.pofs) */
    arb_zero(c);
    for (j = 0; j < d; j++)
        arb_submul(c, S->p + j, S->pofs + j, prec);
    arb_mul_2exp_si(c, c, 1);
    arb_sin_cos_pi(acb_imagref(ph), acb_realref(ph), c, prec);
    acb_mul(res, res, ph, prec);

    acb_mat_clear(tau);
    _acb_vec_clear(z, d);
    arb_clear(c);
    arb_clear(t);
    acb_clear(ph);
    acb_clear(e);
}

/* lower bound for min t_n over the non-excluded points */
static void
tmin_bound(arb_t tmin, const theta_side_struct * S, slong prec)
{
    slong d = S->d, i, j, k, nb;
    arb_t R;
    arf_t R2;
    arb_ptr vc;
    acb_theta_eld_t E;
    slong * pts;

    arb_init(R);
    arf_init(R2);
    vc = _arb_vec_init(d);
    acb_theta_eld_init(E, d, d);

    /* radius^2: the smallest squared diagonal entry of C (few points) */
    arb_sqr(R, arb_mat_entry(S->C, 0, 0), prec);
    for (j = 1; j < d; j++)
    {
        arb_t t;
        arb_init(t);
        arb_sqr(t, arb_mat_entry(S->C, j, j), prec);
        arb_min(R, R, t, prec);
        arb_clear(t);
    }
    arb_get_lbound_arf(R2, R, prec);
    arb_set_arf(tmin, R2);

    /* |C n - C sg|^2 <= R2 */
    for (i = 0; i < d; i++)
    {
        arb_zero(vc + i);
        for (j = i; j < d; j++)
            arb_submul(vc + i, arb_mat_entry(S->C, i, j), S->sg + j, prec);
    }
    if (acb_theta_eld_set(E, S->C, R2, vc))
    {
        nb = acb_theta_eld_nb_pts(E);
        pts = flint_malloc(FLINT_MAX(nb, 1) * d * sizeof(slong));
        acb_theta_eld_points(pts, E);
        for (k = 0; k < nb; k++)
        {
            const slong * n = pts + k * d;
            slong e;
            int skip = 0;
            arb_t q, s;
            for (e = 0; e < S->nexcl && !skip; e++)
            {
                skip = 1;
                for (j = 0; j < d; j++)
                    skip = skip && (n[j] == S->excl[e * d + j]);
            }
            if (skip)
                continue;
            arb_init(q);
            arb_init(s);
            for (i = 0; i < d; i++)
            {
                arb_zero(s);
                for (j = i; j < d; j++)
                {
                    arb_t v;
                    arb_init(v);
                    arb_sub_si(v, S->sg + j, n[j], prec);
                    arb_addmul(s, arb_mat_entry(S->C, i, j), v, prec);
                    arb_clear(v);
                }
                arb_addmul(q, s, s, prec);
            }
            arb_get_lbound_arf(arb_midref(s), q, prec);
            mag_zero(arb_radref(s));
            if (arf_sgn(arb_midref(s)) < 0)
                arb_zero(s);
            arb_min(tmin, tmin, s, prec);
            arb_clear(q);
            arb_clear(s);
        }
        flint_free(pts);
    }
    else
        arb_zero(tmin);

    arb_clear(R);
    arf_clear(R2);
    _arb_vec_clear(vc, d);
    acb_theta_eld_clear(E);
}

/* log2 Phi_all(r) in double (parameters only) */
static double
d_log2_phi(const theta_side_struct * S, double r)
{
    slong j;
    double s = 0.0;
    for (j = 0; j < S->d; j++)
        s += log2(1.0 + sqrt(3.141592653589793 / r)
            / arf_get_d(arb_midref(arb_mat_entry(S->C, j, j)), ARF_RND_NEAR));
    return s;
}

/* Direct sums.  For large w, Theta(w) has few terms of any size: with
   the points n (not excluded) with t_n <= T0, their t_n and phases
   e(-p.n) precomputed,

     Theta(w) = e(-p.pofs) sum_(t_n <= T) e(-p.n) e^(-w t_n) + tail,
     |tail| <= e^(-(w - 1) T) Phi_all(1)   (w >= 1),

   since e^(-w t) <= e^(-(w - 1) T) e^(-t) for t > T.  A node uses the
   direct sum when the T needed for its accuracy is at most T0; T0 is
   chosen so that there are at most direct_max_points(d, prec) points,
   about where a direct sum (one real exponential per point) stops being
   cheaper than a theta value (measured: about 80 2^d at 1000 bits and
   20 2^d at 4000 bits). */

static slong
direct_max_points(slong d, slong prec)
{
    slong c = 80000 / FLINT_MAX(prec, 1);
    c = FLINT_MAX(c, 16);
    c = FLINT_MIN(c, 80);
    return c * (WORD(1) << d);
}

typedef struct
{
    slong nb;
    arb_ptr t;
    acb_ptr ph;
    double T0;
    mag_t phi1;
}
direct_struct;

static void
direct_init(direct_struct * D, const theta_side_struct * S, double tmind, slong prec)
{
    slong d = S->d, i, j, k, nb, nmax = direct_max_points(S->d, prec), it;
    acb_theta_eld_t E;
    arf_t R2;
    arb_ptr vc, v;
    arb_t q, s, one;
    double T, Tbest = 0.0;
    slong * pts;

    D->nb = 0;
    D->t = NULL;
    D->ph = NULL;
    D->T0 = 0.0;
    mag_init(D->phi1);

    arf_init(R2);
    vc = _arb_vec_init(d);
    v = _arb_vec_init(d);
    arb_init(q);
    arb_init(s);
    arb_init(one);

    arb_one(one);
    phi_all(D->phi1, S->C, one);

    /* |C n - C sg|^2 <= T */
    for (i = 0; i < d; i++)
        for (j = i; j < d; j++)
            arb_submul(vc + i, arb_mat_entry(S->C, i, j), S->sg + j, prec);

    /* the largest T (by factors 1.25) with at most nmax points */
    T = FLINT_MAX(tmind, 1e-3);
    for (it = 0; it < 200; it++)
    {
        acb_theta_eld_init(E, d, d);
        arf_set_d(R2, T);
        if (!acb_theta_eld_set(E, S->C, R2, vc) || acb_theta_eld_nb_pts(E) > nmax)
        {
            acb_theta_eld_clear(E);
            break;
        }
        acb_theta_eld_clear(E);
        Tbest = T;
        T *= 1.25;
    }

    if (Tbest > 0.0)
    {
        acb_theta_eld_init(E, d, d);
        arf_set_d(R2, Tbest);
        if (acb_theta_eld_set(E, S->C, R2, vc))
        {
            nb = acb_theta_eld_nb_pts(E);
            pts = flint_malloc(FLINT_MAX(nb, 1) * d * sizeof(slong));
            acb_theta_eld_points(pts, E);
            D->t = _arb_vec_init(FLINT_MAX(nb, 1));
            D->ph = _acb_vec_init(FLINT_MAX(nb, 1));

            for (k = 0; k < nb; k++)
            {
                const slong * n = pts + k * d;
                slong e;
                int skip = 0;

                for (e = 0; e < S->nexcl && !skip; e++)
                {
                    skip = 1;
                    for (j = 0; j < d; j++)
                        skip = skip && (n[j] == S->excl[e * d + j]);
                }
                if (skip)
                    continue;

                /* t_n = |C (n - sg)|^2 */
                arb_zero(q);
                for (i = 0; i < d; i++)
                {
                    arb_zero(s);
                    for (j = i; j < d; j++)
                    {
                        arb_sub_si(v + j, S->sg + j, n[j], prec);
                        arb_submul(s, arb_mat_entry(S->C, i, j), v + j, prec);
                    }
                    arb_addmul(q, s, s, prec);
                }
                arb_set(D->t + D->nb, q);

                /* e(-p.n) */
                arb_zero(q);
                for (j = 0; j < d; j++)
                    arb_submul_si(q, S->p + j, n[j], prec);
                arb_mul_2exp_si(q, q, 1);
                arb_sin_cos_pi(acb_imagref(D->ph + D->nb), acb_realref(D->ph + D->nb), q, prec);
                D->nb++;
            }

            D->T0 = Tbest;
            flint_free(pts);
        }
        acb_theta_eld_clear(E);
    }

    arf_clear(R2);
    _arb_vec_clear(vc, d);
    _arb_vec_clear(v, d);
    arb_clear(q);
    arb_clear(s);
    arb_clear(one);
}

static void
direct_clear(direct_struct * D)
{
    if (D->t != NULL)
    {
        _arb_vec_clear(D->t, FLINT_MAX(D->nb, 1));
        _acb_vec_clear(D->ph, FLINT_MAX(D->nb, 1));
    }
    mag_clear(D->phi1);
}

/* the T for which the tail e^(-(w - 1) T) Phi_all(1) is about 2^-lerr
   (a choice; w >= 2) */
static double
direct_T(const direct_struct * D, double wd, double lerr)
{
    return (lerr * 0.6931471805599453 + log(FLINT_MAX(mag_get_d(D->phi1), 1.0))) / (wd - 1.0);
}

/* Theta(w) by the direct sum with cutoff T, w >= 2 */
static void
direct_eval(acb_t res, const theta_side_struct * S, const direct_struct * D,
    const arb_t w, double T, slong prec)
{
    slong k, d = S->d, j;
    arb_t e, Tb;
    arf_t lb;
    acb_t ph;
    mag_t tail, m;

    arb_init(e);
    arb_init(Tb);
    arf_init(lb);
    acb_init(ph);
    mag_init(tail);
    mag_init(m);

    acb_zero(res);
    arb_set_d(Tb, T);
    for (k = 0; k < D->nb; k++)
    {
        arb_get_lbound_arf(lb, D->t + k, MAG_BITS);
        if (arf_cmp(lb, arb_midref(Tb)) > 0)
            continue;
        arb_mul(e, D->t + k, w, prec);
        arb_neg(e, e);
        arb_exp(e, e, prec);
        acb_addmul_arb(res, D->ph + k, e, prec);
    }

    /* the tail over the points with t_n > T: e^(-(w - 1) T) Phi_all(1) */
    arb_sub_ui(e, w, 1, MAG_BITS);
    arb_mul(e, e, Tb, MAG_BITS);
    arb_neg(e, e);
    arb_exp(e, e, MAG_BITS);
    arb_get_mag(tail, e);
    mag_mul(tail, tail, D->phi1);
    acb_add_error_mag(res, tail);

    /* the constant phase e(-p.pofs) */
    arb_zero(e);
    for (j = 0; j < d; j++)
        arb_submul(e, S->p + j, S->pofs + j, prec);
    arb_mul_2exp_si(e, e, 1);
    arb_sin_cos_pi(acb_imagref(ph), acb_realref(ph), e, prec);
    acb_mul(res, res, ph, prec);

    arb_clear(e);
    arb_clear(Tb);
    arf_clear(lb);
    acb_clear(ph);
    mag_clear(tail);
    mag_clear(m);
}

/* Gauss-Legendre nodes shared by the intervals and the two sides:
   the nodes of degree n at precision prec, for the most recent few
   degrees */
#define GL_CACHE_SIZE 4

typedef struct
{
    slong num;
    slong n[GL_CACHE_SIZE];
    slong prec[GL_CACHE_SIZE];
    arb_ptr x[GL_CACHE_SIZE];
    arb_ptr w[GL_CACHE_SIZE];
}
gl_cache_struct;

static void
gl_cache_init(gl_cache_struct * G)
{
    G->num = 0;
}

static void
gl_cache_clear(gl_cache_struct * G)
{
    slong i;
    for (i = 0; i < G->num; i++)
    {
        _arb_vec_clear(G->x[i], (G->n[i] + 1) / 2);
        _arb_vec_clear(G->w[i], (G->n[i] + 1) / 2);
    }
}

/* S = int_1^oo w^(a-1) Theta(w) dw to absolute accuracy ~2^-bits; returns
   0 if the bounds are useless (tmin = 0) */
/* Gauss-Legendre nodes x_k, k < ceil(n/2) (x_(n-1-k) = -x_k with the
   same weight), and the weighted integrand values w_k f(y_k), in
   parallel; the values are summed in a fixed order afterwards, so the
   result does not depend on the number of threads */
typedef struct
{
    arb_ptr x, w;
    slong n, prec;
}
gl_nodes_ctx;

static void
gl_nodes_worker(slong k, void * vc)
{
    gl_nodes_ctx * c = vc;
    arb_hypgeom_legendre_p_ui_root(c->x + k, c->w + k, c->n, k, c->prec);
}

/* the degree to use for at least n nodes: a cached degree in
   [n, n + n/8], else n */
static slong
gl_cache_degree(const gl_cache_struct * G, slong n)
{
    slong i, best = n;
    for (i = 0; i < G->num; i++)
        if (G->n[i] >= n && G->n[i] <= n + n / 8 && (best == n || G->n[i] < best))
            best = G->n[i];
    return best;
}

/* the nodes x_k, k < ceil(n/2), and weights of degree n to at least prec
   bits (computed with some extra bits, so that intervals needing
   slightly more precision share them) */
static void
gl_cache_get(arb_srcptr * x, arb_srcptr * w, gl_cache_struct * G, slong n, slong prec)
{
    slong i, j;
    gl_nodes_ctx c;

    for (i = 0; i < G->num; i++)
    {
        if (G->n[i] == n && G->prec[i] >= prec)
        {
            *x = G->x[i];
            *w = G->w[i];
            return;
        }
    }

    if (G->num == GL_CACHE_SIZE)
    {
        _arb_vec_clear(G->x[0], (G->n[0] + 1) / 2);
        _arb_vec_clear(G->w[0], (G->n[0] + 1) / 2);
        for (j = 1; j < G->num; j++)
        {
            G->n[j - 1] = G->n[j];
            G->prec[j - 1] = G->prec[j];
            G->x[j - 1] = G->x[j];
            G->w[j - 1] = G->w[j];
        }
        G->num--;
    }

    i = G->num;
    G->n[i] = n;
    G->prec[i] = prec + 32;
    G->x[i] = _arb_vec_init((n + 1) / 2);
    G->w[i] = _arb_vec_init((n + 1) / 2);
    G->num++;

    c.x = G->x[i];
    c.w = G->w[i];
    c.n = n;
    c.prec = G->prec[i];
    flint_parallel_do(gl_nodes_worker, &c, (n + 1) / 2, flint_get_num_threads(),
        FLINT_PARALLEL_STRIDED);

    *x = G->x[i];
    *w = G->w[i];
}

typedef struct
{
    const theta_side_struct * S;
    const direct_struct * D;
    arb_srcptr gl_x, gl_w, h, yc, tmin;
    acb_ptr vals;
    slong n, wp, bits;
}
gl_eval_ctx;

/* The value w_k f(y_k) of one node, needed to absolute accuracy about
   2^-bits / n.  Its precision is tapered (the enclosures are rigorous at
   any precision): with lmul = log2 |w_k e^(a y_k)|, Theta(w) is needed
   to absolute accuracy 2^-(bits + lmul + log2 n).  A direct sum has
   relative accuracy, so it needs about bits + lmul + log2 Phi(w) bits,
   which decreases rapidly for large w = e^y; acb_theta_one gives an
   absolute accuracy (relative to e^(-w sg^T Q sg) theta), so a theta
   value needs bits + lmul. */
static void
gl_eval_worker(slong k, void * vc)
{
    gl_eval_ctx * c = vc;
    slong kk = FLINT_MIN(k, c->n - 1 - k), wp = c->wp, extra;
    arb_t x, wt, y, w;
    acb_t f, ea;
    double lmul, lphi, wd, yd, T;
    int direct;

    arb_init(x);
    arb_init(wt);
    arb_init(y);
    arb_init(w);
    acb_init(f);
    acb_init(ea);

    /* the size of the multiplier and of Theta (low precision; choices) */
    {
        mag_t m;
        mag_init(m);
        arb_set_round(x, c->gl_x + kk, MAG_BITS);
        if (kk != k)
            arb_neg(x, x);
        arb_mul(y, x, c->h, MAG_BITS);
        arb_add(y, y, c->yc, MAG_BITS);
        arb_exp(w, y, MAG_BITS);
        phi_bound(m, c->S->C, c->tmin, w);
        yd = arf_get_d(arb_midref(y), ARF_RND_NEAR);
        wd = arf_get_d(arb_midref(w), ARF_RND_NEAR);
        lmul = arf_get_d(arb_midref(acb_realref(c->S->a)), ARF_RND_NEAR) * yd * 1.4426950408889634
            + log2(arf_get_d(arb_midref(c->gl_w + kk), ARF_RND_NEAR));
        lphi = mag_get_d_log2_approx(m);
        mag_clear(m);
        if (!(fabs(lmul) < 1e9))
            lmul = 0.0;
        if (!(fabs(lphi) < 1e9))
            lphi = 0.0;
    }
    extra = FLINT_BIT_COUNT(c->n) + 32;

    /* the direct sum when the cutoff it needs has precomputed points */
    T = 0.0;
    if (c->D->nb > 0 && wd >= 2.0)
        T = direct_T(c->D, wd, c->bits + lmul + FLINT_BIT_COUNT(c->n) + 16);
    direct = (T > 0.0 && T <= c->D->T0);

    if (direct)
        wp = FLINT_MIN(wp, (slong) (c->bits + lmul + FLINT_MIN(lphi, 0.0)) + extra);
    else
        wp = FLINT_MIN(wp, (slong) (c->bits + lmul) + extra);
    wp = FLINT_MAX(wp, 32);

    arb_set_round(x, c->gl_x + kk, wp);
    if (kk != k)
        arb_neg(x, x);
    arb_set_round(wt, c->gl_w + kk, wp);
    arb_mul(y, x, c->h, wp);
    arb_add(y, y, c->yc, wp);
    arb_exp(w, y, wp);

    if (direct)
        direct_eval(f, c->S, c->D, w, T, wp);
    else
        theta_eval(f, c->S, w, wp);
    /* e^(a y) */
    acb_mul_arb(ea, c->S->a, y, wp);
    acb_exp(ea, ea, wp);
    acb_mul(f, f, ea, wp);
    acb_mul_arb(c->vals + k, f, wt, wp);

    arb_clear(x);
    arb_clear(wt);
    arb_clear(y);
    arb_clear(w);
    acb_clear(f);
    acb_clear(ea);
}

static int
theta_side_sum(acb_t res, theta_side_struct * S, gl_cache_struct * G, slong bits)
{
    slong prec0 = 64, iv, nint, k;
    arb_t bmin, cosb;
    arb_t tmin, W, Y, ya, yc, r, t, u, x, wt, y, w, h, rho;
    acb_t f, sum, ea;
    mag_t M, err, m, tot_err;
    double tmind, Wd, rhod, are, aim, hd;
    arb_srcptr gl_x, gl_w;
    direct_struct D;
    int ok = 1;

    arb_init(tmin); arb_init(W); arb_init(Y); arb_init(ya); arb_init(yc); arb_init(r);
    arb_init(t); arb_init(u); arb_init(x); arb_init(wt); arb_init(y); arb_init(w);
    arb_init(h); arb_init(rho); arb_init(bmin); arb_init(cosb);
    acb_init(f); acb_init(sum); acb_init(ea);
    mag_init(M); mag_init(err); mag_init(m); mag_init(tot_err);
    D.t = NULL;
    mag_init(D.phi1);

    tmin_bound(tmin, S, prec0);
    tmind = arf_get_d(arb_midref(tmin), ARF_RND_DOWN);
    {
        arf_t ub;
        arf_init(ub);
        arb_get_ubound_arf(ub, acb_realref(S->a), prec0);
        are = arf_get_d(ub, ARF_RND_UP);
        arf_clear(ub);
    }
    aim = fabs(arf_get_d(arb_midref(acb_imagref(S->a)), ARF_RND_UP));
    if (!(tmind > 1e-6) || !(fabs(are) < 1e6) || !(aim < 1e6))
    {
        ok = 0;
        goto cleanup;
    }

    /* the tail: W with Phi_all(W/2) e^(-W tmin/2) W^(Re a - 1) / (tmin - (Re a - 1)/W)
       below 2^-(bits + 4); W found in doubles, the bound in balls */
    Wd = 2.0;
    while (1)
    {
        double lb = (Wd * tmind / 2) * 1.4426950408889634 - (are - 1) * log2(Wd)
            - d_log2_phi(S, Wd / 2) + log2(tmind);
        if (lb > bits + 8 || Wd > 1e12)
            break;
        Wd *= 1.25;
    }
    /* intervals of half-length hd <= HALF_LEN covering [0, log W] exactly */
    nint = (slong) ceil(log(Wd) / (2 * HALF_LEN));
    if (nint < 1)
        nint = 1;
    hd = log(Wd) / (2 * nint);
    hd = FLINT_MAX(hd, 0.25);
    arb_set_d(h, hd);
    arb_mul_ui(Y, h, 2 * nint, prec0);
    arb_exp(W, Y, prec0);
    {
        /* with s >= Re a (an exact upper bound) and Phi(w) <=
           e^(-(w - W/2) tmin) Phi_all(W/2) for w >= W/2:
           int_W^oo |w^(a-1)| Phi(w) dw <= Phi_all(W/2) e^(W tmin/2)
             W^(s - 1) e^(-W tmin) / (tmin - max(0, s - 1)/W), using
           w^(s-1) <= W^(s-1) e^((s-1)(w - W)/W) */
        arb_t c, su;
        mag_t pm;
        arb_init(c);
        arb_init(su);
        mag_init(pm);
        arb_get_ubound_arf(arb_midref(su), acb_realref(S->a), prec0);
        mag_zero(arb_radref(su));
        arb_mul_2exp_si(r, W, -1);
        phi_all(pm, S->C, r);
        arb_mul(t, r, tmin, prec0);
        arb_neg(t, t);
        arb_exp(t, t, prec0);
        arb_sub_ui(c, su, 1, prec0);
        arb_pow(u, W, c, prec0);
        arb_mul(t, t, u, prec0);
        if (arb_is_negative(c))
            arb_zero(c);
        arb_div(c, c, W, prec0);
        arb_sub(c, tmin, c, prec0);
        if (!arb_is_positive(c))
            mag_inf(tot_err);
        else
        {
            arb_div(t, t, c, prec0);
            arb_get_mag(tot_err, t);
            mag_mul(tot_err, tot_err, pm);
            /* |w^(a-1)| <= w^(Re a - 1) for real w */
        }
        arb_clear(c);
        arb_clear(su);
        mag_clear(pm);
    }

    /* the ellipse: semi-minor axis h (rho - 1/rho)/2 = THETA_MAX */
    rhod = THETA_MAX / hd;
    rhod = rhod + sqrt(rhod * rhod + 1.0);
    arb_set_d(rho, rhod);
    /* the actual semi-minor axis b = h (rho - 1/rho)/2 (rho is only a
       parameter), checked to be < pi/2, and cos(b) */
    {
        arb_t pi2;
        arb_init(pi2);
        arb_inv(u, rho, prec0);
        arb_sub(bmin, rho, u, prec0);
        arb_mul(bmin, bmin, h, prec0);
        arb_mul_2exp_si(bmin, bmin, -1);
        arb_const_pi(pi2, prec0);
        arb_mul_2exp_si(pi2, pi2, -1);
        if (!arb_lt(bmin, pi2))
            ok = 0;
        arb_cos(cosb, bmin, prec0);
        arb_clear(pi2);
    }

    /* the points for the direct sums, at the precision of the setup */
    mag_clear(D.phi1);
    direct_init(&D, S, tmind, bits + 64 + (slong) log2(Wd + 2));

    acb_zero(sum);
    for (iv = 0; iv < nint && ok; iv++)
    {
        slong n, wp;
        double Md, rmind;

        /* interval [ya, ya + 2h], center yc */
        arb_mul_ui(ya, h, 2 * iv, prec0);
        arb_add(yc, ya, h, prec0);

        /* M: |e^(a y)| <= e^(Re a Re y + |Im a| |Im y|) on the ellipse,
           times Phi(r_min) */
        arb_add_ui(t, rho, 0, prec0);
        arb_inv(u, rho, prec0);
        arb_add(t, t, u, prec0);
        arb_mul(t, t, h, prec0);
        arb_mul_2exp_si(t, t, -1);          /* semi-major */
        arb_sub(r, yc, t, prec0);
        arb_exp(r, r, prec0);
        arb_mul(r, r, cosb, prec0);         /* r_min (a ball; the bounds are monotone) */
        phi_bound(M, S->C, tmin, r);
        arb_add(u, yc, t, prec0);
        arb_sub(x, yc, t, prec0);
        arb_mul(u, u, acb_realref(S->a), prec0);
        arb_mul(x, x, acb_realref(S->a), prec0);
        arb_max(u, u, x, prec0);
        arb_abs(y, acb_imagref(S->a));
        arb_addmul(u, bmin, y, prec0);
        arb_exp(u, u, prec0);
        arb_get_mag(m, u);
        mag_mul(M, M, m);

        /* n: (64/15) h M rho^-2n / (rho^2 - 1) <= 2^-(bits + 4) / nint */
        Md = mag_get_d_log2_approx(M);
        rmind = log2(rhod);
        n = (slong) ceil((Md + (bits + 4) + log2(64.0 / 15 * hd / (rhod * rhod - 1)) + log2((double) nint)) / (2 * rmind));
        n = FLINT_MAX(n, 1);
        n = ((n + 7) / 8) * 8;
        if (n > 1000000)
        {
            ok = 0;
            break;
        }
        /* a cached degree slightly above n, if any: a few more values
           cost less than new nodes */
        n = gl_cache_degree(G, n);

        /* the error bound */
        arb_sqr(t, rho, prec0);
        arb_sub_ui(t, t, 1, prec0);
        arb_pow_ui(u, rho, 2 * n, prec0);
        arb_mul(t, t, u, prec0);
        arb_mul_ui(u, h, 64, prec0);
        arb_div_ui(u, u, 15, prec0);
        arb_div(u, u, t, prec0);
        arb_get_mag(err, u);
        mag_mul(err, err, M);
        mag_add(tot_err, tot_err, err);

        /* precision: absolute 2^-bits for values of size <= M */
        wp = bits + FLINT_MAX(0, (slong) Md) + FLINT_BIT_COUNT(n) + 16;

        /* the nodes, shared by the intervals and sides with the same n */
        gl_cache_get(&gl_x, &gl_w, G, n, wp);

        {
            gl_eval_ctx c;
            acb_t Sk;
            acb_init(Sk);
            c.S = S;
            c.D = &D;
            c.gl_x = gl_x;
            c.gl_w = gl_w;
            c.h = h;
            c.yc = yc;
            c.tmin = tmin;
            c.n = n;
            c.wp = wp;
            c.bits = bits;
            c.vals = _acb_vec_init(n);
            flint_parallel_do(gl_eval_worker, &c, n, flint_get_num_threads(),
                FLINT_PARALLEL_STRIDED);
            for (k = 0; k < n; k++)
                acb_add(Sk, Sk, c.vals + k, wp);
            _acb_vec_clear(c.vals, n);
            acb_mul_arb(Sk, Sk, h, wp);
            acb_add(sum, sum, Sk, wp);
            acb_clear(Sk);
        }
    }

    if (ok)
    {
        acb_set(res, sum);
        acb_add_error_mag(res, tot_err);
    }

cleanup:
    direct_clear(&D);
    arb_clear(tmin); arb_clear(W); arb_clear(Y); arb_clear(ya); arb_clear(yc); arb_clear(r);
    arb_clear(t); arb_clear(u); arb_clear(x); arb_clear(wt); arb_clear(y); arb_clear(w);
    arb_clear(h); arb_clear(rho); arb_clear(bmin); arb_clear(cosb);
    acb_clear(f); acb_clear(sum); acb_clear(ea);
    mag_clear(M); mag_clear(err); mag_clear(m); mag_clear(tot_err);
    return ok;
}

static void
side_setup(theta_side_struct * S, const _acb_mzeta_epstein_problem_struct * P,
    int side, const arb_mat_t C, slong prec)
{
    slong d = P->d, i, j;
    arb_mat_t CT;

    S->P = P;
    S->side = side;
    S->d = d;
    S->a = side ? P->b : P->a;
    S->Q = flint_malloc(sizeof(arb_mat_struct));
    S->C = flint_malloc(sizeof(arb_mat_struct));
    arb_mat_init(S->Q, d, d);
    arb_mat_init(S->C, d, d);
    arb_mat_init(CT, d, d);
    arb_mat_set(S->C, C);
    arb_mat_transpose(CT, C);
    arb_mat_mul(S->Q, CT, C, prec);
    arb_mat_clear(CT);
    S->sg = _arb_vec_init(d);
    S->p = _arb_vec_init(d);
    S->pofs = _arb_vec_init(d);
    S->Qsg = _arb_vec_init(d);
    arb_init(S->sQs);

    if (side == 0)
    {
        /* t = |Cd (n - u0)|^2, phase e(-w0.n) */
        _arb_vec_set(S->sg, P->u0, d);
        _arb_vec_set(S->p, P->w0, d);
    }
    else
    {
        /* t = |Cr (m + w0)|^2, phase e(-u0.(m + w0)) */
        _arb_vec_neg(S->sg, P->w0, d);
        _arb_vec_set(S->p, P->u0, d);
        _arb_vec_set(S->pofs, P->w0, d);
    }
    for (i = 0; i < d; i++)
    {
        for (j = 0; j < d; j++)
            arb_addmul(S->Qsg + i, arb_mat_entry(S->Q, i, j), S->sg + j, prec);
        arb_addmul(S->sQs, S->Qsg + i, S->sg + i, prec);
    }

    /* the excluded points */
    S->excl = flint_calloc(2 * d, sizeof(slong));
    S->nexcl = 0;
    if ((side == 0 && P->xs == 1) || (side == 1 && P->ys == 1))
        S->nexcl++;    /* n = 0 */
    if (side == 1 && P->reg)
    {
        for (j = 0; j < d; j++)
            S->excl[S->nexcl * d + j] = P->lpt[j];
        S->nexcl++;
    }
}

static void
side_clear(theta_side_struct * S)
{
    slong d = S->d;
    arb_mat_clear(S->Q);
    arb_mat_clear(S->C);
    flint_free(S->Q);
    flint_free(S->C);
    _arb_vec_clear(S->sg, d);
    _arb_vec_clear(S->p, d);
    _arb_vec_clear(S->pofs, d);
    _arb_vec_clear(S->Qsg, d);
    arb_clear(S->sQs);
    flint_free(S->excl);
}

int
acb_mzeta_epstein_theta(acb_t res, const acb_t nu, const arb_mat_t A,
    arb_srcptr x, arb_srcptr y, int flags, acb_mzeta_cache_t cache, slong prec)
{
    slong d = arb_mat_nrows(A), wp;
    acb_mzeta_cache_t tmp_cache;
    const _acb_mzeta_lattice_struct * L;
    _acb_mzeta_epstein_problem_struct P;
    theta_side_struct S0, S1;
    acb_t Sd, Sr;
    gl_cache_struct G;
    int ok;

    if (d < 1 || d > 8 || arb_mat_ncols(A) != d)
    {
        acb_indeterminate(res);
        return 0;
    }

    if (cache == NULL)
    {
        acb_mzeta_cache_init(tmp_cache);
        ok = acb_mzeta_epstein_theta(res, nu, A, x, y, flags, tmp_cache, prec);
        acb_mzeta_cache_clear(tmp_cache);
        return ok;
    }

    /* the errors of the input to the theta function are amplified by
       about w^(Re a) relative to the integrand, w up to ~prec / tmin:
       extra precision for the setup (accuracy only) */
    {
        double ra = arf_get_d(arb_midref(acb_realref(nu)), ARF_RND_NEAR) / 2;
        double mx = FLINT_MAX(0.0, FLINT_MAX(ra, d / 2.0 - ra));
        if (!(mx < 1e6))
            mx = 0.0;
        wp = prec + 32 + (slong) (mx * (log2((double) prec + 2) + 12));
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

    acb_init(Sd);
    acb_init(Sr);
    side_setup(&S0, &P, 0, L->Cd, wp);
    side_setup(&S1, &P, 1, L->Cr, wp);

    gl_cache_init(&G);
    ok = theta_side_sum(Sd, &S0, &G, prec + 10) && theta_side_sum(Sr, &S1, &G, prec + 10);
    gl_cache_clear(&G);
    if (ok)
        _acb_mzeta_epstein_finish(res, Sd, Sr, L, &P, prec, wp);
    else
        acb_indeterminate(res);

    side_clear(&S0);
    side_clear(&S1);
    acb_clear(Sd);
    acb_clear(Sr);
    _acb_mzeta_epstein_problem_clear(&P);
    return 1;
}
