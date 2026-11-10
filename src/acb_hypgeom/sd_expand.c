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
    Saddle-point expansions for the integrals of exp(phi(w)) along the
    steepest-descent contours of sd_plan.c, as an alternative to the
    quadrature of sd_integrate.c when the saddles are well separated.

    Let w_s be a simple saddle, c = phi(w_s), A = phi''(w_s).  Since
    phi''' = phi' + nu and phi'''' = phi'' for both kernels,

        phi(w_s + u) - c = P(u) = A (cosh u - 1) + nu (sinh u - u).

    Along the two descent paths from the saddle, P(u) = -tau^2 with tau
    real, and u(tau) = sum_{k >= 1} b_k tau^k is an analytic function of
    tau near 0 (b_1 = sqrt(-2/A), with the sign of the outgoing path;
    the incoming path is tau < 0).  The integral from one valley to the
    other is

        e^c int e^(-tau^2) u'(tau) dtau ~ e^c sum_{k odd} k b_k Gamma(k/2),

    the expansion of Debye type (for the Bessel functions, this is the
    Debye expansion, summed numerically for the given parameters).  The
    coefficients follow from D u' = -2 tau, D = P'(u) = A sinh u +
    nu (cosh u - 1), together with (sinh u)' = cosh u u' and
    (cosh u)' = sinh u u', by a recurrence with O(n^2) operations for n
    coefficients.

    Error bound.  Write P(u) = (A/2) u^2 (1 + g(u)); for |u| <= r,

        |g(u)| <= G(r) = 2 (cosh r - 1 - r^2/2) / r^2
                         + 2 |nu/A| (sinh r - r) / r^2.

    If G(r) < 1, then s(u) = u sqrt(1 + g(u)) satisfies |s(u) - u| <=
    r (1 - sqrt(1 - G)) on |u| = r, so by Rouche's theorem s takes each
    value sigma with |sigma| < r sqrt(1 - G) exactly once in |u| < r.
    Since tau = sqrt(-A/2) s(u), the function u(tau) is analytic with
    |u| < r in the disk |tau| < rho = sqrt(|A|/2) r sqrt(1 - G(r)), and
    the descent paths stay on this branch there.  By Cauchy's
    inequality |b_k| <= r rho^(-k).  For tau_0 < rho (on either side),

        int_0^{tau_0} e^(-tau^2) |u' - sum_{k <= n} k b_k tau^(k-1)| dtau
            <= r sum_{k > n} min(k Gamma(k/2) / 2, tau_0^k) rho^(-k),

        sum_{k <= n} k |b_k| Gamma(k/2, tau_0^2) / 2
            <= r e^(-tau_0^2) tau_0^(-2) t / (1 - t)^2,   t = tau_0 / rho

    (for tau_0^2 >= n - 2), and the integral over the rest of the
    descent path, beyond u(tau_0), is that along the segment from
    w_s + u(tau_0) (enclosed with the Taylor polynomial and the bound
    r t^(n+1) / (1 - t) for its remainder) to a vertex of the planned
    polygon followed by the rest of the polygon and its tail, which are
    bounded by the majorants of the planner.  The first sum decreases
    like the terms of the expansion until k ~ 2 rho^2, so that the
    attainable accuracy is about e^(-rho^2), where rho^2 is typically
    within a factor 0.6 to 0.9 of the distance from c to the nearest
    other critical value of phi.

    Since the integrand is entire, the result is an enclosure whatever
    the polygon is; only the size of the last bound depends on it.

    This is implemented for the symmetric kernel alpha (e^w - e^-w) -
    nu w with exact parameters, where the saddles, A and c are given
    by square roots and logarithms.  As for the quadrature, the
    tolerance is 2^-prec times the L1 norm of the integrand on the
    whole contour: the expansion of a piece with a small contribution
    needs fewer terms and less precision, and a piece estimated to be
    below the tolerance is left to the quadrature.

    What is rigorous and what is not.  Whether the expansions are used,
    and their parameters (r, the number of terms n, tau_0 and the
    vertices where the polygon takes over), are chosen heuristically in
    double precision (sde_prepare and the functions before it) so that
    the error bound comes out below the tolerance.  The bound itself is
    then computed in arb and mag arithmetic from enclosures of the
    saddle, A and c (sde_bound_radius, sde_bound_truncation,
    sde_bound_rest, and sd_bound.c for the path beyond tau_0), and is
    valid for any choice of the parameters.
*/

#include <math.h>
#include "acb.h"
#include "acb_hypgeom/sd.h"

/* The coefficients b_2, ..., b_n from b_1, where A b_1^2 = -2.  With
   s = sinh u, h = cosh u, D = A s + nu (h - 1) (as power series in tau)
   and [k] denoting a coefficient,

       k [k] s = sum_{i=1}^{k} i b_i [k-i] h,
       k [k] h = sum_{i=1}^{k-1} i b_i [k-i] s,
       sum_{i=1}^{k} [i] D (k - i + 1) b_{k-i+1} = -2 [k = 1],

   where b_k appears in [k] s = S + b_k and [k] D = A S + nu [k] h
   + A b_k with the sum S over i < k, and in the term i = 1 of the last
   sum, so that (k + 1) A b_1 b_k = -(known terms) and
   b_k = (known terms) b_1 / (2 (k + 1)).

   The terms of the expansion decrease by a factor of about 2^(-prec)
   over the n coefficients, so the working precision decreases along
   the recurrence (wp for the first coefficients). */
static void
sde_coeffs(acb_ptr b, const acb_t A, const acb_t nu, slong n, slong prec, slong wp)
{
    acb_ptr ib, s, h, D;
    acb_t S, H, Dn, t;
    slong k, p;

    ib = _acb_vec_init(n + 2); s = _acb_vec_init(n + 2);
    h = _acb_vec_init(n + 2); D = _acb_vec_init(n + 2);
    acb_init(S); acb_init(H); acb_init(Dn); acb_init(t);

    acb_set(ib + 1, b + 1);
    acb_set(s + 1, b + 1);
    acb_one(h);
    acb_mul(D + 1, A, b + 1, wp);

    for (k = 2; k <= n; k++)
    {
        p = wp - (slong) (0.7 * (double) prec * (double) k / (double) n);
        p = FLINT_MAX(p, 40);

        acb_dot(S, NULL, 0, ib + 1, 1, h + k - 1, -1, k - 1, p);
        acb_dot(H, NULL, 0, ib + 1, 1, s + k - 1, -1, k - 1, p);
        acb_div_ui(S, S, k, p);
        acb_div_ui(H, H, k, p);
        acb_mul(Dn, A, S, p);
        acb_addmul(Dn, nu, H, p);
        acb_mul(t, Dn, b + 1, p);
        acb_dot(t, t, 0, D + 2, 1, ib + k - 1, -1, k - 2, p);
        acb_mul(t, t, b + 1, p);
        acb_div_ui(b + k, t, 2 * (k + 1), p);
        acb_mul_ui(ib + k, b + k, k, p);
        acb_add(s + k, S, b + k, p);
        acb_swap(h + k, H);
        acb_addmul(Dn, A, b + k, p);
        acb_swap(D + k, Dn);
    }

    _acb_vec_clear(ib, n + 2); _acb_vec_clear(s, n + 2);
    _acb_vec_clear(h, n + 2); _acb_vec_clear(D, n + 2);
    acb_clear(S); acb_clear(H); acb_clear(Dn); acb_clear(t);
}

/* The plan for one piece.  Heuristic: everything from here to
   sde_prepare is done in double precision to choose the parameters of
   the expansion (the radius r, the number of terms, the points tau_0
   where the expansion is cut off) so that the rigorous error bound,
   which is computed afterwards from these parameters (sde_bound_*),
   meets the tolerance. */
typedef struct
{
    int n;                /* number of coefficients b_1, ..., b_n, or 0
                             for a piece estimated to be below the
                             tolerance, which is left to the quadrature
                             (-1: not yet checked) */
    double r, rho;        /* |u| < r for |tau| < rho (estimate) */
    int j[2];             /* the vertices beyond which the polygon is only bounded */
    double tau0[2];
    int len;              /* the terms of the error bound up to this
                             index are bounded individually */
    double logtol;        /* log of the tolerance, relative to e^c */
    slong prec;           /* relative precision needed for this piece */
    double cost;
}
sde_plan;

/* the factor 1 + 2e-6 allows for the difference between phi'' at the
   double saddle and at the true saddle (checked by the caller) */
#define SDE_SLACK 2e-6

static double
sde_G(double r, double q)
{
    double a, b;

    if (r < 0.3)
    {
        /* 2 (cosh r - 1 - r^2/2) / r^2 = r^2/12 + r^4/360 + ... */
        a = r * r / 12.0 * (1.0 + r * r / 29.0);
        /* 2 (sinh r - r) / r^2 = r/3 + r^3/60 + ... */
        b = r / 3.0 * (1.0 + r * r / 19.0);
    }
    else
    {
        a = 2.0 * (cosh(r) - 1.0 - 0.5 * r * r) / (r * r);
        b = 2.0 * (sinh(r) - r) / (r * r);
    }

    return (a + q * b) * (1.0 + 1e-9);
}

/* the largest rho (rounded down) over a range of r, and its r */
static double
sde_radius(double * rbest, double absA, double q)
{
    double r, G, rho, best = 0.0;

    *rbest = 0.0;
    for (r = 0.01; r < 4.0; r *= 1.2)
    {
        G = sde_G(r, q);
        if (!(G < 1.0))
            break;
        rho = sqrt(0.5 * absA) * r * sqrt(1.0 - G) * (1.0 - 1e-9);
        if (rho > best)
        {
            best = rho;
            *rbest = r;
        }
    }

    return best;
}

/* log of sum_{k > n} min(k Gamma(k/2) / 2, tau0^k) rho^(-k), rounded up,
   for all n < len: tab[n]; requires tau0 < rho */
static void
sde_tail_table(double * tab, int len, double tau0, double rho)
{
    double lt = log(tau0), lr = log(rho), t = tau0 / rho, s, lg[2];
    int k;

    /* beyond len: geometric */
    s = (len + 1) * (lt - lr) - log1p(-t);

    /* log Gamma(k/2 + 1) for k and k - 1, downwards by the recurrence */
    lg[0] = lgamma(0.5 * len + 1.0);
    lg[1] = lgamma(0.5 * (len - 1) + 1.0);

    for (k = len; k >= 1; k--)
    {
        double a = lg[(len - k) & 1], b = k * lt;
        s = _sd_logaddexp(s, fmin(a, b) - k * lr);
        tab[k - 1] = s + 1e-9 * (fabs(s) + 1.0);
        lg[(len - k) & 1] = a - log(0.5 * k);
    }
}

/* number of coefficients such that the truncation error is below
   exp(logtol), or -1; tau0 < 0.9 rho */
static int
sde_choose_n(double * logerr, int * lenout, double r, double rho, double tau0, double logtol, int nmax)
{
    double * tab, lr = log(r);
    int len, n, res = -1;

    /* beyond len, the terms are bounded by the geometric series, which
       should be negligible */
    len = (int) fmin(4.0 * tau0 * tau0 + 64.0, 1e6);
    len = FLINT_MIN(len, nmax + 64);
    len = (int) fmin(len, (logtol - 3.0 - lr + log1p(-tau0 / rho)) / log(tau0 / rho) + 2.0);
    len = FLINT_MAX(len, 2);
    *lenout = len;
    tab = flint_malloc(sizeof(double) * (len + 1));
    sde_tail_table(tab, len, tau0, rho);

    for (n = 1; n < len && n <= nmax; n++)
    {
        if (lr + tab[n] <= logtol)
        {
            res = n;
            *logerr = lr + tab[n];
            break;
        }
    }

    flint_free(tab);
    return res;
}

/* Cost of the expansion for one saddle (microseconds on a typical
   machine, as measured): the saddle and the bounds, and the recurrence
   with 3 n^2 / 2 complex products in dot products at the working
   precision wp.  With real parameters, the coefficients are real or
   imaginary. */
static double
sde_cost(int n, slong wp, int real)
{
    double p = (double) wp, m = (double) n, c2;

    c2 = fmax(fmin(0.09 + 0.0011 * (p - 80.0), 0.2), 0.1 + 0.0005 * p + 5e-8 * p * p);

    return (real ? 6.0 : 8.0) + (real ? 0.06 : 0.08) * p + m * (1.0 + 0.002 * p)
        + (real ? 0.6 : 1.0) * m * m * c2;
}

/* Cost of the quadrature for one piece (two descent paths) with the
   given number of bits in the arithmetic, likewise: in N-fold double
   precision up to N = 4. */
double
_sd_quad_cost(double bits)
{
    static const double dcost[4] = { 30.0, 45.0, 100.0, 200.0 };

    if (bits <= 4 * 53)
        return dcost[(int) ((fmax(bits, 2.0) - 1.0) / 53.0)];
    return 650.0 * pow(bits / 256.0, 2.1);
}

double
_sd_expand_cost(double absA, double q, int real, slong prec)
{
    double r, rho, tau0, logerr;
    int n, len;

    if (!(absA > 0.0) || !(q < 1e300))
        return INFINITY;

    rho = sde_radius(&r, absA * (1.0 - SDE_SLACK), q * (1.0 + SDE_SLACK));
    tau0 = sqrt((prec + 20.0) * SD_LN2);
    if (!(tau0 < 0.85 * rho))
        return INFINITY;

    n = sde_choose_n(&logerr, &len, r, rho, tau0, -(prec + 8.0) * SD_LN2, 100000);
    if (n < 0)
        return INFINITY;
    if (n > tau0 * tau0)
    {
        tau0 = sqrt(n + 2.0);
        if (!(tau0 < 0.85 * rho))
            return INFINITY;
    }

    return sde_cost(n, prec + 40, real);
}

/* Re(phi - top) at a vertex */
static double
sde_rephi(const sd_problem * P, cd w, double top)
{
    cd phi, d1, d2;
    _sd_phi(&phi, &d1, &d2, P, w);
    return phi.re - top;
}

static int
sde_prepare(sde_plan * pl, const sd_problem * P, const sd_piece * pc, double logL,
    slong prec, slong wp, int check, double maxcost)
{
    cd phi, d1, d2, ws;
    double absA, q, tol, lead, T, top = pc->c.re, lE1, lE2, g[2];
    int side, iter, mid = pc->mid, npts = pc->npts, n = 1, real;

    /* A piece below the tolerance (as for the quadrature, 2^-prec times
       the L1 norm of the integrand on the whole contour) is only
       bounded.  This is only checked where the integrand is small at
       the top of the piece. */
    tol = logL - top - (prec + 5.0) * SD_LN2;
    if (tol > -12.0 && !check)
    {
        /* to be checked if the expansions are used */
        pl->n = -1;
        pl->cost = 15.0;
        return 1;
    }
    if (tol > -12.0)
    {
        lE1 = _sd_logaddexp(pc->logtail[0], pc->logtail[1]);
        for (iter = 0; iter + 1 < npts; iter++)
            lE1 = _sd_logaddexp(lE1, _sd_panel(g, P, pc->pts[iter], pc->pts[iter + 1], top));
        if (lE1 <= tol)
        {
            pl->n = 0;
            pl->prec = 2;
            pl->cost = 15.0;
            return 1;
        }
    }

    if (mid < 1 || mid > npts - 2)
        return 0;

    ws = pc->pts[mid];
    _sd_phi(&phi, &d1, &d2, P, ws);
    absA = cd_abs(d2);
    if (!(absA > 1e-300) || !(absA < 1e300))
        return 0;
    q = P->anu / absA;

    pl->rho = sde_radius(&pl->r, absA * (1.0 - SDE_SLACK), q * (1.0 + SDE_SLACK));
    if (!(pl->rho > 0.0))
        return 0;

    /* The tolerance, relative to e^c.  The first term for this piece
       is b_1 sqrt(pi), which gives the relative precision needed for
       it. */
    lead = 0.5 * log(2.0 / absA) + 0.5 * log(SD_PI);
    tol = fmax(lead, logL - top) - (prec + 5.0) * SD_LN2;
    tol = fmin(tol, lead - 4.0 * SD_LN2);
    pl->prec = (slong) ((lead - tol) / SD_LN2) + 1;
    T = -tol + 4.0;
    real = (P->nu.im == 0.0) && (P->al.re == 0.0 || P->al.im == 0.0);

    /* A quick exit if the cost is certainly too high: the bound for
       the first term omitted, r (k/2) Gamma(k/2) rho^-k for k = n + 1
       (as tau_0^2 >= n - 1), must be below the tolerance. */
    if (maxcost < 1e300)
    {
        double lr = log(pl->rho), lr0 = log(pl->r);

        for (n = 1; lgamma(0.5 * (n + 1) + 1.0) - (n + 1) * lr + lr0 > tol; n++)
            if (n > 100000 || sde_cost(n + 1, wp - (prec - pl->prec), real) > maxcost)
                return 0;
    }

    for (iter = 0; iter < 4; iter++)
    {
        double tmax = 0.0, tmin = INFINITY;
        int nn;

        /* the first vertices where the integrand has dropped by e^-T */
        for (side = 0; side < 2; side++)
        {
            int jj, step = side ? 1 : -1;
            double v = 0.0;

            for (jj = mid + step; jj >= 0 && jj < npts; jj += step)
            {
                v = sde_rephi(P, pc->pts[jj], top);
                if (v <= -T)
                    break;
            }
            if (jj < 0 || jj >= npts)
                return 0;

            pl->j[side] = jj;
            pl->tau0[side] = sqrt(-v);
            tmax = fmax(tmax, pl->tau0[side]);
            tmin = fmin(tmin, pl->tau0[side]);
        }

        if (!(tmax < 0.85 * pl->rho))
            return 0;

        /* each of the two sides is bounded with the larger tau0 */
        nn = sde_choose_n(&lE2, &pl->len, pl->r, pl->rho, tmax, tol - log(8.0), 200000);
        if (nn < 0)
            return 0;
        n = nn;

        if (tmin * tmin >= n - 1.0)
            break;
        T = n + 2.0;
        if (iter == 3)
            return 0;
    }

    pl->n = n;
    pl->logtol = tol - log(8.0);

    pl->cost = sde_cost(n, wp - (prec - pl->prec), real);
    return 1;
}

/* The saddle near the double w0, A = phi''(w_s) and c = phi(w_s) for
   the kernel alpha (e^w - e^-w) - nu w: with phi' = alpha (e^w + e^-w)
   - nu = 0 and A = alpha (e^w - e^-w),

       A^2 = nu^2 - 4 alpha^2,   e^w = (nu + A) / (2 alpha)
           = 2 alpha / (nu - A),   c = A - nu w,

   with the sign of A and the branch of the logarithm of the saddle
   closest to w0 (the other saddles are -w0 + 2 pi i k, and both the
   expansion and the choice between them are only valid when A is not
   small, which is checked by the caller).  Returns 0 on failure. */
static int
sde_saddle(acb_t w, acb_t A, acb_t c, const acb_t alpha, const acb_t nu,
    cd w0, cd A0, slong wp)
{
    acb_t t, u;
    arb_t pi;
    double ar, ai, wi, kk;
    int ok = 0;

    acb_init(t); acb_init(u);
    arb_init(pi);

    acb_mul(t, nu, nu, wp);
    acb_mul(u, alpha, alpha, wp);
    acb_mul_2exp_si(u, u, 2);
    acb_sub(t, t, u, wp);
    acb_sqrt(A, t, wp);

    ar = arf_get_d(arb_midref(acb_realref(A)), ARF_RND_NEAR);
    ai = arf_get_d(arb_midref(acb_imagref(A)), ARF_RND_NEAR);
    if (ar * A0.re + ai * A0.im < 0.0)
        acb_neg(A, A);

    /* without cancellation */
    acb_mul_2exp_si(u, alpha, 1);
    acb_add(t, nu, A, wp);
    acb_sub(w, nu, A, wp);
    if (arf_cmpabs(arb_midref(acb_realref(t)), arb_midref(acb_realref(w))) >= 0
        && arf_cmpabs(arb_midref(acb_imagref(t)), arb_midref(acb_imagref(w))) >= 0)
        acb_div(t, t, u, wp);
    else
        acb_div(t, u, w, wp);

    if (!acb_is_finite(t) || acb_contains_zero(t))
        goto cleanup;

    acb_log(w, t, wp);
    wi = arf_get_d(arb_midref(acb_imagref(w)), ARF_RND_NEAR);
    kk = floor((w0.im - wi) / (2.0 * SD_PI) + 0.5);
    if (!(fabs(kk) < 1e15))
        goto cleanup;
    if (kk != 0.0)
    {
        arb_const_pi(pi, wp);
        arb_set_d(acb_realref(t), 2.0 * kk);
        arb_addmul(acb_imagref(w), pi, acb_realref(t), wp);
    }

    /* the saddle of the contour */
    arb_set_d(acb_realref(t), w0.re);
    arb_set_d(acb_imagref(t), w0.im);
    acb_sub(t, t, w, 60);
    acb_mul_2exp_si(t, t, 20);
    arb_set_d(acb_realref(u), 1.0 + fabs(w0.re) + fabs(w0.im));
    arb_zero(acb_imagref(u));
    acb_get_mag(arb_radref(acb_realref(t)), t);
    acb_get_mag_lower(arb_radref(acb_imagref(t)), u);
    if (mag_cmp(arb_radref(acb_realref(t)), arb_radref(acb_imagref(t))) > 0)
        goto cleanup;

    acb_mul(c, nu, w, wp);
    acb_sub(c, A, c, wp);
    ok = acb_is_finite(c);

cleanup:
    acb_clear(t); acb_clear(u);
    arb_clear(pi);
    return ok;
}

/* Rigorous.  rho = a lower bound for sqrt(|A| / 2) r sqrt(1 - G(r)),

       G(r) = 2 (cosh r - 1 - r^2/2) / r^2 + 2 |nu / A| (sinh r - r) / r^2,

   for the exact number r > 0 and the balls A and nu; returns 0 if G(r)
   is not certified to be less than 1. */
static int
sde_bound_radius(arb_t rho, const acb_t A, const acb_t nu, double r)
{
    arb_t R, G, sh, ch, t;
    mag_t q, a;
    int ok = 0;

    arb_init(R); arb_init(G); arb_init(sh); arb_init(ch); arb_init(t);
    mag_init(q); mag_init(a);

    acb_get_mag(q, nu);
    acb_get_mag_lower(a, A);

    if (r > 0.0 && r < 100.0 && !mag_is_zero(a) && mag_is_finite(q))
    {
        mag_div(q, q, a);                       /* |nu / A| */

        arb_set_d(R, r);
        arb_sinh_cosh(sh, ch, R, 128);
        arb_mul(t, R, R, 128);                  /* r^2 */
        arb_sub(sh, sh, R, 128);                /* sinh r - r */
        arb_sub_ui(ch, ch, 1, 128);
        arb_mul_2exp_si(G, t, -1);
        arb_sub(ch, ch, G, 128);                /* cosh r - 1 - r^2 / 2 */
        /* upper bounds, before multiplying by the upper bound of q */
        arb_get_ubound_arf(arb_midref(sh), sh, 128);
        mag_zero(arb_radref(sh));
        arf_set_mag(arb_midref(G), q);
        mag_zero(arb_radref(G));
        arb_mul(G, G, sh, 128);
        arb_add(G, G, ch, 128);
        arb_div(G, G, t, 128);
        arb_mul_2exp_si(G, G, 1);

        arb_sub_ui(G, G, 1, 128);
        arb_neg(G, G);                          /* 1 - G */
        if (arb_is_positive(G))
        {
            arb_sqrt(G, G, 128);
            arf_set_mag(arb_midref(t), a);
            mag_zero(arb_radref(t));
            arb_mul_2exp_si(t, t, -1);
            arb_sqrt(t, t, 128);
            arb_mul(G, G, t, 128);
            arb_mul(G, G, R, 128);
            arb_get_lbound_arf(arb_midref(rho), G, 128);
            mag_zero(arb_radref(rho));
            ok = arb_is_finite(rho) && arb_is_positive(rho);
        }
    }

    arb_clear(R); arb_clear(G); arb_clear(sh); arb_clear(ch); arb_clear(t);
    mag_clear(q); mag_clear(a);
    return ok;
}

/* Rigorous.  For one side, with the exact numbers r, tau_0 > 0 and a
   lower bound rho for the radius of analyticity (t = tau_0 / rho):

       err >= r sum_{k > n} min((k/2) Gamma(k/2), tau_0^k) rho^(-k)
              + r e^(-tau_0^2) tau_0^(-2) t / (1 - t)^2,
       rem >= r t^(n+1) / (1 - t),

   the bounds for the truncation of the expansion (the first sum, with
   the terms up to k = len bounded individually and the rest by the
   geometric series) and of the Taylor polynomial of u at tau_0.
   Returns 0 unless t < 1 and tau_0^2 >= n - 2 are certified. */
static int
sde_bound_truncation(mag_t err, mag_t rem, const arb_t rho, double r, double tau0, slong n, slong len)
{
    arb_t T, t, omt, u;
    mag_struct g[2];
    mag_t pw, tp, tm, irho, sum, term, rr, iomt;
    slong k;
    int ok = 0;

    arb_init(T); arb_init(t); arb_init(omt); arb_init(u);
    mag_init(g); mag_init(g + 1); mag_init(pw); mag_init(tp); mag_init(tm);
    mag_init(irho); mag_init(sum); mag_init(term); mag_init(rr); mag_init(iomt);

    arb_set_d(T, tau0);
    arb_div(t, T, rho, 128);
    arb_sub_ui(omt, t, 1, 128);
    arb_neg(omt, omt);                          /* 1 - t */
    arb_mul(u, T, T, 128);                      /* tau_0^2 */
    arb_sub_si(T, u, n - 2, 128);

    if (tau0 > 0.0 && arb_is_positive(omt) && arb_is_nonnegative(T) && r > 0.0)
    {
        len = FLINT_MAX(len, n + 1);

        mag_set_d(rr, r);
        arb_get_mag(tm, t);
        arb_inv(T, rho, 128);
        arb_get_mag(irho, T);
        arb_inv(T, omt, 128);
        arb_get_mag(iomt, T);

        /* Gamma(k/2 + 1) for k = n + 1, n + 2 */
        for (k = 0; k < 2; k++)
        {
            arb_set_ui(T, n + 1 + k + 2);
            arb_mul_2exp_si(T, T, -1);
            arb_gamma(T, T, 64);
            arb_get_mag(g + k, T);
        }

        mag_pow_ui(pw, irho, n + 1);
        mag_pow_ui(tp, tm, n + 1);

        /* rem */
        mag_mul(rem, tp, iomt);
        mag_mul(rem, rem, rr);

        for (k = n + 1; k <= len; k++)
        {
            int par = (k - (n + 1)) & 1;

            mag_mul(term, g + par, pw);
            mag_min(term, term, tp);
            mag_add(sum, sum, term);

            /* Gamma((k + 2)/2 + 1) = (k + 2)/2 Gamma(k/2 + 1) */
            mag_mul_ui(g + par, g + par, k + 2);
            mag_mul_2exp_si(g + par, g + par, -1);
            mag_mul(pw, pw, irho);
            mag_mul(tp, tp, tm);
        }

        /* sum_{k > len} t^k */
        mag_mul(term, tp, iomt);
        mag_add(sum, sum, term);

        /* e^(-tau_0^2) tau_0^(-2) t / (1 - t)^2 */
        arb_neg(T, u);
        arb_exp(T, T, 128);
        arb_div(T, T, u, 128);
        arb_mul(T, T, t, 128);
        arb_div(T, T, omt, 128);
        arb_div(T, T, omt, 128);
        arb_get_mag(term, T);
        mag_add(sum, sum, term);

        mag_mul(err, sum, rr);
        ok = mag_is_finite(err) && mag_is_finite(rem);
    }

    arb_clear(T); arb_clear(t); arb_clear(omt); arb_clear(u);
    mag_clear(g); mag_clear(g + 1); mag_clear(pw); mag_clear(tp); mag_clear(tm);
    mag_clear(irho); mag_clear(sum); mag_clear(term); mag_clear(rr); mag_clear(iomt);
    return ok;
}

/* Rigorous.  err >= the integral of |exp(phi - c)| over a path from the
   ball W to the valley at the given end of the piece: the segment
   from W to the vertex j, the polygon from there to its end, and the
   ray into the valley.  Returns 0 if the ray is not certified. */
static int
sde_bound_rest(mag_t err, const sd_problem * P, const sd_piece * pc,
    const acb_t alpha, const acb_t nu, const acb_t W, const acb_t c, int j, int end, double logtol)
{
    acb_t b;
    mag_t t;
    int ok;

    acb_init(b);
    mag_init(t);

    arb_set_d(acb_realref(b), pc->pts[j].re);
    arb_set_d(acb_imagref(b), pc->pts[j].im);
    _sd_bound_segment(err, 1, alpha, 1, nu, W, b, c);

    ok = _sd_bound_end(t, alpha, 1, nu, pc, c, j, end, P->lin_dir, P, logtol);
    mag_add(err, err, t);

    acb_clear(b);
    mag_clear(t);
    return ok && mag_is_finite(err);
}

/* The expansion for one piece: res = the integral over the piece (in
   the order of its vertices), without the sign.  The parameters of the
   plan (n, r, tau_0 and the vertices j) were chosen heuristically; the
   error bound is computed here from the enclosures of the saddle, A
   and c, and is valid whatever the plan is. */
static int
sde_piece(acb_t res, const sde_plan * pl, const sd_problem * P, const sd_piece * pc,
    const acb_t alpha, const acb_t nu, slong prec, slong wp)
{
    acb_ptr b;
    acb_t w, A, c, S, t, u;
    arb_t gam, tau, rho;
    mag_t err, e, rem;
    cd dir, phi, d1, d2;
    slong n = pl->n, k;
    int ok = 0, side, mid = pc->mid;

    b = _acb_vec_init(n + 2);
    acb_init(w); acb_init(A); acb_init(c);
    acb_init(S); acb_init(t); acb_init(u);
    arb_init(gam); arb_init(tau); arb_init(rho);
    mag_init(err); mag_init(e); mag_init(rem);

    /* (the sign of A and the branch of the saddle from the doubles) */
    _sd_phi(&phi, &d1, &d2, P, pc->pts[mid]);

    if (!sde_saddle(w, A, c, alpha, nu, pc->pts[mid], d2, wp + 30))
        goto cleanup;

    /* b_1 = sqrt(-2/A), along the outgoing path */
    acb_inv(b + 1, A, wp);
    acb_mul_si(b + 1, b + 1, -2, wp);
    acb_sqrt(b + 1, b + 1, wp);
    dir = cd_sub(pc->pts[mid + 1], pc->pts[mid]);
    {
        double br = arf_get_d(arb_midref(acb_realref(b + 1)), ARF_RND_NEAR);
        double bi = arf_get_d(arb_midref(acb_imagref(b + 1)), ARF_RND_NEAR);
        double dot = br * dir.re + bi * dir.im;

        /* (the polygon leaves the saddle along +-b_1 up to a small
           angle; with the other sign the result would still be an
           enclosure, with a large bound for the rest of the path) */
        if (!(fabs(dot) > 0.5 * hypot(br, bi) * cd_abs(dir)))
            goto cleanup;
        if (dot < 0.0)
            acb_neg(b + 1, b + 1);
    }

    /* the radius of analyticity */
    if (!sde_bound_radius(rho, A, nu, pl->r))
        goto cleanup;

    /* the recurrence only needs the relative accuracy of the result,
       not the bits for the cancellation in c */
    sde_coeffs(b, A, nu, n, prec, FLINT_MIN(wp, prec + 14 + 2 * FLINT_BIT_COUNT(n)));

    /* sum_{k odd} k b_k Gamma(k/2) */
    arb_const_sqrt_pi(gam, wp);
    acb_zero(S);
    for (k = 1; k <= n; k += 2)
    {
        acb_mul_ui(t, b + k, k, wp);
        acb_addmul_arb(S, t, gam, wp);
        arb_mul_ui(gam, gam, k, wp);
        arb_mul_2exp_si(gam, gam, -1);
    }

    if (!acb_is_finite(S))
        goto cleanup;

    /* the error, relative to exp(c), on each side */
    for (side = 0; side < 2; side++)
    {
        double t0 = pl->tau0[side];

        if (!sde_bound_truncation(e, rem, rho, pl->r, t0, n, pl->len))
            goto cleanup;
        mag_add(err, err, e);

        /* u(+-tau_0): the Taylor polynomial and its remainder */
        arb_set_d(tau, side ? t0 : -t0);
        acb_zero(u);
        for (k = n; k >= 1; k--)
        {
            acb_add(u, u, b + k, 64);
            acb_mul_arb(u, u, tau, 64);
        }
        acb_add_error_mag(u, rem);
        acb_add(u, u, w, 64);

        if (!sde_bound_rest(e, P, pc, alpha, nu, u, c, pl->j[side], side, pl->logtol))
            goto cleanup;
        mag_add(err, err, e);
    }

    acb_add_error_mag(S, err);

    acb_exp(t, c, wp);
    acb_mul(res, S, t, wp);
    ok = acb_is_finite(res);

cleanup:
    _acb_vec_clear(b, n + 2);
    acb_clear(w); acb_clear(A); acb_clear(c);
    acb_clear(S); acb_clear(t); acb_clear(u);
    arb_clear(gam); arb_clear(tau); arb_clear(rho);
    mag_clear(err); mag_clear(e); mag_clear(rem);
    return ok;
}

int
_sd_expand(acb_t res, const sd_contour * C, const sd_problem * P,
    const acb_t alpha, int sym, const acb_t nu, slong prec, double maxratio)
{
    sde_plan pl[SD_MAXPIECES];
    acb_t sum, gsum, t;
    double cost = 0.0, qcost = 0.0, mag = 0.0;
    slong wp;
    int i, pass, ok = 1, nsel = 0;
    char sel[SD_MAXPIECES];

    /* only implemented for the symmetric kernel */
    if (!sym || C->npieces == 0 || !acb_is_exact(alpha) || !acb_is_exact(nu)
        || P->has_start)
        return 0;

    /* the cancellation in c = phi(w_s), whose imaginary part is needed
       to an absolute accuracy of 2^-prec */
    for (i = 0; i < C->npieces; i++)
    {
        const sd_piece * pc = C->piece + i;
        if (pc->mid < 0)
            continue;
        mag = fmax(mag, P->aal * exp(fabs(pc->pts[pc->mid].re)) * 2.0
            + P->anu * cd_abs(pc->pts[pc->mid]) + cd_abs(pc->c) + 2.0);
    }
    if (!(mag < 1e300))
        return 0;
    wp = prec + 16 + (slong) log2(mag);

    /* first without checking which pieces are below the tolerance */
    for (pass = 0; pass < 2; pass++)
    {
        cost = qcost = 0.0;

        for (i = 0; i < C->npieces; i++)
        {
            if ((pass == 0 || pl[i].n < 0)
                && !sde_prepare(pl + i, P, C->piece + i, C->logL[0], prec, wp, pass,
                    maxratio * C->npieces * _sd_quad_cost(prec + 17.0 + log2(mag))))
                return 0;
            cost += pl[i].cost;
            qcost += (pl[i].n <= 0) ? 15.0 : _sd_quad_cost(pl[i].prec + 12.0 + log2(mag));
        }

        if (cost > maxratio * qcost)
            return 0;
    }

    acb_init(sum); acb_init(gsum); acb_init(t);

    for (i = 0; i < C->npieces && ok; i++)
    {
        const sd_piece * pc = C->piece + i;

        sel[i] = (pl[i].n == 0);
        nsel += sel[i];
        if (sel[i])
            continue;

        ok = sde_piece(t, pl + i, P, pc, alpha, nu, pl[i].prec,
                wp - (prec - pl[i].prec) + 2 * FLINT_BIT_COUNT(pl[i].n));
        if (pc->sign < 0)
            acb_neg(t, t);
        acb_add(pc->geom ? gsum : sum, pc->geom ? gsum : sum, t, wp);
    }

    /* the pieces with the factor 1 / (1 - q), q = exp(-2 pi i sigma nu) */
    if (ok && !acb_is_zero(gsum))
    {
        acb_const_pi(t, wp);
        acb_mul(t, t, nu, wp);
        acb_mul_onei(t, t);
        acb_mul_si(t, t, (P->lin_dir.im > 0.0) ? -2 : 2, wp);
        acb_exp(t, t, wp);
        /* (the geometric series must converge: |q| < 1) */
        {
            mag_t qm;
            mag_init(qm);
            acb_get_mag(qm, t);
            ok = mag_cmp_2exp_si(qm, 0) < 0;
            mag_clear(qm);
        }
        acb_sub_ui(t, t, 1, wp);
        acb_neg(t, t);
        acb_div(gsum, gsum, t, wp);
        acb_add(sum, sum, gsum, wp);
    }

    /* the pieces estimated to be below the tolerance: by quadrature
       (cheaply, with the tolerance of the whole contour) */
    if (ok && nsel > 0)
    {
        ok = _sd_integrate_pieces(t, C, P, alpha, sym, nu, NULL, prec, sel);
        acb_add(sum, sum, t, wp);
    }

    if (ok)
    {
        acb_set(res, sum);
        ok = acb_is_finite(res);
    }

    acb_clear(sum); acb_clear(gsum); acb_clear(t);
    return ok;
}
