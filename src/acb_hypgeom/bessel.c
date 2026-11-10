/*
    Copyright (C) 2026 Fredrik Johansson
    Developed using Claude Fable 5

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/* Automatic algorithm selection for J, Y, I, K (and the scaled
   e^{-z} I, e^{z} K), choosing between the z = 0 series, the asymptotic
   expansion at infinity, and rigorous numerical integration along
   steepest-descent contours (bessel_integration.c).

   The methods are tried in order of cheapness at a doubling working
   precision, each result verified a posteriori and the best kept:

     * first, unless the integration is expected to be cheaper, the
       series-and-asymptotics implementation
       (acb_hypgeom_bessel_*_nointegration) at the precision raised by
       the estimated cancellation of the series, or for K and Y the
       asymptotic expansion directly where it reaches the precision;
       this settles the easy cases (small orders and arguments in
       particular);

     * for an inexact z and exact nu, evaluation at the midpoint of z
       with the propagated error bounded through the differential
       equation (see midpoint_route_ode), at little more than the cost
       of an evaluation with exact input;

     * the series and the asymptotic expansion, in whichever order the
       cheaper one comes first (the series costs ~|z| terms, the
       expansion ~wp/2), starting at the precision raised by the
       estimated cancellation of the series, which absorbs the rest as
       the precision doubles, the expansion being self-limiting;

     * numerical integration, which evaluates the function on the input
       balls directly with essentially optimal accuracy, at a cost
       nearly independent of the parameters (some hundred microseconds
       at machine precision): it is tried first when its estimated cost
       is below that of the series and expansions, and otherwise as
       soon as it is cheaper than their next round;

     * for other inexact inputs, once the working precision exceeds the
       input accuracy, evaluation at the input midpoints with the
       propagated error added from rigorous derivative bounds,

           |F(p) - F(mid)| <= sup |dF/dz| rz + sup |dF/dnu| rnu

       (bounds from the same steepest-descent contours); since the
       integration gives about the same enclosure without evaluating
       anything at the midpoint, it is preferred here unless much more
       expensive, and its result is then final;

     * finally the series-and-asymptotics implementation at the target
       precision.

   An attempt that fails or returns a poor enclosure is cheap relative
   to the one that succeeds, and the doubling schedule bounds the total
   wasted work by a constant factor of the final attempt. */

#include <math.h>
#include "acb.h"
#include "arb_hypgeom.h"
#include "acb_hypgeom.h"
#include "acb_hypgeom/impl.h"
#include "acb_hypgeom/sd.h"

#define KIND_J ACB_HYPGEOM_KIND_BESSEL_J
#define KIND_Y ACB_HYPGEOM_KIND_BESSEL_Y
#define KIND_I ACB_HYPGEOM_KIND_BESSEL_I
#define KIND_K ACB_HYPGEOM_KIND_BESSEL_K

static void bessel(acb_t res, int kind, int scaled, const acb_t nu,
    const acb_t z, slong prec);
static void
method_series(acb_t res, int kind, int scaled, const acb_t nu, const acb_t z, slong wp)
{
    switch (kind)
    {
        case KIND_J: acb_hypgeom_bessel_j_0f1(res, nu, z, wp); break;
        case KIND_I: acb_hypgeom_bessel_i_0f1(res, nu, z, scaled, wp); break;
        case KIND_K: acb_hypgeom_bessel_k_0f1(res, nu, z, scaled, wp); break;
        default: acb_hypgeom_bessel_y_nointegration(res, nu, z, wp); break;
    }
}

static void
method_asymp(acb_t res, int kind, int scaled, const acb_t nu, const acb_t z, slong wp)
{
    acb_t h1, h2;

    switch (kind)
    {
        case KIND_J: acb_hypgeom_bessel_j_asymp(res, nu, z, wp); break;
        case KIND_I: acb_hypgeom_bessel_i_asymp(res, nu, z, scaled, wp); break;
        case KIND_K: acb_hypgeom_bessel_k_asymp(res, nu, z, scaled, wp); break;
        default:
            acb_init(h1);
            acb_init(h2);
            acb_hypgeom_bessel_h1_asymp(h1, nu, z, wp);
            acb_hypgeom_bessel_h2_asymp(h2, nu, z, wp);
            acb_sub(res, h1, h2, wp);
            acb_div_onei(res, res);
            acb_mul_2exp_si(res, res, -1);
            acb_clear(h1);
            acb_clear(h2);
    }
}

static void
method_nointegration(acb_t res, int kind, int scaled, const acb_t nu, const acb_t z, slong wp)
{
    switch (kind)
    {
        case KIND_J: acb_hypgeom_bessel_j_nointegration(res, nu, z, wp); break;
        case KIND_I: _acb_hypgeom_bessel_i_nointegration(res, nu, z, scaled, wp); break;
        case KIND_K: _acb_hypgeom_bessel_k_nointegration(res, nu, z, scaled, wp); break;
        default: acb_hypgeom_bessel_y_nointegration(res, nu, z, wp); break;
    }
}

/* an upper bound for the absolute value of e^{-z} (scaled I) or e^{z}
   (scaled K) over the ball z */
static void
scale_factor_bound(mag_t bound, int kind, const acb_t z)
{
    arb_t r;

    arb_init(r);
    if (kind == KIND_I)
        arb_neg(r, acb_realref(z));
    else
        arb_set(r, acb_realref(z));
    arb_exp(r, r, MAG_BITS);
    arb_get_mag(bound, r);
    arb_clear(r);
}

/* the radius of an enclosure, as a single magnitude */
static void
enclosure_rad(mag_t r, const acb_t x)
{
    mag_hypot(r, arb_radref(acb_realref(x)), arb_radref(acb_imagref(x)));
}

/* replace best by cand if cand is more informative */
static void
keep_best(acb_t best, slong * best_acc, const acb_t cand)
{
    slong a;

    if (!acb_is_finite(cand))
        return;

    a = acb_rel_accuracy_bits(cand);

    if (a > *best_acc)
    {
        acb_set(best, cand);
        *best_acc = a;
        return;
    }

    /* the relative accuracy is degenerate for an enclosure straddling
       zero; since every candidate encloses the same value, the tighter
       one is more informative: decide on the radius whenever neither
       candidate has a relative digit */
    if (a <= 0 && *best_acc <= 0)
    {
        mag_t rc, rb;
        mag_init(rc);
        mag_init(rb);
        enclosure_rad(rc, cand);
        enclosure_rad(rb, best);
        if (!acb_is_finite(best) || mag_cmp(rc, rb) < 0)
        {
            acb_set(best, cand);
            *best_acc = a;
        }
        mag_clear(rc);
        mag_clear(rb);
    }
}

/* the propagated error sup |dF/dz| rz + sup |dF/dnu| rnu (for the
   scaled functions, e^{-+z} F has z-derivative e^{-+z} (F' -+ F)) */
static void
propagated_error(mag_t err, const mag_t Dz, const mag_t Dnu, const mag_t M,
    int kind, int scaled, const acb_t nu, const acb_t z)
{
    mag_t rr, D;

    mag_init(rr);
    mag_init(D);
    mag_zero(err);

    if (!acb_is_exact(z))
    {
        mag_set(D, Dz);
        if (scaled)
            mag_add(D, D, M);
        enclosure_rad(rr, z);
        mag_mul(rr, rr, D);
        mag_add(err, err, rr);
    }

    if (!acb_is_exact(nu))
    {
        enclosure_rad(rr, nu);
        mag_mul(rr, rr, Dnu);
        mag_add(err, err, rr);
    }

    if (scaled)
    {
        scale_factor_bound(rr, kind, z);
        mag_mul(err, err, rr);
    }

    mag_clear(rr);
    mag_clear(D);
}

/* Evaluation at the input midpoints (at precision mprec) with the
   propagated error sup |dF/dz| rz + sup |dF/dnu| rnu added from the
   bounds along the steepest-descent contour; res is set to the
   enclosure, or to an indeterminate value on failure, and *mid_failed
   is set if the midpoint value has no relative accuracy (then neither
   has anything else). */
static void
midpoint_route(acb_t res, int * mid_failed, int kind, int scaled,
    const acb_t nu, const acb_t z, slong mprec)
{
    mag_t Dz, Dnu, M, err;
    acb_t numid, zmid;

    mag_init(Dz);
    mag_init(Dnu);
    mag_init(M);
    mag_init(err);
    acb_init(numid);
    acb_init(zmid);

    acb_get_mid(numid, nu);
    acb_get_mid(zmid, z);
    bessel(res, kind, scaled, numid, zmid, mprec);

    *mid_failed = !acb_is_finite(res) || acb_rel_accuracy_bits(res) <= 0;

    if (!*mid_failed && _acb_hypgeom_bessel_integration_bounds(M, Dz, Dnu, kind, nu, z))
    {
        propagated_error(err, Dz, Dnu, M, kind, scaled, nu, z);
        acb_add_error_mag(res, err);
    }
    else
        acb_indeterminate(res);

    mag_clear(Dz);
    mag_clear(Dnu);
    mag_clear(M);
    mag_clear(err);
    acb_clear(numid);
    acb_clear(zmid);
}

/* Evaluation at the midpoint of z for exact nu, with the propagated
   error r sup |g'| bounded through the differential equation.  The
   function g (J, Y, I, K, or e^{sigma z} times I or K with sigma = -1,
   +1 for the scaled functions) satisfies

       g'' = a g' - b g,    a = 2 sigma - 1/z,
                            b = 1 - nu^2/z^2             (J, Y)
                            b = -nu^2/z^2 - sigma/z      (I, K)

   so that v = (g, c g') for a constant c > 0 satisfies v' = A v with
   A = [[0, 1/c], [-c b, a]], and along the segment from the midpoint
   m to any point of the ball (of radius r)

       |v(z)| <= |v(m)| exp(r L),   L^2 = 1/c^2 + c^2 sup |b|^2 + sup |a|^2,

   whence |g'(z)| <= sqrt(|g(m)|^2 / c^2 + |g'(m)|^2) exp(r L).  The
   derivative at the midpoint is obtained from one more function value,
   g' = -+g_{nu+1} + (nu/z + sigma) g.  The constant c balances the two
   components (c ~ |b|^(-1/2)), and is increased where |g'| is small
   compared to |g| as far as the exponential factor allows.  This
   costs little more than the evaluation at the midpoint, and the
   bound is within a small factor of |g'(m)| unless the ball is so
   wide that nothing can be said anyway (or unless g' nearly vanishes
   at the midpoint).  Returns 0 if not applicable, and 2 if the bound
   is within a factor four of r |g'(m)|, so that the enclosure is about
   as good as possible; *mid_failed is set if the midpoint value has no
   relative accuracy. */
static int
midpoint_route_ode(acb_t res, int * mid_failed, int kind, int scaled,
    const acb_t nu, const acb_t z, double zn, slong mprec)
{
    acb_t zm, nu1, g1, a, b, t;
    mag_t r, G0, G1, Q, L, cu, cl, u, v;
    double c, c0, c1, cap, Qd, rd, G0d, G1d;
    slong lp;
    int sigma, ok = 0;

    /* (nu + 1 must be exact) */
    if (!acb_is_exact(nu) || acb_bits(nu) > 1000
        || (arf_cmpabs_2exp_si(arb_midref(acb_realref(nu)), -1000) < 0
            && !arb_is_zero(acb_realref(nu)))
        || (arf_cmpabs_2exp_si(arb_midref(acb_imagref(nu)), -1000) < 0
            && !arb_is_zero(acb_imagref(nu)))
        || arf_cmpabs_2exp_si(arb_midref(acb_realref(nu)), 1000) > 0
        || arf_cmpabs_2exp_si(arb_midref(acb_imagref(nu)), 1000) > 0)
        return 0;

    /* the segments must not cross the branch cut */
    if (!(arb_is_positive(acb_realref(z)) || !arb_contains_zero(acb_imagref(z))
        || ((kind == KIND_J || kind == KIND_I) && acb_is_int(nu))))
        return 0;

    acb_init(zm); acb_init(nu1); acb_init(g1);
    acb_init(a); acb_init(b); acb_init(t);
    mag_init(r); mag_init(G0); mag_init(G1); mag_init(Q); mag_init(L);
    mag_init(cu); mag_init(cl); mag_init(u); mag_init(v);

    sigma = scaled ? ((kind == KIND_I) ? -1 : 1) : 0;
    lp = 20 + (slong) log2(2.0 + zn);

    acb_get_mid(zm, z);
    bessel(res, kind, scaled, nu, zm, mprec);
    *mid_failed = !acb_is_finite(res) || acb_rel_accuracy_bits(res) <= 0;

    if (*mid_failed)
        goto cleanup;

    /* g'(m) */
    acb_add_ui(nu1, nu, 1, ARF_PREC_EXACT);
    bessel(g1, kind, scaled, nu1, zm, lp);
    if (kind != KIND_I)
        acb_neg(g1, g1);
    acb_div(t, nu, zm, lp);
    acb_add_si(t, t, sigma, lp);
    acb_addmul(g1, t, res, lp);

    if (!acb_is_finite(g1))
        goto cleanup;

    acb_get_mag(G0, res);
    acb_get_mag(G1, g1);

    /* a and b on the ball */
    acb_inv(t, z, lp);
    acb_set_si(a, 2 * sigma);
    acb_sub(a, a, t, lp);
    acb_mul(b, nu, t, lp);
    acb_mul(b, b, b, lp);
    acb_neg(b, b);
    if (kind == KIND_J || kind == KIND_Y)
        acb_add_ui(b, b, 1, lp);
    else if (sigma == 1)
        acb_sub(b, b, t, lp);
    else if (sigma == -1)
        acb_add(b, b, t, lp);

    acb_get_mag(Q, b);
    enclosure_rad(r, z);

    Qd = mag_get_d(Q);
    rd = mag_get_d(r);
    G0d = mag_get_d(G0);
    G1d = mag_get_d(G1);

    c0 = (Qd > 0.0) ? 1.0 / sqrt(Qd) : 1e100;
    c1 = (G1d > 0.0) ? 8.0 * G0d / G1d : 1e100;
    cap = (rd * Qd > 0.0) ? 0.125 / (rd * Qd) : 1e100;
    c = fmax(c0, fmin(c1, cap));

    if (!(c >= 1e-100 && c <= 1e100) || !mag_is_finite(Q))
        goto cleanup;

    mag_set_d(cu, c);
    mag_set_d_lower(cl, c);

    /* L */
    mag_mul(u, cu, Q);
    mag_mul(L, u, u);
    mag_one(u);
    mag_div(u, u, cl);
    mag_mul(v, u, u);
    mag_add(L, L, v);
    acb_get_mag(v, a);
    mag_mul(v, v, v);
    mag_add(L, L, v);
    mag_sqrt(L, L);

    /* sqrt(G0^2 / c^2 + G1^2) exp(r L) r */
    mag_mul(u, u, G0);
    mag_hypot(u, u, G1);
    mag_mul(L, L, r);
    mag_exp(L, L);
    mag_mul(u, u, L);
    mag_mul(u, u, r);

    if (mag_is_finite(u))
    {
        acb_add_error_mag(res, u);
        acb_get_mag_lower(v, g1);
        mag_mul_lower(v, v, r);
        mag_mul_2exp_si(v, v, 2);
        ok = (mag_cmp(u, v) <= 0) ? 2 : 1;
    }

cleanup:
    acb_clear(zm); acb_clear(nu1); acb_clear(g1);
    acb_clear(a); acb_clear(b); acb_clear(t);
    mag_clear(r); mag_clear(G0); mag_clear(G1); mag_clear(Q); mag_clear(L);
    mag_clear(cu); mag_clear(cl); mag_clear(u); mag_clear(v);

    return ok;
}

/* Rough cost estimates (in microseconds on a typical machine; only the
   ratios matter) of the series at the working precision wp (~1.4 |z| +
   wp/10 terms, plus the limit formula at integer order with ~|nu|
   terms), of the asymptotic expansions (~wp/2 terms, see below), and of the
   integration at precision prec: per descent path from a saddle, in
   N-fold double precision up to N = 4, and from there on growing
   quadratically. */
static double
series_cost(int kind, double num, double zn, int intorder, slong wp)
{
    double terms = 1.4 * zn + (double) wp / 10.0 + 10.0, c;

    /* Y and K through J_{-nu} and I_{-nu} need ~|nu| terms before the
       terms start decreasing when |z| is not small compared to |nu|
       (and the limit formula at integer order has a sum of |nu| terms) */
    if ((kind == KIND_Y || kind == KIND_K) && zn > 0.5 * num)
        terms += 1.5 * num;
    c = terms * (0.3 + 0.001 * (double) wp);
    if (intorder)
        c += 10.0 + 0.1 * (double) wp + 2.0 * (1.4 * zn + 10.0) * (0.3 + 0.001 * (double) wp);
    return c;
}

/* (for large orders, the terms of the expansion first grow, to about
   e^(|nu|^2 / (2 |z|)) at index |nu|^2 / (2 |z|)) */
static double
asymp_cost(double num, double zn, slong wp)
{
    double t = fmin(num * num / zn, 1e9);
    return 30.0 + (0.4 * (double) wp + 1.4 * t) * (0.3 + 0.001 * ((double) wp + 0.72 * t));
}

/* |phi''| = |sqrt(nu^2 -+ z^2)| at the saddles of the integrals for exact
   nu and z, where the integrals can be computed by the saddle-point
   expansions (sd_expand.c), and 0 otherwise */
static double
saddle_curvature(int kind, const acb_t nu, const acb_t z)
{
    double nr, ni, zr, zi, ar, ai;

    if (!acb_is_exact(nu) || !acb_is_exact(z))
        return 0.0;

    nr = arf_get_d(arb_midref(acb_realref(nu)), ARF_RND_NEAR);
    ni = arf_get_d(arb_midref(acb_imagref(nu)), ARF_RND_NEAR);
    zr = arf_get_d(arb_midref(acb_realref(z)), ARF_RND_NEAR);
    zi = arf_get_d(arb_midref(acb_imagref(z)), ARF_RND_NEAR);

    if (!(fabs(nr) + fabs(ni) + fabs(zr) + fabs(zi) < 1e150))
        return 0.0;

    ar = (zr - zi) * (zr + zi);
    ai = 2.0 * zr * zi;
    if (kind == KIND_J || kind == KIND_Y)
    {
        ar = -ar;
        ai = -ai;
    }
    ar += (nr - ni) * (nr + ni);
    ai += 2.0 * nr * ni;

    return sqrt(hypot(ar, ai));
}

/* absA: saddle_curvature() */
static double
integration_cost(int kind, double num, double zn, int real, double absA, slong prec)
{
    static const double dcost[4] = { 30.0, 45.0, 100.0, 200.0 };
    double bits = (double) prec + 12.0 + log2(2.0 + zn + num), c;
    int npaths;

    if (bits <= 4 * 53)
        c = dcost[(int) ((bits - 1.0) / 53.0)];
    else
        c = 600.0 * (bits / 256.0) * (bits / 256.0);

    /* per saddle (two paths): the planning of the contour and the
       expansion in place of the quadrature, where that matters */
    if (absA > 0.0 && bits > 4 * 53)
        c = fmin(c, 0.5 * (45.0 + 0.1 * bits + _sd_expand_cost(absA, num / absA, real, prec)));

    /* one saddle for K and for J, I in the monotone regime; two for
       the oscillatory functions (one with real variables, where the
       contributions are complex conjugates) */
    if (kind == KIND_K)
        npaths = 2;
    else if (kind == KIND_Y || zn > num)
        npaths = real ? 2 : 4;
    else
        npaths = 2;

    return npaths * c;
}

/* The cancellation (in bits) which the series must absorb through its
   working precision: the terms grow to about e^|z| while J and Y are of
   size e^|Im z|, I of size e^|Re z| and K of size e^(-Re z); for J and
   I, no more than the largest term (z/2)^(2k) / (k! (nu+1)_k), at
   k ~ nser, exceeds the first one (for large orders, the terms hardly
   grow). */
static double
series_cancellation(int kind, double num, double zn, double nser, const acb_t z)
{
    double zre = fabs(arf_get_d(arb_midref(acb_realref(z)), ARF_RND_NEAR));
    double zim = fabs(arf_get_d(arb_midref(acb_imagref(z)), ARF_RND_NEAR));
    double cancel, growth;

    if (kind == KIND_K)
        return 2.9 * zre;

    cancel = 1.45 * (zn - ((kind == KIND_I) ? zre : zim));

    if (kind != KIND_Y && nser >= 1.0 && num < 1e15)
    {
        growth = 2.0 * nser * log(0.5 * zn) - lgamma(nser + 1.0)
            - (lgamma(nser + num + 1.0) - lgamma(num + 1.0));
        cancel = fmin(cancel, 1.45 * growth + 10.0);
    }

    return fmax(cancel, 0.0);
}

/* whether the asymptotic expansion reaches wp bits: its error is about
   e^(-2 |z|) for |z| beyond the turning point region (empirically, from
   |z| ~ 1.9 |nu|) */
static int
asymp_converges(double num, double zn, slong wp)
{
    return zn > 1.9 * num + 2.0 && 2.885 * zn >= (double) wp + 3.0;
}

/* the expected cost of the series and expansions until success at the
   goal precision: the asymptotic expansion where it converges, else the
   series at the precision raised by its cancellation; infinite where
   neither applies */
static double
first_cost(int kind, double num, double zn, double cancel,
    int intorder, int series_plausible, slong goal)
{
    if (asymp_converges(num, zn, goal))
        return asymp_cost(num, zn, goal + 30);
    if (series_plausible)
        return 1.2 * series_cost(kind, num, zn, intorder, goal + 30 + (slong) cancel);
    return INFINITY;
}

/* whether the integration can be attempted at all (the contours are
   planned in double precision) */
static int
integration_applicable(const acb_t nu, const acb_t z)
{
    return arf_cmpabs_2exp_si(arb_midref(acb_realref(z)), 44) < 0
        && arf_cmpabs_2exp_si(arb_midref(acb_imagref(z)), 44) < 0
        && arf_cmpabs_2exp_si(arb_midref(acb_realref(nu)), 44) < 0
        && arf_cmpabs_2exp_si(arb_midref(acb_imagref(nu)), 44) < 0
        && !(arf_cmpabs_2exp_si(arb_midref(acb_realref(z)), -900) < 0
             && arf_cmpabs_2exp_si(arb_midref(acb_imagref(z)), -900) < 0);
}

/* Heuristic estimate of the relative accuracy attainable from inexact
   input balls: the relative condition numbers are roughly
   |z| + |nu| + 1 with respect to z and |nu| (1 + |log(2(|nu|+1)/|z|)|)
   with respect to nu (for small |nu|, the accuracy of nu is counted
   relative to 1: the functions are not sensitive to the relative
   error of an order close to zero).  Near a negative integer -n at distance d (except
   for K, which is even in nu), J_nu ~ (-1)^n J_n + O(d) Y_n and similarly
   for Y and I, so the condition number with respect to nu can grow by a
   factor up to min(1/d, |Y_n / J_n|), where |Y_n / J_n| is roughly
   (2n / (e|z|))^(2n) for n > |z|.  This is only used to decide when to
   stop trying harder; it does not affect the correctness of the
   enclosure. */
static slong
attainable_accuracy(int kind, const acb_t nu, const acb_t z, double num,
    double zn, slong prec)
{
    slong acc = prec, t;

    if (!acb_is_exact(z))
    {
        t = acb_rel_accuracy_bits(z);
        if (t > 0)
            t -= (slong) FLINT_MIN(log2(1.0 + zn + num), 1e6);
        acc = FLINT_MIN(acc, t);
    }

    if (!acb_is_exact(nu))
    {
        double lc;

        lc = log2(1.0 + num * (1.0 + fabs(log(2.0 * (num + 1.0))
                                          - log(fmax(zn, 1e-300)))));

        if (kind != KIND_K && arf_cmp_si(arb_midref(acb_realref(nu)), 0) < 0)
        {
            double n, d, lpole;
            arf_t u, v;

            arf_init(u);
            arf_init(v);

            /* distance of the midpoint to the nearest integer, exactly */
            arf_set_round(u, arb_midref(acb_realref(nu)), 53, ARF_RND_NEAR);
            n = -arf_get_d(u, ARF_RND_NEAR);
            arf_floor(u, arb_midref(acb_realref(nu)));
            arf_sub(v, arb_midref(acb_realref(nu)), u, 53, ARF_RND_NEAR);
            d = arf_get_d(v, ARF_RND_NEAR);
            d = fmin(d, 1.0 - d);
            d = fmax(d, fabs(arf_get_d(arb_midref(acb_imagref(nu)), ARF_RND_NEAR)));

            arf_clear(u);
            arf_clear(v);

            n = floor(n + 0.5);
            lpole = 2.0 * n * log2(fmax(2.0 * n / (2.718281828 * fmax(zn, 1e-300)), 1.0));
            if (d > 0.0)
                lpole = fmin(lpole, -log2(d));
            lpole = fmax(lpole, 0.0);
            lc = fmax(lc, log2(num) + lpole) + 1.0;
        }

        t = acb_rel_one_accuracy_bits(nu);
        if (t > 0)
            t -= (slong) FLINT_MIN(lc, 1e6);
        acc = FLINT_MIN(acc, t);
    }

    return acc;
}

/* the function is real for real nu on the nonnegative real z-axis, and
   for integer nu on the whole real axis (J and I); the computed
   imaginary part is then a ball around the exact value 0 */
static int
bessel_real_case(int kind, const acb_t nu, const acb_t z)
{
    if (!acb_is_real(nu) || !acb_is_real(z))
        return 0;
    return arb_is_nonnegative(acb_realref(z))
        || ((kind == KIND_J || kind == KIND_I) && acb_is_int(nu));
}

/* kind is J, Y, I or K; scaled is only allowed with I and K */
static void
bessel(acb_t res, int kind, int scaled, const acb_t nu, const acb_t z, slong prec)
{
    double zn, num, nser, nlim, cancel, c_int, c_first;
    slong wp, acc_in, goal, target, accept, best_acc = -WORD_MAX;
    acb_t best, t;
    int inexact, real_case, round, series_plausible, intorder;
    int mid_tried = 0, mid_failed = 0, int_ok, int_tried = 0, int_done = 0, int_first;
    slong mprec;

    if (res == nu || res == z)
    {
        acb_init(t);
        bessel(t, kind, scaled, nu, z, prec);
        acb_swap(res, t);
        acb_clear(t);
        return;
    }

    /* reflection formulas: K_{-nu} = K_nu, and for integer n,
       J_{-n} = (-1)^n J_n, Y_{-n} = (-1)^n Y_n, I_{-n} = I_n; the
       methods below work better with Re(nu) >= 0 */
    if (acb_is_finite(nu) && arf_sgn(arb_midref(acb_realref(nu))) < 0
        && (kind == KIND_K || (acb_is_int(nu) && acb_is_exact(nu))))
    {
        acb_init(t);
        acb_neg(t, nu);
        bessel(res, kind, scaled, t, z, prec);
        if ((kind == KIND_J || kind == KIND_Y)
            && !arf_is_int_2exp_si(arb_midref(acb_realref(nu)), 1))
            acb_neg(res, res);
        acb_clear(t);
        return;
    }

    real_case = bessel_real_case(kind, nu, z);

    if (!acb_is_finite(nu) || !acb_is_finite(z) || acb_contains_zero(z))
    {
        method_nointegration(res, kind, scaled, nu, z, prec);
        goto finish;
    }

    zn = hypot(arf_get_d(arb_midref(acb_realref(z)), ARF_RND_NEAR),
               arf_get_d(arb_midref(acb_imagref(z)), ARF_RND_NEAR));
    num = hypot(arf_get_d(arb_midref(acb_realref(nu)), ARF_RND_NEAR),
                arf_get_d(arb_midref(acb_imagref(nu)), ARF_RND_NEAR));

    /* inputs not faithfully representable as doubles: the cost
       estimate below is meaningless there, and the series and
       asymptotic expansions handle all such inputs quickly */
    if (!isfinite(zn) || !isfinite(num)
        || zn > 1e300 || num > 1e300 || zn < 1e-300)
    {
        method_nointegration(res, kind, scaled, nu, z, prec);
        goto finish;
    }

    inexact = !acb_is_exact(nu) || !acb_is_exact(z);

    /* the goal is capped by the relative accuracy of the inputs (that
       of nu relative to max(1, |nu|)):
       nothing more is attainable, and hopeless input balls then
       terminate after the first cheap attempt; the target also
       accounts (heuristically) for the conditioning */
    acc_in = FLINT_MIN(acb_rel_one_accuracy_bits(nu),
                       acb_rel_accuracy_bits(z));
    acc_in = FLINT_MAX(acc_in, 0);
    goal = FLINT_MIN(prec, acc_in);
    target = FLINT_MIN(goal, attainable_accuracy(kind, nu, z, num, zn, prec));

    /* Accept results this close to the target.  For inexact input the
       target is only a heuristic estimate of what is attainable; a
       result short of it is compared with the integration on the
       input balls. */
    accept = target - 5;

    /* nothing is attainable: direct interval propagation is as good as
       anything, quickly */
    if (goal <= 0 || target < 4)
    {
        method_nointegration(res, kind, scaled, nu, z, prec);
        goto finish;
    }

    /* number of terms of the 0F1 series (the connection formulas for Y
       and K cost a small multiple); at integer order, Y and K are
       computed from a limit formula which also includes a finite sum of
       |nu| terms (this is feasible for moderately large |nu|) */
    nser = 0.5 * (-(num + 1.0) + sqrt((num + 1.0) * (num + 1.0) + zn * zn));
    intorder = (kind == KIND_Y || kind == KIND_K) && acb_contains_int(nu);
    nlim = intorder ? num : 0.0;

    series_plausible = (nser <= 10.0 * (16.0 * (double) goal + 480.0))
                        && (nlim <= 20.0 * (double) goal + 20000.0);

    /* precision for the midpoint evaluation */
    mprec = (acc_in < prec - 30) ? acc_in + 30 : prec;

    /* the integration goes first when it is expected to be cheaper
       than the series (at the precision raised by the cancellation it
       must absorb) or the asymptotic expansion, where that converges */
    int_ok = integration_applicable(nu, z);
    c_int = integration_cost(kind, num, zn, real_case,
                saddle_curvature(kind, nu, z), FLINT_MIN(prec, mprec));
    cancel = series_cancellation(kind, num, zn, nser, z);
    c_first = first_cost(kind, num, zn, cancel, intorder, series_plausible, goal);
    int_first = int_ok && (c_int < c_first);

    acb_init(best);
    acb_init(t);
    acb_indeterminate(best);

    /* The initial working precision: with the cancellation of the
       series, when that is the method to be used.  With inexact input,
       the goal may be below prec; the series can still lose up to
       ~1.45 |z| bits to cancellation (which is only recovered by
       increasing the working precision for exact input). */
    wp = goal;
    if (goal < prec)
        wp = FLINT_MIN(prec, goal + (slong) FLINT_MIN((double) goal, 1.45 * zn));
    else if (series_plausible && !asymp_converges(num, zn, goal))
        wp = goal + (slong) FLINT_MIN(cancel, 4.0 * (double) goal);

    for (wp = wp + 30, round = 0; wp <= 16 * goal + 480;
            wp = 2 * wp + 30, round++)
    {
        int m, series_first, series_done = 0;
        slong iter_start_acc = best_acc;

        /* the series needs ~nser terms, the expansion ~wp/2; where the
           expansion is known to reach the working precision (its error
           is ~ e^{-2 |z|}), it comes first since
           it is self-limiting (and the limit formula for Y and K at
           integer order is expensive) */
        series_first = (FLINT_MAX(nser, nlim) <= 0.5 * (double) wp)
            && !asymp_converges(num, zn, wp);

        /* First, unless the integration is expected to be cheaper, the
           series-and-asymptotics implementation, which is tuned for
           the easy cases, or the asymptotic expansion.  Only the
           implementation of J raises the precision by the cancellation
           of the series. */
        if (round == 0 && !int_first)
        {
            slong p0 = prec + 10;
            int conv, direct;

            /* K and Y by the expansion as soon as it converges (their
               series cancel heavily or are expensive); the
               implementations of J and I switch from the series (which
               is cheaper up to about |z| = prec / 2) by themselves,
               but for inexact input the expansion is preferable as
               there is no cancellation to inflate the input radii */
            conv = asymp_converges(num, zn, goal);
            direct = conv && (kind == KIND_K || inexact)
                && !(kind == KIND_Y && !arb_is_positive(acb_realref(z)));

            if (kind != KIND_J && !direct
                && !(conv && 2.0 * zn >= (double) p0))
                p0 += (slong) FLINT_MIN(cancel, 4.0 * (double) goal);

            /* (with an inexact z, the cancellation of the series
               inflates the radius regardless of the precision, and so
               does the expansion close to where it stops converging;
               the evaluation at the midpoint is then better) */
            if (!(inexact && acb_is_exact(nu) && cancel >= 4.0
                    && !asymp_converges(num, zn, goal + 12)))
            {
                if (direct)
                    method_asymp(t, kind, scaled, nu, z, inexact ? wp : p0);
                else
                    method_nointegration(t, kind, scaled, nu, z, inexact ? wp : p0);
                keep_best(best, &best_acc, t);
                /* (if only z is inexact, the evaluation at the midpoint
                   is cheap and nearly optimal: less slack) */
                if (best_acc >= ((inexact && acb_is_exact(nu)) ? target - 2 : accept))
                    goto done;
            }

            /* inexact z only: the evaluation at the midpoint with the
               derivative bounded through the differential equation */
            if (inexact)
            {
                int ode = midpoint_route_ode(t, &mid_failed, kind, scaled, nu, z, zn, mprec);

                if (ode)
                {
                    keep_best(best, &best_acc, t);
                    if (best_acc >= accept || ode == 2)
                        goto done;
                }
            }
        }

        /* the integration, first when expected to be cheapest overall,
           or once it is cheaper than the next round of the series and
           expansions */
        if (int_ok && !int_tried
            && (int_first || c_int <= (series_first ? series_cost(kind, num, zn, intorder, wp) : asymp_cost(num, zn, wp))))
        {
            int_tried = 1;
            if (_acb_hypgeom_bessel_integration(t, NULL, kind, nu, z, scaled, prec))
            {
                int_done = 1;
                /* (this is as good as it gets for inexact input) */
                keep_best(best, &best_acc, t);
                if (best_acc >= accept || (inexact && acb_rel_accuracy_bits(t) > 0))
                    goto done;
            }
        }

        /* inexact inputs: once the working precision exceeds the
           input accuracy, direct interval propagation cannot gain
           more; evaluate at the midpoint with certified derivative
           bounds for the input radii instead */
        if (inexact && !mid_tried && acc_in <= wp - 30)
        {
            mid_tried = 1;

            /* the integration evaluates on the input balls directly,
               and is preferred unless much more expensive */
            if (int_ok && !int_tried && c_int <= 4.0 * (c_first + 50.0))
            {
                int_tried = 1;
                if (_acb_hypgeom_bessel_integration(t, NULL, kind, nu, z, scaled, prec))
                {
                    int_done = 1;
                    /* (this is as good as it gets for inexact input) */
                    keep_best(best, &best_acc, t);
                    if (best_acc >= accept || acb_rel_accuracy_bits(t) > 0)
                        goto done;
                }
            }

            midpoint_route(t, &mid_failed, kind, scaled, nu, z, mprec);
            keep_best(best, &best_acc, t);

            /* the midpoint result is as good as this route gets; if it
               falls far short (the derivative bounds can be pessimistic
               for wide balls), integration with the input balls may do
               better, at a relatively high cost */
            if (best_acc >= FLINT_MIN(accept, target / 2))
                goto done;
            if (best_acc >= 2)
                break;
        }

        for (m = 0; m < 2; m++)
        {
            if ((m == 0) == series_first)
            {
                if (nser <= 4.0 * (double) prec + 1000.0 && series_plausible
                    && !(inexact && round == 0 && kind == KIND_Y))
                {
                    method_series(t, kind, scaled, nu, z, wp);
                    keep_best(best, &best_acc, t);
                    series_done = 1;
                }
            }
            else if (kind == KIND_Y && !arb_is_positive(acb_realref(z)))
            {
                /* the Hankel expansions need Re(z) > 0; the series
                   method for Y also uses the asymptotic expansion of J
                   where appropriate */
                if (!series_done)
                {
                    acb_hypgeom_bessel_y_nointegration(t, nu, z, wp);
                    keep_best(best, &best_acc, t);
                }
            }
            else
            {
                method_asymp(t, kind, scaled, nu, z, wp);
                keep_best(best, &best_acc, t);
            }
            if (best_acc >= accept)
                goto done;
        }

        /* diminishing returns: cancellation is recovered at about one
           bit per extra bit of working precision; stop doubling when
           the accuracy is not recovering.  While nothing useful has
           been found, keep going (the method choices inside the series
           code can change at a higher precision, and the series may
           need more terms than a lower precision allows), unless the
           series is hopelessly long. */
        if (round >= 1 && best_acc < iter_start_acc + wp / 8
            && (best_acc > 0 || !series_plausible))
            break;
    }

    /* inexact inputs which are more accurate than the precision can
       still spoil the series (near a singularity in nu of the formulas
       used); the midpoint evaluation is cheap */
    if (best_acc < accept && inexact && !mid_tried)
    {
        midpoint_route(t, &mid_failed, kind, scaled, nu, z, mprec);
        keep_best(best, &best_acc, t);
    }

    if (best_acc < accept && int_ok && !int_tried && !mid_failed
        && _acb_hypgeom_bessel_integration(t, NULL, kind, nu, z, scaled, prec))
    {
        int_done = 1;
        keep_best(best, &best_acc, t);
    }

    /* the real-variable integral representations for I and K (for
       parameters too large for the contours to be planned) */
    if (best_acc < accept && !int_done
        && (kind == KIND_I || kind == KIND_K)
        && acb_is_real(nu) && acb_is_real(z)
        && arf_cmp_d(arb_midref(acb_realref(nu)), 0.51) > 0
        && arf_cmpabs_2exp_si(arb_midref(acb_realref(nu)), 60) < 0
        && arf_cmp_2exp_si(arb_midref(acb_realref(z)), -16) > 0
        && arf_cmpabs_2exp_si(arb_midref(acb_realref(z)), 60) < 0)
    {
        acb_zero(t);
        if (kind == KIND_I)
            arb_hypgeom_bessel_i_integration(acb_realref(t), acb_realref(nu),
                acb_realref(z), scaled, prec + 30);
        else
            arb_hypgeom_bessel_k_integration(acb_realref(t), acb_realref(nu),
                acb_realref(z), scaled, prec + 30);
        keep_best(best, &best_acc, t);
    }

    if (best_acc < accept)
    {
        method_nointegration(t, kind, scaled, nu, z, prec);
        keep_best(best, &best_acc, t);
    }

done:
    acb_set_round(res, best, prec);
    acb_clear(best);
    acb_clear(t);

finish:
    if (real_case && acb_is_finite(res))
        arb_zero(acb_imagref(res));
}

void
acb_hypgeom_bessel_j(acb_t res, const acb_t nu, const acb_t z, slong prec)
{
    bessel(res, KIND_J, 0, nu, z, prec);
}

void
acb_hypgeom_bessel_y(acb_t res, const acb_t nu, const acb_t z, slong prec)
{
    bessel(res, KIND_Y, 0, nu, z, prec);
}

void
acb_hypgeom_bessel_i(acb_t res, const acb_t nu, const acb_t z, slong prec)
{
    bessel(res, KIND_I, 0, nu, z, prec);
}

void
acb_hypgeom_bessel_i_scaled(acb_t res, const acb_t nu, const acb_t z, slong prec)
{
    bessel(res, KIND_I, 1, nu, z, prec);
}

void
acb_hypgeom_bessel_k(acb_t res, const acb_t nu, const acb_t z, slong prec)
{
    bessel(res, KIND_K, 0, nu, z, prec);
}

void
acb_hypgeom_bessel_k_scaled(acb_t res, const acb_t nu, const acb_t z, slong prec)
{
    bessel(res, KIND_K, 1, nu, z, prec);
}

void
acb_hypgeom_bessel_jy(acb_t res1, acb_t res2, const acb_t nu, const acb_t z,
    slong prec)
{
    double zn, num;
    slong acc_in, goal, target, accept;
    acb_t t, u;
    int classic, hopeless, done = 0;

    zn = hypot(arf_get_d(arb_midref(acb_realref(z)), ARF_RND_NEAR),
               arf_get_d(arb_midref(acb_imagref(z)), ARF_RND_NEAR));
    num = hypot(arf_get_d(arb_midref(acb_realref(nu)), ARF_RND_NEAR),
                arf_get_d(arb_midref(acb_imagref(nu)), ARF_RND_NEAR));

    acb_init(t);
    acb_init(u);

    acc_in = FLINT_MIN(acb_rel_one_accuracy_bits(nu), acb_rel_accuracy_bits(z));
    acc_in = FLINT_MAX(acc_in, 0);
    goal = FLINT_MIN(prec, acc_in);
    target = FLINT_MIN(goal, attainable_accuracy(KIND_Y, nu, z, num, zn, prec));
    target = FLINT_MIN(target, attainable_accuracy(KIND_J, nu, z, num, zn, prec));
    accept = target - 5;

    hopeless = goal <= 0 || target < 4 || !acb_is_finite(nu) || !acb_is_finite(z)
        || !isfinite(zn) || !isfinite(num);

    /* where the series is short or the asymptotic expansion converges,
       computing J and Y together through the series and asymptotic
       expansions shares work; the result is kept if it is accurate */
    classic = hopeless || (zn < 24.0 && num < 64.0)
        || _acb_hypgeom_bessel_asymp_converges(nu, z, prec + 10);

    if (classic)
    {
        acb_hypgeom_bessel_jy_nointegration((res1 != NULL) ? t : NULL, u, nu, z, prec);

        if ((res1 == NULL || bessel_real_case(KIND_J, nu, z)) && acb_is_finite(t))
            arb_zero(acb_imagref(t));
        if (bessel_real_case(KIND_Y, nu, z) && acb_is_finite(u))
            arb_zero(acb_imagref(u));

        done = hopeless || (acb_rel_accuracy_bits(u) >= accept
            && (res1 == NULL || acb_rel_accuracy_bits(t) >= accept));
    }

    if (!done && res1 == NULL)
    {
        bessel(u, KIND_Y, 0, nu, z, prec);
        done = 1;
    }

    /* J and Y together from the two Hankel integrals along the
       steepest-descent contours, when that is expected to be
       cheaper than the series for Y alone */
    if (!done && integration_applicable(nu, z))
    {
        double nser = 0.5 * (-(num + 1.0) + sqrt((num + 1.0) * (num + 1.0) + zn * zn));
        int intorder = acb_contains_int(nu);
        int series_plausible = (nser <= 10.0 * (16.0 * (double) goal + 480.0))
                                && (!intorder || num <= 20.0 * (double) goal + 20000.0);

        if (integration_cost(KIND_Y, num, zn, bessel_real_case(KIND_Y, nu, z),
                    saddle_curvature(KIND_Y, nu, z), goal)
            < first_cost(KIND_Y, num, zn, series_cancellation(KIND_Y, num, zn, nser, z),
                    intorder, series_plausible, goal)
            && _acb_hypgeom_bessel_integration(t, u, ACB_HYPGEOM_KIND_BESSEL_JY, nu, z, 0, prec)
            && acb_rel_accuracy_bits(t) >= accept
            && acb_rel_accuracy_bits(u) >= accept)
        {
            if (bessel_real_case(KIND_J, nu, z))
                arb_zero(acb_imagref(t));
            if (bessel_real_case(KIND_Y, nu, z))
                arb_zero(acb_imagref(u));
            done = 1;
        }
    }

    if (!done)
    {
        bessel(t, KIND_J, 0, nu, z, prec);
        bessel(u, KIND_Y, 0, nu, z, prec);
    }

    if (res1 != NULL)
        acb_swap(res1, t);
    acb_swap(res2, u);
    acb_clear(t);
    acb_clear(u);
}
