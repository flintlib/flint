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

   For each method, the number of guard bits it needs (its
   cancellation) and its cost are estimated in double precision (see
   bessel_model); the cheapest is tried first, with its guard bits.
   Every result is a rigorous enclosure and is verified a posteriori,
   so that the estimates only affect the running time; the best result
   is kept.

     * The first attempt, unless the integration is expected to be
       cheaper: the series at the precision raised by its cancellation
       (for J and Y through the implementations
       acb_hypgeom_bessel_*_nointegration, which do that by themselves),
       or the asymptotic expansion with guard bits for its own
       cancellation where it reaches the precision and is cheaper.

     * With exact input, an attempt which falls short by a known number
       of bits is repeated with that many more bits.  Otherwise the
       working precision doubles from one round to the next, and in
       each round the series and the expansion are tried in the order
       of their cost (the expansion only where it is expected to get
       there).

     * Numerical integration evaluates the function on the input balls
       directly with essentially optimal accuracy, at a cost nearly
       independent of the parameters (some hundred microseconds at
       machine precision): it is tried first when its estimated cost is
       below that of the series and expansions, and otherwise before
       any attempt which is expected to cost more.

     * For an inexact z and exact nu, evaluation at the midpoint of z
       with the propagated error bounded through the differential
       equation (see midpoint_route_ode), at little more than the cost
       of an evaluation with exact input.

     * For other inexact inputs, once the working precision exceeds the
       input accuracy, evaluation at the input midpoints with the
       propagated error added from rigorous derivative bounds,

           |F(p) - F(mid)| <= sup |dF/dz| rz + sup |dF/dnu| rnu

       (bounds from the same steepest-descent contours); the
       integration gives about the same enclosure without evaluating
       anything at the midpoint, and is preferred unless more
       expensive; its result is then final.

     * Finally the series-and-asymptotics implementation at the target
       precision.

   An attempt that fails or returns a poor enclosure is cheap relative
   to the one that succeeds, and the schedule bounds the total wasted
   work by a constant factor of the final attempt. */

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

/* (in the main loop of bessel(): an attempt at the working precision wp
   which has fallen short of the accuracy to accept by a known number
   of bits) */
/* (the integration, if it has not been tried and is not expected to
   cost more than the attempt which would come next) */
#define TRY_INTEGRATION_IF(next_cost) \
    do { \
        if (int_ok && !int_tried && !inexact && c_int <= (next_cost)) \
        { \
            int_tried = 1; \
            if (_acb_hypgeom_bessel_integration(t, NULL, kind, nu, z, scaled, prec)) \
            { \
                int_done = 1; \
                keep_best(best, &best_acc, t); \
                if (best_acc >= accept) \
                    goto done; \
            } \
        } \
    } while (0)

#define NOTE_ATTEMPT() \
    do { \
        if (!inexact && acb_is_finite(t)) \
        { \
            slong a__ = acb_rel_accuracy_bits(t); \
            if (a__ > 0 && a__ < accept) \
                need_wp = FLINT_MIN(need_wp, wp + (accept - a__) + 20); \
        } \
    } while (0)
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

/* The estimates on which the choice of method is based.  All of this
   is heuristic, in double precision: it decides which method is tried
   first and with how many guard bits, and nothing else; every result
   is verified a posteriori. */
typedef struct
{
    double num, zn;     /* |nu|, |z| */
    double nr, ni, zr, zi;
    double lz, az;      /* log |z/2|, arg z */
    double nser;        /* index of the largest term of the series for J and I */
    double xi;          /* Re(xi) of the Debye approximations (see below; 0 for small orders) */
    double cancel;      /* bits lost by the series to cancellation */
    double cancel_pos, cancel_neg;  /* for Y at noninteger order (see series_cancellation) */
    double areach;      /* bits which the asymptotic expansion can reach */
    double aloss;       /* bits lost by the asymptotic expansion */
    double apen;        /* bits by which its error bound exceeds the first omitted term */
    double apeak;       /* size of its largest term, in bits */
    int real;
    int intorder;
    int hp;             /* which cost model (see series_cost) */
}
bessel_model;

/* log Gamma(x) for x > 0, roughly (these estimates are evaluated on
   every call of the functions, so that their cost matters at low
   precision) */
static double
d_lgamma(double x)
{
    if (x < 8.0)
        return lgamma(x);
    return (x - 0.5) * log(x) - x + 0.9189385332046727 + 1.0 / (12.0 * x);
}

/* Re log Gamma(x + iy), roughly: Stirling's formula after shifting the
   argument to the right.  For x < 1/2, the reflection formula without
   the factor pi / sin(pi (x + iy)) unless regularize is zero: with
   regularize set, this is the size of Gamma away from its poles and
   zeros, which is what matters for the terms of a series at (or near)
   an integer order. */
static double
d_lgamma_re(double x, double y, int regularize)
{
    double s = 0.0, r2;

    if (x < 0.5)
    {
        s = -d_lgamma_re(1.0 - x, -y, 0);
        if (!regularize || fabs(y) > 0.5)
        {
            double sn = sin(3.14159265358979323846 * x);
            double ay = 3.14159265358979323846 * fabs(y);
            /* log |sin(pi (x + iy))| */
            s += 1.1447298858494002
                - ((ay > 20.0) ? ay - 0.6931471805599453
                               : 0.5 * log(fmax(sn * sn + sinh(ay) * sinh(ay), 1e-300)));
        }
        return s;
    }

    if (y == 0.0)
        return d_lgamma(x);

    while (x < 8.0)
    {
        s -= 0.5 * log(x * x + y * y);
        x += 1.0;
    }

    r2 = x * x + y * y;
    return s + (x - 0.5) * 0.5 * log(r2) - y * atan2(y, x) - x
        + 0.9189385332046727 + x / (12.0 * r2);
}

/* Re(xi), xi = s - nu log((nu + s) / w), s = sqrt(nu^2 - w^2) on the
   branch continuous with s = nu at w = 0: the exponent of the
   leading-order (Debye) approximations, J_nu(w) ~ e^xi for |nu| > |w|
   and of size e^|Re xi| in the oscillatory region, Y_nu(w) of size
   e^|Re xi|, and I_nu(z) ~ e^xi, K_nu(z) ~ e^-xi with w = iz. */
double
_acb_hypgeom_bessel_debye_re(double nr, double ni, double wr, double wi)
{
    double qr, qi, m, sr, si, ur, ui, d;

    qr = (nr - ni) * (nr + ni) - (wr - wi) * (wr + wi);
    qi = 2.0 * (nr * ni - wr * wi);
    m = hypot(qr, qi);
    sr = sqrt(0.5 * fmax(m + qr, 0.0));
    si = sqrt(0.5 * fmax(m - qr, 0.0));
    if (qi < 0.0)
        si = -si;
    if (sr * nr + si * ni < 0.0)
    {
        sr = -sr;
        si = -si;
    }

    ur = nr + sr;
    ui = ni + si;
    m = hypot(ur, ui);
    if (!(m > 0.0))
        return sr;

    d = atan2(ui, ur) - atan2(wi, wr);
    if (d > 3.14159265358979323846)
        d -= 6.283185307179586;
    if (d < -3.14159265358979323846)
        d += 6.283185307179586;

    return sr - (nr * (log(m) - log(hypot(wr, wi))) - ni * d);
}

/* log of the absolute value of the k-th term (z/2)^(2k +- nu) /
   (k! Gamma(k +- nu + 1)) of the series for J and I (sign = 1) or of
   the series with -nu which Y and K also need (sign = -1) */
static double
series_log_term(const bessel_model * B, double k, int sign)
{
    return 2.0 * k * B->lz + sign * (B->nr * B->lz - B->ni * B->az) - d_lgamma(k + 1.0)
        - d_lgamma_re(k + 1.0 + sign * B->nr, sign * B->ni, 1);
}

/* The cancellation (in bits) which the series must absorb through its
   working precision: the largest term relative to the function, whose
   size is estimated by the Debye approximations (so that this also
   accounts for J being exponentially small beyond the turning point,
   where the alternating series loses twice the growth of its terms,
   and for the series of K and Y with large orders, which do not cancel
   as long as their terms decrease). */
static double
series_cancellation(int kind, bessel_model * B)
{
    double zre = fabs(B->zr), zim = fabs(B->zi), num = B->num, zn = B->zn;
    double cancel, lmax, lfun, x;

    /* the classical estimates: the terms grow to about e^|z|, while J
       and Y are of size e^|Im z|, I of size e^|Re z| and K of size
       e^(-Re z) */
    if (kind == KIND_K)
        cancel = 1.45 * (zn + B->zr);
    else
        cancel = 1.45 * (zn - ((kind == KIND_I) ? zre : zim));

    B->cancel_pos = B->cancel_neg = cancel;

    /* (the terms of the series for J and I decrease from the first one) */
    if ((kind == KIND_J || kind == KIND_I) && B->nser < 1.0 && B->nr >= 0.0)
        return 0.0;

    /* (these are accurate for small orders; the estimates below are
       for Re(nu) >= 0, and have only been validated for orders close
       to the real axis) */
    if (num <= 2.0 || B->nr < 0.0 || fabs(B->ni) > 0.05 * num
        || !(num < 1e100 && zn < 1e100 && zn > 1e-100))
    {
        if ((kind == KIND_J || kind == KIND_I) && B->nser >= 1.0 && num < 1e15)
        {
            double growth = 2.0 * B->nser * B->lz - d_lgamma(B->nser + 1.0)
                - (d_lgamma(B->nser + num + 1.0) - d_lgamma(num + 1.0));
            cancel = fmin(cancel, 1.45 * growth + 10.0);
        }
        return fmax(cancel, 0.0);
    }

    x = B->xi;

    /* the largest term of the series with +nu: the first one or the one
       at the peak */
    lmax = series_log_term(B, 0.0, 1);
    if (B->nser >= 1.0)
        lmax = fmax(lmax, series_log_term(B, B->nser, 1));

    if (kind == KIND_J || kind == KIND_I)
    {
        lfun = (num > zn) ? x : fabs(x);
    }
    else
    {
        /* the series with -nu: the terms decrease from the first one
           while k (|nu| - k) > |z|^2 / 4, and peak again beyond k = |nu| */
        double s = sqrt(num * num + zn * zn);

        lmax = fmax(lmax, series_log_term(B, 0.0, -1));
        lmax = fmax(lmax, series_log_term(B, 0.5 * (num + s), -1));
        if (zn < num)
            lmax = fmax(lmax, series_log_term(B, 0.5 * (num - sqrt(num * num - zn * zn)), -1));

        if (kind == KIND_K)
            lfun = (B->zr >= 0.0) ? -x : fabs(x);
        else
            lfun = fabs(x);

        /* (the combination of the two series is divided by sin(pi nu),
           which is large for complex orders: that is not cancellation) */
        lfun += fmax(3.14159265358979323846 * fabs(B->ni) - 0.6931471805599453, 0.0);

        /* Y at noninteger order is computed from J_nu and J_{-nu} by an
           implementation which raises the precision of each series by
           its own cancellation, and there is none between the two:
           nothing to add to the working precision (the two amounts
           enter the cost) */
        if (kind == KIND_Y && !B->intorder)
        {
            double lpos = series_log_term(B, 0.0, 1);
            if (B->nser >= 1.0)
                lpos = fmax(lpos, series_log_term(B, B->nser, 1));
            B->cancel_pos = fmax(1.45 * (lpos - ((num > zn) ? x : fabs(x))) + 6.0, 0.0);
            B->cancel_neg = fmax(1.45 * (lmax - lfun) + 6.0, 0.0);
            return 0.0;
        }
    }

    /* (the estimates are of leading order: a margin of one percent
       costs less than a second attempt) */
    return fmax(1.465 * (lmax - lfun) + 6.0, 0.0);
}

/* The number of terms of the series for J and I needed at the working
   precision wp: beyond the peak at nser, until the terms have decreased
   by wp bits.  The ratio of consecutive terms at the index k is
   |z|^2 / (4 k (k + |nu|)), and the terms decrease by about half the
   log of the ratio at the last index per term on average. */
static double
series_terms(const bessel_model * B, slong wp)
{
    double num = B->num, zn = B->zn, N, r, bits = 0.6931471805599453 * (double) wp;
    int i;

    if (!(num < 1e15 && zn < 1e15 && zn > 1e-300))
        return 1.4 * zn + (double) wp / 10.0 + 10.0;

    N = B->nser + 10.0 + 0.25 * (double) wp;
    for (i = 0; i < 2; i++)
    {
        r = log(4.0 * N * (N + num) / (zn * zn));
        N = B->nser + 2.0 + 2.0 * bits / fmax(r, 0.05);
    }

    return N + 5.0;
}

/* Rough cost estimates (in microseconds on a typical machine; only the
   ratios matter) of the series at the working precision wp, of the
   asymptotic expansions, and of the integration at precision prec: per
   descent path from a saddle, in N-fold double precision up to N = 4,
   and from there on growing with the exponent 1.6 up to about 1000
   bits and 2.7 beyond (measured up to 10000 bits; the number of nodes
   grows with the precision).  The cost per term of the series grows
   slowly with the precision (measured: like its cube root for real and
   its square root for complex arguments, up to some thousand bits). */
static double
series_term_cost(const bessel_model * B, double wp)
{
    return B->real ? 0.07 * cbrt(wp) : 0.03 * sqrt(wp);
}

static double
series_cost(int kind, const bessel_model * B, slong wp)
{
    double terms, c, per_term;

    /* Up to about 200 bits (where the integration runs in N-fold double
       precision), simpler estimates which were tuned there: ~1.4 |z| +
       wp/10 terms, plus ~|nu| for Y and K when |z| is not small
       compared to |nu| and for the limit formula at integer order. */
    if (!B->hp)
    {
        double p = (double) wp + ((kind == KIND_Y && !B->intorder) ? B->cancel_pos : 0.0);

        terms = 1.4 * B->zn + p / 10.0 + 10.0;
        if ((kind == KIND_Y || kind == KIND_K) && B->zn > 0.5 * B->num)
            terms += 1.5 * B->num;
        c = terms * (0.3 + 0.001 * p);
        if (B->intorder)
            c += 10.0 + 0.1 * p + 2.0 * (1.4 * B->zn + 10.0) * (0.3 + 0.001 * p);
        return c;
    }

    terms = series_terms(B, wp);
    per_term = series_term_cost(B, (double) wp);

    if (kind == KIND_Y && !B->intorder)
    {
        /* the series for J_nu and J_{-nu}, each at the precision raised
           by its cancellation; the terms of the second start decreasing
           at about (|nu| + sqrt(|nu|^2 + |z|^2)) / 2 if |z| is not
           small compared to |nu| */
        double tneg;

        tneg = (terms - B->nser) + ((B->zn > 0.5 * B->num)
                ? 0.5 * (B->num + sqrt(B->num * B->num + B->zn * B->zn)) : 0.0);
        return 4.0 + 0.02 * (double) wp * log2(2.0 + B->num)
            + terms * series_term_cost(B, (double) wp + B->cancel_pos)
            + tneg * series_term_cost(B, (double) wp + B->cancel_neg);
    }

    /* Y and K need the series with -nu as well, which has ~|nu| more
       terms before they start decreasing when |z| is not small compared
       to |nu| (the limit formula at integer order has a sum of |nu| terms) */
    if (kind == KIND_Y || kind == KIND_K)
    {
        terms *= 2.0;
        if (B->zn > 0.5 * B->num)
            terms += 1.5 * B->num;
    }

    /* (the prefactor: a power and a gamma function) */
    c = 2.0 + 0.01 * (double) wp * log2(2.0 + B->num) + terms * per_term;
    /* (the limit formula: a finite sum of |nu| terms as well) */
    if (B->intorder)
        c += 10.0 + 0.1 * (double) wp + 2.0 * (1.4 * B->zn + 10.0 + B->num) * per_term;
    return c;
}

/* log of the absolute value of the k-th term
   prod_{j <= k} (4 nu^2 - (2j-1)^2) / (k! (8z)^k) of the asymptotic
   expansions */
static double
asymp_log_term(double k, double num, double zn)
{
    double v;

    if (k <= 0.0)
        return 0.0;
    if (k < num)
        v = d_lgamma(num + k + 0.5) - d_lgamma(num - k + 0.5);
    else
        v = d_lgamma(num + k + 0.5) + d_lgamma(k - num + 0.5);
    return v - d_lgamma(k + 1.0) - k * log(2.0 * zn);
}

/* The asymptotic expansions sum_k a_k / (c z)^k: for large orders, the
   terms first grow, up to the index k1 = sqrt(|z|^2 + |nu|^2) - |z| (to
   about e^(|nu|^2 / (2 |z|)) if |nu| is small compared to |z|); they
   are smallest (about e^(-2|z|) for small orders) at the index
   k2 = sqrt(|z|^2 + |nu|^2) + |z|.  The sum is the function divided by
   its leading behaviour e^(+-iz) or e^(+-z), of size e^(-+d) with
   d = |Re xi| - |Im z| (|Re z| for I): 0 on the real axis, where the
   expansions lose the growth of their terms to cancellation, and more
   than that off the axis for the component with the small sum.  For K,
   the sum is of size e^(Re z - Re xi), which on the positive real axis
   is the size of the largest term (the terms have the same sign: no
   cancellation).  The error bounds of the implementation exceed the
   first omitted term by a factor of about e^(|nu|^2 / |z|), and need
   |z| > |nu| for K; the number of terms is limited to 50 + 10 times
   the working precision. */
static void
asymp_model(bessel_model * B, int kind)
{
    double num = B->num, zn = B->zn, s, peak, low, lsum;

    B->apen = 0.0;
    B->apeak = 0.0;

    /* (of no use for small |z|, and not needed where the simpler
       criterion is used and fails) */
    if (zn < 4.0 || (!B->hp && !(zn > 1.9 * num + 2.0)))
    {
        B->aloss = 0.0;
        B->areach = 0.0;
        return;
    }

    if (num <= 2.0)
    {
        B->aloss = 0.0;
        B->areach = 2.885 * zn;
        return;
    }

    if (!(num < 1e15 && zn < 1e15) || (kind == KIND_K && zn < num))
    {
        B->aloss = 0.0;
        B->areach = (zn > 1e15 * num) ? 2.885 * zn : 0.0;
        return;
    }

    /* (complex orders: the simpler criterion and estimates) */
    if (fabs(B->ni) > 0.05 * num)
    {
        B->apen = 1.45 * num * num / zn;
        B->aloss = B->apeak = 0.72 * num * num / zn;
        B->areach = (zn > 1.9 * num + 2.0) ? 2.885 * zn - B->apen : 0.0;
        return;
    }

    B->apen = 1.45 * num * num / zn;

    s = sqrt(zn * zn + num * num);
    peak = fmax(asymp_log_term(num * num / (s + zn), num, zn), 0.0);
    low = asymp_log_term(s + zn, num, zn);

    if (kind == KIND_K)
        lsum = fmin(B->zr - B->xi, peak);
    else
        lsum = -fabs(fabs(B->xi) - fabs((kind == KIND_I) ? B->zr : B->zi));

    B->apeak = 1.45 * peak;
    B->aloss = 1.45 * fmax(peak - lsum, 0.0);
    B->areach = 1.45 * (lsum - low) - B->apen;
}

/* whether the asymptotic expansion reaches wp bits, with a number of
   guard bits which is not out of proportion */
static int
asymp_converges(const bessel_model * B, slong wp)
{
    /* (the simpler criterion: beyond the turning point region,
       empirically from |z| ~ 1.9 |nu|, where the error is about
       e^(-2 |z|)) */
    if (!B->hp)
        return B->zn > 1.9 * B->num + 2.0 && 2.885 * B->zn >= (double) wp + 3.0;

    return B->areach >= (double) wp + 3.0 && B->aloss <= 3.0 * (double) wp + 100.0;
}

/* the number of terms of the expansion for wp bits: the implementation
   adds terms until they are smaller than the unit roundoff at its
   working precision (wp + aloss bits) divided by the factor in the
   error bound, that is, from the peak until the terms have decreased by
   apeak + wp + aloss + apen bits (by bisection) */
static double
asymp_terms(const bessel_model * B, slong wp)
{
    double num = B->num, zn = B->zn, bits = (double) wp + B->aloss + B->apen + B->apeak;
    double s, lo, hi, target;
    int i;

    if (num <= 2.0 || !(num < 1e15 && zn < 1e15))
    {
        /* n! / (2 |z|)^n = 2^(-bits): at most 0.4 terms per bit, fewer
           for large |z| */
        double n = 0.4 * bits;

        if (zn < 1e15)
            for (i = 0; i < 3; i++)
                n = fmin(0.6931471805599453 * bits / fmax(1.0 - log((n + 1.0) / (2.0 * zn)), 1.0), 0.4 * bits);
        return n;
    }

    s = sqrt(zn * zn + num * num);
    lo = num * num / (s + zn);
    hi = s + zn;
    target = fmax(asymp_log_term(lo, num, zn), 0.0) - 0.6931471805599453 * bits;

    lo = fmax(lo, 1.0);
    for (i = 0; i < 7; i++)
    {
        double mid = sqrt(lo * hi);
        if (asymp_log_term(mid, num, zn) > target)
            lo = mid;
        else
            hi = mid;
    }

    return hi;
}

/* the precision at which to evaluate the expansion for wp bits: with
   guard bits for its cancellation, and large enough for the number of
   terms to be allowed */
static slong
asymp_prec(const bessel_model * B, slong wp)
{
    double p = (double) wp + fmin(B->aloss, 3.0 * (double) wp + 100.0);

    if (B->num > 2.0)
        p = fmax(p, fmin(0.11 * asymp_terms(B, wp) + 20.0, 4.0 * (double) wp + 100.0));

    return (slong) p;
}

/* the cost of the expansion for wp bits (two expansions except for K) */
static double
asymp_cost(int kind, const bessel_model * B, slong wp)
{
    double N, wpe;

    if (!B->hp)
    {
        double t = fmin(B->num * B->num / B->zn, 1e9);
        return 30.0 + (0.4 * (double) wp + 1.4 * t) * (0.3 + 0.001 * ((double) wp + 0.72 * t));
    }

    N = asymp_terms(B, wp);
    wpe = fmax((double) wp + B->aloss, 0.11 * N);

    /* (measured up to 30000 bits: the cost per term grows linearly with
       the precision up to some hundred bits, like its square root up to
       some thousand, and slowly from there) */
    return 5.0 + (N + 20.0) * ((wpe < 350.0) ? 0.0024 * wpe
            : (wpe < 3500.0) ? 0.045 * sqrt(wpe) : 2.66 * pow(wpe / 3500.0, 0.15))
        * ((kind == KIND_K) ? 0.45 : 1.0);
}

static void
bessel_model_init(bessel_model * B, int kind, const acb_t nu, const acb_t z,
    int real_case, int intorder, slong goal)
{
    B->hp = (goal > 200);
    B->nr = arf_get_d(arb_midref(acb_realref(nu)), ARF_RND_NEAR);
    B->ni = arf_get_d(arb_midref(acb_imagref(nu)), ARF_RND_NEAR);
    B->zr = arf_get_d(arb_midref(acb_realref(z)), ARF_RND_NEAR);
    B->zi = arf_get_d(arb_midref(acb_imagref(z)), ARF_RND_NEAR);
    B->num = hypot(B->nr, B->ni);
    B->zn = hypot(B->zr, B->zi);
    B->real = real_case;
    B->intorder = intorder;
    B->nser = 0.5 * (-(B->num + 1.0) + sqrt((B->num + 1.0) * (B->num + 1.0) + B->zn * B->zn));
    B->lz = log(0.5 * B->zn);
    B->az = (B->zi == 0.0 && B->zr > 0.0) ? 0.0 : atan2(B->zi, B->zr);
    B->xi = 0.0;
    if (B->num > 2.0 && B->num < 1e100 && B->zn < 1e100 && B->zn > 1e-100
        && !((kind == KIND_J || kind == KIND_I) && B->nser < 1.0 && B->zn < 4.0))
    {
        if (kind == KIND_J || kind == KIND_Y)
            B->xi = _acb_hypgeom_bessel_debye_re(B->nr, B->ni, B->zr, B->zi);
        else
            B->xi = _acb_hypgeom_bessel_debye_re(B->nr, B->ni, -B->zi, B->zr);
    }
    B->cancel = series_cancellation(kind, B);

    /* Y at a noninteger order close to an integer:
       (J_nu cos(pi nu) - J_{-nu}) / sin(pi nu) loses the bits by which
       sin(pi nu) is small */
    if (kind == KIND_Y && !intorder && fabs(B->ni) < 0.25 && B->num < 1e15)
    {
        double d = fabs(B->nr - floor(B->nr + 0.5));
        d = hypot(d, B->ni);
        if (d < 0.25)
            B->cancel += fmin(-log2(fmax(3.2 * d, 1e-300)), 1000.0);
    }

    asymp_model(B, kind);
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
integration_cost(int kind, const bessel_model * B, double absA, slong prec, double other)
{
    double num = B->num, zn = B->zn;
    int real = B->real;
    static const double dcost[4] = { 30.0, 45.0, 100.0, 200.0 };
    double bits = (double) prec + 12.0 + log2(2.0 + zn + num), c;
    int npaths;

    if (bits <= 4 * 53)
        c = dcost[(int) ((bits - 1.0) / 53.0)];
    else if (bits <= 1024.0)
        c = 550.0 * pow(bits / 256.0, 1.6);
    else
        c = 550.0 * 9.19 * pow(bits / 1024.0, 2.7);

    /* one saddle for K and for J, I in the monotone regime; two for
       the oscillatory functions (one with real variables, where the
       contributions are complex conjugates, and one for J and I where
       the contribution e^(-|Re xi|) of the second is negligible
       compared to e^(|Re xi|)) */
    if (kind == KIND_K)
        npaths = 2;
    else if (kind == KIND_Y || zn > num)
        npaths = (real || (B->hp && kind != KIND_Y && 2.9 * fabs(B->xi) > bits + 20.0)) ? 2 : 4;
    else
        npaths = 2;

    /* per saddle (two paths): the planning of the contour and the
       expansion in place of the quadrature, where that matters (its
       estimate takes from a few microseconds to some hundred at
       thousands of bits, and is only computed if a cheap lower bound
       for it leaves the integration a chance against the other
       methods, of cost other) */
    if (absA > 0.0 && bits > 4 * 53
        && npaths * 0.5 * (45.0 + 0.1 * bits + _sd_expand_cost_lower(absA, real, prec))
            < fmin(other, npaths * c))
        c = fmin(c, 0.5 * (45.0 + 0.1 * bits + _sd_expand_cost(absA, num / absA, real, prec)));

    return npaths * c;
}

/* the expected cost of the series and expansions until success at the
   goal precision: the cheaper of the asymptotic expansion, where it
   converges, and the series at the precision raised by its
   cancellation; infinite where neither applies.  Also sets
   *asymp_first to whether the expansion is the cheaper one. */
static double
first_cost(int * asymp_first, int kind, const bessel_model * B,
    int series_plausible, slong goal)
{
    double ca = INFINITY, cs = INFINITY;
    slong swp = goal + 30 + (slong) fmin(B->cancel, 1e9);
    int conv = asymp_converges(B, goal);

    /* (the series has at least nser terms and the expansion at least
       0.35 per bit: there is often no need to count them) */
    if (B->hp && series_plausible && conv
        && 1.2 * B->nser * series_term_cost(B, (double) swp)
            > 4.0 * (5.0 + (0.4 * ((double) goal + 30.0 + B->aloss + B->apen + B->apeak) + 20.0)
                * 0.0024 * ((double) goal + 30.0 + B->aloss)))
        series_plausible = 0;

    if (series_plausible)
        cs = 1.2 * series_cost(kind, B, swp);
    if (conv && (!B->hp || 5.0 + (0.35 * ((double) goal + 30.0) + 20.0) * 0.02 * sqrt((double) goal + 30.0) < cs))
        ca = asymp_cost(kind, B, goal + 30);

    *asymp_first = (ca <= cs) && (ca < INFINITY);
    return fmin(ca, cs);
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
    bessel_model B;
    int asymp_first;
    slong wp_next, need_wp;
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
    intorder = (kind == KIND_Y || kind == KIND_K) && acb_contains_int(nu);
    nlim = intorder ? num : 0.0;
    bessel_model_init(&B, kind, nu, z, real_case, intorder, goal);
    nser = B.nser;

    series_plausible = (nser <= 10.0 * (16.0 * (double) goal + 480.0))
                        && (nlim <= 20.0 * (double) goal + 20000.0);

    /* precision for the midpoint evaluation */
    mprec = (acc_in < prec - 30) ? acc_in + 30 : prec;

    /* the integration goes first when it is expected to be cheaper
       than the series (at the precision raised by the cancellation it
       must absorb) or the asymptotic expansion, where that converges */
    int_ok = integration_applicable(nu, z);
    cancel = fmin(B.cancel, 1e9);
    c_first = first_cost(&asymp_first, kind, &B, series_plausible, goal);
    c_int = integration_cost(kind, &B, saddle_curvature(kind, nu, z), FLINT_MIN(prec, mprec), c_first);
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
    else if (series_plausible && !asymp_first)
        wp = goal + (slong) FLINT_MIN(cancel, 4.0 * (double) goal);

    for (wp = wp + 30, round = 0; wp <= 16 * goal + 480;
            wp = wp_next, round++)
    {
        int m, series_first, series_done = 0, asymp_done = 0;
        slong iter_start_acc = best_acc, awp, swp;

        need_wp = WORD_MAX;

        /* the cheaper of the series and the expansion comes first */
        series_first = !asymp_first;

        /* the expansion with guard bits for its own cancellation, and
           the series with those for its cancellation (which the
           working precision includes if the series is the first choice) */
        awp = inexact ? wp + (slong) FLINT_MIN(B.aloss, (double) wp + 100.0) : asymp_prec(&B, wp);
        swp = wp;
        if (!inexact && asymp_first && goal >= prec)
            swp += (slong) FLINT_MIN(cancel, 4.0 * (double) goal);

        /* First, unless the integration is expected to be cheaper, the
           series-and-asymptotics implementation, which is tuned for
           the easy cases, or the asymptotic expansion.  Only the
           implementation of J raises the precision by the cancellation
           of the series. */
        if (round == 0 && !int_first)
        {
            /* (guard bits for the rounding errors: a few more at high
               precision, where the sums are longer) */
            slong p0 = prec + 10 + FLINT_MIN(prec / 40, 15), p0a;
            int conv, direct, series0;

            /* The expansion directly: for K as soon as it converges
               (its series cancels heavily); where it is the cheaper
               method by the estimates (at low precision, where the
               simpler estimates are used, the implementations of J and
               I switch from the series by themselves instead); and for
               inexact input well beyond the turning point, where there
               is no cancellation to inflate the input radii. */
            conv = asymp_converges(&B, goal);
            direct = conv && (kind == KIND_K || (B.hp && asymp_first)
                    || (inexact && zn > 1.9 * num + 2.0))
                && !(kind == KIND_Y && !arb_is_positive(acb_realref(z)));
            p0a = inexact ? wp + (slong) FLINT_MIN(B.aloss, (double) goal + 100.0)
                          : asymp_prec(&B, p0);

            /* With exact input, when the series is the choice: the
               series itself for I and K (whose implementations only
               choose between the series and the expansion, by a simpler
               criterion and without guard bits for the latter), and
               also for J where its implementation would choose the
               expansion.  Otherwise the implementation of J raises the
               precision of the series by its cancellation. */
            series0 = !inexact && !direct && !asymp_first
                && (kind == KIND_I || kind == KIND_K
                    || (kind == KIND_J && zn >= 8.0 && 2.0 * zn >= (double) p0
                        && _acb_hypgeom_bessel_asymp_converges(nu, z, p0)));

            if ((kind != KIND_J || series0) && !direct
                && (series0 || !(conv && 2.0 * zn >= (double) p0)))
                p0 += (slong) FLINT_MIN(cancel, 4.0 * (double) goal);

            /* (with an inexact z, the cancellation of the series
               inflates the radius regardless of the precision, and so
               does the expansion close to where it stops converging;
               the evaluation at the midpoint is then better) */
            if (!(inexact && cancel >= (acb_is_exact(nu) ? 4.0 : 16.0)
                    && !(asymp_converges(&B, goal + 12) && zn > 1.9 * num + 2.0)))
            {
                if (direct)
                {
                    method_asymp(t, kind, scaled, nu, z, p0a);
                    asymp_done = 1;
                }
                else if (series0)
                {
                    method_series(t, kind, scaled, nu, z, p0);
                    series_done = 1;
                }
                else
                {
                    method_nointegration(t, kind, scaled, nu, z, inexact ? wp : p0);
                    /* (this is the series, for exact input) */
                    series_done = !inexact;
                }
                keep_best(best, &best_acc, t);
                /* (if only z is inexact, the evaluation at the midpoint
                   is cheap and nearly optimal: less slack) */
                if (best_acc >= ((inexact && acb_is_exact(nu)) ? target - 2 : accept))
                    goto done;

                /* a near miss with exact input is a matter of guard
                   bits: the accuracy is recovered bit for bit */
                if (!inexact && best_acc > 0)
                {
                    wp += (accept - best_acc) + 20;
                    asymp_done = series_done = 0;
                    awp = asymp_prec(&B, wp);
                    swp = wp;
                    if (asymp_first && goal >= prec)
                        swp += (slong) FLINT_MIN(cancel, 4.0 * (double) goal);
                }
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
            && (int_first || c_int <= (series_first ? series_cost(kind, &B, wp) : asymp_cost(kind, &B, wp))))
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
            if (int_ok && !int_tried && c_int <= 1.5 * (c_first + 150.0))
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
                    && !(inexact && round == 0 && kind == KIND_Y)
                    && !series_done)
                {
                    /* (the integration first if the series has become
                       the more expensive of the two) */
                    TRY_INTEGRATION_IF(series_cost(kind, &B, swp));
                    method_series(t, kind, scaled, nu, z, swp);
                    keep_best(best, &best_acc, t);
                    series_done = 1;
                    NOTE_ATTEMPT();
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
                    NOTE_ATTEMPT();
                }
            }
            else if (!asymp_done && (inexact || round >= 2 || asymp_converges(&B, accept)))
            {
                /* (with exact input, not where the expansion is not
                   expected to get there, unless nothing else has) */
                TRY_INTEGRATION_IF(asymp_cost(kind, &B, wp));
                method_asymp(t, kind, scaled, nu, z, awp);
                keep_best(best, &best_acc, t);
                NOTE_ATTEMPT();
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
            && (best_acc > 0 || !series_plausible)
            && need_wp == WORD_MAX)
            break;

        /* the next working precision: where an attempt of this round
           has fallen short by a known number of bits (with exact
           input), that many more; otherwise doubled */
        wp_next = 2 * wp + 30;
        if (need_wp < wp_next)
            wp_next = FLINT_MAX(need_wp, wp + 30);
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

        bessel_model B;
        int asymp_first;

        bessel_model_init(&B, KIND_Y, nu, z, bessel_real_case(KIND_Y, nu, z), intorder, goal);

        double c_first = first_cost(&asymp_first, KIND_Y, &B, series_plausible, goal);

        if (integration_cost(KIND_Y, &B, saddle_curvature(KIND_Y, nu, z), goal, c_first) < c_first
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
