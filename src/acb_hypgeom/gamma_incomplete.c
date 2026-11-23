/*
    Copyright (C) 2026 Fredrik Johansson
    Developed using Claude Fable 5

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/* Automatic algorithm selection for the incomplete gamma functions with
   complex parameter and argument, combining the convergent and
   asymptotic series (reliable and fast except in the large complex
   transition region |s| ~ |z|) with rigorous numerical integration
   along steepest-descent contours (gamma_incomplete_integration.c).
   This follows the same scheme as for the Bessel functions (bessel.c):
   the cheapest method by estimates of the costs is tried first, at the
   precision raised by the estimated cancellation; each result is
   verified a posteriori and the best kept; with exact input, an
   attempt which falls short by a known number of bits is repeated with
   that many more, and otherwise the working precision doubles; the
   integration, which evaluates the function on the input balls
   directly, is tried before any attempt which is expected to cost
   more; and for inexact inputs there is also the evaluation at the
   midpoints with the propagated error added from derivative bounds
   (the z-derivative -+z^(s-1) e^(-z) in closed form, which makes this
   the first choice when only z is inexact; the s-derivative bounded
   along the steepest-descent contour).  For real s and z > 0 where the
   function is monotone, the values at two corners of the input box
   give its range (see monotone_route), and an inexact s is otherwise
   first tried with Cauchy's estimates and a crude bound for the
   function (see midpoint_route_cauchy).  The easiest cases (a small
   |z|, or the upper function far out) are computed by the series and
   the expansions before any estimates; a ball z which meets the cut
   from both sides is left to them as well; and for input which is as
   accurate as the precision, up to 1000 bits, a first result which is
   a few more bits short of the target is accepted. */

#include <math.h>

#include "arb.h"
#include "acb.h"
#include "arb_hypgeom.h"
#include "acb_hypgeom.h"
#include "acb_hypgeom/impl.h"

/* the number of bits of s if it is large: acb_gamma and acb_pow lose
   them for complex s */
static slong
large_bits(const acb_t s)
{
    slong e, f;

    if (!acb_is_finite(s))
        return 0;
    e = arf_abs_bound_lt_2exp_si(arb_midref(acb_realref(s)));
    f = arf_abs_bound_lt_2exp_si(arb_midref(acb_imagref(s)));
    e = FLINT_MAX(e, f);
    return (e >= 44) ? FLINT_MIN(e, 4000) : 0;
}

/* the radius of an enclosure, as a single magnitude */
static void
enclosure_rad(mag_t r, const acb_t x)
{
    mag_hypot(r, arb_radref(acb_realref(x)), arb_radref(acb_imagref(x)));
}

static void
keep_best(acb_t best, slong * best_acc, const acb_t cand)
{
    slong a;

    if (!acb_is_finite(cand))
        return;

    /* unclamped: deeply negative accuracies carry the recovery
       signal that drives the precision-doubling loop through large
       cancellation */
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

/* Heuristic estimate of the relative accuracy attainable from inexact
   input balls: the relative condition numbers are roughly |z| + |s| + 1
   with respect to z and |s| (1 + |log(2(|s|+1)/|z|)|) with respect
   to s.  This is only used to decide when to stop trying harder; it
   does not affect the correctness of the enclosure. */
static slong
attainable_accuracy(const acb_t s, const acb_t z, int upper, double sn, double zn,
    slong prec)
{
    slong acc = prec, t;

    if (!acb_is_exact(z))
    {
        t = acb_rel_accuracy_bits(z);
        if (t > 0)
        {
            /* the relative condition number in z, |z^s e^(-z) / F|, is
               taken to be at most |z| + |s|; for real s >= 1 and z > 0
               it is at most z for the upper function (Gamma(s, z) >=
               z^(s-1) e^(-z)) and s for the lower one (gamma(s, z) >=
               z^s e^(-z) / s), and on the side of the median of the
               function which holds most of Gamma(s) (z <= s - 1 for
               the upper one, z >= s for the lower one) at most
               2 z^s e^(-z) / Gamma(s): Gamma(1000, z) at z = 905
               rounded has 0.1, where |z| + |s| would take 11 bits.
               (Not for input as accurate as the precision up to 1000
               bits, of which a first result may miss a few bits more:
               a target closer to the precision would only add the
               value at the midpoint there.) */
            double lc = log2(1.0 + zn + sn);

            if ((prec > 1000 || t < prec - 16)
                && acb_is_real(s) && acb_is_real(z) && arb_is_positive(acb_realref(z)))
            {
                double sr = arf_get_d(arb_midref(acb_realref(s)), ARF_RND_NEAR);
                double zr = arf_get_d(arb_midref(acb_realref(z)), ARF_RND_NEAR);

                if (sr >= 1.0 && sr < 1e15 && zr > 0.0 && zr < 1e15)
                {
                    double lb = 1.0 + (sr * log(zr) - zr - lgamma(sr)) * 1.4426950408889634;
                    double lcr = upper ? log2(zr) : log2(sr);

                    if (upper ? (zr <= sr - 1.0) : (zr >= sr))
                        lcr = FLINT_MIN(lcr, lb);

                    lc = fmin(lc, log2(1.0 + exp2(fmin(lcr, 1e6))));
                }
            }

            t -= (slong) FLINT_MIN(lc, 1e6);
        }
        acc = FLINT_MIN(acc, t);
    }

    if (!acb_is_exact(s))
    {
        t = upper ? acb_rel_one_accuracy_bits(s) : acb_rel_accuracy_bits(s);
        if (t > 0)
            t -= (slong) fmin(log2(1.0 + sn * (1.0 + fabs(
                log(2.0 * (sn + 1.0)) - log(fmax(zn, 1e-300))))), 1e6);
        acc = FLINT_MIN(acc, t);
    }

    return acc;
}

/* rigorous |d/dz| bound over the input balls: |z^(s-1) e^(-z)|; also
   returns whether the bound is within a factor two of the smallest
   value over the balls */
static int
dz_bound(mag_t Dz, const acb_t s, const acb_t z)
{
    acb_t t, u;
    mag_t lo;
    int tight;

    acb_init(t);
    acb_init(u);
    mag_init(lo);
    acb_sub_ui(t, s, 1, 32);
    acb_pow(u, z, t, 32);
    acb_neg(t, z);
    acb_exp(t, t, 32);
    acb_mul(u, u, t, 32);
    acb_get_mag(Dz, u);
    acb_get_mag_lower(lo, u);
    mag_mul_2exp_si(lo, lo, 1);
    tight = mag_cmp(Dz, lo) <= 0;
    acb_clear(t);
    acb_clear(u);
    mag_clear(lo);
    return tight;
}

/* whether the series route for the upper function goes through the
   series for the lower function (see method_series) */
static int
upper_by_lower_series(const acb_t s)
{
    return arf_cmp_si(arb_midref(acb_realref(s)), 1) >= 0;
}

/* gamma(s, z) = s^-1 z^s 1F1(s, s+1, -z) = s^-1 z^s e^-z 1F1(1, s+1, z)
   by the series itself (acb_hypgeom_m would switch to the asymptotic
   expansion for large |z| regardless of s, which fails and is
   expensive for |s| comparable to |z|; the expansion is the other
   route).  The terms of the first series grow to e^|z| and its sum is
   of size e^(-Re z); the terms z^k / (s+1)_k of the second do not grow
   if |z| < |s|, and otherwise grow to at most e^|z| with a sum of size
   e^(Re z).  The form with the smaller cancellation is used (see
   series_cancellation for the working precision this takes). */
static void
lower_series(acb_t res, const acb_t s, const acb_t z, slong wp)
{
    acb_t s1, w, t;
    double zre, zim, zn, sn;
    int kummer;

    zre = arf_get_d(arb_midref(acb_realref(z)), ARF_RND_NEAR);
    zim = arf_get_d(arb_midref(acb_imagref(z)), ARF_RND_NEAR);
    zn = hypot(zre, zim);
    sn = hypot(arf_get_d(arb_midref(acb_realref(s)), ARF_RND_NEAR),
               arf_get_d(arb_midref(acb_imagref(s)), ARF_RND_NEAR));
    kummer = (zre > 0.0 || zn < sn);

    acb_init(s1);
    acb_init(w);
    acb_init(t);
    acb_add_ui(s1, s, 1, wp);

    if (kummer)
    {
        /* sum of z^k / (s+1)_k */
        acb_hypgeom_pfq_direct(t, NULL, 0, s1, 1, z, -1, wp);
        acb_neg(w, z);
        acb_exp(w, w, wp);
        acb_mul(t, t, w, wp);
    }
    else
    {
        /* sum of (s)_k (-z)^k / ((s+1)_k k!) */
        acb_struct bb[2];
        bb[0] = *s1;
        acb_init(bb + 1);
        acb_one(bb + 1);
        acb_neg(w, z);
        acb_hypgeom_pfq_direct(t, s, 1, bb, 2, w, -1, wp);
        acb_clear(bb + 1);
    }

    acb_pow(res, z, s, wp);
    acb_mul(t, t, res, wp);
    acb_div(res, t, s, wp);
    acb_clear(s1);
    acb_clear(w);
    acb_clear(t);
}

/* The series-class route, at the working precision wp (which must
   include the cancellation, see series_cancellation): the series for
   the lower function, and for the upper function with Re(s) >= 1 that
   series subtracted from Gamma(s); otherwise the implementation for
   the upper function, which bundles the convergent series with their
   elementary and singular special cases. */
static void
method_series(acb_t res, const acb_t s, const acb_t z, int upper,
    slong wp)
{
    acb_t g;

    if (upper && !upper_by_lower_series(s))
    {
        acb_hypgeom_gamma_upper_nointegration(res, s, z, 0, wp);
        return;
    }

    acb_init(g);
    lower_series(res, s, z, wp);

    if (upper)
    {
        if (acb_is_finite(res))
        {
            acb_gamma(g, s, wp);
            acb_sub(res, g, res, wp);
        }
        else
            acb_hypgeom_gamma_upper_nointegration(res, s, z, 0, wp);
    }
    else if (!acb_is_finite(res))
    {
        /* the lower series can fail where the upper series
           converges; reflect through Gamma(s) */
        acb_hypgeom_gamma_upper_nointegration(res, s, z, 0, wp);
        acb_gamma(g, s, wp);
        acb_sub(res, g, res, wp);
    }

    acb_clear(g);
}

/* the asymptotic route: the (self-limiting, divergent) 2F0 expansion
   of the upper gamma; the lower gamma by reflection, whose
   cancellation at small |z| is absorbed by the doubling working
   precision and is complementary to the series route */
static void
method_asymp(acb_t res, const acb_t s, const acb_t z, int upper,
    slong wp, slong lsmall)
{
    /* (for the lower function, Gamma(s,z) is smaller than Gamma(s) by
       about lsmall bits, and needs that many fewer) */
    acb_hypgeom_gamma_upper_asymp(res, s, z, 0,
        upper ? wp : FLINT_MAX(32, wp - lsmall + 10));
    if (!upper)
    {
        acb_t g;
        acb_init(g);
        acb_gamma(g, s, wp);
        acb_sub(res, g, res, wp);
        acb_clear(g);
    }
}

static void
method_nointegration(acb_t res, const acb_t s, const acb_t z, int upper,
    int regularized, slong wp)
{
    if (upper)
        acb_hypgeom_gamma_upper_nointegration(res, s, z, regularized, wp);
    else
        acb_hypgeom_gamma_lower_nointegration(res, s, z, regularized, wp);
}

static void gamma_incomplete_core(acb_t res, const acb_t s, const acb_t z, int upper, slong prec);

/* (in the main loop of gamma_incomplete_core(): the integration, if it
   has not been tried and is not expected to cost more than the attempt
   which would come next) */
#define TRY_INTEGRATION_IF(next_cost) \
    do { \
        if (int_ok && !int_tried && !inexact && c_int <= (next_cost)) \
        { \
            int_tried = 1; \
            if (_acb_hypgeom_gamma_incomplete_integration(t, s, z, upper, prec)) \
            { \
                int_done = 1; \
                keep_best(best, &best_acc, t); \
                if (best_acc >= accept) \
                    goto done; \
            } \
        } \
    } while (0)

/* (an attempt at the working precision wp which has fallen short of
   the accuracy to accept by a known number of bits) */
#define NOTE_ATTEMPT() \
    do { \
        if (!inexact && acb_is_finite(t)) \
        { \
            slong a__ = acb_rel_accuracy_bits(t); \
            if (a__ > 0 && a__ < accept) \
                need_wp = FLINT_MIN(need_wp, wp + (accept - a__) + 20); \
        } \
    } while (0)

/* (the bits by which a first result on rounded input may fall short,
   see gamma_incomplete_core) */
#define LAZY_SLACK 4

static int
z_meets_cut(const acb_t s, const acb_t z)
{
    return arb_contains_zero(acb_imagref(z)) && !arb_is_exact(acb_imagref(z))
        && !arb_is_positive(acb_realref(z))
        && !(acb_is_int(s) && arb_is_positive(acb_realref(s)));
}

/* Evaluation at the input midpoints with the propagated error
   sup |dF/dz| rz + sup |dF/ds| rs added: the z-derivative is
   -+z^(s-1) e^(-z), and the s-derivative is bounded along the
   steepest-descent contour.  The result is indeterminate on failure
   (in particular, for a ball z crossing the branch cut); returns 2 if
   the error bound is about as good as possible (s exact, and the
   derivative bound within a factor two of the derivative), and
   otherwise 1 on success. */
static int
gamma_incomplete_midpoint(acb_t res, const acb_t s, const acb_t z, int upper, slong prec)
{
    mag_t D, err, rr;
    acb_t sm, zm;
    int status = 1;

    if (z_meets_cut(s, z))
    {
        acb_indeterminate(res);
        return 0;
    }

    mag_init(D);
    mag_init(err);
    mag_init(rr);
    acb_init(sm);
    acb_init(zm);

    if (!acb_is_exact(z))
    {
        if (dz_bound(D, s, z))
            status = 2;
        enclosure_rad(rr, z);
        mag_mul(err, rr, D);
    }

    if (!acb_is_exact(s))
    {
        status = 1;
        if (!_acb_hypgeom_gamma_incomplete_integration_bounds(rr, D, s, z, upper))
            mag_inf(D);
        enclosure_rad(rr, s);
        mag_addmul(err, rr, D);
    }

    if (mag_is_finite(err))
    {
        acb_get_mid(sm, s);
        acb_get_mid(zm, z);
        gamma_incomplete_core(res, sm, zm, upper, prec);
        acb_add_error_mag(res, err);
    }
    else
    {
        acb_indeterminate(res);
        status = 0;
    }

    mag_clear(D);
    mag_clear(err);
    mag_clear(rr);
    acb_clear(sm);
    acb_clear(zm);

    return status;
}

static void gamma_incomplete(acb_t res, const acb_t s, const acb_t z, int upper,
    int regularized, slong prec);

/* Real s and positive z: where the function is monotone in both
   variables (or in the one which is inexact), its range over the input
   box is the interval between the values at two corners, which are
   computed with exact input: the optimal enclosure, whatever the size
   of the box (and for the regularized functions without the dependency
   between the incomplete and the complete gamma function).  With
   x > 0:

   - Gamma(s,x) = int_x^inf t^(s-1) e^(-t) dt is decreasing in x, and
     increasing in s for x >= 1 (the derivative has the factor log t);
   - gamma(s,x) = int_0^x (s > 0) is increasing in x, and decreasing
     in s for x <= 1;
   - Q(s,x) = Gamma(s,x) / Gamma(s) (s > 0), the probability that a
     gamma variable with the shape s exceeds x, is decreasing in x and
     increasing in s; P = 1 - Q;
   - Gamma(s,x) x^(-s) = int_1^inf u^(s-1) e^(-x u) du is decreasing
     in x and increasing in s;
   - gamma^*(s,x) = e^(-x) sum_k x^k / Gamma(s+k+1) (s > 0) is
     decreasing in x (it is int_0^1 u^(s-1) e^(-x u) du / Gamma(s)),
     and decreasing in s for s >= 1/2, where the gamma functions
     increase.

   dirs: +1 (increasing), -1 (decreasing), 0 (unknown) in s and x.
   Returns 0 if the corners do not apply. */
#define MONOTONE_MAX_ACC 64

static int
monotone_dirs(int * ds, int * dx, const acb_t s, const acb_t z, int upper, int regularized)
{
    const arb_struct * sr = acb_realref(s), * zr = acb_realref(z);
    arb_t t;
    int ok;

    if (!acb_is_real(s) || !acb_is_real(z) || !arb_is_positive(zr))
        return 0;
    if (regularized < 0 || regularized > 2)
        return 0;
    if (!(upper && regularized != 1) && !arb_is_positive(sr))
        return 0;

    arb_init(t);
    *ds = *dx = 0;

    if (upper)
    {
        *dx = -1;
        if (regularized == 0)
        {
            arb_sub_ui(t, zr, 1, 64);
            *ds = arb_is_nonnegative(t) ? 1 : 0;
        }
        else
            *ds = 1;
    }
    else if (regularized == 0)
    {
        *dx = 1;
        arb_sub_ui(t, zr, 1, 64);
        *ds = arb_is_nonpositive(t) ? -1 : 0;
    }
    else if (regularized == 1)
    {
        *dx = 1;
        *ds = -1;
    }
    else
    {
        *dx = -1;
        arb_set_d(t, 0.5);
        arb_sub(t, sr, t, 64);
        *ds = arb_is_nonnegative(t) ? -1 : 0;
    }

    ok = (*ds != 0 || acb_is_exact(s)) && (*dx != 0 || acb_is_exact(z));
    arb_clear(t);
    return ok;
}

static int
monotone_route(acb_t res, const acb_t s, const acb_t z, int upper, int regularized,
    int ds, int dx, slong mprec, slong prec)
{
    arf_t sl, sh, zl, zh;
    acb_t a, b, lo, hi;
    mag_t t;
    int ok, pass;

    arf_init(sl); arf_init(sh); arf_init(zl); arf_init(zh);
    acb_init(a); acb_init(b); acb_init(lo); acb_init(hi);

    /* (the corners, rounded outwards; exact inputs as they are) */
    if (acb_is_exact(s))
    {
        arf_set(sl, arb_midref(acb_realref(s)));
        arf_set(sh, sl);
    }
    else
    {
        arb_get_lbound_arf(sl, acb_realref(s), mprec);
        arb_get_ubound_arf(sh, acb_realref(s), mprec);
    }
    if (acb_is_exact(z))
    {
        arf_set(zl, arb_midref(acb_realref(z)));
        arf_set(zh, zl);
    }
    else
    {
        arb_get_lbound_arf(zl, acb_realref(z), mprec);
        arb_get_ubound_arf(zh, acb_realref(z), mprec);
    }

    /* (still in the domain where the directions hold: the conditions
       s > 0, s >= 1/2, x >= 1, x <= 1 are closed or hold for the
       rounded corners too, except at the boundary of s > 0, x > 0) */
    ok = arf_sgn(zl) > 0 && arf_is_finite(zh) && arf_is_finite(sl) && arf_is_finite(sh)
        && ((upper && regularized != 1) || arf_sgn(sl) > 0);
    if (ok && !upper && regularized == 2 && ds != 0)
        ok = arf_cmp_d(sl, 0.5) >= 0;
    if (ok && regularized == 0 && ds != 0)
        ok = upper ? (arf_cmp_si(zl, 1) >= 0) : (arf_cmp_si(zh, 1) <= 0);

    mag_init(t);

    /* (once more with the full precision if the range is narrower
       than the accuracy of the input suggests) */
    for (pass = 0; pass < 2 && ok; pass++)
    {
        arb_set_arf(acb_realref(a), (ds >= 0) ? sh : sl);
        arb_set_arf(acb_realref(b), (dx >= 0) ? zh : zl);
        gamma_incomplete(hi, a, b, upper, regularized, mprec);
        arb_set_arf(acb_realref(a), (ds >= 0) ? sl : sh);
        arb_set_arf(acb_realref(b), (dx >= 0) ? zl : zh);
        gamma_incomplete(lo, a, b, upper, regularized, mprec);

        ok = acb_is_finite(lo) && acb_is_finite(hi);
        if (ok)
        {
            arb_union(acb_realref(res), acb_realref(lo), acb_realref(hi), mprec);
            arb_zero(acb_imagref(res));
        }

        mag_add(t, arb_radref(acb_realref(lo)), arb_radref(acb_realref(hi)));
        mag_mul_2exp_si(t, t, 3);
        if (mprec >= prec || !ok || mag_cmp(t, arb_radref(acb_realref(res))) <= 0)
            break;
        mprec = prec;
    }

    mag_clear(t);

    arf_clear(sl); arf_clear(sh); arf_clear(zl); arf_clear(zh);
    acb_clear(a); acb_clear(b); acb_clear(lo); acb_clear(hi);
    return ok;
}

/* An inexact s, cheaply (see midpoint_route_cauchy in bessel.c): with M
   a bound for |F(s', z)| on the disk |s' - m| <= R around the midpoint
   m of s, for z in its ball, and for s in its ball of radius r,

       |F(s, z) - F(m, z)| <= r |F'(m, z)| + r^2 M / (R - r)^2,
       |F'(m, z)| <= |F(m + h, z) - F(m, z)| / h + h M / (R - h)^2,

   with F(m, z) and F(m + h, z) enclosed on the ball z by the methods
   for an exact s.  The bounds for Re(z) > 0, with s = sigma + i tau:
   on the ray t = z + u, u >= 0, |t^(s-1)| <= e^(pi |tau| / 2) |t|^(sigma-1)
   and Re(z) <= |t| <= |z| + u, so that

       |Gamma(s,z)| <= e^(pi |tau| / 2) e^(|z| - Re z) Gamma(sigma, |z|)

   if sigma >= 1 and |z| >= 1 on the disk and the ball, and otherwise,
   with p = max(sigma - 1, 0) and q = max(1 - sigma, 0),

       |Gamma(s,z)| <= e^(pi |tau| / 2) max(1, 1 / Re z)^q e^(1 + |z| - Re z) Gamma(p + 1, 1 + |z|);

   for the lower function (Re(s) > 0 on the disk),
   gamma(s,z) = z^s int_0^1 u^(s-1) e^(-z u) du gives |gamma(s,z)| <= |z^s| / sigma,
   and also |gamma(s,z)| <= |Gamma(s)| + |Gamma(s,z)|. */
static int
midpoint_route_cauchy(acb_t res, const acb_t s, const acb_t z, int upper,
    double num, double zn, slong mprec)
{
    mag_t r, R, M, A, t, u, h, hlo, zmax, zmin;
    acb_t sm, sb, w, dd;
    arb_t x, y, g;
    arf_t xl, smax, smin, e;
    double Rd;
    slong k = 12;
    int status = 0;
    int ok, second = 0;

    acb_indeterminate(res);

    if (acb_is_exact(s) || !arb_is_positive(acb_realref(z)))
        return 0;

    mag_init(r); mag_init(R); mag_init(M); mag_init(A); mag_init(t); mag_init(u); mag_init(h);
    mag_init(zmax); mag_init(zmin); mag_init(hlo);
    acb_init(sm); acb_init(sb); acb_init(w); acb_init(dd);
    arb_init(x); arb_init(y); arb_init(g);
    arf_init(xl); arf_init(smax); arf_init(smin); arf_init(e);

    acb_get_mid(sm, s);

    Rd = fmin(0.25, 1.0 / (1.0 + fabs(log(fmax(fmax(zn, num), 1e-300)))));
    if (!upper)
        Rd = fmin(Rd, 0.7 * arf_get_d(arb_midref(acb_realref(sm)), ARF_RND_DOWN));
    ok = (Rd > 1e-30);
    if (ok)
        mag_set_ui_2exp_si(R, 1, (slong) floor(log2(Rd)));
    enclosure_rad(r, s);
    mag_mul_2exp_si(t, r, 2);
    ok = ok && (mag_cmp(t, R) <= 0);

    /* the disk, in a box */
    acb_set(sb, sm);
    arb_add_error_mag(acb_realref(sb), R);
    arb_add_error_mag(acb_imagref(sb), R);
    arb_get_lbound_arf(smin, acb_realref(sb), 30);
    arb_get_ubound_arf(smax, acb_realref(sb), 30);
    arb_get_lbound_arf(xl, acb_realref(z), 30);
    acb_get_mag(zmax, z);
    acb_get_mag_lower(zmin, z);

    ok = ok && arf_sgn(xl) > 0 && (upper || arf_sgn(smin) > 0);

    if (ok)
    {
        /* the upper function */
        arf_set_mag(e, zmax);
        if (arf_cmp_si(smin, 1) >= 0 && mag_cmp_2exp_si(zmin, 0) >= 0)
        {
            /* e^(|z| - Re z) Gamma(sigma, |z|) */
            arb_set_arf(x, smax);
            arb_set_arf(y, e);
            arb_hypgeom_gamma_upper(g, x, y, 0, 30);
            arb_sub_arf(y, y, xl, 30);
        }
        else
        {
            /* max(1, 1 / Re z)^q e^(1 + |z| - Re z) Gamma(p + 1, 1 + |z|) */
            arb_set_arf(x, smax);
            arb_sub_ui(x, x, 1, 30);
            if (arb_is_negative(x))
                arb_zero(x);
            arb_add_ui(x, x, 1, 30);
            arb_set_arf(y, e);
            arb_add_ui(y, y, 1, 30);
            arb_hypgeom_gamma_upper(g, x, y, 0, 30);
            arb_sub_arf(y, y, xl, 30);
            if (arf_cmp_si(smin, 1) < 0 && arf_cmp_si(xl, 1) < 0)
            {
                arb_t q, v;
                arb_init(q); arb_init(v);
                arb_set_arf(q, smin);
                arb_sub_ui(q, q, 1, 30);
                arb_neg(q, q);
                arb_set_arf(v, xl);
                arb_inv(v, v, 30);
                arb_pow(v, v, q, 30);
                arb_mul(g, g, v, 30);
                arb_clear(q); arb_clear(v);
            }
        }
        arb_exp(y, y, 30);
        arb_mul(g, g, y, 30);
        /* e^(pi |tau| / 2) */
        arb_get_mag(t, acb_imagref(sb));
        arf_set_mag(e, t);
        arb_const_pi(y, 30);
        arb_mul_arf(y, y, e, 30);
        arb_mul_2exp_si(y, y, -1);
        arb_exp(y, y, 30);
        arb_mul(g, g, y, 30);
        arb_get_mag(M, g);
        ok = mag_is_finite(M);

        if (ok && !upper)
        {
            /* |Gamma(s)| + |Gamma(s,z)|, or |z^s| / sigma */
            acb_gamma(w, sb, 30);
            acb_get_mag(t, w);
            mag_add(M, M, t);

            acb_pow(w, z, sb, 30);
            acb_get_mag(A, w);
            arb_set_arf(x, smin);
            arb_get_mag_lower(t, x);
            mag_div(A, A, t);
            mag_min(M, M, A);
            ok = mag_is_finite(M);
        }
    }

    if (ok)
    {
        gamma_incomplete_core(res, sm, z, upper, mprec);

        /* the step: see bessel.c */
        if (acb_is_finite(res))
        {
            k = 12;

            acb_get_mag_lower(t, res);
            if (!mag_is_zero(t))
            {
                mag_div(u, M, t);
                mag_div(u, u, R);
                arf_set_mag(e, u);
                k = 4 + arf_abs_bound_lt_2exp_si(e);
                k = FLINT_MAX(k, 12);
                k = FLINT_MIN(k, 4 * mprec + 1000);
            }
            mag_mul_2exp_si(h, R, -k);
            mag_mul_2exp_si(t, r, 3);
            if (mag_cmp(t, h) > 0)
            {
                arf_set_mag(e, t);
                mag_set_ui_2exp_si(h, 1, arf_abs_bound_lt_2exp_si(e));
            }
            mag_mul_2exp_si(t, R, -2);
            second = (mag_cmp(h, t) <= 0);
        }

        if (second)
        {
            /* the second point m' = m + h, rounded: with the bounds
               hlo <= |m' - m| <= h */
            slong wp2 = mprec + k + 200;

            arf_set_mag(e, h);
            arb_add_arf(acb_realref(w), acb_realref(sm), e, wp2);
            arb_set(acb_imagref(w), acb_imagref(sm));
            acb_get_mid(w, w);
            acb_sub(dd, w, sm, wp2);
            acb_get_mag(h, dd);
            acb_get_mag_lower(hlo, dd);
            mag_mul_2exp_si(t, R, -2);
            second = !mag_is_zero(hlo) && (mag_cmp(h, t) <= 0);
            if (second)
                acb_set(sm, w);
        }

        if (second)
        {
            /* (its error is multiplied by r / h: fewer bits will do) */
            {
                slong p2 = mprec;

                if (!mag_is_zero(r))
                {
                    arf_set_mag(e, r);
                    p2 = arf_abs_bound_lt_2exp_si(e);
                    arf_set_mag(e, hlo);
                    p2 = mprec - (arf_abs_bound_lt_2exp_si(e) - p2) + 20;
                    p2 = FLINT_MAX(p2, 40);
                    p2 = FLINT_MIN(p2, mprec);
                }
                gamma_incomplete_core(w, sm, z, upper, p2);
            }
            acb_sub(w, w, res, mprec);
            acb_get_mag(A, w);
            mag_div(A, A, hlo);
            mag_mul(A, A, r);                 /* the term of the first order */
            mag_sub_lower(t, R, h);
            mag_mul(t, t, t);
            mag_div(u, M, t);
            mag_mul(u, u, h);
            mag_mul(u, u, r);
            mag_sub_lower(t, R, r);
            mag_mul(t, t, t);
            mag_div(t, M, t);
            mag_mul(t, t, r);
            mag_mul(t, t, r);
            mag_add(u, u, t);                 /* the terms with M */
            mag_add(M, A, u);
            /* (about as good as possible: the terms with M are small) */
            mag_mul_2exp_si(u, u, 2);
            status = (mag_cmp(u, A) <= 0) ? 2 : 1;
        }
        else
        {
            /* r M / (R - r) */
            mag_sub_lower(t, R, r);
            mag_mul(M, M, r);
            mag_div(M, M, t);
            status = 1;
        }

        if (acb_is_finite(res) && mag_is_finite(M))
            acb_add_error_mag(res, M);
        else
        {
            acb_indeterminate(res);
            status = 0;
        }
    }

    mag_clear(r); mag_clear(R); mag_clear(M); mag_clear(A); mag_clear(t); mag_clear(u); mag_clear(h);
    mag_clear(zmax); mag_clear(zmin); mag_clear(hlo);
    acb_clear(sm); acb_clear(sb); acb_clear(w); acb_clear(dd);
    arb_clear(x); arb_clear(y); arb_clear(g);
    arf_clear(xl); arf_clear(smax); arf_clear(smin); arf_clear(e);
    return status;
}

/* Rough cost estimates (in microseconds on a typical machine; only the
   ratios matter) of the series at the working precision wp, of the
   asymptotic expansion (see asymp_terms), and of the integration at
   precision prec (see bessel.c).  The terms z^k / (s+1)_k of the
   series decrease from k ~ |z| - Re(s), by a factor |z| / |s + k|. */
static double
series_terms(double sre, double zn, slong wp)
{
    double s1 = fmax(sre, 0.0) + 1.0, k0, k1, b, m, L, f, lb, bits;
    int i;

    /* (an upper bound: the terms decrease by a factor of two or more
       from k0 on) */
    k0 = fmax(2.0 * zn - s1 + 1.0, 0.0);
    m = k0 + (double) wp / fmax(1.0, log2((s1 + k0) / (zn + 1e-300)));

    if (!(zn > 1e-300 && zn < 1e15 && s1 < 1e15))
        return m + 5.0;

    /* The terms grow up to k1 = |z| - s1 (if that is positive), and
       m terms further they have decreased by
       sum_{j < m} log((b + j) / |z|) ~ (b + m) log((b + m) / |z|)
       - b log(b / |z|) - m, with b = max(s1, |z|): the m for which
       this is wp bits, by Newton's method from above (the function is
       convex). */
    k1 = fmax(zn - s1, 0.0);
    b = fmax(s1, zn);
    lb = b * log(b / zn);
    bits = 0.6931471805599453 * (double) wp;
    m = fmax(m - k1, 1.0);

    for (i = 0; i < 4; i++)
    {
        L = log((b + m) / zn);
        f = (b + m) * L - lb - m - bits;
        if (!(f > 0.0) || !(L > 1e-6))
            break;
        m = fmax(m - f / L, 0.0);
    }

    return k1 + m + 5.0;
}

static double
series_cost(double sre, double zn, slong wp, int intpar)
{
    double terms = series_terms(sre, zn, wp), p = (double) wp, c;

    /* (at higher precision, the powers and the gamma function cost as
       much as some hundred terms, and the terms are cheap, in particular
       with an integer parameter) */
    if (wp <= 128)
        return 5.0 + terms * (0.12 + 0.0005 * p);

    c = (0.05 + 0.00012 * p) * ((intpar & 1) ? 0.4 : 1.0);

    /* (the series 1F1(s, s+1, -z), which lower_series sums for z in the
       left half-plane with |z| >= |s|: 0.35 to 0.45 microseconds per
       term with a short s, for 100 to 4000 terms from 350 to 3400
       bits; taken into account for an integer s, where the terminating
       expansion is the alternative) */
    if (intpar == 3)
        c = fmax(c, 0.44);

    /* (a complex z of full length is multiplied in at every term: with
       s = 10^5, z = s (1 + i) / sqrt(2) rounded, 0.54, 1.9 and 12.7
       microseconds per term at 333, 3333 and 33333 bits) */
    if (intpar & 4)
        c = fmax(c, 0.4 + 3.5e-4 * p);

    return 5.0 + 0.2 * (p - 128.0) + terms * c;
}

/* The number of terms (1-s)_k / z^k of the asymptotic expansion needed
   for wp bits: the k-th ratio is about hypot(|s|, k) / |z|, so that n
   terms gain about n log(|z| / hypot(|s|, n/e)) (n! / z^n for small
   |s|); infinite if the terms start growing before that. */
static double
asymp_terms(double num, double zn, slong wp, double nterm)
{
    double n = 0.0, L;
    int i;

    for (i = 0; i < 4; i++)
    {
        L = log(zn / hypot(num, 0.368 * n + 1.0));
        if (!(L > 0.01))
        {
            n = INFINITY;
            break;
        }
        n = 0.693 * (double) wp / L;
    }

    if (!(hypot(num, n) < zn))
        n = INFINITY;

    /* (for a positive integer s, the expansion terminates after s
       terms: Gamma(s, z) = (s-1)! e^(-z) sum_{k<s} z^k / k!) */
    return (nterm > 0.0) ? fmin(n, nterm) : n;
}

static double
asymp_cost(double num, double zn, slong wp, double nterm)
{
    return 4.0 + asymp_terms(num, zn, wp, nterm) * (0.15 + 0.0005 * (double) wp);
}

static double
integration_cost(double num, double zn, int real, int upper, slong prec)
{
    static const double dcost[4] = { 30.0, 45.0, 100.0, 200.0 };
    double bits = (double) prec + 12.0 + log2(2.0 + zn + num);

    /* (measured at 64 to 256 bits, with the nodes at hand: 30 to 330
       microseconds for the upper function and 50 to 600 for the lower
       one at 64 bits, for real and complex s and z of sizes up to 1000;
       two to four times that at 128 bits) */
    if (bits <= 4 * 53)
        return (upper ? 4.0 : 8.0) * dcost[(int) ((bits - 1.0) / 53.0)];
    /* (measured: growing with the exponent 1.6 up to about 1000 bits and
       2.7 beyond; a third of that with real variables) */
    if (bits <= 1024.0)
        return (real ? 0.3 : 0.75) * 550.0 * pow(bits / 256.0, 1.6);
    return (real ? 0.25 : 0.75) * 550.0 * 9.19 * pow(bits / 1024.0, 2.7);
}

/* The cancellation (in bits) which the series route must absorb through
   its working precision: that of the series for the lower function in
   the form used (see lower_series), and for the upper function computed
   from it also that of Gamma(s) - gamma(s,z), where
   |Gamma(s,z)| ~ |z|^(Re(s)-1) e^(-Re z) for |z| beyond Re(s).  For the
   implementation of the upper function used otherwise, the terms grow
   to e^|z| while the function is of size e^(-Re z) or more. */
static double
series_cancellation(int upper, int by_lower, double sre, double sn, double zn, double zre)
{
    double c;

    if (upper && !by_lower)
    {
        c = 1.45 * (zn - fabs(zre));
        if (zre > 0.0 && zn > sre)
            c += fmax(0.0, 1.45 * (lgamma(fmax(sre, 1.0)) - (sre - 1.0) * log(zn) + zre));
        return c;
    }

    if (zre > 0.0 || zn < sn)
    {
        c = (zn < sn) ? 0.0 : 1.45 * (zn - zre);

        /* (The sum of z^k / (s+1)_k is s e^z z^(-s) gamma(s, z).  For
           real s > 0 its size e^(Re z) holds where the integral for
           gamma(s, z) is dominated by its saddle point at t = s - 1,
           which the ray to z passes for Re(z) > s - 1; for 0 < Re(z) < s
           the end point z dominates, the sum is about s / |s - z|, and
           the cancellation is the growth of the largest term,
           e^((|z| - s) - s log(|z| / s)) for |z| > s, with the
           number of terms of that size, up to sqrt(s): with
           s = 10^5, z = s (1 + i) / sqrt(2), a few bits where
           1.45 (|z| - Re(z)) would take 42500.) */
        if (zn >= sn && zre > 0.0 && zre < sre && sn < 1.001 * sre)
        {
            double m = (zn > sn) ? (zn - sn) - sn * log(zn / sn) : 0.0;

            c = fmin(c, 1.45 * m + log2(1.0 + sqrt(sn)) + 4.0);
        }
    }
    else
        c = 1.45 * (zn + zre);

    if (upper && zn > sre)
        c += fmax(0.0, 1.45 * (lgamma(fmax(sre, 1.0)) - (sre - 1.0) * log(zn) + zre));

    return c;
}

static void gamma_incomplete_core_u(acb_t res, const acb_t s, const acb_t z,
    int upper, slong prec, slong uprec);

static void
gamma_incomplete_core(acb_t res, const acb_t s, const acb_t z, int upper,
    slong prec)
{
    gamma_incomplete_core_u(res, s, z, upper, prec, prec);
}

/* (uprec: the precision of the caller, where prec has been raised for
   the regularization) */
static void
gamma_incomplete_core_u(acb_t res, const acb_t s, const acb_t z, int upper,
    slong prec, slong uprec)
{
    double num, zn, sre, c_int, c_first, cancel, nterm = 0.0, aloss = 0.0, afac = 0.0;
    int asymp_ok, asymp_first, intpar;
    slong wp_next, need_wp, awpl, lsm = 0, amin = 0;
    slong wp, acc_in, goal, target, accept, best_acc = -WORD_MAX;
    acb_t best, t;
    int inexact, round, series_plausible, mid_tried = 0, int_ok, int_tried = 0, int_done = 0, int_first;
    int cheap_tried = 0, lazy, easy = 0;

    if (!acb_is_finite(s) || !acb_is_finite(z))
    {
        method_nointegration(res, s, z, upper, 0, prec);
        return;
    }

    zn = hypot(arf_get_d(arb_midref(acb_realref(z)), ARF_RND_NEAR),
               arf_get_d(arb_midref(acb_imagref(z)), ARF_RND_NEAR));
    num = hypot(arf_get_d(arb_midref(acb_realref(s)), ARF_RND_NEAR),
                arf_get_d(arb_midref(acb_imagref(s)), ARF_RND_NEAR));

    /* inputs not faithfully representable as doubles: the cost
       estimate below is meaningless, and the series and asymptotic
       expansions handle all such inputs quickly */
    if (!isfinite(zn) || !isfinite(num) || zn > 1e300 || num > 1e300)
    {
        method_nointegration(res, s, z, upper, 0, prec);
        return;
    }

    inexact = !acb_is_exact(s) || !acb_is_exact(z);

    /* a ball z which meets the cut (-inf, 0] from both sides (unless s
       is a positive integer, for which there is no cut): no method
       gives more than the jump of the function */
    if (z_meets_cut(s, z))
    {
        method_nointegration(res, s, z, upper, 0, prec);
        return;
    }

    /* the goal is capped by the relative accuracy of the inputs:
       nothing more is attainable, and hopeless input balls then
       terminate after the first cheap attempt; the target also
       accounts (heuristically) for the conditioning */
    /* (Gamma(s,z) is entire in s, while gamma(s,z) has a pole at 0) */
    acc_in = upper ? acb_rel_one_accuracy_bits(s) : acb_rel_accuracy_bits(s);
    target = acb_rel_accuracy_bits(z);
    acc_in = FLINT_MIN(acc_in, target);
    acc_in = FLINT_MAX(acc_in, 0);
    goal = FLINT_MIN(prec, acc_in);
    target = attainable_accuracy(s, z, upper, num, zn, prec);
    target = FLINT_MIN(goal, target);

    /* Accept results this close to the target.  For inexact input the
       target is only a heuristic estimate of what is attainable; a
       result short of it is compared with the integration on the
       input balls. */
    accept = target - (inexact ? 8 : 5);

    /* Input which is as accurate as the precision (the result of
       rounding at this precision, as a rule): up to 1000 bits, a first
       result on the input balls which is a few more bits short is
       accepted, since the function at the midpoints with a bound for
       the variation costs several times as much there, and the caller
       has more bits for less by raising the precision. */
    lazy = inexact && uprec <= 1000 && acc_in >= uprec - 16;

    /* nothing is attainable: direct interval propagation is as good as
       anything, quickly */
    if (goal <= 0 || target < 4)
    {
        /* (unless the function is monotone: the corners of the box) */
        int ds, dx;

        if (!(inexact && monotone_dirs(&ds, &dx, s, z, upper, 0)
                && monotone_route(res, s, z, upper, 0, ds, dx, FLINT_MIN(prec, 64), FLINT_MIN(prec, 64))))
            method_nointegration(res, s, z, upper, 0, prec);
        return;
    }

    /* The easiest cases, before any estimates (which cost a
       microsecond, as much as these functions at low precision): the
       first attempt as it would be made below.  If the result is not
       accurate, everything follows as usual, with this result for the
       first attempt.

       A small |z|: a few terms of a series.

       The upper function with |z| >= 4 (|s| + prec): the terms of its
       expansion decrease by a factor of four or more each, down to
       the precision. */
    {
        slong p = prec + 10 + FLINT_MIN(prec / 40, 15);

        if (zn <= 0.0625)
        {
            if (inexact)
                method_nointegration(res, s, z, upper, 0, FLINT_MIN(prec, goal) + 30);
            else
                method_series(res, s, z, upper, p);
            easy = 1;
        }
        else if (upper && !inexact && zn >= 4.0 * (num + (double) prec + 30.0))
        {
            method_asymp(res, s, z, upper, p, 0);
            easy = 2;
        }

        if (easy && acb_is_finite(res)
            && acb_rel_accuracy_bits(res) >= ((inexact && acb_is_exact(s)) ? target - 2 : accept - (lazy ? LAZY_SLACK : 0)))
        {
            acb_set_round(res, res, prec);
            return;
        }
    }

    /* with s near the negative real axis and large, the terms of the
       convergent series do not decrease before the Pochhammer symbol
       has passed the negative integers (the other parameter regions
       are handled by the method selection in the series code) */
    {
        double sim = arf_get_d(arb_midref(acb_imagref(s)), ARF_RND_NEAR);
        double nterms = 0.0;

        sre = arf_get_d(arb_midref(acb_realref(s)), ARF_RND_NEAR);

        if (sre < 0.0 && fabs(sim) < -sre)
            nterms = -sre;

        series_plausible = (nterms <= 10.0 * (16.0 * (double) goal + 480.0));
    }

    /* The integration (planned in double precision) goes first when it
       is expected to be cheaper than the asymptotic expansion, where
       that converges (|z| well beyond |s|), or the series at the
       precision raised by its cancellation of ~1.45 (|z| - |Re z|)
       bits. */
    int_ok = arf_cmpabs_2exp_si(arb_midref(acb_realref(z)), 860) < 0
        && arf_cmpabs_2exp_si(arb_midref(acb_imagref(z)), 860) < 0
        && arf_cmpabs_2exp_si(arb_midref(acb_realref(s)), 860) < 0
        && arf_cmpabs_2exp_si(arb_midref(acb_imagref(s)), 860) < 0
        && zn > 1e-270 && !acb_contains_zero(z);
    c_int = integration_cost(num, zn, acb_is_real(s) && acb_is_real(z), upper, (acc_in < prec - 30) ? acc_in + 30 : prec);
    intpar = acb_is_exact(s) && acb_is_int(s);
    /* (and the form of the series, see series_cost) */
    if ((!upper || upper_by_lower_series(s))
        && !(arf_sgn(arb_midref(acb_realref(z))) > 0 || zn < num))
        intpar |= 2;
    if (!acb_is_real(z) && acb_bits(z) > 64)
        intpar |= 4;

    /* A positive integer s: the expansion of the upper function
       terminates, Gamma(s, z) = (s-1)! e^(-z) sum_{j<s} z^j / j!.  The
       terms of this partial sum of the exponential series are largest
       at j = min(s - 1, |z|), and the sum is of the size of e^z or of
       its last term, whichever is larger: for s well beyond |z| it
       loses about 1.45 (|z| - Re z) bits, and nothing for s <= |z|.
       Where the expansion does not converge, all the terms are summed,
       which the implementation does at a precision of more than s
       bits.  (Timings for s = 1000: 0.20 microseconds per term at 1100
       bits and 0.23 at 3333 bits for a short z; 0.42 and 1.7 for a
       complex z of full length.) */
    if (acb_is_exact(s) && acb_is_int(s) && sre >= 1.0 && sre <= 1e9)
    {
        double zre = arf_get_d(arb_midref(acb_realref(z)), ARF_RND_NEAR);

        nterm = sre;
        if (!isfinite(asymp_terms(num, zn, goal + 8, 0.0)))
        {
            /* (none for s - 1 <= |z|, where the largest term is the last
               one: computed as the difference of two equal expressions,
               it can be a rounding error on x87, where one of them may
               be kept with extended precision) */
            if (zn > 1e-300 && sre - 1.0 > zn)
            {
                double lzn = log(zn);

                aloss = 1.45 * fmax(zn * lzn - lgamma(zn + 1.0)
                    - fmax(zre, (sre - 1.0) * lzn - lgamma(sre)), 0.0);
            }

            amin = (slong) sre + 2;

            /* (the measured cost of the terms for s <= |z| + 1, where
               they do not grow; the general estimate otherwise) */
            if (sre - 1.0 <= zn)
                afac = (acb_bits(z) <= 128) ? 0.03 : acb_is_real(z) ? 0.6 : 1.0;
        }
    }

/* the cost of the expansion at w bits, and its precision */
#define ACOST(w) ((afac > 0.0) \
    ? 4.0 + nterm * (0.19 + (0.1 + 0.0004 * (double) FLINT_MAX((w), amin)) * afac) \
    : asymp_cost(num, zn, FLINT_MAX((w), amin), nterm))
#define AWP(w) FLINT_MAX((w), amin)

    asymp_ok = isfinite(asymp_terms(num, zn, goal + 8, nterm))
        && (upper || aloss == 0.0);

    /* The lower function from the expansion, as Gamma(s) - Gamma(s,z):
       where Gamma(s,z) is small compared to Gamma(s), by lsmall bits,
       that many fewer bits of it are needed (and the expansion need not
       converge to more). */
    awpl = goal + 8;
    if (!upper)
    {
        double zre = arf_get_d(arb_midref(acb_realref(z)), ARF_RND_NEAR);

        if (zn > sre && zre > 0.0 && sre < 1e15 && zn < 1e15)
        {
            double lsmall = 1.45 * (lgamma(fmax(sre, 1.0)) - (sre - 1.0) * log(zn) + zre);

            if (lsmall > 0.0)
            {
                awpl = goal + 8 - (slong) fmin(lsmall, 1e9);
                awpl = FLINT_MAX(awpl, 16);
                asymp_ok = asymp_ok || isfinite(asymp_terms(num, zn, awpl, nterm));
                /* (only for exact input: the bits are needed otherwise) */
                if (asymp_ok && !inexact)
                    lsm = (slong) fmin(lsmall, 1e9);
            }
        }
    }
    /* (For |z| below Re(s) the lower function is small compared to
       Gamma(s), by about Gamma(s + 1) e^z / z^s: the difference
       Gamma(s) - Gamma(s,z) cancels, and the series, whose terms
       decrease from the first one, is the method.) */
    if (!upper && asymp_ok && sre >= 1.0 && zn < sre && zn > 1e-300 && sre < 1e15)
    {
        double zre = arf_get_d(arb_midref(acb_realref(z)), ARF_RND_NEAR);

        if (1.45 * (lgamma(sre + 1.0) + zre - sre * log(zn)) > 4.0)
            asymp_ok = 0;
    }

    cancel = series_cancellation(upper, upper_by_lower_series(s), sre, num, zn,
                arf_get_d(arb_midref(acb_realref(z)), ARF_RND_NEAR));
    cancel = FLINT_MIN(cancel, 4.0 * (double) goal);
    /* the cheaper of the expansion, where it converges, and the series */
    c_first = 1.5 * series_cost(sre, zn, goal + 30 + (slong) cancel, intpar);
    asymp_first = 0;
    if (asymp_ok)
    {
        double ca = ACOST(awpl + (slong) aloss);
        if (ca <= c_first)
        {
            c_first = ca;
            asymp_first = 1;
        }
    }
    int_first = int_ok && (c_int < c_first);
    /* (elementary closed forms) */
    if (acb_is_exact(s) && acb_is_int(s)
        && arf_cmp_si(arb_midref(acb_realref(s)), 1) >= 0
        && arf_cmp_si(arb_midref(acb_realref(s)), 3) <= 0)
        int_first = 0;

    acb_init(best);
    acb_init(t);
    acb_indeterminate(best);

    /* the algorithms in order of cheapness, with doubling working
       precision: the series route absorbs any cancellation given
       enough extra precision, and the asymptotic expansion is
       self-limiting.  With inexact input, the goal may be below prec;
       the series can still lose up to ~1.44 |z| bits to cancellation
       (which is only recovered by increasing the working precision for
       exact input) */
    wp = goal;
    if (goal < prec)
        wp = FLINT_MIN(prec, goal + (slong) FLINT_MIN((double) goal, 1.45 * zn));
    else if (!asymp_first)
        wp = goal + (slong) cancel;

    for (wp = wp + 30, round = 0; wp <= 16 * goal + 480;
            wp = wp_next, round++)
    {
        int m, series_first, series_done = 0, asymp_done = 0;
        slong iter_start_acc = best_acc, awp, swp;

        need_wp = WORD_MAX;

        /* the cheaper of the series and the expansion comes first; the
           expansion with guard bits for its cancellation, and the
           series with those for its own (which the working precision
           includes if the series is the first choice) */
        series_first = !asymp_first;
        awp = AWP(wp + (slong) aloss);
        swp = wp;
        if (!inexact && asymp_first && goal >= prec)
            swp += (slong) cancel;

        /* First, unless the integration is expected to be cheaper, the
           series-and-asymptotics implementation, which is tuned for
           the easy cases. */
        if (round == 0 && (!int_first || easy))
        {
            /* (with an inexact z, the cancellation of the series
               inflates the radius regardless of the precision) */
            if (easy || !(inexact && acb_is_exact(s) && cancel >= 4.0 && !asymp_ok))
            {
                /* (guard bits for the rounding errors: a few more at
                   high precision, where the sums are longer) */
                slong p0 = prec + 10 + FLINT_MIN(prec / 40, 15);

                /* (the expansion directly where it is the choice: the
                   thresholds of the implementations differ, and for
                   inexact input with no more bits than the input is
                   good for; with exact input, the series route with its
                   guard bits, and otherwise the implementations, which
                   limit the precision by the accuracy of the input) */
                if (easy)
                {
                    /* (the attempt before the estimates) */
                    acb_set(t, res);
                    asymp_done = (easy == 2);
                    series_done = (easy == 1) && !inexact;
                }
                else if (asymp_first)
                {
                    method_asymp(t, s, z, upper, AWP((inexact ? FLINT_MIN(wp, goal + 20) : p0) + (slong) aloss), lsm);
                    asymp_done = 1;
                }
                else if (!inexact)
                {
                    method_series(t, s, z, upper, p0 + (slong) cancel);
                    series_done = 1;
                }
                else
                    method_nointegration(t, s, z, upper, 0, wp);
                keep_best(best, &best_acc, t);
                /* (if only z is inexact, the evaluation at the midpoint
                   is cheap and nearly optimal: less slack) */
                if (best_acc >= ((inexact && acb_is_exact(s)) ? target - 2 : accept - (lazy ? LAZY_SLACK : 0)))
                    goto done;

                /* a near miss with exact input is a matter of guard
                   bits: the accuracy is recovered bit for bit */
                if (!inexact && best_acc > 0)
                {
                    wp += (accept - best_acc) + 20;
                    awp = AWP(wp + (slong) aloss);
                    swp = wp + ((asymp_first && goal >= prec) ? (slong) cancel : 0);
                    asymp_done = series_done = 0;
                }
            }

            /* inexact z only: the evaluation at the midpoint with the
               derivative in closed form is cheap and nearly optimal */
            if (inexact && acb_is_exact(s))
            {
                int status = gamma_incomplete_midpoint(t, s, z, upper,
                                (acc_in < prec - 30) ? acc_in + 30 : prec);

                mid_tried = 1;
                keep_best(best, &best_acc, t);
                if (best_acc >= accept || status == 2)
                    goto done;
            }
        }

        /* an inexact s which the first attempt has not settled: the
           corners of the box where the function is monotone (the
           optimal enclosure), and otherwise the midpoint of s with a
           cheap bound for the derivative; both at the cost of two
           evaluations with an exact s, unless the integration on the
           input balls costs less (the evaluations cost more than
           estimated for complex arguments) */
        if (inexact && !acb_is_exact(s) && !int_first && !cheap_tried)
        {
            int ds, dx;
            slong mprec = (acc_in < prec - 30) ? acc_in + 30 : prec;

            cheap_tried = 1;

            /* (the corners for wide boxes, where two evaluations at the
               accuracy of the input are cheap; for narrow balls at a
               high precision, one evaluation and a cheap radius are
               better) */
            if (acc_in <= MONOTONE_MAX_ACC && monotone_dirs(&ds, &dx, s, z, upper, 0))
            {
                if (monotone_route(t, s, z, upper, 0, ds, dx, mprec, prec))
                {
                    keep_best(best, &best_acc, t);
                    goto done;
                }
            }
            else if (2.0 * ((acb_is_real(s) && acb_is_real(z)) ? 1.0 : 3.0) * c_first < c_int)
            {
                int status = midpoint_route_cauchy(t, s, z, upper, num, zn, mprec);

                /* (final if the bound is about as good as possible) */
                keep_best(best, &best_acc, t);
                if (best_acc >= accept || status == 2)
                    goto done;
            }
        }

        /* the integration, first when expected to be cheapest overall,
           or once it is cheaper than the next round of the series and
           expansions */
        if (int_ok && !int_tried
            && (int_first || c_int <= (series_first ? series_cost(sre, zn, swp, intpar) : ACOST(awp))))
        {
            int_tried = 1;
            if (_acb_hypgeom_gamma_incomplete_integration(t, s, z, upper, prec))
            {
                int_done = 1;
                /* (this is as good as it gets for inexact input) */
                keep_best(best, &best_acc, t);
                if (best_acc >= accept || (inexact && acb_rel_accuracy_bits(t) > 0))
                    goto done;
            }
        }

        for (m = 0; m < 2; m++)
        {
            if ((m == 0) == series_first)
            {
                if (!series_done)
                {
                    /* (the integration first if the series has become
                       the more expensive of the two) */
                    TRY_INTEGRATION_IF(series_cost(sre, zn, swp, intpar));
                    method_series(t, s, z, upper, swp);
                    keep_best(best, &best_acc, t);
                    NOTE_ATTEMPT();
                }
            }
            else if (!asymp_done && (inexact || round >= 2 || asymp_ok))
            {
                /* (with exact input, not where the expansion is not
                   expected to get there, unless nothing else has) */
                TRY_INTEGRATION_IF(ACOST(awp));
                method_asymp(t, s, z, upper, awp, lsm);
                keep_best(best, &best_acc, t);
                NOTE_ATTEMPT();
            }
            if (best_acc >= accept)
                goto done;
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
                if (_acb_hypgeom_gamma_incomplete_integration(t, s, z, upper, prec))
                {
                    int_done = 1;
                    /* (this is as good as it gets for inexact input) */
                    keep_best(best, &best_acc, t);
                    if (best_acc >= accept || acb_rel_accuracy_bits(t) > 0)
                        goto done;
                }
            }

            gamma_incomplete_midpoint(t, s, z, upper, (acc_in < prec - 30) ? acc_in + 30 : prec);
            keep_best(best, &best_acc, t);

            /* the midpoint result is as good as this route gets; if it
               falls far short (the derivative bounds can be pessimistic),
               integration with the input balls may do better, at a
               relatively high cost */
            if (best_acc >= FLINT_MIN(accept, target / 2))
                goto done;

            /* (more precision does not improve the results on the
               input balls) */
            break;
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
       still spoil the series; the midpoint evaluation is cheap */
    if (best_acc < accept && inexact && !mid_tried)
    {
        gamma_incomplete_midpoint(t, s, z, upper, prec);
        keep_best(best, &best_acc, t);
    }

    if (best_acc < accept && int_ok && !int_tried
        && _acb_hypgeom_gamma_incomplete_integration(t, s, z, upper, prec))
    {
        int_done = 1;
        keep_best(best, &best_acc, t);
    }

    /* the real-variable integral representation (for parameters too
       large for the contours to be planned) */
    if (best_acc < accept && !int_done
        && upper && acb_is_real(s) && acb_is_real(z)
        && arf_cmp_2exp_si(arb_midref(acb_realref(z)), -16) > 0
        && arf_cmpabs_2exp_si(arb_midref(acb_realref(s)), 60) < 0
        && arf_cmpabs_2exp_si(arb_midref(acb_realref(z)), 60) < 0)
    {
        acb_zero(t);
        arb_hypgeom_gamma_upper_integration(acb_realref(t), acb_realref(s),
            acb_realref(z), 0, prec + 30);
        keep_best(best, &best_acc, t);
    }

    /* last resort: the series at the target precision */
    if (best_acc < accept)
    {
        method_nointegration(t, s, z, upper, 0, prec);
        keep_best(best, &best_acc, t);
    }

done:
    acb_set_round(res, best, prec);
    acb_clear(best);
    acb_clear(t);
}

/* the incomplete gamma functions are real for real s on the positive
   real z-axis (at z = 0 the finite cases are real as well), and for
   positive integer s on the whole real axis; the computed imaginary
   part is then a ball around the exact value 0 */
static int
gamma_real_case(const acb_t s, const acb_t z)
{
    if (!acb_is_real(s) || !acb_is_real(z))
        return 0;
    return arb_is_nonnegative(acb_realref(z))
        || (acb_is_int(s) && arb_is_positive(acb_realref(s)));
}

/* plain: the function on the input balls from the unregularized one,
   without the complementary function and the corners */
static void
gamma_incomplete_x(acb_t res, const acb_t s, const acb_t z, int upper,
    int regularized, slong prec, int plain, slong uprec)
{
    acb_t g;
    slong wp;
    int real_case;

    if (res == s || res == z)
    {
        acb_init(g);
        gamma_incomplete_x(g, s, z, upper, regularized, prec, plain, uprec);
        acb_swap(res, g);
        acb_clear(g);
        return;
    }

    real_case = gamma_real_case(s, z);

    /* special values at z = 0 and non-finite input are handled by
       the series code */
    if (!acb_is_finite(s) || !acb_is_finite(z) || acb_contains_zero(z))
    {
        method_nointegration(res, s, z, upper, regularized, prec);
        goto finish;
    }

    /* Gamma(s, z) for s = 1, 2, 3: the elementary closed forms
       e^(-z) (1, z + 1, z^2 + 2 z + 2), before any estimates (unless
       the polynomial has lost bits close to one of its roots). */
    if (upper && acb_is_exact(s) && acb_is_int(s)
        && arf_cmp_si(arb_midref(acb_realref(s)), 1) >= 0
        && arf_cmp_si(arb_midref(acb_realref(s)), 3) <= 0)
    {
        slong zacc = acb_is_exact(z) ? prec + 4 : acb_rel_accuracy_bits(z);

        acb_hypgeom_gamma_upper_nointegration(res, s, z, regularized, prec + 4);

        /* (on a ball z: within a few bits of its accuracy) */
        if (acb_rel_accuracy_bits(res) >= FLINT_MIN(prec, zacc) - (acb_is_exact(z) ? 2 : 6))
        {
            acb_set_round(res, res, prec);
            goto finish;
        }
    }

    /* compensation for the log-scale of Gamma(s) and z^s; capped, since
       no accuracy is attainable when the exponent is huge anyway */
    wp = prec;
    if (regularized)
    {
        mag_t msn;
        slong cap = 4 * prec + 64, e = 0;

        mag_init(msn);
        acb_get_mag(msn, s);
        if (mag_cmp_2exp_si(msn, 1) > 0)
            e = (mag_cmp_2exp_si(msn, cap) >= 0) ? cap : MAG_EXP(msn);
        wp = prec + 15 + e;
        mag_clear(msn);
    }

    /* Q = 1 - P and P = 1 - Q for real s > 0 and z > 0, where the
       other function is small (the lower one for z < s): it is needed
       with fewer bits, by its size, which is estimated from the leading
       terms P ~ z^s e^(-z) / Gamma(s + 1) and Q ~ z^(s-1) e^(-z) / Gamma(s). */
    if (!plain && regularized == 1 && real_case && arb_is_positive(acb_realref(s))
        && arb_is_positive(acb_realref(z)))
    {
        double sd = arf_get_d(arb_midref(acb_realref(s)), ARF_RND_NEAR);
        double zd = arf_get_d(arb_midref(acb_realref(z)), ARF_RND_NEAR);
        double lsm = 0.0, lsum = 0.0;

        if (sd < 1e15 && zd < 1e15 && zd > 1e-300 && ((zd < sd) == (upper != 0)))
        {
            if (upper)
            {
                lsm = sd * log(zd) - zd - lgamma(sd + 1.0);
                lsum = log((sd + 1.0) / fmax(sd + 1.0 - zd, sqrt(sd)));
            }
            else
            {
                lsm = (sd - 1.0) * log(zd) - zd - lgamma(sd);
                lsum = log(zd / fmax(zd - sd + 1.0, sqrt(sd)));
            }
            lsm *= 1.4427;
            lsum *= 1.4427;
        }

        /* (The other function is that small where the leading term is,
           times the factor lsum of the sum of the terms, (s + 1) /
           (s + 1 - z) for P and z / (z - s + 1) for Q, which is about
           sqrt(s) at z = s: there neither function is small.  With an
           exact s the other function is computed anyway, since it costs
           less there for large s (P(10^4, 10^4): 0.2 ms instead of 0.4
           at 64 bits); with an inexact s the function itself has the
           classical implementation below, which costs less than the
           other one on the balls: P(s, z) at s, z = 1000 (1 +/- 2^-33325)
           took 124 ms instead of 47 at 33333 bits.) */
        if (lsm < -3.0 && (acb_is_exact(s) || lsm + lsum < -3.0))
        {
            slong p2 = prec + (slong) fmax(lsm, -1e9) + 30, acc, acc2;

            acb_init(g);
            gamma_incomplete_x(g, s, z, !upper, 1, FLINT_MAX(p2, 32), 1, uprec);

            if (acb_is_finite(g))
            {
                acb_sub_ui(g, g, 1, prec + 10);
                acb_neg(g, g);
                acc = acb_rel_accuracy_bits(s);
                acc2 = acb_rel_accuracy_bits(z);
                acc = FLINT_MIN(acc, acc2);
                acc = FLINT_MIN(prec, FLINT_MAX(acc, 0));
                /* (as good as it gets for the input) */
                if (acb_rel_accuracy_bits(g) >= acc - 8)
                {
                    acb_set_round(res, g, prec);
                    acb_clear(g);
                    goto finish;
                }
            }

            acb_clear(g);
        }
    }

    /* The regularized functions with an inexact s.  First the classical
       implementation, which handles the regularization by itself; if
       its result falls short of the accuracy which the input permits
       (for its conditioning with respect to z) by more than a few
       bits, for real s and z > 0 where the function is monotone, the
       corners of the box (for wide boxes: see gamma_incomplete_core),
       where Gamma(s,z) and Gamma(s) are not taken on the same ball. */
    if (!plain && regularized && !acb_is_exact(s))
    {
        int ds, dx, mono;
        slong acc0, acc;

        acc0 = acb_rel_accuracy_bits(s);
        acc = acb_rel_accuracy_bits(z);
        acc0 = FLINT_MIN(acc0, acc);
        acc0 = FLINT_MAX(acc0, 0);
        acc = attainable_accuracy(s, z, upper,
            fabs(arf_get_d(arb_midref(acb_realref(s)), ARF_RND_NEAR)),
            fabs(arf_get_d(arb_midref(acb_realref(z)), ARF_RND_NEAR)), prec);

        mono = (acc0 <= MONOTONE_MAX_ACC) && monotone_dirs(&ds, &dx, s, z, upper, regularized);

        /* (for a wide box the estimate of what the input permits is
           poor, and the corners are cheap) */
        if (!(mono && acc0 <= 24))
        {
            /* (with some more bits, for the rounding errors; a few
               bits less are accepted for input which is as accurate as
               the precision: see gamma_incomplete_core) */
            int lazy = (uprec <= 1000 && acc0 >= uprec - 16);

            method_nointegration(res, s, z, upper, regularized, prec + 30);
            acb_set_round(res, res, prec);

            if (acb_is_finite(res)
                && acb_rel_accuracy_bits(res) >= FLINT_MIN(prec, acc) - 8 - (lazy ? LAZY_SLACK : 0))
                goto finish;
        }
        else
            acb_indeterminate(res);

        if (mono)
        {
            acb_init(g);
            if (monotone_route(g, s, z, upper, regularized, ds, dx,
                    (acc0 < prec - 40) ? acc0 + 40 : prec, prec)
                && (!acb_is_finite(res)
                    || acb_rel_accuracy_bits(g) > acb_rel_accuracy_bits(res)))
                acb_set_round(res, g, prec);
            acb_clear(g);
            if (acb_is_finite(res))
                goto finish;
        }
    }

    acb_init(g);
    gamma_incomplete_core_u(g, s, z, upper, wp, uprec);

    if (!acb_is_finite(g))
    {
        /* poles in s (the regularized functions may still be finite) */
        method_nointegration(res, s, z, upper, regularized, prec);
    }
    else if (regularized == 0)
    {
        acb_swap(res, g);
    }
    else if (regularized == 1)
    {
        acb_t rg;
        acb_init(rg);
        acb_rgamma(rg, s, wp + large_bits(s));
        acb_mul(res, g, rg, prec);
        acb_clear(rg);
    }
    else  /* regularized == 2: upper Gamma(s,z) z^-s;
              lower gamma^*(s,z) = z^-s gamma(s,z) / Gamma(s) */
    {
        acb_t p;
        acb_init(p);
        acb_pow(p, z, s, wp + large_bits(s));
        acb_div(res, g, p, wp);
        if (!upper)
        {
            acb_rgamma(p, s, wp + large_bits(s));
            acb_mul(res, res, p, wp);
        }
        acb_set_round(res, res, prec);
        acb_clear(p);
    }

    acb_clear(g);

finish:
    if (real_case && acb_is_finite(res))
        arb_zero(acb_imagref(res));
}

static void
gamma_incomplete(acb_t res, const acb_t s, const acb_t z, int upper,
    int regularized, slong prec)
{
    gamma_incomplete_x(res, s, z, upper, regularized, prec, 0, prec);
}

void
acb_hypgeom_gamma_upper(acb_t res, const acb_t s, const acb_t z,
    int regularized, slong prec)
{
    gamma_incomplete(res, s, z, 1, regularized, prec);
}

void
acb_hypgeom_gamma_lower(acb_t res, const acb_t s, const acb_t z,
    int regularized, slong prec)
{
    gamma_incomplete(res, s, z, 0, regularized, prec);
}
