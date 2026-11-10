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
   the methods are tried in order of cheapness at a doubling working
   precision, each result verified a posteriori and the best kept; the
   integration, which evaluates the function on the input balls
   directly, is tried as soon as its estimated cost is below that of
   the series and expansions; and for inexact inputs there is also the
   evaluation at the midpoints with the propagated error added from
   derivative bounds (the z-derivative -+z^(s-1) e^(-z) in closed form,
   which makes this the first choice when only z is inexact; the
   s-derivative bounded along the steepest-descent contour). */

#include <math.h>
#include "arb.h"
#include "acb.h"
#include "arb_hypgeom.h"
#include "acb_hypgeom.h"
#include "acb_hypgeom/impl.h"

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
            t -= (slong) FLINT_MIN(log2(1.0 + zn + sn), 1e6);
        acc = FLINT_MIN(acc, t);
    }

    if (!acb_is_exact(s))
    {
        t = upper ? acb_rel_one_accuracy_bits(s) : acb_rel_accuracy_bits(s);
        if (t > 0)
            t -= (slong) FLINT_MIN(log2(1.0 + sn * (1.0 + fabs(
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

/* the series-class route */
static void
method_series(acb_t res, const acb_t s, const acb_t z, int upper,
    slong wp)
{
    /* the _nointegration functions bundle the convergent series with
       their elementary and singular special cases; any cancellation is
       absorbed by the caller's doubling working precision */
    if (upper)
        acb_hypgeom_gamma_upper_nointegration(res, s, z, 0, wp);
    else
    {
        acb_t s1, w, t;
        double zre, zim, zn, sn, cancel;
        slong wp2;
        int kummer;

        /* gamma(s, z) = s^-1 z^s 1F1(s, s+1, -z)
                       = s^-1 z^s e^-z 1F1(1, s+1, z)
           by the series itself (acb_hypgeom_m would switch to the
           asymptotic expansion for large |z| regardless of s, which
           fails and is expensive for |s| comparable to |z|; the
           expansion is the other route).  The terms of the first
           series grow to e^|z| and its sum is of size e^(-Re z); the
           terms z^k / (s+1)_k of the second do not grow if |z| < |s|,
           and otherwise grow to at most e^|z| with a sum of size
           e^(Re z).  The form with the smaller cancellation is used,
           with the cancellation absorbed by the working precision. */
        zre = arf_get_d(arb_midref(acb_realref(z)), ARF_RND_NEAR);
        zim = arf_get_d(arb_midref(acb_imagref(z)), ARF_RND_NEAR);
        zn = hypot(zre, zim);
        sn = hypot(arf_get_d(arb_midref(acb_realref(s)), ARF_RND_NEAR),
                   arf_get_d(arb_midref(acb_imagref(s)), ARF_RND_NEAR));
        kummer = (zre > 0.0 || zn < sn);
        cancel = 1.45 * (kummer ? ((zn < sn) ? 0.0 : zn - zre) : zn + zre);

        wp2 = wp;
        if (zn < 1e15)
            wp2 += (slong) FLINT_MIN(cancel, 4.0 * (double) wp);
        /* (no more than the input accuracy allows, as in acb_hypgeom_m) */
        if (!acb_is_exact(s) || !acb_is_exact(z))
        {
            slong acc = FLINT_MIN(acb_rel_one_accuracy_bits(s), acb_rel_one_accuracy_bits(z));
            wp2 = FLINT_MAX(2, FLINT_MIN(wp2, FLINT_MAX(acc, 0) + 20));
        }

        acb_init(s1);
        acb_init(w);
        acb_init(t);
        acb_add_ui(s1, s, 1, wp2);

        if (kummer)
        {
            /* sum of z^k / (s+1)_k */
            acb_hypgeom_pfq_direct(t, NULL, 0, s1, 1, z, -1, wp2);
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
            acb_hypgeom_pfq_direct(t, s, 1, bb, 2, w, -1, wp2);
            acb_clear(bb + 1);
        }

        acb_pow(res, z, s, wp);
        acb_mul(t, t, res, wp);
        acb_div(res, t, s, wp);
        acb_clear(s1);
        acb_clear(w);
        acb_clear(t);

        /* the lower series can fail where the upper series
           converges; reflect through Gamma(s), the cancellation again
           being absorbed by the doubling working precision */
        if (!acb_is_finite(res))
        {
            acb_t g;
            acb_init(g);
            acb_hypgeom_gamma_upper_nointegration(res, s, z, 0, wp);
            acb_gamma(g, s, wp);
            acb_sub(res, g, res, wp);
            acb_clear(g);
        }
    }
}

/* the asymptotic route: the (self-limiting, divergent) 2F0 expansion
   of the upper gamma; the lower gamma by reflection, whose
   cancellation at small |z| is absorbed by the doubling working
   precision and is complementary to the series route */
static void
method_asymp(acb_t res, const acb_t s, const acb_t z, int upper,
    slong wp)
{
    acb_hypgeom_gamma_upper_asymp(res, s, z, 0, wp);
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

    if (arb_contains_zero(acb_imagref(z)) && !arb_is_exact(acb_imagref(z))
        && !arb_is_positive(acb_realref(z))
        && !(acb_is_int(s) && arb_is_positive(acb_realref(s))))
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

/* Rough cost estimates (in microseconds on a typical machine; only the
   ratios matter) of the series at the working precision wp, of the
   asymptotic expansion (see asymp_terms), and of the integration at
   precision prec (see bessel.c).  The terms z^k / (s+1)_k of the
   series decrease from k ~ |z| - Re(s), by a factor |z| / |s + k|. */
static double
series_terms(double sre, double zn, slong wp)
{
    double k0 = fmax(2.0 * zn - fmax(sre, 0.0), 0.0);
    return k0 + (double) wp / fmax(1.0, log2((fmax(sre, 0.0) + k0 + 1.0) / (zn + 1e-300))) + 5.0;
}

static double
series_cost(double sre, double zn, slong wp)
{
    return 5.0 + series_terms(sre, zn, wp) * (0.12 + 0.0005 * (double) wp);
}

/* The number of terms (1-s)_k / z^k of the asymptotic expansion needed
   for wp bits: the k-th ratio is about hypot(|s|, k) / |z|, so that n
   terms gain about n log(|z| / hypot(|s|, n/e)) (n! / z^n for small
   |s|); infinite if the terms start growing before that. */
static double
asymp_terms(double num, double zn, slong wp)
{
    double n = 0.0, L;
    int i;

    for (i = 0; i < 4; i++)
    {
        L = log(zn / hypot(num, 0.368 * n + 1.0));
        if (!(L > 0.01))
            return INFINITY;
        n = 0.693 * (double) wp / L;
    }

    return (hypot(num, n) < zn) ? n : INFINITY;
}

static double
asymp_cost(double num, double zn, slong wp)
{
    return 20.0 + asymp_terms(num, zn, wp) * (0.3 + 0.001 * (double) wp);
}

static double
integration_cost(double num, double zn, slong prec)
{
    static const double dcost[4] = { 30.0, 45.0, 100.0, 200.0 };
    double bits = (double) prec + 12.0 + log2(2.0 + zn + num);

    if (bits <= 4 * 53)
        return 2.0 * dcost[(int) ((bits - 1.0) / 53.0)];
    return 2.0 * 600.0 * (bits / 256.0) * (bits / 256.0);
}

/* The cancellation (in bits) which the series must absorb through its
   working precision: that of the lower series in the form used (see
   method_series), and for the upper function with Re(z) > 0 also that
   of Gamma(s) - gamma(s,z), where |Gamma(s,z)| ~ |z|^(Re(s)-1) e^(-Re z)
   for |z| beyond Re(s). */
static double
series_cancellation(int upper, double sre, double sn, double zn, double zre)
{
    double c;

    if (!upper)
    {
        if (zre > 0.0 || zn < sn)
            return (zn < sn) ? 0.0 : 1.45 * (zn - zre);
        return 1.45 * (zn + zre);
    }

    c = 1.45 * (zn - fabs(zre));
    if (zre > 0.0 && zn > sre)
        c += fmax(0.0, 1.45 * (lgamma(fmax(sre, 1.0)) - (sre - 1.0) * log(zn) + zre));
    return c;
}

static void
gamma_incomplete_core(acb_t res, const acb_t s, const acb_t z, int upper,
    slong prec)
{
    double num, zn, sre, c_int, c_first, cancel;
    int asymp_ok;
    slong wp, acc_in, goal, target, accept, best_acc = -WORD_MAX;
    acb_t best, t;
    int inexact, round, series_plausible, mid_tried = 0, int_ok, int_tried = 0, int_done = 0, int_first;

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

    /* the goal is capped by the relative accuracy of the inputs:
       nothing more is attainable, and hopeless input balls then
       terminate after the first cheap attempt; the target also
       accounts (heuristically) for the conditioning */
    /* (Gamma(s,z) is entire in s, while gamma(s,z) has a pole at 0) */
    acc_in = FLINT_MIN(upper ? acb_rel_one_accuracy_bits(s) : acb_rel_accuracy_bits(s),
                       acb_rel_accuracy_bits(z));
    acc_in = FLINT_MAX(acc_in, 0);
    goal = FLINT_MIN(prec, acc_in);
    target = FLINT_MIN(goal, attainable_accuracy(s, z, upper, num, zn, prec));

    /* Accept results this close to the target.  For inexact input the
       target is only a heuristic estimate of what is attainable; a
       result short of it is compared with the integration on the
       input balls. */
    accept = target - 5;

    /* nothing is attainable: direct interval propagation is as good as
       anything, quickly */
    if (goal <= 0 || target < 4)
    {
        method_nointegration(res, s, z, upper, 0, prec);
        return;
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
    int_ok = arf_cmpabs_2exp_si(arb_midref(acb_realref(z)), 44) < 0
        && arf_cmpabs_2exp_si(arb_midref(acb_imagref(z)), 44) < 0
        && arf_cmpabs_2exp_si(arb_midref(acb_realref(s)), 44) < 0
        && arf_cmpabs_2exp_si(arb_midref(acb_imagref(s)), 44) < 0
        && zn > 1e-270 && !acb_contains_zero(z);
    c_int = integration_cost(num, zn, (acc_in < prec - 30) ? acc_in + 30 : prec);
    asymp_ok = isfinite(asymp_terms(num, zn, goal + 8));
    cancel = series_cancellation(upper, sre, num, zn,
                arf_get_d(arb_midref(acb_realref(z)), ARF_RND_NEAR));
    cancel = FLINT_MIN(cancel, 4.0 * (double) goal);
    if (asymp_ok)
        c_first = asymp_cost(num, zn, goal + 8);
    else
        c_first = 1.5 * series_cost(sre, zn, goal + 30 + (slong) cancel);
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
    else if (!asymp_ok)
        wp = goal + (slong) cancel;

    for (wp = wp + 30, round = 0; wp <= 16 * goal + 480;
            wp = 2 * wp + 30, round++)
    {
        int m, series_first;
        slong iter_start_acc = best_acc;

        series_first = (series_terms(sre, zn, wp) <= asymp_terms(num, zn, wp));

        /* First, unless the integration is expected to be cheaper, the
           series-and-asymptotics implementation, which is tuned for
           the easy cases. */
        if (round == 0 && !int_first)
        {
            /* (with an inexact z, the cancellation of the series
               inflates the radius regardless of the precision) */
            if (!(inexact && acb_is_exact(s) && cancel >= 4.0 && !asymp_ok))
            {
                /* (directly where the expansion converges: the
                   thresholds of the implementations differ) */
                if (asymp_ok)
                    method_asymp(t, s, z, upper, inexact ? wp : prec + 10);
                else
                    method_nointegration(t, s, z, upper, 0, inexact ? wp
                        : prec + 10 + (slong) cancel);
                keep_best(best, &best_acc, t);
                /* (if only z is inexact, the evaluation at the midpoint
                   is cheap and nearly optimal: less slack) */
                if (best_acc >= ((inexact && acb_is_exact(s)) ? target - 2 : accept))
                    goto done;
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

        /* the integration, first when expected to be cheapest overall,
           or once it is cheaper than the next round of the series and
           expansions */
        if (int_ok && !int_tried
            && (int_first || c_int <= (series_first ? series_cost(sre, zn, wp) : asymp_cost(num, zn, wp))))
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
                /* no term-count gate: the bundled implementations
                   select their own methods and abort cheaply when
                   nothing converges; the doubling then finds the
                   cancellation-compensated precision */
                method_series(t, s, z, upper, wp);
                keep_best(best, &best_acc, t);
            }
            else
            {
                method_asymp(t, s, z, upper, wp);
                keep_best(best, &best_acc, t);
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
            if (int_ok && !int_tried && c_int <= 4.0 * (c_first + 50.0))
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
            if (best_acc >= 2)
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
            && (best_acc > 0 || !series_plausible))
            break;
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

static void
gamma_incomplete(acb_t res, const acb_t s, const acb_t z, int upper,
    int regularized, slong prec)
{
    acb_t g;
    slong wp;
    int real_case;

    if (res == s || res == z)
    {
        acb_init(g);
        gamma_incomplete(g, s, z, upper, regularized, prec);
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

    acb_init(g);
    gamma_incomplete_core(g, s, z, upper, wp);

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
        acb_rgamma(rg, s, wp);
        acb_mul(res, g, rg, prec);
        acb_clear(rg);
    }
    else  /* regularized == 2: upper Gamma(s,z) z^-s;
              lower gamma^*(s,z) = z^-s gamma(s,z) / Gamma(s) */
    {
        acb_t p;
        acb_init(p);
        acb_pow(p, z, s, wp);
        acb_div(res, g, p, wp);
        if (!upper)
        {
            acb_rgamma(p, s, wp);
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
