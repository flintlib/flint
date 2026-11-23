/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <math.h>
#include "acb.h"
#include "acb_hypgeom.h"
#include "acb_hypgeom/impl.h"

/* Enclosures for inexact input from the value at the midpoints.

   Let f be analytic on the closed polydisc with centre m = (m_i) and
   radii R_i, and |f| <= M there.  The Taylor coefficients at m satisfy
   |c_k| <= M / R^k (Cauchy), so that for |d_i| <= r_i < R_i, with
   rho_i = r_i / R_i and s = sum rho_i < 1,

       |f(m + d) - f(m)| <= M (prod 1/(1 - rho_i) - 1) <= M s / (1 - s),

       |f(m + d) - f(m) - sum_i D_i f(m) d_i| <= M s^2 / (1 - s),

   using prod (1 - rho_i) >= 1 - s.  The first bound is good to a small
   factor if M is close to |f(m)| and the R_i are about the scales on
   which f varies.  If M is a crude bound (too large by the cancellation
   of a series, say), the second bound is used: the remainder has the
   additional factor s, and the derivatives at m are bounded by divided
   differences,

       |D_i f(m)| <= |f(m + h e_i) - f(m)| / h + M h / (R_i^2 (1 - h / R_i)),

   for 0 < h < R_i (the Taylor series in the variable i alone), with h so
   small that the second term does not matter.  This costs one value of
   f per inexact variable, at a precision of the bits by which M is too
   large plus a few: the difference must be resolved, not f.  It also
   improves the first bound where M is close to |f(m)|: that bound is
   too large by M / |f(m)|, by the number of variables and by the
   factors by which the R_i are smaller than the scales of f.

   The bounds are as good with radii smaller than the scales of f by a
   factor 2^shrink (given by the caller) if M is replaced by M 2^shrink
   in the choice of h and of the precision, and in what counts as
   close. */

static void
_mag_set_d_lower(mag_t res, double x)
{
    /* x > 0, a double: exact as an arf */
    arf_t t;
    arf_init(t);
    arf_set_d(t, x);
    arf_get_mag_lower(res, t);
    arf_clear(t);
}

int
_acb_hypgeom_inexact_enclosure(acb_t res, slong * need_prec, const acb_t v, acb_srcptr args,
    slong n, const double * R, const mag_t M, slong shrink,
    void (*f)(acb_t, acb_srcptr, slong, void *), void * ctx,
    int real, slong max_prec)
{
    mag_t s, t, u, e0, e1, vlo, rad;
    slong i, K, nin, logn, p2;
    int ok, use_fd, good = 0;

    *need_prec = 0;

    if (!acb_is_finite(v) || !mag_is_finite(M))
        return 0;

    mag_init(s);
    mag_init(t);
    mag_init(u);
    mag_init(e0);
    mag_init(e1);
    mag_init(vlo);
    mag_init(rad);

    ok = 1;
    nin = 0;
    mag_zero(s);

    /* s = sum r_i / R_i */
    for (i = 0; i < n && ok; i++)
    {
        if (acb_is_exact(args + i))
            continue;

        if (!(R[i] > 0.0) || !(R[i] < 1e300))
        {
            ok = 0;
            break;
        }

        nin++;
        mag_hypot(t, arb_radref(acb_realref(args + i)), arb_radref(acb_imagref(args + i)));
        _mag_set_d_lower(u, R[i]);
        mag_div(t, t, u);
        mag_add(s, s, t);
    }

    /* (the boxes must be well inside the polydisc) */
    ok = ok && nin > 0 && mag_cmp_2exp_si(s, -2) < 0;

    acb_get_mag_lower(vlo, v);
    ok = ok && !mag_is_zero(vlo);

    if (ok)
    {
        /* M / |f(m)| < 2^K */
        mag_div(t, M, vlo);
        if (mag_cmp_2exp_si(t, 0) <= 0)
            K = 0;
        else if (mag_cmp_2exp_si(t, WORD_MAX / 8) >= 0)
            K = WORD_MAX / 8;
        else
            K = MAG_EXP(t);

        /* Radii smaller than the scales of f by 2^shrink: the bounds
           are larger by this factor than with a bound M 2^shrink for
           radii that are not. */
        K += shrink;

        /* e0 = M s / (1 - s) */
        mag_one(t);
        mag_sub_lower(u, t, s);
        mag_mul(e0, M, s);
        mag_div(e0, e0, u);

        logn = FLINT_BIT_COUNT(nin);
        /* (the differences are resolved with steps smaller by another
           factor 2^shrink) */
        p2 = K + shrink + 26 + 2 * logn;

        /* The divided differences are used if f is given (the caller
           decides whether they are worth their cost), if the value at
           the midpoint is accurate enough to resolve them, and if the
           remainder M s^2 is small then. */
        use_fd = (f != NULL) && (p2 <= max_prec)
            && (mag_cmp_2exp_si(s, -(K + 4)) < 0);

        /* (the points m + h e_i must be exact numbers of a reasonable
           size: not for midpoints with huge or tiny exponents) */
        for (i = 0; i < n && use_fd; i++)
        {
            arf_srcptr x = arb_midref(acb_realref(args + i));
            slong lim = 4 * max_prec + 2000;

            if (!acb_is_exact(args + i) && !arf_is_zero(x)
                && (arf_cmpabs_2exp_si(x, lim) >= 0 || arf_cmpabs_2exp_si(x, -lim) <= 0))
                use_fd = 0;
        }

        if (use_fd && acb_rel_accuracy_bits(v) < p2)
        {
            /* (the caller can provide a more accurate v) */
            *need_prec = p2 + 10;
            use_fd = 0;
        }

        good = (K <= 8);

        if (use_fd)
        {
            acb_ptr w;
            acb_t fw;
            arf_t h;

            w = _acb_vec_init(n);
            acb_init(fw);
            arf_init(h);

            /* e1 = M s^2 / (1 - s) + sum r_i (bound for |D_i f(m)|) */
            mag_mul(e1, M, s);
            mag_mul(e1, e1, s);
            mag_div(e1, e1, u);

            for (i = 0; i < n && mag_is_finite(e1); i++)
            {
                slong j;

                if (acb_is_exact(args + i))
                    continue;

                /* h = R_i 2^-(K + 4 + log n) */
                arf_set_d(h, R[i]);
                arf_mul_2exp_si(h, h, -(K + 4 + logn));

                for (j = 0; j < n; j++)
                    acb_get_mid(w + j, args + j);
                arb_add_arf(acb_realref(w + i), acb_realref(w + i), h, ARF_PREC_EXACT);

                f(fw, w, p2, ctx);

                if (!acb_is_finite(fw))
                {
                    mag_inf(e1);
                    break;
                }

                /* |f(m + h e_i) - f(m)| / h */
                acb_sub(fw, fw, v, p2 + 30);
                acb_get_mag(t, fw);
                arf_get_mag_lower(u, h);
                mag_div(t, t, u);

                /* M h / (R_i^2 (1 - h / R_i)), h / R_i <= 2^-4 */
                arf_get_mag(u, h);
                mag_mul(u, u, M);
                _mag_set_d_lower(rad, R[i]);
                mag_div(u, u, rad);
                mag_div(u, u, rad);
                mag_mul_2exp_si(u, u, 1);
                mag_add(t, t, u);

                mag_hypot(u, arb_radref(acb_realref(args + i)), arb_radref(acb_imagref(args + i)));
                mag_mul(t, t, u);
                mag_add(e1, e1, t);
            }

            if (mag_cmp(e1, e0) < 0)
            {
                mag_set(e0, e1);
                good = 1;
            }

            _acb_vec_clear(w, n);
            acb_clear(fw);
            arf_clear(h);
        }

        ok = mag_is_finite(e0);

        if (ok)
        {
            acb_set(res, v);
            if (real)
                arb_add_error_mag(acb_realref(res), e0);
            else
                acb_add_error_mag(res, e0);
        }
    }

    mag_clear(s);
    mag_clear(t);
    mag_clear(u);
    mag_clear(e0);
    mag_clear(e1);
    mag_clear(vlo);
    mag_clear(rad);

    return ok ? 1 + good : 0;
}

/* A bound for |pFq(a, b, z)| on balls, p <= q + 1 <= 3, from the sum of
   bounds for the absolute values of the terms (with the denominator
   k! as the last lower parameter 1):

       sum_k t_k,    t_0 = 1,    t_(k+1) = t_k |z| prod_i sup |a_i + k|
                                               / (prod_i inf |b_i + k| (k + 1)).

   This is what the series gives in ball arithmetic (too large by the
   cancellation), and close to the function where the terms have equal
   signs.  It is computed with doubles for the factors, rounded up
   (the midpoints are below 10^15; their rounding to doubles is
   covered by 4e-16 times their absolute values added to the radii,
   and the other rounding errors by the relative error 10^-10), and
   the terms in mag arithmetic: 50 to 100 ns per term.

   For j >= k, with U_i = |Re a_i| + |Im a_i| + rad(a_i) and L_i =
   Re(b_i) - rad(b_i) > -k,

       sup |a_i + j| / inf |b_i + j| <= (U_i + j) / (L_i + j)
                                     <= max(1, (U_i + k) / (L_i + k))

   (the quotient is monotone in j), and similarly with a missing a_i
   replaced by 1 / (L_i + j) <= 1 / (L_i + k): the terms decrease by a
   factor q below the product of these bounds and |z|, and the tail is
   at most t_k q / (1 - q).  The a_i can be paired with the b_i (and
   with the factor k + 1) in either order; the smaller q is used (for
   2F1(1, 1000; 2000; z), the pairs (1000, 2000) and (1, k + 1) give
   q = |z| from the first term on, and the other pairs only after
   10^5 terms).  The sum stops when q <= 1/2, or when the tail bound
   is at most 1/8 of the sum (for p = q + 1, where q stays above |z|),
   or after max_terms terms if q < 1 there.  It is infinite if a ball
   b_i contains a nonpositive integer, or if q >= 1 after max_terms
   terms. */
slong
_acb_hypgeom_inexact_terms_bound(mag_t M, acb_srcptr a, slong p, acb_srcptr b, slong q,
    const acb_t z, slong max_terms, const mag_t limit)
{
    double ar[2], ai[2], arad[2], br[2], bi[2], brad[2], zd, ratio, qk, qp;
    mag_t t, r, u, S;
    slong i, k;
    int ok, perm;

    mag_init(t);
    mag_init(r);
    mag_init(u);
    mag_init(S);

    ok = (p <= 2 && q <= 1 && p <= q + 1 && acb_is_finite(z));

    for (i = 0; i < p && ok; i++)
    {
        ok = acb_is_finite(a + i);
        ar[i] = arf_get_d(arb_midref(acb_realref(a + i)), ARF_RND_NEAR);
        ai[i] = arf_get_d(arb_midref(acb_imagref(a + i)), ARF_RND_NEAR);
        mag_hypot(r, arb_radref(acb_realref(a + i)), arb_radref(acb_imagref(a + i)));
        arad[i] = mag_get_d(r) * (1.0 + 1e-10) + 4e-16 * (fabs(ar[i]) + fabs(ai[i]) + 1.0);
        ok = ok && (hypot(ar[i], ai[i]) < 1e15) && (arad[i] < 1e15);
    }

    for (i = 0; i < q && ok; i++)
    {
        ok = acb_is_finite(b + i);
        br[i] = arf_get_d(arb_midref(acb_realref(b + i)), ARF_RND_NEAR);
        bi[i] = arf_get_d(arb_midref(acb_imagref(b + i)), ARF_RND_NEAR);
        mag_hypot(r, arb_radref(acb_realref(b + i)), arb_radref(acb_imagref(b + i)));
        brad[i] = mag_get_d(r) * (1.0 + 1e-10) + 4e-16 * (fabs(br[i]) + fabs(bi[i]) + 1.0);
        ok = ok && (hypot(br[i], bi[i]) < 1e15) && (brad[i] < 1e15);
    }

    zd = 0.0;
    if (ok)
    {
        acb_get_mag(r, z);
        zd = mag_get_d(r) * (1.0 + 1e-10);
        ok = (zd < 1e15) && (p <= q || zd < 0.999);
    }

    mag_one(t);
    mag_one(S);

    for (k = 0; ok; k++)
    {
        /* the bound q_k for the ratios from k on: a_i with b_i, or
           in the other order */
        qk = 1e300;
        for (perm = 0; perm < 2; perm++)
        {
            qp = zd;
            for (i = 0; i < q + 1 && qp < 1e300; i++)
            {
                double L = (i < q) ? br[i] - brad[i] : 1.0;
                slong j = perm ? q - i : i;

                if (!(L + k > 0.0))
                    qp = 1e300;
                else if (j < p)
                    qp *= FLINT_MAX(1.0, (fabs(ar[j]) + fabs(ai[j]) + arad[j] + k) * (1.0 + 1e-10) / (L + k));
                else
                    qp /= (L + k);
            }

            qk = FLINT_MIN(qk, qp);
        }

        qk *= (1.0 + 1e-10);

        if (qk <= 0.999999)
        {
            int stop;

            /* t_k q / (1 - q) */
            mag_set_d(r, qk / (1.0 - qk) * (1.0 + 1e-10));
            mag_mul(r, t, r);

            stop = (qk <= 0.5) || (k >= max_terms);

            if (!stop)
            {
                mag_mul_2exp_si(u, r, 3);
                stop = (mag_cmp(u, S) <= 0);
            }

            if (stop)
            {
                mag_add(S, S, r);
                break;
            }
        }

        if (k >= max_terms)
        {
            ok = 0;
            break;
        }

        /* t_(k+1) / t_k */
        ratio = zd / (k + 1.0);
        /* (the sums ar + k and br + k are rounded, k < 2^53) */
        for (i = 0; i < p; i++)
            ratio *= (hypot(ar[i] + k, ai[i]) + arad[i] + 4e-16 * k) * (1.0 + 1e-10);
        for (i = 0; i < q; i++)
        {
            double l = hypot(br[i] + k, bi[i]) * (1.0 - 1e-10) - brad[i] - 4e-16 * k;

            if (!(l > 0.0))
            {
                ok = 0;
                break;
            }
            ratio /= l;
        }

        if (!ok || !(ratio < 1e300))
        {
            ok = 0;
            break;
        }

        mag_set_d(r, ratio * (1.0 + 1e-10));
        mag_mul(t, t, r);
        mag_add(S, S, t);

        if (limit != NULL && mag_cmp(S, limit) > 0)
        {
            ok = 0;
            break;
        }
    }

    if (ok)
        mag_set(M, S);
    else
        mag_inf(M);

    mag_clear(t);
    mag_clear(r);
    mag_clear(u);
    mag_clear(S);

    return k;
}

/* The polydisc as balls: the midpoint of each inexact argument with
   radius R_i in the real and the imaginary part (a square containing
   the disc); exact arguments are copied. */
void
_acb_hypgeom_inexact_polydisc(acb_ptr balls, acb_srcptr args, slong n, const double * R)
{
    slong i;
    mag_t r;

    mag_init(r);

    for (i = 0; i < n; i++)
    {
        acb_get_mid(balls + i, args + i);

        if (!acb_is_exact(args + i))
        {
            mag_set_d(r, R[i]);
            arb_add_error_mag(acb_realref(balls + i), r);
            arb_add_error_mag(acb_imagref(balls + i), r);
        }
    }

    mag_clear(r);
}

/* A bound on a smaller polydisc.  The integrals are not computed for
   wide balls, and the series on balls can fail to converge or be far
   too large for small parameters: the smallest one of the bounds
   majorant(stage) for stage0 <= stage <= stage1 is tried with radii
   smaller by 2^step and 2^(2 step) (and with the given ones first if
   first is set).  The enclosure is as tight with such radii
   if divided differences are used; shrink counts the bits by which the
   radii are smaller than the scales of f.  The radii and M are
   replaced if this improves M 2^shrink (the bound for the previous
   radii remains valid otherwise); returns whether they are. */
static int
_inexact_smaller(mag_t M, double * R, slong * shrink, acb_ptr balls, acb_srcptr args, slong n,
    void (*majorant)(mag_t, acb_srcptr, int, void *), void * ctx, int stage0, int stage1,
    int step, int first)
{
    mag_t M2, M3;
    double R2[ACB_HYPGEOM_INEXACT_MAX_ARGS];
    slong i;
    int k, j, improved = 0;

    mag_init(M2);
    mag_init(M3);

    for (i = 0; i < n; i++)
        R2[i] = R[i];

    for (k = first ? 0 : 1; k < 3; k++)
    {
        if (k > 0)
        {
            for (i = 0; i < n; i++)
                R2[i] = ldexp(R2[i], -step);
            _acb_hypgeom_inexact_polydisc(balls, args, n, R2);
        }

        mag_inf(M2);
        for (j = stage0; j <= stage1; j++)
        {
            majorant(M3, balls, j, ctx);
            mag_min(M2, M2, M3);
        }

        if (mag_is_finite(M2))
        {
            mag_mul_2exp_si(M3, M2, step * k);
            if (mag_cmp(M3, M) < 0)
            {
                mag_set(M, M2);
                *shrink += step * k;
                for (i = 0; i < n; i++)
                    R[i] = R2[i];
                improved = 1;
            }
            break;
        }
    }

    _acb_hypgeom_inexact_polydisc(balls, args, n, R);

    mag_clear(M2);
    mag_clear(M3);
    return improved;
}

/* The function on inexact input: its value at the midpoints, and the
   bound for the variation with a bound for |f| on the polydisc.

   - majorant(M, balls, 0, ctx): a bound that is close to |f| where it
     applies (from a series with positive terms, or from an integral
     with a positive integrand).  The radii are meant to be lower
     bounds for the scales on which f varies, and the enclosure is
     tighter with large radii as long as the bound stays close to the
     value: it is tried with four times the radii first if the bound
     costs little compared to the value (stage0_cheap, or if it is a
     value of f at low precision and that is cheap), and with twice
     the radii otherwise.
   - majorant(M, balls, k, ctx) for 1 <= k <= ncrude: bounds that cost
     little, typically enclosures of f in ball arithmetic which are too
     large by the cancellation in a series.
   - majorant(M, balls, -1, ctx): a bound that is close to |f| but
     expensive (numerical integration at low precision); exp_ratio is
     an estimate of its cost divided by that of the value of f.

   A majorant that does not apply sets M to infinity.  The cost of a
   value of f at a low precision p, divided by that of the value at
   the precision wp, is estimated as the ratio of the costs of a term
   of a series at the precisions p + c and wp + c, where c is the
   larger one of cancel (an estimate of the bits that the series for f
   lose to cancellation) and of the bits by which the bound from the
   series is too large.  (The number of terms depends little on the
   precision: a value at 40 bits costs a quarter to a half of one at
   a few hundred bits, and a tenth of one at a few thousand.)

   With a bound M within 2^8 of |f(m)|, the first-order bound M s is
   used, with divided differences if they add at most half to the
   cost.  Otherwise, with M = 2^K |f(m)|:

   - Divided differences give an enclosure close to the range of f if
     the input is accurate to more than K bits, at the cost of one
     value of f per inexact variable at a precision of about K bits.
   - The expensive bound gives one at the cost of the integrals.
   - The first-order bound M s loses K bits, like the series in ball
     arithmetic.

   The cheaper one of the first two is used.  If lazy is set (the input
   is as accurate as the precision: the caller can increase both), only
   if it costs no more than twice the value of f itself (less than the
   function at twice the precision): otherwise the result is the
   first-order bound, rather than several milliseconds are spent on
   integrals where the series take microseconds.  A caller who
   increases the precision gets an accurate result from that bound, or
   from divided differences when the precision is large compared to
   K.

   Returns 0 if there is no result, 1 if the result is a valid
   enclosure which may be far from the range of f, and 2 if it is
   close. */
static double
_inexact_low_cost(double p, double wp, double c, int real)
{
    return _acb_hypgeom_est_series_term_cost(p + c, p + c, real)
        / _acb_hypgeom_est_series_term_cost(wp + c, wp + c, real);
}

int
_acb_hypgeom_inexact(acb_t res, acb_srcptr args, slong n, const double * R0,
    void (*f)(acb_t, acb_srcptr, slong, void *),
    void (*majorant)(mag_t, acb_srcptr, int, void *), int ncrude,
    void * ctx, int real, slong prec, slong goal, double cancel, double exp_ratio, int lazy,
    int stage0_cheap)
{
    acb_ptr balls, mids;
    acb_t v, res2;
    mag_t M, M2, vlo, s, t;
    double R[ACB_HYPGEOM_INEXACT_MAX_ARGS], fd_cost, limit, c;
    slong i, wp, need, max_prec, shrink = 0, K, p2, nin, logn;
    int status, k, feasible;

    if (n > ACB_HYPGEOM_INEXACT_MAX_ARGS)
        return 0;

    balls = _acb_vec_init(n);
    mids = _acb_vec_init(n);
    acb_init(v);
    acb_init(res2);
    mag_init(M);
    mag_init(M2);
    mag_init(vlo);
    mag_init(s);
    mag_init(t);

    nin = 0;
    for (i = 0; i < n; i++)
    {
        acb_get_mid(mids + i, args + i);
        nin += !acb_is_exact(args + i);
    }
    logn = FLINT_BIT_COUNT(nin);

    wp = FLINT_MIN(prec, goal + 30);
    max_prec = goal + 30;
    status = 0;
    cancel = FLINT_MIN(FLINT_MAX(cancel, 0.0), 1e12);

    /* Boxes that are too large for the largest polydisc: nothing is
       computed. */
    k = 1;
    mag_zero(s);
    for (i = 0; i < n; i++)
    {
        if (!acb_is_exact(args + i))
        {
            if (!(R0[i] > 0.0) || !(R0[i] < 1e300))
            {
                k = 0;
                break;
            }

            mag_hypot(t, arb_radref(acb_realref(args + i)), arb_radref(acb_imagref(args + i)));
            _mag_set_d_lower(M2, 4.0 * R0[i]);
            mag_div(t, t, M2);
            mag_add(s, s, t);
        }
    }

    if (k == 0 || nin == 0 || mag_cmp_2exp_si(s, -2) >= 0)
    {
        _acb_vec_clear(balls, n);
        _acb_vec_clear(mids, n);
        acb_clear(v);
        acb_clear(res2);
        mag_clear(M);
        mag_clear(M2);
        mag_clear(vlo);
        mag_clear(s);
        mag_clear(t);
        return 0;
    }

    f(v, mids, wp, ctx);
    acb_get_mag_lower(vlo, v);

    if (acb_is_finite(v) && !mag_is_zero(vlo))
    {
        if (stage0_cheap || _inexact_low_cost(40.0, wp, cancel, real) <= 0.25)
        {
            for (i = 0; i < n; i++)
                R[i] = 4.0 * R0[i];
            _acb_hypgeom_inexact_polydisc(balls, args, n, R);
            majorant(M, balls, 0, ctx);

            mag_mul_2exp_si(M2, vlo, 2);
            if (mag_cmp(M, M2) > 0)
            {
                for (i = 0; i < n; i++)
                    R[i] = R0[i];
                _acb_hypgeom_inexact_polydisc(balls, args, n, R);
                majorant(M, balls, 0, ctx);
            }
        }
        else
        {
            for (i = 0; i < n; i++)
                R[i] = 2.0 * R0[i];
            _acb_hypgeom_inexact_polydisc(balls, args, n, R);
            majorant(M, balls, 0, ctx);
        }

        for (k = 1; k <= ncrude; k++)
        {
            mag_mul_2exp_si(M2, vlo, 8);
            if (mag_cmp(M, M2) <= 0)
                break;

            majorant(M2, balls, k, ctx);
            mag_min(M, M, M2);
        }

        /* no bound, or one that is too large by more than 2^24: the
           cheap bounds on a smaller polydisc */
        mag_mul_2exp_si(M2, vlo, 24);
        if (ncrude > 0 && mag_cmp(M, M2) > 0)
            _inexact_smaller(M, R, &shrink, balls, args, n, majorant, ctx, 1, ncrude, 4, 0);

        limit = lazy ? 2.0 : 1e300;
        fd_cost = 1e300;

        /* no bound at all: the expensive one */
        if (!mag_is_finite(M) && exp_ratio <= limit)
            _inexact_smaller(M, R, &shrink, balls, args, n, majorant, ctx, -1, -1, 8, 1);

        /* the first-order bound */
        status = _acb_hypgeom_inexact_enclosure(res, &need, v, args, n, R, M, shrink,
                    NULL, NULL, real, max_prec);

        for (k = 0; k < 2 && status != 0; k++)
        {
            /* M 2^shrink < 2^K |f(m)|, s = sum r_i / R_i */
            mag_div(t, M, vlo);
            K = (mag_cmp_2exp_si(t, 0) <= 0) ? 0
                : (mag_cmp_2exp_si(t, WORD_MAX / 8) >= 0) ? WORD_MAX / 8 : MAG_EXP(t);
            K += shrink;

            mag_zero(s);
            for (i = 0; i < n; i++)
            {
                if (!acb_is_exact(args + i))
                {
                    mag_hypot(t, arb_radref(acb_realref(args + i)), arb_radref(acb_imagref(args + i)));
                    _mag_set_d_lower(M2, R[i]);
                    mag_div(t, t, M2);
                    mag_add(s, s, t);
                }
            }

            /* the divided differences: nin values of f at the
               precision p2 */
            p2 = K + shrink + 26 + 2 * logn;
            feasible = (K < WORD_MAX / 16) && (p2 <= max_prec)
                && (mag_cmp_2exp_si(s, -(K + 4)) < 0);
            c = FLINT_MAX(cancel, (double) K);
            fd_cost = feasible ? nin * _inexact_low_cost(p2, wp, c, real) : 1e300;

            /* the expensive bound, if the bound is not close, and if
               the divided differences are not cheaper */
            if (k == 1 || status == 2 || fd_cost <= 1.0
                || !(exp_ratio <= limit) || !(exp_ratio < fd_cost)
                || !_inexact_smaller(M, R, &shrink, balls, args, n, majorant, ctx, -1, -1, 8, 1))
                break;

            status = _acb_hypgeom_inexact_enclosure(res, &need, v, args, n, R, M, shrink,
                        NULL, NULL, real, max_prec);
        }

        /* the divided differences: with a bound that is close, if
           they add little to the cost */
        if (status != 0 && fd_cost <= ((status == 2) ? 0.5 : limit))
        {
            int status2;

            status2 = _acb_hypgeom_inexact_enclosure(res2, &need, v, args, n, R, M, shrink,
                        f, ctx, real, max_prec);

            if (need != 0 && status != 2)
            {
                f(v, mids, FLINT_MAX(wp, need), ctx);
                status2 = _acb_hypgeom_inexact_enclosure(res2, &need, v, args, n, R, M, shrink,
                            f, ctx, real, max_prec);
            }

            if (status2 != 0)
            {
                acb_swap(res, res2);
                status = FLINT_MAX(status, status2);
            }
        }
    }

    _acb_vec_clear(balls, n);
    _acb_vec_clear(mids, n);
    acb_clear(v);
    acb_clear(res2);
    mag_clear(M);
    mag_clear(M2);
    mag_clear(vlo);
    mag_clear(s);
    mag_clear(t);

    return status;
}

int
_acb_hypgeom_inexact_is_coarse(const acb_t x, slong prec)
{
    return !acb_is_exact(x) && acb_rel_accuracy_bits(x) < prec - 16;
}

/* Whether remainders which continue to decrease as from r0 to r1 fall
   below the target within the given number of further orders (with a
   margin of a factor 2 per order) */
static int
_taylor_in_reach(const mag_t r0, const mag_t r1, const mag_t target, slong orders)
{
    double l0, l1, lt, slope;

    if (mag_is_zero(r1))
        return 1;
    if (mag_is_inf(r0) || mag_is_inf(r1) || mag_is_zero(target))
        return 0;

    l0 = mag_get_d_log2_approx(r0);
    l1 = mag_get_d_log2_approx(r1);
    lt = mag_get_d_log2_approx(target);
    slope = l1 - l0;

    if (slope > -2.0)
        return 0;

    return l1 + (slope - 1.0) * orders <= lt;
}

int
_acb_hypgeom_inexact_taylor(acb_t res, slong * order, const acb_t z,
    void (*coeff)(acb_t, slong, const acb_t, slong, int, void *),
    void (*recur)(acb_t, const acb_t, const acb_t, slong, const acb_t, slong, void *),
    void * ctx, int real, slong prec, slong goal, slong max_order, slong max_values)
{
    acb_t m, v, u0, u1, t;
    mag_t r, rj, target, vlo, u, err;
    mag_ptr rem;
    slong J, Jg, j, wp, have, values;
    int status = 0;

    *order = 0;

    if (!acb_is_finite(z) || max_order < 1)
        return 0;

    acb_init(m);
    acb_init(v);
    acb_init(u0);
    acb_init(u1);
    acb_init(t);
    mag_init(r);
    mag_init(rj);
    mag_init(target);
    mag_init(vlo);
    mag_init(u);
    mag_init(err);
    rem = _mag_vec_init(max_order + 2);

    acb_get_mid(m, z);
    mag_hypot(r, arb_radref(acb_realref(z)), arb_radref(acb_imagref(z)));

    /* The remainders rem[J] = sup |c_J| r^J with the crude bounds, as
       far as they are needed (have is the last one) */
#define TAYLOR_REMAINDER(JJ) \
    do { \
        while (have < (JJ) && have >= 0) \
        { \
            have++; \
            coeff(t, have, z, 64, 1, ctx); \
            if (!acb_is_finite(t)) \
            { \
                have = -1; \
                break; \
            } \
            acb_get_mag(u, t); \
            mag_mul(rem + have, u, rj); \
            mag_mul(rj, rj, r); \
        } \
    } while (0)

    have = 0;
    mag_set(rj, r);

    TAYLOR_REMAINDER(1);
    if (have < 0)
        goto cleanup;

    /* Before anything expensive is computed: the order which would be
       needed if f(m) and the first-order term were as large as their
       crude bounds (a lower bound for the order) */
    coeff(t, 0, m, 64, 1, ctx);
    if (acb_is_finite(t))
    {
        acb_get_mag(u, t);
        mag_mul_2exp_si(target, u, -(goal + 4));
        mag_mul_2exp_si(u, rem + 1, -3);
        mag_max(target, target, u);

        for (Jg = 1; Jg <= max_order; Jg++)
        {
            TAYLOR_REMAINDER(Jg);
            if (have < 0)
                goto cleanup;
            if (mag_cmp(rem + Jg, target) <= 0)
                break;

            /* (the remainders decrease by a factor which changes
               little with the order: stop if the largest order will
               not do, or if they do not decrease) */
            if (Jg >= 2 && !_taylor_in_reach(rem + Jg - 1, rem + Jg, target, max_order - Jg))
            {
                Jg = max_order;
            }
        }

        if (Jg > max_order)
        {
            *order = max_order + 1;
            goto cleanup;
        }
    }
    else
    {
        Jg = 1;
    }

    /* the value at the midpoint */
    wp = FLINT_MIN(prec, goal + 30) + 10;
    coeff(v, 0, m, wp, 0, ctx);
    acb_get_mag_lower(vlo, v);

    if (!acb_is_finite(v) || mag_is_zero(vlo))
        goto cleanup;

    mag_mul_2exp_si(target, vlo, -(goal + 2));
    J = 1;
    mag_zero(err);

    /* (the first order with the crude bound for f' may do) */
    if (mag_cmp(rem + 1, target) > 0)
    {
        if (max_order < 2)
        {
            *order = max_order + 1;
            goto cleanup;
        }

        /* the first-order term, which is the variation of f over the
           box if the remainder is small */
        coeff(u1, 1, m, (recur != NULL) ? wp : 30, 0, ctx);
        if (!acb_is_finite(u1))
            goto cleanup;

        acb_get_mag(u, u1);
        mag_mul(err, u, r);

        mag_mul_2exp_si(target, vlo, -(goal + 4));
        mag_mul_2exp_si(u, err, -3);
        mag_max(target, target, u);

        for (J = FLINT_MAX(Jg, 2); J <= max_order; J++)
        {
            TAYLOR_REMAINDER(J);
            if (have < 0)
                goto cleanup;
            if (mag_cmp(rem + J, target) <= 0)
                break;
            if (!_taylor_in_reach(rem + J - 1, rem + J, target, max_order - J))
                J = max_order;
        }

        if (J > max_order)
        {
            *order = max_order + 1;
            goto cleanup;
        }

        /* The terms of the orders 2, ..., J - 1: from the recurrence
           which the differential equation gives for the coefficients,
           in ball arithmetic, and from values of the derivatives where
           there is none or it has lost its accuracy. */
        acb_set(u0, v);
        mag_mul(rj, r, r);
        values = 2;

        for (j = 2; j < J; j++)
        {
            if (recur != NULL)
                recur(t, u0, u1, j - 2, m, wp, ctx);

            if (recur == NULL || !acb_is_finite(t)
                || (!acb_is_zero(t) && acb_rel_accuracy_bits(t) < 4))
            {
                if (values >= max_values)
                {
                    *order = max_order + 1;
                    goto cleanup;
                }

                values++;
                coeff(t, j, m, (recur != NULL) ? wp : 30, 0, ctx);
                if (!acb_is_finite(t))
                    goto cleanup;
            }

            acb_get_mag(u, t);
            mag_addmul(err, u, rj);
            mag_mul(rj, rj, r);

            acb_swap(u0, u1);
            acb_swap(u1, t);
        }
    }

    *order = J;
    mag_add(err, err, rem + J);

    if (real)
        arb_add_error_mag(acb_realref(v), err);
    else
        acb_add_error_mag(v, err);

    acb_set_round(res, v, prec);
    status = 1;

cleanup:
#undef TAYLOR_REMAINDER
    acb_clear(m);
    acb_clear(v);
    acb_clear(u0);
    acb_clear(u1);
    acb_clear(t);
    mag_clear(r);
    mag_clear(rj);
    mag_clear(target);
    mag_clear(vlo);
    mag_clear(u);
    mag_clear(err);
    _mag_vec_clear(rem, max_order + 2);

    return status;
}
