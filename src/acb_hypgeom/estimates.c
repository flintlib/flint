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

/* Estimates in double precision for the choice of algorithm in the
   hypergeometric functions 1F1, U and 2F1: sizes of terms, numbers of terms,
   cancellation and costs.  They decide what is tried and with how much
   working precision; the accuracy of each result is checked afterwards,
   so that nothing here needs to be rigorous. */

#define M_INVLOG2 1.4426950408889634
#define M_NEGINF (-0x1p996)  /* (exact also in extended precision, as on x87: -1e300 is not) */

typedef _acb_hypgeom_cd m_cd;

m_cd
_acb_hypgeom_cd_mul(m_cd x, m_cd y)
{
    m_cd r;
    r.re = x.re * y.re - x.im * y.im;
    r.im = x.re * y.im + x.im * y.re;
    return r;
}

m_cd
_acb_hypgeom_cd_div(m_cd x, m_cd y)
{
    m_cd r;
    double d = y.re * y.re + y.im * y.im;
    r.re = (x.re * y.re + x.im * y.im) / d;
    r.im = (x.im * y.re - x.re * y.im) / d;
    return r;
}

m_cd
_acb_hypgeom_cd_sqrt(m_cd x)
{
    m_cd r;
    double m = hypot(x.re, x.im);
    r.re = sqrt(0.5 * (m + fabs(x.re)));
    r.im = (r.re == 0.0) ? 0.0 : 0.5 * fabs(x.im) / r.re;
    if (x.re < 0.0)
    {
        double t = r.re;
        r.re = r.im;
        r.im = t;
    }
    if (x.im < 0.0)
        r.im = -r.im;
    return r;
}

m_cd
_acb_hypgeom_cd_log(m_cd x)
{
    m_cd r;
    r.re = log(hypot(x.re, x.im));
    r.im = atan2(x.im, x.re);
    return r;
}

double
_acb_hypgeom_est_lgamma_re(m_cd x)
{
    return _acb_hypgeom_d_lgamma_re(x.re, x.im, 0);
}

/* log |(a)_k|, with -infinity for a vanishing product */
double
_acb_hypgeom_est_lrise(m_cd a, double k)
{
    m_cd u, v;

    if (k <= 0.0)
        return 0.0;

    if (a.im == 0.0 && a.re <= 0.0 && a.re == floor(a.re) && k > -a.re)
        return M_NEGINF;

    if (a.re >= 0.5)
    {
        u.re = a.re + k; u.im = a.im;
        return _acb_hypgeom_est_lgamma_re(u) - _acb_hypgeom_est_lgamma_re(a);
    }

    /* (a)_j = (-1)^j Gamma(1 - a) / Gamma(1 - a - j) as far as to the
       left of 1/2 */
    u.re = 1.0 - a.re; u.im = -a.im;

    if (a.re + k <= 0.5)
    {
        v.re = u.re - k; v.im = u.im;
        return _acb_hypgeom_est_lgamma_re(u) - _acb_hypgeom_est_lgamma_re(v);
    }
    else
    {
        double j = ceil(0.5 - a.re), s;
        v.re = u.re - j; v.im = u.im;
        s = _acb_hypgeom_est_lgamma_re(u) - _acb_hypgeom_est_lgamma_re(v);
        u.re = a.re + j; u.im = a.im;
        v.re = a.re + k; v.im = a.im;
        return s + _acb_hypgeom_est_lgamma_re(v) - _acb_hypgeom_est_lgamma_re(u);
    }
}

/* log of the absolute value of term k of 1F1(a, b, z) */
static double
m_logterm(m_cd a, m_cd b, double lz, double k)
{
    double u = _acb_hypgeom_est_lrise(a, k), w;

    if (u <= M_NEGINF)
        return M_NEGINF;

    /* (Beyond a pole of the lower parameter, b = -m: without the
       factor zero, m! (k - m - 1)!.  This is the size of the terms of
       the regularized series times Gamma(m + 1), which is what the
       connection formula of U has in the limit of an integer b.) */
    w = _acb_hypgeom_est_lrise(b, k);
    if (w <= M_NEGINF)
        w = _acb_hypgeom_d_lgamma_re(1.0 - b.re, 0.0, 0) + _acb_hypgeom_d_lgamma_re(k + b.re, 0.0, 0);

    return u - w + k * lz - _acb_hypgeom_d_lgamma_re(k + 1.0, 0.0, 0);
}

/* the last index of a terminating series, or a huge number */
static double
m_last_index(m_cd a)
{
    if (a.im == 0.0 && a.re <= 0.0 && a.re == floor(a.re))
        return -a.re;
    return 1e300;
}

/* The largest term of 1F1(a, b, z): its logarithm, and its index in kpk. */
double
_acb_hypgeom_est_1f1_peak(double * kpk, m_cd a, m_cd b, m_cd z)
{
    double lz, best = 0.0, kb = 0.0, k, kprev, lo, hi, v, kmax, klast;
    int i;

    lz = log(hypot(z.re, z.im) + 1e-300);
    klast = m_last_index(a);
    kmax = 4.0 * (hypot(z.re, z.im) + hypot(a.re, a.im) + hypot(b.re, b.im)) + 64.0;
    kmax = FLINT_MIN(kmax, klast);

    lo = 0.0;
    hi = FLINT_MIN(kmax, 2.0);
    kprev = 0.0;
    for (k = 1.0; k <= kmax; k = ceil(k * 1.5))
    {
        v = m_logterm(a, b, lz, k);
        if (v > best)
        {
            best = v;
            kb = k;
            lo = kprev;
            hi = fmin(kmax, ceil(k * 1.5));
        }
        kprev = k;
    }

    /* the terms are large where b + k is close to zero */
    for (i = -1; i <= 2; i++)
    {
        k = floor(-b.re) + i;
        if (k >= 1.0 && k <= kmax)
        {
            v = m_logterm(a, b, lz, k);
            if (v > best)
            {
                best = v;
                kb = k;
                lo = hi = k;
            }
        }
    }

    while (hi - lo > 2.0)
    {
        double m1 = floor(lo + (hi - lo) / 3.0), m2 = ceil(hi - (hi - lo) / 3.0);
        if (m2 <= m1)
            break;
        if (m_logterm(a, b, lz, m1) < m_logterm(a, b, lz, m2))
            lo = m1;
        else
            hi = m2;
    }

    for (k = lo; k <= hi; k += 1.0)
    {
        v = m_logterm(a, b, lz, k);
        if (v > best)
        {
            best = v;
            kb = k;
        }
    }

    *kpk = kb;
    return best;
}

/* The logarithm of the width of a saddle point with the second and third
   derivatives d2 and d3 (absolute values) of the exponent: that of the
   Gaussian, sqrt(2 pi / d2), but no more than that of a double saddle
   point, 2 Gamma(4/3) (6 / d3)^(1/3), which two saddle points approach
   where they merge (a turning point of the function). */
static double
m_saddle_width(double d2, double d3)
{
    double w2 = 0.9189385332046727 - 0.5 * log(d2 + 1e-300);
    double w3 = 0.5798 - log(d3 / 6.0 + 1e-300) / 3.0;

    return FLINT_MIN(w2, w3);
}

/* The number of terms of 1F1(a, b, z) until the terms are smaller than
   exp(target). */
double
_acb_hypgeom_est_1f1_terms(m_cd a, m_cd b, m_cd z, double kpk, double target)
{
    double lz, lo, hi, klast;

    lz = log(hypot(z.re, z.im) + 1e-300);
    klast = m_last_index(a);
    lo = kpk;
    hi = kpk + 16.0;

    while (hi < klast && hi < 1e15 && m_logterm(a, b, lz, hi) > target)
    {
        lo = hi;
        hi = 2.0 * hi + 16.0;
    }

    if (hi >= klast)
    {
        if (m_logterm(a, b, lz, klast) > target)
            return klast + 1.0;
        hi = klast;
    }

    while (hi - lo > 1.0 + 0.02 * lo)
    {
        double mid = floor(0.5 * (lo + hi));
        if (m_logterm(a, b, lz, mid) > target)
            lo = mid;
        else
            hi = mid;
    }

    /* A lower parameter with negative real part: the terms grow again
       where b + k passes zero.  If they matter there, the series is
       summed past that point. */
    if (b.re < 0.0 && floor(-b.re) + 2.0 > hi && floor(-b.re) + 2.0 < klast)
    {
        double kp = floor(-b.re) + 1.0, d, last = -1.0, step = 1.0;

        for (d = 0.0; kp + d < klast && d < 4.0 * (kp + hypot(z.re, z.im)) + 64.0; d = 2.0 * d + 1.0)
        {
            if (m_logterm(a, b, lz, kp + d) > target)
            {
                last = kp + d;
                step = d + 1.0;
            }
        }

        if (last >= 0.0)
        {
            lo = last;
            hi = FLINT_MIN(last + step + 1.0, klast);
            while (hi - lo > 1.0 + 0.02 * lo)
            {
                double mid = floor(0.5 * (lo + hi));
                if (m_logterm(a, b, lz, mid) > target)
                    lo = mid;
                else
                    hi = mid;
            }
        }
    }

    return hi + 1.0;
}

/* Estimates of log |1F1(a, b, z)|, from the saddle points of
   Gamma(b) / (Gamma(a) Gamma(b - a)) int_0^1 exp(phi(t)) dt,
   phi(t) = z t + (a - 1) log(t) + (b - a - 1) log(1 - t)
   (for parameters outside of the domain of the integral, by analytic
   continuation of the approximation).  There are two saddle points, and
   the function is of the size of the integrand at one of them (or at
   both).  With real parameters and argument, which one is known in most
   cases:

   * a > 0, b - a > 0: the maximum of the integrand on (0, 1)
   * a < 0, b - a > 0: the integral over a loop from 1 around 0 crosses
     the negative real axis at a minimum of phi along the axis
   * a > 0, b - a < 0: the same after t -> 1 - t, at a minimum on (1, inf)

   A pair of complex conjugate saddle points gives the same size.  With
   complex numbers, the same conditions on the real parts select the
   saddle point that the real one is continued to, which is right as
   long as the other one does not take over.  Where no saddle point is
   selected, the first estimate is from the higher one and the second
   from the lower one.  Both estimates are
   at most cap, the size of the largest terms of the series. */
void
_acb_hypgeom_est_1f1_logabs(double * first, double * second, m_cd a, m_cd b, m_cd z, double cap)
{
    m_cd a1, c1, ba, p, q, disc, t, u, phi, d2, d3, one_t, z2;
    double C, h, vhi = M_NEGINF, vlo = 1e300, chosen = INFINITY;
    int i, real;

    *first = *second = cap;

    if (hypot(z.re, z.im) < 1e-300)
        return;

    real = (a.im == 0.0 && b.im == 0.0 && z.im == 0.0);

    a1.re = a.re - 1.0; a1.im = a.im;
    ba.re = b.re - a.re; ba.im = b.im - a.im;
    c1.re = ba.re - 1.0; c1.im = ba.im;

    C = _acb_hypgeom_est_lgamma_re(b) - _acb_hypgeom_est_lgamma_re(a) - _acb_hypgeom_est_lgamma_re(ba);

    /* z t^2 - (z - b + 2) t - (a - 1) = 0 */
    p.re = z.re - b.re + 2.0; p.im = z.im - b.im;
    q.re = 4.0 * a1.re; q.im = 4.0 * a1.im;
    disc = _acb_hypgeom_cd_mul(p, p);
    q = _acb_hypgeom_cd_mul(q, z);
    disc.re += q.re; disc.im += q.im;
    disc = _acb_hypgeom_cd_sqrt(disc);
    z2.re = 2.0 * z.re; z2.im = 2.0 * z.im;

    for (i = 0; i < 2; i++)
    {
        u.re = p.re + (i ? -disc.re : disc.re);
        u.im = p.im + (i ? -disc.im : disc.im);
        t = _acb_hypgeom_cd_div(u, z2);
        one_t.re = 1.0 - t.re; one_t.im = -t.im;

        if (!(hypot(t.re, t.im) > 1e-12) || !(hypot(one_t.re, one_t.im) > 1e-12))
            continue;

        phi = _acb_hypgeom_cd_mul(z, t);
        u = _acb_hypgeom_cd_mul(a1, _acb_hypgeom_cd_log(t));
        phi.re += u.re;
        u = _acb_hypgeom_cd_mul(c1, _acb_hypgeom_cd_log(one_t));
        phi.re += u.re;

        /* phi'' = -(a - 1) / t^2 - (b - a - 1) / (1 - t)^2 */
        d2 = _acb_hypgeom_cd_div(a1, _acb_hypgeom_cd_mul(t, t));
        u = _acb_hypgeom_cd_div(c1, _acb_hypgeom_cd_mul(one_t, one_t));
        d2.re = -(d2.re + u.re); d2.im = -(d2.im + u.im);

        /* phi''' = 2 (a - 1) / t^3 - 2 (b - a - 1) / (1 - t)^3 */
        d3 = _acb_hypgeom_cd_div(a1, _acb_hypgeom_cd_mul(t, _acb_hypgeom_cd_mul(t, t)));
        u = _acb_hypgeom_cd_div(c1, _acb_hypgeom_cd_mul(one_t, _acb_hypgeom_cd_mul(one_t, one_t)));
        d3.re = 2.0 * (d3.re - u.re); d3.im = 2.0 * (d3.im - u.im);

        h = C + phi.re + m_saddle_width(hypot(d2.re, d2.im), hypot(d3.re, d3.im));

        if (h != h)
            continue;

        vhi = FLINT_MAX(vhi, h);
        vlo = FLINT_MIN(vlo, h);

        /* (for complex numbers: the saddle point that these are
           continued to) */
        if (real ? (t.im == 0.0) : 1)
        {
            if (a.re > 0.0 && ba.re > 0.0)
            {
                if (t.re > 0.0 && t.re < 1.0)
                    chosen = h;
            }
            else if (a.re < 0.0 && ba.re > 0.0)
            {
                if (t.re < 0.0 && d2.re > 0.0)
                    chosen = h;
            }
            else if (a.re > 0.0 && ba.re < 0.0)
            {
                if (t.re > 1.0 && d2.re > 0.0)
                    chosen = h;
            }
        }
    }

    if (vhi > M_NEGINF)
    {
        if (chosen < INFINITY)
        {
            *first = FLINT_MIN(chosen, cap);
            *second = FLINT_MIN((chosen == vhi) ? vlo : vhi, cap);
        }
        else
        {
            *first = FLINT_MIN(vhi, cap);
            *second = FLINT_MIN(vlo, cap);
        }
    }
}

/* The asymptotic series sum_k (a)_k (a-b+1)_k / (k! (-z)^k): the logarithm
   of its smallest term, and in *n the number of terms down to exp(target)
   (or to the smallest term). */
static double
m_asymp_logterm(m_cd a, m_cd c, double lz, double k)
{
    double u = _acb_hypgeom_est_lrise(a, k), v = _acb_hypgeom_est_lrise(c, k);
    if (u <= M_NEGINF || v <= M_NEGINF)
        return M_NEGINF;
    return u + v - _acb_hypgeom_d_lgamma_re(k + 1.0, 0.0, 0) - k * lz;
}

/* Sets *n to the number of terms and *hump to the logarithm of the
   largest term among them.  The terms need not decrease from the start:
   with large parameters they can grow first (and with parameters of
   large negative real part, the series almost terminates near
   k = -Re(a) or k = -Re(a - b + 1)).  The sum is then affected by
   cancellation of the size of the hump. */
static double
m_asymp_min(double * n, double * hump, m_cd a, m_cd b, m_cd z, double target)
{
    m_cd c;
    double lz, best = 0.0, kb = 0.0, top = 0.0, ktop = 0.0, k, v, kmax, lo, hi;
    int i, j;

    c.re = a.re - b.re + 1.0; c.im = a.im - b.im;
    lz = log(hypot(z.re, z.im) + 1e-300);
    kmax = 2.0 * (hypot(z.re, z.im) + hypot(a.re, a.im) + hypot(c.re, c.im)) + 16.0;

    /* the smallest term */
    for (k = 1.0; k <= kmax; k = ceil(k * 1.3))
    {
        v = m_asymp_logterm(a, c, lz, k);
        if (v < best)
        {
            best = v;
            kb = k;
        }
        if (v <= target)
            break;
    }

    /* (where the series almost terminates) */
    if (best > target)
    {
        for (i = 0; i < 2; i++)
        {
            m_cd x = i ? c : a;

            for (j = 0; j < 3 && x.re < 0.0; j++)
            {
                k = floor(-x.re) + j;
                v = m_asymp_logterm(a, c, lz, k);
                if (v < best)
                {
                    best = v;
                    kb = k;
                }
            }
        }
    }

    /* the largest term before */
    for (k = 1.0; k < kb; k = ceil(k * 1.3))
    {
        v = m_asymp_logterm(a, c, lz, k);
        if (v > top)
        {
            top = v;
            ktop = k;
        }
    }

    /* the first term below the target after it */
    lo = ktop;
    hi = kb;
    if (best <= target)
    {
        while (hi - lo > 1.0 + 0.02 * lo)
        {
            double mid = floor(0.5 * (lo + hi));
            if (m_asymp_logterm(a, c, lz, mid) > target)
                lo = mid;
            else
                hi = mid;
        }
    }

    *hump = top;
    *n = hi + 1.0;
    return best;
}

/* The asymptotic series z^(-a) 2F0(a, a-b+1, -1/z) of U(a, b, z), as
   evaluated by acb_hypgeom_u_asymp: whether the error bound can be
   smaller than exp(-need) times the first term.  The bound is the first
   omitted term times exp(more), where more comes from the factors of
   the bound; the truncation index is limited by the precision.  Sets
   *n to the number of terms, *more to that exponent and *cancel to the
   logarithm of the cancellation in the sum (the terms can grow before
   they decrease; no cancellation if they have the same sign). */
int
_acb_hypgeom_est_u_asymp(double * n, double * more, double * cancel,
    m_cd a, m_cd b, m_cd z, double need, int real)
{
    m_cd u, v;
    double zabs, r, rho0, sigma, nu, sp, alpha, C1, target, lmin, hump;
    int R;

    *n = 1.0;
    *more = 0.0;
    *cancel = 0.0;

    zabs = hypot(z.re, z.im);

    /* |b - 2a| and |2a^2 - 2ab + b| / 2 */
    r = hypot(b.re - 2.0 * a.re, b.im - 2.0 * a.im);
    u.re = 2.0 * (a.re - b.re); u.im = 2.0 * (a.im - b.im);
    v = _acb_hypgeom_cd_mul(u, a);
    rho0 = 0.5 * hypot(v.re + b.re, v.im + b.im);

    if (zabs < 1.0 || zabs <= r)
        return 0;

    if (z.re >= r)
        R = 1;
    else if (fabs(z.im) >= r || z.re >= 0.0)
        R = 2;
    else if (zabs >= 2.0 * r)
        R = 3;
    else
        return 0;

    sigma = r / zabs;
    nu = 1.0 + 2.0 * sigma * sigma;
    sp = (R == 3) ? sigma * nu : sigma;

    if (sp >= 0.99)
        return 0;

    alpha = 1.0 / (1.0 - sp);
    C1 = (R == 1) ? 1.0 : 3.1415926535897932;
    if (R == 3)
        C1 = (C1 + nu * nu * sigma) * nu;

    *more = 2.0 * alpha * C1 * (rho0 + sp * (1.0 + 0.25 * sp) * alpha * alpha) / zabs
        + log(2.0 * alpha) + 3.0;

    /* (with a margin: near the smallest term, the number of terms and
       the accuracy are sensitive to the estimates) */
    target = 1.1 * (-need - *more) - 10.0;
    lmin = m_asymp_min(n, &hump, a, b, z, target);

    /* (if the terms grow first, the sum can be smaller than the largest
       term by about as much as that is larger than 1) */
    if (hump > 0.0 && lmin <= target)
    {
        target -= hump;
        lmin = m_asymp_min(n, &hump, a, b, z, target);
    }

    if (R == 3)
        lmin += *n * log(nu);

    if (lmin > target)
        return 0;

    if (hump > 0.0)
    {
        /* the sign of (a + k) (a - b + 1 + k) / (-z) */
        double k = floor(0.5 * *n);

        if (real && (a.re + k) * (a.re - b.re + 1.0 + k) * (-z.re) > 0.0)
            *cancel = 0.0;
        else
            *cancel = (real ? 2.0 : 1.0) * hump;
    }

    return 1;
}

/* Estimates of log |U(a, b, z)| from the saddle points of
   1 / Gamma(a) int_0^inf exp(phi(t)) dt,
   phi(t) = -z t + (a - 1) log(t) + (b - a - 1) log(1 + t),
   after U(a, b, z) = z^(1-b) U(a-b+1, 2-b, z) if this makes the real
   part of the first parameter positive.  For a > 1 and z > 0, the
   integrand has its maximum at the positive saddle point; with complex
   numbers, a saddle point with Re(z t) > 0 is preferred.  Otherwise the
   first estimate is from the higher saddle point and the second from
   the lower one.  Both are at most cap. */
void
_acb_hypgeom_est_u_logabs(double * first, double * second, m_cd a, m_cd b, m_cd z, double cap)
{
    m_cd a1, c1, p, q, disc, t, u, phi, d2, d3, one_t, z2, lz, lt;
    double C, h, vhi = M_NEGINF, vlo = 1e300, chosen = INFINITY, shift = 0.0, argz;
    int i;

    *first = *second = cap;

    if (hypot(z.re, z.im) < 1e-300)
        return;

    argz = atan2(z.im, z.re);

    if (a.re < 1.0 && a.re - b.re + 1.0 >= 1.0)
    {
        lz = _acb_hypgeom_cd_log(z);
        shift = (1.0 - b.re) * lz.re + b.im * lz.im;
        u.re = a.re - b.re + 1.0; u.im = a.im - b.im;
        b.re = 2.0 - b.re; b.im = -b.im;
        a = u;
    }

    a1.re = a.re - 1.0; a1.im = a.im;
    c1.re = b.re - a.re - 1.0; c1.im = b.im - a.im;

    C = -_acb_hypgeom_est_lgamma_re(a) + shift;

    /* z t^2 + (z - b + 2) t - (a - 1) = 0 */
    p.re = z.re - b.re + 2.0; p.im = z.im - b.im;
    q.re = 4.0 * a1.re; q.im = 4.0 * a1.im;
    disc = _acb_hypgeom_cd_mul(p, p);
    q = _acb_hypgeom_cd_mul(q, z);
    disc.re += q.re; disc.im += q.im;
    disc = _acb_hypgeom_cd_sqrt(disc);
    z2.re = 2.0 * z.re; z2.im = 2.0 * z.im;

    for (i = 0; i < 2; i++)
    {
        u.re = -p.re + (i ? -disc.re : disc.re);
        u.im = -p.im + (i ? -disc.im : disc.im);
        t = _acb_hypgeom_cd_div(u, z2);
        one_t.re = 1.0 + t.re; one_t.im = t.im;

        if (!(hypot(t.re, t.im) > 1e-12) || !(hypot(one_t.re, one_t.im) > 1e-12))
            continue;

        phi = _acb_hypgeom_cd_mul(z, t);
        phi.re = -phi.re;

        /* (the determination of log t which continues the one on the
           ray from 0 on which exp(-z t) decreases, arg t = -arg z) */
        lt = _acb_hypgeom_cd_log(t);
        if (lt.im + argz > 3.141592653589793)
            lt.im -= 2.0 * 3.141592653589793;
        else if (lt.im + argz < -3.141592653589793)
            lt.im += 2.0 * 3.141592653589793;

        u = _acb_hypgeom_cd_mul(a1, lt);
        phi.re += u.re;
        u = _acb_hypgeom_cd_mul(c1, _acb_hypgeom_cd_log(one_t));
        phi.re += u.re;

        d2 = _acb_hypgeom_cd_div(a1, _acb_hypgeom_cd_mul(t, t));
        u = _acb_hypgeom_cd_div(c1, _acb_hypgeom_cd_mul(one_t, one_t));
        d2.re = -(d2.re + u.re); d2.im = -(d2.im + u.im);

        /* phi''' = 2 (a - 1) / t^3 + 2 (b - a - 1) / (1 + t)^3 */
        d3 = _acb_hypgeom_cd_div(a1, _acb_hypgeom_cd_mul(t, _acb_hypgeom_cd_mul(t, t)));
        u = _acb_hypgeom_cd_div(c1, _acb_hypgeom_cd_mul(one_t, _acb_hypgeom_cd_mul(one_t, one_t)));
        d3.re = 2.0 * (d3.re + u.re); d3.im = 2.0 * (d3.im + u.im);

        h = C + phi.re + m_saddle_width(hypot(d2.re, d2.im), hypot(d3.re, d3.im));

        if (h != h)
            continue;

        vhi = FLINT_MAX(vhi, h);
        vlo = FLINT_MIN(vlo, h);

        u = _acb_hypgeom_cd_mul(z, t);
        if (a.re > 0.0 && u.re > 0.0 && (chosen == INFINITY || h > chosen))
            chosen = h;
    }

    if (vhi > M_NEGINF)
    {
        if (chosen < INFINITY)
        {
            *first = FLINT_MIN(chosen, cap);
            *second = FLINT_MIN((chosen == vhi) ? vlo : vhi, cap);
        }
        else
        {
            *first = FLINT_MIN(vhi, cap);
            *second = FLINT_MIN(vlo, cap);
        }

        /* (An estimate far above the cap is that of a saddle point
           which does not contribute: the other one first.) */
        if (*first >= cap && *second < cap
            && ((chosen < INFINITY) ? chosen : vhi) > cap + 10.0)
        {
            *first = *second;
            *second = cap;
        }
    }
}

/* The same for the terms (a)_k (b)_k / ((c)_k k!) z^k of 2F1, |z| < 1:
   the logarithm of the largest term, and the number of terms down to
   exp(target). */
static double
f_logterm(m_cd a, m_cd b, m_cd c, double lz, double k)
{
    double u = _acb_hypgeom_est_lrise(a, k), v = _acb_hypgeom_est_lrise(b, k), w;

    if (u <= M_NEGINF || v <= M_NEGINF)
        return M_NEGINF;

    /* (Beyond a pole of the lower parameter, c = -m: without the
       factor zero, m! (k - m - 1)!.  This is the size of the terms of
       the regularized series times Gamma(m + 1), which is what a
       connection formula has in a limit case.) */
    w = _acb_hypgeom_est_lrise(c, k);
    if (w <= M_NEGINF)
        w = _acb_hypgeom_d_lgamma_re(1.0 - c.re, 0.0, 0) + _acb_hypgeom_d_lgamma_re(k + c.re, 0.0, 0);

    return u + v - w + k * lz - _acb_hypgeom_d_lgamma_re(k + 1.0, 0.0, 0);
}

/* For parameters of moderate size, the terms themselves in double
   precision (a few nanoseconds per term): the squares of their absolute
   values with the exponent kept apart.  With terms = NULL, the
   logarithm of the largest term (its index in kpk); otherwise the
   number of terms down to exp(target) is stored in terms, or 10^15
   if there are more than 10 per bit below the largest term.  Returns
   1e300 if this does not apply (large parameters, or too many
   terms). */
#define F_SCAN_MAX_PARAMETER 256.0
#define F_SCAN_MAX_TERMS 200000.0

static double
f_scan(double * kpk, double * terms, m_cd a, m_cd b, m_cd c, double lz, double target, double kmax)
{
    double x, pm, k, ksafe, z2, num, den, kp, target2, r, k2 = -1.0;
    slong e, pe, cur2;
    int ex, real;

    if (!(hypot(a.re, a.im) <= F_SCAN_MAX_PARAMETER) || !(hypot(b.re, b.im) <= F_SCAN_MAX_PARAMETER)
        || !(hypot(c.re, c.im) <= F_SCAN_MAX_PARAMETER) || !(lz < -0.001) || !(lz > -700.0))
        return 1e300;

    z2 = exp(2.0 * lz);
    /* (beyond the parameters with negative real parts, the ratios of
       the terms are monotone) */
    ksafe = FLINT_MAX(FLINT_MAX(-a.re, -b.re), FLINT_MAX(-c.re, 0.0)) + 2.0;
    target2 = 2.0 * target * M_INVLOG2;

    /* With real parameters, the ratio of the terms beyond ksafe is
       |z| (a + k) (b + k) / ((c + k) (k + 1)), and this minus 1 has the
       sign of (|z| - 1) k^2 + (|z| (a + b) - c - 1) k + (|z| a b - c):
       the terms increase between its roots and decrease beyond the
       larger one, k2.  The largest term is then found without the
       terms between ksafe and k2 (there are (a + b - c - 1) / (1 - |z|)
       of them). */
    real = (a.im == 0.0 && b.im == 0.0 && c.im == 0.0);
    if (real && terms == NULL)
    {
        double zz = exp(lz), q2 = zz - 1.0, q1 = zz * (a.re + b.re) - c.re - 1.0,
            q0 = zz * a.re * b.re - c.re, disc = q1 * q1 - 4.0 * q2 * q0;

        if (disc > 0.0)
            k2 = (-q1 - sqrt(disc)) / (2.0 * q2);
    }

    x = 0.5; e = 1;      /* |t_0|^2 = 1 */
    pm = 0.5; pe = 1; kp = 0.0;

    for (k = 0.0; k < F_SCAN_MAX_TERMS; k += 1.0)
    {
        /* (more terms than the caller is interested in, also while
           they increase) */
        if (terms != NULL && k > kmax)
        {
            *terms = 1e15;
            return 0.0;
        }

        num = ((a.re + k) * (a.re + k) + a.im * a.im) * ((b.re + k) * (b.re + k) + b.im * b.im) * z2;
        den = ((c.re + k) * (c.re + k) + c.im * c.im) * (k + 1.0) * (k + 1.0);

        /* the last term */
        if (num == 0.0)
            break;

        /* (a pole of the lower parameter: without the factor zero, as
           in f_logterm) */
        if (den == 0.0)
            den = (k + 1.0) * (k + 1.0);

        r = num / den;
        x = frexp(x * r, &ex);
        e += ex;
        cur2 = e;

        if (e > pe || (e == pe && x > pm))
        {
            pe = e;
            pm = x;
            kp = k + 1.0;
        }
        else if (terms == NULL)
        {
            if (r < 1.0 && k >= ksafe && cur2 < pe - 60)
                break;
        }

        else if (r < 1.0 && k >= ksafe && cur2 < target2)
        {
            break;
        }
        else if (k - kp > 5.0 * FLINT_MAX(pe - target2, 0.0) + 2.0 * ksafe + 64.0)
        {
            /* More than 10 terms per bit beyond the largest term (the
               squares are counted): more than a series is allowed to
               have. */
            *terms = 1e15;
            return 0.0;
        }

        /* (real parameters: the largest term is this one, or the one
           at the larger root) */
        if (terms == NULL && real && k >= ksafe)
        {
            if (k2 > k + 1.0 && k2 < 1e15)
            {
                double kk = floor(k2) + 1.0, v = f_logterm(a, b, c, lz, kk);

                if (v > 0.5 * (pe + log(pm) * M_INVLOG2) / M_INVLOG2)
                {
                    *kpk = kk;
                    return v;
                }
            }

            break;
        }
    }

    if (k >= F_SCAN_MAX_TERMS)
        return 1e300;

    if (terms != NULL)
    {
        *terms = k + 2.0;
        return 0.0;
    }

    *kpk = kp;
    return 0.5 * (pe + log(pm) * M_INVLOG2) / M_INVLOG2;
}

double
_acb_hypgeom_est_2f1_peak(double * kpk, m_cd a, m_cd b, m_cd c, double lz)
{
    double best = 0.0, kb = 0.0, k, kprev, lo, hi, v, kmax, klast;
    int i;

    v = f_scan(kpk, NULL, a, b, c, lz, 0.0, 1e300);
    if (v < 1e300)
        return v;

    klast = fmin(m_last_index(a), m_last_index(b));
    kmax = 4.0 * (hypot(a.re, a.im) + hypot(b.re, b.im) + hypot(c.re, c.im)) / FLINT_MAX(-lz, 1e-6) + 64.0;
    kmax = FLINT_MIN(kmax, 1e12);
    kmax = FLINT_MIN(kmax, klast);

    lo = 0.0;
    hi = FLINT_MIN(kmax, 2.0);
    kprev = 0.0;
    for (k = 1.0; k <= kmax; k = ceil(k * 1.5))
    {
        v = f_logterm(a, b, c, lz, k);
        if (v > best)
        {
            best = v;
            kb = k;
            lo = kprev;
            hi = fmin(kmax, ceil(k * 1.5));
        }
        kprev = k;
    }

    /* the terms are large where c + k is close to zero */
    for (i = -1; i <= 2; i++)
    {
        k = floor(-c.re) + i;
        if (k >= 1.0 && k <= kmax)
        {
            v = f_logterm(a, b, c, lz, k);
            if (v > best)
            {
                best = v;
                kb = k;
                lo = hi = k;
            }
        }
    }

    while (hi - lo > 2.0)
    {
        double m1 = floor(lo + (hi - lo) / 3.0), m2 = ceil(hi - (hi - lo) / 3.0);
        if (m2 <= m1)
            break;
        if (f_logterm(a, b, c, lz, m1) < f_logterm(a, b, c, lz, m2))
            lo = m1;
        else
            hi = m2;
    }

    for (k = lo; k <= hi; k += 1.0)
    {
        v = f_logterm(a, b, c, lz, k);
        if (v > best)
        {
            best = v;
            kb = k;
        }
    }

    *kpk = kb;
    return best;
}

double
_acb_hypgeom_est_2f1_terms(m_cd a, m_cd b, m_cd c, double lz, double kpk, double target)
{
    return _acb_hypgeom_est_2f1_terms_max(a, b, c, lz, kpk, target, 1e300);
}

/* (the same; 10^15 if there are more than kmax terms and that is found
   out early) */
double
_acb_hypgeom_est_2f1_terms_max(m_cd a, m_cd b, m_cd c, double lz, double kpk, double target, double kmax)
{
    double lo, hi, klast;

    if (f_scan(&kpk, &lo, a, b, c, lz, target, kmax) < 1e300)
        return lo;

    klast = fmin(m_last_index(a), m_last_index(b));
    lo = kpk;
    hi = kpk + 16.0;

    while (hi < klast && hi < 1e15 && f_logterm(a, b, c, lz, hi) > target)
    {
        lo = hi;
        hi = 2.0 * hi + 16.0;
    }

    if (hi >= klast)
    {
        if (f_logterm(a, b, c, lz, klast) > target)
            return klast + 1.0;
        hi = klast;
    }

    while (hi - lo > 1.0 + 0.02 * lo)
    {
        double mid = floor(0.5 * (lo + hi));
        if (f_logterm(a, b, c, lz, mid) > target)
            lo = mid;
        else
            hi = mid;
    }

    /* (past a pole of the lower parameter, if the terms matter there) */
    if (c.re < 0.0 && floor(-c.re) + 2.0 > hi && floor(-c.re) + 2.0 < klast)
    {
        double kp = floor(-c.re) + 1.0, d, last = -1.0, step = 1.0;

        for (d = 0.0; kp + d < klast && d < 64.0 * (kp + 64.0); d = 2.0 * d + 1.0)
        {
            if (f_logterm(a, b, c, lz, kp + d) > target)
            {
                last = kp + d;
                step = d + 1.0;
            }
        }

        if (last >= 0.0)
        {
            lo = last;
            hi = FLINT_MIN(last + step + 1.0, klast);
            while (hi - lo > 1.0 + 0.02 * lo)
            {
                double mid = floor(0.5 * (lo + hi));
                if (f_logterm(a, b, c, lz, mid) > target)
                    lo = mid;
                else
                    hi = mid;
            }
        }
    }

    return hi + 1.0;
}

/* Estimates of log |2F1(a, b, c, z)| from the saddle points of
   Gamma(c) / (Gamma(b) Gamma(c - b)) int_0^1 exp(phi(t)) dt,
   phi(t) = (b - 1) log(t) + (c - b - 1) log(1 - t) - a log(1 - z t)
   (with a and b exchanged if that integral converges and this one does
   not).  If the integral converges, a saddle point with real part in
   (0, 1) is preferred: for real parameters and z < 1, the integrand has
   its maximum there.  Otherwise the first estimate is from the higher
   saddle point and the second from the lower one.  Both are at most
   cap. */
void
_acb_hypgeom_est_2f1_logabs(double * first, double * second,
    m_cd a, m_cd b, m_cd c, m_cd z, double cap)
{
    m_cd b1, c1, cb, q2, q1, disc, t[2], u, v, phi, d2, d3, one_t, one_zt;
    double C, h, vhi = M_NEGINF, vlo = 1e300, chosen = INFINITY;
    int i, nroots, converges;

    *first = *second = cap;

    cb.re = c.re - b.re; cb.im = c.im - b.im;

    if (!(b.re > 0.0 && cb.re > 0.0) && a.re > 0.0 && c.re - a.re > 0.0)
    {
        u = a; a = b; b = u;
        cb.re = c.re - b.re; cb.im = c.im - b.im;
    }

    converges = (b.re > 0.0 && cb.re > 0.0);

    b1.re = b.re - 1.0; b1.im = b.im;
    c1.re = cb.re - 1.0; c1.im = cb.im;

    C = _acb_hypgeom_est_lgamma_re(c) - _acb_hypgeom_est_lgamma_re(b) - _acb_hypgeom_est_lgamma_re(cb);

    /* z (c - a - 2) t^2 + ((2 - c) + z (a - b + 1)) t + (b - 1) = 0 */
    u.re = c.re - a.re - 2.0; u.im = c.im - a.im;
    q2 = _acb_hypgeom_cd_mul(z, u);
    u.re = a.re - b.re + 1.0; u.im = a.im - b.im;
    q1 = _acb_hypgeom_cd_mul(z, u);
    q1.re += 2.0 - c.re; q1.im += -c.im;

    if (hypot(q2.re, q2.im) < 1e-12 * (hypot(q1.re, q1.im) + hypot(b1.re, b1.im)))
    {
        if (hypot(q1.re, q1.im) < 1e-300)
            return;
        u.re = -b1.re; u.im = -b1.im;
        t[0] = _acb_hypgeom_cd_div(u, q1);
        nroots = 1;
    }
    else
    {
        disc = _acb_hypgeom_cd_mul(q1, q1);
        u = _acb_hypgeom_cd_mul(q2, b1);
        disc.re -= 4.0 * u.re; disc.im -= 4.0 * u.im;
        disc = _acb_hypgeom_cd_sqrt(disc);
        v.re = 2.0 * q2.re; v.im = 2.0 * q2.im;
        u.re = -q1.re + disc.re; u.im = -q1.im + disc.im;
        t[0] = _acb_hypgeom_cd_div(u, v);
        u.re = -q1.re - disc.re; u.im = -q1.im - disc.im;
        t[1] = _acb_hypgeom_cd_div(u, v);
        nroots = 2;
    }

    for (i = 0; i < nroots; i++)
    {
        one_t.re = 1.0 - t[i].re; one_t.im = -t[i].im;
        u = _acb_hypgeom_cd_mul(z, t[i]);
        one_zt.re = 1.0 - u.re; one_zt.im = -u.im;

        if (!(hypot(t[i].re, t[i].im) > 1e-12) || !(hypot(one_t.re, one_t.im) > 1e-12)
            || !(hypot(one_zt.re, one_zt.im) > 1e-12))
            continue;

        u = _acb_hypgeom_cd_mul(b1, _acb_hypgeom_cd_log(t[i]));
        phi.re = u.re;
        u = _acb_hypgeom_cd_mul(c1, _acb_hypgeom_cd_log(one_t));
        phi.re += u.re;
        u = _acb_hypgeom_cd_mul(a, _acb_hypgeom_cd_log(one_zt));
        phi.re -= u.re;

        /* phi'' = -(b-1)/t^2 - (c-b-1)/(1-t)^2 + a z^2/(1-zt)^2 */
        d2 = _acb_hypgeom_cd_div(b1, _acb_hypgeom_cd_mul(t[i], t[i]));
        u = _acb_hypgeom_cd_div(c1, _acb_hypgeom_cd_mul(one_t, one_t));
        d2.re = -(d2.re + u.re); d2.im = -(d2.im + u.im);
        u = _acb_hypgeom_cd_div(z, one_zt);
        u = _acb_hypgeom_cd_mul(u, u);
        u = _acb_hypgeom_cd_mul(u, a);
        d2.re += u.re; d2.im += u.im;

        /* phi''' = 2 (b-1)/t^3 - 2 (c-b-1)/(1-t)^3 + 2 a z^3/(1-zt)^3 */
        d3 = _acb_hypgeom_cd_div(b1, _acb_hypgeom_cd_mul(t[i], _acb_hypgeom_cd_mul(t[i], t[i])));
        u = _acb_hypgeom_cd_div(c1, _acb_hypgeom_cd_mul(one_t, _acb_hypgeom_cd_mul(one_t, one_t)));
        d3.re -= u.re; d3.im -= u.im;
        u = _acb_hypgeom_cd_div(z, one_zt);
        u = _acb_hypgeom_cd_mul(u, _acb_hypgeom_cd_mul(u, u));
        u = _acb_hypgeom_cd_mul(u, a);
        d3.re = 2.0 * (d3.re + u.re); d3.im = 2.0 * (d3.im + u.im);

        h = C + phi.re + m_saddle_width(hypot(d2.re, d2.im), hypot(d3.re, d3.im));

        if (h != h)
            continue;

        vhi = FLINT_MAX(vhi, h);
        vlo = FLINT_MIN(vlo, h);

        if (converges && t[i].re > 0.0 && t[i].re < 1.0 && (chosen == INFINITY || h > chosen))
            chosen = h;
    }

    if (vhi > M_NEGINF)
    {
        if (chosen < INFINITY)
        {
            *first = FLINT_MIN(chosen, cap);
            *second = FLINT_MIN((chosen == vhi) ? vlo : vhi, cap);
        }
        else
        {
            *first = FLINT_MIN(vhi, cap);
            *second = FLINT_MIN(vlo, cap);
        }
    }

    /* No saddle point between the end points of a convergent integral:
       the integrand can have its maximum at an end point, where it
       decreases like t^(b-1) exp(-lambda t) with
       lambda = (c - b - 1) - a z, or like s^(c-b-1) exp(-mu s),
       s = 1 - t, with mu = (b - 1) + a z / (1 - z); the integral is then
       about Gamma(b) lambda^(-b), or Gamma(c - b) mu^(b-c) (1 - z)^(-a).
       (For example a = 1000, b = 1, c = 2000, z = 0.99: the saddle
       points are beyond 1, where the integrand is 2^-4677, and the
       function is 1.98.) */
    if (converges && chosen == INFINITY)
    {
        double e = M_NEGINF;
        m_cd lam, mu, one_z;

        u = _acb_hypgeom_cd_mul(a, z);
        lam.re = c1.re - u.re; lam.im = c1.im - u.im;

        if (lam.re > 0.0)
        {
            v = _acb_hypgeom_cd_mul(b, _acb_hypgeom_cd_log(lam));
            h = _acb_hypgeom_est_lgamma_re(c) - _acb_hypgeom_est_lgamma_re(cb) - v.re;
            if (h == h)
                e = FLINT_MAX(e, h);
        }

        one_z.re = 1.0 - z.re; one_z.im = -z.im;

        if (hypot(one_z.re, one_z.im) > 1e-12)
        {
            mu = _acb_hypgeom_cd_div(u, one_z);
            mu.re += b1.re; mu.im += b1.im;

            if (mu.re > 0.0)
            {
                v = _acb_hypgeom_cd_mul(cb, _acb_hypgeom_cd_log(mu));
                u = _acb_hypgeom_cd_mul(a, _acb_hypgeom_cd_log(one_z));
                h = _acb_hypgeom_est_lgamma_re(c) - _acb_hypgeom_est_lgamma_re(b) - v.re - u.re;
                if (h == h)
                    e = FLINT_MAX(e, h);
            }
        }

        if (e > M_NEGINF && (vhi == M_NEGINF || e > vhi))
        {
            *second = *first;
            *first = FLINT_MIN(e, cap);
        }
    }
}

/* The connection formulas of 2F1 (which = 2, ..., 5 as in
   acb_hypgeom_2f1_transform: the arguments 1/z, 1/(1 - z), 1 - z and
   1 - 1/z), each a sum of two series in w with factors of gamma
   functions and powers (DLMF 15.8.2 to 15.8.5): the parameters of the
   two series, the logarithms of their factors and of their largest
   terms, and Lam, the logarithm of the largest term of the formula.
   In a limit case (b - a an integer for the first two, c - a - b for
   the others) the two terms have poles which cancel; what remains has
   the size of the terms of the two series with the gamma function at
   its pole replaced by its size away from the poles, and with the
   Pochhammer symbol of a lower parameter beyond its pole without the
   factor zero (times a logarithm).
   Returns 0 if the formula does not apply (|w| is not below 0.99, or
   the estimate is not finite). */
static double
f_lgamma(m_cd x, int limit)
{
    if (x.im == 0.0 && x.re <= 0.0 && x.re == floor(x.re) && !limit)
        return 1e300;

    return _acb_hypgeom_d_lgamma_re(x.re, x.im, limit);
}

static double
f_re_mul(m_cd x, m_cd y)
{
    return x.re * y.re - x.im * y.im;
}

int
_acb_hypgeom_est_2f1_conn(_acb_hypgeom_est_2f1_conn_struct * E, int which,
    m_cd a, m_cd b, m_cd c, m_cd z)
{
    m_cd d, ca, cb, cab, one_z, mz, lz, l1z, lmz, num[2], den0[2], den1[2], t;
    double g[2];
    int j, limit;

    one_z.re = 1.0 - z.re; one_z.im = -z.im;
    mz.re = -z.re; mz.im = -z.im;

    if (!(hypot(z.re, z.im) > 1e-300) || !(hypot(one_z.re, one_z.im) > 1e-300))
        return 0;

    lz = _acb_hypgeom_cd_log(z);
    l1z = _acb_hypgeom_cd_log(one_z);
    lmz = _acb_hypgeom_cd_log(mz);

    ca.re = c.re - a.re; ca.im = c.im - a.im;
    cb.re = c.re - b.re; cb.im = c.im - b.im;
    cab.re = ca.re - b.re; cab.im = ca.im - b.im;

    if (which == 2)
        E->lw = -lz.re;
    else if (which == 3)
        E->lw = -l1z.re;
    else if (which == 4)
        E->lw = l1z.re;
    else
        E->lw = l1z.re - lz.re;

    /* (not close to the circle of convergence: log(0.99)) */
    if (!(E->lw < -0.01))
        return 0;

    if (which <= 3)
    {
        /* d = b - a; Gamma(c) Gamma(b - a) / (Gamma(b) Gamma(c - a)) and
           Gamma(c) Gamma(a - b) / (Gamma(a) Gamma(c - b)) */
        d.re = b.re - a.re; d.im = b.im - a.im;
        num[0] = d;
        num[1].re = -d.re; num[1].im = -d.im;
        den0[0] = b; den1[0] = ca;
        den0[1] = a; den1[1] = cb;

        E->a[0] = a; E->a[1] = b;
        E->c[0].re = 1.0 - d.re; E->c[0].im = -d.im;
        E->c[1].re = 1.0 + d.re; E->c[1].im = d.im;

        if (which == 2)
        {
            /* (-z)^(-a) F(a, a - c + 1, a - b + 1, 1/z) and the same
               with a and b exchanged */
            E->b[0].re = 1.0 - ca.re; E->b[0].im = -ca.im;
            E->b[1].re = 1.0 - cb.re; E->b[1].im = -cb.im;
            E->pref[0] = -f_re_mul(a, lmz);
            E->pref[1] = -f_re_mul(b, lmz);
        }
        else
        {
            /* (1 - z)^(-a) F(a, c - b, a - b + 1, 1/(1 - z)) and the
               same with a and b exchanged */
            E->b[0] = cb;
            E->b[1] = ca;
            E->pref[0] = -f_re_mul(a, l1z);
            E->pref[1] = -f_re_mul(b, l1z);
        }
    }
    else
    {
        /* d = c - a - b; Gamma(c) Gamma(c - a - b) / (Gamma(c - a)
           Gamma(c - b)) and Gamma(c) Gamma(a + b - c) / (Gamma(a)
           Gamma(b)) */
        d = cab;
        num[0] = d;
        num[1].re = -d.re; num[1].im = -d.im;
        den0[0] = ca; den1[0] = cb;
        den0[1] = a; den1[1] = b;

        E->c[0].re = 1.0 - d.re; E->c[0].im = -d.im;
        E->c[1].re = 1.0 + d.re; E->c[1].im = d.im;
        E->a[0] = a;
        E->a[1] = ca;

        if (which == 4)
        {
            /* F(a, b, a + b - c + 1, 1 - z) and (1 - z)^(c - a - b)
               F(c - a, c - b, c - a - b + 1, 1 - z) */
            E->b[0] = b;
            E->b[1] = cb;
            E->pref[0] = 0.0;
            E->pref[1] = f_re_mul(cab, l1z);
        }
        else
        {
            /* z^(-a) F(a, a - c + 1, a + b - c + 1, 1 - 1/z) and
               (1 - z)^(c - a - b) z^(a - c) F(c - a, 1 - a,
               c - a - b + 1, 1 - 1/z) */
            E->b[0].re = 1.0 - ca.re; E->b[0].im = -ca.im;
            E->b[1].re = 1.0 - a.re; E->b[1].im = -a.im;
            E->pref[0] = -f_re_mul(a, lz);
            E->pref[1] = f_re_mul(cab, l1z) - f_re_mul(ca, lz);
        }
    }

    limit = (d.im == 0.0 && d.re == floor(d.re));
    E->limit = limit;

    g[0] = _acb_hypgeom_d_lgamma_re(c.re, c.im, 0);
    E->Lam = M_NEGINF;

    for (j = 0; j < 2; j++)
    {
        double u = f_lgamma(den0[j], 0), v = f_lgamma(den1[j], 0);

        E->pk[j] = 0.0;
        E->kpk[j] = 0.0;

        /* (a term which vanishes) */
        if (u >= 1e300 || v >= 1e300)
        {
            E->pref[j] = M_NEGINF;
            continue;
        }

        E->pref[j] += g[0] + f_lgamma(num[j], limit) - u - v;

        if (limit)
        {
            /* (the logarithm and the digamma functions of the limit) */
            E->pref[j] += 1.0 + log(1.0 + fabs(E->lw) + fabs(lz.re));

            /* The limit is a derivative with respect to the parameters:
               a series which terminates continues with terms of the
               size (-m)_m (k - m - 1)!, which is about that of
               (1/2 - m)_k. */
            if (E->a[j].im == 0.0 && E->a[j].re <= 0.0 && E->a[j].re == floor(E->a[j].re))
                E->a[j].re += 0.5;
            if (E->b[j].im == 0.0 && E->b[j].re <= 0.0 && E->b[j].re == floor(E->b[j].re))
                E->b[j].re += 0.5;
        }

        t = E->c[j];
        E->pk[j] = _acb_hypgeom_est_2f1_peak(E->kpk + j, E->a[j], E->b[j], t, E->lw);

        if (E->pref[j] != E->pref[j] || E->pk[j] != E->pk[j])
            return 0;

        E->Lam = FLINT_MAX(E->Lam, E->pref[j] + E->pk[j]);
    }

    return E->Lam > M_NEGINF && E->Lam < 1e299;
}

/* The number of terms of the two series down to exp(level) (the
   logarithm of an absolute value, as Lam), and in tper the larger one
   of the two numbers (10^15 if a series has more than kmax terms) */
double
_acb_hypgeom_est_2f1_conn_terms(double * tper, const _acb_hypgeom_est_2f1_conn_struct * E, double level, double kmax)
{
    double T[2];
    int j;

    for (j = 0; j < 2; j++)
    {
        T[j] = 0.0;

        if (E->pref[j] <= M_NEGINF)
            continue;

        T[j] = _acb_hypgeom_est_2f1_terms_max(E->a[j], E->b[j], E->c[j], E->lw, E->kpk[j],
                level - E->pref[j], kmax);

        /* (In a limit case, the series with the lower parameter 1 - m
           is summed past the pole at k = m, where the terms of the
           limit continue with the factor w^m: at least m + 1 terms.) */
        if (E->limit && E->c[j].re <= 0.25 && T[j] < 2.0 - E->c[j].re)
        {
            T[j] = 2.0 - E->c[j].re;

            if (T[j] > kmax)
                T[j] = 1e15;
        }
    }

    *tper = FLINT_MAX(T[0], T[1]);

    return T[0] + T[1];
}

/* Costs in microseconds, measured on one machine; only the ratios matter. */

/* One term of the 1F1 series at precision wp with parameters of the
   given number of bits. */
double
_acb_hypgeom_est_series_term_cost(double wp, double bits, int real)
{
    double extra;

    bits = FLINT_MIN(bits, wp);
    if (bits <= 64.0)
        bits = 0.0;

    /* (a working precision above that of the parameters, for series
       that cancel: products of numbers of both sizes; 3363-bit
       parameters at 8000, 16000 and 28000 bits take 6, 11 and 35
       microseconds) */
    extra = 1.5e-7 * bits * (wp - bits);
    if (wp > 16000.0)
        extra *= (wp / 16000.0) * (wp / 16000.0);

    /* (at most three multiplications at the working precision, which
       cost about 6e-5 wp log2(wp) microseconds each from 16000 bits on:
       the fit above grows like wp^3 beyond the measurements, and the
       series must remain cheaper than an integral, whose cost grows
       like wp^2.7, at a high enough precision) */
    extra = fmin(extra, 1.8e-4 * wp * log2(wp + 2.0));

    if (real)
        return 0.25 + 0.006 * sqrt(wp) + 1.5e-5 * wp + 7e-4 * bits + extra;
    else
        return 1.1 + 0.018 * sqrt(wp) + 4e-5 * wp + 1.8e-3 * bits + 2.2 * extra;
}

/* Whether the two saddle points of exp(z t) t^(b-a-1) (1 + t)^(a-1),
   the roots of z t^2 + (z + b - 2) t + (b - a - 1), are not both close
   to the real line: the two integrals behind 1F1 and U for Re(a) < 1
   (_acb_hypgeom_2u_integration) are of the size of the function then
   (it oscillates for real arguments), and cancel otherwise. */
int
_acb_hypgeom_est_2u_oscillatory(_acb_hypgeom_cd a, _acb_hypgeom_cd b, _acb_hypgeom_cd z)
{
    _acb_hypgeom_cd p, q, D;

    p.re = z.re + b.re - 2.0; p.im = z.im + b.im;
    q.re = b.re - a.re - 1.0; q.im = b.im - a.im;
    D = _acb_hypgeom_cd_mul(p, p);
    q = _acb_hypgeom_cd_mul(z, q);
    D.re -= 4.0 * q.re;
    D.im -= 4.0 * q.im;

    return (D.re < 0.0) || (fabs(D.im) > fabs(D.re));
}

/* Numerical integration at precision prec. */
double
_acb_hypgeom_est_integration_cost(double prec, int real)
{
    /* (with the nodes in the cache) */
    static const double tab_real[9] = { 80.0, 125.0, 520.0, 1500.0, 5400.0, 39000.0, 3.2e5, 1.8e6, 1.6e7 };
    static const double tab_cplx[9] = { 410.0, 790.0, 2600.0, 5600.0, 22000.0, 1.1e5, 8.1e5, 6.2e6, 6.0e7 };
    const double * tab = real ? tab_real : tab_cplx;
    double x;
    int i;

    /* the table is for prec = 64 * 2^i */
    x = log(FLINT_MAX(prec, 64.0) / 64.0) * M_INVLOG2;
    if (x >= 8.0)
        x = tab[8] * pow(9.0, x - 8.0);
    else
    {
        i = (int) x;
        x = tab[i] * pow(tab[i + 1] / tab[i], x - i);
    }

    /* The quadrature nodes are computed on the first call at a given
       precision, at several times the cost of an integral, and cached.
       A part of this is counted: the series are preferred for a single
       evaluation, and integration if it is much faster. */
    if (prec >= 256.0)
        x += 0.1 * 60000.0 * pow(prec / 1024.0, 2.6);

    return x;
}

/* The relative accuracy that an inexact parameter x is good for in a
   function f with d log(f) / dx of order 1: its relative accuracy
   minus log2 |x|.  (WORD_MAX for an exact x; 0 from |x| = 2^60, the
   same on 32-bit and 64-bit systems.) */
slong
_acb_hypgeom_est_input_accuracy(const acb_t x)
{
    slong acc;
    mag_t m;

    if (acb_is_exact(x))
        return WORD_MAX;

    acc = acb_rel_accuracy_bits(x);
    acc = FLINT_MAX(acc, 0);

    mag_init(m);
    acb_get_mag(m, x);
    if (mag_cmp_2exp_si(m, 0) > 0)
    {
        if (mag_cmp_2exp_si(m, 60) < 0)
            acc -= FLINT_MIN(MAG_EXP(m), acc);
        else
            acc = 0;
    }
    mag_clear(m);

    return acc;
}

/* Whether the two saddle points of t^(b-1) (1 - t)^(a-c) (1 - x t)^(-a),
   x = 1 - z, the roots of x (b - c) t^2 + ((a - b) x - (a + b - c)) t + b
   (for large parameters), are not both close to the real line: the two
   integrals behind 2F1 for Re(z) < 0 (_acb_hypgeom_2f1_integration_left)
   are of the size of the function then (it oscillates for real
   arguments), and are not found or cancel otherwise.  The condition is
   symmetric in a and b. */
int
_acb_hypgeom_est_2f1_left_oscillatory(_acb_hypgeom_cd a, _acb_hypgeom_cd b, _acb_hypgeom_cd c,
    _acb_hypgeom_cd z)
{
    _acb_hypgeom_cd x, p, q, D;

    x.re = 1.0 - z.re; x.im = -z.im;
    p.re = a.re - b.re; p.im = a.im - b.im;
    p = _acb_hypgeom_cd_mul(p, x);
    p.re -= a.re + b.re - c.re;
    p.im -= a.im + b.im - c.im;
    q.re = b.re - c.re; q.im = b.im - c.im;
    q = _acb_hypgeom_cd_mul(_acb_hypgeom_cd_mul(x, b), q);
    D = _acb_hypgeom_cd_mul(p, p);
    D.re -= 4.0 * q.re;
    D.im -= 4.0 * q.im;

    return (D.re < 0.0) || (fabs(D.im) > fabs(D.re));
}
