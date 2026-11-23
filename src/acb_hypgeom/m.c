/*
    Copyright (C) 2015 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <math.h>
#include "acb.h"
#include "arb_hypgeom.h"
#include "acb_hypgeom.h"
#include "acb_hypgeom/impl.h"

void
acb_hypgeom_m_asymp(acb_t res, const acb_t a, const acb_t b, const acb_t z, int regularized, slong prec)
{
    acb_t t, u, v, c;

    acb_init(t);
    acb_init(u);
    acb_init(v);
    acb_init(c);

    acb_sub(c, b, a, prec);
    acb_neg(v, z);

    acb_hypgeom_u_asymp(t, a, b, z, -1, prec);
    acb_hypgeom_u_asymp(u, c, b, v, -1, prec);

    /* gamma(b-a) */
    acb_rgamma(v, c, prec);
    acb_mul(t, t, v, prec);

    /* z^(a-b) */
    acb_neg(c, c);
    acb_pow(v, z, c, prec);
    acb_mul(u, u, v, prec);

    /* gamma(a) */
    acb_rgamma(v, a, prec);
    acb_mul(u, u, v, prec);

    /* exp(z) */
    acb_exp(v, z, prec);
    acb_mul(u, u, v, prec);

    /* (-z)^(-a) */
    acb_neg(c, a);
    acb_neg(v, z);
    acb_pow(v, v, c, prec);
    acb_mul(t, t, v, prec);

    acb_add(t, t, u, prec);

    if (!regularized)
    {
        acb_gamma(v, b, prec);

        if (acb_is_finite(v))
            acb_mul(t, t, v, prec);
        else
            acb_indeterminate(t);
    }

    if (acb_is_real(a) && acb_is_real(b) && acb_is_real(z))
    {
        arb_zero(acb_imagref(t));
    }

    acb_swap(res, t);

    acb_clear(t);
    acb_clear(u);
    acb_clear(v);
    acb_clear(c);
}

static void
_acb_hypgeom_m_1f1(acb_t res, const acb_t a, const acb_t b, const acb_t z,
    int regularized, slong prec, slong gamma_prec, int kummer, slong n_max)
{
    if (regularized)
    {
        /* Remove singularity */
        if (acb_is_int(b) && arb_is_nonpositive(acb_realref(b)) &&
            arf_cmpabs_2exp_si(arb_midref(acb_realref(b)), 30) < 0)
        {
            acb_t c, d, t, u;
            slong n;

            n = arf_get_si(arb_midref(acb_realref(b)), ARF_RND_DOWN);

            acb_init(c);
            acb_init(d);
            acb_init(t);
            acb_init(u);

            acb_sub(c, a, b, prec);
            acb_add_ui(c, c, 1, prec);

            acb_neg(d, b);
            acb_add_ui(d, d, 2, prec);

            _acb_hypgeom_m_1f1(t, c, d, z, 0, prec, gamma_prec, kummer, n_max);

            acb_pow_ui(u, z, 1 - n, prec);
            acb_mul(t, t, u, prec);

            acb_rising_ui(u, a, 1 - n, prec);
            acb_mul(t, t, u, prec);

            arb_fac_ui(acb_realref(u), 1 - n, prec);
            acb_div_arb(res, t, acb_realref(u), prec);

            acb_clear(c);
            acb_clear(d);
            acb_clear(t);
            acb_clear(u);
        }
        else
        {
            acb_t t;
            acb_init(t);
            acb_rgamma(t, b, gamma_prec);
            _acb_hypgeom_m_1f1(res, a, b, z, 0, prec, gamma_prec, kummer, n_max);
            acb_mul(res, res, t, prec);
            acb_clear(t);
        }
        return;
    }

    /* Kummer's transformation */
    if (kummer)
    {
        acb_t u, v;
        acb_init(u);
        acb_init(v);

        acb_sub(u, b, a, prec);
        acb_neg(v, z);

        _acb_hypgeom_m_1f1(u, u, b, v, regularized, prec, gamma_prec, 0, n_max);
        acb_exp(v, z, prec);
        acb_mul(res, u, v, prec);

        acb_clear(u);
        acb_clear(v);
        return;
    }

    /* (n_max: the largest number of terms; by default, a multiple of
       the precision) */
    n_max = FLINT_MAX(n_max, 50 + 10 * prec);

    if (acb_is_one(a))
    {
        _acb_hypgeom_pfq_direct_max(res, NULL, 0, b, 1, z, prec, n_max);
    }
    else
    {
        acb_struct c[3];
        c[0] = *a;
        c[1] = *b;

        acb_init(c + 2);
        acb_one(c + 2);

        _acb_hypgeom_pfq_direct_max(res, c, 1, c + 1, 2, z, prec, n_max);

        acb_clear(c + 2);
    }
}

void
acb_hypgeom_m_1f1(acb_t res, const acb_t a, const acb_t b, const acb_t z, int regularized, slong prec)
{
    if (arf_sgn(arb_midref(acb_realref(z))) >= 0
        || (acb_is_int(a) && arb_is_nonpositive(acb_realref(a))))
    {
        _acb_hypgeom_m_1f1(res, a, b, z, regularized, prec, prec, 0, 0);
    }
    else
    {
        _acb_hypgeom_m_1f1(res, a, b, z, regularized, prec, prec, 1, 0);
    }
}

static double
hypotmx(double x, double y)
{
    if (x > 0.0 && x > 1e6 * fabs(y))
        return y * y / (2.0 * x);
    else
        return sqrt(x * x + y * y) - x;
}

static void
acb_hypgeom_m_choose(int * asymp, int * kummer, slong * wp,
    const acb_t a, const acb_t b, const acb_t z, int regularized, slong prec)
{
    double x, y, t, cancellation;
    double input_accuracy, direct_accuracy, asymp_accuracy;
    slong m = WORD_MAX;
    slong n = WORD_MAX;

    if (acb_is_int(a) &&
            arf_cmpabs_2exp_si(arb_midref(acb_realref(a)), 30) < 0)
    {
        m = arf_get_si(arb_midref(acb_realref(a)), ARF_RND_DOWN);
    }

    if (acb_is_int(b) &&
            arf_cmpabs_2exp_si(arb_midref(acb_realref(b)), 30) < 0)
    {
        n = arf_get_si(arb_midref(acb_realref(b)), ARF_RND_DOWN);
    }

    *asymp = 0;
    *kummer = 0;
    *wp = prec;

    /* The 1F1 series terminates. */
    /* TODO: for large m, estimate extra precision here. */
    if (m <= 0 && m < n && m > -10 * prec && (n > 0 || !regularized))
    {
        *asymp = 0;
        return;
    }

    /* The 1F1 series terminates with the Kummer transform. */
    /* TODO: for large m, estimate extra precision here. */
    if (m >= 1 && n >= 1 && m < 0.1 * prec && n < 0.1 * prec && n <= m)
    {
        *asymp = 0;
        *kummer = 1;
        return;
    }

    input_accuracy = acb_rel_one_accuracy_bits(z);
    t = acb_rel_one_accuracy_bits(a);
    input_accuracy = FLINT_MIN(input_accuracy, t);
    t = acb_rel_one_accuracy_bits(b);
    input_accuracy = FLINT_MIN(input_accuracy, t);
    input_accuracy = FLINT_MAX(input_accuracy, 0.0);

    /* From here we ignore the values of a, b. Taking them into account is
       a possible future improvement... */

    /* Tiny |z|. */
    if ((arf_cmpabs_2exp_si(arb_midref(acb_realref(z)), 2) < 0 &&
         arf_cmpabs_2exp_si(arb_midref(acb_imagref(z)), 2) < 0))
    {
        *asymp = 0;
        *wp = FLINT_MAX(2, FLINT_MIN(input_accuracy + 20, prec));
        return;
    }

    /* Huge |z|. */
    if ((arf_cmpabs_2exp_si(arb_midref(acb_realref(z)), 64) > 0 ||
         arf_cmpabs_2exp_si(arb_midref(acb_imagref(z)), 64) > 0))
    {
        *asymp = 1;
        *wp = FLINT_MAX(2, FLINT_MIN(input_accuracy + 20, prec));
        return;
    }

    x = arf_get_d(arb_midref(acb_realref(z)), ARF_RND_DOWN);
    y = arf_get_d(arb_midref(acb_imagref(z)), ARF_RND_DOWN);

    asymp_accuracy = sqrt(x * x + y * y) * 1.44269504088896 - 5.0;

    /* The Kummer transformation gives less cancellation with the 1F1 series. */
    if (x < 0.0)
    {
        *kummer = 1;
        x = -x;
    }

    if (asymp_accuracy >= prec)
    {
        *asymp = 1;
        *wp = FLINT_MAX(2, FLINT_MIN(input_accuracy + 20, prec));
        return;
    }

    cancellation = hypotmx(x, y) * 1.44269504088896;

    direct_accuracy = input_accuracy - cancellation;

    if (direct_accuracy > asymp_accuracy)
    {
        *asymp = 0;
        *wp = FLINT_MAX(2, FLINT_MIN(input_accuracy + 20, prec + cancellation));
    }
    else
    {
        *asymp = 1;
        *wp = FLINT_MAX(2, FLINT_MIN(input_accuracy + 20, prec));
    }
}

void
_acb_hypgeom_m_nointegration(acb_t res, const acb_t a, const acb_t b, const acb_t z, int regularized, slong prec)
{
    int asymp, kummer;
    slong wp;

    acb_hypgeom_m_choose(&asymp, &kummer, &wp, a, b, z, regularized, prec);

    if (asymp)
    {
        acb_hypgeom_m_asymp(res, a, b, z, regularized, wp);
    }
    else
    {
        _acb_hypgeom_m_1f1(res, a, b, z, regularized, wp, FLINT_MIN(wp, prec), kummer, 0);
    }

    acb_set_round(res, res, prec);
}

/* ------------------------------------------------------------------------- */
/* Choice of algorithm from estimates of the costs                            */
/* ------------------------------------------------------------------------- */

/* The estimates of sizes, numbers of terms and costs (estimates.c) are
   computed in double precision from the midpoints of the input.  They
   decide what is tried and with how much working precision; the accuracy
   of each result is checked afterwards. */

#define M_INVLOG2 1.4426950408889634

typedef _acb_hypgeom_cd m_cd;

#define m_cd_log _acb_hypgeom_cd_log
#define m_lgre _acb_hypgeom_est_lgamma_re
#define m_peak _acb_hypgeom_est_1f1_peak
#define m_num_terms _acb_hypgeom_est_1f1_terms
#define m_logabs _acb_hypgeom_est_1f1_logabs
#define m_series_term_cost _acb_hypgeom_est_series_term_cost
#define m_integration_cost _acb_hypgeom_est_integration_cost

static double
m_arb_bits(const arb_t x)
{
    return (double) arf_bits(arb_midref(x));
}

/* ------------------------------------------------------------------------- */
/* Inexact input                                                               */
/* ------------------------------------------------------------------------- */

/* Bound for |1F1| (regularized or not) on balls A, B, Z with
   Re(A) > 0 and Re(B - A) > 0, from the Euler integral:

   |1F1(a, b, z)| <= G 1F1(Re a, Re b, Re z),
   G = |Gamma(b)| Gamma(Re a) Gamma(Re(b - a))
          / (|Gamma(a)| |Gamma(b - a)| Gamma(Re b)),

   since |exp(z t) t^(a-1) (1-t)^(b-a-1)| is the integrand at the real
   parts.  Here |Gamma(b)| <= Gamma(Re b), and

   Gamma(x) / |Gamma(x + iy)| = prod_k (1 + y^2 / (x + k)^2)^(1/2)
                              <= exp(y^2 psi'(x) / 2),
   psi'(x) = sum_k (x + k)^(-2) <= 1 / x + 1 / x^2.

   1F1(al, be, x) is the mean of exp(x t) for the beta distribution with
   parameters (al, be - al), which increases stochastically with al and
   decreases with be: 1F1 increases with x, and for x >= 0 it increases
   with al and decreases with be (the opposite for x <= 0).  Its largest
   value on the box of real parts is at a corner.  This bound is close
   to |1F1| unless z or the parameters are far from real. */
static void
m_majorant_euler(mag_t M, const acb_t A, const acb_t B, const acb_t Z, int regularized)
{
    acb_t BA, t, al, be, x;
    arb_t g, u, w;
    arf_t c;
    slong prec = 40;
    int i;

    acb_init(BA); acb_init(t); acb_init(al); acb_init(be); acb_init(x);
    arb_init(g); arb_init(u); arb_init(w);
    arf_init(c);

    acb_sub(BA, B, A, prec);

    if (arb_is_positive(acb_realref(A)) && arb_is_positive(acb_realref(BA))
        && acb_is_finite(A) && acb_is_finite(B) && acb_is_finite(Z))
    {
        /* log G <= sum over a and b - a of y^2 (1 / x + 1 / x^2) / 2 */
        arb_zero(g);
        for (i = 0; i < 2; i++)
        {
            const acb_struct * p = i ? BA : A;

            arb_get_lbound_arf(c, acb_realref(p), prec);
            arb_set_arf(u, c);
            arb_inv(u, u, prec);
            arb_mul(w, u, u, prec);
            arb_add(u, u, w, prec);
            arb_mul(w, acb_imagref(p), acb_imagref(p), prec);
            arb_mul(u, u, w, prec);
            arb_mul_2exp_si(u, u, -1);
            arb_add(g, g, u, prec);
        }
        arb_exp(g, g, prec);

        /* 1 / Gamma(Re b) for the regularized function */
        if (regularized)
        {
            arb_rgamma(u, acb_realref(B), prec);
            arb_mul(g, g, u, prec);
        }

        /* the corner */
        arb_get_ubound_arf(c, acb_realref(Z), prec);
        acb_zero(x); arb_set_arf(acb_realref(x), c);

        if (arf_sgn(c) >= 0)
        {
            arb_get_ubound_arf(c, acb_realref(A), prec);
            acb_zero(al); arb_set_arf(acb_realref(al), c);
            arb_get_lbound_arf(c, acb_realref(B), prec);
            acb_zero(be); arb_set_arf(acb_realref(be), c);
        }
        else
        {
            arb_get_lbound_arf(c, acb_realref(A), prec);
            acb_zero(al); arb_set_arf(acb_realref(al), c);
            arb_get_ubound_arf(c, acb_realref(B), prec);
            acb_zero(be); arb_set_arf(acb_realref(be), c);
        }

        acb_hypgeom_m(t, al, be, x, 0, prec);
        acb_mul_arb(t, t, g, prec);
        acb_get_mag(M, t);
    }
    else
    {
        mag_inf(M);
    }

    acb_clear(BA); acb_clear(t); acb_clear(al); acb_clear(be); acb_clear(x);
    arb_clear(g); arb_clear(u); arb_clear(w);
    arf_clear(c);
}

/* Bound for |1F1| on balls from the absolute values of the terms of
   the series, and of the series of exp(z) 1F1(b - a, b, -z): close to
   |1F1| where the terms of one of them have equal signs (the domain of
   the Euler integral with real arguments), and too large by the
   cancellation otherwise. */
static void
m_majorant_terms(mag_t M, const acb_t A, const acb_t B, const acb_t Z, int regularized)
{
    acb_t t, u;
    mag_t M2, e, lim;
    slong k;

    acb_init(t);
    acb_init(u);
    mag_init(M2);
    mag_init(e);
    mag_init(lim);

    k = _acb_hypgeom_inexact_terms_bound(M, A, 1, B, 1, Z, 30000, NULL);

    /* exp(sup Re(z)) */
    arb_exp(acb_realref(t), acb_realref(Z), 40);
    arb_get_mag(e, acb_realref(t));

    /* (the second sum only as far as it can give a smaller bound, and
       with no more than a few times the terms of the first: a longer
       sum has larger terms) */
    if (mag_is_finite(M) && !mag_is_zero(e) && mag_is_finite(e))
        mag_div_lower(lim, M, e);
    else
        mag_inf(lim);

    acb_sub(t, B, A, 40);
    acb_neg(u, Z);
    _acb_hypgeom_inexact_terms_bound(M2, t, 1, B, 1, u,
        mag_is_finite(M) ? FLINT_MIN(30000, 4 * k + 100) : 30000, mag_is_finite(lim) ? lim : NULL);

    if (mag_is_finite(M2))
        mag_mul(M2, M2, e);

    mag_min(M, M, M2);

    if (regularized && mag_is_finite(M))
    {
        acb_rgamma(t, B, 40);
        acb_get_mag(e, t);
        mag_mul(M, M, e);
    }

    acb_clear(t);
    acb_clear(u);
    mag_clear(M2);
    mag_clear(e);
    mag_clear(lim);
}

/* The asymptotic series on balls */
static void
m_majorant_asymp(mag_t M, const acb_t A, const acb_t B, const acb_t Z, int regularized)
{
    acb_t t;
    acb_init(t);
    acb_hypgeom_m_asymp(t, A, B, Z, regularized, 64);
    acb_get_mag(M, t);
    acb_clear(t);
}

static void
m_eval_exact(acb_t res, acb_srcptr w, slong prec, void * ctx)
{
    acb_hypgeom_m(res, w, w + 1, w + 2, *((int *) ctx), prec);
}

/* The bounds for |1F1| on balls (see _acb_hypgeom_inexact) */
static void
m_majorant(mag_t M, acb_srcptr w, int stage, void * ctx)
{
    int regularized = *((int *) ctx);

    if (stage == 0)
    {
        m_majorant_terms(M, w, w + 1, w + 2, regularized);
    }
    else if (stage == 1)
    {
        m_majorant_euler(M, w, w + 1, w + 2, regularized);
    }
    else if (stage == 2)
    {
        mag_t t;
        mag_init(t);
        acb_get_mag_lower(t, w + 2);
        if (mag_cmp_2exp_si(t, 5) >= 0)
            m_majorant_asymp(M, w, w + 1, w + 2, regularized);
        else
            mag_inf(M);
        mag_clear(t);
    }
    else
    {
        /* the integrals along paths of steepest descent: close to the
           integral of the absolute value of the integrand */
        acb_t t;
        acb_init(t);
        acb_hypgeom_1f1_integration(t, w, w + 1, w + 2, regularized, 30);
        acb_get_mag(M, t);
        acb_clear(t);
    }
}

/* The function on inexact input from its value at the midpoints and a
   bound for the variation (inexact.c).  Returns 0 if there is no
   result, 2 if the enclosure is close to the range of the function. */
static int
m_inexact(acb_t res, const acb_t a, const acb_t b, const acb_t z, int regularized,
    slong prec, slong goal, double cancel, double exp_ratio, int lazy)
{
    acb_struct args[3];
    double R[3], sa, sb, sz, spar, na, nb, nz;
    int n, real;

    args[0] = *a;
    args[1] = *b;
    args[2] = *z;

    n = !acb_is_exact(a) + !acb_is_exact(b) + !acb_is_exact(z);

    na = arf_get_d(arb_midref(acb_realref(a)), ARF_RND_NEAR);
    sa = arf_get_d(arb_midref(acb_imagref(a)), ARF_RND_NEAR);
    nb = arf_get_d(arb_midref(acb_realref(b)), ARF_RND_NEAR);
    sb = arf_get_d(arb_midref(acb_imagref(b)), ARF_RND_NEAR);
    nz = arf_get_d(arb_midref(acb_realref(z)), ARF_RND_NEAR);
    sz = arf_get_d(arb_midref(acb_imagref(z)), ARF_RND_NEAR);

    if (!(hypot(na, sa) < 1e15) || !(hypot(nb, sb) < 1e15) || !(hypot(nz, sz) < 1e15))
        return 0;

    /* The radii of the polydisc: about the scales on which the function
       varies by a factor e.  The logarithmic derivatives in a and b are
       of the size of log |parameters|, and in z of size 1 or |a / b|. */
    spar = 2.0 + log(2.0 + hypot(na, sa) + hypot(nb, sb) + hypot(nz, sz));
    R[0] = R[1] = 1.0 / (n * spar);
    R[2] = 1.0 / (n * (1.0 + hypot(na, sa) / FLINT_MAX(hypot(nb, sb), 1.0)));

    /* (inside of the domains of the integrals, if the midpoint is) */
    if (nb - na > 0.0)
    {
        R[0] = FLINT_MIN(R[0], 0.05 * (nb - na));
        R[1] = FLINT_MIN(R[1], 0.05 * (nb - na));
    }
    if (na > 0.0)
        R[0] = FLINT_MIN(R[0], 0.05 * na);

    real = acb_is_real(a) && acb_is_real(b) && acb_is_real(z);

    return _acb_hypgeom_inexact(res, args, 3, R, m_eval_exact, m_majorant, 2,
                &regularized, real, prec, goal, cancel, exp_ratio, lazy, 1);
}

#define M_SERIES 0      /* and M_SERIES + 1 for the Kummer transformation */
#define M_ASYMP 2
#define M_INTEGRATION 3

static int
m_is_nonpositive_int(const acb_t x)
{
    return acb_is_int(x) && arb_is_nonpositive(acb_realref(x));
}

/* For integers 0 < a < b both series of the asymptotic expansion
   terminate, and with b terms in all

       1F1(a, b, z) = Gamma(b) / Gamma(b - a) (-z)^(-a)
                          sum_{k < b - a} (a)_k (a - b + 1)_k / (k! (-z)^k)
                    + Gamma(b) / Gamma(a) e^z z^(a - b)
                          sum_{k < a} (b - a)_k (1 - a)_k / (k! z^k)

   exactly (acb_hypgeom_m_asymp).  For small |z| the terms grow and the
   two parts cancel.  The estimated number of bits by which the largest
   term exceeds the function (which is about the larger part for large
   |z|, and at most e^max(0, Re z) in size), or a negative number if
   this is not the case, or if the sum is not shorter than the series,
   which has at least |z| + prec / 8 terms. */
static double
m_terminating_loss(const acb_t a, const acb_t b, m_cd A, m_cd B, m_cd Z, slong prec)
{
    double n1, n2, zabs, lz, t, r, h1, h2, p1, p2, lm, loss, k;

    if (!acb_is_int(a) || !acb_is_int(b) || !(A.re >= 1.0) || !(B.re > A.re) || !(B.re < 1e9))
        return -1.0;

    zabs = hypot(Z.re, Z.im);
    if (!(zabs > 1e-300) || !(B.re <= zabs + 0.125 * prec))
        return -1.0;

    lz = log(zabs);
    n1 = B.re - A.re;
    n2 = A.re;

    /* the largest terms of the two sums (the ratios decrease) */
    h1 = t = 0.0;
    for (k = 0.0; k < n1 - 1.0; k += 1.0)
    {
        r = (A.re + k) * (n1 - 1.0 - k) / ((k + 1.0) * zabs);
        if (r <= 1.0)
            break;
        t += log(r);
        h1 = t;
    }

    h2 = t = 0.0;
    for (k = 0.0; k < n2 - 1.0; k += 1.0)
    {
        r = (n1 + k) * (n2 - 1.0 - k) / ((k + 1.0) * zabs);
        if (r <= 1.0)
            break;
        t += log(r);
        h2 = t;
    }

    p1 = m_lgre(B) - _acb_hypgeom_d_lgamma_re(n1, 0.0, 0) - A.re * lz;
    p2 = m_lgre(B) - m_lgre(A) + Z.re - n1 * lz;
    lm = FLINT_MIN(FLINT_MAX(p1, p2), FLINT_MAX(0.0, Z.re));
    loss = (FLINT_MAX(p1 + h1, p2 + h2) - lm) * M_INVLOG2;
    loss = FLINT_MAX(loss, 0.0);

    if (loss > 2.0 * prec + 100.0 || B.re * (prec + loss) > (zabs + 0.125 * prec) * prec)
        return -1.0;

    return loss;
}

/* The Taylor coefficients in z for _acb_hypgeom_inexact_taylor:
   (a)_j / ((b)_j j!) 1F1(a + j, b + j, z), without (b)_j for the
   regularized function. */
typedef struct
{
    const acb_struct * a;
    const acb_struct * b;
    int regularized;
}
m_taylor_ctx;

#define M_TAYLOR_MAX_ORDER 60
#define M_TAYLOR_MAX_VALUES 8

static void
m_taylor_coeff(acb_t res, slong j, const acb_t z, slong prec, int crude, void * vctx)
{
    m_taylor_ctx * ctx = (m_taylor_ctx *) vctx;
    acb_t a, b, t, u;
    arb_t f;
    slong wp = prec + 10;

    acb_init(a);
    acb_init(b);
    acb_init(t);
    acb_init(u);
    arb_init(f);

    acb_rising_ui(t, ctx->a, j, wp);

    if (acb_is_zero(t))
    {
        /* (a polynomial of degree less than j) */
        acb_zero(res);
    }
    else
    {
        if (!ctx->regularized)
        {
            acb_rising_ui(u, ctx->b, j, wp);
            acb_div(t, t, u, wp);
        }

        arb_fac_ui(f, j, wp);
        acb_div_arb(t, t, f, wp);

        /* (exact, unless the parameters have more bits than that) */
        acb_add_ui(a, ctx->a, j, wp + 200);
        acb_add_ui(b, ctx->b, j, wp + 200);

        /* (crude: the series, or the one of Kummer's transformation,
           on the box) */
        if (!acb_is_exact(a) || !acb_is_exact(b))
            acb_indeterminate(res);
        else if (crude)
            acb_hypgeom_m_1f1(res, a, b, z, ctx->regularized, prec);
        else
            acb_hypgeom_m(res, a, b, z, ctx->regularized, prec);

        acb_mul(res, res, t, wp);
    }

    acb_clear(a);
    acb_clear(b);
    acb_clear(t);
    acb_clear(u);
    arb_clear(f);
}


/* c_{j+2} from c_j and c_{j+1} at m, for a solution of Kummer's
   equation z w'' + (b - z) w' - a w = 0:
   m (j+2)(j+1) c_{j+2} = (j + a) c_j - (j+1)(j + b - m) c_{j+1}. */
static void
m_taylor_recur(acb_t res, const acb_t c0, const acb_t c1, slong j, const acb_t m,
    slong prec, void * vctx)
{
    m_taylor_ctx * ctx = (m_taylor_ctx *) vctx;
    acb_t s, t;

    acb_init(s);
    acb_init(t);

    acb_add_si(s, ctx->a, j, prec);
    acb_mul(s, s, c0, prec);

    acb_add_si(t, ctx->b, j, prec);
    acb_sub(t, t, m, prec);
    acb_mul(t, t, c1, prec);
    acb_mul_si(t, t, j + 1, prec);
    acb_sub(s, s, t, prec);

    acb_mul_si(t, m, (j + 2) * (j + 1), prec);
    acb_div(res, s, t, prec);

    acb_clear(s);
    acb_clear(t);
}

/* Numerical integration, replacing best if more accurate (the complex
   integral if the real one does not apply, or is less accurate than
   enough). */
static void
m_integrate(acb_t best, slong * bestacc, const acb_t a, const acb_t b, const acb_t z,
    int regularized, slong iprec, slong enough)
{
    acb_t t;
    slong acc, racc = -WORD_MAX;

    acb_init(t);

    /* real arguments: on the real line (and for large parameters along
       paths of steepest descent) */
    if (acb_is_real(a) && acb_is_real(b) && acb_is_real(z))
    {
        arb_hypgeom_1f1_integration(acb_realref(t),
                acb_realref(a), acb_realref(b), acb_realref(z), regularized, iprec);
        arb_zero(acb_imagref(t));
        acc = racc = acb_rel_accuracy_bits(t);

        if (acc > *bestacc || (acb_is_finite(t) && !acb_is_finite(best)))
        {
            acb_swap(best, t);
            *bestacc = acc;
        }
    }

    /* complex arguments (or real ones outside of the domain of the real
       integral): along paths of steepest descent */
    if (racc < enough)
    {
        acb_hypgeom_1f1_integration(t, a, b, z, regularized, iprec);
        acc = acb_rel_accuracy_bits(t);

        if (acc > *bestacc || (acb_is_finite(t) && !acb_is_finite(best)))
        {
            acb_swap(best, t);
            *bestacc = acc;
        }
    }

    acb_clear(t);
}

void
acb_hypgeom_m(acb_t res, const acb_t a, const acb_t b, const acb_t z, int regularized, slong prec)
{
    acb_t best, t;
    m_cd A, B, Z, BA;
    double Lam[2], kpk[2], lm, lm2, pbits, c_asymp, n_asymp[2], asymp_wp;
    slong goal, accept, acc, bestacc, in_acc, in_min, guard;
    double budget, spent, cref = 1e300, spare;
    int real, exact, int_ok, tried_asymp, tried_int, iter, have_model, loop_first = 0, lazy_input;
    int moderate;
    int tried_taylor = 0, poly_z = 0;
    int series_runs[2];

    if (!acb_is_finite(a) || !acb_is_finite(b) || !acb_is_finite(z))
    {
        _acb_hypgeom_m_nointegration(res, a, b, z, regularized, prec);
        return;
    }

    acb_init(best);
    acb_init(t);

    real = acb_is_real(a) && acb_is_real(b) && acb_is_real(z);
    exact = acb_is_exact(a) && acb_is_exact(b) && acb_is_exact(z);
    lazy_input = !exact && !(_acb_hypgeom_inexact_is_coarse(a, prec)
        || _acb_hypgeom_inexact_is_coarse(b, prec) || _acb_hypgeom_inexact_is_coarse(z, prec));

    /* The accuracy that the input is good for: the relative condition
       number for an inexact z, a or b is about its absolute value. */
    in_min = _acb_hypgeom_est_input_accuracy(z);
    in_acc = _acb_hypgeom_est_input_accuracy(a);
    in_min = FLINT_MIN(in_min, in_acc);
    in_acc = _acb_hypgeom_est_input_accuracy(b);
    in_min = FLINT_MIN(in_min, in_acc);
    goal = FLINT_MIN(prec, in_min);
    goal = FLINT_MAX(goal, 2);

    /* The bits by which inexact input is more accurate than the
       precision: series which cancel by fewer bits are evaluated on
       the input with more precision, as for exact input. */
    spare = exact ? 1e300 : (double) in_min - (double) prec;

    /* (such input is not the result of rounding at this precision:
       there is no reason to expect a caller to provide more) */
    if (spare >= 30.0)
        lazy_input = 0;

    /* Nearly full accuracy, or what the input is good for. */
    accept = goal - (exact ? 10 : 12);
    bestacc = -WORD_MAX;
    acb_indeterminate(best);

    /* A polynomial in z (times e^z if b - a is the nonpositive integer)
       of which only the argument is inexact */
    if (!exact && acb_is_exact(a) && acb_is_exact(b))
    {
        acb_sub(t, b, a, prec);
        poly_z = m_is_nonpositive_int(a) || (acb_is_exact(t) && m_is_nonpositive_int(t));
    }

    A.re = arf_get_d(arb_midref(acb_realref(a)), ARF_RND_NEAR);
    A.im = arf_get_d(arb_midref(acb_imagref(a)), ARF_RND_NEAR);
    B.re = arf_get_d(arb_midref(acb_realref(b)), ARF_RND_NEAR);
    B.im = arf_get_d(arb_midref(acb_imagref(b)), ARF_RND_NEAR);
    Z.re = arf_get_d(arb_midref(acb_realref(z)), ARF_RND_NEAR);
    Z.im = arf_get_d(arb_midref(acb_imagref(z)), ARF_RND_NEAR);

    /* Integers 0 < a < b: the asymptotic expansion is a finite sum. */
    if (exact)
    {
        double loss = m_terminating_loss(a, b, A, B, Z, prec);

        if (loss >= 0.0)
        {
            acb_hypgeom_m_asymp(best, a, b, z, regularized,
                prec + (slong) loss + 12 + 2 * FLINT_BIT_COUNT((ulong) B.re));
            bestacc = acb_rel_accuracy_bits(best);

            if (bestacc >= accept)
            {
                acb_set_round(res, best, prec);
                acb_clear(best);
                acb_clear(t);
                return;
            }
        }
    }

    /* The estimates are not available: for huge numbers, and at the poles
       (where only a terminating series can make sense). */
    have_model = hypot(A.re, A.im) < 1e12 && hypot(B.re, B.im) < 1e12 && hypot(Z.re, Z.im) < 1e12
        && !(m_is_nonpositive_int(b) && !regularized);

    /* Parameters of moderate size: the choice by the size of z alone is
       usually right, and the estimates would cost as much as the
       function at low precision. */
    moderate = hypot(A.re, A.im) <= 30.0 && (hypot(B.re, B.im) <= 30.0 || B.re >= fabs(B.im));

    /* (Not with a large b for which that choice is the asymptotic
       expansion and its terms do not get small enough, as for
       1F1(1, 1000 + 1/2, 1000): the attempt would cost several times
       as much as the series.) */
    if (have_model && moderate && hypot(B.re, B.im) > 30.0)
    {
        int asymp, kummer;
        slong cwp;

        acb_hypgeom_m_choose(&asymp, &kummer, &cwp, a, b, z, regularized, prec);

        if (asymp)
        {
            m_cd Zn, BA0;
            double n0, more, cancel, need = (goal + 3.0) / M_INVLOG2;

            Zn.re = -Z.re; Zn.im = -Z.im;
            BA0.re = B.re - A.re; BA0.im = B.im - A.im;

            if (!_acb_hypgeom_est_u_asymp(&n0, &more, &cancel, A, B, Z, need, real)
                || !_acb_hypgeom_est_u_asymp(&n0, &more, &cancel, BA0, B, Zn, need, real))
                moderate = 0;
        }
    }

    if (!have_model || moderate)
    {
        _acb_hypgeom_m_nointegration(best, a, b, z, regularized, prec);
        bestacc = acb_rel_accuracy_bits(best);

        /* (for inexact input that is as accurate as the precision, a
           few bits less are accepted up to 1000 bits: the function at
           the midpoints with a bound for the variation costs several
           times as much there; not for the polynomials, which have
           their Taylor expansion) */
        if (bestacc >= accept - ((lazy_input && prec <= 1000 && !poly_z) ? 18 : 0) || !have_model)
        {
            /* (without estimates, integration if less than half of the
               bits are accurate) */
            if (bestacc < accept && bestacc < goal / 2 && goal > 2)
                m_integrate(best, &bestacc, a, b, z, regularized,
                    (goal < prec - 30) ? goal + 30 : prec, goal / 2);

            acb_set_round(res, best, prec);
            acb_clear(best);
            acb_clear(t);
            return;
        }
    }

    /* The series behind the regularized function at a pole of Gamma(b) */
    if (regularized && m_is_nonpositive_int(b))
    {
        A.re = A.re - B.re + 1.0; A.im = A.im - B.im;
        B.re = 2.0 - B.re; B.im = -B.im;
    }

    BA.re = B.re - A.re; BA.im = B.im - A.im;

    {
        m_cd Zn;
        Zn.re = -Z.re; Zn.im = -Z.im;
        Lam[0] = m_peak(kpk, A, B, Z);
        Lam[1] = m_peak(kpk + 1, BA, B, Zn);
    }

    m_logabs(&lm, &lm2, A, B, Z, FLINT_MIN(Lam[0], Z.re + Lam[1]) + 2.0);
    Lam[1] += Z.re;

    pbits = FLINT_MAX(m_arb_bits(acb_realref(a)), m_arb_bits(acb_imagref(a)));
    pbits = FLINT_MAX(pbits, m_arb_bits(acb_realref(b)));
    pbits = FLINT_MAX(pbits, m_arb_bits(acb_imagref(b)));
    pbits = FLINT_MAX(pbits, m_arb_bits(acb_realref(z)));
    pbits = FLINT_MAX(pbits, m_arb_bits(acb_imagref(z)));

    /* The asymptotic expansion: the two series
       (-z)^(-a) / Gamma(b-a) 2F0(a, a-b+1, -1/z) and
       exp(z) z^(a-b) / Gamma(a) 2F0(b-a, 1-a, 1/z) must reach exp(-goal)
       relative to the larger part, and the working precision must cover
       the cancellation if the terms grow before they decrease. */
    {
        m_cd Zn, lz, lzn;
        double w[2], wmax, zabs, base;
        int i;

        Zn.re = -Z.re; Zn.im = -Z.im;
        lz = m_cd_log(Z);
        lzn = m_cd_log(Zn);
        zabs = hypot(Z.re, Z.im);

        w[0] = -(A.re * lzn.re - A.im * lzn.im) - m_lgre(BA);
        w[1] = Z.re - (BA.re * lz.re - BA.im * lz.im) - m_lgre(A);
        wmax = FLINT_MAX(w[0], w[1]);

        base = FLINT_MIN(prec, goal + 20) + 15.0
            + 1.5 * log(1.0 + hypot(A.re, A.im) + hypot(B.re, B.im) + zabs) * M_INVLOG2;
        asymp_wp = base;
        c_asymp = (zabs < 4.0) ? 1e300 : 30.0;

        for (i = 0; i < 2 && c_asymp < 1e300; i++)
        {
            double more, cancel, need;

            need = (goal + 3.0) / M_INVLOG2 - (wmax - w[i]);

            if (!_acb_hypgeom_est_u_asymp(n_asymp + i, &more, &cancel, i ? BA : A, B, i ? Zn : Z, need, real))
            {
                c_asymp = 1e300;
                break;
            }

            /* (the number of terms, and the factor of the error bound
               that the truncation accounts for, are limited by the
               precision) */
            asymp_wp = FLINT_MAX(asymp_wp, n_asymp[i] / 9.5);
            asymp_wp = FLINT_MAX(asymp_wp, (need + more) * M_INVLOG2 / 10.5 + 4.0);

            /* (the cancellation matters as far as this part of the
               function does) */
            cancel -= wmax - w[i];
            if (cancel > 0.0)
                asymp_wp = FLINT_MAX(asymp_wp, base + cancel * M_INVLOG2);
        }

        if (asymp_wp > 1e8)
            c_asymp = 1e300;

        if (c_asymp < 1e300)
            c_asymp += (n_asymp[0] + n_asymp[1]) * 1.5 * m_series_term_cost(asymp_wp, pbits, real)
                + 4.0 * 7e-5 * pow(asymp_wp, 1.6);
    }

    /* What the attempts may cost in total: without a limit, the series
       with millions of terms or bits would be attempted for huge
       parameters.  (The rest is cut down after an attempt that gives
       no information.) */
    budget = 5e4 + 5.0 * m_integration_cost(FLINT_MIN(prec, goal + 30), real);
    spent = 0.0;
    tried_asymp = tried_int = 0;
    series_runs[0] = series_runs[1] = 0;

    /* The integrals: the Euler integral for Re(a) >= 1 and
       Re(b - a) >= 1 (int_ok = 1); two integrals for U where only one
       of the two holds, for z in the right half-plane (Re(b - a) >= 1)
       or for exp(z) 1F1(b - a, b, -z) and z in the left half-plane
       (Re(a) >= 1), if the function oscillates (int_ok = 2). */
    {
        m_cd A0, B0, BA0, Zn;
        int ok_a, ok_ba;

        A0.re = arf_get_d(arb_midref(acb_realref(a)), ARF_RND_NEAR);
        A0.im = arf_get_d(arb_midref(acb_imagref(a)), ARF_RND_NEAR);
        B0.re = arf_get_d(arb_midref(acb_realref(b)), ARF_RND_NEAR);
        B0.im = arf_get_d(arb_midref(acb_imagref(b)), ARF_RND_NEAR);
        BA0.re = B0.re - A0.re; BA0.im = B0.im - A0.im;
        Zn.re = -Z.re; Zn.im = -Z.im;

        acb_sub_ui(t, a, 1, prec);
        ok_a = arb_is_nonnegative(acb_realref(t));
        acb_sub(t, b, a, prec);
        acb_sub_ui(t, t, 1, prec);
        ok_ba = arb_is_nonnegative(acb_realref(t));

        int_ok = 0;
        if (ok_a && ok_ba)
            int_ok = 1;
        else if (have_model && ok_ba && arb_is_positive(acb_realref(z)))
        {
            if (_acb_hypgeom_est_2u_oscillatory(A0, B0, Z))
                int_ok = 2;
        }
        else if (have_model && ok_a && arb_is_negative(acb_realref(z)))
        {
            if (_acb_hypgeom_est_2u_oscillatory(BA0, B0, Zn))
                int_ok = 2;
        }
    }

    /* Inexact input: the error bounds of the series on balls are
       inflated by the cancellation, which more precision does not help.
       The function is evaluated at the midpoints, with a bound for its
       variation (see inexact.c for what is attempted at which cost).
       What follows after this is for the case that there is no result
       (a ball containing a zero of the function, or huge numbers). */
inexact_input:

    if (!exact)
    {
        double cmain, cint, cexp, cmin = 1e300;
        slong wp = FLINT_MIN(prec, goal + 30);
        int status, lazy, nothing, done = 0, s;

        /* the cost of the value: a series, the asymptotic series or
           the integral */
        cmain = c_asymp;
        for (s = 0; s < 2; s++)
        {
            m_cd Zn;
            double c, tw, T;

            Zn.re = -Z.re; Zn.im = -Z.im;
            c = (Lam[s] - lm) * M_INVLOG2;
            c = FLINT_MIN(FLINT_MAX(c, 0.0), 1e9);
            cmin = FLINT_MIN(cmin, c);
            tw = wp + c + 16.0 + (wp + c) / 64.0;
            T = (s == 0) ? m_num_terms(A, B, Z, kpk[0], Lam[0] - tw / M_INVLOG2)
                         : m_num_terms(BA, B, Zn, kpk[1], Lam[1] - Z.re - tw / M_INVLOG2);
            if (tw < 1e8 && T < 1e9)
                cmain = FLINT_MIN(cmain, T * m_series_term_cost(tw, pbits, real) + 7e-5 * pow(tw, 1.6));
        }

        cint = (int_ok == 0) ? 1e300
            : (int_ok == 1) ? m_integration_cost(wp, real)
            : (real ? 1.0 : 2.0) * m_integration_cost(wp, 0);

        /* (the reference for the cost of an integral below) */
        cref = cmain;

        /* At low precision, a series that loses little (or the
           asymptotic series): on the input, as for exact input, which
           saves the bound.  (If that is not accurate, this is reached
           once more.) */
        if (loop_first == 0 && wp <= 100 && (cmin <= 8.0 || c_asymp < 1e300) && cmain < cint)
        {
            loop_first = 1;
            goto methods;
        }

        /* Input which is more accurate than the precision by more
           bits than a series loses: on the input with more precision,
           as for exact input.  (If that is not accurate, this is
           reached once more.) */
        if (loop_first == 0 && spare >= 20.0 + ((c_asymp < 1e300) ? 0.0 : cmin) && cmain < cint)
        {
            loop_first = 1;
            goto methods;
        }

        /* A polynomial in z (times e^z if b - a is the nonpositive
           integer) of which only the argument is inexact: the Taylor
           expansion at its midpoint, whose terms are polynomials of
           the same kind with exact arguments. */
        if (!tried_taylor && poly_z && cmain < 1e300 && cmain <= 4.0 * cint)
        {
            tried_taylor = 1;

            {
                m_taylor_ctx ctx;
                slong order;

                ctx.a = a;
                ctx.b = b;
                ctx.regularized = regularized;

                if (_acb_hypgeom_inexact_taylor(t, &order, z, m_taylor_coeff,
                        m_taylor_recur, &ctx, real, prec, goal,
                        M_TAYLOR_MAX_ORDER, M_TAYLOR_MAX_VALUES))
                {
                    acc = acb_rel_accuracy_bits(t);

                    if (acc > bestacc || (acb_is_finite(t) && !acb_is_finite(best)))
                    {
                        acb_swap(best, t);
                        bestacc = acc;
                    }

                    acb_set_round(res, best, prec);
                    acb_clear(best);
                    acb_clear(t);
                    return;
                }
            }
        }

        /* the bound from the integrals on complex balls */
        cexp = 5.0 * m_integration_cost(64, 0);

        /* (less accurate than the precision: more precision does not
           improve the result) */
        lazy = lazy_input;

        /* The series cancel by more bits than the input has, and the
           bound is out of proportion: the series give no result that
           is worth its cost; the caller can increase the precision.
           (Not in the domain of the Euler integral, where there is a
           bound that costs little.) */
        nothing = lazy && cmin > goal + 11.0 && c_asymp >= 1e300
            && cexp > 2.0 * FLINT_MIN(cmain, cint)
            && !(arb_is_positive(acb_realref(a)) && int_ok == 1);

        /* Where the value would be computed by an integral: the
           integral on the input, which encloses the function without
           the loss of the series.  Where the series on the input
           give nothing (no result at all, or one from the midpoints
           with bounds that are computed whatever they cost), also an
           integral of up to four times the cost of the series at the
           precision at which they succeed. */
        if (int_ok != 0 && cint <= ((nothing || (!lazy && cmin > goal + 11.0)) ? 4.0 : 1.0) * cmain)
        {
            m_integrate(best, &bestacc, a, b, z, regularized,
                (goal < prec - 30) ? goal + 30 : prec, FLINT_MIN(accept, goal - 20));

            /* (the value and a bound would cost two more integrals) */
            done = (bestacc >= accept - 18);
        }

        if (nothing || done)
        {
            acb_set_round(res, best, prec);
            acb_clear(best);
            acb_clear(t);
            return;
        }

        cmain = FLINT_MIN(cmain, cint);

        status = m_inexact(t, a, b, z, regularized, prec, goal, cmin, cexp / cmain, lazy);

        if (status != 0)
        {
            acc = acb_rel_accuracy_bits(t);

            if (acc > bestacc || (acb_is_finite(t) && !acb_is_finite(best)))
            {
                acb_swap(best, t);
                bestacc = acc;
            }

            /* a bound that is not close to the function: the integral
               on the input, if that is not out of proportion */
            if (status == 1 && bestacc < accept && int_ok && (!lazy || cint <= 1.5 * cmain))
                m_integrate(best, &bestacc, a, b, z, regularized,
                    (goal < prec - 30) ? goal + 30 : prec, FLINT_MIN(accept, goal - 20));

            acb_set_round(res, best, prec);
            acb_clear(best);
            acb_clear(t);
            return;
        }
    }

methods:

    for (iter = 0; iter < 8 && bestacc < accept && loop_first != 2; iter++)
    {
        double cost[4], T[2], cbest;
        slong wps[2], cancel[2];
        int s, method;

        for (s = 0; s < 2; s++)
        {
            m_cd Zn;
            double c, base, tw;

            Zn.re = -Z.re; Zn.im = -Z.im;

            c = (Lam[s] - lm) * M_INVLOG2;
            c = FLINT_MIN(FLINT_MAX(c, 0.0), 1e9);
            cancel[s] = (slong) c;

            /* no more precision than an inexact input is good for */
            base = FLINT_MIN(prec, goal + 20);
            guard = 16 + (slong) ((base + c) / 64.0);
            tw = base + c + guard;

            T[s] = (s == 0) ? m_num_terms(A, B, Z, kpk[0], Lam[0] - tw / M_INVLOG2)
                            : m_num_terms(BA, B, Zn, kpk[1], Lam[1] - Z.re - tw / M_INVLOG2);

            if (tw > 1e8 || T[s] > 1e9 || series_runs[s] >= 4
                || (!exact && spare < c + 20.0 && c > ((int_ok && !tried_int)
                        ? FLINT_MAX(12.0, goal / 8.0) : goal / 2.0 + 10.0)))
            {
                cost[M_SERIES + s] = 1e300;
                wps[s] = 0;
            }
            else
            {
                wps[s] = (slong) tw;
                /* (and the exponential or gamma function) */
                cost[M_SERIES + s] = T[s] * m_series_term_cost(tw, pbits, real)
                    + 7e-5 * pow(tw, 1.6);
            }

            if (cost[M_SERIES + s] > budget - spent)
                cost[M_SERIES + s] = 1e300;
        }

        cost[M_ASYMP] = (tried_asymp >= 2 || c_asymp > budget - spent) ? 1e300 : c_asymp;
        cost[M_INTEGRATION] = (tried_int || !int_ok) ? 1e300
                    : (int_ok == 1) ? m_integration_cost(FLINT_MIN(prec, goal + 30), real)
                    : (real ? 1.0 : 2.0) * m_integration_cost(FLINT_MIN(prec, goal + 30), 0);

        /* (for inexact input that is as accurate as the precision:
           no integral at many times the cost of the other methods,
           whose results are limited by the input) */
        if (iter == 0 && !(cref < 1e300))
            cref = FLINT_MIN(FLINT_MIN(cost[0], cost[1]), cost[M_ASYMP]);
        if (lazy_input && cref < 1e300 && cost[M_INTEGRATION] > 4.0 * cref)
            cost[M_INTEGRATION] = 1e300;

        method = 0;
        cbest = cost[0];
        for (s = 1; s < 4; s++)
        {
            if (cost[s] < cbest)
            {
                cbest = cost[s];
                method = s;
            }
        }

        if (cbest >= 1e300)
            break;

        spent += cbest;

        if (method == M_ASYMP)
        {
            acb_hypgeom_m_asymp(t, a, b, z, regularized, (slong) asymp_wp);
            acc = acb_rel_accuracy_bits(t);
            tried_asymp++;

            /* a near miss: once more with the missing bits */
            if (acc < accept && acc > 4 && tried_asymp == 1 && exact)
            {
                asymp_wp += accept - acc + 10;
                c_asymp *= 1.0 + (accept - acc + 10) / asymp_wp;
            }
            else
                tried_asymp = 2;
        }
        else if (method == M_INTEGRATION)
        {
            /* (no more precision than the input is good for) */
            m_integrate(best, &bestacc, a, b, z, regularized,
                (goal < prec - 30) ? goal + 30 : prec, FLINT_MIN(accept, goal - 20));
            tried_int = 1;

            /* (the integral encloses the function on inexact input
               without the loss of the series) */
            if (!exact && acb_is_finite(best))
                break;

            continue;
        }
        else
        {
            slong macc;

            s = method - M_SERIES;
            _acb_hypgeom_m_1f1(t, a, b, z, regularized, wps[s], FLINT_MIN(wps[s], prec + 10), s,
                (slong) (1.25 * T[s] + 100.0));
            acc = macc = acb_rel_accuracy_bits(t);
            series_runs[s]++;

            if (!exact && spare < cancel[s] + 20.0)
            {
                /* (more precision does not improve the error bounds
                   for inexact parameters) */
                series_runs[s] = 4;
            }
            else if (macc < accept)
            {
                /* The cancellation was underestimated.  If the result
                   is not all noise, this shows by how much. */
                double loss;
                if (macc >= 4)
                {
                    loss = (double) cancel[s] + (accept - macc) + 10.0;
                    lm = FLINT_MIN(lm, Lam[s] - loss / M_INVLOG2);
                }
                else if (lm2 < lm - 20.0)
                {
                    /* the other saddle point */
                    lm = lm2;
                    lm2 = 1e300;
                    budget = spent + (budget - spent) / 3.0;
                }
                else
                {
                    /* twice the precision */
                    loss = cancel[s] + wps[s];
                    lm = FLINT_MIN(lm, Lam[s] - loss / M_INVLOG2);
                    budget = spent + (budget - spent) / 3.0;
                }
            }
            else if (acc < accept)
            {
                /* accurate at the midpoints, with error bounds for
                   inexact parameters that more precision does not
                   improve */
                series_runs[s] = 4;
            }
        }

        if (acc > bestacc || (acb_is_finite(t) && !acb_is_finite(best)))
        {
            acb_swap(best, t);
            bestacc = acc;
        }
    }

    if (loop_first == 1 && bestacc < accept)
    {
        loop_first = 2;
        goto inexact_input;
    }

    if (bestacc == -WORD_MAX)
    {
        _acb_hypgeom_m_nointegration(best, a, b, z, regularized, prec);
        bestacc = acb_rel_accuracy_bits(best);
    }

    acb_set_round(res, best, prec);

    acb_clear(best);
    acb_clear(t);
}

void
acb_hypgeom_1f1(acb_t res, const acb_t a, const acb_t b, const acb_t z, int regularized, slong prec)
{
    acb_hypgeom_m(res, a, b, z, regularized, prec);
}
