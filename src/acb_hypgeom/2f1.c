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
#include "acb_hypgeom.h"
#include "acb_hypgeom/impl.h"
#include "arb_hypgeom.h"

static void
_acb_hypgeom_2f1r_reduced(acb_t res,
    const acb_t b, const acb_t c, const acb_t z, slong prec)
{
    acb_t t, u;
    acb_init(t);
    acb_init(u);
    acb_sub_ui(t, z, 1, prec);
    acb_neg(t, t);
    acb_neg(u, b);
    acb_pow(t, t, u, prec);
    acb_rgamma(u, c, prec);
    acb_mul(t, t, u, prec);
    acb_set(res, t);
    acb_clear(t);
    acb_clear(u);
    return;
}

void
_acb_hypgeom_2f1_nointegration(acb_t res, const acb_t a, const acb_t b,
        const acb_t c, const acb_t z, int flags, slong prec)
{
    int algorithm, regularized;

    regularized = flags & ACB_HYPGEOM_2F1_REGULARIZED;

    if (!acb_is_finite(a) || !acb_is_finite(b) || !acb_is_finite(c) || !acb_is_finite(z))
    {
        acb_indeterminate(res);
        return;
    }

    if (acb_is_zero(z))
    {
        if (regularized)
            acb_rgamma(res, c, prec);
        else
            acb_one(res);
        return;
    }

    if (regularized && acb_is_int(c) && arb_is_nonpositive(acb_realref(c)))
    {
        if ((acb_is_int(a) && arb_is_nonpositive(acb_realref(a)) &&
            arf_cmp(arb_midref(acb_realref(a)), arb_midref(acb_realref(c))) >= 0) ||
            (acb_is_int(b) && arb_is_nonpositive(acb_realref(b)) &&
            arf_cmp(arb_midref(acb_realref(b)), arb_midref(acb_realref(c))) >= 0))
        {
            acb_zero(res);
            return;
        }
    }

    if (regularized && acb_eq(a, c))
    {
        _acb_hypgeom_2f1r_reduced(res, b, c, z, prec);
        return;
    }

    if (regularized && acb_eq(b, c))
    {
        _acb_hypgeom_2f1r_reduced(res, a, c, z, prec);
        return;
    }

    /* polynomial */
    if (acb_is_int(a) && arf_sgn(arb_midref(acb_realref(a))) <= 0 &&
         arf_cmpabs_ui(arb_midref(acb_realref(a)), prec) < 0)
    {
        acb_hypgeom_2f1_direct(res, a, b, c, z, regularized, prec);
        return;
    }

    /* polynomial */
    if (acb_is_int(b) && arf_sgn(arb_midref(acb_realref(b))) <= 0 &&
         arf_cmpabs_ui(arb_midref(acb_realref(b)), prec) < 0)
    {
        acb_hypgeom_2f1_direct(res, a, b, c, z, regularized, prec);
        return;
    }

    /* Try to reduce to a polynomial case using the Pfaff transformation */
    /* TODO: look at flags for integer c-b, c-a here, even when c is nonexact */
    if (acb_is_exact(c))
    {
        acb_t t;
        acb_init(t);

        acb_sub(t, c, b, prec);

        if (acb_is_int(t) && arb_is_nonpositive(acb_realref(t)))
        {
            acb_hypgeom_2f1_transform(res, a, b, c, z, flags, 1, prec);
            acb_clear(t);
            return;
        }

        acb_sub(t, c, a, prec);

        if (acb_is_int(t) && arb_is_nonpositive(acb_realref(t)))
        {
            int f1, f2;

            /* When swapping a, b, also swap the flags. */
            f1 = flags & ACB_HYPGEOM_2F1_AC;
            f2 = flags & ACB_HYPGEOM_2F1_BC;

            flags &= ~ACB_HYPGEOM_2F1_AC;
            flags &= ~ACB_HYPGEOM_2F1_BC;

            if (f1) flags |= ACB_HYPGEOM_2F1_BC;
            if (f2) flags |= ACB_HYPGEOM_2F1_AC;

            acb_hypgeom_2f1_transform(res, b, a, c, z, flags, 1, prec);
            acb_clear(t);
            return;
        }

        acb_clear(t);
    }

    /* special value at z = 1 */
    if (acb_is_one(z))
    {
        acb_t t, u, v;

        acb_init(t);
        acb_init(u);
        acb_init(v);

        acb_sub(t, c, a, prec);
        acb_sub(u, c, b, prec);
        acb_sub(v, t, b, prec);

        if (arb_is_positive(acb_realref(v)))
        {
            acb_rgamma(t, t, prec);
            acb_rgamma(u, u, prec);
            acb_mul(t, t, u, prec);
            acb_gamma(v, v, prec);
            acb_mul(t, t, v, prec);

            if (!regularized)
            {
                acb_gamma(v, c, prec);
                acb_mul(t, t, v, prec);
            }

            acb_set(res, t);
        }
        else
        {
            acb_indeterminate(res);
        }

        acb_clear(t);
        acb_clear(u);
        acb_clear(v);

        return;
    }

    algorithm = acb_hypgeom_2f1_choose(z);

    if (algorithm == 0)
    {
        acb_hypgeom_2f1_direct(res, a, b, c, z, regularized, prec);
    }
    else if (algorithm >= 1 && algorithm <= 5)
    {
        acb_hypgeom_2f1_transform(res, a, b, c, z, flags, algorithm, prec);
    }
    else
    {
        acb_hypgeom_2f1_corner(res, a, b, c, z, regularized, prec);
    }
}

/* ------------------------------------------------------------------------- */
/* Choice of algorithm from estimates of the costs                            */
/* ------------------------------------------------------------------------- */

/* See the comments for acb_hypgeom_m; the estimates are in estimates.c.
   Four series are compared for |z| < 1 or |z / (z - 1)| < 1:

   0:  2F1(a, b, c, z)
   1:  (1-z)^(c-a-b) 2F1(c-a, c-b, c, z)
   2:  (1-z)^(-a) 2F1(a, c-b, c, z/(z-1))
   3:  (1-z)^(-b) 2F1(c-a, b, c, z/(z-1))

   They differ in the cancellation and in the number of terms.  For z
   where a connection formula (at 1/z, 1 - z, ...) is chosen, that is
   tried first, and its precision is raised once by the accuracy that is
   found missing. */

/* (the largest order of the Taylor expansion in z for polynomials, and
   the largest number of values of the function and its derivatives at
   the midpoint) */
#define F_TAYLOR_MAX_ORDER 60
#define F_TAYLOR_MAX_VALUES 8

#define F_INVLOG2 1.4426950408889634

typedef _acb_hypgeom_cd f_cd;

static int
f_swap_flags(int flags)
{
    int f1, f2;

    f1 = flags & ACB_HYPGEOM_2F1_AC;
    f2 = flags & ACB_HYPGEOM_2F1_BC;
    flags &= ~ACB_HYPGEOM_2F1_AC;
    flags &= ~ACB_HYPGEOM_2F1_BC;
    if (f1) flags |= ACB_HYPGEOM_2F1_BC;
    if (f2) flags |= ACB_HYPGEOM_2F1_AC;
    return flags;
}

/* one of the four series at precision wp */
static void
f_series_form(acb_t res, int form, const acb_t a, const acb_t b, const acb_t c,
    const acb_t z, int flags, slong wp)
{
    int regularized = flags & ACB_HYPGEOM_2F1_REGULARIZED;

    if (form == 0)
    {
        acb_hypgeom_2f1_direct(res, a, b, c, z, regularized, wp);
    }
    else if (form == 1)
    {
        acb_t t, u, v;

        acb_init(t);
        acb_init(u);
        acb_init(v);

        acb_sub(t, c, a, wp);
        acb_sub(u, c, b, wp);
        acb_hypgeom_2f1_direct(res, t, u, c, z, regularized, wp);

        /* (1-z)^(c-a-b) */
        acb_sub(t, t, b, wp);
        acb_sub_ui(v, z, 1, wp);
        acb_neg(v, v);
        acb_pow(v, v, t, wp);
        acb_mul(res, res, v, wp);

        acb_clear(t);
        acb_clear(u);
        acb_clear(v);
    }
    else if (form == 2)
    {
        acb_hypgeom_2f1_transform(res, a, b, c, z, flags, 1, wp);
    }
    else
    {
        acb_hypgeom_2f1_transform(res, b, a, c, z, f_swap_flags(flags), 1, wp);
    }
}

/* ------------------------------------------------------------------------- */
/* Inexact input                                                              */
/* ------------------------------------------------------------------------- */

/* Bound for |2F1(a, b, c, z)| on balls with Re(B) > 0, Re(C - B) > 0 and
   Re(Z) < 1, from the Euler integral

       2F1 = Gamma(c) / (Gamma(b) Gamma(c - b))
                 int_0^1 t^(b-1) (1 - t)^(c-b-1) (1 - z t)^(-a) dt.

   With the real parts al, be, ga of a, b, c, and T distributed as
   Beta(be, ga - be),

       |2F1(a, b, c, z)| <= G E[|(1 - z T)^(-a)|],

   G = |Gamma(c)| Gamma(be) Gamma(ga - be) / (Gamma(ga) |Gamma(b)| |Gamma(c - b)|)
     <= exp(Im(b)^2 psi'(be) / 2 + Im(c - b)^2 psi'(ga - be) / 2)

   (see the bound for 1F1).  Here |(1 - z t)^(-a)| = |1 - z t|^(-al)
   exp(Im(a) arg(1 - z t)), |arg(1 - z t)| <= |arg(1 - z)| (the argument
   is monotone on a line segment), and

       |1 - z t|^(-al) <= (1 - x t)^(-al),    al >= 0,  x >= Re(z),
       |1 - z t|^(-al) <= (1 + r t)^(-al),    al <= 0,  r >= |z|,

   so that the expectation is at most 2F1(al, be, ga, x) or
   2F1(al, be, ga, -r).  These are increasing functions of T for
   x >= 0 and al >= 0 or for the second bound, and decreasing for x <= 0
   and al >= 0; the distribution of T increases stochastically with be
   and decreases with ga - be.  The largest value on the box of real
   parts is therefore at a corner in (be, ga) (at the largest be and the
   smallest ga, or the opposite), with the largest al (the smallest for
   x <= 0, the smallest negative one for the second bound).  This is
   close to |2F1| unless z or the parameters are far from real. */
static void
f_majorant_euler(mag_t M, const acb_t A, const acb_t B, const acb_t C, const acb_t Z)
{
    acb_t CB, t, al, be, ga, x;
    arb_t g, u, w;
    arf_t alo, ahi, x1, c;
    mag_t F1, F2, zm;
    slong prec = 40;
    int i;

    acb_init(CB); acb_init(t); acb_init(al); acb_init(be); acb_init(ga); acb_init(x);
    arb_init(g); arb_init(u); arb_init(w);
    arf_init(alo); arf_init(ahi); arf_init(x1); arf_init(c);
    mag_init(F1); mag_init(F2); mag_init(zm);

    acb_sub(CB, C, B, prec);
    arb_get_ubound_arf(x1, acb_realref(Z), prec);

    if (arb_is_positive(acb_realref(B)) && arb_is_positive(acb_realref(CB))
        && arf_cmp_si(x1, 1) < 0
        && acb_is_finite(A) && acb_is_finite(B) && acb_is_finite(C) && acb_is_finite(Z))
    {
        /* log G */
        arb_zero(g);
        for (i = 0; i < 2; i++)
        {
            const acb_struct * p = i ? CB : B;

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

        /* |Im(a)| sup |arg(1 - z)| */
        acb_sub_ui(t, Z, 1, prec);
        acb_neg(t, t);
        acb_arg(u, t, prec);
        arb_mul(u, u, acb_imagref(A), prec);
        arb_abs(u, u);
        arb_add(g, g, u, prec);
        arb_exp(g, g, prec);

        arb_get_lbound_arf(alo, acb_realref(A), prec);
        arb_get_ubound_arf(ahi, acb_realref(A), prec);

        mag_zero(F1);
        mag_zero(F2);

        /* members with al >= 0 */
        if (arf_sgn(ahi) >= 0)
        {
            acb_zero(x); arb_set_arf(acb_realref(x), x1);
            acb_zero(al); acb_zero(be); acb_zero(ga);

            if (arf_sgn(x1) >= 0)
            {
                arb_set_arf(acb_realref(al), ahi);
                arb_get_ubound_arf(c, acb_realref(B), prec);
                arb_set_arf(acb_realref(be), c);
                arb_get_lbound_arf(c, acb_realref(C), prec);
                arb_set_arf(acb_realref(ga), c);
            }
            else
            {
                if (arf_sgn(alo) > 0)
                    arb_set_arf(acb_realref(al), alo);
                arb_get_lbound_arf(c, acb_realref(B), prec);
                arb_set_arf(acb_realref(be), c);
                arb_get_ubound_arf(c, acb_realref(C), prec);
                arb_set_arf(acb_realref(ga), c);
            }

            acb_hypgeom_2f1(t, al, be, ga, x, 0, prec);
            acb_get_mag(F1, t);
        }

        /* members with al < 0 */
        if (arf_sgn(alo) < 0)
        {
            acb_get_mag(zm, Z);
            acb_zero(x);
            arf_set_mag(arb_midref(acb_realref(x)), zm);
            acb_neg(x, x);
            acb_zero(al); acb_zero(be); acb_zero(ga);
            arb_set_arf(acb_realref(al), alo);
            arb_get_ubound_arf(c, acb_realref(B), prec);
            arb_set_arf(acb_realref(be), c);
            arb_get_lbound_arf(c, acb_realref(C), prec);
            arb_set_arf(acb_realref(ga), c);

            acb_hypgeom_2f1(t, al, be, ga, x, 0, prec);
            acb_get_mag(F2, t);
        }

        mag_max(F1, F1, F2);
        arb_get_mag(M, g);
        mag_mul(M, M, F1);
    }
    else
    {
        mag_inf(M);
    }

    acb_clear(CB); acb_clear(t); acb_clear(al); acb_clear(be); acb_clear(ga); acb_clear(x);
    arb_clear(g); arb_clear(u); arb_clear(w);
    arf_clear(alo); arf_clear(ahi); arf_clear(x1); arf_clear(c);
    mag_clear(F1); mag_clear(F2); mag_clear(zm);
}

static void
f_eval_exact(acb_t res, acb_srcptr w, slong prec, void * ctx)
{
    acb_hypgeom_2f1(res, w, w + 1, w + 2, w + 3, *((int *) ctx), prec);
}

/* Bound for |2F1| on balls from the absolute values of the terms of
   the four series in z and in z / (z - 1), with their factors: close
   to |2F1| where the terms of one of them have equal signs (the domain
   of the Euler integral with real arguments), and too large by the
   cancellation otherwise. */
static void
f_majorant_terms(mag_t M, const acb_t A, const acb_t B, const acb_t C, const acb_t Z, int flags)
{
    acb_struct ab[2];
    acb_t CA, CB, W, one_z, t, e;
    mag_t M2, f, lim;
    slong k, nmax = 30000;
    int i;

    acb_init(ab); acb_init(ab + 1);
    acb_init(CA); acb_init(CB); acb_init(W); acb_init(one_z); acb_init(t); acb_init(e);
    mag_init(M2); mag_init(f); mag_init(lim);

    acb_sub(CA, C, A, 40);
    acb_sub(CB, C, B, 40);
    acb_sub_ui(one_z, Z, 1, 40);
    acb_div(W, Z, one_z, 40);
    acb_neg(one_z, one_z);

    mag_inf(M);

    for (i = 0; i < 4; i++)
    {
        acb_set(ab, (i == 0 || i == 2) ? A : CA);
        acb_set(ab + 1, (i == 0 || i == 3) ? B : CB);

        /* (1 - z)^(c - a - b), (1 - z)^(-a), (1 - z)^(-b) */
        mag_one(f);
        if (i != 0)
        {
            if (i == 1)
                acb_sub(e, CA, B, 40);
            else if (i == 2)
                acb_neg(e, A);
            else
                acb_neg(e, B);

            acb_pow(t, one_z, e, 40);
            acb_get_mag(f, t);

            if (!mag_is_finite(f))
                continue;
        }

        /* (the sum only as far as it can give a smaller bound, and
           with no more than a few times the terms of the shortest sum
           so far: a longer one has larger terms) */
        if (mag_is_finite(M) && !mag_is_zero(f))
        {
            mag_div_lower(lim, M, f);
            k = _acb_hypgeom_inexact_terms_bound(M2, ab, 2, C, 1, (i < 2) ? Z : W, nmax, lim);
        }
        else
            k = _acb_hypgeom_inexact_terms_bound(M2, ab, 2, C, 1, (i < 2) ? Z : W, nmax, NULL);

        if (!mag_is_finite(M2))
            continue;

        nmax = FLINT_MIN(nmax, 4 * k + 100);

        mag_mul(M2, M2, f);
        mag_min(M, M, M2);
    }

    if ((flags & ACB_HYPGEOM_2F1_REGULARIZED) && mag_is_finite(M))
    {
        acb_rgamma(t, C, 40);
        acb_get_mag(f, t);
        mag_mul(M, M, f);
    }

    acb_clear(ab); acb_clear(ab + 1);
    acb_clear(CA); acb_clear(CB); acb_clear(W); acb_clear(one_z); acb_clear(t); acb_clear(e);
    mag_clear(M2); mag_clear(f); mag_clear(lim);
}

/* The bounds for |2F1| on balls (see _acb_hypgeom_inexact) */
static void
f_majorant(mag_t M, acb_srcptr w, int stage, void * ctx)
{
    int flags = *((int *) ctx);
    acb_t t;

    acb_init(t);

    if (stage == 0)
    {
        f_majorant_terms(M, w, w + 1, w + 2, w + 3, flags);
    }
    else if (stage == 1)
    {
        /* from the integral in terms of b, or of a */
        mag_t M2;
        mag_init(M2);

        f_majorant_euler(M, w, w + 1, w + 2, w + 3);
        f_majorant_euler(M2, w + 1, w, w + 2, w + 3);
        mag_min(M, M, M2);

        if ((flags & ACB_HYPGEOM_2F1_REGULARIZED) && mag_is_finite(M))
        {
            acb_rgamma(t, w + 2, 40);
            acb_get_mag(M2, t);
            mag_mul(M, M, M2);
        }

        mag_clear(M2);
    }
    else if (stage == 2)
    {
        /* the algorithm chosen by z on the balls */
        _acb_hypgeom_2f1_nointegration(t, w, w + 1, w + 2, w + 3, flags, 64);
        acb_get_mag(M, t);
    }
    else
    {
        /* the integral along a path of steepest descent */
        _acb_hypgeom_2f1_integration_paths(t, w, w + 1, w + 2, w + 3, flags, 30);
        acb_get_mag(M, t);
    }

    acb_clear(t);
}

/* The function on inexact input from its value at the midpoints and a
   bound for the variation (inexact.c).  Returns 0 if there is no
   result, 2 if the enclosure is close to the range of the function. */
static int
f_inexact(acb_t res, const acb_t a, const acb_t b, const acb_t c, const acb_t z, int flags,
    slong prec, slong goal, double cancel, double exp_ratio, int lazy)
{
    acb_struct args[4];
    double R[4], sa, sb, sc, sz, spar, na, nb, nc, nz, d1, dcut, q;
    int n, real;

    args[0] = *a;
    args[1] = *b;
    args[2] = *c;
    args[3] = *z;

    n = !acb_is_exact(a) + !acb_is_exact(b) + !acb_is_exact(c) + !acb_is_exact(z);

    na = arf_get_d(arb_midref(acb_realref(a)), ARF_RND_NEAR);
    sa = arf_get_d(arb_midref(acb_imagref(a)), ARF_RND_NEAR);
    nb = arf_get_d(arb_midref(acb_realref(b)), ARF_RND_NEAR);
    sb = arf_get_d(arb_midref(acb_imagref(b)), ARF_RND_NEAR);
    nc = arf_get_d(arb_midref(acb_realref(c)), ARF_RND_NEAR);
    sc = arf_get_d(arb_midref(acb_imagref(c)), ARF_RND_NEAR);
    nz = arf_get_d(arb_midref(acb_realref(z)), ARF_RND_NEAR);
    sz = arf_get_d(arb_midref(acb_imagref(z)), ARF_RND_NEAR);
    d1 = hypot(1.0 - nz, sz);

    if (!(hypot(na, sa) < 1e15) || !(hypot(nb, sb) < 1e15) || !(hypot(nc, sc) < 1e15)
        || !(hypot(nz, sz) < 1e15) || !(d1 > 1e-15))
        return 0;

    /* The radii of the polydisc: about the scales on which the function
       varies by a factor e.  The logarithmic derivatives in the
       parameters are of the size of log |parameters| or log |1 - z|,
       and in z of the size min(|a|, |b|, |a b / c|) / min(|1 - z|, 1).
       The disc for z must not meet the cut. */
    spar = 2.0 + log(2.0 + hypot(na, sa) + hypot(nb, sb) + hypot(nc, sc)) + fabs(log(d1));
    R[0] = R[1] = R[2] = 1.0 / (n * spar);
    q = fmin(hypot(na, sa), hypot(nb, sb));
    q = fmin(q, hypot(na, sa) * hypot(nb, sb) / fmax(hypot(nc, sc), 1.0));
    R[3] = FLINT_MIN(d1, 1.0) / (n * (1.0 + q));
    dcut = (nz > 1.0) ? fabs(sz) : d1;
    R[3] = FLINT_MIN(R[3], 0.05 * dcut);

    /* (inside of the domains of the integrals, if the midpoint is) */
    if (nb > 0.0)
        R[1] = FLINT_MIN(R[1], 0.05 * nb);
    if (na > 0.0)
        R[0] = FLINT_MIN(R[0], 0.05 * na);
    q = nc - nb;
    if (q > 0.0)
    {
        R[1] = FLINT_MIN(R[1], 0.025 * q);
        R[2] = FLINT_MIN(R[2], 0.025 * q);
    }
    q = nc - na;
    if (q > 0.0)
    {
        R[0] = FLINT_MIN(R[0], 0.025 * q);
        R[2] = FLINT_MIN(R[2], 0.025 * q);
    }

    /* (real for z < 1) */
    {
        arb_t one;
        arb_init(one);
        arb_one(one);
        real = acb_is_real(a) && acb_is_real(b) && acb_is_real(c) && acb_is_real(z)
            && arb_lt(acb_realref(z), one);
        arb_clear(one);
    }

    return _acb_hypgeom_inexact(res, args, 4, R, f_eval_exact, f_majorant, 2,
                &flags, real, prec, goal, cancel, exp_ratio, lazy, 1);
}

/* The Taylor coefficients in z for _acb_hypgeom_inexact_taylor:
   (a)_j (b)_j / ((c)_j j!) 2F1(a + j, b + j, c + j, z), without (c)_j
   for the regularized function. */
typedef struct
{
    const acb_struct * a;
    const acb_struct * b;
    const acb_struct * c;
    int flags;
}
f_taylor_ctx;

static void
f_taylor_coeff(acb_t res, slong j, const acb_t z, slong prec, int crude, void * vctx)
{
    f_taylor_ctx * ctx = (f_taylor_ctx *) vctx;
    acb_t a, b, c, t, u;
    arb_t f;
    slong wp = prec + 10;

    acb_init(a);
    acb_init(b);
    acb_init(c);
    acb_init(t);
    acb_init(u);
    arb_init(f);

    acb_rising_ui(t, ctx->a, j, wp);
    acb_rising_ui(u, ctx->b, j, wp);
    acb_mul(t, t, u, wp);

    if (acb_is_zero(t))
    {
        /* (a polynomial of degree less than j) */
        acb_zero(res);
    }
    else
    {
        if (!(ctx->flags & ACB_HYPGEOM_2F1_REGULARIZED))
        {
            acb_rising_ui(u, ctx->c, j, wp);
            acb_div(t, t, u, wp);
        }

        arb_fac_ui(f, j, wp);
        acb_div_arb(t, t, f, wp);

        /* (exact, unless the parameters have more bits than that) */
        acb_add_ui(a, ctx->a, j, wp + 200);
        acb_add_ui(b, ctx->b, j, wp + 200);
        acb_add_ui(c, ctx->c, j, wp + 200);

        if (!acb_is_exact(a) || !acb_is_exact(b) || !acb_is_exact(c))
            acb_indeterminate(res);
        else if (crude)
            _acb_hypgeom_2f1_nointegration(res, a, b, c, z, ctx->flags, prec);
        else
            acb_hypgeom_2f1(res, a, b, c, z, ctx->flags, prec);

        acb_mul(res, res, t, wp);
    }

    acb_clear(a);
    acb_clear(b);
    acb_clear(c);
    acb_clear(t);
    acb_clear(u);
    arb_clear(f);
}

/* c_{j+2} from c_j and c_{j+1} at m, for a solution of the
   hypergeometric equation z (1 - z) w'' + (c - (a + b + 1) z) w' - a b w = 0:
   m (1 - m) (j+2)(j+1) c_{j+2} = (j + a)(j + b) c_j
       - (j+1) ((1 - 2m) j + c - (a + b + 1) m) c_{j+1}. */
static void
f_taylor_recur(acb_t res, const acb_t c0, const acb_t c1, slong j, const acb_t m,
    slong prec, void * vctx)
{
    f_taylor_ctx * ctx = (f_taylor_ctx *) vctx;
    acb_t s, t, u;

    acb_init(s);
    acb_init(t);
    acb_init(u);

    acb_add_si(s, ctx->a, j, prec);
    acb_add_si(t, ctx->b, j, prec);
    acb_mul(s, s, t, prec);
    acb_mul(s, s, c0, prec);

    /* (1 - 2m) j + c - (a + b + 1) m */
    acb_mul_2exp_si(t, m, 1);
    acb_sub_ui(t, t, 1, prec);
    acb_neg(t, t);
    acb_mul_si(t, t, j, prec);
    acb_add(t, t, ctx->c, prec);
    acb_add(u, ctx->a, ctx->b, prec);
    acb_add_ui(u, u, 1, prec);
    acb_submul(t, u, m, prec);

    acb_mul(t, t, c1, prec);
    acb_mul_si(t, t, j + 1, prec);
    acb_sub(s, s, t, prec);

    acb_sub_ui(t, m, 1, prec);
    acb_neg(t, t);
    acb_mul(t, t, m, prec);
    acb_mul_si(t, t, (j + 2) * (j + 1), prec);
    acb_div(res, s, t, prec);

    acb_clear(s);
    acb_clear(t);
    acb_clear(u);
}

/* Numerical integration, replacing best if more accurate (the complex
   integral if the real one does not apply, or is less accurate than
   enough). */
static void
f_integrate(acb_t best, slong * bestacc, const acb_t a, const acb_t b, const acb_t c,
    const acb_t z, int flags, slong iprec, slong enough)
{
    acb_t t;
    slong acc, racc = -WORD_MAX;

    acb_init(t);

    /* real arguments: on the real line (and for large parameters along
       paths of steepest descent) */
    if (acb_is_real(a) && acb_is_real(b) && acb_is_real(c) && acb_is_real(z))
    {
        arb_hypgeom_2f1_integration(acb_realref(t),
                acb_realref(a), acb_realref(b), acb_realref(c), acb_realref(z), flags, iprec);
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
        _acb_hypgeom_2f1_integration_paths(t, a, b, c, z, flags, iprec);
        acc = acb_rel_accuracy_bits(t);

        if (acc > *bestacc || (acb_is_finite(t) && !acb_is_finite(best)))
        {
            acb_swap(best, t);
            *bestacc = acc;
        }
    }

    acb_clear(t);
}

static int
f_is_nonpositive_int(const acb_t x)
{
    return acb_is_int(x) && arf_sgn(arb_midref(acb_realref(x))) <= 0;
}

/* Whether a function which is a polynomial for exact parameters (a, b,
   c - a or c - b a nonpositive integer) is analytic on the box z: a
   polynomial in z, or (1 - z)^(c - a - b) times one, which is rational
   if c - a - b is an integer and has the cut [1, inf) otherwise. */
static int
f_poly_analytic(const acb_t a, const acb_t b, const acb_t c, const acb_t z, slong prec)
{
    acb_t t;
    int res;

    if (f_is_nonpositive_int(a) || f_is_nonpositive_int(b))
        return 1;

    acb_init(t);

    acb_sub_ui(t, z, 1, prec + 64);

    if (acb_contains_zero(t))
    {
        res = 0;
    }
    else if (arb_is_negative(acb_realref(t)) || !arb_contains_zero(acb_imagref(t)))
    {
        res = 1;
    }
    else
    {
        /* (an integer which is computed exactly) */
        acb_sub(t, c, a, prec + 200);
        acb_sub(t, t, b, prec + 200);
        res = acb_is_int(t);
    }

    acb_clear(t);
    return res;
}

/* The integrals that apply: 1 for the Euler integral (Re(b) >= 1 and
   Re(c - b) >= 1, or the same for a), 2 for the two integrals past the
   point 1 / (1 - z) (Re(z) < 0, Re(b) >= 1 and Re(a - c) >= 0, or the
   same with a and b exchanged, where the function oscillates), 3 for
   the loop around 1 (|z - 1/2| < 1/2, Re(b) >= 1 and Re(a) <= 0, or
   the same with a and b exchanged, where the function oscillates: the
   two integrals for b - c + 1, b, a + b - c + 1 and 1 - 1 / z), 4 for
   the Euler integral past the point 1 / z, for a real z > 1 (with an
   imaginary part which is zero exactly), 0 for none. */
static int
f_integral_domain(const acb_t a, const acb_t b, const acb_t c, const acb_t z,
    f_cd A, f_cd B, f_cd C, f_cd Z, slong prec)
{
    acb_t t;
    int i, ok, cut, res = 0;

    acb_init(t);

    /* (the segment [0, 1] is not a path for z on the cut [1, inf): see 4) */
    acb_sub_ui(t, z, 1, prec);
    cut = arb_contains_zero(acb_imagref(z)) && !arb_is_negative(acb_realref(t));

    for (i = 0; i < 2 && res == 0 && !cut; i++)
    {
        const acb_struct * p = i ? a : b;

        acb_sub_ui(t, p, 1, prec);
        ok = arb_is_nonnegative(acb_realref(t));
        acb_sub(t, c, p, prec);
        acb_sub_ui(t, t, 1, prec);
        if (ok && arb_is_nonnegative(acb_realref(t)))
            res = 1;
    }

    for (i = 0; i < 2 && res == 0 && arb_is_negative(acb_realref(z))
        && _acb_hypgeom_est_2f1_left_oscillatory(A, B, C, Z); i++)
    {
        const acb_struct * p = i ? a : b;
        const acb_struct * q = i ? b : a;

        acb_sub_ui(t, p, 1, prec);
        ok = arb_is_nonnegative(acb_realref(t));
        acb_sub(t, q, c, prec);
        if (ok && arb_is_nonnegative(acb_realref(t)))
            res = 2;
    }

    if (res == 0 && !acb_contains_zero(z))
    {
        /* (Re(1 / z) > 1) */
        acb_inv(t, z, prec);
        acb_sub_ui(t, t, 1, prec);

        if (arb_is_positive(acb_realref(t)))
        {
            f_cd A2, C2, Z2;

            Z2.re = -arf_get_d(arb_midref(acb_realref(t)), ARF_RND_NEAR);
            Z2.im = -arf_get_d(arb_midref(acb_imagref(t)), ARF_RND_NEAR);

            for (i = 0; i < 2 && res == 0; i++)
            {
                const acb_struct * p = i ? a : b;
                const acb_struct * q = i ? b : a;
                f_cd P = i ? A : B, Q = i ? B : A;

                A2.re = P.re - C.re + 1.0; A2.im = P.im - C.im;
                C2.re = Q.re + A2.re; C2.im = Q.im + A2.im;

                acb_sub_ui(t, p, 1, prec);
                if (arb_is_nonnegative(acb_realref(t)) && arb_is_nonpositive(acb_realref(q))
                    && _acb_hypgeom_est_2f1_left_oscillatory(A2, P, C2, Z2))
                    res = 3;
            }
        }
    }

    if (res == 0 && cut && arb_is_zero(acb_imagref(z)))
    {
        acb_sub_ui(t, z, 1, prec);

        for (i = 0; i < 2 && res == 0 && arb_is_positive(acb_realref(t)); i++)
        {
            const acb_struct * p = i ? a : b;
            acb_t u;

            acb_init(u);
            acb_sub_ui(u, p, 1, prec);
            ok = arb_is_nonnegative(acb_realref(u));
            acb_sub(u, c, p, prec);
            acb_sub_ui(u, u, 1, prec);
            if (ok && arb_is_nonnegative(acb_realref(u)))
                res = 4;
            acb_clear(u);
        }
    }

    acb_clear(t);
    return res;
}

/* (the cost of the integrals, in microseconds) */
static double
f_integral_cost(int domain, slong prec, int real)
{
    if (domain == 1)
        return 1.3 * _acb_hypgeom_est_integration_cost(prec, real);
    if (domain == 2 || domain == 3)
        return (real ? 1.3 : 2.6) * _acb_hypgeom_est_integration_cost(prec, 0);
    if (domain == 4)
        return 1.3 * _acb_hypgeom_est_integration_cost(prec, 0);
    return 1e300;
}

/* The number of terms of the shortest of the series that terminate (a,
   b, c - a or c - b a nonpositive integer), or 0. */
static double
f_poly_terms(const acb_t a, const acb_t b, const acb_t c, slong prec)
{
    double n = 0.0, m;
    acb_t t;
    int i;

    acb_init(t);

    for (i = 0; i < 4; i++)
    {
        if (i == 0)
            acb_set(t, a);
        else if (i == 1)
            acb_set(t, b);
        else if (!acb_is_exact(c))
            break;
        else
            acb_sub(t, c, (i == 2) ? a : b, prec);

        if (f_is_nonpositive_int(t))
        {
            m = 1.0 - arf_get_d(arb_midref(acb_realref(t)), ARF_RND_UP);
            if (n == 0.0 || m < n)
                n = m;
        }
    }

    acb_clear(t);
    return n;
}

/* A polynomial with exact parameters and argument whose terms cancel
   by more bits than the precision provides: the same computation with
   the missing bits, or with twice the precision when the result is all
   noise, while the cost of nterms terms stays within the budget (in
   microseconds; see f_poly_budget).  The precision is limited as well,
   since a polynomial which vanishes at z is not accurate at any
   precision (and by maxwp, for inexact input which is that accurate). */
static double
f_poly_cost(const acb_t a, const acb_t b, const acb_t c, const acb_t z, double wp,
    double nterms, int real)
{
    double pbits, zterm;

    pbits = acb_bits(a);
    pbits = fmax(pbits, acb_bits(b));
    pbits = fmax(pbits, acb_bits(c));

    /* (an argument of many bits is multiplied in by rectangular
       splitting, at about 1.3e-4 wp microseconds per term for real
       data: 2F1(-10^4, 10^4 + 1; 1; z) with z of 3363 or 33333 bits
       takes 9 to 10 microseconds per term at 66666 bits and 17 to 18
       at 133333; the cost of a term with parameters of that many bits,
       which are multiplied in at every term, grows faster) */
    zterm = (acb_bits(z) > 64) ? (real ? 1.3e-4 : 4e-4) * wp : 0.0;

    return nterms * 1.6 * (_acb_hypgeom_est_series_term_cost(wp, pbits, real) + zterm)
        + 20.0 + 0.004 * pow(wp, 1.5);
}

/* The budget for the polynomial with more precision: 50 ms, or half
   the estimated cost of the integral cint which follows if the
   polynomial fails, and where there is none, 16 times the cost of the
   polynomial at the precision prec.  (At high precision an integral
   costs seconds to minutes, and the polynomial with all the bits that
   it loses is usually much cheaper: 2F1(-10^4, 10^4 + 1; 1; 1/2) at
   33333 bits by the series at 53000 bits in 7 ms.) */
static double
f_poly_budget(double cint, double cpoly)
{
    double budget = (cint < 1e300) ? 0.5 * cint : 16.0 * cpoly;

    return FLINT_MAX(budget, 5e4);
}

static void
f_poly_refine(acb_t best, slong * bestacc, const acb_t a, const acb_t b, const acb_t c,
    const acb_t z, int flags, slong prec, slong accept, double nterms, int real, double maxwp,
    double budget)
{
    acb_t t;
    slong wp = prec, wp2, acc = *bestacc;
    double limit, spent = 0.0, cost;
    int iter;

    acb_init(t);

    limit = 16.0 * prec + 8.0 * nterms + 1000.0;
    limit = FLINT_MIN(limit, maxwp);

    for (iter = 0; iter < 30 && *bestacc < accept; iter++)
    {
        wp2 = (acc >= 4) ? wp + (accept - acc) + 10 + wp / 64 : 2 * wp;

        cost = f_poly_cost(a, b, c, z, wp2, nterms, real);

        if (wp2 > limit || spent + cost > budget)
            break;

        spent += cost;
        wp = wp2;
        _acb_hypgeom_2f1_nointegration(t, a, b, c, z, flags, wp);
        acc = acb_rel_accuracy_bits(t);

        if (acc > *bestacc || (acb_is_finite(t) && !acb_is_finite(best)))
        {
            acb_swap(best, t);
            *bestacc = acc;
        }
    }

    acb_clear(t);
}

/* ------------------------------------------------------------------------- */
/* The candidates of the model                                                */
/* ------------------------------------------------------------------------- */

/* 0, ..., 3: the series at z and at z / (z - 1), each in two forms;
   4, ..., 7: the connection formulas at 1/z, 1/(1 - z), 1 - z and
   1 - 1/z (acb_hypgeom_2f1_transform with which = 2, ..., 5); 8: the
   last one with a and b exchanged, which has other series; 9: the
   corner algorithm (two series at |z| = 0.73 and two steps of
   analytic continuation), for z near exp(+-pi i/3). */
#define F_NC 10
#define F_MAX_LOG_ARGUMENT (-0.01)
#define F_CONN 4
#define F_CORNER 9
/* (the algorithm chosen by z alone once more, and the integrals) */
#define F_RETRY F_NC
#define F_INTEGRAL (F_NC + 1)

/* log |z1| for the first point of the corner algorithm, 0.375 + 0.625i */
#define F_CORNER_LZ (-0.3162615)

/* a connection formula chosen by z without estimates: the largest
   parameter, and the guard bits */
#define F_EASY_MAX_PARAMETER 12.0
#define F_EASY_GUARD 10

typedef struct
{
    f_cd fa[4];
    f_cd fb[4];
    f_cd C;
    double flz[4];
    double pk[4];
    double kpk[4];
    _acb_hypgeom_est_2f1_conn_struct E[5];
    double corner_pk;
    double corner_kpk;
    double pbits;       /* bits of the parameters and of z */
    double wpbits;      /* bits of the parameters and of z / (z - 1) */
    double parbits;     /* bits of the parameters */
    double maxpar;      /* the largest absolute value of a parameter */
    int real;           /* real parameters and argument */
    int cpar;           /* a parameter which is not real */
    int intpar;         /* a, b and c are integers */
    int moderate;       /* parameters below 256: costs per term from timings there */
}
f_model;

/* The number of terms of the candidate i with its terms down to
   2^-tw of the largest one; in tper, the number for one series (10^15
   if that is more than kmax, as far as this is found out early). */
static double
f_model_terms(double * tper, const f_model * M, int i, double tw, double kmax)
{
    double T;

    if (i < F_CONN)
    {
        T = _acb_hypgeom_est_2f1_terms_max(M->fa[i], M->fb[i], M->C, M->flz[i], M->kpk[i],
                M->pk[i] - tw / F_INVLOG2, kmax);
        *tper = T;
        return T;
    }

    if (i == F_CORNER)
    {
        T = _acb_hypgeom_est_2f1_terms_max(M->fa[0], M->fb[0], M->C, F_CORNER_LZ, M->corner_kpk,
                M->corner_pk - tw / F_INVLOG2, kmax);
        *tper = T;
        return T;
    }

    return _acb_hypgeom_est_2f1_conn_terms(tper, M->E + (i - F_CONN),
        M->E[i - F_CONN].Lam - tw / F_INVLOG2, kmax);
}

/* Bits which the candidate loses besides the cancellation of its
   terms, from measurements: the limit cases of the connection formulas
   (the derivative with respect to a parameter, with a tail bound which
   is less tight) about 0.6 bits per unit of the largest parameter
   above 12, up to 60 (beyond the parameters that were measured, a
   result which is not accurate shows what is missing), and 2 bits per
   unit if a, b and c are integers (85 bits for 2F1(32, 14; 48; z));
   the other connection formulas with complex parameters 0.03
   bits per term; the corner algorithm, in its continuation steps, 2.5
   bits per unit of the largest parameter (10 to 20 bits for
   parameters below 5, up to 60 below 20, and more than 150 at 60). */
static double
f_model_extra(const f_model * M, int i, double tper)
{
    double x = 0.0;

    if (i == F_CORNER)
    {
        x = 2.5 * M->maxpar;
    }
    else if (i >= F_CONN)
    {
        if (M->E[i - F_CONN].limit)
        {
            x = 0.6 * FLINT_MAX(M->maxpar - 12.0, 0.0);
            x = FLINT_MIN(x, 60.0);

            if (M->intpar)
            {
                x = 2.0 * FLINT_MAX(M->maxpar - 12.0, 0.0);
                x = FLINT_MIN(x, 150.0);
            }
        }
        else if (M->cpar)
        {
            x = FLINT_MIN(0.03 * tper, 60.0);
        }
    }

    return x;
}

/* The cost in microseconds of the candidate i with T terms at the
   working precision tw.  For moderate parameters, measured per term
   (0.4 to 1 microseconds at 64 and 333 bits, 0.9 to 2.7 at 3333 bits;
   the series and the connection formulas alike, the limit cases 3.2
   times as much), with the gamma functions and powers of the
   connection formulas apart.  Parameters with many bits (1/3 rounded
   at a high precision) cost a multiplication of that size for each
   term: 2.3 instead of 0.45 microseconds at 3333 bits. */
static double
f_model_term_cost(const f_model * M, int i, double tw)
{
    double u;

    if (M->moderate)
    {
        u = M->real ? 0.36 + 1.5e-4 * tw : M->cpar ? 0.78 + 5.6e-4 * tw : 0.5 + 2.8e-4 * tw;
        u += 0.8 * (_acb_hypgeom_est_series_term_cost(tw, M->parbits, M->real)
                    - _acb_hypgeom_est_series_term_cost(tw, 0.0, M->real));
    }
    else if (i < F_CONN && M->parbits <= 0.125 * tw)
    {
        /* (A series at z or z / (z - 1) with short parameters and a
           long z: z enters by its powers in the rectangular splitting,
           and a term costs a multiplication by a small number.
           Measured with a = 1000, b = 10, c = 3000 at 333, 1000 and
           3333 bits: 0.49, 0.53 and 0.57 microseconds per term for a
           short z, 0.55, 0.72 and 1.1 for a real z of full length, and
           0.8, 0.9 and 1.7 for a complex one.) */
        /* (The argument of the series at z / (z - 1) has all the bits
           of the working precision unless that is a short dyadic
           number: z = -1/4 has two bits and z / (z - 1) = 1/5 all of
           them.  2F1(-1000 - 1/2, 1001 + 1/2; 1; -1/4) at 33333 bits:
           0.44 microseconds per term at z, 3.8 at 1/5.) */
        double zb = (i < 2) ? M->pbits : M->wpbits;

        u = _acb_hypgeom_est_series_term_cost(tw, M->parbits, 1);
        if (zb > 64.0)
            u += 1e-4 * FLINT_MIN(zb, tw);
        u *= (M->real ? 1.6 : M->cpar ? 3.5 : 2.4);
    }
    else
    {
        /* (the argument of a connection formula enters by its powers
           in the rectangular splitting: few multiplications) */
        u = 1.6 * _acb_hypgeom_est_series_term_cost(tw,
            (i < 2) ? M->pbits : (i < F_CONN) ? M->wpbits : M->parbits, M->real);
    }

    /* (two series for the function and its derivative, and the two
       continuation steps, whose Taylor series get longer with the
       parameters: 2.2 to 2.8 microseconds per term of one series for
       parameters below 5, 3 to 6 below 20, 5 to 20 below 60) */
    if (i == F_CORNER)
        u *= 4.5 * (1.0 + M->maxpar / 12.0);
    else if (i >= F_CONN && M->E[i - F_CONN].limit)
        u *= 3.2;

    return u;
}

static double
f_model_cost(const f_model * M, int i, double T, double tw)
{
    double u = f_model_term_cost(M, i, tw), beta;

    /* (the tail bound and the set-up of a series: 1 to 5 microseconds
       at 64 bits, 8 to 14 at 333 bits, 50 (real) or 500 at 3333 bits) */
    if (i < F_CONN)
        return T * u + 3.0 + (M->real ? 0.0003 : 0.0025) * pow(tw, 1.5);

    if (i == F_CORNER)
        return T * u + 30.0;

    /* (the gamma functions and powers of a connection formula: 13 to
       16 microseconds at 64 bits, 45 to 50 at 333 bits, 700 to 1100 at
       3333 bits; with complex parameters 43, 170 and 4100) */
    beta = 0.12 * tw * (1.0 + tw / 3000.0) * (M->cpar ? 3.3 : 1.0) + 5.0;

    if (M->E[i - F_CONN].limit)
        beta = 1.5 * beta + 5.0;

    return T * u + beta;
}

/* The working precision of the candidate i for base bits with the
   cancellation cc: with the bits which it loses otherwise, and with as
   many bits as the terms of a series require (their number is limited
   to 10 times the precision); 0 if there is none below 10^8, or if
   the candidate costs more than maxcost (the terms are counted only as
   far as that).  The number of terms is stored in T. */
static double
f_model_precision(double * T, const f_model * M, int i, double base, double cc, double maxcost)
{
    double tw, twmin, tper, need, kmax;
    int j;

    /* (guard bits: with moderate parameters, the candidates lose 4 to
       13 bits more than the estimated cancellation, the limit cases up
       to 19) */
    twmin = base + cc + floor((base + cc) / 64.0);
    if (!M->moderate)
        twmin += 20.0;
    else if (i >= F_CONN && i != F_CORNER && M->E[i - F_CONN].limit)
        twmin += 18.0;
    else
        twmin += 12.0;

    tw = twmin;
    *T = 0.0;
    kmax = (maxcost < 1e300) ? maxcost / f_model_term_cost(M, i, twmin) + 16.0 : 1e300;

    for (j = 0; j < 8; j++)
    {
        *T = f_model_terms(&tper, M, i, tw, kmax);

        /* (more terms than a series is allowed to have at any
           precision) */
        if (tper >= 1e15)
            return 0.0;

        need = twmin + f_model_extra(M, i, tper);
        need = FLINT_MAX(need, tper / 9.0);

        if (need <= tw + 1.0)
            break;

        tw = FLINT_MAX(1.1 * tw, need);

        /* (the terms do not get fewer than 9 per bit) */
        if (tw > 16.0 * twmin + 1000.0)
            return 0.0;
    }

    if (j == 8 || tw > 1e8 || *T > 1e9)
        return 0.0;

    return tw;
}

/* the candidate i at the precision wp */
static void
f_model_evaluate(acb_t res, int i, const acb_t a, const acb_t b, const acb_t c,
    const acb_t z, int flags, slong wp)
{
    if (i < F_CONN)
        f_series_form(res, i, a, b, c, z, flags, wp);
    else if (i == F_CORNER)
        acb_hypgeom_2f1_corner(res, a, b, c, z, flags & ACB_HYPGEOM_2F1_REGULARIZED, wp);
    else if (i == F_CONN + 4)
        acb_hypgeom_2f1_transform(res, b, a, c, z, f_swap_flags(flags), 5, wp);
    else
        acb_hypgeom_2f1_transform(res, a, b, c, z, flags, i - F_CONN + 2, wp);
}

/* The costs of the candidates for results with base bits: the working
   precisions wps, the numbers of terms T and the cancellations cancel
   (in bits, from the largest terms Lam and the size lf of the function).
   The candidates are taken in the order of a first guess of their costs
   (the number of terms from the argument alone), so that the terms of
   the expensive ones are counted only as far as the cost of the best
   one, and no further than limit.  Candidates which cancel by more than
   ccmax bits are left out.  Returns the least cost (1e300 if there is no
   candidate). */
static double
f_model_costs(double * cost, double * T, slong * wps, slong * cancel,
    const f_model * M, const int * valid, const int * runs, const double * Lam,
    double lf, double base, double ccmax, double limit)
{
    double proxy[F_NC], csofar = 1e300;
    int i, j, order;

    for (i = 0; i < F_NC; i++)
    {
        double cc, lw;

        cost[i] = 1e300;
        wps[i] = 0;
        cancel[i] = 0;
        T[i] = 0.0;
        proxy[i] = 1e300;

        if (!valid[i] || runs[i] >= 4)
            continue;

        cc = (Lam[i] - lf) * F_INVLOG2;
        cc = FLINT_MIN(FLINT_MAX(cc, 0.0), 1e9);
        cancel[i] = (slong) cc;

        if (cc > ccmax)
            continue;

        lw = (i < F_CONN) ? M->flz[i] : (i == F_CORNER) ? F_CORNER_LZ : M->E[i - F_CONN].lw;
        proxy[i] = (base + cc + 20.0) / (-lw) * ((i < F_CONN) ? 1.0 : 2.0)
            * f_model_term_cost(M, i, base + cc);
    }

    for (order = 0; order < F_NC; order++)
    {
        double tw;

        i = -1;
        for (j = 0; j < F_NC; j++)
        {
            if (proxy[j] < 1e300 && (i < 0 || proxy[j] < proxy[i]))
                i = j;
        }

        if (i < 0)
            break;

        proxy[i] = 1e300;
        tw = f_model_precision(T + i, M, i, base, (double) cancel[i], FLINT_MIN(csofar, limit));

        if (tw == 0.0)
            continue;

        wps[i] = (slong) tw;
        cost[i] = f_model_cost(M, i, T[i], tw);
        csofar = FLINT_MIN(csofar, cost[i]);

        if (cost[i] > limit)
            cost[i] = 1e300;
    }

    return csofar;
}

void
acb_hypgeom_2f1(acb_t res, const acb_t a, const acb_t b,
        const acb_t c, const acb_t z, int flags, slong prec)
{
    acb_t best, t;
    f_cd A, B, C, Z, W, CA, CB, one_z, l1z, u;
    f_model M;
    double Lam[F_NC], pref[4], lf, lf2, budget, spent, cap;
    slong goal, accept, acc, bestacc;
    int valid[F_NC], runs[F_NC];
    int real, exact, int_ok, tried_int, tried_retry, special, have_model, connection, iter, i;
    int loop_first = 0, lazy_input;
    double cref = 1e300, spare;
    double cost[F_NC + 2], T[F_NC];
    slong wps[F_NC], cancel[F_NC];
    slong std_acc = -WORD_MAX, in_min;
    int have_costs = 0, small, easy, tried_std = 0, std_cand, zchoice;

    if (!acb_is_finite(a) || !acb_is_finite(b) || !acb_is_finite(c) || !acb_is_finite(z))
    {
        acb_indeterminate(res);
        return;
    }

    acb_init(best);
    acb_init(t);

    real = acb_is_real(a) && acb_is_real(b) && acb_is_real(c) && acb_is_real(z);
    exact = acb_is_exact(a) && acb_is_exact(b) && acb_is_exact(c) && acb_is_exact(z);
    lazy_input = !exact && !(_acb_hypgeom_inexact_is_coarse(a, prec)
        || _acb_hypgeom_inexact_is_coarse(b, prec) || _acb_hypgeom_inexact_is_coarse(c, prec)
        || _acb_hypgeom_inexact_is_coarse(z, prec));

    /* the accuracy that the input is good for (the relative condition
       number for an inexact parameter is about its absolute value) */
    in_min = acb_rel_accuracy_bits(z);
    acc = _acb_hypgeom_est_input_accuracy(a);
    in_min = FLINT_MIN(in_min, acc);
    acc = _acb_hypgeom_est_input_accuracy(b);
    in_min = FLINT_MIN(in_min, acc);
    acc = _acb_hypgeom_est_input_accuracy(c);
    in_min = FLINT_MIN(in_min, acc);

    /* Close to z = 1 the function depends on 1 - z, which has fewer
       accurate bits than z: with s = c - a - b it behaves like
       (1 - z)^s for Re(s) < 0, and its derivative like (1 - z)^(s - 1)
       for Re(s) < 1, with a factor s (or 1 / log(1 - z) for small |s|).
       The bits by which 1 - z is less accurate than z are lost for
       Re(s) <= 0, and their part 1 - Re(s) for Re(s) up to 1. */
    if (!acb_is_exact(z) && arf_cmpabs_2exp_si(arb_midref(acb_imagref(z)), -1) < 0
        && arf_cmp_d(arb_midref(acb_realref(z)), 0.5) > 0
        && arf_cmp_d(arb_midref(acb_realref(z)), 1.5) < 0)
    {
        slong zacc = acb_rel_accuracy_bits(z), lost;

        acb_sub_ui(t, z, 1, FLINT_MAX(prec, zacc) + 64);
        lost = zacc - acb_rel_accuracy_bits(t);

        if (lost > 0 && lost < WORD_MAX / 4)
        {
            double sr, si, w, f;

            sr = arf_get_d(arb_midref(acb_realref(c)), ARF_RND_NEAR)
                - arf_get_d(arb_midref(acb_realref(a)), ARF_RND_NEAR)
                - arf_get_d(arb_midref(acb_realref(b)), ARF_RND_NEAR);
            si = arf_get_d(arb_midref(acb_imagref(c)), ARF_RND_NEAR)
                - arf_get_d(arb_midref(acb_imagref(a)), ARF_RND_NEAR)
                - arf_get_d(arb_midref(acb_imagref(b)), ARF_RND_NEAR);

            w = FLINT_MIN(FLINT_MAX(1.0 - sr, 0.0), 1.0);
            f = fmin(fmax(hypot(sr, si), 1.0 / (1.0 + 0.69 * (double) lost)), 1.0);
            w = w * (double) lost + log(f) * F_INVLOG2;

            if (w > 0.0)
                in_min = FLINT_MIN(in_min, zacc - (slong) w);
        }
    }

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

    A.re = arf_get_d(arb_midref(acb_realref(a)), ARF_RND_NEAR);
    A.im = arf_get_d(arb_midref(acb_imagref(a)), ARF_RND_NEAR);
    B.re = arf_get_d(arb_midref(acb_realref(b)), ARF_RND_NEAR);
    B.im = arf_get_d(arb_midref(acb_imagref(b)), ARF_RND_NEAR);
    C.re = arf_get_d(arb_midref(acb_realref(c)), ARF_RND_NEAR);
    C.im = arf_get_d(arb_midref(acb_imagref(c)), ARF_RND_NEAR);
    Z.re = arf_get_d(arb_midref(acb_realref(z)), ARF_RND_NEAR);
    Z.im = arf_get_d(arb_midref(acb_imagref(z)), ARF_RND_NEAR);

    /* cases with formulas of their own: polynomials, z = 0 or 1, poles */
    special = acb_is_zero(z) || acb_is_one(z) || f_is_nonpositive_int(a) || f_is_nonpositive_int(b)
        || f_is_nonpositive_int(c) || acb_eq(a, c) || acb_eq(b, c);
    if (!special && acb_is_exact(c))
    {
        acb_sub(t, c, a, prec);
        special = f_is_nonpositive_int(t);
        acb_sub(t, c, b, prec);
        special = special || f_is_nonpositive_int(t);
    }

    have_model = !special && hypot(A.re, A.im) < 1e12 && hypot(B.re, B.im) < 1e12
        && hypot(C.re, C.im) < 1e12 && hypot(Z.re, Z.im) < 1e12;

    /* (a ball for z which meets the cut [1, inf) from both sides: no
       formula gives more than the jump of the function) */
    if (have_model && !arb_is_exact(acb_imagref(z)) && arb_contains_zero(acb_imagref(z)))
    {
        acb_sub_ui(t, z, 1, prec);
        if (!arb_is_negative(acb_realref(t)))
            have_model = 0;
    }

    /* the four series */
    one_z.re = 1.0 - Z.re; one_z.im = -Z.im;
    CA.re = C.re - A.re; CA.im = C.im - A.im;
    CB.re = C.re - B.re; CB.im = C.im - B.im;

    if (hypot(one_z.re, one_z.im) > 1e-12)
    {
        u.re = -Z.re; u.im = -Z.im;
        W = _acb_hypgeom_cd_div(u, one_z);
        l1z = _acb_hypgeom_cd_log(one_z);
    }
    else
    {
        W.re = 1e300; W.im = 0.0;
        l1z.re = l1z.im = 0.0;
    }

    M.fa[0] = A;  M.fb[0] = B;  M.flz[0] = log(hypot(Z.re, Z.im) + 1e-300);
    M.fa[1] = CA; M.fb[1] = CB; M.flz[1] = M.flz[0];
    M.fa[2] = A;  M.fb[2] = CB; M.flz[2] = log(hypot(W.re, W.im) + 1e-300);
    M.fa[3] = CA; M.fb[3] = B;  M.flz[3] = M.flz[2];
    M.C = C;

    /* Parameters of moderate size, or z where a connection formula
       (at 1/z, 1 - z, ...) converges faster than the four series: the
       choice by z alone. */
    zchoice = acb_hypgeom_2f1_choose(z);
    connection = (zchoice >= 2);

    /* (the domain of the integrals: found when it is needed) */
    int_ok = -1;

    small = hypot(A.re, A.im) <= 30.0 && hypot(B.re, B.im) <= 30.0 && hypot(C.re, C.im) <= 30.0;

    /* (Also larger parameters with the series at z, if its terms
       decrease by a factor of four or more each until they reach the
       precision: their sum is then within a third of the first term,
       and nothing is lost.) */
    if (!small && have_model && zchoice == 0)
    {
        double zabs = hypot(Z.re, Z.im);
        double bits = 0.0, k, r;

        for (k = 0.0; bits < prec + 10.0 && k < prec; k += 1.0)
        {
            r = zabs * hypot(A.re + k, A.im) * hypot(B.re + k, B.im)
                / (hypot(C.re + k, C.im) * (k + 1.0));
            if (!(r <= 0.25))
                break;
            bits -= log(r + 1e-300) * F_INVLOG2;
        }

        small = (bits >= prec + 10.0);
    }

    /* (With parameters up to 12, also a connection formula whose
       argument is at most 1/2 and which is not a limit case: no other
       one has many fewer terms, and
       the estimates cost a part of what the series do.  If the result
       is not accurate, the model follows.  In an experiment with
       random parameters and z, this is 10 to 20 percent faster than
       the model in total; with parameters between 12 and 30 it is 20
       percent slower, as the integer ones lose more bits than that.) */
    easy = 0;
    if (connection && have_model && hypot(A.re, A.im) <= F_EASY_MAX_PARAMETER
        && hypot(B.re, B.im) <= F_EASY_MAX_PARAMETER && hypot(C.re, C.im) <= F_EASY_MAX_PARAMETER)
    {
        double w1 = hypot(Z.re, Z.im), w2 = hypot(one_z.re, one_z.im);

        w1 = FLINT_MIN(1.0 / w1, w2 / w1);
        w2 = FLINT_MIN(1.0 / w2, w2);
        easy = (FLINT_MIN(w1, w2) <= 0.5);

        /* (not the limit cases, which cost three times as much per
           term and lose more bits: a series at z or z / (z - 1) with
           twice the terms is faster) */
        if (easy)
        {
            double dr = (w1 < w2) ? B.re - A.re : C.re - A.re - B.re;
            double di = (w1 < w2) ? B.im - A.im : C.im - A.im - B.im;

            if (di == 0.0 && dr == floor(dr))
                easy = 0;
        }
    }

    if (!have_model || (!connection && small) || easy)
    {
        _acb_hypgeom_2f1_nointegration(best, a, b, c, z, flags, prec + (easy ? F_EASY_GUARD : 0));
        tried_std = 1;
        bestacc = std_acc = acb_rel_accuracy_bits(best);

        /* (for inexact input that is as accurate as the precision, a
           few bits less are accepted up to 1000 bits: the function at
           the midpoints with a bound for the variation costs several
           times as much there) */
        if (bestacc >= accept - ((lazy_input && prec <= 1000) ? 18 : 0) || !have_model)
        {
            if (bestacc < accept && goal > 2)
            {
                double nterms = exact ? f_poly_terms(a, b, c, prec) : 0.0;

                int_ok = f_integral_domain(a, b, c, z, A, B, C, Z, prec);

                slong iprec = (goal < prec - 30) ? goal + 30 : prec;

                if (nterms > 0.0)
                {
                    /* A polynomial which cancels, with exact
                       parameters and argument: an integral, if it
                       costs less than a few attempts with more
                       precision; the polynomial with more precision;
                       the integral if that has failed. */
                    double cint = f_integral_cost(int_ok, iprec, real);
                    int tried = 0;

                    if (bestacc < goal / 2
                        && cint < 4.0 * f_poly_cost(a, b, c, z, 2.0 * prec, nterms, real))
                    {
                        f_integrate(best, &bestacc, a, b, c, z, flags, iprec, goal / 2);
                        tried = 1;
                    }

                    if (bestacc < accept)
                        f_poly_refine(best, &bestacc, a, b, c, z, flags, prec, accept, nterms, real, 1e300,
                            f_poly_budget(tried ? 1e300 : cint, f_poly_cost(a, b, c, z, prec, nterms, real)));

                    if (bestacc < goal / 2 && !tried && cint < 1e300)
                        f_integrate(best, &bestacc, a, b, c, z, flags, iprec, goal / 2);
                }
                else
                {
                    nterms = f_poly_terms(a, b, c, prec);

                    if (nterms > 0.0)
                    {
                        /* A polynomial with an inexact argument or
                           inexact other parameters: the integral on
                           the input, if it costs less than a few
                           attempts at the polynomial; otherwise at the
                           midpoints, with a bound for the variation
                           (see inexact.c). */
                        double cexp = int_ok ? 2.5 * _acb_hypgeom_est_integration_cost(64, 0) : 1e300;
                        double cpoly = f_poly_cost(a, b, c, z, 2.0 * prec, nterms, real);
                        int status, taylor;

                        taylor = acb_is_exact(a) && acb_is_exact(b) && acb_is_exact(c)
                            && f_poly_analytic(a, b, c, z, prec);

                        /* (input which is more accurate than the
                           precision: first the polynomial with more
                           precision, as far as the input allows) */
                        if (spare >= 30.0)
                        {
                            f_poly_refine(best, &bestacc, a, b, c, z, flags, prec, accept,
                                nterms, real, prec + spare - 10.0, 5e4);

                            if (bestacc >= accept)
                            {
                                acb_set_round(res, best, prec);
                                acb_clear(best);
                                acb_clear(t);
                                return;
                            }
                        }

                        if (bestacc < goal / 2 && f_integral_cost(int_ok, iprec, real) < 4.0 * cpoly)
                        {
                            f_integrate(best, &bestacc, a, b, c, z, flags, iprec,
                                FLINT_MIN(accept, goal - 20));

                            /* (the Taylor expansion would cost two
                               more integrals or polynomials) */
                            if (bestacc >= accept - 18)
                            {
                                acb_set_round(res, best, prec);
                                acb_clear(best);
                                acb_clear(t);
                                return;
                            }
                        }

                        /* Only the argument is inexact: the Taylor
                           expansion at its midpoint, whose terms are
                           polynomials of the same kind with exact
                           arguments. */
                        if (taylor)
                        {
                            f_taylor_ctx ctx;
                            slong order;

                            ctx.a = a;
                            ctx.b = b;
                            ctx.c = c;
                            ctx.flags = flags;

                            if (_acb_hypgeom_inexact_taylor(t, &order, z, f_taylor_coeff,
                                    f_taylor_recur, &ctx, real, prec, goal,
                                    F_TAYLOR_MAX_ORDER, F_TAYLOR_MAX_VALUES))
                            {
                                acc = acb_rel_accuracy_bits(t);

                                if (acc > bestacc || (acb_is_finite(t) && !acb_is_finite(best)))
                                {
                                    acb_swap(best, t);
                                    bestacc = acc;
                                }

                                /* (the first-order term is the variation
                                   of the function over the input) */
                                acb_set_round(res, best, prec);
                                acb_clear(best);
                                acb_clear(t);
                                return;
                            }
                        }

                        status = f_inexact(t, a, b, c, z, flags, prec, goal, 0.0,
                            cexp / cpoly, lazy_input);

                        if (status != 0)
                        {
                            acc = acb_rel_accuracy_bits(t);

                            if (acc > bestacc || (acb_is_finite(t) && !acb_is_finite(best)))
                            {
                                acb_swap(best, t);
                                bestacc = acc;
                            }
                        }
                    }

                    /* (without estimates, integration if less than
                       half of the bits are accurate) */
                    if (bestacc < goal / 2)
                        f_integrate(best, &bestacc, a, b, c, z, flags, iprec, goal / 2);
                }
            }

            acb_set_round(res, best, prec);
            acb_clear(best);
            acb_clear(t);
            return;
        }
    }

    if (int_ok < 0)
        int_ok = f_integral_domain(a, b, c, z, A, B, C, Z, prec);

    M.parbits = fmax(fmax(acb_bits(a), acb_bits(b)), acb_bits(c));
    M.pbits = fmax(M.parbits, acb_bits(z));

    /* (z / (z - 1), as the series at it gets it) */
    acb_sub_ui(t, z, 1, prec);
    acb_div(t, z, t, prec);
    M.wpbits = fmax(M.parbits, acb_bits(t));
    M.maxpar = fmax(hypot(A.re, A.im), fmax(hypot(B.re, B.im), hypot(C.re, C.im)));
    M.real = real;
    M.cpar = (A.im != 0.0 || B.im != 0.0 || C.im != 0.0);
    M.intpar = !M.cpar && A.re == floor(A.re) && B.re == floor(B.re) && C.re == floor(C.re);
    M.moderate = (M.maxpar <= 256.0);

    pref[0] = 0.0;
    pref[1] = (CA.re - B.re) * l1z.re - (CA.im - B.im) * l1z.im;
    pref[2] = -(A.re * l1z.re - A.im * l1z.im);
    pref[3] = -(B.re * l1z.re - B.im * l1z.im);

    cap = 1e300;
    for (i = 0; i < F_NC; i++)
    {
        valid[i] = 0;
        runs[i] = 0;
        Lam[i] = 1e300;
    }

    for (i = 0; i < 4; i++)
    {
        /* (The series are candidates for arguments below 0.99.  Their
           number of terms is limited to 10 per bit of the working
           precision, which is reached at 0.933 unless the parameters
           make the terms decrease faster: see f_model_precision.) */
        valid[i] = (M.flz[i] < F_MAX_LOG_ARGUMENT);

        if (valid[i])
        {
            M.pk[i] = _acb_hypgeom_est_2f1_peak(M.kpk + i, M.fa[i], M.fb[i], C, M.flz[i]);
            Lam[i] = pref[i] + M.pk[i];
            cap = fmin(cap, Lam[i] + log(M.kpk[i] + 1.0 / (-M.flz[i]) + 6.0));
        }
    }

    /* The connection formulas, where their argument is below 0.99:
       the largest term of their two series with the gamma factors and
       powers. */
    {
        for (i = F_CONN; i < F_CONN + 5; i++)
        {
            int which = (i == F_CONN + 4) ? 5 : i - F_CONN + 2;

            if (i == F_CONN + 4)
                valid[i] = _acb_hypgeom_est_2f1_conn(M.E + (i - F_CONN), which, B, A, C, Z);
            else
                valid[i] = _acb_hypgeom_est_2f1_conn(M.E + (i - F_CONN), which, A, B, C, Z);

            valid[i] = valid[i] && (M.E[i - F_CONN].lw < F_MAX_LOG_ARGUMENT);

            /* (the limit is computed only if the difference b - a or
               c - a - b is an exact integer, or if the caller's flags
               say that it is an integer) */
            if (valid[i] && M.E[i - F_CONN].limit
                && !(flags & ((which <= 3) ? ACB_HYPGEOM_2F1_AB : ACB_HYPGEOM_2F1_ABC)))
            {
                if (which <= 3)
                {
                    acb_sub(t, b, a, prec);
                }
                else
                {
                    acb_sub(t, c, a, prec);
                    acb_sub(t, t, b, prec);
                }

                if (!acb_is_int(t))
                    valid[i] = 0;
            }

            if (valid[i])
            {
                Lam[i] = M.E[i - F_CONN].Lam;
                cap = fmin(cap, Lam[i] + log(2.0 * (FLINT_MAX(M.E[i - F_CONN].kpk[0],
                    M.E[i - F_CONN].kpk[1]) + 1.0 / (-M.E[i - F_CONN].lw) + 6.0)));
            }
        }

        /* (a and b exchanged give the same series if a - c + 1 = 1 - b) */
        if (valid[F_CONN + 3] && valid[F_CONN + 4] && A.re == B.re && A.im == B.im)
            valid[F_CONN + 4] = 0;
    }

    /* the corner algorithm, where it would be chosen by z */
    if (acb_hypgeom_2f1_choose(z) == 6 && Z.im != 0.0)
    {
        valid[F_CORNER] = 1;
        M.corner_pk = _acb_hypgeom_est_2f1_peak(&M.corner_kpk, A, B, C, F_CORNER_LZ);
        Lam[F_CORNER] = M.corner_pk;
    }

    /* (The function is at most the sum of the terms of any of these:
       the largest term times the number of terms of that size, which
       is about the index of the largest term plus the number of terms
       by which the argument alone reduces them by a factor e.) */
    if (cap < 1e300)
        _acb_hypgeom_est_2f1_logabs(&lf, &lf2, A, B, C, Z, cap);
    else
        lf = lf2 = 0.0;

    /* (The algorithm chosen by z, if it was evaluated above and its
       result is not accurate: the bits which it has lost show the size
       of the function, as for a candidate below.) */
    std_cand = -1;
    if (tried_std && connection)
    {
        int which = acb_hypgeom_2f1_choose(z);

        std_cand = (which == 6) ? F_CORNER : F_CONN + which - 2;

        if (!valid[std_cand])
            std_cand = -1;
        else if (std_acc >= 4)
            lf = FLINT_MIN(lf, Lam[std_cand] - (double) (prec + (easy ? F_EASY_GUARD : 0) - std_acc) / F_INVLOG2);
    }

    budget = 5e4 + 5.0 * _acb_hypgeom_est_integration_cost(FLINT_MIN(prec, goal + 30), real);
    spent = 0.0;
    tried_int = tried_retry = 0;

    /* Inexact input: the error bounds of the series on balls are
       inflated by the cancellation, which more precision does not help.
       The function is evaluated at the midpoints, with a bound for its
       variation (see inexact.c for what is attempted at which cost).
       What follows after this is for the case that there is no result
       (z on the cut, a ball containing a zero of the function, or huge
       numbers). */
inexact_input:

    if (!exact)
    {
        double cmain, cint, cexp, cmin = 1e300;
        slong wp = FLINT_MIN(prec, goal + 30);
        int status, lazy, euler, nothing, done = 0;

        /* the cost of the value: a series (or the algorithm chosen by
           z, which is at most as expensive as its series) or the
           integral */
        {
            f_model_costs(cost, T, wps, cancel, &M, valid, runs, Lam, lf,
                (double) FLINT_MIN(prec, goal + 20), 1e300, budget);
            have_costs = 1;

            /* (the cancellation of the cheapest one) */
            cmain = 1e300;
            for (i = 0; i < F_NC; i++)
            {
                if (cost[i] < cmain)
                {
                    cmain = cost[i];
                    cmin = (double) cancel[i];
                }
            }
        }

        if (cmain >= 1e300)
            cmain = 9.0 * wp * 1.6 * _acb_hypgeom_est_series_term_cost(wp, M.pbits, real);

        cint = f_integral_cost(int_ok, wp, real);

        /* (the reference for the cost of an integral below) */
        cref = cmain;

        /* A series that loses little: on the input, as for exact
           input, which saves the bound (up to 8 bits at low precision,
           up to 3 bits up to 1000 bits, where a result within 18 bits
           of what the input allows is accepted: a series with positive
           terms on rounded input costs half as much as the function at
           the midpoints with a bound).  (If that is not accurate, this
           is reached once more.) */
        if (loop_first == 0 && cmin <= ((wp <= 100) ? 8.0 : (wp <= 1000) ? 3.0 : -1.0) && cmain < cint)
        {
            loop_first = 1;
            goto methods;
        }

        /* Input which is more accurate than the precision by more
           bits than a series loses (or than the algorithm chosen by z
           has lost): on the input with more precision, as for exact
           input.  (If that is not accurate, this is reached once
           more.) */
        if (loop_first == 0 && ((cmin < 1e299 && spare >= cmin + 20.0 && cmain < cint)
            || (connection && std_acc >= 4 && spare >= (accept - std_acc) + 30.0)))
        {
            loop_first = 1;
            goto methods;
        }

        /* the bound from the integral on complex balls */
        cexp = int_ok ? 2.5 * _acb_hypgeom_est_integration_cost(64, 0) : 1e300;

        /* (less accurate than the precision: more precision does not
           improve the result) */
        lazy = lazy_input;

        /* The series cancel by more bits than the input has, and the
           bound is out of proportion: the series give no result that
           is worth its cost; the caller can increase the precision.
           (Not in the domain of the Euler integral with z < 1, where
           there is a bound that costs little, or where a connection
           formula is used.) */
        euler = (Z.re < 1.0) && ((B.re > 0.0 && CB.re > 0.0) || (A.re > 0.0 && CA.re > 0.0));
        nothing = lazy && !connection && cmin < 1e299 && cmin > goal + 11.0
            && cexp > 2.0 * FLINT_MIN(cmain, cint) && !euler;

        /* Where the value would be computed by the Euler integral:
           the integral on the input, which encloses the function
           without the loss of the series.  Where the series on the
           input give nothing (no result at all, or one from the
           midpoints with bounds that are computed whatever they cost),
           also an integral of up to four times the cost of the series
           at the precision at which they succeed. */
        if (int_ok && cint <= ((nothing || (!lazy && cmin < 1e299 && cmin > goal + 11.0)) ? 4.0 : 1.0) * cmain)
        {
            f_integrate(best, &bestacc, a, b, c, z, flags,
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

        status = f_inexact(t, a, b, c, z, flags, prec, goal, (cmin < 1e299) ? cmin : 0.0,
                    cexp / cmain, lazy);

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
                f_integrate(best, &bestacc, a, b, c, z, flags,
                    (goal < prec - 30) ? goal + 30 : prec, FLINT_MIN(accept, goal - 20));

            acb_set_round(res, best, prec);
            acb_clear(best);
            acb_clear(t);
            return;
        }
    }

methods:

    for (iter = 0; iter < 10 && bestacc < accept && loop_first != 2; iter++)
    {
        double cbest, csofar, ccmax;
        int method;

        /* No more precision than an inexact input is good for (and as
           much as the terms of the series require).  With input that is
           not more accurate than a series loses, only the series that
           lose little are candidates. */
        ccmax = 1e300;
        if (!exact)
        {
            ccmax = (int_ok && !tried_int) ? FLINT_MAX(12.0, goal / 8.0) : goal / 2.0 + 10.0;
            ccmax = FLINT_MAX(ccmax, spare - 20.0);
        }

        if (have_costs)
        {
            /* (those of the block for inexact input) */
            csofar = 1e300;
            for (i = 0; i < F_NC; i++)
            {
                if ((double) cancel[i] > ccmax)
                    cost[i] = 1e300;
                csofar = FLINT_MIN(csofar, cost[i]);
            }
            have_costs = 0;
        }
        else
        {
            csofar = f_model_costs(cost, T, wps, cancel, &M, valid, runs, Lam, lf,
                (double) FLINT_MIN(prec, goal + 20), ccmax, budget - spent);
        }

        /* The algorithm chosen by z alone, with the missing bits (for
           the connection formulas) */
        cost[F_RETRY] = 1e300;
        if (!tried_retry && std_acc >= 4 && std_acc < accept && connection && std_cand < 0
            && (exact || spare >= (accept - std_acc) + 30.0))
            cost[F_RETRY] = 0.0;

        cost[F_INTEGRAL] = tried_int ? 1e300 : f_integral_cost(int_ok, FLINT_MIN(prec, goal + 30), real);

        /* (for inexact input that is as accurate as the precision:
           no integral at many times the cost of the series, whose
           results are limited by the input) */
        if (iter == 0 && !(cref < 1e300))
        {
            for (i = 0; i < F_NC; i++)
                cref = FLINT_MIN(cref, cost[i]);
        }
        if (lazy_input && cref < 1e300 && cost[F_INTEGRAL] > 4.0 * cref)
            cost[F_INTEGRAL] = 1e300;

        method = 0;
        cbest = cost[0];
        for (i = 1; i < F_NC + 2; i++)
        {
            if (cost[i] < cbest)
            {
                cbest = cost[i];
                method = i;
            }
        }

        /* No candidate with an estimate where a connection formula
           would be chosen from z: the algorithm chosen by z, at the
           precision. */
        if (iter == 0 && connection && !tried_std && csofar >= 1e300)
        {
            _acb_hypgeom_2f1_nointegration(t, a, b, c, z, flags, prec);
            tried_std = 1;
            std_acc = acc = acb_rel_accuracy_bits(t);
            if (acc > bestacc || (acb_is_finite(t) && !acb_is_finite(best)))
            {
                acb_swap(best, t);
                bestacc = acc;
            }
            continue;
        }

        if (cbest >= 1e300)
        {
            /* nothing else: the algorithm chosen by z alone */
            if (!tried_std)
            {
                _acb_hypgeom_2f1_nointegration(t, a, b, c, z, flags, prec);
                tried_std = 1;
                std_acc = acc = acb_rel_accuracy_bits(t);
                if (acc > bestacc || (acb_is_finite(t) && !acb_is_finite(best)))
                {
                    acb_swap(best, t);
                    bestacc = acc;
                }
                continue;
            }
            break;
        }

        spent += cbest;

        if (method == F_INTEGRAL)
        {
            f_integrate(best, &bestacc, a, b, c, z, flags,
                (goal < prec - 30) ? goal + 30 : prec, FLINT_MIN(accept, goal - 20));
            tried_int = 1;

            /* (the integral encloses the function on inexact input
               without the loss of the series) */
            if (!exact && acb_is_finite(best))
                break;

            continue;
        }
        else if (method == F_RETRY)
        {
            _acb_hypgeom_2f1_nointegration(t, a, b, c, z, flags, prec + (accept - std_acc) + 20);
            acc = acb_rel_accuracy_bits(t);
            tried_retry = 1;
        }
        else
        {
            slong macc;

            i = method;

            /* (the input limits the result of a series that cancels:
               no more precision than the input has) */
            if (!exact && wps[i] > FLINT_MAX(in_min, 32) + 30)
                wps[i] = FLINT_MAX(in_min, 32) + 30;

            f_model_evaluate(t, i, a, b, c, z, flags, wps[i]);
            acc = macc = acb_rel_accuracy_bits(t);
            runs[i]++;

            if (!exact && spare < cancel[i] + 20.0)
            {
                int j;

                /* (more precision does not improve the error bounds
                   for inexact parameters, nor does a candidate which
                   cancels as much) */
                for (j = 0; j < F_NC; j++)
                {
                    if (j == i || cancel[j] + 8 >= cancel[i])
                        runs[j] = 4;
                }
            }
            else if (macc < accept)
            {
                /* The cancellation was underestimated.  If the result
                   is not all noise, this shows by how much. */
                double loss;

                if (macc >= 4)
                {
                    loss = (double) cancel[i] + (accept - macc) + 10.0;
                    lf = FLINT_MIN(lf, Lam[i] - loss / F_INVLOG2);
                }
                else if (lf2 < lf - 20.0)
                {
                    /* the other saddle point */
                    lf = lf2;
                    lf2 = 1e300;
                    budget = spent + (budget - spent) / 3.0;
                }
                else
                {
                    /* twice the precision */
                    loss = cancel[i] + wps[i];
                    lf = FLINT_MIN(lf, Lam[i] - loss / F_INVLOG2);
                    budget = spent + (budget - spent) / 3.0;
                }
            }
            else if (acc < accept)
            {
                runs[i] = 4;
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

    if (bestacc == -WORD_MAX && !tried_std)
    {
        _acb_hypgeom_2f1_nointegration(best, a, b, c, z, flags, prec);
        bestacc = acb_rel_accuracy_bits(best);
    }

    acb_set_round(res, best, prec);

    acb_clear(best);
    acb_clear(t);
}
