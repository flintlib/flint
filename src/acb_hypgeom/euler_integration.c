/*
    Copyright (C) 2026 Fredrik Johansson
    Developed using Claude Fable 5.2

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/* The Euler integrals of the confluent and Gauss hypergeometric
   functions for complex parameters and argument, along paths of
   steepest descent (exp_log_integration.c). */

#include <math.h>
#include "acb.h"
#include "acb_hypgeom.h"
#include "acb_hypgeom/impl.h"

int
_acb_hypgeom_1f1_integral(acb_t res, const acb_t a1, const acb_t ba1, const acb_t z, slong prec)
{
    acb_struct mu[2];
    acb_t al[2], be[2], dir;
    int ok;

    acb_init(al[0]); acb_init(al[1]); acb_init(be[0]); acb_init(be[1]); acb_init(dir);

    /* t^a1 (1 - t)^ba1 */
    acb_one(be[0]);
    acb_one(al[1]);
    acb_set_si(be[1], -1);
    mu[0] = *a1;
    mu[1] = *ba1;

    {
        acb_struct alv[2], bev[2];
        alv[0] = *al[0]; alv[1] = *al[1];
        bev[0] = *be[0]; bev[1] = *be[1];
        ok = _acb_hypgeom_integrate_exp_log(res, dir, z, mu, alv, bev, 2, 0, 1, 0, -INFINITY, prec);
    }

    acb_clear(al[0]); acb_clear(al[1]); acb_clear(be[0]); acb_clear(be[1]); acb_clear(dir);
    return ok;
}

/* the bits of the largest parameter beyond 2^16: extra precision for
   the gamma functions */
static slong
euler_extra_bits(const acb_t a, const acb_t b)
{
    mag_t m, t;
    slong e;

    mag_init(m);
    mag_init(t);
    acb_get_mag(m, a);
    acb_get_mag(t, b);
    mag_max(m, m, t);
    e = (mag_is_zero(m) || mag_is_inf(m)) ? 0 : MAG_EXP(m);
    mag_clear(m);
    mag_clear(t);
    return FLINT_MIN(FLINT_MAX(e - 16, 0), 2000);
}

/* the bits of the logarithm of a product of gamma functions and powers
   with the parameters a, b, c (of the order of x log(x) for the largest
   one): extra precision for its exponential, and for an argument which
   is computed and is raised to such a power */
static slong
euler_log_bits(const acb_t a, const acb_t b, const acb_t c)
{
    mag_t m, t;
    slong e;

    mag_init(m);
    mag_init(t);
    acb_get_mag(m, a);
    acb_get_mag(t, b);
    mag_max(m, m, t);
    acb_get_mag(t, c);
    mag_max(m, m, t);
    e = (mag_is_zero(m) || mag_is_inf(m)) ? 0 : MAG_EXP(m);
    mag_clear(m);
    mag_clear(t);
    return FLINT_MIN(FLINT_MAX(e, 0), 2000) + 6;
}

/* For -pi < theta < pi and |theta + arg(z)| < pi / 2,

       Gamma(a) U(a, b, z) = int_0^(inf e^(i theta)) e^(-z t) t^(a-1) (1 + t)^(b-a-1) dt

   (the rotation of the contour from the positive real axis does not
   cross the singular point -1), where arg(z) is not restricted to the
   principal range: this computes Gamma(a) U(a, b, z e^(rot pi i)) for
   rot = 0, 1, -1.  The direction is chosen among Gaussian integers,
   closest to the direction of the steepest decay.  On the cut of U,
   z is taken to have the argument pi. */
static int
u_integral(acb_t res, const acb_t a1, const acb_t ba1, const acb_t z, int rot, slong prec)
{
    static const int dirs[15][2] = { { 1, 0 }, { 1, 1 }, { 0, 1 }, { -1, 1 }, { 1, -1 }, { 0, -1 }, { -1, -1 },
        { 2, 1 }, { 1, 2 }, { -1, 2 }, { -2, 1 }, { 2, -1 }, { 1, -2 }, { -1, -2 }, { -2, -1 } };
    const double pi = 3.14159265358979323846;
    acb_struct mu[2], alv[2], bev[2];
    acb_t al[2], be[2], dir, mz;
    arb_t t, u;
    double zr, zi, az, best = 0.3;
    int ok, i, bi = -1;

    if (!acb_is_finite(z))
        return 0;

    zr = arf_get_d(arb_midref(acb_realref(z)), ARF_RND_NEAR);
    zi = arf_get_d(arb_midref(acb_imagref(z)), ARF_RND_NEAR);
    az = atan2(zi, zr) + rot * pi;
    for (i = 0; i < (rot == 0 ? 7 : 15); i++)
    {
        double s = az + atan2(dirs[i][1], dirs[i][0]);
        if (fabs(s) < 0.5 * pi && cos(s) > best + 1e-9)
        {
            best = cos(s);
            bi = i;
        }
    }
    if (bi < 0)
        return 0;

    acb_init(al[0]); acb_init(al[1]); acb_init(be[0]); acb_init(be[1]); acb_init(dir); acb_init(mz);
    arb_init(t); arb_init(u);

    acb_set_si_si(dir, dirs[bi][0], dirs[bi][1]);

    /* |theta + arg(z) + rot pi| < pi / 2 */
    acb_arg(t, z, 64);
    acb_arg(u, dir, 64);
    arb_add(t, t, u, 64);
    arb_const_pi(u, 64);
    arb_mul_si(u, u, rot, 64);
    arb_add(t, t, u, 64);
    arb_abs(t, t);
    arb_const_pi(u, 64);
    arb_mul_2exp_si(u, u, -1);
    ok = arb_lt(t, u);

    if (ok)
    {
        acb_one(be[0]);
        acb_one(al[1]);
        acb_one(be[1]);
        mu[0] = *a1;
        mu[1] = *ba1;
        alv[0] = *al[0]; alv[1] = *al[1];
        bev[0] = *be[0]; bev[1] = *be[1];
        /* exp(-z e^(rot pi i) t) */
        if (rot == 0)
            acb_neg(mz, z);
        else
            acb_set(mz, z);
        ok = _acb_hypgeom_integrate_exp_log(res, dir, mz, mu, alv, bev, 2, 0, -1, 1, -INFINITY, prec);
    }

    acb_clear(al[0]); acb_clear(al[1]); acb_clear(be[0]); acb_clear(be[1]); acb_clear(dir); acb_clear(mz);
    arb_clear(t); arb_clear(u);
    return ok;
}

int
_acb_hypgeom_u_integral(acb_t res, const acb_t a1, const acb_t ba1, const acb_t z, slong prec)
{
    return u_integral(res, a1, ba1, z, 0, prec);
}

/* The functions with Re(b - a) >= 1 and z in the right half-plane, for
   any a, from two integrals: with

       I(+-) = Gamma(b - a) U(b - a, b, z e^(+-pi i))
             = int e^(z t) t^(b-a-1) (1 + t)^(a-1) dt

   over rays from 0 to infinity passing below and above -1,

       M(a, b, z) / Gamma(b) = e^z Gamma(1 - a) / (2 pi i Gamma(b - a))
                                 (e^(pi i b) I(+) - e^(-pi i b) I(-)),

       U(a, b, z) = e^z Gamma(1 - a) / (2 pi i)
                        (e^(pi i (b - a)) I(+) - e^(-pi i (b - a)) I(-))

   (the connection formulas [DLMF 13.2.41] with both signs, without the
   term U(a, b, z) and without the left-hand side respectively).  The
   integral of the first kind does not exist for Re(a) <= 0.  The two
   terms are of the size of the function where it oscillates (they are
   conjugates for real a, b and z > 0); they cancel where the function
   is recessive, for U with z > 0 beyond the turning point. */
static int
_acb_hypgeom_2u_integration(acb_t res, const acb_t a, const acb_t b, const acb_t z,
    int which, int regularized, slong prec)
{
    acb_t t, u, a1, ba1, Ip, Im;
    slong wp = prec + 10 + euler_extra_bits(a, b);
    int ok, real;

    acb_init(t);
    acb_init(u);
    acb_init(a1);
    acb_init(ba1);
    acb_init(Ip);
    acb_init(Im);

    /* the exponents of t and of 1 + t */
    acb_sub(a1, b, a, wp);
    acb_sub_ui(a1, a1, 1, wp);
    acb_sub_ui(ba1, a, 1, wp);

    real = acb_is_real(a) && acb_is_real(b) && acb_is_real(z);

    ok = acb_is_finite(z) && acb_is_finite(a1) && acb_is_finite(ba1)
        && arb_is_nonnegative(acb_realref(a1)) && arb_is_positive(acb_realref(z));
    ok = ok && u_integral(Ip, a1, ba1, z, 1, prec);

    if (ok)
    {
        if (real)
            acb_conj(Im, Ip);
        else
            ok = u_integral(Im, a1, ba1, z, -1, prec);
    }

    if (ok)
    {
        /* e^(pi i s) I(+) - e^(-pi i s) I(-), s = b or b - a */
        if (which == 0)
            acb_set(u, b);
        else
            acb_sub(u, b, a, wp);

        acb_exp_pi_i(t, u, wp);
        acb_mul(Ip, Ip, t, wp);
        acb_neg(u, u);
        acb_exp_pi_i(t, u, wp);
        acb_mul(Im, Im, t, wp);
        acb_sub(Ip, Ip, Im, wp);

        /* e^z Gamma(1 - a) / (2 pi i) */
        acb_exp(t, z, wp);
        acb_mul(Ip, Ip, t, wp);
        acb_sub_ui(t, a, 1, wp);
        acb_neg(t, t);
        acb_gamma(t, t, wp);
        acb_mul(Ip, Ip, t, wp);
        acb_const_pi(t, wp);
        acb_mul_onei(t, t);
        acb_mul_2exp_si(t, t, 1);
        acb_div(Ip, Ip, t, wp);

        if (which == 0)
        {
            acb_sub(t, b, a, wp);
            acb_rgamma(t, t, wp);
            acb_mul(Ip, Ip, t, wp);

            if (!regularized)
            {
                acb_gamma(t, b, wp);
                acb_mul(Ip, Ip, t, wp);
            }
        }

        /* (real on the real line) */
        if (real)
            arb_zero(acb_imagref(Ip));

        ok = acb_is_finite(Ip);
        if (ok)
            acb_set(res, Ip);
    }

    acb_clear(t);
    acb_clear(u);
    acb_clear(a1);
    acb_clear(ba1);
    acb_clear(Ip);
    acb_clear(Im);
    return ok;
}

int
_acb_hypgeom_2f1_integral(acb_t res, const acb_t b1, const acb_t cb1, const acb_t nega,
    const acb_t z, slong prec)
{
    acb_struct mu[3], alv[3], bev[3];
    acb_t al[3], be[3], dir, zero;
    int ok, j;

    acb_init(zero);
    for (j = 0; j < 3; j++)
    {
        acb_init(al[j]);
        acb_init(be[j]);
    }
    acb_init(dir);

    /* t^b1 (1 - t)^cb1 (1 - z t)^nega */
    acb_one(be[0]);
    acb_one(al[1]);
    acb_set_si(be[1], -1);
    acb_one(al[2]);
    acb_neg(be[2], z);
    mu[0] = *b1;
    mu[1] = *cb1;
    mu[2] = *nega;
    for (j = 0; j < 3; j++)
    {
        alv[j] = *al[j];
        bev[j] = *be[j];
    }

    ok = _acb_hypgeom_integrate_exp_log(res, dir, zero, mu, alv, bev, 3, 0, 1, 0, -INFINITY, prec);
    acb_clear(zero);

    for (j = 0; j < 3; j++)
    {
        acb_clear(al[j]);
        acb_clear(be[j]);
    }
    acb_clear(dir);
    return ok;
}

void
acb_hypgeom_1f1_integration(acb_t res, const acb_t a, const acb_t b, const acb_t z,
    int regularized, slong prec)
{
    acb_t t, a1, ba1, I;
    slong wp = prec + 10 + euler_extra_bits(a, b);
    int ok;

    acb_init(t);
    acb_init(a1);
    acb_init(ba1);
    acb_init(I);

    acb_sub_ui(a1, a, 1, wp);
    acb_sub(ba1, b, a, wp);
    acb_sub_ui(ba1, ba1, 1, wp);

    ok = acb_is_finite(z) && acb_is_finite(a1) && acb_is_finite(ba1);

    if (ok && !(arb_is_nonnegative(acb_realref(a1)) && arb_is_nonnegative(acb_realref(ba1))))
    {
        /* Outside of the domain of the Euler integral: from two
           integrals for U, also for exp(z) 1F1(b - a, b, -z). */
        if (arb_is_nonnegative(acb_realref(ba1)))
        {
            ok = _acb_hypgeom_2u_integration(I, a, b, z, 0, regularized, prec);
        }
        else if (arb_is_nonnegative(acb_realref(a1)))
        {
            acb_sub(ba1, b, a, wp);
            acb_neg(a1, z);
            ok = _acb_hypgeom_2u_integration(I, ba1, b, a1, 0, regularized, prec);
            if (ok)
            {
                acb_exp(t, z, wp);
                acb_mul(I, I, t, wp);
            }
        }
        else
            ok = 0;

        if (ok)
            acb_set_round(res, I, prec);
        else
            acb_indeterminate(res);

        acb_clear(t);
        acb_clear(a1);
        acb_clear(ba1);
        acb_clear(I);
        return;
    }

    ok = ok && _acb_hypgeom_1f1_integral(I, a1, ba1, z, prec);

    if (ok)
    {
        /* Gamma(b) / (Gamma(a) Gamma(b - a)) */
        if (!regularized)
        {
            acb_gamma(t, b, wp);
            acb_mul(I, I, t, wp);
        }
        acb_rgamma(t, a, wp);
        acb_mul(I, I, t, wp);
        acb_sub(t, b, a, wp);
        acb_rgamma(t, t, wp);
        acb_mul(I, I, t, wp);
        acb_set_round(res, I, prec);
    }
    else
        acb_indeterminate(res);

    acb_clear(t);
    acb_clear(a1);
    acb_clear(ba1);
    acb_clear(I);
}

static int
_acb_hypgeom_u_integration(acb_t res, const acb_t a, const acb_t b, const acb_t z, slong prec)
{
    acb_t t, a1, ba1, I;
    slong wp = prec + 10 + euler_extra_bits(a, b);
    int ok;

    acb_init(t);
    acb_init(a1);
    acb_init(ba1);
    acb_init(I);

    acb_sub_ui(a1, a, 1, wp);
    acb_sub(ba1, b, a, wp);
    acb_sub_ui(ba1, ba1, 1, wp);

    ok = acb_is_finite(z) && acb_is_finite(a1) && acb_is_finite(ba1)
        && arb_is_nonnegative(acb_realref(a1));
    ok = ok && _acb_hypgeom_u_integral(I, a1, ba1, z, prec);

    if (ok)
    {
        acb_rgamma(t, a, wp);
        acb_mul(res, I, t, wp);
    }

    acb_clear(t);
    acb_clear(a1);
    acb_clear(ba1);
    acb_clear(I);
    return ok;
}

void
acb_hypgeom_u_integration(acb_t res, const acb_t a, const acb_t b, const acb_t z, slong prec)
{
    acb_t res2;
    acb_init(res2);

    if (!_acb_hypgeom_u_integration(res2, a, b, z, prec))
    {
        /* U(a, b, z) = z^(1 - b) U(a - b + 1, 2 - b, z) */
        acb_t c, d;
        acb_init(c);
        acb_init(d);
        acb_sub(c, a, b, prec + 20);
        acb_add_ui(c, c, 1, prec + 20);
        acb_sub_ui(d, b, 2, prec + 20);
        acb_neg(d, d);

        if (_acb_hypgeom_u_integration(res2, c, d, z, prec))
        {
            acb_sub_ui(c, b, 1, prec + 20);
            acb_neg(c, c);
            acb_pow(c, z, c, prec + 10 + euler_extra_bits(b, b));
            acb_mul(res2, res2, c, prec + 10);
        }
        else if (_acb_hypgeom_2u_integration(res2, a, b, z, 1, 0, prec))
        {
            /* (Re(a) < 1 and Re(b - a) >= 1: from two integrals) */
        }
        else if (_acb_hypgeom_2u_integration(res2, c, d, z, 1, 0, prec))
        {
            acb_sub_ui(c, b, 1, prec + 20);
            acb_neg(c, c);
            acb_pow(c, z, c, prec + 10 + euler_extra_bits(b, b));
            acb_mul(res2, res2, c, prec + 10);
        }
        else
            acb_indeterminate(res2);

        acb_clear(c);
        acb_clear(d);
    }

    acb_set_round(res, res2, prec);
    acb_clear(res2);
}

static int
_acb_hypgeom_2f1_integration(acb_t res, const acb_t a, const acb_t b, const acb_t c,
    const acb_t z, int regularized, slong prec)
{
    acb_t t, b1, cb1, nega, I;
    slong wp = prec + 10 + euler_extra_bits(b, c);
    int ok;

    acb_init(t);
    acb_init(b1);
    acb_init(cb1);
    acb_init(nega);
    acb_init(I);

    acb_sub_ui(b1, b, 1, wp);
    acb_sub(cb1, c, b, wp);
    acb_sub_ui(cb1, cb1, 1, wp);
    acb_neg(nega, a);

    ok = acb_is_finite(z) && acb_is_finite(a) && acb_is_finite(b1) && acb_is_finite(cb1)
        && arb_is_nonnegative(acb_realref(b1)) && arb_is_nonnegative(acb_realref(cb1));
    ok = ok && _acb_hypgeom_2f1_integral(I, b1, cb1, nega, z, prec);

    if (ok)
    {
        /* Gamma(c) / (Gamma(b) Gamma(c - b)) */
        if (!(regularized & ACB_HYPGEOM_2F1_REGULARIZED))
        {
            acb_gamma(t, c, wp);
            acb_mul(I, I, t, wp);
        }
        acb_rgamma(t, b, wp);
        acb_mul(I, I, t, wp);
        acb_sub(t, c, b, wp);
        acb_rgamma(t, t, wp);
        acb_mul(res, I, t, wp);
    }

    acb_clear(t);
    acb_clear(b1);
    acb_clear(cb1);
    acb_clear(nega);
    acb_clear(I);
    return ok;
}

/* 2F1(a, b; c; z) for Re(z) < 0 where the Euler integral does not exist
   at 1 (Re(c - b) < 1), if Re(b) >= 1 and Re(a - c) >= 0.  With
   x = 1 - z, and

       J = int_0^1 t^(b-1) (1 - t)^(a-c) (1 - x t)^(-a) dt

   taken from 0 to 1 past the point 1 / x on the lower side (J_-) and
   on the upper side (J_+),

       F(a, b; c; z) = Gamma(c) Gamma(b - c + 1) / (2 pi i Gamma(b))
                         (e^(pi i c) J_- - e^(-pi i c) J_+).

   For z < 0, J_-+ are the Euler integrals of F(a, b; a + b - c + 1; x -+ i0),
   the values on the two sides of its branch cut; for other z, one of
   them is the integral over [0, 1] and the other its continuation
   across the cut.  (The integrals of t^(b-1) (1 - t)^(c-b-1) (1 - z t)^(-a)
   from 0 to infinity above and below the point 1 are e^(+-pi i b) J_-+
   by t -> t / (t - 1); a combination of them does not depend on the
   part beyond 1, and continues the Euler integral of F(a, b; c; z) in
   c.)  The two integrals are of the size of the function where it
   oscillates (z < 0 for real parameters) and around it.

   The contours are the polygons through 1/2 -+ i Y, with the term
   1 - x t given as -+i (1 - x t), whose principal logarithm is
   continuous along them; this is undone by the factors
   e^(-+ pi i a / 2).  For real parameters and z, J_+ is the conjugate
   of J_-.

   This computes (e^(pi i theta) J_- - e^(-pi i theta) J_+) / (2 pi i)
   for a given theta (c for the formula above). */
static int
_acb_hypgeom_2f1_two_integrals(acb_t res, const acb_t a, const acb_t b, const acb_t c,
    const acb_t z, const acb_t theta, slong prec)
{
    acb_struct mu[3], al[3], be[3];
    acb_t x, via, zero, t, u, J[2];
    slong wp = prec + 10 + euler_log_bits(a, b, c);
    double sr, si, Y;
    int j, side, real, ok;

    if (!acb_is_finite(a) || !acb_is_finite(b) || !acb_is_finite(c) || !acb_is_finite(z)
        || !acb_is_finite(theta) || !arb_is_negative(acb_realref(z)))
        return 0;

    acb_init(x); acb_init(via); acb_init(zero); acb_init(t); acb_init(u);
    acb_init(J[0]); acb_init(J[1]);
    for (j = 0; j < 3; j++)
    {
        acb_init(mu + j);
        acb_init(al + j);
        acb_init(be + j);
    }

    /* t^(b-1) (1 - t)^(a-c) */
    acb_sub_ui(mu, b, 1, wp);
    acb_sub(mu + 1, a, c, wp);
    acb_neg(mu + 2, a);
    acb_one(be);
    acb_one(al + 1);
    acb_set_si(be + 1, -1);

    acb_sub_ui(x, z, 1, wp);
    acb_neg(x, x);

    ok = arb_is_nonnegative(acb_realref(mu)) && arb_is_nonnegative(acb_realref(mu + 1));

    /* The height of the polygons: the point 1 / x = sr + i si, which is
       inside of the circle through 0 and 1, well inside of the
       triangle (0, 1/2 + i Y sgn(si), 1). */
    acb_inv(u, x, 64);
    sr = arf_get_d(arb_midref(acb_realref(u)), ARF_RND_NEAR);
    si = fabs(arf_get_d(arb_midref(acb_imagref(u)), ARF_RND_NEAR));
    Y = 1.0;
    if (!(sr > 0.0 && sr < 1.0))
        ok = 0;
    else
    {
        while (Y < 1048576.0 && (Y * sr < 2.0 * si || Y * (1.0 - sr) < 2.0 * si))
            Y *= 2.0;
        if (Y >= 1048576.0)
            ok = 0;
    }

    real = acb_is_real(a) && acb_is_real(b) && acb_is_real(c) && acb_is_real(z)
        && acb_is_real(theta);

    for (side = 0; side < 2 && ok; side++)
    {
        if (side == 1 && real)
        {
            acb_conj(J[1], J[0]);
            break;
        }

        /* -+i (1 - x t), by way of 1/2 -+ i Y */
        acb_onei(al + 2);
        if (side == 0)
            acb_neg(al + 2, al + 2);
        acb_mul(be + 2, al + 2, x, wp);
        acb_neg(be + 2, be + 2);
        acb_set_d_d(via, 0.5, (side == 0) ? -Y : Y);

        ok = _acb_hypgeom_integrate_exp_log_via(J[side], via, zero, mu, al, be, 3, 0, 1,
                -INFINITY, 0, prec);

        /* e^(-+ pi i a / 2) e^(+- pi i theta) */
        acb_mul_2exp_si(t, a, -1);
        acb_sub(t, theta, t, wp);
        if (side == 1)
            acb_neg(t, t);
        acb_exp_pi_i(t, t, wp);
        acb_mul(J[side], J[side], t, wp);
    }

    if (ok)
    {
        if (real)
        {
            /* (J_- - conj(J_-)) / (2 i) */
            arb_set(acb_realref(t), acb_imagref(J[0]));
            arb_zero(acb_imagref(t));
        }
        else
        {
            acb_sub(t, J[0], J[1], wp);
            acb_div_onei(t, t);
            acb_mul_2exp_si(t, t, -1);
        }

        acb_const_pi(u, wp);
        acb_div(t, t, u, wp);

        ok = acb_is_finite(t);
        if (ok)
            acb_set(res, t);
    }

    acb_clear(x); acb_clear(via); acb_clear(zero); acb_clear(t); acb_clear(u);
    acb_clear(J[0]); acb_clear(J[1]);
    for (j = 0; j < 3; j++)
    {
        acb_clear(mu + j);
        acb_clear(al + j);
        acb_clear(be + j);
    }

    return ok;
}

static int
_acb_hypgeom_2f1_integration_left(acb_t res, const acb_t a, const acb_t b, const acb_t c,
    const acb_t z, int regularized, slong prec)
{
    acb_t t, u, x;
    slong wp;
    int ok, real;

    if (!acb_is_finite(a) || !acb_is_finite(b) || !acb_is_finite(c) || !acb_is_finite(z))
        return 0;

    wp = prec + 10 + euler_log_bits(a, b, c);

    acb_init(t); acb_init(u); acb_init(x);

    /* Gamma(b - c + 1) */
    acb_sub(t, b, c, wp);
    acb_add_ui(t, t, 1, wp);
    ok = !(acb_contains_int(t) && !arb_is_positive(acb_realref(t)));

    ok = ok && _acb_hypgeom_2f1_two_integrals(t, a, b, c, z, c, prec);

    if (ok)
    {
        real = acb_is_real(a) && acb_is_real(b) && acb_is_real(c) && acb_is_real(z);

        /* Gamma(c) Gamma(b - c + 1) / Gamma(b) */
        acb_sub(u, b, c, wp);
        acb_add_ui(u, u, 1, wp);
        acb_lgamma(u, u, wp);
        acb_lgamma(x, b, wp);
        acb_sub(u, u, x, wp);
        if (!(regularized & ACB_HYPGEOM_2F1_REGULARIZED))
        {
            acb_lgamma(x, c, wp);
            acb_add(u, u, x, wp);
        }
        acb_exp(u, u, wp);
        acb_mul(t, t, u, wp);

        if (real)
            arb_zero(acb_imagref(t));

        ok = acb_is_finite(t);
        if (ok)
            acb_set(res, t);
    }

    acb_clear(t); acb_clear(u); acb_clear(x);
    return ok;
}

/* 2F1(a, b; c; z) for Re(a) <= 0 and Re(b) >= 1 inside of the circle
   |z - 1/2| = 1/2 (Re(1/z) > 1), where the Euler integral does not
   exist at 1 for c < b + 1 (the functions of the type of the Jacobi
   polynomials, a < 0 < c < b).  The integral over the loop from 0
   around 1 in the positive sense [DLMF 15.6.3],

       F(a, b; c; z) = Gamma(c) Gamma(1 + b - c) / (2 pi i Gamma(b))
                         int_0^(1+) t^(b-1) (t - 1)^(c-b-1) (1 - z t)^(-a) dt,

   has an integrand which vanishes at 1 / z for Re(a) < 0: the loop is
   the path from 0 to 1 / z below the point 1 minus the path above it,
   with (t - 1)^(c-b-1) = e^(-+ pi i (c-b-1)) (1 - t)^(c-b-1) on them.
   With t = s / z these are z^(-b) times the integrals J_-+ above for
   the parameters A = b - c + 1, B = b, C = a + b - c + 1 at 1 - 1/z:

       F(a, b; c; z) = Gamma(c) Gamma(A) / (2 pi i Gamma(b)) z^(-b)
                         (e^(pi i A) J_- - e^(-pi i A) J_+),

   which differs from the formula for F(A, B; C; 1 - 1/z) in the phases
   (A in place of C) and in the factors.  The two integrals are of the
   size of the function where it oscillates (0 < z < 1 for real
   parameters, between the turning points). */
static int
_acb_hypgeom_2f1_integration_loop(acb_t res, const acb_t a, const acb_t b, const acb_t c,
    const acb_t z, int regularized, slong prec)
{
    acb_t A, C, Z, t, u;
    slong wp;
    int ok, real;

    if (!acb_is_finite(a) || !acb_is_finite(b) || !acb_is_finite(c) || !acb_is_finite(z)
        || acb_contains_zero(z))
        return 0;

    /* (1 - 1/z is raised to the power -a in the integrals) */
    wp = prec + 10 + euler_log_bits(a, b, c);

    acb_init(A); acb_init(C); acb_init(Z); acb_init(t); acb_init(u);

    /* A = b - c + 1 (and Gamma(A)), C = a + A, Z = 1 - 1/z */
    acb_sub(A, b, c, wp);
    acb_add_ui(A, A, 1, wp);
    acb_add(C, A, a, wp);
    acb_inv(Z, z, wp);
    acb_sub_ui(Z, Z, 1, wp);
    acb_neg(Z, Z);

    ok = !(acb_contains_int(A) && !arb_is_positive(acb_realref(A)))
        && arb_is_nonpositive(acb_realref(a));

    ok = ok && _acb_hypgeom_2f1_two_integrals(t, A, b, C, Z, A, prec);

    if (ok)
    {
        real = acb_is_real(a) && acb_is_real(b) && acb_is_real(c) && acb_is_real(z);

        /* Gamma(c) Gamma(A) / Gamma(b) z^(-b) */
        acb_lgamma(u, A, wp);
        acb_lgamma(C, b, wp);
        acb_sub(u, u, C, wp);
        if (!(regularized & ACB_HYPGEOM_2F1_REGULARIZED))
        {
            acb_lgamma(C, c, wp);
            acb_add(u, u, C, wp);
        }
        acb_log(C, z, wp);
        acb_submul(u, C, b, wp);
        acb_exp(u, u, wp);
        acb_mul(t, t, u, wp);

        if (real)
            arb_zero(acb_imagref(t));

        ok = acb_is_finite(t);
        if (ok)
            acb_set(res, t);
    }

    acb_clear(A); acb_clear(C); acb_clear(Z); acb_clear(t); acb_clear(u);
    return ok;
}

/* The Euler integral for Re(z) > 1, where the singular point 1 / z of
   (1 - z t)^(-a) is on the segment [0, 1] or close to it: over the
   polygon from 0 to 1 by way of 1/2 - i Y for Im(z) <= 0, which passes
   below the point 1 / z (for z on the cut (1, inf) this is the
   principal value, the limit from the lower half-plane), and by way of
   1/2 + i Y for Im(z) > 0.  The term 1 - z t is given as -+i (1 - z t),
   whose principal logarithm is continuous along the polygon.  (For
   large real parameters the saddle points can both be real, with the
   line of descent of the first running into the second: the path is
   then planned with the singular point moved off the real axis.) */
static int
_acb_hypgeom_2f1_integration_cut(acb_t res, const acb_t a, const acb_t b, const acb_t c,
    const acb_t z, int regularized, slong prec)
{
    acb_struct mu[3], al[3], be[3];
    acb_t via, zero, t, u, J;
    slong wp;
    double sr, si, Y;
    int j, ok, lower;

    if (!acb_is_finite(a) || !acb_is_finite(b) || !acb_is_finite(c) || !acb_is_finite(z))
        return 0;

    wp = prec + 10 + euler_log_bits(a, b, c);

    /* (the side: z on the real line exactly, or off it) */
    if (arb_is_zero(acb_imagref(z)) || arb_is_negative(acb_imagref(z)))
        lower = 1;
    else if (arb_is_positive(acb_imagref(z)))
        lower = 0;
    else
        return 0;

    acb_init(via); acb_init(zero); acb_init(t); acb_init(u); acb_init(J);
    for (j = 0; j < 3; j++)
    {
        acb_init(mu + j);
        acb_init(al + j);
        acb_init(be + j);
    }

    /* t^(b-1) (1 - t)^(c-b-1) */
    acb_sub_ui(mu, b, 1, wp);
    acb_sub(mu + 1, c, b, wp);
    acb_sub_ui(mu + 1, mu + 1, 1, wp);
    acb_neg(mu + 2, a);
    acb_one(be);
    acb_one(al + 1);
    acb_set_si(be + 1, -1);

    ok = arb_is_nonnegative(acb_realref(mu)) && arb_is_nonnegative(acb_realref(mu + 1));

    /* the point 1 / z = sr + i si, with 0 < sr < 1, and the height of
       the polygon as above */
    acb_inv(u, z, 64);
    sr = arf_get_d(arb_midref(acb_realref(u)), ARF_RND_NEAR);
    si = fabs(arf_get_d(arb_midref(acb_imagref(u)), ARF_RND_NEAR));
    Y = 1.0;
    acb_sub_ui(t, z, 1, wp);
    if (!arb_is_positive(acb_realref(t)) || !(sr > 0.0 && sr < 1.0))
        ok = 0;
    else
    {
        while (Y < 1048576.0 && (Y * sr < 2.0 * si || Y * (1.0 - sr) < 2.0 * si))
            Y *= 2.0;
        if (Y >= 1048576.0)
            ok = 0;
    }

    if (ok)
    {
        /* -+i (1 - z t), by way of 1/2 -+ i Y */
        acb_onei(al + 2);
        if (lower)
            acb_neg(al + 2, al + 2);
        acb_mul(be + 2, al + 2, z, wp);
        acb_neg(be + 2, be + 2);
        acb_set_d_d(via, 0.5, lower ? -Y : Y);

        /* (the tilted plan first if 1 / z is on the segment or close
           to it) */
        ok = _acb_hypgeom_integrate_exp_log_via(J, via, zero, mu, al, be, 3, 0, 1, -INFINITY,
                si <= 1e-6 * FLINT_MIN(sr, 1.0 - sr), prec);
    }

    if (ok)
    {
        /* e^(-+ pi i a / 2) */
        acb_mul_2exp_si(t, a, -1);
        if (lower)
            acb_neg(t, t);
        acb_exp_pi_i(t, t, wp);
        acb_mul(J, J, t, wp);

        /* Gamma(c) / (Gamma(b) Gamma(c - b)) */
        acb_lgamma(u, b, wp);
        acb_neg(u, u);
        acb_sub(t, c, b, wp);
        acb_lgamma(t, t, wp);
        acb_sub(u, u, t, wp);
        if (!(regularized & ACB_HYPGEOM_2F1_REGULARIZED))
        {
            acb_lgamma(t, c, wp);
            acb_add(u, u, t, wp);
        }
        acb_exp(u, u, wp);
        acb_mul(J, J, u, wp);

        ok = acb_is_finite(J);
        if (ok)
            acb_set(res, J);
    }

    acb_clear(via); acb_clear(zero); acb_clear(t); acb_clear(u); acb_clear(J);
    for (j = 0; j < 3; j++)
    {
        acb_clear(mu + j);
        acb_clear(al + j);
        acb_clear(be + j);
    }

    return ok;
}

void
acb_hypgeom_2f1_integration(acb_t res, const acb_t a, const acb_t b, const acb_t c,
    const acb_t z, int regularized, slong prec)
{
    acb_t res2;
    acb_init(res2);

    /* The Euler integral for the larger one of a and b, or for the
       other; over a polygon past the point 1 / z for Re(z) > 1 (first
       if z is on the cut or close to it, where the segment [0, 1] meets
       that point or the plan of the path depends on rounding errors).
       For Re(z) > 1 also the integrals of

           F(a, b; c; z) = (1 - z)^(c-a-b) F(c - a, c - b; c; z),

       which exist for the same parameters, and those first whose
       integrand has a zero at 1 / z, not a pole: with large a > c and
       two real saddle points outside of [0, 1] the path of the first
       integral leads around the pole by way of both, while the other
       one is the sum of two integrals over real segments.  Then the
       two integrals for Re(z) < 0, and the loop. */
    {
        acb_struct par[4][2];
        acb_t ca, cb, t, u;
        int tr[4], zero[4], n = 0, i, j, ok = 0, near_cut, right;
        double zr, zi;
        const acb_struct * p, * q;
        slong wp = prec + 10 + euler_log_bits(a, b, c);

        acb_init(ca); acb_init(cb); acb_init(t); acb_init(u);

        if (arf_cmp(arb_midref(acb_realref(a)), arb_midref(acb_realref(b))) < 0)
        {
            p = a; q = b;
        }
        else
        {
            p = b; q = a;
        }

        zr = arf_get_d(arb_midref(acb_realref(z)), ARF_RND_NEAR);
        zi = arf_get_d(arb_midref(acb_imagref(z)), ARF_RND_NEAR);
        right = (zr > 1.0);
        near_cut = right && (fabs(zi) <= 1e-6 * (zr - 1.0));

        /* (the pairs: the parameter of the term 1 - z t, and the one of
           the integral) */
        par[n][0] = *p; par[n][1] = *q; tr[n++] = 0;
        par[n][0] = *q; par[n][1] = *p; tr[n++] = 0;

        if (right && acb_is_finite(a) && acb_is_finite(b) && acb_is_finite(c))
        {
            acb_sub(ca, c, a, wp + 200);
            acb_sub(cb, c, b, wp + 200);

            if (arf_cmp(arb_midref(acb_realref(ca)), arb_midref(acb_realref(cb))) < 0)
            {
                par[n][0] = *ca; par[n][1] = *cb; tr[n++] = 1;
                par[n][0] = *cb; par[n][1] = *ca; tr[n++] = 1;
            }
            else
            {
                par[n][0] = *cb; par[n][1] = *ca; tr[n++] = 1;
                par[n][0] = *ca; par[n][1] = *cb; tr[n++] = 1;
            }
        }

        for (i = 0; i < n; i++)
            zero[i] = right && arb_is_nonpositive(acb_realref(par[i]));

        for (j = 1; j >= 0 && !ok; j--)
        {
            for (i = 0; i < n && !ok; i++)
            {
                if (zero[i] != j)
                    continue;

                if (near_cut)
                    ok = _acb_hypgeom_2f1_integration_cut(res2, par[i], par[i] + 1, c, z, regularized, prec)
                        || _acb_hypgeom_2f1_integration(res2, par[i], par[i] + 1, c, z, regularized, prec);
                else
                    ok = _acb_hypgeom_2f1_integration(res2, par[i], par[i] + 1, c, z, regularized, prec)
                        || _acb_hypgeom_2f1_integration_cut(res2, par[i], par[i] + 1, c, z, regularized, prec);

                if (ok && tr[i])
                {
                    /* (1 - z)^(c-a-b); on the cut, 1 - z is negative
                       with the argument pi, the limit from below */
                    acb_sub(t, ca, b, wp);
                    acb_sub_ui(u, z, 1, wp);
                    acb_neg(u, u);
                    acb_pow(u, u, t, wp);
                    acb_mul(res2, res2, u, wp);
                    ok = acb_is_finite(res2);
                }
            }
        }

        ok = ok || _acb_hypgeom_2f1_integration_left(res2, p, q, c, z, regularized, prec)
            || _acb_hypgeom_2f1_integration_left(res2, q, p, c, z, regularized, prec)
            || _acb_hypgeom_2f1_integration_loop(res2, p, q, c, z, regularized, prec)
            || _acb_hypgeom_2f1_integration_loop(res2, q, p, c, z, regularized, prec);

        if (!ok)
            acb_indeterminate(res2);

        acb_clear(ca); acb_clear(cb); acb_clear(t); acb_clear(u);
    }

    acb_set_round(res, res2, prec);
    acb_clear(res2);
}
