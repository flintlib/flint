/*
    Copyright (C) 2017, 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <math.h>
#include "arb.h"
#include "acb.h"
#include "arb_hypgeom.h"
#include "impl.h"

/* The approximation x = cos(psi) of K. Petras to x_k, 0 <= k < n/2, as a
   ball for cos(psi) (without the approximation error), and the factors
   of the error bound |x - x_k| <= errc errd, errc = phi^2,
   errd = 22 y^4 (1 + 2 y^2), y = 1/((4k+3) pi); a Newton step from x
   squares errd [Pet1999]. */
void
_arb_hypgeom_legendre_p_ui_root_petras(arb_t res, mag_t errc, mag_t errd, ulong n, ulong k, slong prec)
{
    arb_t phi, psi, s, c, t, u;
    mag_t err;

    arb_init(phi);
    arb_init(psi);
    arb_init(s);
    arb_init(c);
    arb_init(t);
    arb_init(u);
    mag_init(err);

    /* Petras numbering starts from 1 */
    k++;
    if (k > n / 2)
        flint_throw(FLINT_ERROR, "(%s)\n", __func__);

    /* t = 4n+2 */
    arb_set_ui(t, n);
    arb_mul_2exp_si(t, t, 2);
    arb_add_ui(t, t, 2, prec);

    /* u = (4k-1) pi */
    arb_set_ui(u, k);
    arb_mul_2exp_si(u, u, 2);
    arb_sub_ui(u, u, 1, prec);
    arb_const_pi(phi, prec);
    arb_mul(u, u, phi, prec);

    /* phi = ((4k-1)/(4n+2)) pi */
    arb_div(phi, u, t, prec);

    /* errc = phi^2 */
    arb_get_mag(errc, phi);
    mag_mul(errc, errc, errc);

    /* errd = 22*y^4 * (1 + 2*y^2), y = 1/u */
    /* err = y^2 */
    arb_get_mag_lower(err, u);
    mag_one(errd);
    mag_div(err, errd, err);
    mag_mul(err, err, err);
    /* errd = 1+2y^2 */
    mag_mul_2exp_si(errd, err, 1);
    mag_add_ui(errd, errd, 1);
    /* err = y^4 */
    mag_mul(err, err, err);
    /* errd *= 22y^4 */
    mag_mul(errd, errd, err);
    mag_mul_ui(errd, errd, 22);

    /* s, c = sin(phi), cos(phi) */
    arb_sin_cos(s, c, phi, prec);

    /* psi = phi + 2 cos(phi) / (t^2 sin(phi)) (1 - 11/[t^2 sin(phi)^2]) */
    arb_mul(psi, t, s, prec);
    arb_mul(psi, psi, psi, prec);
    arb_ui_div(psi, 11, psi, prec);
    arb_sub_ui(psi, psi, 1, prec);
    arb_neg(psi, psi);

    /* (destroying t) */
    arb_mul(t, t, t, prec);
    arb_mul(t, t, s, prec);
    arb_div(t, c, t, prec);
    arb_mul_2exp_si(t, t, 1);
    arb_mul(psi, psi, t, prec);
    arb_add(psi, psi, phi, prec);

    arb_cos(res, psi, prec);

    arb_clear(phi);
    arb_clear(psi);
    arb_clear(s);
    arb_clear(c);
    arb_clear(t);
    arb_clear(u);
    mag_clear(err);
}

/* Petras' enclosure of x_k */
static void
_petras_enclosure(arb_t res, ulong n, ulong k, slong prec)
{
    mag_t errc, errd;
    mag_init(errc);
    mag_init(errd);
    _arb_hypgeom_legendre_p_ui_root_petras(res, errc, errd, n, k, prec);
    mag_mul(errc, errc, errd);
    arb_add_error_mag(res, errc);
    mag_clear(errc);
    mag_clear(errd);
}

/* An enclosure of x_k (0 <= k < n/2) that contains no other root: the
   enclosure of Petras, checked to lie in (0, 1) and to be disjoint from
   those of the neighbouring positive roots (the next root below the
   smallest positive one is 0 or negative).  Returns 0 if the check
   fails. */
int
_arb_hypgeom_legendre_p_ui_root_isolate(arb_t res, ulong n, ulong k, slong prec)
{
    arb_t t;
    int ok;

    arb_init(t);

    _petras_enclosure(res, n, k, prec);
    arb_one(t);
    ok = arb_is_positive(res) && arb_lt(res, t);

    if (ok && k > 0)
    {
        _petras_enclosure(t, n, k - 1, prec);
        ok = arb_lt(res, t);
    }

    if (ok && k + 1 < n / 2)
    {
        _petras_enclosure(t, n, k + 1, prec);
        ok = arb_gt(res, t);
    }

    arb_clear(t);
    return ok;
}

/* Bounds for the Taylor series of P_n at a real point m.  For z in the
   disk |z - m| <= r, the Laplace integral
       P_n(z) = (1/pi) int_0^pi (z + sqrt(z^2 - 1) cos phi)^n dphi
   gives |P_n(z)| <= M = max rho(z)^n, rho(z) = max |z +- i sqrt(1 - z^2)|
   (the integrand is a convex combination of the two endpoints), bounded
   here on the box enclosing the disk.  The Cauchy estimates
   |c_j| <= M / r^j give, for |t| <= tau = q r with q < 1/2,
       T0 = M q^(R+1) / (1 - q)
          >= sum_(j > R) |c_j| tau^j,
       T1 = (M / r) (R + 1) q^R / (1 - q)^2
          >= sum_(j > R) j |c_j| tau^(j-1),
   bounds for the tails of S(t) = sum c_j t^j and S'(t) beyond degree R.
   Returns 0 (and infinite bounds) if the box is not in the half-plane
   Re(1 - z^2) > 0 where the square root is analytic, or if q >= 1/2. */
int
_arb_hypgeom_legendre_p_ui_taylor_tail(mag_t T0, mag_t T1, ulong n,
    const arb_t m, const mag_t r, const mag_t tau, slong R)
{
    acb_t z, w, t;
    mag_t rho, rho2, M, q, u;
    slong tp = MAG_BITS + 2 * FLINT_BIT_COUNT(n) + 16;
    int ok;

    acb_init(z);
    acb_init(w);
    acb_init(t);
    mag_init(rho);
    mag_init(rho2);
    mag_init(M);
    mag_init(q);
    mag_init(u);

    arb_set(acb_realref(z), m);
    arb_add_error_mag(acb_realref(z), r);
    arb_zero(acb_imagref(z));
    arb_add_error_mag(acb_imagref(z), r);

    /* w = i sqrt(1 - z^2); near x = 1, 1 - z^2 is about 1/n^2 */
    acb_mul(w, z, z, tp);
    acb_neg(w, w);
    acb_add_ui(w, w, 1, tp);
    ok = arb_is_positive(acb_realref(w));
    acb_sqrt(w, w, tp);
    acb_mul_onei(w, w);

    acb_add(t, z, w, tp);
    acb_get_mag(rho, t);
    acb_sub(t, z, w, tp);
    acb_get_mag(rho2, t);
    mag_max(rho, rho, rho2);
    mag_pow_ui(M, rho, n);

    mag_div(q, tau, r);
    ok = ok && mag_cmp_2exp_si(q, -1) < 0;

    if (ok)
    {
        /* u <= 1 - q */
        mag_one(u);
        mag_sub_lower(u, u, q);

        mag_pow_ui(T0, q, R + 1);
        mag_mul(T0, T0, M);
        mag_div(T0, T0, u);

        mag_pow_ui(T1, q, R);
        mag_mul(T1, T1, M);
        mag_div(T1, T1, r);
        mag_mul_ui(T1, T1, R + 1);
        mag_div(T1, T1, u);
        mag_div(T1, T1, u);
    }
    else
    {
        mag_inf(T0);
        mag_inf(T1);
    }

    acb_clear(z);
    acb_clear(w);
    acb_clear(t);
    mag_clear(rho);
    mag_clear(rho2);
    mag_clear(M);
    mag_clear(q);
    mag_clear(u);

    return ok;
}

/* A disk radius for the Taylor bounds at m in (0, 1): about the distance
   1/n over which P_n changes by a bounded factor, and well inside the
   region away from the singular point 1.  Any positive value is valid;
   this only affects how tight the bounds are. */
void
_arb_hypgeom_legendre_p_ui_disk_radius(mag_t r, const arb_t m, ulong n)
{
    double md, sd;

    md = fabs(arf_get_d(arb_midref(m), ARF_RND_NEAR));
    md = FLINT_MIN(md, 1.0);
    sd = sqrt(1.0 - md * md);
    /* zero (and the bounds infinite) at |m| >= 1 */
    mag_set_d_lower(r, FLINT_MIN(sd / n, (1.0 - md) / 4.0) * 0.99);
}

/* Petras' approximation followed by Newton steps with the a priori error
   bound, for 0 <= k < n/2. */
static void
_arb_hypgeom_legendre_p_ui_root_initial(arb_t res, ulong n, ulong k, slong prec)
{
    arb_t t, u;
    mag_t err, errc, errd;
    slong i, tol;

    arb_init(t);
    arb_init(u);
    mag_init(err);
    mag_init(errc);
    mag_init(errd);

    tol = -prec;
    /* We need slightly higher precision since the Newton iteration
       *arithmetic* error bounds are not self-correcting. */
    prec = prec * 1.2 + 10;

    _arb_hypgeom_legendre_p_ui_root_petras(res, errc, errd, n, k, prec);
    mag_mul(err, errc, errd);

    for (i = 0; i < FLINT_BITS; i++)
    {
        if (mag_cmp_2exp_si(err, tol) < 0)
            break;

        arb_hypgeom_legendre_p_ui(t, u, n, res, prec);

        arb_div(t, t, u, prec);
        arb_sub(res, res, t, prec);
        mag_mul(errd, errd, errd);
        mag_mul(err, errc, errd);
    }

    arb_add_error_mag(res, err);

    arb_clear(t);
    arb_clear(u);
    mag_clear(err);
    mag_clear(errc);
    mag_clear(errd);
}

/* x_k for 0 <= k < n/2 by Newton iteration with evaluations of P_n:
   the above at low precision, then interval Newton steps at doubling
   precision with arb_hypgeom_legendre_p_ui. */
static void
_arb_hypgeom_legendre_p_ui_root_basecase(arb_t res, ulong n, ulong k, slong prec)
{
    slong padding, initial_prec, step, wp;
    slong steps[FLINT_BITS];
    arb_t t, u, v;
    mag_t err, pb, p2b;

    arb_init(t);
    arb_init(u);
    arb_init(v);
    mag_init(err);
    mag_init(pb);
    mag_init(p2b);

    padding = 8 + 2 * FLINT_BIT_COUNT(n);
    initial_prec = 40 + padding;

    if (initial_prec > prec / 2)
    {
        _arb_hypgeom_legendre_p_ui_root_initial(res, n, k, prec + padding);
    }
    else
    {
        step = 0;
        steps[step] = prec + padding;

        while (step < FLINT_BITS - 1 && (steps[step] / 2) > initial_prec)
        {
            steps[step + 1] = (steps[step] / 2);
            step++;
        }

        wp = steps[step] + padding;
        _arb_hypgeom_legendre_p_ui_root_initial(res, n, k, wp);
        step--;

        arb_mul(t, res, res, wp);
        arb_sub_ui(t, t, 1, wp);
        arb_hypgeom_legendre_p_ui_deriv_bound(pb, p2b, n, res, t);

        for ( ; step >= 0; step--)
        {
            wp = steps[step] + padding;

            /* Interval Newton update: mid(x) - f(mid(x)) / f'(x) */
            /* We compute f'(mid(x)) and use the bound on f'' to get f'(x) */
            arb_set(v, res);
            mag_mul(err, p2b, arb_radref(v));
            mag_zero(arb_radref(v));
            arb_hypgeom_legendre_p_ui(t, u, n, v, wp);
            arb_add_error_mag(u, err);
            arb_div(t, t, u, wp);
            arb_sub(v, v, t, wp);

            if (mag_cmp(arb_radref(v), arb_radref(res)) >= 0)
                break;

            arb_set(res, v);
        }
    }

    arb_clear(t);
    arb_clear(u);
    arb_clear(v);
    mag_clear(err);
    mag_clear(pb);
    mag_clear(p2b);
}

/* The weight 2 / ((1 - x^2) P_n'(x)^2) at an enclosure x of a root */
static void
_arb_hypgeom_legendre_p_ui_root_weight(arb_t weight, ulong n, const arb_t x, slong prec)
{
    arb_t t, u;
    slong wp;

    arb_init(t);
    arb_init(u);

    wp = FLINT_MAX(prec, 40) + 8 + 2 * FLINT_BIT_COUNT(n);
    arb_hypgeom_legendre_p_ui(NULL, t, n, x, wp);
    arb_mul(t, t, t, wp);
    arb_mul(u, x, x, wp);
    arb_sub_ui(u, u, 1, wp);
    arb_neg(u, u);
    arb_mul(t, t, u, wp);
    arb_ui_div(weight, 2, t, prec);

    arb_clear(t);
    arb_clear(u);
}

/* Methods for the positive roots: Newton iteration with evaluations of
   P_n (small n or precision, and the fallback), bit-burst Newton steps
   with Taylor transfers (precision comparable to n), and Newton-Taylor
   steps on the Taylor series of P_n at a lower-precision approximation
   (precision much larger or much smaller than n). */
#define METHOD_BASECASE 0
#define METHOD_BITBURST 1
#define METHOD_TAYLOR 2

static int
_arb_hypgeom_legendre_p_ui_root_method(ulong n, slong prec)
{
    if (n < ARB_HYPGEOM_LEGENDRE_ROOT_MIN_N || n > ARB_HYPGEOM_LEGENDRE_ROOT_MAX_N
        || prec < ARB_HYPGEOM_LEGENDRE_ROOT_MIN_PREC)
        return METHOD_BASECASE;

    if ((double) prec >= 0.3 * n && (double) prec <= 15.0 * n)
        return METHOD_BITBURST;

    return METHOD_TAYLOR;
}

void
_arb_hypgeom_legendre_p_ui_root_enclosure(arb_t res, arb_t weight, ulong n, ulong k, slong prec)
{
    int method, ok = 0;

    method = _arb_hypgeom_legendre_p_ui_root_method(n, prec);

    if (method == METHOD_BITBURST)
        ok = _arb_hypgeom_legendre_p_ui_root_bitburst(res, weight, n, k, prec);
    else if (method == METHOD_TAYLOR)
        ok = _arb_hypgeom_legendre_p_ui_root_taylor(res, weight, n, k, prec);

    if (!ok)
    {
        _arb_hypgeom_legendre_p_ui_root_basecase(res, n, k, prec);
        if (weight != NULL)
            _arb_hypgeom_legendre_p_ui_root_weight(weight, n, res, prec);
    }
}

void
arb_hypgeom_legendre_p_ui_root(arb_t res, arb_t weight, ulong n, ulong k, slong prec)
{
    int sign;

    if (k >= n)
    {
        flint_throw(FLINT_ERROR, "require n > 0 and a root index 0 <= k < n\n");
    }

    sign = 1;

    if (n % 2 == 1 && k == n / 2)
    {
        sign = 0;
    }
    else if (k >= n / 2)
    {
        k = n - k - 1;
        sign = -1;
    }

    if (sign == 0)
    {
        arb_zero(res);
        if (weight != NULL)
            _arb_hypgeom_legendre_p_ui_root_weight(weight, n, res, prec);
        return;
    }

    _arb_hypgeom_legendre_p_ui_root_enclosure(res, weight, n, k, prec);

    if (sign == -1)
        arb_neg(res, res);

    arb_set_round(res, res, prec);
}
