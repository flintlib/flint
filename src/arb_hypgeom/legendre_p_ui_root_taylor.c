/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb.h"
#include "arb_hypgeom.h"
#include "impl.h"

/*
    The root x_k of P_n from an enclosure X_1 = m +- tau at a fraction
    1/r of the precision (computed recursively by the best method):

    1. One evaluation of P_n(m) and P_n'(m) at the full precision with
       arb_hypgeom_legendre_p_ui, and the Taylor coefficients c_j of
       P_n(m + t) from the recurrence of the Legendre equation, with
       precision tapered by tau^j.
    2. The fixed-point iteration d <- -(c_0 + sum_(j>=2) c_j d^j) / c_1
       on midpoints, which gains about prec/r bits per step.
    3. An interval Newton step on the Taylor polynomial of degree R with
       the Cauchy bound for its tail (_arb_hypgeom_legendre_p_ui_taylor_tail),
       R chosen so that the tail is below 2^-prec at |t| <= tau, and
       containment in the isolating enclosure of Petras.

    This pays off when the precision is much larger than n (the single
    full-precision evaluation dominates) or much smaller (the
    evaluation uses the asymptotic expansion).
*/

#define ORDER 4
#define RMAX 64

int
_arb_hypgeom_legendre_p_ui_root_taylor(arb_t res, arb_t weight, ulong n, ulong k, slong prec)
{
    slong wp, R, j, it, step, pj, pt;
    ulong lambda;
    arb_t X0, X1, m, inv, ic1, d, t, s, D, N, Sd, Spd;
    arb_ptr c;
    mag_t r, tau, T0, T1, eps, aD;
    int ok = 0;

    if (n < 2 || k >= n / 2 || n > ARB_HYPGEOM_LEGENDRE_ROOT_MAX_N || prec < 64)
        return 0;

    lambda = n * (n + 1);
    wp = prec + 32 + FLINT_BIT_COUNT(n);
    step = prec / ORDER;

    arb_init(X0);
    arb_init(X1);
    arb_init(m);
    arb_init(inv);
    arb_init(ic1);
    arb_init(d);
    arb_init(t);
    arb_init(s);
    arb_init(D);
    arb_init(N);
    arb_init(Sd);
    arb_init(Spd);
    mag_init(r);
    mag_init(tau);
    mag_init(T0);
    mag_init(T1);
    mag_init(eps);
    mag_init(aD);
    c = _arb_vec_init(RMAX + 1);

    /* an enclosure of x_k containing no other root */
    if (!_arb_hypgeom_legendre_p_ui_root_isolate(X0, n, k, 64 + 2 * FLINT_BIT_COUNT(n)))
        goto cleanup;

    /* the starting enclosure m +- tau */
    _arb_hypgeom_legendre_p_ui_root_enclosure(X1, NULL, n, k, step + 32);
    arf_set(arb_midref(m), arb_midref(X1));
    mag_set(tau, arb_radref(X1));
    if (mag_is_zero(tau))
        mag_set_ui_2exp_si(tau, 1, -(step + 32));

    /* the degree: the tail at |t| <= tau below 2^-wp */
    _arb_hypgeom_legendre_p_ui_disk_radius(r, m, n);

    for (R = ORDER + 1; R <= RMAX; R++)
    {
        if (!_arb_hypgeom_legendre_p_ui_taylor_tail(T0, T1, n, m, r, tau, R))
            goto cleanup;
        if (mag_cmp_2exp_si(T0, -wp) < 0)
            break;
    }

    if (R > RMAX)
        goto cleanup;

    /* c_0 = P_n(m), c_1 = P_n'(m) */
    arb_hypgeom_legendre_p_ui(c + 0, c + 1, n, m, wp);

    /* (1 - m^2) (j+1)(j+2) c_(j+2) = 2 m (j+1)^2 c_(j+1) + (j(j+1) - lambda) c_j,
       c_j needed to about 2^-wp / tau^j */
    arb_mul(inv, m, m, wp);
    arb_sub_ui(inv, inv, 1, wp);
    arb_neg(inv, inv);
    arb_inv(inv, inv, wp);
    for (j = 0; j + 2 <= R; j++)
    {
        pj = FLINT_MAX(wp - (j + 1) * step, 64);
        arb_mul(t, m, c + j + 1, pj);
        arb_mul_ui(t, t, 2 * (ulong) (j + 1) * (j + 1), pj);
        arb_mul_si(s, c + j, (slong) ((ulong) j * (j + 1)) - (slong) lambda, pj);
        arb_add(t, t, s, pj);
        arb_mul(t, t, inv, pj);
        arb_div_ui(c + j + 2, t, (ulong) (j + 1) * (j + 2), pj);
    }

    /* the fixed-point iteration on midpoints */
    arb_inv(ic1, c + 1, wp);
    arb_mul(d, c + 0, ic1, wp);
    arb_neg(d, d);
    mag_zero(arb_radref(d));
    pt = wp - step + 16;
    for (it = 0; it < 2 * ORDER; it++)
    {
        arb_set_round(s, c + R, pt);
        for (j = R - 1; j >= 2; j--)
        {
            arb_mul(s, s, d, pt);
            arb_add(s, s, c + j, pt);
        }
        arb_mul(s, s, d, pt);
        arb_mul(s, s, d, pt);
        arb_add(s, s, c + 0, wp);
        arb_mul(s, s, ic1, wp);
        arb_neg(s, s);
        arb_sub(t, s, d, wp);
        arf_set(arb_midref(d), arb_midref(s));
        if (arf_is_zero(arb_midref(t)) || arf_abs_bound_lt_2exp_si(arb_midref(t)) < -wp + 8)
            break;
    }

    /* interval Newton step within D = d +- eps, |D| <= tau */
    mag_set_ui_2exp_si(eps, 1, -(prec + 24));
    arb_set(D, d);
    arb_add_error_mag(D, eps);
    arb_get_mag(aD, D);
    if (mag_cmp(aD, tau) > 0)
        goto cleanup;

    /* S(d) and S'(D), each with its tail */
    arb_set(Sd, c + R);
    for (j = R - 1; j >= 0; j--)
    {
        arb_mul(Sd, Sd, d, wp);
        arb_add(Sd, Sd, c + j, wp);
    }
    arb_add_error_mag(Sd, T0);

    arb_mul_ui(Spd, c + R, R, wp);
    for (j = R - 1; j >= 1; j--)
    {
        arb_mul(Spd, Spd, D, wp);
        arb_addmul_ui(Spd, c + j, j, wp);
    }
    arb_add_error_mag(Spd, T1);

    arb_div(N, Sd, Spd, wp);
    arb_sub(N, d, N, wp);
    if (!arb_contains(D, N))
        goto cleanup;

    arb_add(res, m, N, wp);
    if (!arb_contains(X0, res))
        goto cleanup;

    if (weight != NULL)
    {
        /* P_n'(x_k) = S'(N) +- T1 */
        arb_mul_ui(Spd, c + R, R, wp);
        for (j = R - 1; j >= 1; j--)
        {
            arb_mul(Spd, Spd, N, wp);
            arb_addmul_ui(Spd, c + j, j, wp);
        }
        arb_add_error_mag(Spd, T1);
        arb_mul(t, res, res, wp);
        arb_sub_ui(t, t, 1, wp);
        arb_neg(t, t);
        arb_mul(s, Spd, Spd, wp);
        arb_mul(s, s, t, wp);
        arb_ui_div(weight, 2, s, prec);
    }

    arb_set_round(res, res, prec);
    ok = 1;

cleanup:
    arb_clear(X0);
    arb_clear(X1);
    arb_clear(m);
    arb_clear(inv);
    arb_clear(ic1);
    arb_clear(d);
    arb_clear(t);
    arb_clear(s);
    arb_clear(D);
    arb_clear(N);
    arb_clear(Sd);
    arb_clear(Spd);
    mag_clear(r);
    mag_clear(tau);
    mag_clear(T0);
    mag_clear(T1);
    mag_clear(eps);
    mag_clear(aD);
    _arb_vec_clear(c, RMAX + 1);

    return ok;
}
