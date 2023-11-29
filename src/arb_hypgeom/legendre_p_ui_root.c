/*
    Copyright (C) 2017 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb_hypgeom.h"

/* Compute initial isolating interval, following K. Petras. */
void
arb_hypgeom_legendre_p_ui_root_initial(arb_t res, ulong n, ulong k, slong prec)
{
    arb_t phi, psi, s, c, t, u;
    mag_t err, errc, errd;
    slong i, tol;

    arb_init(phi);
    arb_init(psi);
    arb_init(s);
    arb_init(c);
    arb_init(t);
    arb_init(u);
    mag_init(err);
    mag_init(errc);
    mag_init(errd);

    /* Petras numbering starts from 1 */
    k++;
    if (k > n / 2)
        flint_throw(FLINT_ERROR, "(%s)\n", __func__);

    tol = -prec;
    /* We need slightly higher precision since the Newton iteration
       *arithmetic* error bounds are not self-correcting. */
    prec = prec * 1.2 + 10;

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
    arb_get_mag_lower(errc, phi);
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

    arb_clear(phi);
    arb_clear(psi);
    arb_clear(s);
    arb_clear(c);
    arb_clear(t);
    arb_clear(u);
    mag_clear(err);
    mag_clear(errc);
    mag_clear(errd);
}

void
arb_hypgeom_legendre_p_ui_root(arb_t res, arb_t weight, ulong n, ulong k, slong prec)
{
    slong padding, initial_prec, step, wp;
    slong steps[FLINT_BITS];
    arb_t t, u, v, v0;
    mag_t err, err2, pb, p2b;
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

    arb_init(t);
    arb_init(u);
    arb_init(v);
    arb_init(v0);
    mag_init(err);
    mag_init(err2);
    mag_init(pb);
    mag_init(p2b);

    padding = 8 + 2 * FLINT_BIT_COUNT(n);
    initial_prec = 40 + padding;

    if (sign == 0)
    {
        arb_zero(res);
    }
    else if (initial_prec > prec / 2)
    {
        arb_hypgeom_legendre_p_ui_root_initial(res, n, k, prec + padding);
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
        arb_hypgeom_legendre_p_ui_root_initial(res, n, k, wp);
        step--;

        arb_mul(t, res, res, wp);
        arb_sub_ui(t, t, 1, wp);
        arb_hypgeom_legendre_p_ui_deriv_bound(pb, p2b, n, res, t);

        arb_set(v0, res);

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
            {
                /* flint_printf("unexpected Newton iteration failure...\n"); */
                break;
            }

            arb_set(res, v);
        }
    }

    if (weight != NULL)
    {
        wp = FLINT_MAX(prec, 40) + padding;
        arb_hypgeom_legendre_p_ui(NULL, t, n, res, wp);
        arb_mul(t, t, t, wp);
        arb_mul(u, res, res, wp);
        arb_sub_ui(u, u, 1, wp);
        arb_neg(u, u);
        arb_mul(t, t, u, wp);
        arb_ui_div(weight, 2, t, prec);
    }

    if (sign == -1)
        arb_neg(res, res);

    arb_set_round(res, res, prec);

    arb_clear(t);
    arb_clear(u);
    arb_clear(v);
    arb_clear(v0);
    mag_clear(err);
    mag_clear(err2);
    mag_clear(pb);
    mag_clear(p2b);
}

