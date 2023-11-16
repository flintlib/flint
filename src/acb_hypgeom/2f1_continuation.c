/*
    Copyright (C) 2015 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb.h"
#include "acb_hypgeom.h"

/* Differential equation for F(a,b,c,y+z):

   (y+z)(y-1+z) F''(z) + ((y+z)(a+b+1) - c) F'(z) + a b F(z) = 0

   Coefficients in the Taylor series are bounded by

       A * binomial(N+k, k) * nu^k

   using the Cauchy-Kovalevskaya majorant method.
   See J. van der Hoeven, "Fast evaluation of holonomic functions near
   and in regular singularities"
*/
static void
bound(mag_t A, mag_t nu, mag_t N,
    const acb_t a, const acb_t b, const acb_t c, const acb_t y,
    const acb_t f0, const acb_t f1)
{
    mag_t M0, M1, t, u;
    acb_t d;

    acb_init(d);
    mag_init(M0);
    mag_init(M1);
    mag_init(t);
    mag_init(u);

    /* nu = max(1/|y-1|, 1/|y|) = 1/min(|y-1|, |y|) */
    acb_get_mag_lower(t, y);
    acb_sub_ui(d, y, 1, MAG_BITS);
    acb_get_mag_lower(u, d);
    mag_min(t, t, u);
    mag_one(u);
    mag_div(nu, u, t);

    /* M0 = 2 nu |ab| */
    acb_get_mag(t, a);
    acb_get_mag(u, b);
    mag_mul(M0, t, u);
    mag_mul(M0, M0, nu);
    mag_mul_2exp_si(M0, M0, 1);

    /* M1 = nu |a+b+1| + 2|c| */
    acb_add(d, a, b, MAG_BITS);
    acb_add_ui(d, d, 1, MAG_BITS);
    acb_get_mag(t, d);
    mag_mul(t, t, nu);

    acb_get_mag(u, c);
    mag_mul_2exp_si(u, u, 1);
    mag_add(M1, t, u);

    /* N = max(sqrt(2 M0), 2 M1) / nu */
    mag_mul_2exp_si(M0, M0, 1);
    mag_sqrt(M0, M0);
    mag_mul_2exp_si(M1, M1, 1);
    mag_max(N, M0, M1);
    mag_div(N, N, nu);

    /* A = max(|f0|, |f1| / (nu (N+1)) */
    acb_get_mag(t, f0);
    acb_get_mag(u, f1);
    mag_div(u, u, nu);
    mag_div(u, u, N);  /* upper bound for dividing by N+1 */
    mag_max(A, t, u);

    acb_clear(d);
    mag_clear(M0);
    mag_clear(M1);
    mag_clear(t);
    mag_clear(u);
}

/*
   F(x)  = c0   +     c1 x    +     c2 x^2   +     c3 x^3    +    [...]
   F'(x) = c1   +   2 c2 x    +   3 c3 x^2   +   4 c4 x^3    +    [...]
*/
static void
evaluate_sum(acb_t res, acb_t res1,
    const acb_t a, const acb_t b, const acb_t c, const acb_t y,
    const acb_t x, const acb_t f0, const acb_t f1, slong num, slong prec)
{
    acb_t s, s2, w, d, e, xpow, ck, cknext;
    slong k;

    acb_init(s);
    acb_init(s2);
    acb_init(w);
    acb_init(d);
    acb_init(e);
    acb_init(xpow);
    acb_init(ck);
    acb_init(cknext);

    /* d = (y-1)*y */
    acb_sub_ui(d, y, 1, prec);
    acb_mul(d, d, y, prec);
    acb_one(xpow);

    for (k = 0; k < num; k++)
    {
        if (k == 0)
        {
            acb_set(ck, f0);
            acb_set(cknext, f1);
        }
        else
        {
            acb_add_ui(w, b, k-1, prec);
            acb_mul(w, w, ck, prec);
            acb_add_ui(e, a, k-1, prec);
            acb_mul(w, w, e, prec);

            acb_add(e, a, b, prec);
            acb_add_ui(e, e, 2*(k+1)-3, prec);
            acb_mul(e, e, y, prec);
            acb_sub(e, e, c, prec);
            acb_sub_ui(e, e, k-1, prec);
            acb_mul_ui(e, e, k, prec);
            acb_addmul(w, e, cknext, prec);

            acb_mul_ui(e, d, k+1, prec);
            acb_mul_ui(e, e, k, prec);
            acb_div(w, w, e, prec);
            acb_neg(w, w);

            acb_set(ck, cknext);
            acb_set(cknext, w);
        }

        acb_addmul(s, ck, xpow, prec);
        acb_mul_ui(w, cknext, k+1, prec);
        acb_addmul(s2, w, xpow, prec);

        acb_mul(xpow, xpow, x, prec);
    }

    acb_set(res, s);
    acb_set(res1, s2);

    acb_clear(s);
    acb_clear(s2);
    acb_clear(w);
    acb_clear(d);
    acb_clear(e);
    acb_clear(xpow);
    acb_clear(ck);
    acb_clear(cknext);
}

void
acb_hypgeom_2f1_continuation(acb_t res, acb_t res1,
    const acb_t a, const acb_t b, const acb_t c, const acb_t y,
    const acb_t z, const acb_t f0, const acb_t f1, slong prec)
{
    mag_t A, nu, N, w, err, err1, R, T, goal;
    acb_t x;
    slong j, k;

    mag_init(A);
    mag_init(nu);
    mag_init(N);
    mag_init(err);
    mag_init(err1);
    mag_init(w);
    mag_init(R);
    mag_init(T);
    mag_init(goal);
    acb_init(x);

    bound(A, nu, N, a, b, c, y, f0, f1);

    acb_sub(x, z, y, prec);

    /* |T(k)| <= A * binomial(N+k, k) * nu^k * |x|^k */
    acb_get_mag(w, x);
    mag_mul(w, w, nu); /* w = nu |x| */
    mag_mul_2exp_si(goal, A, -prec-2);

    /* bound for T(0) */
    mag_set(T, A);
    mag_inf(R);

    for (k = 1; k < 100 * prec; k++)
    {
        /* T(k) = T(k) * R(k), R(k) = (N+k)/k * w = (1 + N/k) w */
        mag_div_ui(R, N, k);
        mag_add_ui(R, R, 1);
        mag_mul(R, R, w);

        /* T(k) */
        mag_mul(T, T, R);

        if (mag_cmp(T, goal) <= 0 && mag_cmp_2exp_si(R, 0) < 0)
            break;
    }

    /* T(k) [1 + R + R^2 + R^3 + ...] */
    mag_geom_series(err, R, 0);
    mag_mul(err, T, err);

    /* Now compute T, R for the derivative */
    /* Coefficients are A * (k+1) * binomial(N+k+1, k+1) */
    mag_add_ui(T, N, 1);
    mag_mul(T, T, A);
    mag_inf(R);

    for (j = 1; j <= k; j++)
    {
        mag_add_ui(R, N, k + 1);
        mag_div_ui(R, R, k);
        mag_mul(R, R, w);
        mag_mul(T, T, R);
    }

    mag_geom_series(err1, R, 0);
    mag_mul(err1, T, err1);

    if (mag_is_inf(err))
    {
        acb_indeterminate(res);
        acb_indeterminate(res1);
    }
    else
    {
        evaluate_sum(res, res1, a, b, c, y, x, f0, f1, k, prec);

        acb_add_error_mag(res, err);
        acb_add_error_mag(res1, err1);
    }

    mag_clear(A);
    mag_clear(nu);
    mag_clear(N);
    mag_clear(err);
    mag_clear(err1);
    mag_clear(w);
    mag_clear(R);
    mag_clear(T);
    mag_clear(goal);
    acb_clear(x);
}

