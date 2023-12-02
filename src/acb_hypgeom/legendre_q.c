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

/* invalid in (-1,0) */
int
_acb_hypgeom_legendre_q_single_valid(const acb_t z)
{
    arb_t t;
    int ok;

    if (!arb_contains_zero(acb_imagref(z)))
        return 1;

    if (arb_is_positive(acb_imagref(z)))
        return 1;

    arb_init(t);
    arb_one(t);
    arb_neg(t, t);
    ok = arb_lt(acb_realref(z), t);
    arb_clear(t);
    return ok;
}

void
_acb_hypgeom_legendre_q_double(acb_t res, const acb_t n, const acb_t m,
    const acb_t z, slong prec)
{
    acb_t t, u, v;

    acb_init(t);
    acb_init(u);
    acb_init(v);

    if (acb_is_int(m))
    {
        acb_sub_ui(t, z, 1, prec);
        acb_mul_2exp_si(u, m, -1);
        acb_pow(v, t, u, prec);
        acb_neg(t, t);
        acb_neg(u, u);
        acb_pow(t, t, u, prec);
        acb_mul(t, t, v, prec);

        acb_hypgeom_legendre_q(u, n, m, z, 0, prec);
        acb_mul(t, t, u, prec);

        acb_mul_2exp_si(u, m, -1);
        if (!acb_is_int(u))
            acb_neg(t, t);

        acb_sub_ui(u, z, 1, prec);
        acb_sqrt(u, u, prec);
        acb_sub_ui(v, z, 1, prec);
        acb_neg(v, v);
        acb_rsqrt(v, v, prec);
        acb_mul(u, u, v, prec);
        acb_hypgeom_legendre_p(v, n, m, z, 1, prec);
        acb_mul(u, u, v, prec);
        acb_const_pi(v, prec);
        acb_mul(u, u, v, prec);
        acb_mul_2exp_si(u, u, -1);

        acb_sub(res, t, u, prec);
    }
    else
    {
        acb_sub(t, n, m, prec);
        acb_add_ui(t, t, 1, prec);
        acb_mul_2exp_si(u, m, 1);
        acb_rising(t, t, u, prec);
        acb_neg(u, m);
        acb_hypgeom_legendre_p(u, n, u, z, 1, prec);
        acb_mul(t, t, u, prec);

        acb_hypgeom_legendre_p(u, n, m, z, 1, prec);
        acb_sub(t, u, t, prec);

        acb_exp_pi_i(u, m, prec);
        acb_mul(t, t, u, prec);

        acb_sin_pi(u, m, prec);
        acb_div(t, t, u, prec);
        acb_const_pi(u, prec);
        acb_mul(t, t, u, prec);
        acb_mul_2exp_si(t, t, -1);

        acb_set(res, t);
    }

    acb_clear(t);
    acb_clear(u);
    acb_clear(v);
}

void
_acb_hypgeom_legendre_q_single(acb_t res, const acb_t n, const acb_t m,
    const acb_t z, slong prec)
{
    acb_t a, b, c, z2, t, u;

    acb_init(a);
    acb_init(b);
    acb_init(c);
    acb_init(z2);
    acb_init(t);
    acb_init(u);

    /* invalid in (-1,0) */
    if (!_acb_hypgeom_legendre_q_single_valid(z))
    {
        acb_indeterminate(res);
        return;
    }

    acb_pow_si(z2, z, -2, prec); /* z2 = 1/z^2 */

    /* t = 2F1r((m+n+1)/2, (m+n)/2+1, n+3/2, 1/z^2) */
    acb_add(b, m, n, prec);
    acb_add_ui(a, b, 1, prec);
    acb_mul_2exp_si(a, a, -1);
    acb_mul_2exp_si(b, b, -1);
    acb_add_ui(b, b, 1, prec);
    acb_set_ui(c, 3);
    acb_mul_2exp_si(c, c, -1);
    acb_add(c, c, n, prec);
    acb_hypgeom_2f1(t, a, b, c, z2, 1, prec);

    /* prefactor sqrt(pi) 2^-n (z+1)^(m/2) (z-1)^(m/2) exp(i pi m) */
    /*           (1/2) gamma(m+n+1) z^(-m-n-1) */
    if (!acb_is_zero(m))
    {
        acb_add_ui(z2, z, 1, prec);
        acb_mul_2exp_si(c, m, -1);
        acb_pow(z2, z2, c, prec);
        acb_mul(t, t, z2, prec);

        acb_sub_ui(z2, z, 1, prec);
        acb_mul_2exp_si(c, m, -1);
        acb_pow(z2, z2, c, prec);
        acb_mul(t, t, z2, prec);

        acb_exp_pi_i(z2, m, prec);
        acb_mul(t, t, z2, prec);
    }

    acb_set_ui(z2, 2);
    acb_neg(c, n);
    acb_pow(z2, z2, c, prec);
    acb_mul(t, t, z2, prec);

    acb_add(c, m, n, prec);
    acb_add_ui(c, c, 1, prec);
    acb_gamma(z2, c, prec);
    acb_mul(t, t, z2, prec);

    acb_neg(c, c);
    acb_pow(z2, z, c, prec);
    acb_mul(t, t, z2, prec);

    acb_mul_2exp_si(t, t, -1);

    arb_const_sqrt_pi(acb_realref(u), prec);
    acb_mul_arb(t, t, acb_realref(u), prec);

    acb_set(res, t);

    acb_clear(a);
    acb_clear(b);
    acb_clear(c);
    acb_clear(z2);
    acb_clear(t);
    acb_clear(u);
}

void
acb_hypgeom_legendre_q(acb_t res, const acb_t n, const acb_t m,
    const acb_t z, int type, slong prec)
{
    if (type == 0)
    {
        /* http://functions.wolfram.com/07.11.26.0033.01 */
        /* todo: simplify the gamma quotients and the sqrt pi factor... */
        acb_t a, b, c, z2, mn, nm, t, u;

        acb_init(a);
        acb_init(b);
        acb_init(c);
        acb_init(z2);
        acb_init(mn);
        acb_init(nm);
        acb_init(t);
        acb_init(u);

        acb_add(mn, m, n, prec); /* mn = m + n */
        acb_sub(nm, n, m, prec); /* nm = n - m */
        acb_mul(z2, z, z, prec); /* z2 = z^2 */

        /* t = 2F1((1-m-n)/2, (n-m)/2+1, 3/2, z^2) */
        acb_sub_ui(a, mn, 1, prec);
        acb_neg(a, a);
        acb_mul_2exp_si(a, a, -1);
        acb_mul_2exp_si(b, nm, -1);
        acb_add_ui(b, b, 1, prec);
        acb_set_ui(c, 3);
        acb_mul_2exp_si(c, c, -1);
        acb_hypgeom_2f1(t, a, b, c, z2, 0, prec);

        /* u = 2F1(-(m+n)/2, (n-m+1)/2, 1/2, z^2) */
        acb_neg(a, mn);
        acb_mul_2exp_si(a, a, -1);
        acb_add_ui(b, nm, 1, prec);
        acb_mul_2exp_si(b, b, -1);
        acb_one(c);
        acb_mul_2exp_si(c, c, -1);
        acb_hypgeom_2f1(u, a, b, c, z2, 0, prec);

        /* a = cospi((m+n)/2) gamma((m+n)/2+1) rgamma((n-m+1)/2) z */
        /* b = sinpi((m+n)/2) gamma((m+n+1)/2) rgamma((n-m)/2+1) / 2 */
        acb_mul_2exp_si(a, mn, -1);
        acb_sin_cos_pi(b, a, a, prec);

        acb_mul_2exp_si(c, mn, -1);
        acb_add_ui(c, c, 1, prec);
        acb_gamma(c, c, prec);
        acb_mul(a, a, c, prec);
        acb_add_ui(c, nm, 1, prec);
        acb_mul_2exp_si(c, c, -1);
        acb_rgamma(c, c, prec);
        acb_mul(a, a, c, prec);
        acb_mul(a, a, z, prec);

        acb_add_ui(c, mn, 1, prec);
        acb_mul_2exp_si(c, c, -1);
        acb_gamma(c, c, prec);
        acb_mul(b, b, c, prec);
        acb_mul_2exp_si(c, nm, -1);
        acb_add_ui(c, c, 1, prec);
        acb_rgamma(c, c, prec);
        acb_mul(b, b, c, prec);
        acb_mul_2exp_si(b, b, -1);

        /* at - bu */
        acb_mul(t, t, a, prec);
        acb_mul(u, u, b, prec);
        acb_sub(t, t, u, prec);

        /* prefactor sqrt(pi) 2^m (1-z^2)^(-m/2) */
        if (!acb_is_zero(m))
        {
            acb_sub_ui(u, z2, 1, prec);
            acb_neg(u, u);
            acb_neg(c, m);
            acb_mul_2exp_si(c, c, -1);
            acb_pow(u, u, c, prec);
            acb_set_ui(c, 2);
            acb_pow(c, c, m, prec);
            acb_mul(u, u, c, prec);
            acb_mul(t, t, u, prec);
        }

        arb_const_sqrt_pi(acb_realref(u), prec);
        acb_mul_arb(t, t, acb_realref(u), prec);

        acb_set(res, t);

        acb_clear(a);
        acb_clear(b);
        acb_clear(c);
        acb_clear(z2);
        acb_clear(mn);
        acb_clear(nm);
        acb_clear(t);
        acb_clear(u);
    }
    else if (type == 1)
    {
        if ((arf_cmpabs_2exp_si(arb_midref(acb_realref(z)), -2) < 0 &&
             arf_cmpabs_2exp_si(arb_midref(acb_imagref(z)), -2) < 0) ||
            !_acb_hypgeom_legendre_q_single_valid(z))
        {
            _acb_hypgeom_legendre_q_double(res, n, m, z, prec);
        }
        else
        {
            _acb_hypgeom_legendre_q_single(res, n, m, z, prec);
        }
    }
    else
    {
        flint_throw(FLINT_ERROR, "unsupported 'type' %d for legendre q\n", type);
    }
}

