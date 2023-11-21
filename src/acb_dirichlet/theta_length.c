/*
    Copyright (C) 2016 Pascal Molin

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb_dirichlet.h"
#include <math.h>
#define PI   3.14159265358
#define LOG2 0.69314718055

ulong
acb_dirichlet_theta_length_d(ulong q, double t, slong prec)
{
    double a, la;
    a = PI / (double)q * t * t;
    la = (a < .3) ? -log(2*a*(1-a)) : .8;
    la = ((double)prec * LOG2 + la) / a;
    return ceil(sqrt(la)+.5);
}

ulong
acb_dirichlet_theta_length(ulong q, const arb_t t, slong prec)
{
    double dt;
    ulong len;
    arf_t at;
    arf_init(at);
    arb_get_lbound_arf(at, t, 53);
    dt = arf_get_d(at, ARF_RND_DOWN);
    len = acb_dirichlet_theta_length_d(q, dt, prec);
    arf_clear(at);
    return len;
}

/* bound for sum_{k>n} k*exp(-a k^2) */
void
mag_tail_kexpk2_arb(mag_t res, const arb_t a, ulong n)
{
    mag_t m;
    mag_init(m);
    arb_get_mag_lower(m, a);
    /* a < 1/4 */
    if (mag_cmp_2exp_si(m, -2) <= 0)
    {
        mag_t c;
        mag_init(c);
        mag_mul_ui_lower(res, m, n*n-n+1);  /* todo: possible overflow */
        mag_expinv(res, res);
        /* c = 2a(1+2a) */
        mag_mul_2exp_si(m, m, 1);
        mag_one(c);
        mag_add_lower(c, m, c);
        mag_mul_lower(c, m, c);
        mag_div(res, res, c);
        mag_clear(c);
    }
    else
    {
        mag_mul_ui_lower(res, m, n*n-n-1);  /* todo: possible overflow */
        mag_expinv(res, res);
        mag_mul_ui(res, res, 2);
    }
    mag_clear(m);
}

/* a(t) = Pi / q * t^2 */
void
_acb_dirichlet_theta_argument_at_arb(arb_t xt, ulong q, const arb_t t, slong prec)
{
    arb_const_pi(xt, prec);
    arb_div_ui(xt, xt, q, prec);
    arb_mul(xt, xt, t, prec);
    arb_mul(xt, xt, t, prec);
}
