/*
    Copyright (C) 2023 Jean Kieffer

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "acb_theta.h"

void
acb_theta_agm_ext_conv_rate(arf_t c1, arf_t c2, arf_t r, const arf_t eps,
                            const arf_t m, const arf_t M, slong prec)
{
    arb_t M0, m0, M1, u0, U0, u1, U1;
    arb_t t, c;

    arb_init(M0);
    arb_init(m0);
    arb_init(M1);
    arb_init(u0);
    arb_init(U0);
    arb_init(u1);
    arb_init(U1);
    arb_init(t);
    arb_init(c);

    /* Get convergence rate of regular Borchardt; set ui, Mi */
    acb_theta_agm_conv_rate(c1, r, eps, prec);
    arb_set_arf(M0, M);
    arb_set_arf(m0, m);
    
    arb_set_arf(U0, c1);
    arb_mul_arf(U0, U0, r, prec);
    arb_mul_arf(U1, U0, r, prec);
    arb_neg(u0, U0);
    arb_neg(u1, U1);
    
    arb_add_si(U0, U0, 1, prec);
    arb_add_si(u0, u0, 1, prec);
    arb_add_si(U1, U1, 1, prec);
    arb_add_si(u1, u1, 1, prec);
    arb_mul(M1, M0, U0, prec);
    arb_sqrt(M1, M1, prec);

    /* Get c such that |u_0^(n+1) - v_0^n t_0^n| <= x_n:= c r^(2^(n-1)) */
    arb_sqrt(c, M1, prec);
    arb_mul_2exp_si(c, c, -1);
    arb_mul_arf(c, c, r, prec);
    
    arb_mul(t, M0, U0, prec);
    arb_div(t, t, m0, prec);
    arb_sqrt(t, t, prec);
    arb_add(c, c, t, prec);

    arb_mul(c, c, c1, prec);
    arb_mul(c, c, U0, prec);
    arb_sqrt(t, u0, prec);
    arb_div(c, c, t, prec);

    /* Get c' such that x_n^2 + 2 x_n v_0^(n) t_0^(n) <= c' r^(2^(n-1)) */
    arb_mul(t, M1, U0, prec);
    arb_sqrt(t, t, prec);
    arb_mul_si(t, t, prec);
    arb_mul(t, t, c, prec);
    
    arb_sqr(c, c, prec);
    arb_mul_arf(c, c, r, prec);
    arb_add(c, c, t, prec);
    
    /* Get c2 such that |q_{n+1} - 1| <= c2 r^(2^(n-1)) */
    arb_div(c, c, u0);
    arb_div(c, c, u0);

    arb_set_arf(t, r);
    arb_sqr(t, t);
    arb_sub_si(t, t, 1, prec);
    arb_neg(t, t);
    arb_inv(t, t, prec);
    arb_mul(t, t, U1, prec);
    arb_div(t, t, u1, prec);
    arb_mul_arf(t, t, r, prec);
    arb_mul_arf(t, t, c1, prec);
    arb_add(c, c, t, prec);
    
    arb_get_ubound_arf(c2, c, prec);

    arb_clear(M0);
    arb_clear(m0);
    arb_init(M1);
    arb_init(u0);
    arb_init(U0);
    arb_init(u1);
    arb_init(U1);
    arb_clear(t);
    arb_clear(c);
}
