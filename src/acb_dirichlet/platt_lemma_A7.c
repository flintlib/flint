/*
    Copyright (C) 2019 D.H.J Polymath

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb.h"
#include "acb_dirichlet.h"
#include "acb_dirichlet/impl.h"
#include "arb_poly.h"
#include "arb_hypgeom.h"

static void
_platt_lemma_A7_S(arb_t out, slong sigma,
        const arb_t t0, const arb_t h, slong k, slong A, slong prec)
{
    slong l;
    arb_t total, summand;
    arb_t pi, half;
    arb_t a;
    arb_t l_factorial, kd2, t02;
    arb_t x1, x2, x3, x4, x5;

    arb_init(total);
    arb_init(summand);
    arb_init(pi);
    arb_init(half);
    arb_init(a);
    arb_init(l_factorial);
    arb_init(kd2);
    arb_init(t02);
    arb_init(x1);
    arb_init(x2);
    arb_init(x3);
    arb_init(x4);
    arb_init(x5);

    arb_one(half);
    arb_mul_2exp_si(half, half, -1);
    arb_const_pi(pi, prec);
    arb_one(l_factorial);
    arb_set_si(kd2, k);
    arb_mul_2exp_si(kd2, kd2, -1);
    arb_sqr(t02, t0, prec);

    for (l=0; l<=(sigma-1)/2; l++)
    {
        if (l > 1)
        {
            arb_mul_si(l_factorial, l_factorial, l, prec);
        }

        arb_mul_si(a, pi, 4*l+1, prec);
        arb_mul_si(a, a, A, prec);

        arb_inv(x1, a, prec);
        arb_add_ui(x1, x1, 1, prec);

        arb_add_si(x2, half, 2*l, prec);
        arb_sqr(x2, x2, prec);
        arb_add(x2, x2, t02, prec);
        arb_pow(x2, x2, kd2, prec);
        arb_div(x2, x2, l_factorial, prec);

        arb_set_si(x3, 4*l + 1);
        arb_div(x3, x3, h, prec);
        arb_sqr(x3, x3, prec);
        arb_mul_2exp_si(x3, x3, -3);

        arb_mul_2exp_si(x4, a, -1);

        arb_sub(x5, x3, x4, prec);
        arb_exp(x5, x5, prec);

        arb_mul(summand, x1, x2, prec);
        arb_mul(summand, summand, x5, prec);

        arb_add(total, total, summand, prec);
    }

    arb_set(out, total);

    arb_clear(total);
    arb_clear(summand);
    arb_clear(pi);
    arb_clear(half);
    arb_clear(a);
    arb_clear(l_factorial);
    arb_clear(kd2);
    arb_clear(t02);
    arb_clear(x1);
    arb_clear(x2);
    arb_clear(x3);
    arb_clear(x4);
    arb_clear(x5);
}

void
acb_dirichlet_platt_lemma_A7(arb_t out, slong sigma,
        const arb_t t0, const arb_t h, slong k, slong A, slong prec)
{
    arb_t S, C;
    arb_t pi, a;
    arb_t x1, x2;
    arb_t y1, y2, y3, y4;
    arb_t z1, z2;

    if (sigma % 2 == 0 || sigma < 3)
    {
        arb_zero_pm_inf(out);
        return;
    }

    arb_init(S);
    arb_init(C);
    arb_init(pi);
    arb_init(a);
    arb_init(x1);
    arb_init(x2);
    arb_init(y1);
    arb_init(y2);
    arb_init(y3);
    arb_init(y4);
    arb_init(z1);
    arb_init(z2);

    arb_const_pi(pi, prec);

    arb_pow_ui(x1, pi, (ulong) (k+1), prec);
    arb_mul_2exp_si(x1, x1, k+3);

    arb_div(x2, t0, h, prec);
    arb_sqr(x2, x2, prec);
    arb_mul_2exp_si(x2, x2, -1);
    arb_neg(x2, x2);
    arb_exp(x2, x2, prec);

    _platt_lemma_A7_S(S, sigma, t0, h, k, A, prec);

    arb_mul(z1, x1, x2, prec);
    arb_mul(z1, z1, S, prec);

    arb_mul_si(a, pi, 2*sigma-1, prec);
    arb_mul_si(a, a, A, prec);

    arb_inv(y1, a, prec);
    arb_add_ui(y1, y1, 1, prec);

    arb_set_si(y2, 2*sigma + 1);
    arb_div(y2, y2, h, prec);
    arb_sqr(y2, y2, prec);
    arb_mul_2exp_si(y2, y2, -3);

    arb_mul_2exp_si(y3, a, -1);

    arb_sub(y4, y2, y3, prec);
    arb_exp(y4, y4, prec);

    acb_dirichlet_platt_c_bound(C, sigma, t0, h, k, prec);

    arb_mul(z2, y1, y4, prec);
    arb_mul(z2, z2, C, prec);
    arb_mul_2exp_si(z2, z2, 1);

    arb_add(out, z1, z2, prec);

    arb_clear(S);
    arb_clear(C);
    arb_clear(pi);
    arb_clear(a);
    arb_clear(x1);
    arb_clear(x2);
    arb_clear(y1);
    arb_clear(y2);
    arb_clear(y3);
    arb_clear(y4);
    arb_clear(z1);
    arb_clear(z2);
}


/* acb_dirichlet_platt_c_bound for k = 0, ..., K - 1: the same formula
   (see platt_c_bound.c), with the incomplete gamma values
   Gamma((k+l+1)/2, (sigma + 1/2)^2 / (2 h^2)), which depend on k + l
   only, computed once for all k (K + (sigma-1)/2 values instead of
   K ((sigma-1)/2 + 1)) */
static void
_platt_c_bound_vec(arb_ptr res, slong sigma, const arb_t t0, const arb_t h,
    slong K, slong prec)
{
    slong len = (sigma - 1)/2 + 1, k, l, m;
    arb_ptr G, p, cl;
    arb_t lhs, pi, two, u, e, v, base, b, x1, x2, Xa, Xb, f, t;

    if (!arb_is_positive(h))
        flint_throw(FLINT_ERROR, "requires positive h\n");

    arb_init(lhs);
    arb_set_si(lhs, 2*sigma + 1);
    arb_mul_2exp_si(lhs, lhs, -1);          /* sigma + 1/2, exact */
    if (!arb_le(lhs, t0))
    {
        for (k = 0; k < K; k++)
            arb_zero_pm_inf(res + k);
        arb_clear(lhs);
        return;
    }

    G = _arb_vec_init(K + len);
    p = _arb_vec_init(len);
    cl = _arb_vec_init(len);
    arb_init(pi); arb_init(two); arb_init(u); arb_init(e); arb_init(v);
    arb_init(base); arb_init(b); arb_init(x1); arb_init(x2); arb_init(Xa);
    arb_init(Xb); arb_init(f); arb_init(t);

    arb_const_pi(pi, prec);
    arb_set_si(two, 2);

    /* u = (sigma + 1/2)^2 / (2 h^2); G[m] = Gamma(m/2, u), 1 <= m < K + len */
    arb_div(u, lhs, h, prec);
    arb_sqr(u, u, prec);
    arb_mul_2exp_si(u, u, -1);
    for (m = 1; m < K + len; m++)
    {
        arb_set_si(x1, m);
        arb_mul_2exp_si(x1, x1, -1);
        arb_hypgeom_gamma_upper(G + m, x1, u, 0, prec);
    }

    /* e = exp((1 + sqrt(8))/(6 t0)), v = (sigma + 1/2 + t0)^((sigma - 1)/2) */
    arb_sqrt_ui(e, 8, prec);
    arb_add_ui(e, e, 1, prec);
    arb_div_ui(e, e, 6, prec);
    arb_div(e, e, t0, prec);
    arb_exp(e, e, prec);
    arb_add(v, t0, lhs, prec);
    arb_pow_ui(v, v, (sigma - 1)/2, prec);

    /* cl[l] = binom(len - 1, l) (sqrt(2) h)^l */
    arb_sqrt_ui(base, 2, prec);
    arb_mul(base, base, h, prec);
    arb_one(b);
    for (l = 0; l < len; l++)
    {
        arb_bin_uiui(x1, len - 1, l, prec);
        arb_mul(cl + l, x1, b, prec);
        arb_mul(b, b, base, prec);
    }

    for (k = 0; k < K; k++)
    {
        /* Xa = 2^((6k+5-sigma)/4) pi^k (sigma + 1/2)^k h */
        arb_set_si(x1, 6*k + 5 - sigma);
        arb_mul_2exp_si(x1, x1, -2);
        arb_pow(x1, two, x1, prec);
        arb_mul(x2, lhs, pi, prec);
        arb_pow_ui(x2, x2, k, prec);
        arb_mul(Xa, x1, x2, prec);
        arb_mul(Xa, Xa, h, prec);

        /* Xb = 2^((6k+7-sigma)/4) pi^(k - 1/2) */
        arb_set_si(x1, 6*k + 7 - sigma);
        arb_mul_2exp_si(x1, x1, -2);
        arb_pow(x1, two, x1, prec);
        arb_set_si(x2, 2*k - 1);
        arb_mul_2exp_si(x2, x2, -1);        /* k - 1/2, exact */
        arb_pow(x2, pi, x2, prec);
        arb_mul(Xb, x1, x2, prec);

        /* p[len-1-l] = binom(len-1, l) 2^((k-1)/2) h^(k+1) (sqrt(2) h)^l
           Gamma((k+l+1)/2, u) */
        arb_set_si(f, k - 1);
        arb_mul_2exp_si(f, f, -1);
        arb_pow(f, two, f, prec);
        arb_pow_ui(t, h, k + 1, prec);
        arb_mul(f, f, t, prec);
        for (l = 0; l < len; l++)
        {
            arb_mul(x1, cl + l, f, prec);
            arb_mul(p + (len - 1 - l), x1, G + (k + l + 1), prec);
        }

        /* res = e (v Xa + Xb p(t0)) */
        _arb_poly_evaluate(res + k, p, len, t0, prec);
        arb_mul(res + k, res + k, Xb, prec);
        arb_addmul(res + k, v, Xa, prec);
        arb_mul(res + k, res + k, e, prec);
    }

    _arb_vec_clear(G, K + len);
    _arb_vec_clear(p, len);
    _arb_vec_clear(cl, len);
    arb_clear(lhs); arb_clear(pi); arb_clear(two); arb_clear(u); arb_clear(e);
    arb_clear(v); arb_clear(base); arb_clear(b); arb_clear(x1); arb_clear(x2);
    arb_clear(Xa); arb_clear(Xb); arb_clear(f); arb_clear(t);
}

/* The bounds of acb_dirichlet_platt_lemma_A7 for k = 0, ..., K - 1 at
   once: in the sum over l, the factor ((1/2 + 2l)^2 + t0^2)^(k/2) is
   the k-th power of one square root, and the rest of the summand does
   not depend on k, so the K sums cost one pass over l with K
   multiplications per l instead of K passes with a power, an
   exponential and an inversion per term. */
void
_acb_dirichlet_platt_lemma_A7_vec(arb_ptr out, slong sigma,
        const arb_t t0, const arb_t h, slong K, slong A, slong prec)
{
    slong l, k;
    arb_ptr S, pw, Cv;
    arb_t pi, half, a, l_factorial, t02, x1, x3, x4, x5, w, q;
    arb_t z1, x2, C, y1, y2, y3, y4, z2, p;

    if (K <= 0)
        return;

    if (sigma % 2 == 0 || sigma < 3)
    {
        for (k = 0; k < K; k++)
            arb_zero_pm_inf(out + k);
        return;
    }

    S = _arb_vec_init(K);
    pw = _arb_vec_init(K);
    Cv = _arb_vec_init(K);
    arb_init(pi); arb_init(half); arb_init(a); arb_init(l_factorial);
    arb_init(t02); arb_init(x1); arb_init(x3); arb_init(x4); arb_init(x5);
    arb_init(w); arb_init(q); arb_init(z1); arb_init(x2); arb_init(C);
    arb_init(y1); arb_init(y2); arb_init(y3); arb_init(y4); arb_init(z2);
    arb_init(p);

    arb_one(half);
    arb_mul_2exp_si(half, half, -1);
    arb_const_pi(pi, prec);
    arb_one(l_factorial);
    arb_sqr(t02, t0, prec);

    /* S_k = sum_l w_l q_l^k, q_l = sqrt((1/2 + 2l)^2 + t0^2),
       w_l = (1 + 1/a) exp((4l+1)^2/(8h^2) - a/2) / l!, a = pi (4l+1) A */
    for (l = 0; l <= (sigma - 1) / 2; l++)
    {
        if (l > 1)
            arb_mul_si(l_factorial, l_factorial, l, prec);

        arb_mul_si(a, pi, 4*l+1, prec);
        arb_mul_si(a, a, A, prec);

        arb_inv(x1, a, prec);
        arb_add_ui(x1, x1, 1, prec);

        arb_set_si(x3, 4*l + 1);
        arb_div(x3, x3, h, prec);
        arb_sqr(x3, x3, prec);
        arb_mul_2exp_si(x3, x3, -3);
        arb_mul_2exp_si(x4, a, -1);
        arb_sub(x5, x3, x4, prec);
        arb_exp(x5, x5, prec);

        arb_mul(w, x1, x5, prec);
        arb_div(w, w, l_factorial, prec);

        arb_add_si(q, half, 2*l, prec);
        arb_sqr(q, q, prec);
        arb_add(q, q, t02, prec);
        arb_sqrt(q, q, prec);

        arb_set(p, w);
        for (k = 0; k < K; k++)
        {
            arb_add(S + k, S + k, p, prec);
            if (k + 1 < K)
                arb_mul(p, p, q, prec);
        }
    }

    /* exp(-(t0/h)^2 / 2), and the second part as in the scalar version */
    arb_div(x2, t0, h, prec);
    arb_sqr(x2, x2, prec);
    arb_mul_2exp_si(x2, x2, -1);
    arb_neg(x2, x2);
    arb_exp(x2, x2, prec);

    arb_mul_si(a, pi, 2*sigma-1, prec);
    arb_mul_si(a, a, A, prec);
    arb_inv(y1, a, prec);
    arb_add_ui(y1, y1, 1, prec);
    arb_set_si(y2, 2*sigma + 1);
    arb_div(y2, y2, h, prec);
    arb_sqr(y2, y2, prec);
    arb_mul_2exp_si(y2, y2, -3);
    arb_mul_2exp_si(y3, a, -1);
    arb_sub(y4, y2, y3, prec);
    arb_exp(y4, y4, prec);
    arb_mul(y4, y4, y1, prec);

    /* the c bounds, all k at once */
    _platt_c_bound_vec(Cv, sigma, t0, h, K, prec);

    for (k = 0; k < K; k++)
    {
        /* 2^(k+3) pi^(k+1) exp(-(t0/h)^2/2) S_k */
        arb_pow_ui(z1, pi, (ulong) (k+1), prec);
        arb_mul_2exp_si(z1, z1, k+3);
        arb_mul(z1, z1, x2, prec);
        arb_mul(z1, z1, S + k, prec);

        arb_mul(z2, y4, Cv + k, prec);
        arb_mul_2exp_si(z2, z2, 1);

        arb_add(out + k, z1, z2, prec);
    }

    _arb_vec_clear(S, K);
    _arb_vec_clear(pw, K);
    _arb_vec_clear(Cv, K);
    arb_clear(pi); arb_clear(half); arb_clear(a); arb_clear(l_factorial);
    arb_clear(t02); arb_clear(x1); arb_clear(x3); arb_clear(x4); arb_clear(x5);
    arb_clear(w); arb_clear(q); arb_clear(z1); arb_clear(x2); arb_clear(C);
    arb_clear(y1); arb_clear(y2); arb_clear(y3); arb_clear(y4); arb_clear(z2);
    arb_clear(p);
}
