/*
    Copyright (C) 2026 Fredrik Johansson
    Developed using Claude Fable 5.2

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/* The Gauss-Legendre rule with n nodes for int_a^b exp(phi(w) - c) dw in
   dNb ball arithmetic (included from sd_integrate.c for N = 1, 2, 3, 4
   with DB, DBV, QUAD defined): with w = m + h t, m = (a + b)/2,
   h = (b - a)/2, the nodes x_k and the weights w_k (as balls),

       res = h sum_k w_k exp(alpha E(m + h x_k) - nu (m + h x_k) - c),

   E(w) = e^w - e^-w if sym and E(w) = e^w otherwise, everything in real
   vector operations (the complex scalar products with the parameters as
   four real ones, so that complex parameters cost the same as real
   ones). */

/* the parameters as real balls */
typedef struct
{
    DB(struct) ar, ai, nr, ni;
    int ai_zero, ni_zero;
}
QUAD(params);

static void
QUAD(set_params)(QUAD(params) * par, const acb_t alpha, const acb_t nu)
{
    DB(set_arb)(&par->ar, acb_realref(alpha));
    DB(set_arb)(&par->ai, acb_imagref(alpha));
    DB(set_arb)(&par->nr, acb_realref(nu));
    DB(set_arb)(&par->ni, acb_imagref(nu));
    par->ai_zero = DB(is_zero)(&par->ai);
    par->ni_zero = DB(is_zero)(&par->ni);
}

/* res = (b - a) / 2, vec = (a + b) / 2 + ((b - a) / 2) x */
static void
QUAD(affine)(arb_t res, DB(ptr) vec, double a, double b, const DB(struct) * x, slong n)
{
    DB(t) A, B, M, H;
    slong k;

    DB(set_d)(A, a);
    DB(set_d)(B, b);
    DB(add)(M, A, B);
    DB(mul_2exp_si)(M, M, -1);
    DB(sub)(H, B, A);
    DB(mul_2exp_si)(H, H, -1);
    for (k = 0; k < n; k++)
        vec[k] = *M;
    DBV(vec_addmul_scalar)(vec, x, n, H);
    DB(get_arb)(res, H);
}

static void
QUAD(panel)(acb_t res, const QUAD(params) * par, int sym, cd c, cd a, cd b,
    const DB(struct) * x, const DB(struct) * wts, slong n, void * scratch)
{
    DB(t) C, S;
    DB(ptr) xr, xi, e1, e2, sn, cs, sh, ch, t1, t2, pr, pi;
    acb_t sum;
    slong k;

    xr = (DB(ptr)) scratch;
    xi = xr + n; e1 = xi + n; e2 = e1 + n; sn = e2 + n; cs = sn + n;
    sh = cs + n; ch = sh + n; t1 = ch + n; t2 = t1 + n; pr = t2 + n; pi = pr + n;

    /* w_k = m + h x_k; x = Re w, y = Im w */
    QUAD(affine)(acb_realref(res), xr, a.re, b.re, x, n);
    QUAD(affine)(acb_imagref(res), xi, a.im, b.im, x, n);

    /* E = t1 + i t2 = (e^x -+ e^-x) cos y + i (e^x +- e^-x) sin y */
    DBV(vec_exp)(e1, xr, n);
    DBV(vec_sin_cos)(sn, cs, xi, n);
    if (sym)
    {
        for (k = 0; k < n; k++)
            DB(neg)(t1 + k, xr + k);
        DBV(vec_exp)(e2, t1, n);
        DBV(vec_sub)(sh, e1, e2, n);
        DBV(vec_add)(ch, e1, e2, n);
        DBV(vec_mul)(t1, sh, cs, n);
        DBV(vec_mul)(t2, ch, sn, n);
    }
    else
    {
        DBV(vec_mul)(t1, e1, cs, n);
        DBV(vec_mul)(t2, e1, sn, n);
    }

    /* psi = alpha E - nu w - c */
    DB(set_d)(C, -c.re);
    for (k = 0; k < n; k++)
        pr[k] = *C;
    DB(set_d)(C, -c.im);
    for (k = 0; k < n; k++)
        pi[k] = *C;
    DBV(vec_addmul_scalar)(pr, t1, n, &par->ar);
    DBV(vec_addmul_scalar)(pi, t2, n, &par->ar);
    if (!par->ai_zero)
    {
        DBV(vec_submul_scalar)(pr, t2, n, &par->ai);
        DBV(vec_addmul_scalar)(pi, t1, n, &par->ai);
    }
    DBV(vec_submul_scalar)(pr, xr, n, &par->nr);
    DBV(vec_submul_scalar)(pi, xi, n, &par->nr);
    if (!par->ni_zero)
    {
        DBV(vec_addmul_scalar)(pr, xi, n, &par->ni);
        DBV(vec_submul_scalar)(pi, xr, n, &par->ni);
    }

    /* exp(psi) = e^(Re psi) (cos(Im psi) + i sin(Im psi)) */
    DBV(vec_exp)(e1, pr, n);
    DBV(vec_sin_cos)(sn, cs, pi, n);
    DBV(vec_mul)(t1, e1, cs, n);
    DBV(vec_mul)(t2, e1, sn, n);

    /* the sum, times h */
    acb_init(sum);
    DBV(vec_dot)(S, NULL, 0, wts, t1, n);
    DB(get_arb)(acb_realref(sum), S);
    DBV(vec_dot)(S, NULL, 0, wts, t2, n);
    DB(get_arb)(acb_imagref(sum), S);
    acb_mul(res, res, sum, 53 * N + 10);
    acb_clear(sum);
}

/* the nodes and weights of the rule with n nodes, from the arb values
   (the nonnegative half in decreasing order, and their weights) */
static void
QUAD(nodes)(DB(struct) * x, DB(struct) * w, arb_srcptr xa, arb_srcptr wa, slong n)
{
    slong k, kk;
    for (k = 0; k < n; k++)
    {
        kk = (2 * k < n) ? k : n - 1 - k;
        DB(set_arb)(x + k, xa + kk);
        if (2 * k >= n)
            DB(neg)(x + k, x + k);
        DB(set_arb)(w + k, wa + kk);
    }
}
