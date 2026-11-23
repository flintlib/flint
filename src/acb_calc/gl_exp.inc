/*
    Copyright (C) 2026 Fredrik Johansson
    Developed using Claude Fable 5.2

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/* Gauss-Legendre rules for integrals of exp(p(w)) in dNb ball arithmetic
   (included from gl_exp.c for N = 1, 2, 3, 4 with DB, DBV, GLE
   defined), everything in real vector operations. */

/* the nodes and weights of the rule with n nodes, from the arb values
   (the nonnegative half in decreasing order, and their weights) */
static void
GLE(nodes)(DB(struct) * x, DB(struct) * w, arb_srcptr xa, arb_srcptr wa, slong n)
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

/* res = (b - a) / 2, vec = (a + b) / 2 + ((b - a) / 2) x */
static void
GLE(affine)(arb_t res, DB(ptr) vec, double a, double b, const DB(struct) * x, slong n)
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

/* The exponent is the polynomial sum_{j < nc} c_j w^j, evaluated by the
   Horner scheme at the nodes of the chord [a, b] between two complex
   doubles. */
static void
GLE(exp_poly)(acb_t res, acb_srcptr coef, slong nc, double are, double aim,
    double bre, double bim, const DB(struct) * x, const DB(struct) * wts, slong n, void * scratch)
{
    DB(t) Cr, Ci, S;
    DB(ptr) xr, xi, e1, sn, cs, t1, t2, pr, pi, tp;
    acb_t sum;
    slong k, j;

    xr = (DB(ptr)) scratch;
    xi = xr + n; e1 = xi + n; sn = e1 + n; cs = sn + n;
    t1 = cs + n; t2 = t1 + n; pr = t2 + n; pi = pr + n;

    GLE(affine)(acb_realref(res), xr, are, bre, x, n);
    GLE(affine)(acb_imagref(res), xi, aim, bim, x, n);

    DB(set_arb)(Cr, acb_realref(coef + nc - 1));
    DB(set_arb)(Ci, acb_imagref(coef + nc - 1));
    for (k = 0; k < n; k++)
    {
        pr[k] = *Cr;
        pi[k] = *Ci;
    }

    for (j = nc - 2; j >= 0; j--)
    {
        /* (pr + i pi) (xr + i xi) + c_j */
        DBV(vec_mul)(t1, pr, xr, n);
        DBV(vec_mul)(t2, pi, xi, n);
        DBV(vec_sub)(e1, t1, t2, n);
        DBV(vec_mul)(t1, pr, xi, n);
        DBV(vec_mul)(t2, pi, xr, n);
        DBV(vec_add)(pi, t1, t2, n);
        tp = pr; pr = e1; e1 = tp;

        DB(set_arb)(Cr, acb_realref(coef + j));
        DB(set_arb)(Ci, acb_imagref(coef + j));
        for (k = 0; k < n; k++)
        {
            DB(add)(pr + k, pr + k, Cr);
            DB(add)(pi + k, pi + k, Ci);
        }
    }

    /* exp(p) = e^(Re p) (cos(Im p) + i sin(Im p)) */
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
