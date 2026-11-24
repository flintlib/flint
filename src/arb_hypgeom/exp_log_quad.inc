/*
    Copyright (C) 2026 Fredrik Johansson
    Developed using Claude Fable 5.2

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/* The Gauss-Legendre rule with n nodes for the integral of
   exp(lambda t + sum_j mu_j log(al_j + be_j t) - c) over [m - h, m + h]
   in dNb ball arithmetic (included from exp_log_integration.c for
   N = 1, 2, 3, 4), evaluated directly.  Returns 0 if a logarithm is
   not defined on the balls. */
static int
ELQ(direct)(arb_t res, const el_kernel * K, const arb_t c, const arb_t m, const arb_t h,
    const DB(struct) * x, const DB(struct) * wts, slong n, void * scratch)
{
    DB(t) M, H, C, S;
    DB(ptr) t, u, lg, ps, e1;
    slong k;
    int j;

    t = (DB(ptr)) scratch;
    u = t + n; lg = u + n; ps = lg + n; e1 = ps + n;

    DB(set_arb)(M, m);
    DB(set_arb)(H, h);
    for (k = 0; k < n; k++)
        t[k] = *M;
    DBV(vec_addmul_scalar)(t, x, n, H);

    DB(set_arb)(C, c);
    DB(neg)(C, C);
    for (k = 0; k < n; k++)
        ps[k] = *C;
    DB(set_arb)(C, K->lam);
    DBV(vec_addmul_scalar)(ps, t, n, C);

    for (j = 0; j < K->m; j++)
    {
        if (arb_is_zero(K->mu[j]))
            continue;

        DB(set_arb)(C, K->al[j]);
        for (k = 0; k < n; k++)
            u[k] = *C;
        DB(set_arb)(C, K->be[j]);
        DBV(vec_addmul_scalar)(u, t, n, C);
        /* (the status: zero on success) */
        if (DBV(vec_log)(lg, u, n) != 0)
            return 0;
        DB(set_arb)(C, K->mu[j]);
        DBV(vec_addmul_scalar)(ps, lg, n, C);
    }

    DBV(vec_exp)(e1, ps, n);
    DBV(vec_dot)(S, NULL, 0, wts, e1, n);
    DB(get_arb)(res, S);
    arb_mul(res, res, h, 53 * N + 10);
    return arb_is_finite(res);
}
