/*
    Copyright (C) 2026 Fredrik Johansson
    Developed using Claude Fable 5.2

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/* The sum of a Gauss-Legendre rule with n nodes x_k (real, in [-1, 1])
   and weights w_k in dNb ball arithmetic (included from
   exp_log_integration.c for N = 1, 2, 3, 4) for the integrand on a
   chord m + h x in the local form

       exp(E + (lambda h) x + sum_j mu_j log(1 + q_j x)),

   E = phi(m) - c, q_j = be_j h / (al_j + be_j m) with |q_j| < 1, so that
   1 + q_j x has a positive real part:

       log(1 + q x) = log((1 + Re(q) x)^2 + (Im(q) x)^2) / 2
                      + i atan(Im(q) x / (1 + Re(q) x)).

   Everything is done in real vector operations.  The sum is not
   multiplied by h.  Returns 0 if a logarithm is not defined on the
   balls. */
static int
CELQ(direct)(acb_t res, const acb_t E, const acb_t lh, acb_srcptr mu, acb_srcptr q, int nt,
    const DB(struct) * x, const DB(struct) * wts, slong n, void * scratch)
{
    DB(t) C, S;
    DB(ptr) er, ei, a, b, t1, t2, lg, th;
    slong k;
    int j;

    er = (DB(ptr)) scratch;
    ei = er + n; a = ei + n; b = a + n; t1 = b + n; t2 = t1 + n; lg = t2 + n; th = lg + n;

    DB(set_arb)(C, acb_realref(E));
    for (k = 0; k < n; k++)
        er[k] = *C;
    DB(set_arb)(C, acb_imagref(E));
    for (k = 0; k < n; k++)
        ei[k] = *C;

    DB(set_arb)(C, acb_realref(lh));
    DBV(vec_addmul_scalar)(er, x, n, C);
    DB(set_arb)(C, acb_imagref(lh));
    DBV(vec_addmul_scalar)(ei, x, n, C);

    for (j = 0; j < nt; j++)
    {
        int real_q = arb_is_zero(acb_imagref(q + j));

        /* a + b i = 1 + q x */
        DB(one)(C);
        for (k = 0; k < n; k++)
            a[k] = *C;
        DB(set_arb)(C, acb_realref(q + j));
        DBV(vec_addmul_scalar)(a, x, n, C);

        if (real_q)
        {
            /* log(a) */
            if (DBV(vec_log)(lg, a, n) != 0)
                return 0;
            DB(set_arb)(C, acb_realref(mu + j));
            DBV(vec_addmul_scalar)(er, lg, n, C);
            if (!arb_is_zero(acb_imagref(mu + j)))
            {
                DB(set_arb)(C, acb_imagref(mu + j));
                DBV(vec_addmul_scalar)(ei, lg, n, C);
            }
        }
        else
        {
            DB(set_arb)(C, acb_imagref(q + j));
            DBV(vec_mul_scalar)(b, x, n, C);

            /* lg = log(a^2 + b^2) (twice the real part), th = atan(b / a) */
            DBV(vec_mul)(t1, a, a, n);
            DBV(vec_mul)(t2, b, b, n);
            DBV(vec_add)(t1, t1, t2, n);
            if (DBV(vec_log)(lg, t1, n) != 0)
                return 0;
            DBV(vec_div)(t2, b, a, n);
            DBV(vec_atan)(th, t2, n);

            /* mu (lg / 2 + i th) */
            DB(set_arb)(C, acb_realref(mu + j));
            DBV(vec_addmul_scalar)(ei, th, n, C);
            DB(mul_2exp_si)(C, C, -1);
            DBV(vec_addmul_scalar)(er, lg, n, C);
            if (!arb_is_zero(acb_imagref(mu + j)))
            {
                DB(set_arb)(C, acb_imagref(mu + j));
                DBV(vec_submul_scalar)(er, th, n, C);
                DB(mul_2exp_si)(C, C, -1);
                DBV(vec_addmul_scalar)(ei, lg, n, C);
            }
        }
    }

    /* exp(er) (cos(ei) + i sin(ei)) */
    DBV(vec_exp)(t1, er, n);
    DBV(vec_sin_cos)(a, b, ei, n);
    DBV(vec_mul)(t2, t1, b, n);
    DBV(vec_dot)(S, NULL, 0, wts, t2, n);
    DB(get_arb)(acb_realref(res), S);
    DBV(vec_mul)(t2, t1, a, n);
    DBV(vec_dot)(S, NULL, 0, wts, t2, n);
    DB(get_arb)(acb_imagref(res), S);

    return acb_is_finite(res);
}
