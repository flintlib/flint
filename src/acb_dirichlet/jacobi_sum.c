/*
    Copyright (C) 2016 Pascal Molin

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb_dirichlet.h"

/* J_N(1,a) = sum on x = 1 mod some p | q */
ulong
jacobi_one_prime(ulong p, ulong e, ulong pe, ulong cond)
{
    if (e > 1 && cond % (p*p) == 0)
    {
        return 0;
    }
    else
    {
        slong r = (e > 1) ? pe / p : 1;
        if (cond % p)
            return r * (p - 2);
        else
            return -r;
    }
}

static ulong
jacobi_one(const dirichlet_group_t G, ulong cond)
{
    slong k, r = 1;

    for (k = 0; k < G->num; k++)
        r *= jacobi_one_prime(G->P[k].p, G->P[k].e,
                G->P[k].pe.n, cond);
    return r;
}

void
acb_dirichlet_jacobi_sum(acb_t res, const dirichlet_group_t G, const dirichlet_char_t chi1, const dirichlet_char_t chi2, slong prec)
{
    if (G->q_even > 1)
    {
        acb_zero(res);
    }
    else if (chi1->n == 1 || chi2->n == 1)
    {
        ulong cond = (chi1->n == 1) ? dirichlet_conductor_char(G, chi2) : dirichlet_conductor_char(G, chi1);
        acb_set_si(res, jacobi_one(G, cond));
    }
    else if (nmod_mul(chi1->n, chi2->n, G->mod) == 1)
    {
        ulong n;
        n = jacobi_one(G, dirichlet_conductor_char(G, chi1));
        if (dirichlet_parity_char(G, chi1))
            acb_set_si(res, -n);
        else
            acb_set_si(res, n);
    }
    else
    {
        if (G->q <= 150)
            acb_dirichlet_jacobi_sum_naive(res, G, chi1, chi2, prec);
        else if (G->num > 1)
            acb_dirichlet_jacobi_sum_factor(res, G, chi1, chi2, prec);
        else if (G->P[0].e > 1)
            acb_dirichlet_jacobi_sum_naive(res, G, chi1, chi2, prec);
        else
            acb_dirichlet_jacobi_sum_gauss(res, G, chi1, chi2, prec);
    }
}

void
acb_dirichlet_jacobi_sum_ui(acb_t res, const dirichlet_group_t G, ulong a, ulong b, slong prec)
{
    if (G->q_even > 1)
    {
        acb_zero(res);
    }
    else if (a == 1 || b == 1)
    {
        ulong cond = (a == 1) ? dirichlet_conductor_ui(G, b) : dirichlet_conductor_ui(G, a);
        acb_set_si(res, jacobi_one(G, cond));
    }
    else if (nmod_mul(a, b, G->mod) == 1)
    {
        ulong n;
        n = jacobi_one(G, dirichlet_conductor_ui(G, a));
        if (dirichlet_parity_ui(G, a))
            acb_set_si(res, -n);
        else
            acb_set_si(res, n);
    }
    else
    {
        dirichlet_char_t chi1, chi2;
        dirichlet_char_init(chi1, G);
        dirichlet_char_init(chi2, G);
        dirichlet_char_log(chi1, G, a);
        dirichlet_char_log(chi2, G, b);
        acb_dirichlet_jacobi_sum(res, G, chi1, chi2, prec);
        dirichlet_char_clear(chi1);
        dirichlet_char_clear(chi2);
    }
}
