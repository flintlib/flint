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

/*

[S(k+1)] = [ R(k)  0   ] [S(k)]
[T(k+1)]   [ 1     1   ] [T(k)]

[S(k+1)] = [ P(k) / Q(k)  0   ] [S(k)]
[T(k+1)]   [ 1            1   ] [T(k)]


  1  [ P(k)         ]
---- [              ]
Q(k) [ Q(k)   Q(k)  ]

[[A2 0] [B2 C2]] . [[A1 0] [B1 C1]] = [[A1 A2 0] [A1 B2 + B1 C2   C1 C2]

A1 B2 + B1 B2 = B2 (A1 + B1) -- use to save time?

*/

static void
factor(acb_t A, acb_t tmp, acb_srcptr a, slong p, const acb_t z, slong k, slong prec)
{
    slong i;

    if (p == 0)
    {
        if (z == NULL)
            acb_one(A);
        else
            acb_set(A, z);
    }
    else
    {
        acb_add_ui(A, a, k, prec);

        for (i = 1; i < p; i++)
        {
            acb_add_ui(tmp, a + i, k, prec);
            acb_mul(A, A, tmp, prec);
        }

        if (z != NULL)
            acb_mul(A, A, z, prec);
    }
}

static void
bsplit(acb_t A1, acb_t B1, acb_t C1,
        acb_srcptr a, slong p,
        acb_srcptr b, slong q,
        const acb_t z,
        slong aa,
        slong bb,
        slong prec,
        int invz)
{
    if (bb - aa == 1)
    {
        factor(A1, B1, a, p, invz ? NULL : z, aa, prec);
        factor(C1, B1, b, q, invz ? z : NULL, aa, prec);
        /* acb_set(B1, C1);   but we skip this */
    }
    else
    {
        slong m;

        acb_t A2, B2, C2;

        acb_init(A2);
        acb_init(B2);
        acb_init(C2);

        m = aa + (bb - aa) / 2;

        bsplit(A1, B1, C1, a, p, b, q, z, aa, m, prec, invz);
        bsplit(A2, B2, C2, a, p, b, q, z, m, bb, prec, invz);

        if (bb - m == 1)  /* B2 = C2 */
        {
            if (m - aa == 1)
                acb_add(B2, A1, C1, prec);
            else
                acb_add(B2, A1, B1, prec);

            acb_mul(B1, B2, C2, prec);
        }
        else
        {
            if (m - aa == 1)
                acb_mul(B1, C1, C2, prec);
            else
                acb_mul(B1, B1, C2, prec);

            acb_addmul(B1, A1, B2, prec);
        }

        acb_mul(A1, A1, A2, prec);
        acb_mul(C1, C1, C2, prec);

        acb_clear(A2);
        acb_clear(B2);
        acb_clear(C2);
    }
}

void
acb_hypgeom_pfq_sum_bs(acb_t s, acb_t t,
    acb_srcptr a, slong p, acb_srcptr b, slong q, const acb_t z, slong n, slong prec)
{
    acb_t u, v, w, tmp;

    if (n < 4)
    {
        acb_hypgeom_pfq_sum_forward(s, t, a, p, b, q, z, n, prec);
        return;
    }

    acb_init(u);
    acb_init(v);
    acb_init(w);
    acb_init(tmp);

    /* we compute to n-1 instead of n to avoid dividing by 0 in the
       denominator when computing a hypergeometric polynomial
       that terminates right before a pole */
    bsplit(u, v, w, a, p, b, q, z, 0, n - 1, prec, 0);

    acb_add(s, u, v, prec); /* s = s + t */
    acb_div(s, s, w, prec);

    /* split off last factor */
    factor(t, tmp, a, p, z, n - 1, prec);
    acb_mul(u, u, t, prec);
    factor(t, tmp, b, q, NULL, n - 1, prec);
    acb_mul(w, w, t, prec);
    acb_div(t, u, w, prec);

    acb_clear(u);
    acb_clear(v);
    acb_clear(w);
    acb_clear(tmp);
}

void
acb_hypgeom_pfq_sum_bs_invz(acb_t s, acb_t t,
    acb_srcptr a, slong p, acb_srcptr b, slong q, const acb_t z, slong n, slong prec)
{
    acb_t u, v, w, tmp;

    if (n < 4)
    {
        acb_init(u);
        acb_inv(u, z, prec);
        acb_hypgeom_pfq_sum_forward(s, t, a, p, b, q, u, n, prec);
        acb_clear(u);
        return;
    }

    acb_init(u);
    acb_init(v);
    acb_init(w);
    acb_init(tmp);

    /* we compute to n-1 instead of n to avoid dividing by 0 in the
       denominator when computing a hypergeometric polynomial
       that terminates right before a pole */
    bsplit(u, v, w, a, p, b, q, z, 0, n - 1, prec, 1);

    acb_add(s, u, v, prec); /* s = s + t */
    acb_div(s, s, w, prec);

    /* split off last factor */
    factor(t, tmp, a, p, NULL, n - 1, prec);
    acb_mul(u, u, t, prec);
    factor(t, tmp, b, q, z, n - 1, prec);
    acb_mul(w, w, t, prec);
    acb_div(t, u, w, prec);

    acb_clear(u);
    acb_clear(v);
    acb_clear(w);
    acb_clear(tmp);
}

