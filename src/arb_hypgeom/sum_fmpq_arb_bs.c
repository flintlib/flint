/*
    Copyright (C) 2021 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb_hypgeom.h"

static void
factor(arb_t A, const fmpq * a, slong alen, const fmpq * b, slong blen, const fmpz_t bden, const arb_t z, slong k, slong prec)
{
    slong i;

    fmpz_t t, u;
    fmpz_init(t);
    fmpz_init(u);

    if (alen == 0)
    {
        if (z == NULL)
            arb_set_fmpz(A, bden);
        else if (fmpz_is_one(bden))
            arb_set(A, z);
        else
            arb_mul_fmpz(A, z, bden, prec);
    }
    else
    {
        /* product of a_i + k = p/q + k = (p + k*q)/q */
        fmpz_mul_ui(t, fmpq_denref(a + 0), k);
        fmpz_add(t, t, fmpq_numref(a + 0));

        for (i = 1; i < alen; i++)
        {
            fmpz_mul_ui(u, fmpq_denref(a + i), k);
            fmpz_add(u, u, fmpq_numref(a + i));
            fmpz_mul(t, t, u);
        }

        if (!fmpz_is_one(bden))
            fmpz_mul(t, t, bden);

        if (z == NULL)
            arb_set_fmpz(A, t);
        else
            arb_mul_fmpz(A, z, t, prec);
    }

    fmpz_clear(t);
    fmpz_clear(u);
}

static void
bsplit(arb_t A1, arb_t B1, arb_t C1,
        const fmpq * a, slong alen, const fmpz_t aden,
        const fmpq * b, slong blen, const fmpz_t bden,
        const arb_t z, int reciprocal,
        slong aa,
        slong bb,
        slong prec)
{
    if (bb - aa == 1)
    {
        factor(A1, a, alen, b, blen, bden, reciprocal ? NULL : z, aa, prec);
        factor(C1, b, blen, a, alen, aden, reciprocal ? z : NULL, aa, prec);
        /* arb_set(B1, C1);   but we skip this */
    }
    else
    {
        slong m;

        arb_t A2, B2, C2;

        arb_init(A2);
        arb_init(B2);
        arb_init(C2);

        m = aa + (bb - aa) / 2;

        bsplit(A1, B1, C1, a, alen, aden, b, blen, bden, z, reciprocal, aa, m, prec);
        bsplit(A2, B2, C2, a, alen, aden, b, blen, bden, z, reciprocal, m, bb, prec);

        if (bb - m == 1)  /* B2 = C2 */
        {
            if (m - aa == 1)
                arb_add(B2, A1, C1, prec);
            else
                arb_add(B2, A1, B1, prec);

            arb_mul(B1, B2, C2, prec);
        }
        else
        {
            if (m - aa == 1)
                arb_mul(B1, C1, C2, prec);
            else
                arb_mul(B1, B1, C2, prec);

            arb_addmul(B1, A1, B2, prec);
        }

        arb_mul(A1, A1, A2, prec);
        arb_mul(C1, C1, C2, prec);

        arb_clear(A2);
        arb_clear(B2);
        arb_clear(C2);
    }
}

void
arb_hypgeom_sum_fmpq_arb_bs(arb_t res, const fmpq * a, slong alen, const fmpq * b, slong blen, const arb_t z, int reciprocal, slong N, slong prec)
{
    arb_t u, v, w;
    fmpz_t aden, bden;
    slong i;

    if (N <= 3)
    {
        arb_hypgeom_sum_fmpq_arb_forward(res, a, alen, b, blen, z, reciprocal, N, prec);
        return;
    }

    arb_init(u);
    arb_init(v);
    arb_init(w);

    fmpz_init(aden);
    fmpz_init(bden);

    fmpz_one(aden);
    fmpz_one(bden);

    for (i = 0; i < alen; i++)
        fmpz_mul(aden, aden, fmpq_denref(a + i));
    for (i = 0; i < blen; i++)
        fmpz_mul(bden, bden, fmpq_denref(b + i));

    /* we compute to N-1 instead of N to avoid dividing by 0 in the
       denominator when computing a hypergeometric polynomial
       that terminates right before a pole */
    bsplit(u, v, w, a, alen, aden, b, blen, bden, z, reciprocal, 0, N - 1, prec);

    arb_add(res, u, v, prec); /* s = s + t */
    arb_div(res, res, w, prec);

    arb_clear(u);
    arb_clear(v);
    arb_clear(w);

    fmpz_clear(aden);
    fmpz_clear(bden);
}

