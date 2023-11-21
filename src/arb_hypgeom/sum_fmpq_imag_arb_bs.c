/*
    Copyright (C) 2021 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb.h"
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
bsplit(acb_t A1, acb_t B1, acb_t C1,
        const fmpq * a, slong alen, const fmpz_t aden,
        const fmpq * b, slong blen, const fmpz_t bden,
        const arb_t z, int reciprocal,
        slong aa,
        slong bb,
        slong prec)
{
    if (bb - aa == 1)
    {
        factor(acb_realref(A1), a, alen, b, blen, bden, reciprocal ? NULL : z, aa, prec);
        factor(acb_realref(C1), b, blen, a, alen, aden, reciprocal ? z : NULL, aa, prec);
        arb_zero(acb_imagref(A1));
        arb_zero(acb_imagref(C1));

        if (reciprocal)
            acb_div_onei(C1, C1);
        else
            acb_mul_onei(A1, A1);

        /* arb_set(B1, C1);   but we skip this */
    }
    else
    {
        slong m;

        acb_t A2, B2, C2;

        acb_init(A2);
        acb_init(B2);
        acb_init(C2);

        m = aa + (bb - aa) / 2;

        bsplit(A1, B1, C1, a, alen, aden, b, blen, bden, z, reciprocal, aa, m, prec);
        bsplit(A2, B2, C2, a, alen, aden, b, blen, bden, z, reciprocal, m, bb, prec);

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
arb_hypgeom_sum_fmpq_imag_arb_bs(arb_t res_real, arb_t res_imag, const fmpq * a, slong alen, const fmpq * b, slong blen, const arb_t z, int reciprocal, slong N, slong prec)
{
    acb_t u, v, w;
    fmpz_t aden, bden;
    slong i;

    if (N <= 3)
    {
        arb_hypgeom_sum_fmpq_imag_arb_forward(res_real, res_imag, a, alen, b, blen, z, reciprocal, N, prec);
        return;
    }

    acb_init(u);
    acb_init(v);
    acb_init(w);

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

    acb_add(u, u, v, prec); /* s = s + t */
    acb_div(u, u, w, prec);

    if (!acb_is_finite(u))
        acb_indeterminate(u);

    arb_swap(res_real, acb_realref(u));
    arb_swap(res_imag, acb_imagref(u));

    acb_clear(u);
    acb_clear(v);
    acb_clear(w);

    fmpz_clear(aden);
    fmpz_clear(bden);
}
