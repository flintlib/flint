/*
    Copyright (C) 2012 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb_mat.h"

void
arb_mat_pow_ui(arb_mat_t B, const arb_mat_t A, ulong exp, slong prec)
{
    slong d = arb_mat_nrows(A);

    if (exp <= 2 || d <= 1)
    {
        if (exp == 0 || d == 0)
        {
            arb_mat_one(B);
        }
        else if (d == 1)
        {
            arb_pow_ui(arb_mat_entry(B, 0, 0),
                 arb_mat_entry(A, 0, 0), exp, prec);
        }
        else if (exp == 1)
        {
            arb_mat_set(B, A);
        }
        else if (exp == 2)
        {
            arb_mat_sqr(B, A, prec);
        }
    }
    else
    {
        arb_mat_t T, U;
        slong i;

        arb_mat_init(T, d, d);
        arb_mat_set(T, A);
        arb_mat_init(U, d, d);

        for (i = ((slong) FLINT_BIT_COUNT(exp)) - 2; i >= 0; i--)
        {
            arb_mat_sqr(U, T, prec);

            if (exp & (WORD(1) << i))
                arb_mat_mul(T, U, A, prec);
            else
                arb_mat_swap(T, U);
        }

        arb_mat_swap(B, T);
        arb_mat_clear(T);
        arb_mat_clear(U);
    }
}

