/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ca_mat.h"

void
ca_mat_pow_ui_binexp(ca_mat_t B, const ca_mat_t A, ulong exp, ca_ctx_t ctx)
{
    slong d = ca_mat_nrows(A);

    if (exp <= 2 || d <= 1)
    {
        if (exp == 0 || d == 0)
        {
            ca_mat_one(B, ctx);
        }
        else if (d == 1)
        {
            ca_pow_ui(ca_mat_entry(B, 0, 0),
                 ca_mat_entry(A, 0, 0), exp, ctx);
        }
        else if (exp == 1)
        {
            ca_mat_set(B, A, ctx);
        }
        else if (exp == 2)
        {
            ca_mat_sqr(B, A, ctx);
        }
    }
    else
    {
        ca_mat_t T, U;
        slong i;

        ca_mat_init(T, d, d, ctx);
        ca_mat_set(T, A, ctx);
        ca_mat_init(U, d, d, ctx);

        for (i = ((slong) FLINT_BIT_COUNT(exp)) - 2; i >= 0; i--)
        {
            ca_mat_sqr(U, T, ctx);

            if (exp & (WORD(1) << i))
                ca_mat_mul(T, U, A, ctx);
            else
                ca_mat_swap(T, U, ctx);
        }

        ca_mat_swap(B, T, ctx);
        ca_mat_clear(T, ctx);
        ca_mat_clear(U, ctx);
    }
}
