/*
    Copyright (C) 2014 Abhinav Baid

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "double_extras.h"
#include "d_mat.h"

void
d_mat_gso(d_mat_t B, const d_mat_t A)
{
    slong i, j, k;
    int flag;
    double t, s;

    if (B->r != A->r || B->c != A->c)
    {
        flint_throw(FLINT_ERROR, "Exception (d_mat_gso). Incompatible dimensions.\n");
    }

    if (B == A)
    {
        d_mat_t t;
        d_mat_init(t, A->r, A->c);
        d_mat_gso(t, A);
        d_mat_swap_entrywise(B, t);
        d_mat_clear(t);
        return;
    }

    if (A->r == 0)
    {
        return;
    }

    for (k = 0; k < A->c; k++)
    {
        for (j = 0; j < A->r; j++)
        {
            d_mat_entry(B, j, k) = d_mat_entry(A, j, k);
        }
        flag = 1;
        while (flag)
        {
            t = 0;
            for (i = 0; i < k; i++)
            {
                s = 0;
                for (j = 0; j < A->r; j++)
                {
                    s += d_mat_entry(B, j, i) * d_mat_entry(B, j, k);
                }
                t += s * s;
                for (j = 0; j < A->r; j++)
                {
                    d_mat_entry(B, j, k) -= s * d_mat_entry(B, j, i);
                }
            }
            s = 0;
            for (j = 0; j < A->r; j++)
            {
                s += d_mat_entry(B, j, k) * d_mat_entry(B, j, k);
            }
            t += s;
            flag = 0;
            if (s < t)
            {
                if (fabs(s * D_EPS) < 1.0e-308)
                    s = 0;
                else
                    flag = 1;
            }
        }
        s = sqrt(s);
        if (s != 0)
            s = 1 / s;
        for (j = 0; j < A->r; j++)
        {
            d_mat_entry(B, j, k) *= s;
        }
    }
}
