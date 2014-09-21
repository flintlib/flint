/*=============================================================================

    This file is part of FLINT.

    FLINT is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    FLINT is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with FLINT; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301 USA

=============================================================================*/
/******************************************************************************

    Copyright (C) 2014 Abhinav Baid

******************************************************************************/

#include "d_mat.h"

void
d_mat_gso(d_mat_t B, const d_mat_t A)
{
    slong i, j, k, flag;
    double t, s;

    if (B->r != A->r || B->c != A->c)
    {
        flint_printf("Exception (d_mat_gso). Incompatible dimensions.\n");
        abort();
    }

    if (B == A)
    {
        d_mat_t t;
        d_mat_init(t, A->r, A->c);
        d_mat_gso(t, A);
        d_mat_swap(B, t);
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
                if (s * D_EPS == 0)
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
