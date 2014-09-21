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
d_mat_qr(d_mat_t Q, d_mat_t R, const d_mat_t A)
{
    slong i, j, k, flag, orig;
    double t, s;

    if (Q->r != A->r || Q->c != A->c || R->r != A->c || R->c != A->c)
    {
        flint_printf("Exception (d_mat_qr). Incompatible dimensions.\n");
        abort();
    }

    if (Q == A)
    {
        d_mat_t t;
        d_mat_init(t, A->r, A->c);
        d_mat_qr(t, R, A);
        d_mat_swap(Q, t);
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
            d_mat_entry(Q, j, k) = d_mat_entry(A, j, k);
        }
        orig = flag = 1;
        while (flag)
        {
            t = 0;
            for (i = 0; i < k; i++)
            {
                s = 0;
                for (j = 0; j < A->r; j++)
                {
                    s += d_mat_entry(Q, j, i) * d_mat_entry(Q, j, k);
                }
                if (orig)
                {
                    d_mat_entry(R, i, k) = s;
                }
                else
                {
                    d_mat_entry(R, i, k) += s;
                }
                t += s * s;
                for (j = 0; j < A->r; j++)
                {
                    d_mat_entry(Q, j, k) -= s * d_mat_entry(Q, j, i);
                }
            }
            s = 0;
            for (j = 0; j < A->r; j++)
            {
                s += d_mat_entry(Q, j, k) * d_mat_entry(Q, j, k);
            }
            t += s;
            flag = 0;
            if (s < t)
            {
                orig = 0;
                if (s * D_EPS == 0)
                    s = 0;
                else
                    flag = 1;
            }
        }
        d_mat_entry(R, k, k) = s = sqrt(s);
        if (s != 0)
            s = 1 / s;
        for (j = 0; j < A->r; j++)
        {
            d_mat_entry(Q, j, k) *= s;
        }
    }
}
