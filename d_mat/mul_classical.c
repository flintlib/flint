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

    Copyright (C) 2010,2011 Fredrik Johansson
    Copyright (C) 2014 Abhinav Baid

******************************************************************************/

#include "d_mat.h"

void
d_mat_mul_classical(d_mat_t C, const d_mat_t A, const d_mat_t B)
{
    slong ar, bc, br;
    slong jj, kk, i, j, k, blocksize;
    double temp;
    d_mat_t Bt;

    ar = A->r;
    br = B->r;
    bc = B->c;
    blocksize = 64 / sizeof(double);

    if (C == A || C == B)
    {
        d_mat_t t;
        d_mat_init(t, ar, bc);
        d_mat_mul_classical(t, A, B);
        d_mat_swap(C, t);
        d_mat_clear(t);
        return;
    }

    if (C->r != ar || C->c != bc)
    {
        flint_printf
            ("Exception (d_mat_mul_classical). Incompatible dimensions.\n");
        abort();
    }

    if (br == 0)
    {
        d_mat_zero(C);
        return;
    }

    d_mat_init(Bt, bc, br);
    d_mat_transpose(Bt, B);
    d_mat_zero(C);

    for (jj = 0; jj < bc; jj += blocksize)
    {
        for (kk = 0; kk < br; kk += blocksize)
        {
            for (i = 0; i < ar; i++)
            {
                for (j = jj; j < FLINT_MIN(jj + blocksize, bc); j++)
                {
                    temp = 0;

                    for (k = kk; k < FLINT_MIN(kk + blocksize, br); k++)
                    {
                        temp += d_mat_entry(A, i, k) * d_mat_entry(Bt, j, k);
                    }
                    d_mat_entry(C, i, j) += temp;
                }
            }
        }
    }
    d_mat_clear(Bt);
}
