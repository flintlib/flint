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

    Copyright (C) 2010 Fredrik Johansson
    Copyright (C) 2014 Abhinav Baid

******************************************************************************/

#include "d_mat.h"

void
d_mat_transpose(d_mat_t B, const d_mat_t A)
{
    slong ii, jj, i, j, blocksize;
    blocksize = 64 / sizeof(double);

    if (B->r != A->c || B->c != A->r)
    {
        flint_printf
            ("Exception (d_mat_transpose). Incompatible dimensions.\n");
        abort();
    }

    if (B == A)
    {
        d_mat_t t;
        d_mat_init(t, A->r, A->c);
        d_mat_transpose(t, A);
        d_mat_swap(B, t);
        d_mat_clear(t);
        return;
    }

    for (ii = 0; ii < B->r; ii += blocksize)
        for (jj = 0; jj < B->c; jj += blocksize)
            for (i = ii; i < FLINT_MIN(ii + blocksize, B->r); i++)
                for (j = jj; j < FLINT_MIN(jj + blocksize, B->c); j++)
                    d_mat_entry(B, i, j) = d_mat_entry(A, j, i);

}
