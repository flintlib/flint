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

#include "mpf_mat.h"

void
mpf_mat_mul(mpf_mat_t C, const mpf_mat_t A, const mpf_mat_t B)
{
    slong ar, bc, br;
    slong i, j, k;
    mpf_t tmp;

    ar = A->r;
    br = B->r;
    bc = B->c;

    if (C == A || C == B)
    {
        mpf_mat_t t;
        mpf_mat_init(t, ar, bc, C->prec);
        mpf_mat_mul(t, A, B);
        mpf_mat_swap(C, t);
        mpf_mat_clear(t);
        return;
    }

    if (C->r != ar || C->c != bc)
    {
        flint_printf("Exception (mpf_mat_mul). Incompatible dimensions.\n");
        abort();
    }

    if (br == 0)
    {
        mpf_mat_zero(C);
        return;
    }

    mpf_init2(tmp, C->prec);

    for (i = 0; i < ar; i++)
    {
        for (j = 0; j < bc; j++)
        {
            mpf_mul(mpf_mat_entry(C, i, j), mpf_mat_entry(A, i, 0),
                    mpf_mat_entry(B, 0, j));

            for (k = 1; k < br; k++)
            {
                mpf_mul(tmp, mpf_mat_entry(A, i, k), mpf_mat_entry(B, k, j));
                mpf_add(mpf_mat_entry(C, i, j), mpf_mat_entry(C, i, j), tmp);
            }
        }
    }

    mpf_clear(tmp);
}
