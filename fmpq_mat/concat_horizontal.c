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

    Copyright (C) 2015 Elena Sergeicheva

******************************************************************************/




#include "fmpq_mat.h"

void
fmpq_mat_concat_horizontal(fmpq_mat_t res, const fmpq_mat_t mat1, const fmpq_mat_t mat2)
{
    slong i, j;
    slong r1 = mat1->r;
    slong c1 = mat1->c;
    slong r2 = mat2->r;
    slong c2 = mat2->c;
    
    for (i = 0; i < r1; i++)
    {
        for (j = 0; j < c1; j++)
        {
            fmpq_set(fmpq_mat_entry(res, i, j), fmpq_mat_entry(mat1, i, j));
        }
    }

    for (i = 0; i < r2; i++)
    {
        for (j = 0; j < c2; j++)
        {
            fmpq_set(fmpq_mat_entry(res, i, j + c1), fmpq_mat_entry(mat2, i, j));
        }
    }
}
