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

    Copyright (C) 2011 Fredrik Johansson

******************************************************************************/

#include "fmpq_mat.h"

void fmpq_mat_print(const fmpq_mat_t mat)
{
    slong i, j;

    flint_printf("<%wd x %wd matrix over Q>\n", mat->r, mat->c);

    for (i = 0; i < mat->r; i++)
    {
        flint_printf("[");
        for (j = 0; j < mat->c; j++)
        {
            fmpq_print(fmpq_mat_entry(mat, i, j));
            if (j + 1 < mat->c)
                flint_printf(", ");
        }
        flint_printf("]\n");
    }
    flint_printf("\n");
}
