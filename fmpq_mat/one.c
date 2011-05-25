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

    Copyright (C) 2011 Sebastian Pancratz

******************************************************************************/

#include <mpir.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpq.h"
#include "fmpq_mat.h"

void fmpq_mat_one(fmpq_mat_t mat)
{
    long i, j;

    if (mat->r != mat->c)
    {
        printf("Exception:  matrix not square in fmpq_mat_one.\n");
        abort();
    }

    for (i = 0; i < mat->r; i++)
        for (j = 0; j < mat->c; j++)
            fmpq_zero(fmpq_mat_entry(mat, i, j));
    for (i = 0; i < mat->r; i++)
        fmpq_one(fmpq_mat_entry(mat, i, i));
}

