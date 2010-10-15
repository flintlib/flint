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

    Copyright (C) 2008-2009 William Hart
    Copyright (C) 2010 Fredrik Johansson

******************************************************************************/

#include <stdlib.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_mat.h"

void
fmpz_mat_copy(fmpz_mat_t mat1, fmpz_mat_t mat2)
{
    if (mat1 != mat2)
    {
        if (mat1->r != mat2->r || mat1->c != mat2->c)
        {
            printf("exception: fmpz_mat_copy: incompatible dimensions\n");
            abort();
        }

        _fmpz_vec_copy(mat1->entries, mat2->entries, mat2->r * mat2->c);
    }
}
