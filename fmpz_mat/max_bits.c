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

******************************************************************************/

#include <stdlib.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_mat.h"

len_t
fmpz_mat_max_bits(const fmpz_mat_t mat)
{
    len_t i, bits, row_bits, sign;

    sign = 1;
    bits = 0;

    if (mat->r == 0 || mat->c == 0)
        return 0;

    for (i = 0; i < mat->r; i++)
    {
        row_bits = _fmpz_vec_max_bits(mat->rows[i], mat->c);
        if (row_bits < 0)
        {
            row_bits = -row_bits;
            sign = -1;
        }
        bits = FLINT_MAX(bits, row_bits);
    }

    return bits * sign;
}
