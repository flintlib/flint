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

#include <mpir.h>
#include "flint.h"
#include "nmod_vec.h"
#include "nmod_mat.h"

int _nmod_mat_pivot(mp_limb_t ** rows, long n, long start_row, long col)
{
    mp_limb_t * tmp;
    long j;

    if (rows[start_row][col] != 0)
        return 1;

    for (j = start_row + 1; j < n; j++)
    {
        if (rows[j][col] != 0)
        {
            tmp = rows[j];
            rows[j] = rows[start_row];
            rows[start_row] = tmp;
            return -1;
        }
    }
    return 0;
}
