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

#include "flint.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_mat.h"

long
fmpz_mat_rref_fraction_free(long * perm, fmpz_mat_t B, fmpz_t den,
    const fmpz_mat_t A)
{
    long i, rank;

    if (B != A)
        fmpz_mat_set(B, A);

    rank = _fmpz_mat_rowreduce(perm, B, ROWREDUCE_FULL);
    rank = FLINT_ABS(rank);

    if (rank == 0)
    {
        fmpz_set_ui(den, 1UL);
    }
    else
    {
        for (i = 0; i < B->c; i++)
        {
            if (!fmpz_is_zero(fmpz_mat_entry(B, 0, i)))
            {
                fmpz_set(den, fmpz_mat_entry(B, 0, i));
                break;
            }
        }
    }

    return rank;
}
