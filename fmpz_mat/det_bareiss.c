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

******************************************************************************/

#include <stdlib.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_mat.h"


void
_fmpz_mat_det_bareiss(fmpz_t det, fmpz_mat_t tmp)
{
    long m = tmp->r;
    long rank;

    rank = _fmpz_mat_rowreduce(tmp, ROWREDUCE_FAST_ABORT);

    if (rank < 0)
    {
        rank = -rank;
        if (rank < m)
            fmpz_zero(det);
        else
            fmpz_neg(det, &tmp->rows[rank-1][rank-1]);
    }
    else
    {
        if (rank < m)
            fmpz_zero(det);
        else
            fmpz_set(det, &tmp->rows[rank-1][rank-1]);
    }
}

void
fmpz_mat_det_bareiss(fmpz_t det, const fmpz_mat_t A)
{
    fmpz_mat_t tmp;

    if (A->r < 1)
    {
        fmpz_set_ui(det, 1UL);
        return;
    }

    fmpz_mat_init_set(tmp, A);
    _fmpz_mat_det_bareiss(det, tmp);
    fmpz_mat_clear(tmp);
}
