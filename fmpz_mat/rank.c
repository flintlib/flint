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


long
fmpz_mat_rank(const fmpz_mat_t A)
{
    long m, n, rank;
    fmpz_mat_t tmp;

    m = A->r;
    n = A->c;

    if (m < 1 || n < 1)
        return 0;

    fmpz_mat_init(tmp, m, n);
    _fmpz_vec_copy(tmp->entries, A->entries, m * n);
    rank = _fmpz_mat_rowreduce(tmp->rows, m, n, ROWREDUCE_ECHELON_FORM);
    fmpz_mat_clear(tmp);
    return FLINT_ABS(rank);
}
