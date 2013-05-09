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
#include "fmpz_mat.h"
#include "fmpz_vec.h"


void
fmpz_mat_randrank(fmpz_mat_t mat, flint_rand_t state, len_t rank,
                  mp_bitcnt_t bits)
{
    len_t i;
    fmpz * diag;

    if (rank < 0 || rank > mat->r || rank > mat->c)
    {
        printf("Exception (fmpz_mat_randrank). Impossible rank.\n");
        abort();
    }

    diag = _fmpz_vec_init(rank);
    for (i = 0; i < rank; i++)
        fmpz_randtest_not_zero(&diag[i], state, bits);

    fmpz_mat_randpermdiag(mat, state, diag, rank);

    _fmpz_vec_clear(diag, rank);
}
