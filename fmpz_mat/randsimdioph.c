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

    Copyright (C) 2005-2009 Damien Stehle
    Copyright (C) 2007 David Cade
    Copyright (C) 2010 William Hart

******************************************************************************/

#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_mat.h"

void
fmpz_mat_randsimdioph(fmpz_mat_t mat, flint_rand_t state, mp_bitcnt_t bits, mp_bitcnt_t bits2)
{
    const long c = mat->c, r = mat->r;

    long i, j;

    if (c != r)
    {
        printf("Exception (fmpz_mat_randsimdioph). Ill-formed matrix.\n");
        abort();
    }

    fmpz_one(mat->rows[0]);
    fmpz_mul_2exp(mat->rows[0], mat->rows[0], bits2);
    for (j = 1; j < c; j++)
        fmpz_randbits(mat->rows[0] + j, state, bits);
    for (i = 1; i < r; i++)
    {
        for (j = 0; j < i; j++)
            fmpz_zero(mat->rows[i] + j);
        fmpz_one(mat->rows[i] + i);
        fmpz_mul_2exp(mat->rows[i] + i, mat->rows[i] + i, bits);
        for (j = i + 1; j < c; j++)
            fmpz_zero(mat->rows[i] + j);
    }
}
