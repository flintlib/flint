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
fmpz_mat_randintrel(fmpz_mat_t mat, flint_rand_t state, mp_bitcnt_t bits)
{
    const len_t c = mat->c, r = mat->r;

    len_t i, j;

    if (c != r + 1)
    {
        printf("Exception (fmpz_mat_randintrel).  c != r + 1.\n");
        abort();
    }

    for (i = 0; i < r; i++)
    {
        fmpz_randbits(mat->rows[i], state, bits);
        for (j = 1; j <= i; j++)
            fmpz_zero(mat->rows[i] + j);
        fmpz_one(mat->rows[i] + i + 1);
        for (j = i + 2; j < c; j++)
            fmpz_zero(mat->rows[i] + j);
    }
}
