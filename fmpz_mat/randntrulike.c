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
#include "fmpz_vec.h"
#include "fmpz_mat.h"

void
fmpz_mat_randntrulike(fmpz_mat_t mat, flint_rand_t state, mp_bitcnt_t bits, ulong q)
{
    const long c = mat->c, r = mat->r, d = r / 2;

    long i, j, k;
    fmpz *h;

    if ((c != r) || (c != 2 * d))
    {
        printf("Exception (fmpz_mat_randntrulike). Ill-formed matrix.\n");
        abort();
    }

    h = _fmpz_vec_init(d);

    for (i = 0; i < d; i++)
        fmpz_randbits(h + i, state, bits);

    for (i = 0; i < d; i++)
    {
        for (j = 0; j < i; j++)
            fmpz_zero(mat->rows[i] + j);
        fmpz_one(mat->rows[i] + i);
        for (j = i + 1; j < d; j++)
            fmpz_zero(mat->rows[i] + j);
    }

    for (i = d; i < r; i++)
        for (j = 0; j < d; j++)
            fmpz_zero(mat->rows[i] + j);

    for (i = d; i < r; i++)
    {
        for (j = d; j < i; j++)
            fmpz_zero(mat->rows[i] + j);
        fmpz_set_ui(mat->rows[i] + i, q);
        for (j = i + 1; j < c; j++)
            fmpz_zero(mat->rows[i] + j);
    }

    for (i = 0; i < d; i++)
    {
        for (j = d; j < c; j++)
        {
            k = j + i;
            while (k >= d)
                k -= d;
            fmpz_set(mat->rows[i] + j, h + k);
        }
    }

    _fmpz_vec_clear(h, d);
}
