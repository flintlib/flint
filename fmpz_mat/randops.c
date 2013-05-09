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
#include "ulong_extras.h"


void
fmpz_mat_randops(fmpz_mat_t mat, flint_rand_t state, len_t count)
{
    len_t c, i, j, k;
    len_t m = mat->r;
    len_t n = mat->c;

    if (mat->r == 0 || mat->c == 0)
        return;

    for (c = 0; c < count; c++)
    {
        if (n_randint(state, 2))
        {
            if ((i = n_randint(state, m)) == (j = n_randint(state, m)))
                continue;
            if (n_randint(state, 2))
                for (k = 0; k < n; k++)
                    fmpz_add(&mat->rows[j][k], &mat->rows[j][k],
                        &mat->rows[i][k]);
            else
                for (k = 0; k < n; k++)
                    fmpz_sub(&mat->rows[j][k], &mat->rows[j][k],
                        &mat->rows[i][k]);
        }
        else
        {
            if ((i = n_randint(state, n)) == (j = n_randint(state, n)))
                continue;
            if (n_randint(state, 2))
                for (k = 0; k < m; k++)
                    fmpz_add(&mat->rows[k][j], &mat->rows[k][j],
                        &mat->rows[k][i]);
            else
                for (k = 0; k < m; k++)
                    fmpz_sub(&mat->rows[k][j], &mat->rows[k][j],
                        &mat->rows[k][i]);
        }
    }
}
