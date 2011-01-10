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
#include "ulong_extras.h"
#include "nmod_mat.h"
#include "nmod_vec.h"



void
nmod_mat_randops(nmod_mat_t mat, long count)
{
    long c, i, j, k;
    long m = mat->r;
    long n = mat->c;

    if (mat->r == 0 || mat->c == 0)
        return;

    for (c = 0; c < count; c++)
    {
        if (n_randint(2, NULL))
        {
            if ((i = n_randint(m, NULL)) == (j = n_randint(m, NULL)))
                continue;

            if (n_randint(2, NULL))
                for (k = 0; k < n; k++)
                    mat->rows[j][k] = nmod_add(mat->rows[j][k],
                        mat->rows[i][k], mat->mod);
            else
                for (k = 0; k < n; k++)
                    mat->rows[j][k] = nmod_sub(mat->rows[j][k],
                        mat->rows[i][k], mat->mod);
        }
        else
        {
            if ((i = n_randint(n, NULL)) == (j = n_randint(n, NULL)))
                continue;
            if (n_randint(2, NULL))
                for (k = 0; k < m; k++)
                    mat->rows[k][j] = nmod_add(mat->rows[k][j],
                        mat->rows[k][i], mat->mod);
            else
                for (k = 0; k < m; k++)
                    mat->rows[k][j] = nmod_sub(mat->rows[k][j],
                        mat->rows[k][i], mat->mod);
        }
    }
}
