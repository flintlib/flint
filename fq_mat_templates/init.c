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

    Copyright (C) 2010 William Hart
    Copyright (C) 2013 Mike Hansen

******************************************************************************/


#ifdef T

#include "templates.h"

void
TEMPLATE(T, mat_init) (TEMPLATE(T, mat_t) mat, slong rows, slong cols,
                       const TEMPLATE(T, ctx_t) ctx)
{
    if ((rows) && (cols))       /* Allocate space for r*c small entries */
    {
        slong i, j;
        mat->entries = flint_malloc(rows * cols * sizeof(TEMPLATE(T, struct)));
        mat->rows = flint_malloc(rows * sizeof(TEMPLATE(T, struct) *)); /* Initialise rows */

        for (i = 0; i < rows; i++)
        {
            mat->rows[i] = mat->entries + i * cols;
            for (j = 0; j < cols; j++)
            {
                TEMPLATE(T, init) (mat->rows[i] + j, ctx);
            }
        }
    }
    else
        mat->entries = NULL;

    mat->r = rows;
    mat->c = cols;
}


#endif
