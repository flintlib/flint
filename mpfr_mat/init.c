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

******************************************************************************/

#include <stdlib.h>
#include <gmp.h>
#include <mpfr.h>
#include "flint.h"
#include "mpfr_mat.h"

void
mpfr_mat_init(mpfr_mat_t mat, long rows, long cols, mpfr_prec_t prec)
{

    if ((rows) && (cols))       /* Allocate space for r*c small entries */
    {
        long i;
        mat->entries =
            (__mpfr_struct *) flint_malloc(rows * cols * sizeof(__mpfr_struct));
        mat->rows = (__mpfr_struct **) flint_malloc(rows * sizeof(__mpfr_struct *));  /* Initialise rows */

        for (i = 0; i < rows * cols; i++)
            mpfr_init2(mat->entries + i, prec);
        for (i = 0; i < rows; i++)
            mat->rows[i] = mat->entries + i * cols;
    }
    else
        mat->entries = NULL;

    mat->r = rows;
    mat->c = cols;
    mat->prec = prec;
}
