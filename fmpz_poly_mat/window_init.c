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

    Copyright (C) 2014 Abhinav Baid
    Copyright (C) 2015 Elena Sergeicheva

******************************************************************************/

#include "fmpz_poly_mat.h"

void
fmpz_poly_mat_window_init(fmpz_poly_mat_t window, const fmpz_poly_mat_t mat, slong r1,
                     slong c1, slong r2, slong c2)
{
    slong i;
    window->entries = NULL;

    if (r2 - r1)
        window->rows = flint_malloc((r2 - r1) * sizeof(fmpz_poly_t));

    for (i = 0; i < r2 - r1; i++)
        window->rows[i] = mat->rows[r1 + i] + c1;

    window->r = r2 - r1;
    window->c = c2 - c1;
}
