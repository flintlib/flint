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

    Copyright (C) 2011 Fredrik Johansson

******************************************************************************/

#include "flint.h"
#include "nmod_poly.h"
#include "nmod_poly_mat.h"

void
nmod_poly_mat_trace(nmod_poly_t trace, const nmod_poly_mat_t mat)
{
    len_t i, n = nmod_poly_mat_nrows(mat);

    if (n == 0)
        nmod_poly_zero(trace);
    else
    {
        nmod_poly_set(trace, nmod_poly_mat_entry(mat, 0, 0));
        for (i = 1; i < n; i++)
            nmod_poly_add(trace, trace, nmod_poly_mat_entry(mat, i, i));
    }
}
