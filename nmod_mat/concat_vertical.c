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

    Copyright (C) 2015 Elena Sergeicheva

******************************************************************************/


#include "nmod_mat.h"

void
nmod_mat_concat_vertical(nmod_mat_t res, const nmod_mat_t mat1, const nmod_mat_t mat2)
{
    slong i;
    slong r1 = mat1->r;
    slong c1 = mat1->c;
    slong r2 = mat2->r;

    for (i = 0; i < r1; i++)
    	flint_mpn_copyi(res->rows[i], mat1->rows[i], c1);

    for (i = 0; i < r2; i++)
    	flint_mpn_copyi(res->rows[i + r1], mat2->rows[i], c1);

}
