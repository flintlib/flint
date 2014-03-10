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

    Copyright (C) 2014 Ashish Kedia

******************************************************************************/

#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "nmod_mat.h"

void
nmod_mat_swap(nmod_mat_t mat1, nmod_mat_t mat2)
{
    nmod_mat_t temp;
    *temp = *mat1;
    *mat1 = *mat2;
    *mat2 = *temp;
}
