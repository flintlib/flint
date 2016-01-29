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

    Copyright (C) 2016 Aaditya Thakkar

******************************************************************************/

#include<stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "fmpz_mat.h"
#include "fmpz.h"
#include "fmpz_vec.h"

void
fmpz_mat_addmul(fmpz_mat_t D, const fmpz_mat_t C,
                                const fmpz_mat_t A, const fmpz_mat_t B)
{
    slong m, k, n;

    m = A->r;
    k = A->c;
    n = B->c;

    fmpz_mat_t tmp;
    fmpz_mat_init(tmp, m, n);
    fmpz_mat_mul_strassen(tmp, A, B);
    fmpz_mat_add(D, C, tmp); 
    fmpz_mat_clear(tmp);
    
}
