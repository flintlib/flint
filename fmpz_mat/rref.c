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

    Copyright (C) 2011-2012 Fredrik Johansson
    Copyright (C) 2014 Alex J. Best

******************************************************************************/

#include "fmpz_mat.h"

slong
fmpz_mat_rref(fmpz_mat_t R, fmpz_t den, const fmpz_mat_t A)
{
    if (FLINT_MIN(A->c, A->r) <= 20)
        return fmpz_mat_rref_fflu(R, den, A);
    else if (A->r <= 105 && A->c >= 1.4 * A->r)
        return fmpz_mat_rref_fflu(R, den, A);
    else
        return fmpz_mat_rref_mul(R, den, A);
}

