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

    Copyright (C) 2014 Alex J. Best

******************************************************************************/

#include "fmpz_mat.h"

void
fmpz_mat_hnf(fmpz_mat_t H, const fmpz_mat_t A)
{
    slong m = A->r, b = fmpz_mat_max_bits(A), cutoff = 3;

    if (b <= 2)
        cutoff = 69;
    else if (b <= 4)
        cutoff = 49;
    else if (b <= 8)
        cutoff = 34;
    else if (b <= 16)
        cutoff = 19;
    else if (b <= 32)
        cutoff = 11;
    else if (b <= 64)
        cutoff = 9;
    else if (b <= 128)
        cutoff = 7;
    else if (b <= 256)
        cutoff = 5;

    if (m < cutoff)
        fmpz_mat_hnf_classical(H, A);
    else
        fmpz_mat_hnf_pernet_stein(H, A);
}
